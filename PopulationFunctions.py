import distances
from scipy import interpolate
import pickle
import numpy as np
import math
import indexTricks as iT
from astropy.io import fits
import pandas as pd


class RedshiftDependentRelation():
    def __init__(self, D=None, #reset=False,
                 cosmo=[0.3,0.7,0.7]):
        self.beginRedshiftDependentRelation(D, cosmo=cosmo)

    def beginRedshiftDependentRelation(self, D, zmax=20, cosmo=[0.3,0.7,0.7]): # NH: increased zmax from 10 to 20
        self.zmax = zmax
        self.zbins, self.dz = np.linspace(0, self.zmax, 401, retstep=True)
        self.z2bins, self.dz2 = np.linspace(0, self.zmax, 201, retstep=True)
        if D==None:
            import distances
            D=distances.Distance(cosmo=cosmo)
        self.D=D

        splinedump = open("redshiftsplines.pkl","rb")
        self.Da_spline, self.Dmod_spline, self.volume_spline, self.Da_bispline = pickle.load(splinedump,encoding='iso-8859-1')

    def Volume(self,z1,z2=None):
        if z2==None:
            return self.splev(z1,self.volume_spline)
        else:
            z1,z2=self.biassert(z1,z2)
            return self.splev(z2,self.volume_spline)-self.splev(z1,self.volume_spline)

    def Da(self,z1,z2=None,units="Mpc"):
        if units=="kpc":
            corfrac=1000
        elif units=="Mpc":
            corfrac=1
        else:
            print("don't know those units yet")
        if z2 is None:
            return self.splev(z1,self.Da_spline)*corfrac
        else:
            z1,z2=self.biassert(z1,z2)
            return self.Da_bispline.ev(z1,z2)*corfrac

    def Dmod(self,z):
        return self.splev(z,self.Dmod_spline)

    def splev(self,x,f_of_x_as_spline):
        return interpolate.splev(x,f_of_x_as_spline)

    def bisplev(self,x,y,f_ofxy_as_bispline):
        return interpolate.bisplev(x,y,f_ofxy_as_bispline)

    def biassert(self,z1,z2):
            try: len(z1)
            except TypeError:z1=[z1]
            try: len(z2)
            except TypeError:z2=[z2]
            if len(z1)==1 and len(z2)!=1:z1=np.ones(len(z2))*z1[0]
            if len(z2)==1 and len(z1)!=1:z2=np.ones(len(z1))*z2[0]
            assert len(z1)==len(z2),"get it together"
            return z1,z2

#====================================================================================


class EinsteinRadiusTools(RedshiftDependentRelation):
    def  __init__(self,D=None):
        self.beginRedshiftDependentRelation(D)
        self.c=299792

    def sie_sig(self,rein,zl,zs):
        self.c=299792
        ds=self.Da(zs)
        print('ds')
        print(ds)
        dls=self.Da(zl,zs)
        sig=(rein*(ds*self.c**2)/(206265*4*math.pi*dls))**0.5
        return sig
    def sie_rein(self,sig,zl,zs):
        self.c=299792
        ds=self.Da(zs)
        dls=self.Da(zl,zs)
        rein=sig**2*((ds*self.c**2)/(206265*4*math.pi*dls))**-1
        rein[rein<0]=0
        return rein


#====================================================================================
class Population(RedshiftDependentRelation):
    def  __init__(self):
        pass

    def draw_apparent_magnitude(self,M,z,band=None,colours=None):
        if band!=None:
            colours=self.colour(z,band)
        if colours is None:
            colours=0
            print("warning no k-correction")
        Dmods=self.Dmod(z)
        ml = M - colours + Dmods
        return ml

    def draw_apparent_size(self,r_phys,z):
        rl = r_phys/(self.Da(z,units="kpc"))
        rl *= 206264
        return rl

#====================================================================================

class LensPopulation_(Population):
    def  __init__(self, zlmax=3, sigfloor=100, D=None, #reset=True,
                  bands=[#'F814W_ACS','g_SDSS','r_SDSS','i_SDSS','z_SDSS','Y_UKIRT','VIS',
                  'JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W'],
                  cosmo=[0.3,0.7,0.7]
                  ): #sadface
        self.sigfloor=sigfloor
        self.zlmax=zlmax
        self.bands=bands

        self.beginRedshiftDependentRelation(D)
        self.beginLensPopulation(D)

    def beginLensPopulation(self,D):
        self.lenspopfunctions()


    def lenspopfunctions(self):
        self.Psigzspline()
        self.Colourspline()
        self.lensPopSplineDump()

    def Psigzspline(self):
        # drawing from a 2d pdf is a pain; should probably make this into its own module

        self.zlbins, self.dzl = np.linspace(0, self.zlmax, 201, retstep=True) # get an array of redshifts and return the step size too (dzl)

        sigmas = np.linspace(self.sigfloor, 400, 401) # get an array of velocity dispersions

        # self.sigbins = sigmas # NH: never used?

        # initialise empty arrays for filling below
        Csiggivenz = np.zeros((sigmas.size, self.zlbins.size))

        CDFbins = np.linspace(0, 1, 1001)

        siggivenCz = np.zeros((CDFbins.size,self.zlbins.size))

        dNdz = self.zlbins*0

        for i in range(len(self.zlbins)): # for each of the redshift bins

            z = self.zlbins[i] # get the current redshift

            dphidsiggivenz = self.phi(sigmas, z) # compute d phi/ d sigma for this redshift

            phisigspline = interpolate.splrep(sigmas, dphidsiggivenz) # interpolate d phi/ d sigma as a function of sigma

            tot = interpolate.splint(self.sigfloor, 500, phisigspline) # integrate the d phi/ d sigma spline between the limits sigfloor and 500

            Csiggivenz[:,i] = np.cumsum(dphidsiggivenz)/np.sum(dphidsiggivenz) # write the cumulative sum of d phi/ d sigma divided by its sum (why this quantity idk) to the empty array

            Csiggivenzspline = interpolate.splrep(Csiggivenz[:,i], sigmas) # interpolate the cumulative sum as a function of sigma

            siggivenCz[:,i] = interpolate.splev(CDFbins, Csiggivenzspline) # splev evaluates the spline Csiggivenzspline at the points given by CDFbins, which are then written to the empty array

            if z!=0:
                dNdz[i] = tot*(self.Volume(z) - self.Volume(z-self.dzl))/self.dzl # multiply the total number (from the integral of d phi/ d sigma) by the comoving volume at the given redshift

        Nofzcdf = np.cumsum(dNdz)/np.sum(dNdz) # get the cumulative distribution function for N(z)

        self.cdfdNdzasspline = interpolate.splrep(Nofzcdf, self.zlbins) # interpolate the N(z) CDF

        self.dNdzspline = interpolate.splrep(self.zlbins, dNdz) # interpolate dN/dz

        # N = interpolate.splint(0, self.zlmax, self.dNdzspline) # integrate the dN/dz spine to get N, the total number of deflectors

        self.cdfdsigdzasspline = interpolate.RectBivariateSpline(CDFbins, self.zlbins, siggivenCz)

        dphidsiggivenz0 = self.phi(sigmas, sigmas*0)

        cdfdNdsigz0 = dphidsiggivenz0.cumsum()/dphidsiggivenz0.sum()

        self.cdfdNdsigz0asspline = interpolate.splrep(cdfdNdsigz0, sigmas) # here

    def Colourspline(self):
        from stellarpop import tools
        sed = tools.getSED('BC_Z=1.0_age=10.00gyr')
        #different SEDs don't change things much

        rband=tools.filterfromfile('r_SDSS')
        z=self.zlbins
        self.colourspline={}
        for band in self.bands:
          if band!="VIS":
            c=z*0
            Cband=tools.filterfromfile(band)
            for i in range(len(z)):
                c[i] = - (tools.ABFM(Cband,sed,z[i]) - tools.ABFM(rband,sed,0))
            self.colourspline[band]=interpolate.splrep(z,c)

    def lensPopSplineDump(self):
        splinedump=open("lenspopsplines.pkl","wb")
        pickle.dump([self.cdfdNdzasspline,self.cdfdNdsigz0asspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,self.zlmax,self.sigfloor,self.colourspline,self.bands],splinedump,2)

    def draw_z(self,N):
        return interpolate.splev(np.random.random(N),self.cdfdNdzasspline)

    def draw_sigma(self,z):
        try: len(z)
        except TypeError:z=[z]
        if self.nozdependence:
            sigs = interpolate.splev(np.random.random(len(z)),self.cdfdNdsigz0asspline)
            return sigs
        else:
            print("Warning: drawing from 2dpdf is low accuracy")
            return self.cdfdsigdzasspline.ev(np.random.random(len(z)),z)

    def draw_zsig(self, N):
        z=self.draw_z(N)
        sig=self.draw_sigma(z)
        return z, sig

    def EarlyTypeRelations(self,sigma,z=None,scatter=True,band=None):#z dependence not encoded currently
        #Hyde and Bernardi, M = r band absolute magnitude.
        V=np.log10(sigma)
        Mr=(-0.37+(0.37**2-(4*(0.006)*(2.97+V)))**0.5)/(2*0.006)
        if scatter:
            Mr+=np.random.randn(len(Mr))*(0.15/2.4)

        #R=4.72+0.63*Mr+0.02*Mr**2 #rest-frame R_band size.
        R=2.46-2.79*V+0.84*V**2
        if scatter:
            R+=np.random.randn(len(R))*0.11

        #convert to observed r band size;
        r_phys = 10**R

        return Mr,r_phys

    def colour(self,z,band):
        return interpolate.splev(z, self.colourspline[band])

    def Ndeflectors(self, z, zmin=0, fsky=1):
        if zmin > z: # make sure the redshifts are passed in the right order
            z, zmin = zmin, z
        N = interpolate.splint(zmin, z, self.dNdzspline)
        N *= fsky
        return N

    def phi(self, sigma, z):
        # equation 3 of Tom's paper
        sigma[sigma==0]+=1e-6
        phi_star=(8*10**-3)*self.D.h**3
        alpha=2.32
        beta=2.67
        sigst=161
        phi=phi_star * \
            ((sigma*1./sigst)**alpha)*\
            np.exp(-(sigma*1./sigst)**beta)*beta/\
            math.gamma(alpha*1./beta)/\
            (1.*sigma)

        phi*=(1+z)**(-2.5)
        return phi

    def draw_flattening_lenspop(self,sigma,z=None):
        # equation 4 of Tom's paper
        x=sigma
        y=0.378-0.000572*x
        e=np.random.rayleigh(y)
        q=1-e
        #dont like ultraflattened masses:
        while len(q[q<0.2])>0 or len(q[q>1])>0:
            q[q<0.2]=1-np.random.rayleigh(y[q<0.2])
            q[q>1]=1-np.random.rayleigh(y[q>1])
        return q

    def drawLensPopulation(self,number):
        self.zl,self.sigl=self.draw_zsig(number)
        self.ql=self.draw_flattening_lenspop(self.sigl)
        self.Mr,self.r_phys_nocol=self.EarlyTypeRelations(self.sigl,self.zl,scatter=True)
        self.ml={}
        self.rl={}
        self.r_phys={}
        for band in self.bands:
            self.r_phys[band]=self.r_phys_nocol#could add a colorfunc here
            if band !="VIS":
                self.ml[band]=self.draw_apparent_magnitude(self.Mr,self.zl,band)
            else: pass
            self.rl[band]=self.draw_apparent_size(self.r_phys[band],self.zl)
        return self.zl,self.sigl,self.ml,self.rl,self.ql

#====================================================================================

class SourcePopulation_(Population):
    def  __init__(self,
                  D=None,
                  bands= ['JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W'],
                  cosmo=[0.3,0.7,0.7],
                  population='jaguar'
                  ):

        self.bands = bands

        self.beginRedshiftDependentRelation(D)

        self.loadjaguar()

    def flux_to_mag(self, f):
        '''
        Converts flux in nJy to AB mag
        First checks if a value is zero and replaces it with the median
        So as not to get infs from the log
        '''
        median_flux = np.median(f[f > 0])
        f[f == 0] = median_flux
        flux_in_Jy = f*1e-9
        mag = -2.5*np.log10(flux_in_Jy) + 8.9
        return mag

    def loadjaguar(self):

        # source population simulated using a sky catalogue made for LSST
        # catalogue was generated by ray-tracing through the Millennium simulation
        # final catalogue is a 4.5x4.5 deg^2 footprint on the sky with halo masses down to 2.5 x 10^9 M_sun
        # resolution of the simulation means the catalogue is complete down to i ~ 27.5
        # this may not be deep enough to perfectly reconstruct the very faint lens population of LSST
        # is this lack of faint sources something we need to consider for COSMOS-Web?

        self.population = 'jaguar'

        self.data_type = 'holloway' # or sf_and_q

        print('loading {} data!'.format(self.data_type))

        #hdul_q = fits.open(r'/pbs/home/n/nhogg/git_lenspop/jaguar/JADES_Q_mock_r1_v1.2.fits') # CC-IN2P3
        #hdul_sf = fits.open(r'/pbs/home/n/nhogg/git_lenspop/jaguar/JADES_SF_mock_r1_v1.2.fits')

        # Load standard jaguar catalogues

        hdul_q = fits.open(r'/home/nataliehogg/Documents/Projects/cosmos_web/lenspop/jaguar/JADES_Q_mock_r1_v1.2.fits') # JADES catalogue made using JAGUAR sim; quiescent galaxies only
        hdul_sf = fits.open(r'/home/nataliehogg/Documents/Projects/cosmos_web/lenspop/jaguar/JADES_SF_mock_r1_v1.2.fits') # JADES catalogue made using JAGUAR sim; star-forming galaxies only

        data_q = hdul_q[1].data  # assume the first extension is a table
        data_sf = hdul_sf[1].data

        # get the number of sources in a single realisation
        single_realisation_source_number = len(list(data_q['redshift']) + list(data_sf['redshift']))

        # choose whether to use the standard catalogue or the Holloway modified one
        if self.data_type == 'sf_and_q':
            self.zc = np.array(list(data_q['redshift']) + list(data_sf['redshift'])) # kind of hacky way to join the two catalogues but whatever
            self.m = {}
            # these are fluxes in nJy
            m1 = np.array(list(data_q['NRC_F115W_fnu']) + list(data_sf['NRC_F115W_fnu']))
            m2 = np.array(list(data_q['NRC_F150W_fnu']) + list(data_sf['NRC_F150W_fnu']))
            m3 = np.array(list(data_q['NRC_F277W_fnu']) + list(data_sf['NRC_F277W_fnu']))
            m4 = np.array(list(data_q['NRC_F444W_fnu']) + list(data_sf['NRC_F444W_fnu']))
            # convert the fluxes to AB magnitudes
            self.m["JWST_NIRCam_F115W"] = self.flux_to_mag(m1)
            self.m["JWST_NIRCam_F150W"] = self.flux_to_mag(m2)
            self.m["JWST_NIRCam_F277W"] = self.flux_to_mag(m3)
            self.m["JWST_NIRCam_F444W"] = self.flux_to_mag(m4)
            self.mstar = np.array(list(data_q['mStar']) + list(data_sf['mStar'])) # log10 stellar mass
            self.r_eff = np.array(list(data_q['Re_circ']) + list(data_sf['Re_circ'])) # this is the effective *circularised* physical radius in kpc
            self.q = np.array(list(data_q['axis_ratio']) + list(data_sf['axis_ratio'])) # axis ratio
        elif self.data_type == 'holloway':
            # or load modified jaguar catalogue from Holloway et al, provided as a csv
            data = pd.read_csv(r'/home/nataliehogg/Documents/Projects/cosmos_web/lenspop/jaguar/holloway_data/Adapted_JAGUAR_Parent_Catalogue.csv')
            # Holloway catalogue is the full 10 realisations of 11x11 arcmin;
            # to match with the standard one used above (1 realisation of 11x11arcmin) we select that number from the Holloway catalogue
            # but the Holloway catalogue is sorted by redshift, so we cannot just take the first n lines as we will only get low z sources
            # so, we generate n indices at random and pick using those
            random_source_indices = np.random.randint(0, len(data['z'])-1, size=single_realisation_source_number)
            self.zc = np.array(data['z'][random_source_indices])
            self.m = {}
            # in the Holloway catalogue the fluxes have already been converted to magnitudes
            self.m["JWST_NIRCam_F115W"] = np.array(data['NRC_F115W_fnu'][random_source_indices])
            self.m["JWST_NIRCam_F150W"] = np.array(data['NRC_F150W_fnu'][random_source_indices])
            self.m["JWST_NIRCam_F277W"] = np.array(data['NRC_F277W_fnu'][random_source_indices])
            self.m["JWST_NIRCam_F444W"] = np.array(data['NRC_F444W_fnu'][random_source_indices])
            self.mstar = np.array(data['Log(M_star)'][random_source_indices])
            self.r_eff = np.array(data['R_eff (kpc)'][random_source_indices])
            # the Holloway catalogue did not preserve this info, reading it from the standard catalogue
            # thanks to the above random selection, the previous arrays will match in size to self.q
            self.q = np.array(list(data_q['axis_ratio']) + list(data_sf['axis_ratio']))
        else:
            print('I don\'t know that data type.')

    def RofMz(self, M, z, scatter=True, band=None):
        #band independent so far
        # print('I am in RofMz')
        # equation 5 of Tom's paper
        #{mosleh et al}, {Huang, Ferguson et al.}, Newton SLACS XI.
        # warning that this sometimes returns zeros
        r_phys=((M/-19.5)**-0.22)*((1.+z)/5.)**(-1.2) # equation 5 1507.02657
        # is the same as
        R=-(M+18.)/4.
        r_phys=(10**R)*((1.+z)/1.6)**(-1.2)

        if scatter!=False:
            if scatter==True:scatter=0.35 #dex
            self.scattered=10**(np.random.randn(len(r_phys))*scatter)
            r_phys*=self.scattered

        return r_phys

    # def draw_flattening_sourcepop(self, N):
    # for jaguar we don't need this; it's in the mock already (axis ratio)
    #     # print('I am in draw_flattening_sourcepop')
    #     # y=np.ones(N*1.5)*0.3
    #     y=np.ones(int(N*1.5))*0.3
    #     e=np.random.rayleigh(y)
    #     q=1-e
    #     q=q[q>0.2]
    #     q=q[:N]
    #
    #     return q

    def drawSourcePopulation(self, number, sourceplaneoverdensity=1, returnmasses=False):

        source_index = np.random.randint(0, len(self.zc), number*3)

        source_index=source_index[:number]

        self.zs=self.zc[source_index]

        self.r_phys=self.r_eff[source_index]

        self.ms={}
        for band in self.bands:
            self.ms[band]=self.m[band][source_index]

        self.rs=self.draw_apparent_size(self.r_phys, self.zs) # converts physical size to angular size in arcsec

        self.qs=self.q[source_index]

        self.ps=np.random.random_sample(number)*180

        # lsst sim has a source density of ~0.06 per square arcsecond
        # fac=(0.06)**-0.5
        # a=fac*(sourceplaneoverdensity)**-.5

        # there are 302,515 SF galaxies and 7464 quiescent galaxies in the 1st realisation of the JADES mock loaded above (len(self.zc))
        # the area of the mock is 11 arcmin^2
        # that is 660 arcseconds^2
        # the source density is therefore ~0.71 (309979/660^2)
        # if I try to compute the LSST source density in the same way I get ~0.03, and it's not clear where the factor 2 went
        # density = len(self.zc)/(660**2.)

        density = 0.71

        fac=(density)**-0.5
        a=fac*(sourceplaneoverdensity)**-0.5

        self.xs=(np.random.random_sample(number)-0.5)*a
        self.ys=(np.random.random_sample(number)-0.5)*a

        # the stuff above basically determines the SL cross-section
        # since the lens==True criterion is
        # b**2 > x**2 + y**2
        # so the bigger x and y are, the bigger the Einstein radii have to be to meet the criteria
        # since the Einstein radii are generated physically and in the same way each time
        # larger x and y means fewer lenses in the sample

        return self.zs, self.ms, self.xs, self.ys, self.qs, self.ps, self.rs

if __name__=="__main__":

    S2 = SourcePopulation_(population = "jaguar")
