import distances
from scipy import interpolate
# import cPickle
# import pickle
import pickle5 as pickle
import numpy
import math
import indexTricks as iT

# to correctly load pickled LSST file
# import sys
# from importlib import reload
# reload(sys)
# sys.setdefaultencoding('utf-8')

#====================================================================================

class RedshiftDependentRelation():
    def __init__(self,D=None,reset=False,cosmo=[0.3,0.7,0.7]):
        self.beginRedshiftDependentRelation(D,reset=reset,cosmo=cosmo)

    def beginRedshiftDependentRelation(self,D,reset,zmax=10,cosmo=[0.3,0.7,0.7]):
        self.zmax=zmax
        self.zbins,self.dz=numpy.linspace(0,self.zmax,401,retstep=True)
        self.z2bins,self.dz2=numpy.linspace(0,self.zmax,201,retstep=True)
        if D==None:
            import distances
            D=distances.Distance(cosmo=cosmo)
        self.D=D

        if reset!=True:
            try:
            #load useful redshift splines
                splinedump=open("redshiftsplines.pkl","rb")
                # self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline=cPickle.load(splinedump)
                self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline=pickle.load(splinedump)
            except IOError or EOFError:
                self.redshiftfunctions()
        else:
            self.redshiftfunctions()

    def redshiftfunctions(self):
        D=self.D
        zbins=self.zbins
        z2bins=self.z2bins
        Dabins=zbins*0.0
        Dmodbins=zbins*0.0
        Da2bins=numpy.zeros((z2bins.size,z2bins.size))
        volumebins=zbins*0.0
        for i in range(zbins.size):
            Dabins[i]=D.Da(zbins[i])
            Dmodbins[i]=D.distance_modulus(zbins[i])
            volumebins[i]=D.volume(zbins[i])
        for i in range(z2bins.size):
            for j in range(z2bins.size):
                if j>i:
                    Da2bins[i,j]=D.Da(z2bins[i],z2bins[j])

        self.Da_spline=interpolate.splrep(zbins,Dabins)
        self.Dmod_spline=interpolate.splrep(zbins,Dmodbins)

        self.volume_spline=interpolate.splrep(zbins,volumebins)

        z2d=iT.coords((z2bins.size,z2bins.size))*self.dz2
        self.Da_bispline=interpolate.RectBivariateSpline(z2bins,z2bins,Da2bins)

        #pickle the splines
        splinedump=open("redshiftsplines.pkl","wb")
        # cPickle.dump([self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline],splinedump,2)
        pickle.dump([self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline],splinedump,2)

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
            # print "don't know those units yet"
            print("don't know those units yet")
        # if z2==None:
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
            if len(z1)==1 and len(z2)!=1:z1=numpy.ones(len(z2))*z1[0]
            if len(z2)==1 and len(z1)!=1:z2=numpy.ones(len(z1))*z2[0]
            assert len(z1)==len(z2),"get it together"
            return z1,z2

#====================================================================================


class EinsteinRadiusTools(RedshiftDependentRelation):
    def  __init__(self,D=None,reset=False):
        self.beginRedshiftDependentRelation(D,reset)
        self.c=299792

    def sie_sig(self,rein,zl,zs):
        self.c=299792
        ds=self.Da(zs)
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
        # if colours==None:
        if colours is None:
            colours=0
            # print "warning no k-correction"
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
    def  __init__(self,zlmax=2,sigfloor=100,D=None,reset=True,
                  bands=['F814W_ACS','g_SDSS','r_SDSS','i_SDSS','z_SDSS','Y_UKIRT','VIS'],cosmo=[0.3,0.7,0.7]
                  ): #sadface
        self.sigfloor=sigfloor
        self.zlmax=zlmax
        self.bands=bands

        self.beginRedshiftDependentRelation(D,reset)
        self.beginLensPopulation(D,reset)


    def beginLensPopulation(self,D,reset):
        reset=True
        if reset!=True:
            try:
            #load Lens-population splines
                splinedump=open("lenspopsplines.pkl","rb")
                # self.cdfdNdzasspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,zlmax,sigfloor,self.colourspline,bands=cPickle.load(splinedump)
                self.cdfdNdzasspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,zlmax,sigfloor,self.colourspline,bands=pickle.load(splinedump)
            except IOError or EOFError or ValueError:
                self.lenspopfunctions()
            #check sigfloor and zlmax are same as requested
            if zlmax!=self.zlmax or self.sigfloor!=sigfloor:
                self.lenspopfunctions()
            #check all the necessary colours are included
            redocolours=False
            for band in self.bands:
                if band not in bands:redocolours=True
            if redocolours:
                self.Colourspline()
                self.lensPopSplineDump()
        else:
            self.lenspopfunctions()

    def lenspopfunctions(self):
        self.Psigzspline()
        self.Colourspline()
        self.lensPopSplineDump()

    def Psigzspline(self):
        #"""
        #drawing from a 2d pdf is a pain; should probably make this into its own module
        self.zlbins,self.dzl=numpy.linspace(0,self.zlmax,201,retstep=True)
        sigmas=numpy.linspace(self.sigfloor,400,401)
        self.sigbins=sigmas
        dNdz=self.zlbins*0
        Csiggivenz=numpy.zeros((sigmas.size,self.zlbins.size))
        CDFbins=numpy.linspace(0,1,1001)
        siggivenCz=numpy.zeros((CDFbins.size,self.zlbins.size))
        for i in range(len(self.zlbins)):
            z=self.zlbins[i]
            dphidsiggivenz=self.phi(sigmas,z)
            phisigspline=interpolate.splrep(sigmas,dphidsiggivenz)
            tot=interpolate.splint(self.sigfloor,500,phisigspline)
            Csiggivenz[:,i]=numpy.cumsum(dphidsiggivenz)/numpy.sum(dphidsiggivenz)
            Csiggivenzspline=interpolate.splrep(Csiggivenz[:,i],sigmas)
            siggivenCz[:,i]=interpolate.splev(CDFbins,Csiggivenzspline)
            if z!=0:
                dNdz[i]=tot*(self.Volume(z)-self.Volume(z-self.dzl))/self.dzl

        Nofzcdf=numpy.cumsum(dNdz)/numpy.sum(dNdz)
        #import pylab as plt
        #plt.plot(self.zlbins,Nofzcdf)
        #plt.show()
        #exit()
        self.cdfdNdzasspline=interpolate.splrep(Nofzcdf,self.zlbins)

        self.dNdzspline=interpolate.splrep(self.zlbins,dNdz)
        N=interpolate.splint(0,self.zlmax,self.dNdzspline)

        self.cdfdsigdzasspline=interpolate.RectBivariateSpline(\
            CDFbins,self.zlbins,siggivenCz)

        dphidsiggivenz0=self.phi(sigmas,sigmas*0)
        cdfdNdsigz0=dphidsiggivenz0.cumsum()/dphidsiggivenz0.sum()
        self.cdfdNdsigz0asspline=interpolate.splrep(cdfdNdsigz0,sigmas)


        #"""
        #phi is redshift independant.



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
        # cPickle.dump([self.cdfdNdzasspline,self.cdfdNdsigz0asspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,self.zlmax,self.sigfloor,self.colourspline,self.bands],splinedump,2)
        pickle.dump([self.cdfdNdzasspline,self.cdfdNdsigz0asspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,self.zlmax,self.sigfloor,self.colourspline,self.bands],splinedump,2)

    def draw_z(self,N):
        return interpolate.splev(numpy.random.random(N),self.cdfdNdzasspline)

    def draw_sigma(self,z):
        try: len(z)
        except TypeError:z=[z]
        if self.nozdependence:
            sigs =interpolate.splev(numpy.random.random(len(z)),self.cdfdNdsigz0asspline)
            return sigs
        else:
            # print "Warning: drawing from 2dpdf is low accuracy"
            print("Warning: drawing from 2dpdf is low accuracy")
            return self.cdfdsigdzasspline.ev(numpy.random.random(len(z)),z)

    def draw_zsig(self,N):
        z=self.draw_z(N)
        sig=self.draw_sigma(z)
        return z,sig

    def EarlyTypeRelations(self,sigma,z=None,scatter=True,band=None):#z dependence not encoded currently
        #Hyde and Bernardi, M = r band absolute magnitude.
        V=numpy.log10(sigma)
        Mr=(-0.37+(0.37**2-(4*(0.006)*(2.97+V)))**0.5)/(2*0.006)
        if scatter:
            Mr+=numpy.random.randn(len(Mr))*(0.15/2.4)

        #R=4.72+0.63*Mr+0.02*Mr**2 #rest-frame R_band size.
        R=2.46-2.79*V+0.84*V**2
        if scatter:
            R+=numpy.random.randn(len(R))*0.11

        #convert to observed r band size;
        r_phys = 10**R

        return Mr,r_phys

    def colour(self,z,band):
        return interpolate.splev(z,self.colourspline[band])

    def Ndeflectors(self,z,zmin=0,fsky=1):
        if zmin>z:
            z,zmin=zmin,z
        N=interpolate.splint(zmin,z,self.dNdzspline)
        N*=fsky
        return N

    def phi(self,sigma,z):
        sigma[sigma==0]+=1e-6
        phi_star=(8*10**-3)*self.D.h**3
        alpha=2.32
        beta=2.67
        sigst=161
        phi=phi_star * \
            ((sigma*1./sigst)**alpha)*\
            numpy.exp(-(sigma*1./sigst)**beta)*beta/\
            math.gamma(alpha*1./beta)/\
            (1.*sigma)

        phi*=(1+z)**(-2.5)
        return phi

    def draw_flattening_lenspop(self,sigma,z=None):
        x=sigma
        y=0.378-0.000572*x
        e=numpy.random.rayleigh(y)
        q=1-e
        #dont like ultraflattened masses:
        while len(q[q<0.2])>0 or len(q[q>1])>0:
            q[q<0.2]=1-numpy.random.rayleigh(y[q<0.2])
            q[q>1]=1-numpy.random.rayleigh(y[q>1])
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
                  reset=False,
                  bands=['F814W_ACS','g_SDSS','r_SDSS','i_SDSS','z_SDSS','Y_UKIRT','VIS'],
                  cosmo=[0.3,0.7,0.7],
                  population='lsst'
                  ):

        self.bands = bands

        self.beginRedshiftDependentRelation(D, reset)

        self.loadlsst()

    def loadlsst(self):

        self.population = 'lsst'

        print("new lsst catalogue")

        data = numpy.genfromtxt('lsst.1sqdegree_catalog2.txt', delimiter=',') # NH: added the delimiter type so it reads the file correctly

        self.zc = data[:,2]

        self.m = {}

        #print data[:,0].max()-data[:,0].min()
        #print data[:,1].max()-data[:,1].min()

        self.m["g_SDSS"]=data[:,3]
        self.m["r_SDSS"]=data[:,4]
        self.m["i_SDSS"]=data[:,5]
        self.m["z_SDSS"]=data[:,6]
        self.m["F814W_ACS"]=data[:,5] # we'll make do with F814==i
        self.m["Y_UKIRT"]=data[:,6]*99 #there is no Y band data atm
        self.mstar=data[:,12]
        self.mhalo=data[:,13]
        self.m["VIS"]=(self.m["r_SDSS"]+self.m["i_SDSS"]+self.m["z_SDSS"])/3
        self.Mv=data[:,7]

    def RofMz(self,M,z,scatter=True,band=None):#band independent so far
    #{mosleh et al}, {Huang, Ferguson et al.}, Newton SLACS XI.
        r_phys=((M/-19.5)**-0.22)*((1.+z)/5.)**(-1.2) # equation 5 1507.02657
        # is the same as
        R=-(M+18.)/4.
        r_phys=(10**R)*((1.+z)/1.6)**(-1.2)

        if scatter!=False:
            if scatter==True:scatter=0.35 #dex
            self.scattered=10**(numpy.random.randn(len(r_phys))*scatter)
            r_phys*=self.scattered

        return r_phys


    def draw_flattening_sourcepop(self,N):
        # y=numpy.ones(N*1.5)*0.3
        y=numpy.ones(int(N*1.5))*0.3
        e=numpy.random.rayleigh(y)
        q=1-e
        q=q[q>0.2]
        q=q[:N]

        return q

    def drawSourcePopulation(self,number,sourceplaneoverdensity=10,returnmasses=False):
        source_index=numpy.random.randint(0,len(self.zc),number*3)
        #source_index=source_index[((self.zc[source_index]<10) & (self.zc[source_index]>0.05))]
        source_index=source_index[:number]
        self.zs=self.zc[source_index]
        self.Mvs=self.Mv[source_index]
        self.ms={}
        for band in self.bands:
            if band !="VIS":
                self.ms[band]=self.m[band][source_index]
            else:
                self.ms[band]=(self.m["r_SDSS"][source_index]+self.m["i_SDSS"][source_index]+self.m["z_SDSS"][source_index])/3.

        self.r_phys=self.RofMz(self.Mvs,self.zs,scatter=True)
        self.rs=self.draw_apparent_size(self.r_phys,self.zs)
        self.qs=self.draw_flattening_sourcepop(number)

        self.ps=numpy.random.random_sample(number )*180

        # #cosmos has a source density of ~0.015 per square arcsecond
        # if self.population=="cosmos":
        #     fac=(0.015)**-0.5
        #     a=fac*(sourceplaneoverdensity)**-.5
        # #lsst sim has a source density of ~0.06 per square arcsecond
        # elif self.population=="lsst":
        #     fac=(0.06)**-0.5
        #     a=fac*(sourceplaneoverdensity)**-.5

        #lsst sim has a source density of ~0.06 per square arcsecond
        fac=(0.06)**-0.5
        a=fac*(sourceplaneoverdensity)**-.5

        # else:
        #     pass

        self.xs=(numpy.random.random_sample(number)-0.5)*a
        self.ys=(numpy.random.random_sample(number)-0.5)*a

        if returnmasses:
            self.mstar_src=self.mstar[source_index]
            self.mhalo_src=self.mhalo[source_index]
            return self.zs,self.ms,self.xs,self.ys,self.qs,self.ps,self.rs,self.mstar_src,self.mhalo_src

        return self.zs,self.ms,self.xs,self.ys,self.qs,self.ps,self.rs


class AnalyticSourcePopulation_(SourcePopulation_):
    def  __init__(self,D=None,reset=False,
                  bands=['F814W_ACS','g_SDSS','r_SDSS','i_SDSS','z_SDSS','Y_UKIRT'],cosmo=[0.3,0.7,0.7]
                  ):
        self.bands=bands
        self.beginRedshiftDependentRelation(D,reset)
        # print "not written!"
        print("not written!")



if __name__=="__main__":
    #RedshiftDependentRelation(reset=True)

    #L=LensPopulation_(reset=True,sigfloor=100)

    # S=SourcePopulation_(reset=False,population="cosmos")
    S2=SourcePopulation_(reset=False,population="lsst")


    # print numpy.median(S.Mv[S.m["i_SDSS"]<25])-numpy.median(S2.Mv[S2.m["i_SDSS"]<25])
    # print len(S.Mv[S.m["i_SDSS"]<25])*1./(len(S2.Mv[S2.m["i_SDSS"]<25])*100)
    # print len(S.Mv)/(60.**2)/2.
    # print len(S2.Mv[S2.m["i_SDSS"]<25])/(0.2**2)/(60.**2)
    # print len(S2.Mv)/(0.2**2)/(60.**2)

    # print(numpy.median(S.Mv[S.m["i_SDSS"]<25])-numpy.median(S2.Mv[S2.m["i_SDSS"]<25]))
    # print(len(S.Mv[S.m["i_SDSS"]<25])*1./(len(S2.Mv[S2.m["i_SDSS"]<25])*100))
    # print(len(S.Mv)/(60.**2)/2.)
    print(len(S2.Mv[S2.m["i_SDSS"]<25])/(0.2**2)/(60.**2))
    print(len(S2.Mv)/(0.2**2)/(60.**2))

    #print EarlyTypeRelations(self,100,z=None,scatter=True,band=None)
