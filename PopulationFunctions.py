import distances
from scipy import interpolate
import pickle
import numpy as np
import math
import indexTricks as iT
from astropy.io import fits
import pandas as pd
import os
from pathlib import Path


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

def find_jaguar_catalog(filename='JADES_Q_mock_r1_v1.2.fits'):
    """
    Auto-detect the path to JAGUAR catalog files.
    Searches in multiple common locations.

    Parameters
    ----------
    filename : str
        Name of the catalog file to find

    Returns
    -------
    str
        Full path to the catalog directory

    Raises
    ------
    FileNotFoundError
        If catalog cannot be found in any common location
    """
    # Get script directory
    script_dir = Path(__file__).parent.absolute()

    # List of paths to search (in order of priority)
    search_paths = [
        # Environment variable (highest priority)
        os.environ.get('JAGUAR_DATA_PATH'),

        # Relative to current working directory
        Path.cwd() / 'jaguar',
        Path.cwd() / 'data' / 'jaguar',

        # Relative to script directory
        script_dir / 'jaguar',
        script_dir / 'data' / 'jaguar',

        # Parent directories
        script_dir.parent / 'jaguar',
        script_dir.parent / 'lenspop' / 'jaguar',
        script_dir.parent / 'cosmos_web' / 'lenspop' / 'jaguar',

        # Common absolute paths
        Path.home() / 'Documents' / 'Projects' / 'cosmos_web' / 'lenspop' / 'jaguar',
        Path('/home/nataliehogg/Documents/Projects/cosmos_web/lenspop/jaguar'),
        Path('/pbs/home/n/nhogg/git_lenspop/jaguar'),  # CC-IN2P3
    ]

    # Filter out None values and search
    for path in search_paths:
        if path is None:
            continue

        path = Path(path)
        catalog_file = path / filename

        if catalog_file.exists():
            print(f'Found JAGUAR catalog at: {path}')
            return str(path)

    # If not found, provide helpful error message
    error_msg = f"""
    Could not find JAGUAR catalog file '{filename}' in any of the following locations:
    """
    for path in search_paths:
        if path is not None:
            error_msg += f"\n  - {path}"

    error_msg += f"""

    Please either:
    1. Place the JAGUAR catalog files in one of the above locations, or
    2. Set the JAGUAR_DATA_PATH environment variable:
       export JAGUAR_DATA_PATH=/path/to/jaguar/directory

    The catalog files needed are:
    - JADES_Q_mock_r1_v1.2.fits
    - JADES_SF_mock_r1_v1.2.fits
    """

    raise FileNotFoundError(error_msg)

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

        # Load the appropriate source catalogue
        if population == 'lsst':
            self.loadlsst()
        elif population == 'cosmos':
            self.loadcosmos()
        elif population == 'jaguar':
            self.loadjaguar()
        elif population == 'flagship':
            self.loadflagship()
        else:
            raise ValueError(f"Unknown source population: {population}. Use 'lsst', 'jaguar', 'flagship', or 'cosmos'.")

    def flux_to_mag(self, f, waveband):
        '''
        Converts flux in nJy to AB mag
        First checks if a value is zero and replaces it with the median
        So as not to get infs from the log
        '''
        median_flux = np.median(f[f > 0])
        f[f == 0] = median_flux
        pixel_sr_dict = {"F115W": 2.12313551253570e-14,
                         "F150W": 2.12313551253570e-14,
                         "F277W": 8.46363636363636e-14,
                         "F444W": 8.46363636363636e-14}
        pixar_sr = pixel_sr_dict[waveband]
        zero_point = -6.10 - 2.5*np.log10(pixar_sr)
        magnitude_ab_pix = zero_point - 2.5*np.log10(f)
        return magnitude_ab_pix

    def loadlsst(self):
        """
        Load LSST 1 sq degree catalog
        Format: ra, dec, redshift, g_ab, r_ab, i_ab, z_ab, absmag_r_total,
                bulge_n, disk_n, BulgeHalfLightRadius, gal_type, mass_stellar, mass_halo
        """
        self.population = 'lsst'

        # Load pickled LSST catalog
        print('loading LSST catalogue!')
        f = open('lsst.1sqdegree_catalog2.pkl', 'rb')
        data = pickle.load(f, encoding='latin1')
        f.close()

        # Extract columns
        self.zc = data[:, 2]  # redshift

        # Magnitudes (g, r, i, z bands)
        self.m = {}
        self.m["g_SDSS"] = data[:, 3]
        self.m["r_SDSS"] = data[:, 4]
        self.m["i_SDSS"] = data[:, 5]
        self.m["z_SDSS"] = data[:, 6]

        # For JWST bands - use SDSS as proxy (will need K-correction)
        # This is approximate - ideally would compute proper K-corrections
        if 'JWST_NIRCam_F115W' in self.bands:
            self.m["JWST_NIRCam_F115W"] = data[:, 5]  # ~i-band
            self.m["JWST_NIRCam_F150W"] = data[:, 6]  # ~z-band
            self.m["JWST_NIRCam_F277W"] = data[:, 6]  # ~z-band
            self.m["JWST_NIRCam_F444W"] = data[:, 6]  # ~z-band

        # VIS band for Euclid (average of r, i, z)
        self.m["VIS"] = (data[:, 4] + data[:, 5] + data[:, 6]) / 3.0

        # Effective radius in kpc (column 10)
        self.r_eff = data[:, 10]

        # Masses
        self.mstar = data[:, 12]  # stellar mass

        # LSST catalog doesn't have axis ratio and position angle
        # Generate them randomly following typical distributions
        N = len(self.zc)

        # Axis ratio: Rayleigh distribution with sigma=0.3
        y = np.ones(N) * 0.3
        e = np.random.rayleigh(y)
        self.q = 1 - e
        self.q[self.q < 0.2] = 0.2  # Minimum axis ratio
        self.q[self.q > 1.0] = 0.9  # Maximum axis ratio

        # Position angle: uniform 0-180 degrees
        self.p = np.random.random_sample(N) * 180.0

    def loadflagship(self):
        """
        Load Euclid Flagship mock catalog
        Paper: https://arxiv.org/abs/2405.13495

        Expected columns:
        - redshift_observed (or redshift): with peculiar velocity
        - VIS_flux_obs_continuum: observed flux in erg/cm²/s/Hz
        - r_bulge_arcsec: bulge half-light radius in arcsec
        - r_disk_arcsec: disk half-light radius in arcsec
        - axis_ratio: b/a ratio (medium/major axis)
        """
        self.population = 'flagship'

        print('Loading Euclid Flagship catalogue!')

        # Load catalog - try multiple locations and formats
        catalog_loaded = False
        data = None

        # Search paths in order of priority
        search_paths = [
            'euclid_flagship/',  # Subdirectory
            './',                 # Current directory
            'flagship/',          # Alternative subdirectory
        ]

        filenames = [
            'euclid_flagship_catalogue.fits',  # British spelling
            'euclid_flagship_catalog.fits',
            'flagship_catalog.fits',
            'flagship_catalogue.fits',
            'euclid_flagship.fits',
            'flagship.fits',
        ]

        # Try FITS format in all locations
        for path in search_paths:
            for fname in filenames:
                filepath = os.path.join(path, fname)
                try:
                    hdul = fits.open(filepath)
                    data = hdul[1].data
                    hdul.close()
                    catalog_loaded = True
                    print(f'  Loaded from FITS: {filepath}')
                    break
                except (FileNotFoundError, OSError):
                    continue
            if catalog_loaded:
                break

        # Try pickle format if FITS not found
        if not catalog_loaded:
            for path in search_paths:
                for fname in ['flagship_catalog.pkl', 'flagship.pkl']:
                    filepath = os.path.join(path, fname)
                    try:
                        with open(filepath, 'rb') as f:
                            data = pickle.load(f)
                        catalog_loaded = True
                        print(f'  Loaded from pickle: {filepath}')
                        break
                    except (FileNotFoundError, OSError):
                        continue
                if catalog_loaded:
                    break

        # Try numpy format if others not found
        if not catalog_loaded:
            for path in search_paths:
                for fname in ['flagship_catalog.npy', 'flagship.npy']:
                    filepath = os.path.join(path, fname)
                    try:
                        data = np.load(filepath, allow_pickle=True)
                        catalog_loaded = True
                        print(f'  Loaded from numpy: {filepath}')
                        break
                    except (FileNotFoundError, OSError):
                        continue
                if catalog_loaded:
                    break

        if not catalog_loaded:
            raise FileNotFoundError("""
            Could not find Flagship catalog. Searched for:
              - euclid_flagship/flagship_catalog.fits
              - euclid_flagship/euclid_flagship.fits
              - ./flagship_catalog.fits
              - And other locations/formats

            Please ensure the catalog file is in the euclid_flagship/ directory.
            """)

        # Extract redshifts (try different column names)
        if 'observed_redshift_gal' in data.dtype.names:
            self.zc = np.array(data['observed_redshift_gal'])
            print('  Using observed redshift (with peculiar velocity)')
        elif 'redshift_observed' in data.dtype.names:
            self.zc = np.array(data['redshift_observed'])
            print('  Using observed redshift (with peculiar velocity)')
        elif 'redshift_obs' in data.dtype.names:
            self.zc = np.array(data['redshift_obs'])
            print('  Using observed redshift (with peculiar velocity)')
        elif 'redshift' in data.dtype.names:
            self.zc = np.array(data['redshift'])
            print('  Using redshift (check if this includes peculiar velocity!)')
        else:
            raise KeyError("Could not find redshift column. Expected 'observed_redshift_gal', 'redshift_observed', 'redshift_obs', or 'redshift'")

        # --- VIS MAGNITUDE: Convert from flux ---
        # Try different possible column names
        flux_column_names = ['euclid_vis', 'VIS_flux_obs_continuum', 'flux_VIS_obs', 'VIS_flux', 'flux_vis_continuum']
        flux_vis_erg = None

        for col_name in flux_column_names:
            if col_name in data.dtype.names:
                flux_vis_erg = np.array(data[col_name])
                print(f'  Using VIS flux from column: {col_name}')
                break

        if flux_vis_erg is None:
            raise KeyError(f"Could not find VIS flux column. Tried: {flux_column_names}")

        # Handle zero/negative fluxes
        flux_safe = np.copy(flux_vis_erg)
        n_bad = np.sum(flux_safe <= 0)
        if n_bad > 0:
            median_flux = np.median(flux_vis_erg[flux_vis_erg > 0])
            flux_safe[flux_safe <= 0] = median_flux / 100.0  # Very faint
            print(f'  Warning: {n_bad} sources with flux ≤ 0, set to faint value')

        # Convert flux (erg/cm²/s/Hz) to AB magnitude
        # AB_mag = -2.5 * log10(f_ν [erg/cm²/s/Hz]) - 48.60
        self.m = {}
        self.m["VIS"] = -2.5 * np.log10(flux_safe) - 48.60

        # Also populate SDSS bands for compatibility (use VIS as proxy)
        self.m["r_SDSS"] = self.m["VIS"]
        self.m["i_SDSS"] = self.m["VIS"]
        self.m["z_SDSS"] = self.m["VIS"]
        self.m["g_SDSS"] = self.m["VIS"]

        # --- EFFECTIVE RADIUS: Combine bulge + disk ---
        # Try different possible column names
        r_bulge_names = ['bulge_r50', 'r_bulge_arcsec', 'bulge_half_light_radius', 'r_bulge']
        r_disk_names = ['disk_r50', 'r_disk_arcsec', 'disk_half_light_radius', 'r_disk']

        r_bulge_arcsec = None
        r_disk_arcsec = None

        for col_name in r_bulge_names:
            if col_name in data.dtype.names:
                r_bulge_arcsec = np.array(data[col_name])
                print(f'  Found bulge radius: {col_name}')
                break

        for col_name in r_disk_names:
            if col_name in data.dtype.names:
                r_disk_arcsec = np.array(data[col_name])
                print(f'  Found disk radius: {col_name}')
                break

        if r_bulge_arcsec is None or r_disk_arcsec is None:
            raise KeyError(f"Could not find both bulge and disk radius columns.\n  Tried bulge: {r_bulge_names}\n  Tried disk: {r_disk_names}")

        # Use maximum of bulge and disk (conservative - captures full extent)
        r_eff_arcsec = np.maximum(r_bulge_arcsec, r_disk_arcsec)

        # Convert from angular (arcsec) to physical (kpc)
        # r_kpc = r_arcsec * Da(z) / 206264
        # Da(z) returns angular diameter distance in kpc
        r_eff_kpc = r_eff_arcsec * self.Da(self.zc, units="kpc") / 206264.0

        self.r_eff = r_eff_kpc

        # Sanity check for sizes
        n_zero = np.sum(self.r_eff <= 0)
        if n_zero > 0:
            print(f'  Warning: {n_zero} sources with r_eff ≤ 0')
            # Set minimum size
            self.r_eff[self.r_eff <= 0] = 0.1  # 0.1 kpc minimum

        # --- AXIS RATIO ---
        axis_ratio_names = ['q_gal', 'axis_ratio', 'b_to_a', 'q', 'axis_ratio_3d']
        self.q = None

        for col_name in axis_ratio_names:
            if col_name in data.dtype.names:
                self.q = np.array(data[col_name])
                print(f'  Found axis ratio: {col_name}')
                break

        if self.q is None:
            # Generate if not in catalog
            print('  Warning: axis_ratio not found, generating randomly')
            N = len(self.zc)
            y = np.ones(N) * 0.3
            e = np.random.rayleigh(y)
            self.q = 1 - e

        # Sanity check: ensure 0.2 < q ≤ 1.0
        self.q[self.q < 0.2] = 0.2
        self.q[self.q > 1.0] = 1.0

        # --- POSITION ANGLE ---
        # Flagship doesn't have this - generate randomly
        N = len(self.zc)
        self.p = np.random.random_sample(N) * 180.0  # Uniform 0-180°
        print('  Position angle: generated randomly (not in catalog)')

        # --- STELLAR MASS (optional, for record-keeping) ---
        mass_names = ['stellar_mass', 'mstar', 'mass_stellar', 'log_mstar']
        self.mstar = None

        for col_name in mass_names:
            if col_name in data.dtype.names:
                self.mstar = np.array(data[col_name])
                print(f'  Found stellar mass: {col_name}')
                break

        if self.mstar is None:
            self.mstar = np.zeros(N)  # Placeholder if not available

        # Print summary
        print(f'\n  Flagship catalog summary:')
        print(f'    Sources loaded: {N}')
        print(f'    Redshift range: {self.zc.min():.3f} - {self.zc.max():.3f}')
        print(f'    VIS mag range: {self.m["VIS"].min():.2f} - {self.m["VIS"].max():.2f}')
        print(f'    r_eff range: {self.r_eff.min():.4f} - {self.r_eff.max():.2f} kpc')
        print(f'    Axis ratio range: {self.q.min():.2f} - {self.q.max():.2f}')

    def loadcosmos(self):
        """
        Load COSMOS catalog
        Note: This method is not yet implemented
        """
        raise NotImplementedError("COSMOS catalog loading not yet implemented. Use 'lsst', 'jaguar', or 'flagship' instead.")

    def loadjaguar(self):

        # source population simulated using a sky catalogue made for LSST
        # catalogue was generated by ray-tracing through the Millennium simulation
        # final catalogue is a 4.5x4.5 deg^2 footprint on the sky with halo masses down to 2.5 x 10^9 M_sun
        # resolution of the simulation means the catalogue is complete down to i ~ 27.5
        # this may not be deep enough to perfectly reconstruct the very faint lens population of LSST
        # is this lack of faint sources something we need to consider for COSMOS-Web?

        self.population = 'jaguar'

        self.data_type = 'sf_and_q' # sf_and_q or holloway

        print('loading {} data!'.format(self.data_type))

        # Auto-detect the path to JAGUAR catalog files
        jaguar_path = find_jaguar_catalog('JADES_Q_mock_r1_v1.2.fits')

        # Load standard jaguar catalogues
        hdul_q = fits.open(os.path.join(jaguar_path, 'JADES_Q_mock_r1_v1.2.fits'))  # Quiescent galaxies
        hdul_sf = fits.open(os.path.join(jaguar_path, 'JADES_SF_mock_r1_v1.2.fits'))  # Star-forming galaxies

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
            self.m["JWST_NIRCam_F115W"] = self.flux_to_mag(m1, 'F115W')
            self.m["JWST_NIRCam_F150W"] = self.flux_to_mag(m2, 'F150W')
            self.m["JWST_NIRCam_F277W"] = self.flux_to_mag(m3, 'F277W')
            self.m["JWST_NIRCam_F444W"] = self.flux_to_mag(m4, 'F444W')
            self.mstar = np.array(list(data_q['mStar']) + list(data_sf['mStar'])) # log10 stellar mass
            self.r_eff = np.array(list(data_q['Re_circ']) + list(data_sf['Re_circ'])) # this is the effective *circularised* physical radius in kpc
            self.q = np.array(list(data_q['axis_ratio']) + list(data_sf['axis_ratio'])) # axis ratio
            self.p = np.array(list(data_q['position_angle']) + list(data_sf['position_angle']))
        elif self.data_type == 'holloway':
            # or load modified jaguar catalogue from Holloway et al, provided as a csv
            holloway_path = os.path.join(jaguar_path, 'holloway_data', 'Adapted_JAGUAR_Parent_Catalogue.csv')
            data = pd.read_csv(holloway_path)
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
            self.p = np.array(list(data_q['position_angle']) + list(data_sf['position_angle']))
        else:
            print('I don\'t know that data type.')

    def RofMz(self, M, z, scatter=True, band=None):
        '''
        effective physical radius of source light profile
        '''
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

        self.ps=self.p[source_index] #np.random.random_sample(number)*180; indeed the position angle can be taken from the catalogues too

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
