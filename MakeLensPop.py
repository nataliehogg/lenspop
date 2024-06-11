import distances
from scipy import interpolate
# import cPickle
import pickle
import numpy
import math
import time
import indexTricks as iT
# import pylab as plt
from PopulationFunctions import *

class LensPopulation(LensPopulation_):
    def  __init__(self, zlmax=3, sigfloor=100, D=None, reset=True,
                  bands=['JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W']
                  ): #sadface
        self.sigfloor=sigfloor
        self.zlmax=zlmax
        self.bands=bands

        self.beginRedshiftDependentRelation(D)
        self.beginLensPopulation(D)

    def phi(self, sigma, z):
    #you can change this, but remember to reset the splines if you do.
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

        # phi*=(1+z)**(-2.5)
        self.nozdependence=True

        return phi


class SourcePopulation(SourcePopulation_):
    def  __init__(self, D=None, reset=False,
                  bands=['JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W'],
                  population="jaguar"
                  ):
        self.bands=bands
        self.beginRedshiftDependentRelation(D, reset)
        if population=="cosmos":
            self.loadcosmos()
        elif population=="lsst":
            self.loadlsst()
        elif population=="jaguar":
            self.loadjaguar()

class LensSample():
    """
    Wrapper for all the other objects so you can just call it, and then run
    Generate_Lens_Pop to get a fairly drawn lens population
    """
    def  __init__(self, D=None,
                  zlmax=3, sigfloor=100,
                  bands=['JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W'],
                  cosmo=[0.3,0.7,0.7], sourcepop="jaguar"
                  ):
        self.sourcepopulation=sourcepop
        if D==None:
            import distances
            D=distances.Distance(cosmo=cosmo)

        reset = False

        self.L=LensPopulation(reset=reset,sigfloor=sigfloor,zlmax=zlmax,bands=bands,D=D)

        self.S=SourcePopulation(reset=reset,bands=bands,D=D,population=sourcepop)

        self.E=EinsteinRadiusTools(D=D)

    def Lenses_on_sky(self):
        self.ndeflectors=self.L.Ndeflectors(self.L.zlmax)
        return self.ndeflectors

    def Generate_Lens_Pop(self, N, firstod=1, nsources=1, prunenonlenses=True, save=True):

        t0 = time.perf_counter()

        if prunenonlenses==False:
            assert N<60000

        self.lens={}
        self.reallens={}
        M=N*1 # initialise M, the number of deflectors computed in PopulationFunctions -> Ndeflectors -> dn/dz, sigma, phi, all that business
        l=-1 # initialise number of lenses (?)
        l2=-1 # initialise number of lenses meeting b^2 > x^ + y^2 criterion

        while M>0: # while the number of deflectors is non-zero... (i.e. this looks at every deflector to see if it's lensing a source?)
            timeleft= 1 # fake number to initialise
            if M!=N:
                tnow = time.perf_counter()
                ti = (tnow-t0)/float(N-M)
                timeleft = ti*M/60.

            print('{:.1f} minutes ({:.1f} hours) left'.format(timeleft, timeleft/60))

            if M>100000: # if the number of deflectors is over 100,000
                n=100000 # start drawing the lens and source pops at 100,000
            else: # otherwise if M < 100,000
                n=M*1 # do the whole lot in one go
            M-=n # then subtract the number computed from the number to be computed

            zl, sigl, ml, rl, ql = self.L.drawLensPopulation(n) # draw n lenses (from the total possible number of deflectors, which is M)

            zs, ms, xs, ys, qs, ps, rs = self.S.drawSourcePopulation(n*nsources, sourceplaneoverdensity=firstod, returnmasses=False) # draw n sources (from the total possible number of sources, which is M)

            zl1=zl*1
            sigl1=sigl*1
            for i in range(nsources-1):
                zl=numpy.concatenate((zl,zl1))
                sigl=numpy.concatenate((sigl,sigl1))

            b=self.E.sie_rein(sigl,zl,zs)
            for i in range(n): # for each of the n lenses and sources, check to see if the nth lens lenses the nth source
                l +=1
                self.lens[l]={}
                if b[i]**2>(xs[i]**2+ys[i]**2):
                    self.lens[l]["lens?"]=True
                else:
                    self.lens[l]["lens?"]=False

                self.lens[l]["b"]={} # Einstein radius
                self.lens[l]["zs"]={} # source redshift
                self.lens[l]["zl"]=zl[i] # lens redshift
                self.lens[l]["sigl"]=sigl[i] # velocity dispersion
                for j in range(nsources):
                    self.lens[l]["zs"][j+1]=zs[i+j*n]
                    self.lens[l]["b"][j+1] =b[i+j*n]

                self.lens[l]["ml"]={}
                self.lens[l]["rl"]={}
                self.lens[l]["ms"]={}

                for band in ml.keys():
                        self.lens[l]["ml"][band]=ml[band][i]
                        self.lens[l]["rl"][band]=rl[band][i]
                self.lens[l]["ql"]=ql[i]

                self.lens[l]["ms"]={}
                self.lens[l]["xs"]={}
                self.lens[l]["ys"]={}
                self.lens[l]["rs"]={}
                self.lens[l]["qs"]={}
                self.lens[l]["ps"]={}
                # self.lens[l]["mstar"]={}
                # self.lens[l]["mhalo"]={}

                for j in range(nsources):
                    self.lens[l]["ms"][j+1]={}
                    for band in ml.keys():
                        self.lens[l]["ms"][j+1][band]=ms[band][i+j*n]
                    self.lens[l]["zs"][j+1]=zs[i+j*n]
                    self.lens[l]["b"][j+1] =b[i+j*n]
                    self.lens[l]["xs"][j+1]=xs[i+j*n]
                    self.lens[l]["ys"][j+1]=ys[i+j*n]
                    self.lens[l]["rs"][j+1]=rs[i+j*n]
                    self.lens[l]["qs"][j+1]=qs[i+j*n]
                    self.lens[l]["ps"][j+1]=ps[i+j*n]
                    # self.lens[l]["mhalo"][j+1]=mstar[i+j*n]
                    # self.lens[l]["mstar"][j+1]=mhalo[i+j*n]


                if self.lens[l]["lens?"]: # NH: what is this doing? for each l it looks to see if lens? is True or False
                    if prunenonlenses == True: # NH: not keen on this happening on the fly... why not save all then remove in a separate step?
                        l2+=1 # adds one to the number of true lenses

                        self.reallens[l2] = self.lens[l].copy()

                        del self.lens # NH: isn't this just killing everything? maybe if we just comment this out and run, everything will be saved ok
                        self.lens={} # NH: no because what should be saves is self.reallens. anyway still not sure why this step is here

                        # if l2%1000==0:
                        #     # print l2
                        #     print(l2)

                        if (l2+1)%10000==0: # if (l2+1)/10000 has no remainder:
                            filename = 'idealisedlenses/lenspopulation_{}_{}.pkl'.format(self.sourcepopulation, l2-10000+1) # NH: also got no clue what all these hardcoded numbers are
                            out_file = open(filename, 'wb')
                            # pickle.dump(self.reallens, out_file, 2) # NH: so the 2 here is the protocol; we ought to increase this, though two is backwards compatible for py3
                            # 5 is fastest, which may be relevant for the amount of data we have
                            # https://stackoverflow.com/questions/23582489/python-pickle-protocol-choice
                            pickle.dump(self.reallens, out_file, protocol=5)
                            out_file.close()

                            del self.reallens
                            self.reallens={}

        if save:
            fn = "idealisedlenses/lenspopulation_%s_residual_%i.pkl"%(self.sourcepopulation, l2)
            f = open(fn, 'wb')
            pickle.dump(self.reallens, f, protocol=5) # NH: update protocol
            f.close()

        if prunenonlenses==False:
          if save:
            f=open("idealisedlenses/nonlenspopulation_%s.pkl"%self.sourcepopulation,'wb')
            pickle.dump(self.lens, f, protocol=5) # NH: update protocol
            f.close()
            print(len(self.lens.keys()))


        self.lens=self.reallens

    def LoadLensPop(self,j=0,sourcepopulation="jaguar"):
        f=open("idealisedlenses/lenspopulation_%s_%i.pkl"%(sourcepopulation,j),'rb')
        self.lens=pickle.load(f)
        f.close()


    def apply_preselection(self,imag,z):
        if imag<15: return False
        if imag>23:return False
        if z<0.05: return False
        return True

if __name__ == "__main__":

    import distances

    D = distances.Distance()

    fsky = 1

    Lpop = LensPopulation(reset=True, sigfloor=100, zlmax=3, D=D)

    Ndeflectors = Lpop.Ndeflectors(z=3, zmin=0, fsky=fsky)

    L = LensSample(sigfloor=100, cosmo=[0.3,0.7,0.7], sourcepop="jaguar")

    L.Generate_Lens_Pop(int(Ndeflectors), firstod=1, nsources=1, prunenonlenses=True)
