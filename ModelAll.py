from __init__ import * # this imports all the below; we can make this more explicit
#from PopulationFunctions import *
#from MakeLensPop import *
#from Surveys import *
#from FastLensSim import *

import pickle
import sys
import time
import os
import fnmatch
import re

sigfloor = 100 # must match that used in MakeLensPop!
cosmo = [0.3, 0.7, 0.7] # omega_m, omega_L, h; should match sourcepop cosmology?

L = LensSample(sigfloor=sigfloor, cosmo=cosmo)

# default settings
experiment = 'COSMOS-Web'
frac       = 0.1 #0.1 # sky fraction i.e. 10%
a          = 20 # SNR threshold
b          = 3 # Magnification threshold
c          = 1000 # second SNR cut
d          = 1000 # third SNR cut

# user defined settings overwrite above defaults
if len(sys.argv)>1:
    experiment=sys.argv[1]
    frac=float(sys.argv[2])
if len(sys.argv)>3:
    a=int(sys.argv[3])
    b=int(sys.argv[4])
    #c=int(sys.argv[5])
    #d=int(sys.argv[6])

print('sky fraction simulated: {}'.format(frac))
if experiment=='COSMOS-Web':
    area = 0.54
    sky_area = 41253
    weighting = area/(sky_area*frac)
    print('weighting to be applied (COSMOS-Web): {}'.format(weighting))

firstod  = 1 # source overdensity
nsources = 1 # initialise nsources


surveys=[]

if experiment=="Euclid":
    surveys+=["Euclid"]
if experiment=="CFHT":
    surveys+=["CFHT"] #full coadd (Gaussianised)
if experiment=="CFHTa":
    surveys+=["CFHTa"] #dummy CFHT
if experiment=="DES":
    surveys+=["DESc"] #Optimal stacking of data
    surveys+=["DESb"] #Best Single epoch image
    surveys+=["DESa"] #full coadd (Gaussianised)
if experiment=="LSST":
    surveys+=["LSSTc"] #Optimal stacking of data
    surveys+=["LSSTb"] #Best Single epoch image
    surveys+=["LSSTa"] #full coadd (Gaussianised)
    #print "only doing LSSTc"
if experiment=="COSMOS-Web":
    surveys+=["COSMOS-Web"]


S={}
n={}
for survey in surveys:
    S[survey]=FastLensSim(survey,fractionofseeing=1)
    S[survey].bfac=float(2)
    S[survey].rfac=float(2)

# this saves having to remember and hardcode the number of idealised lenses each time
for file in os.listdir('/home/nataliehogg/Documents/Projects/cosmos_web/lenspop/idealisedlenses/'):
# for file in os.listdir('/media/nataliehogg/skygate/lenspop_results/old_results/firstjaguar_idealisedlenses'):
    if fnmatch.fnmatch(file, 'lenspopulation_jaguar_residual_*.pkl'):
        num_jag = int(re.findall('\d+', file)[0])

t0 = time.perf_counter()

for sourcepop in ["jaguar"]:
  chunk=0
  Si=0
  SSPL={}
  foundcount={}
  for survey in surveys:
      foundcount[survey]=0

  if sourcepop=="cosmos":
      n_l2=1100000
  elif sourcepop=="lsst":
      n_l2=12530000
  elif sourcepop=="jaguar":
      n_l2=num_jag # this should be the number appended to the _residual pickle file in idealisedlenses/; it's the total number of real lenses for the whole sky
  nall=int(n_l2*frac)

  for i in range(nall):
    if i%10000==0: # if the remainder of i/10000 is zero i.e. if i is a multiple of 10000:
        print("about to load")
        L.LoadLensPop(i,sourcepop)

    if i!=0:
        if i%10000==0 or i==100 or i==300 or i==1000 or i==3000:
            # t1=time.clock()
            t1=time.perf_counter()
            ti=(t1-t0)/float(i)
            tl=(nall-i)*ti
            tl/=60#mins
            hl=numpy.floor(tl/(60))
            ml=tl-(hl*60)
            # print i,"%ih%im left"%(hl,ml)
            # print(i)
            print('{} hours, {} mins left'.format(hl, ml))

    # print(L)
    lenspars=L.lens[i]
    if lenspars["lens?"]==False:
        del L.lens[i]
        continue

    # don't have these bands for COSMOS-Web
    # lenspars["rl"]["VIS"]=(lenspars["rl"]["r_SDSS"]+\
    #                        lenspars["rl"]["i_SDSS"]+lenspars["rl"]["z_SDSS"])/3
    # for mi in [lenspars["ml"],lenspars["ms"][1]]:
    #     mi["VIS"]=(mi["r_SDSS"]+mi["i_SDSS"]+mi["z_SDSS"])/3

    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["SN"]={}
    lenspars["bestband"]={}
    lenspars["pf"]={}
    lenspars["resolved"]={}
    lenspars["poptag"]={}
    lenspars["seeing"]={}
    lenspars["rfpf"]={}
    lenspars["rfsn"]={}

    lastsurvey="non"
    for survey in surveys:

        S[survey].setLensPars(lenspars["ml"],lenspars["rl"],lenspars["ql"])#,reset=True)
        for j in range(nsources):
            S[survey].setSourcePars(lenspars["b"][j+1],lenspars["ms"][j+1],\
                                    lenspars["xs"][j+1],lenspars["ys"][j+1],\
                                    lenspars["qs"][j+1],lenspars["ps"][j+1],\
                                    lenspars["rs"][j+1],sourcenumber=j+1    )

        if survey[:3]+str(i)!=lastsurvey:
            model=S[survey].makeLens(stochasticmode="MP")
            SOdraw=numpy.array(S[survey].SOdraw)
            if type(model)!=type(None):
                lastsurvey=survey[:3]+str(i)
            if S[survey].seeingtest=="Fail":
                lenspars["pf"][survey]={}
                lenspars["rfpf"][survey]={}
                for src in S[survey].sourcenumbers:
                    lenspars["pf"][survey][src]=False
                    lenspars["rfpf"][survey][src]=False
                continue#try next survey
        else:
            S[survey].loadModel(model)
            S[survey].stochasticObserving(mode="MP",SOdraw=SOdraw)
            if S[survey].seeingtest=="Fail":
                lenspars["pf"][survey]={}
                for src in S[survey].sourcenumbers:
                    lenspars["pf"][survey][src]=False
                continue#try next survey
            S[survey].ObserveLens()

        mag,msrc,SN,bestband,pf=S[survey].SourceMetaData(SNcutA=a,magcut=b,SNcutB=[c,d])
        lenspars["SN"][survey]={}
        lenspars["bestband"][survey]={}
        lenspars["pf"][survey]={}
        lenspars["resolved"][survey]={}
        lenspars["poptag"][survey]=i
        lenspars["seeing"][survey]=S[survey].seeing
        rfpf={}
        rfsn={}
        for src in S[survey].sourcenumbers:
            rfpf[src]=False
            rfsn[src]=[0]
            lenspars["mag"][src]=mag[src]
            lenspars["msrc"][src]=msrc[src]
            lenspars["SN"][survey][src]=SN[src]
            lenspars["bestband"][survey][src]=bestband[src]
            lenspars["pf"][survey][src]=pf[src]
            lenspars["resolved"][survey][src]=S[survey].resolved[src]
        if survey!="Euclid":
            if S[survey].seeingtest!="Fail":
                if survey not in ["CFHT","CFHTa"]:
                    S[survey].makeLens(noisy=True,stochasticmode="1P",SOdraw=SOdraw,MakeModel=False)
                    rfpf,rfsn=S[survey].RingFinderSN(SNcutA=a,magcut=b,SNcutB=[c,d],mode="donotcrossconvolve")
                else:
                    rfpf,rfsn=S[survey].RingFinderSN(SNcutA=a,magcut=b,SNcutB=[c,d],mode="crossconvolve")
        lenspars["rfpf"][survey]=rfpf
        lenspars["rfsn"][survey]=rfsn

        L.lens[i]=None #delete used data for memory saving

    accept=False
    for survey in surveys:
        if lenspars["pf"][survey][1]:
            accept=True

    if accept:
        #S[survey].display(band="VIS",bands=["VIS","VIS","VIS"])
        #if Si>100:exit()
        Si+=1
        SSPL[Si]=lenspars.copy()
        if (Si+1)%1000==0:
            f=open("LensStats/%s_%s_Lens_stats_%i.pkl"%(experiment,sourcepop,chunk),"wb")
            # cPickle.dump([frac,SSPL],f,2)
            pickle.dump([frac,SSPL],f,2)
            f.close()
            SSPL={} # reset SSPL or memory fills up
            chunk+=1

    del L.lens[i]

  f=open("LensStats/%s_%s_Lens_stats_%i.pkl"%(experiment,sourcepop,chunk),"wb")
  # cPickle.dump([frac,SSPL],f,2)
  pickle.dump([frac,SSPL],f,2)
  f.close()
  # print Si
  # print(Si)

bl=[]
for j in SSPL.keys():
    try:
        if SSPL[j]["rfpf"][survey][1]:
            bl.append(SSPL[j]["b"][1])
    except KeyError:pass
