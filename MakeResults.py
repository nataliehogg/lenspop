from __init__ import *
import pickle
import sys,os
import pylab as plt
import glob

sourcepops=["jaguar"]

if len(sys.argv)>1:
    experiment=sys.argv[1]

surveystoread=[]
if experiment=="Euclid":
    surveystoread+=["Euclid"]
elif experiment=="CFHT":
    surveystoread+=["CFHT"]
elif experiment=="CFHTa":
    surveystoread+=["CFHTa"]
elif experiment=="DES":
    surveystoread+=["DESc"]
    surveystoread+=["DESb"]
    surveystoread+=["DESa"]
elif experiment=="LSST":
    surveystoread+=["LSSTc"]
    surveystoread+=["LSSTb"]
    surveystoread+=["LSSTa"]
else:
    surveystoread=[str(experiment)]
    # experiment=experiment[:-1]

for survey in surveystoread:
  for sourcepop in sourcepops:
    if survey[-2]=="a":
        surveyname=survey[:-1]+"_full_coadd"
    elif survey[-2]=="b":
        surveyname=survey[:-1]+"_best_epoch"
    elif survey[-2]=="c":
        surveyname=survey[:-1]+"_optimal_coadd"
    else:
        surveyname=survey

    filename = "{}_{}_lists.pkl".format(survey, sourcepop)

    print(filename)

    lensparsfile = "lenses_{}.txt".format(survey)

    print(lensparsfile)

    f = open(lensparsfile,"w")
    # print
    #os.system("rm %s"%filename) #this line resets the read-in
    bl={}
    zs={}
    zl={}
    sigl={}
    ql={}
    rs={}
    ms={}
    mag={}
    weights={}
    for key in ["resolved", "rfpf"]:
        bl[key]=[]
        zs[key]=[]
        rs[key]=[]
        ms[key]=[]
        zl[key]=[]
        sigl[key]=[]
        ql[key]=[]
        mag[key]=[]
        rs[key]=[]
        weights[key]=[]

    if experiment=="CFHT":
      frac=42000.*1./150.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    if experiment=="CFHTa":
      frac=42000.*1./150.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    elif experiment=="Euclid":
      frac=42000.*1./15000.
      bands=["VIS"]

    elif experiment=="DES":
      frac=42000.*1./5000.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    elif experiment=="LSST":
      frac=42000.*1./20000.
      bands=["g_SDSS","r_SDSS","i_SDSS"]

    elif experiment=="COSMOS-Web":
        frac=41253.*1./0.54
        bands=['JWST_NIRCam_F115W', 'JWST_NIRCam_F150W', 'JWST_NIRCam_F277W', 'JWST_NIRCam_F444W']

    filelist=glob.glob("LensStats/{}_{}_Lens_stats_*.pkl".format(experiment,sourcepop))

    chunki=0
    ilist=[]
    # print survey
    print('running {}'.format(survey))
    for chunk in filelist:
        # print chunki
        # print(chunki)
        chunki+=1
        f2=open(chunk,"rb")
        # fracsky,sspl=cPickle.load(f2)
        fracsky,sspl=pickle.load(f2)
        fract=frac*fracsky
        f2.close()
        I=0
        # print(sspl[1])
        for i in sspl.keys():
            if i in ilist:
                continue
            else:
                try:
                    sspl[i]["seeing"][survey]
                except KeyError:
                    continue
                f.write("%.2f "%sspl[i]["zl"])
                f.write("%.2f "%sspl[i]["zs"][1])
                f.write("%.2f "%sspl[i]["b"][1])
                f.write("%.2f "%sspl[i]["sigl"])
                f.write("%.2f "%sspl[i]["ql"])
                #f.write("%.2f "%sspl[i]["rl"]["g_SDSS"])
                for band in bands:
                    f.write("%.2f "%sspl[i]["ml"][band])
                #f.write("%.2f "%sspl[i]["rl"]["g_SDSS"])
                f.write("%.2f "%sspl[i]["xs"][1])
                f.write("%.2f "%sspl[i]["ys"][1])
                f.write("%.2f "%sspl[i]["qs"][1])
                f.write("%.2f "%sspl[i]["ps"][1])
                f.write("%.2f "%sspl[i]["rs"][1])
                f.write("%.2f "%sspl[i]["mag"][1])
                for band in bands:
                    f.write("%.2f "%sspl[i]["seeing"][survey][band])
                    f.write("%.2f "%sspl[i]["SN"][survey][1][band][0])
                if survey!="Euclid":
                    f.write("%.2f "%sspl[i]["rfsn"][survey][1][0])
                f.write("\n")


                ilist.append(str(i))
                if sspl[i]["pf"][survey][1]==False:continue

                try:
                    bb=sspl[i]["bestband"][survey][1]
                    #print sspl[i]["seeing"][survey][bb]
                    #print sspl[i]["mag"][1]*sspl[i]["rs"][1],
                    try:
                        (sspl[i]["b"][1]**2-sspl[i]["rs"][1]**2)**0.5
                    except FloatingPointError: print(0)
                except KeyError:
                  pass
                try:
                  if sspl[i]["resolved"][survey][1][sspl[i]["bestband"][survey][1]]:
                    bb=sspl[i]["bestband"][survey][1]
                    if sspl[i]["mag"][1]<3:continue
                    if sspl[i]["SN"][survey][1][bb][0]<20:continue

                    bl["resolved"].append(sspl[i]["b"][1])
                    weights["resolved"].append(1./fract)
                    zs["resolved"].append(sspl[i]["zs"][1])
                    rs["resolved"].append(sspl[i]["rs"][1])
                    zl["resolved"].append(sspl[i]["zl"])
                    sigl["resolved"].append(sspl[i]["sigl"])
                    ql["resolved"].append(sspl[i]["ql"])
                    mag["resolved"].append(sspl[i]["mag"][1])
                    #ms["resolved"].append(sspl[i]["ms"][1]["g_SDSS"])

                    if sspl[i]["rfpf"][survey][1]:
                        if sspl[i]["rfsn"][survey][1][0]<20:continue
                        if sspl[i]["resolved"][survey][1]["RF"]==False:continue

                        if experiment=="CFHT" or experiment=="CFHTa":
                            if sspl[i]["zl"]>1:continue
                            if sspl[i]["zl"]<0.2:continue
                            if sspl[i]["ml"]["i_SDSS"]<17:continue
                            if sspl[i]["ml"]["i_SDSS"]>22:continue

                        bl["rfpf"].append(sspl[i]["b"][1])
                        weights["rfpf"].append(1./fract)
                        zs["rfpf"].append(sspl[i]["zs"][1])
                        rs["rfpf"].append(sspl[i]["rs"][1])
                        zl["rfpf"].append(sspl[i]["zl"])
                        sigl["rfpf"].append(sspl[i]["sigl"])
                        ql["rfpf"].append(sspl[i]["ql"])
                        mag["rfpf"].append(sspl[i]["mag"][1])
                        #ms["rfpf"].append(sspl[i]["ms"][1]["g_SDSS"])


                except KeyError:
                  pass
    f.close()

    if survey[-2]=="a":
        surveyname=survey[:-1]+" (full coadd)"
    elif survey[-2]=="b":
        surveyname=survey[:-1]+" (best single epoch imaging)"
    elif survey[-2]=="c":
        surveyname=survey[:-1]+" (optimal coadd)"
    else:
        surveyname=survey

    nlenses = numpy.sum(numpy.array(weights["resolved"]).ravel())
    nlenses_gi = numpy.sum(numpy.array(weights["rfpf"]).ravel())

    print('{} will find {} lenses assuming Poisson limited galaxy subtraction in all bands.'.format(survey, nlenses))

    f=open(filename,"wb")
    pickle.dump([weights, bl, zs, rs, ms, zl, sigl, ql, mag], f, 2)
    f.close()
