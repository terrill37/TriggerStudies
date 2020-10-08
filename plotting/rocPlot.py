import ROOT

ROOT.gROOT.SetBatch(True)

import ROOTHelp.FancyROOTStyle
#from JetLevelPlotUtils import getCMSText, getText

from optparse import OptionParser
p = OptionParser()

#change input to txt file with list of input files
p.add_option('--input',  type = 'string', default = "outBTag.FTKBtagging.ttbar.mwt2.All.root", dest = 'inFile', help = 'input File' )

p.add_option('--output', type = 'string', default = "makeRocCurves", dest = 'outDir', help = 'output dir' )

p.add_option('--nameTag', type = 'string', default = "", dest = 'nameTag')

p.add_option('--runEta', action='store_true', dest='runEta', default=False, help='runs eta plot')

p.add_option('--workingPoint', dest='workingPt', action='append', type='float',
            default=[0.2], help='working Point cuts')

p.add_option('--lightDir' , dest='lightDir', type='string', default="matchedJet_L", 
             help="light flavor directory to be used")

#p.add_option('--cmsText', type = 'string', default = "Work in Progress",  help = '' )

#p.add_option('--doCaloJets', action="store_true", help = '' )
(o,a) = p.parse_args()

#check if input file txt file
#make list of input files
file_list=[]
ptCuts = ['0p4', '0p9', '1p5', '2p0', '2p5', '5p0', '7p5', '10p0', '20p0']
infile_pt =[]
PU_type=[]
if (".txt" in o.inFile):
    for line in open(o.inFile, "r").readlines():
        line = line.replace('\n','').strip()
        #print(line)
        if line    == '' : continue
        if line[0] == '#': continue
        file_list.append(line.replace('\n',''))

        #find PU200 or NoPU
        if "PU200" in line:
            PU_type.append("PU200")

        elif "NoPU" or "noPU" in line:
            PU_type.append("NoPU")

        else:
            PU_type.append("none")

        for pt in ptCuts:
            if (pt in line):
                infile_pt.append(pt.replace('p','.'))
    
else:
    file_list.append(o.inFile)

#print(infile_pt)

inFile = []
for file_name in file_list:
   # print(file_name)
    inFile.append(ROOT.TFile(file_name, "READ"))
#print(inFile)

import os
if not os.path.exists(o.outDir):
    os.mkdir(o.outDir)

from rocCurveUtils     import makeRoc
from JetLevelPlotUtils import getCMSText, getText

def computeSigEff(sigNum, sigDen, bkgNum, bkgDen, workingPoint=0.7, debug=False):
    #print(inFile)
    #bkg light flavor
    #sig b-jet
    
    nbins = sigNum.GetNbinsX() #number of bins
    #print("nbins",nbins)

    sigTot=sigDen.Integral(0,nbins+1) #total signal
    #print ("total sig", sigTot)
    bkgTot=bkgDen.Integral(0,nbins+1) #total background
    #print ("total bkg", bkgTot)
   
    for bins in range(nbins, 0, -1):
        if debug: print("bins", bins)

        thisSig=abs(sigNum.Integral(bins, nbins+1))
        thisBkg=abs(bkgNum.Integral(bins, nbins+1))
        
        sigEff= thisSig/sigTot #b-jet
        bkgEff= thisBkg/bkgTot #light flavor

        if bkgEff:
            bkgRej = 1.0/bkgEff
        if bkgEff > workingPoint:
            #print("thisSig: ", thisSig)
            if debug: print("sigEff: ", sigEff,"bkgEff ", bkgEff)
            #cut = inFile.GetName()[:inFile.GetName().find(".")].replace("/","_")
            return sigEff, bkgEff

        #else:
            #return 1, 1

def makeWorkingPointsHists(workingPoint, cuts, sigEffs):
    can = ROOT.TCanvas("can", "B-Jet Eff at L-flavor " + str(workingPoint))
    
    can.SetGrid()
    graph = ROOT.TGraph(len(cuts))
    
    #print('before makeWorkingPointsHists loop')
    for i in range(len(cuts)):
        #print(i)
        graph.SetPoint(i, float(cuts[i]), sigEffs[i])

    graph.GetXaxis().SetTitle("trackPtCut/GeV")
    graph.GetYaxis().SetTitle("B-Jet Efficiency")

    graph.SetLineColor(3)
    graph.SetLineWidth(4)
    graph.SetMarkerColor(5)
    graph.SetMarkerStyle(21)
    graph.SetTitle("B-Jet Eff as L-flavor set to " + str(workingPoint))
    graph.Draw("ALP")

    can.Update()
    can.GetFrame().SetFillColor(21)
    can.GetFrame().SetBorderSize(12)
    can.Modified
    can.Update()
    can.SaveAs(o.outDir+"/workingPoint_"+str(workingPoint).replace('.','p')+\
                o.nameTag+"_"+o.lightDir+".png")
    can.SaveAs(o.outDir+"/workingPoint_"+str(workingPoint).replace('.','p')+\
                o.nameTag+"_"+o.lightDir+".pdf")


def plotSame(name, graphs, colors, styles, pt_name, legend_names, plotCaloJet=False, plotPFJet=False, plotOffJet=False, plotCSV=False, plotDeepCSV=False, workingPts=None, rocType=None):
    
    can = ROOT.TCanvas(name,name)
    can.cd().SetLogy(0)
 #   print(graphs)
    legend = ROOT.TLegend(0.2, 0.9, 0.4, 0.6)
    
    for gItr, g in enumerate(graphs):
        key = legend_names[gItr]+' GeV'
        legend.AddEntry(g, key, 'L')
        legend.SetLineWidth(0)
        
        

    for gItr, g in enumerate(graphs):
        g.SetLineColor(colors[gItr])
        g.SetLineStyle(styles[gItr])
        g.GetXaxis().SetRangeUser(0.6, 1.0)
        g.GetYaxis().SetRangeUser(0.0, 0.6)
        if not gItr:
            g.Draw("AL")
        else:
            g.Draw("L")
    legend.Draw("same")

    #cmsLine1, cmsLine2 = getCMSText(xStart=0.6, yStart=0.875, subtext=o.cmsText)
    #cmsLine1.Draw("same")
    #cmsLine2.Draw("same")

    can.SaveAs(o.outDir+"/roc_combo"+name+\
                o.nameTag+"_"+o.lightDir+".pdf")
    can.SaveAs(o.outDir+"/roc_combo"+name+\
                o.nameTag+"_"+o.lightDir+".png")

#def plotEta(direc, bkg, sig, inFile=[]):
#    #make plots of eta for each ptCut
#    can_L=ROOT.TCanvas()
#    can_B=ROOT.TCanvas()
#    
#    legend_L=ROOT.TLegend(0.2, 0.9, 0.3, 0.6)
#    i=0
#
#    eta_L=[]
#    eta_B=[]
#    colors=[]
#    for files in inFile:
#        eta_L.append(files.Get(direc+"_"+bkg+"/eta"))
#        eta_B.append(files.Get(direc+"_"+sig+"/eta"))
#        colors.append(i)
#        i = i+1
#
#    for j, eta_light in enumerate(eta_L):
#        can_L.cd()
#        eta_light.Draw("h,same")
#        eta_light.SetLineColor(colors[j]+1)
#
#    for k, eta_bjet in enumerate(eta_B):
#        can_B.cd()
#        eta_bjet.Draw("h,same")
#        eta_bjet.SetLineColor(colors[k]+1)
#    
#    print("Nbins of eta_L", eta_L[1].GetNbinsX())
#    print("Nbins of eta_B", eta_B[0].GetNbinsX())
#
#    can_L.SaveAs(o.outDir+"/"+o.nameTag+"light_Eta.png")
#    can_B.SaveAs(o.outDir+"/"+o.nameTag+"b_Eta.png")

name = "PF_deepcsv"
var = "DeepCSV_l"
bkg = o.lightDir 
#"matchedJet_L"  
sig = "matchedJet_B"

varNorm=None

dire="offJets"

color =[]
plotsEff =[]
plotsRej =[]
line_type =[]
legend_names=[]

sigEffs=[]
bkgEffs=[]

sigHists=[]
bkgHists=[]
sigNormHists=[]
bkgNormHists=[]

i = 0
for flea in inFile:
    sigHist = flea.Get(dire+"_"+sig+"/"+var)
    bkgHist = flea.Get(dire+"_"+bkg+"/"+var)

    sigHists.append(sigHist)
    bkgHists.append(bkgHist)
    
    legend_names.append(PU_type[i]+" pt "+str(infile_pt[i]))

    if PU_type[i] == "PU200":
        color.append(i+1)
        line_type.append(ROOT.kSolid)

    elif PU_type[i] == "NoPU":
        color.append(i+1)
        line_type.append(7)

    else:
        color.append(i+1)
        line_type.append(6)

    if varNorm:
        sigNormHist = flea.Get(dire+"_"+sig+"/"+varNorm)
        bkgNormHist = flea.Get(dire+"_"+bkg+"/"+varNorm)

    else:
        sigNormHist = sigHist
        bkgNormHist = bkgHist
    
    sigNormHists.append(sigNormHist)
    bkgNormHists.append(bkgNormHist)

    rocPlots =[]
    #i=0
    for config in [("Rej", 1, 5e4),("Eff",5e-4,1)]:
        #print("before error")
        rocPlots.append(makeRoc(sigHist, sigNormHist, bkgHist, bkgNormHist, doErr=False, bkgMode=config[0], cleanNoCut=True, debug=False))
        can = ROOT.TCanvas(name+"_"+str(i)+config[0], name+"_"+str(i)+config[0])
        can.cd().SetLogy(1)
        rocPlots[-1].SetLineWidth(5)
        rocPlots[-1].GetXaxis().SetTitle("B-Jet Efficiency")
        rocPlots[-1].GetXaxis().SetRangeUser(0.4,1)
        if config[0]=="Rej":
            ytitle = "Light flavor rejection"
        elif config[0]=="Eff":
            ytitle = "Light flavor efficiency"
        rocPlots[-1].GetYaxis().SetTitle(ytitle)
        rocPlots[-1].GetYaxis().SetRangeUser(config[1],config[2])
        rocPlots[-1].Draw("AL, same")
        #can.SaveAs(o.outDir+"/roc_"+name+"_"+config[0]+str(i)+".pdf")
        #print("\n") 
    i=i+1

    plotsEff.append(rocPlots[1])
    plotsRej.append(rocPlots[0])
    #print(flea)

#print(plotsEff)

#print("start plotSame")
#print(plots[1][:])
#print(plots[:][1])
plotSame("PF_deepCSV_"+"Eff",
         plotsEff,
         color,
         line_type,
         pt_name = infile_pt,
         legend_names=legend_names,
         plotOffJet = True,
         plotPFJet = True,
         plotDeepCSV = False,
         rocType = "Eff"
         )
        

for workPt in o.workingPt:
    sigEffsIn=[]
    bkgEffsIn=[]
    for i in range(0, len(sigHists), 1):
        sigEff, bkgEff = computeSigEff(sigNum=sigHists[i], sigDen=sigNormHists[i], \
                                       bkgNum=bkgHists[i], bkgDen=bkgNormHists[i], \
                                       workingPoint=workPt)

        sigEffsIn.append(sigEff) #add to inner list
        bkgEffsIn.append(bkgEff)
    
    sigEffs.append(sigEffsIn) #add list to list
    bkgEffs.append(bkgEffsIn) #may not need

    makeWorkingPointsHists(workingPoint=workPt, cuts=infile_pt, sigEffs=sigEffsIn)        

#if (o.runEta):
#    plotEta(direc=dire, bkg=bkg, sig=sig, inFile=inFile)

