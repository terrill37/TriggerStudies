import ROOT

ROOT.gROOT.SetBatch(True)

import ROOTHelp.FancyROOTStyle

from optparse import OptionParser
p = OptionParser()
p.add_option('--inputMC',  type = 'string', dest = 'inFileMC', help = 'input MC File' )
p.add_option('--output',   type = 'string', dest = 'outDir'  , help = 'output dir')

(o,a) = p.parse_args()

inFileMC = ROOT.TFile(o.inFileMC, "READ")

import os
if not os.path.exists(o.outDir):
    os.mkdir(o.outDir)

from JetLevelPlotUtils import makeEff, drawComp, getHist, drawStackCompRatio, makeStack, makeInverseTurnOn, make2DComp, makeInverseTurnOnAll, plotRatio

#
#  Offline Turnon curves:
#
effBinning = 5

eff_PuppiDeepCSV_MC   = makeEff("DeepCSV_l" , 
    ["offJets_matchedPuppiDeepcsvTag",  "offJets_matchedPuppi"]  ,
    inFileMC, binning=effBinning)

#eff_PuppiDeepCSVJet_MC = makeEff("DeepCSV_l",
#    ["offJets_matchedPuppiDeepcsvTagJet", "offJets_matchedPuppiJet"],
#    inFileMC, binning=effBinning)
#OnlineDeepCSVCutPuppi  = 0.17


reveff_DeepCSV_MC = {}
jetTypes = "Puppi"

reveff_DeepCSV_MC[jetTypes] = {}
for op in ["Loose", "Medium", "Tight"]:
    reveff_DeepCSV_MC [jetTypes][op] = makeEff("DeepCSV_l",
                                                ["offJets"+op+"_matched"+jetTypes+"Jet",
                                                "offJets_matched"+jetTypes+"Jet"],
                                                inFileMC, binning=effBinning)

drawComp("RevEff_DeepCSV_MC", [(reveff_DeepCSV_MC[jetTypes]["Loose"],  "Loose",  ROOT.kBlue),
                               (reveff_DeepCSV_MC[jetTypes]["Medium"], "Medium", ROOT.kSpring),
                               (reveff_DeepCSV_MC[jetTypes]["Tight"],  "Tight",  ROOT.kRed),],
         yTitle="Efficiency", xTitle="DeepCSV Value of Jets", otherText=""+jetTypes+" Jets(MC)", 
         outDir=o.outDir, cmsText="work in progress", lumiText="testing")

drawComp("Eff_Puppi_DeepCSV_MC", [(eff_PuppiDeepCSV_MC, "matched", ROOT.kBlue),],
 #                                 (eff_PuppiDeepCSVJet_MC, "matched Jet", ROOT.kRed),],
         yTitle="Efficiency", xTitle="DeepCSV Value of Puppi", otherText="MC Puppi",
         outDir=o.outDir, cmsText="DeepCSV Cut > 0.17", lumiText="")


#
# Eff_pt
# pt > 30
#

#cut on offline and online in numerator
#cut on offline only in denominator

#      [name,   xStartOther,yStartOther,units,xmin, xmax,cmsText]
plots=[["pt_m", 200,        0.25,      '(GeV)',0,   500, 'pt<0.5'],
       ["eta",  -2,         0.25,      '',    -5,   5,   ''],
       ["phi",  -1,         0.25,      '',    -3.2, 3.2, ''],]
       
for i in range(0,len(plots)):
    eff_Puppi_MC = makeEff(plots[i][0], 
                            ["offJetsMedDeepCSV_matchedPuppiDeepCSV", 
                            "offJetsMedDeepCSV_matchedPuppiJet"],
                            inFileMC, binning=effBinning)
    
    #eff_PuppiJetPt_MC = makeEff("pt_m",
    #                           ["offJets_matchedPuppiDeepcsvTagJet", "offJets_matchedPuppiJet"],
    #                           inFileMC, binning=effBinning)
    
    print("eff pt passed")
    
    drawComp("Eff_Puppi_"+plots[i][0]+"_MC", [(eff_Puppi_MC, "matched", ROOT.kRed),],
     #                            (eff_PuppiJetPt_MC, "matchedJet", ROOT.kBlue),],
            yTitle="Efficiency", xTitle=plots[i][0]+" Value of Puppi "+plots[i][3],
            otherText="#splitline{offline DeepCSV Medium cut (2017) > 0.4941}{online DeepCSV Cut > 0.17}",
            xMin=plots[i][4], xMax=plots[i][5], yMax=1., xStartOther=plots[i][1], yStartOther=plots[i][2],
            outDir=o.outDir, cmsText=plots[i][6], lumiText="")
    

#
# Offline vs Online
#
csvBinning = 3

offDeepCSV_Puppi_MC   = getHist(inFileMC, "offJets_matchedPuppi", "DeepCSV_l", 
                                binning=csvBinning, norm=1)

puppiDeepCSV_MC       = getHist(inFileMC, "offJets_matchedPuppiJet", "DeepCSV_l",
                                binning=csvBinning, norm=1)

offDeepCSV_B_Puppi_MC = getHist(inFileMC, "offJets_matchedPuppi_B", "DeepCSV_l",
                                binning=csvBinning, norm=0)

offDeepCSV_L_Puppi_MC = getHist(inFileMC, "offJets_matchedPuppi_L", "DeepCSV_l",
                                binning=csvBinning, norm=0)

offDeepCSV_C_Puppi_MC = getHist(inFileMC, "offJets_matchedPuppi_C", "DeepCSV_l",
                                binning=csvBinning, norm=0)

drawStackCompRatio("Off_Puppi_DeepCSV_FlavourComp",(puppiDeepCSV_MC,"MC"),
                   [(offDeepCSV_B_Puppi_MC, "B Jets", ROOT.kAzure),
                    (offDeepCSV_L_Puppi_MC, "Light Flavour", ROOT.kGreen),
                    (offDeepCSV_C_Puppi_MC, "Charm Jets", ROOT.kRed)],
                   yTitle="Normalized", xTitle="Offline (Puppi-Jet) DeepCSV",
                   rTitle="MC", setLogy=0, outDir=o.outDir, 
                   cmsText="work in progress", lumiText ="")


#
# make stacks of pt, eta, phi
#

print ("make stacks")
makeStack("OffJet_Pt", "pt_m", "offJets_matchedPuppi_X", binning=2,
          xTitle="Offline Jet Pt", rTitle="MC", logy=0, inFileData='none', inFileMC=inFileMC,
          outDir=o.outDir, min=20, cmsText="work in progress", lumiText="")

makeStack("OffJet_Eta", "eta", "offJets_matchedPuppi_X", binning=2,
          xTitle="Offline Jet Eta", rTitle="MC", logy=0, inFileData='none', inFileMC=inFileMC,
          outDir=o.outDir, x_min=-5, x_max=5, cmsText="work in progress", lumiText="")

makeStack("OffJet_Phi", "phi", "offJets_matchedPuppi_X", binning=2,
          xTitle="Offline Jet Phi", rTitle="MC", logy=0, inFileData='none', inFileMC=inFileMC,
          outDir=o.outDir, cmsText="work in progress", lumiText="")














