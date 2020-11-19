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
effBinning = 10


eff_mass_wDeepCut   = makeEff("mass" , 
    ["L1_wDeepCSVCut",  "deepCSV_noL1"]  ,
   inFileMC, binning=effBinning)

eff_mass_trig1   = makeEff("mass" , 
    ["HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1",  "deepCSV_noL1"]  ,
    inFileMC, binning=effBinning)

eff_mass_trig2   = makeEff("mass" , 
    ["HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1",  "deepCSV_noL1"]  ,
    inFileMC, binning=effBinning)

eff_mass_trig3   = makeEff("mass",
    ["HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",  "deepCSV_noL1"]  ,
    inFileMC, binning=effBinning)

drawComp("mass_Efficiency2",[(eff_mass_wDeepCut, "L1 only", ROOT.kBlue),
                            (eff_mass_trig1, "HLT_QuadPFPuppiJet", ROOT.kRed),
                            (eff_mass_trig2, "TriplePFPuppiBTagDeepCSV", ROOT.kGreen),
                            (eff_mass_trig3, "QuadJetOnly", ROOT.kMagenta)
                            ],
         yTitle="Efficiency", xTitle="m_{hh4b} (GeV)", otherText="DeepCSV > 0.3",
         xMin = -0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
         outDir=o.outDir, cmsText="preliminary", lumiText="")



can = ROOT.TCanvas()

no_cuts = 'noCuts'
L1_wDeepCSVCut = 'L1_wDeepCSVCut'
deepCSV_noL1 = 'deepCSV_noL1'
trig1 = 'HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1'
trig2 = 'HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1'
trig3 = 'HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1'

hist_names=["mass","deepCSV"]

hist2 = inFileMC.Get(L1_wDeepCSVCut+'/mass')
hist= inFileMC.Get(deepCSV_noL1+'/mass')
hist3= inFileMC.Get(trig1+'/mass')
hist4= inFileMC.Get(trig2+'/mass')
hist5= inFileMC.Get(trig3+'/mass')

can.cd()
hist.SetLineColor(ROOT.kBlack)
hist.SetMarkerColor(ROOT.kBlack)

hist2.SetLineColor(ROOT.kBlue)
hist2.SetMarkerColor(ROOT.kBlue)

hist3.SetLineColor(ROOT.kRed)
hist3.SetMarkerColor(ROOT.kRed)

hist4.SetLineColor(ROOT.kGreen)
hist4.SetMarkerColor(ROOT.kGreen)

hist5.SetLineColor(ROOT.kMagenta)
hist5.SetMarkerColor(ROOT.kMagenta)

hist.Draw()
hist2.Draw('same')
hist3.Draw('same')
hist4.Draw('same')
hist5.Draw('same')

hist.GetYaxis().SetTitle("Entries")
hist.GetXaxis().SetTitle("m_{hh4b} (GeV)")

leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(hist, 'no trigger', 'L')
leg.AddEntry(hist2, 'L1 trigger', 'L')
leg.AddEntry(hist3, 'QuadPFPuppiJet', 'L')
leg.AddEntry(hist4, 'TriplePFPuppiBTagDeepCSV', 'L')
leg.AddEntry(hist5, 'Quad Only', 'L')
leg.Draw('same')

can.SaveAs(o.outDir+'/masses2.png')

can.Clear()

hist_deepCSV = inFileMC.Get(no_cuts+'/deepCSV')
hist_deepCSV.Draw()
hist_deepCSV.GetXaxis().SetTitle("DeepCSV score")
can.SetLogy(1)
can.SaveAs(o.outDir+'DeepCSV_precut2.png')
can.Clear()

hist_deepCSV2 = inFileMC.Get('L1_wDeepCSVCut'+'/deepCSV')
hist_deepCSV2.Draw()
hist_deepCSV2.GetXaxis().SetTitle("DeepCSV score")
can.SaveAs(o.outDir+'DeepCSV_trigCut2.png')
can.Clear()


