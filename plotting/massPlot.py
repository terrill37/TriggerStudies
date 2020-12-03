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
plot_names = ["mass"]

#["mass", "Ht", "pt_1",
# "pt_2", "pt_3", "pt_4"]


print("mass efficiency")
for i in range(0,len(plot_names)):    
    eff_mass_wDeepCut   = makeEff(plot_names[i] , 
        ["tagged_L1",  "tagged_noL1"]  ,
       inFileMC, binning=effBinning)
    
    eff_mass_trig1   = makeEff(plot_names[i] , 
        ["tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1",  "tagged_noL1"]  ,
        inFileMC, binning=effBinning)
    
    eff_mass_trig2   = makeEff(plot_names[i] , 
        ["tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1",  "tagged_noL1"]  ,
        inFileMC, binning=effBinning)
    
    eff_mass_trig3   = makeEff(plot_names[i],
        ["tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",  "tagged_noL1"]  ,
        inFileMC, binning=effBinning)
    
    drawComp(plot_names[i]+"_Efficiency",[(eff_mass_wDeepCut, "L1 only", ROOT.kBlue),
                                (eff_mass_trig1, "HLT_QuadPFPuppiJet", ROOT.kRed),
                                (eff_mass_trig2, "TriplePFPuppiBTagDeepCSV", ROOT.kGreen),
                                (eff_mass_trig3, "QuadJetOnly", ROOT.kMagenta)
                                ],
             yTitle="Efficiency", xTitle=plot_names[i]+" (GeV)", otherText="DeepCSV > 0.3",
             xMin = -0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
             outDir=o.outDir, cmsText="preliminary", lumiText="")

plots = ["mass", "Ht", "pt_1",
        "pt_2", "pt_3", "pt_4"]

for i in range(0, len(plots)):
    can = ROOT.TCanvas()
    
    no_cuts = 'noCuts'
    L1_wDeepCSVCut = 'tagged_L1'
    deepCSV_noL1 = 'tagged_noL1'
    trig2 = 'tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1'
    trig3 = 'tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1'
    trig1 = 'tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1'
    
    hist2 = inFileMC.Get(L1_wDeepCSVCut+'/'+plots[i])
    hist= inFileMC.Get(deepCSV_noL1+'/'+plots[i])
    hist3= inFileMC.Get(trig1+'/'+plots[i])
    hist4= inFileMC.Get(trig2+'/'+plots[i])
    hist5= inFileMC.Get(trig3+'/'+plots[i])
    
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
    hist.GetXaxis().SetTitle(plots[i]+" (GeV)")
    
    leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(hist, 'no trigger', 'L')
    leg.AddEntry(hist2, 'L1 trigger', 'L')
    leg.AddEntry(hist3, 'Quad', 'L')
    leg.AddEntry(hist4, 'Ht', 'L')
    leg.AddEntry(hist5, 'triplePuppi', 'L')
    leg.Draw('same')
    
    can.SaveAs(o.outDir+'/'+plots[i]+'.png')
    
    can.Clear()

hist_deepCSV = inFileMC.Get(no_cuts+'/deepCSV')
hist_deepCSV.Draw()
hist_deepCSV.GetXaxis().SetTitle("DeepCSV score")
can.SetLogy(1)
can.SaveAs(o.outDir+'/DeepCSV_precut.png')
can.Clear()

hist_deepCSV2 = inFileMC.Get('tagged_L1'+'/deepCSV')
hist_deepCSV2.Draw()
hist_deepCSV2.GetXaxis().SetTitle("DeepCSV score")
can.SaveAs(o.outDir+'/DeepCSV_trigCut.png')
can.Clear()

trig_plots=['pt_all','pt_cut','pt_initial']

for i in range(0,len(trig_plots)):
    h = inFileMC.Get('trigger/'+trig_plots[i])
    h.Draw()
    h.GetXaxis().SetTitle(trig_plots[i]+' (GeV)')
    can.SetLogy(0)
    can.SaveAs(o.outDir+'/'+trig_plots[i]+'.png')
    can.Clear()


eff_ht   = makeEff("Ht",
    ["tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1","tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1"],
     inFileMC, binning=effBinning)
 
#eff_mass_trig3   = makeEff(plot_names[i],
 #   ["HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",  "deepCSV_noL1"]  ,
  #   inFileMC, binning=effBinning)
 
drawComp("tagged_Ht_Efficiency",[(eff_ht, "Ht to Quad", ROOT.kRed),
                                     ],
             yTitle="Efficiency", xTitle="Ht (GeV)", otherText="DeepCSV > 0.3",
             xMin = -0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
             outDir=o.outDir, cmsText="preliminary", lumiText="")

effs_all=[]
for i in range(1,5,1):
    eff = makeEff("pt_"+str(i),
                  ["tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1", "tagged_L1"],
                  inFileMC, binning=effBinning)

    eff_name = (eff, "pt_"+str(i), 2*i)

    effs_all.append(eff_name)

    drawComp("pt_"+str(i)+"_ratio_Quad_Only", [eff_name],
              yTitle="Efficiency", xTitle="pt (GeV)", otherText="DeepCSV > 0.3",
              xMin=-0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
              outDir=o.outDir, cmsText="preliminary", lumiText="Quad Only")

drawComp("pt_ratio_Quad_Only", effs_all,
         yTitle="Efficiency", xTitle="pt (GeV)", otherText="DeepCSV > 0.3",
         xMin=-0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
         outDir=o.outDir, cmsText="preliminary", lumiText="Quad Only")


L1_untagged = 'L1_untagged'
trig1_untagged = 'HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1'
trig2_untagged = 'HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1'
trig3_untagged = 'HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1'

#hist= inFileMC.Get(deepCSV_noL1+'/'+plots[i])
hist2= inFileMC.Get(L1_untagged +'/Ht')
hist3= inFileMC.Get(trig1_untagged+'/Ht')
hist4= inFileMC.Get(trig2_untagged+'/Ht')
hist5= inFileMC.Get(trig3_untagged+'/Ht')

can.cd()
#hist.SetLineColor(ROOT.kBlack)
#hist.SetMarkerColor(ROOT.kBlack)

hist2.SetLineColor(ROOT.kBlue)
hist2.SetMarkerColor(ROOT.kBlue)

hist3.SetLineColor(ROOT.kRed)
hist3.SetMarkerColor(ROOT.kRed)

hist4.SetLineColor(ROOT.kGreen)
hist4.SetMarkerColor(ROOT.kGreen)

hist5.SetLineColor(ROOT.kMagenta)
hist5.SetMarkerColor(ROOT.kMagenta)

#hist.Draw()
hist2.Draw('')
hist3.Draw('same')
hist4.Draw('same')
hist5.Draw('same')

hist.GetYaxis().SetTitle("Entries")
hist.GetXaxis().SetTitle("Ht (GeV)")

leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
#leg.AddEntry(hist, 'no trigger', 'L')
leg.AddEntry(hist2, 'L1', 'L')
leg.AddEntry(hist3, 'Quad', 'L')
leg.AddEntry(hist4, 'Ht', 'L')
leg.AddEntry(hist5, 'triplePuppi', 'L')
leg.Draw('same')

can.SaveAs(o.outDir+'/Ht_untagged.png')

can.Clear()


