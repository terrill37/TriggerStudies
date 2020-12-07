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

def efficiency_plots(file_name, plot_name, eff_num, eff_denom, effBinning):
    eff_plots = []
    for i in range(0, len(eff_num), 1):
        eff = makeEff(plot_name,
                     [eff_num[i][0], eff_denom],
                     inFileMC, binning=effBinning)

        eff_info = (eff, eff_num[i][1], eff_num[i][2])

        eff_plots.append(eff_info)

    drawComp(plot_name+"_Efficiency"+file_name, eff_plots,
             yTitle="Efficiency", xTitle=plot_name+" (GeV)", otherText="",
             xMin = -0.2, xMax = 1000, yMax = 1.2, xStartOther=600, yStartOther=0.3,
             outDir=o.outDir, cmsText="preliminary", lumiText="")


def plots(plot_name, plot_dir, file_name=""):
    can = ROOT.TCanvas()

    plots=[]
    leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    
    for i in range(0, len(plot_dir), 1):
        hist = inFileMC.Get(plot_dir[i][0] + '/' + plot_name)
        
        hist.SetLineColor(plot_dir[i][2])
        hist.SetMarkerColor(plot_dir[i][2])

        hist.GetYaxis().SetTitle("Entries")
        hist.GetXaxis().SetTitle(plot_name+" (GeV)")

        leg.AddEntry(hist, plot_dir[i][1], 'L')

        plots.append(hist)
    
    can.cd()
    
    plots[0].Draw()
    for j in range(1, len(plots), 1):
        #print(j)
        plots[j].Draw('same')

    leg.Draw('same')

    can.SaveAs(o.outDir+'/'+plot_name+file_name+'.png')
    can.Clear()




def main():
    effBinning = 5
    
    print("mass efficiency")
    tagged_info   = [("tagged_L1","L1", ROOT.kBlack),
                     ("tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1", "QuadOnly",ROOT.kRed),
                     ("tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1", "Quad + HT",ROOT.kMagenta),
                     ("tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1","Quad + Ht + TriplePuppi", ROOT.kBlue),
                     ]
    
    untagged_info = [("noCuts", "noL1", ROOT.kGreen),
                     ("L1_untagged","L1", ROOT.kBlack),
                     ("HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1", "QuadOnly",ROOT.kRed),
                     ("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1", "Quad + HT",ROOT.kMagenta),
                     ("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1","Quad + Ht + TriplePuppi", ROOT.kBlue),
                     ]

    
    eff_plots_den = "tagged_noL1"

    efficiency_plots("_tagged", "mass", tagged_info, eff_plots_den, effBinning)
    
    #print(tagged_info)
    tagged_info.insert(0, ("tagged_noL1", "noL1", ROOT.kGreen))
    #print(tagged_info)
    

    print("tagged plots and untagged plots")
    plot_names = ["mass", "pt_1", "pt_2", "pt_3", "pt_4", "Ht"]
    for i in range(0, len(plot_names), 1):
        plots(plot_names[i], tagged_info,   "_tagged")
        plots(plot_names[i], untagged_info, "_untagged")

    print("pt efficiencies of tagged")
    pts=["pt_1_cut","pt_2_cut","pt_3_cut","pt_4_cut"]
    pt_plot_num = [("tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1", "QuadOnly",ROOT.kRed)]
    pt_plot_den = "tagged_L1"
    for j in range(0, len(pts), 1):
        efficiency_plots("_tagged", pts[j], pt_plot_num, pt_plot_den, effBinning)
   
    



if __name__ == "__main__":
    main()



