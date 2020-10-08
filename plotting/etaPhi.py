import ROOT
from time import sleep
ROOT.gROOT.SetBatch(True)



import ROOTHelp.FancyROOTStyle
#from JetLevelPlotUtils import getCMSText, getText

from optparse import OptionParser
p = OptionParser()

p.add_option('--input', type = 'string', dest = 'input')

p.add_option('--output', type = 'string', default = "makeRocCurves", dest = 'outDir', help = 'output dir' )

p.add_option('--nameTag', type = 'string', default = "", dest = 'nameTag')

(o,a) = p.parse_args()

import os
if not os.path.exists(o.outDir):
    os.mkdir(o.outDir)



def Read_input(inFile=o.input): 
    print("Begin reading files...")
    
    file_list=[]
    ptCuts = ['0p4', '0p9', '1p5', '2p0', '2p5', '5p0', '7p5', '10p0', '20p0']
    infile_pt =[]
    PU_type=[]
    
    if (".txt" in inFile):
        for line in open(inFile, "r").readlines():
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

    inFile_root = []
    for file_name in file_list:
        # print(file_name)
        inFile_root.append(ROOT.TFile(file_name, "READ"))
    
    sleep(3)
    print("files read ...")

    return inFile_root, PU_type, infile_pt

def Eta_Phi(inFile, infile_pt):
    print("Begin Eta-Phi plotting...")
    print("Retrieving ttree ... ")
    direc_1 = "offJets_matchedJet_B_good"
    direc_2 = "offJets_matchedJet_B_bad"
    
    plot_names = ["etaPhi",
                  "matched_dEta_vs_dPhi",
                  "pixHitMap",
                  "innerPixHitMap"]
    
    for i, plot_name in enumerate(plot_names):
        tree_good = inFile.Get(direc_1+"/tracks/"+plot_name)
        tree_bad = inFile.Get(direc_2+"/tracks/"+plot_name)
        
        #tree_good.SetDirectory(0)
        #tree_bad.SetDirectory(0)

        #print("Closing root file ...")
        #inFile[0].Close()
        #sleep(2)
        
        print("Making canvas ... " + plot_name)
        sleep(1)
        can = ROOT.TCanvas("canvas")
        
        can.cd()

        tree_good.SetStats(0)
        tree_bad.SetStats(0)

        tree_good.GetYaxis().SetTitle("Phi")
        tree_good.GetXaxis().SetTitle("Eta")
        tree_bad.GetYaxis().SetTitle("Phi")
        tree_bad.GetXaxis().SetTitle("Eta")
        
        tree_good.SetTitle("track_good")
        tree_bad.SetTitle("track_bad")

        latex_good = ROOT.TLatex()
        latex_bad = ROOT.TLatex()
        #Latex

        print("drawing good ... " + plot_name)
        pad_good = ROOT.TPad("pad_good", "pad_good", 0, 0, 0.5, 1)
        #pad_good.SetTitle('track_good')
        
        pad_good.Draw()
        pad_good.cd()
        tree_good.Draw("colz")
        tree_good.GetYaxis().SetTitleOffset(1.025)
        sleep(1)

        print("drawing bad ... " + plot_name)
        can.cd()
        pad_bad = ROOT.TPad("pad_bad", "pad_bad", 0.5, 0, 1, 1)
        #pad_bad.SetTitle('track_bad')
        pad_bad.Draw()
        pad_bad.cd()
        tree_bad.Draw("colz")
        
        pad_good.SetRightMargin(.2)
        pad_bad.SetRightMargin(.2)
        pad_good.SetTopMargin(0.9)
        pad_bad.SetTopMargin(0.9)
        sleep(1)

        print("printing canvas ... " + plot_name)
        can.SaveAs(o.outDir+"/" + plot_name + "_pt" + str(infile_pt).replace(".","p") + ".png")
        sleep(1)

def main():
    rootFiles, PU_types, infile_pt = Read_input(o.input)
    #print(rootFiles)
    for i, rootFile in enumerate(rootFiles):
        print("Making plots for ... " + str(infile_pt[i]))
        sleep(0.5)
        Eta_Phi(inFile=rootFile, infile_pt=infile_pt[i])

if __name__ == "__main__":
    main()
   

