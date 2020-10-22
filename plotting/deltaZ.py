import ROOT
from time import sleep

from ROOT import TPad

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

def pair_PU200_noPU(pt_check, PU_check, inFile, infile_pt, PU_type):
    #output of second hist that matches from PU200
    for i, pt_match in enumerate(infile_pt):
        
        if (pt_match ==  pt_check and PU_type[i] != "PU200" and PU_check == "PU200"):
            print("Match Found")
            return inFile[i], PU_type[i]

        
    print("ERROR: Match Not Found for " + str(PU_type[i]) + " Pt " + str(pt_check))
    return "NONE","NONE"

def normalize(hist1, hist2):
    #normalize hists and return hists
    if hist1.Integral():
        hist1.Scale(1./hist1.Integral())

    if hist2.Integral():
        hist2.Scale(1./hist2.Integral())

    return hist1, hist2


def delta_z(inFile1, inFile2, infile_pt, PU_type1, PU_type2):
    print("Begin delta-z plotting...")
    print("Retrieving ttree ... ")
    direc = "offJets_matchedJet_B" #online jet
 
    #direc_2 = "offJets_matchedJet_B_bad"
    
    path = "tracks/"
    plots = ["dz_s", "dz_m", "dzSig_m", "dzSig_s", "dzError", "dz_l", "dzSig_l"]
    
    for i, dz in enumerate(plots):
        h_name = direc+"/"+path+dz
        #h_name2 = direc+"/"+path+dz

        hist1 = inFile1.Get(h_name)
        hist2 = inFile2.Get(h_name)

        can = ROOT.TCanvas("can", "canvas", 800, 800)
        can.cd()
        #can.SetLogy(1)

        hist1.SetLineColor(1)
        hist2.SetLineColor(2)
        hist1.SetMarkerColor(1)
        hist2.SetMarkerColor(2)
        
        hist1.SetFillColor(ROOT.kSpring)
        hist1.SetFillStyle(3003)
        
        hist1.SetMinimum(0.00000000001)
        
        #set Pad1
        pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(.01) #joins upper and lower plot
        pad1.Draw()
        
        can.cd() #return to main canvas before defining next pad

        #set Pad2 for ratio
        pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()

        
        #draw upper plots in Pad1
        pad1.cd()
        hist1.Draw("h")
        hist2.Draw("same p")
        hist1.SetMinimum(0.0000001)

        hist1.GetXaxis().SetTitle(plots[i] + " [cm]")
        hist2.GetYaxis().SetTitle("dz/dzError")
        
        #move to pad2 for ratio plot
        pad2.cd()
        
        print("Normalizing hists ...")
        sleep(0.5)
        hist1, hist2 = normalize(hist1, hist2)

        #set info for ratio
        hist3 = hist1.Clone("hist3")
        hist3.SetLineColor(4)
        hist3.SetMarkerStyle(21)
        
        
        
        hist3.Sumw2()
        hist3.SetStats(0)
        hist3.Divide(hist2)
        
        h3Max = hist3.GetMaximum()
        hist3.SetMaximum(h3Max+0.1*h3Max)
        hist3.GetXaxis().SetTitle(plots[i] + " [cm] ")
        
        #draw ratio plot hist3
        hist3.Draw("pe")
        
        maxy = hist1.GetMaximum()
        hist1.GetYaxis().SetRangeUser(0.00000000001, maxy + 0.1*maxy)
        hist1.SetMinimum(0.001)

        pad1.SetLogy(1)
        
        #return to pad1 for legend
        pad1.cd()
        legend = ROOT.TLegend(0.2,0.75,0.45,0.9)
        legend.AddEntry(hist1, PU_type1, 'L')
        legend.AddEntry(hist2, PU_type2, 'L')
        legend.SetLineWidth(0)
        legend.Draw("same")

        can.SaveAs(o.outDir+"/tracks_"+dz+"_pt"+str(infile_pt).replace(".","p")+o.nameTag+".png")
        can.SaveAs(o.outDir+"/tracks_"+dz+"_pt"+str(infile_pt).replace(".","p")+o.nameTag+".pdf")
        
        #close canvas before next iteration through loop
        can.Close()

def main():
    rootFiles, PU_types, infile_pt = Read_input(o.input)
    #print(rootFiles)
    for i, rootFile in enumerate(rootFiles):
        inFile2, PU_type2 = pair_PU200_noPU(pt_check=infile_pt[i], PU_check=PU_types[i], inFile=rootFiles, infile_pt=infile_pt, PU_type=PU_types)
        
        if inFile2 != "NONE":
            if PU_types[i] == "noPU":
                print("Ignoring noPU starting")
                break
            else:
                print("Making plots for ... " + str(infile_pt[i]))
                print(PU_types[i],PU_type2)
                delta_z(inFile1 = rootFile, inFile2 = inFile2, infile_pt=infile_pt[i],
                       PU_type1=PU_types[i], PU_type2=PU_type2)
        
        elif inFile2 =="NONE":
            print("Ignoring Error. Moving to next iteration")
            continue
        
        #print("Making plots for ... " + str(infile_pt[i]))
        #delta_z(inFile=rootFile, infile_pt=infile_pt[i])
        sleep(0.5)
        

if __name__ == "__main__":
    main()
 
