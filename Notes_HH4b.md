

Code is should be fairly self explanatory. 

Big picture:
============

This code processes ntuples that are produced from CMSSW and makes a bunch of histograms.
   >  Actually takes two different inputs: one for HLT and one for Offline

   >  Matches events between the two using Run and Event numbers
     
  Outputs from CMSSW reco:  HLT and Offline are inputs to this code. 
     
 Almost all "analysis" that we do uses these histograms.
   >  eg: ratios of histograms in different folders can give efficiencies or fake rates as a function of any variable plotted

Overview of output ROOT file produced
=======================================

ROOT file has a bunch of directories which contain histograms

####  <ins>Directories:</ins>
  - `*HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1`             - trigger for QuadJet only
  - `*HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1` - trigger for QuadJet and Ht
  - `*HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1` - trigger for QuadJet, Ht, and TriplePFPuppi
  - `L1_untagged` - trigger L1 for untagged bjets
  - `tagged_noL1` - no triggers applied and jets are btagged
  - `tagged_L1` - trigger L1 applied and jets are btagged
  - `noCuts` - noCuts on triggers or tagging
  - `trigger` - useful for debugging

####  <ins>Some more details:</ins>
  - `tagged_*` - passed btagging with set DeepCSV score


#####  <ins>plots:</ins>

  - `deepCSV` - deepCSV score for applied conditions 
  - `mass` - invariant mass of event
  - `Ht` - Ht for applied conditions on events
  - `pt_i` - ith pt stored in event passing given conditions (1 thru 4 available only)
  - `pt_i_cut` - applied turn on conditions from pt_4 down to pt_1
  

