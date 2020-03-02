#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include "TSystem.h"
#include "TChain.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "TriggerStudies/NtupleAna/interface/EventDataOLD.h"
#include "TriggerStudies/NtupleAna/interface/JetDataHandler.h"
#include "TriggerStudies/NtupleAna/interface/LeptonDataHandler.h"
#include "TriggerStudies/NtupleAna/interface/EventDisplayData.h"


#include "TriggerStudies/NtupleAna/interface/TrackHists.h"
#include "TriggerStudies/NtupleAna/interface/JetHists.h"
#include "TriggerStudies/NtupleAna/interface/EventHists.h"
#include "TriggerStudies/NtupleAna/interface/Helpers.h"

using namespace NtupleAna;
using std::cout;
using std::endl;
using std::string;
using optutl::CommandLineParser;

int main(int argc, char * argv[]){
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();


  // initialize command line parser
  CommandLineParser parser ("Analyze FWLite Histograms");

  parser.addOption ("configFile",   CommandLineParser::kString,
		    "configFileName");
  parser.addOption ("puFile",   CommandLineParser::kString,
		    "File with pile-up correction weights");
  // set defaults
  parser.integerValue ("maxEvents"  ) = -1;
  parser.integerValue ("outputEvery") =   10000;
  parser.stringValue  ("outputFile" ) = "TEST.root";


  parser.parseArguments (argc, argv);
  //cout << "Running with "<<parser.integerValue("maxEvents") << " config file " << parser.stringValue("configFile") << endl;
  string configFile     = parser.stringValue("configFile");
  string outputFileName = parser.stringValue("outputFile");
  string puFile         = parser.stringValue("puFile");
  std::vector<string> inputFiles = parser.stringVector("inputFiles");

  cout << "\t Running with configFile: " << configFile << endl;
  cout << "\t Outputting to: " << outputFileName << endl;
  cout << "\t puFile: " << puFile << endl;


  if( !edm::readPSetsFrom(configFile)->existsAs<edm::ParameterSet>("process") ){
    cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << endl; exit(0);
  }

  //
  // Get the config
  //
  // now get each parameter
  const edm::ParameterSet& process = edm::readPSetsFrom(configFile)->getParameter<edm::ParameterSet>("process");
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("procNtupleExample");
  bool makeEventDisplays = ana.getParameter<bool>("MakeEventDisplays");
  bool loadTrkLevel = ana.getParameter<bool>("LoadTrkLevel");
  bool debug = ana.getParameter<bool>("debug");
  //string offJetName = ana.getParameter<string>("offJetName");


  //
  //  Inint Tree
  //
  TChain* tree = new TChain("tree");
  cout << "\t inputFiles are: " << endl;
  for(unsigned int iFile=0; iFile<inputFiles.size(); ++iFile){
    cout << "\t\t " << inputFiles[iFile].c_str() << endl;
    // open input file (can be located on castor)

    //TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    //if( inFile ){
    tree->Add(inputFiles[iFile].c_str());
    //}
  }

  tree->SetBranchStatus("*",0);
  TObjArray *branchList = tree->GetListOfBranches();

  int nBranch     = tree->GetNbranches();
  for(int i=0;i<nBranch;i++){
    //cout << "Turning off " << branchList->At(i)->GetName() << endl;
    tree->SetBranchStatus(branchList->At(i)->GetName(),0);
  }


  //
  // Input Data
  //
  EventDataOLD eventData = EventDataOLD();
  eventData.SetBranchAddress(tree);

  JetDataHandler pfJetsDB = JetDataHandler("pfJets",loadTrkLevel);
  pfJetsDB.SetBranchAddress(tree);

  JetDataHandler caloJetsDB = JetDataHandler("caloJets",loadTrkLevel);
  caloJetsDB.SetBranchAddress(tree);

  JetDataHandler offJetsDB = JetDataHandler("offJets",loadTrkLevel);
  offJetsDB.SetBranchAddress(tree);

  LeptonDataHandler muonDB = LeptonDataHandler("offTightMuons","2017");
  muonDB.SetBranchAddress(tree);

  LeptonDataHandler elecDB = LeptonDataHandler("offTightElectrons","2017");
  elecDB.SetBranchAddress(tree);

  //
  // Make output ntuple/Hists
  //
  fwlite::TFileService fs = fwlite::TFileService(outputFileName);

  TFileDirectory cutFlowDir = fs.mkdir("CutFlow");
  TH1F* hCutFlowHist = cutFlowDir.make<TH1F>("cutflow", "cutflow", 1, 1, 2);
  hCutFlowHist->SetCanExtend(TH1::kAllAxes);
  int cf_bin_all  = hCutFlowHist->GetXaxis()->FindBin("All");
  int cf_bin_trig = hCutFlowHist->GetXaxis()->FindBin("PassTrigger");
  EventHists eventHists     = EventHists("AllEvents", fs);

  TrackHists pfTrackHists            = TrackHists("pfTracks", fs);
  TrackHists pfTrackHists_unmatched  = TrackHists("pfTracks_unmatched", fs);
  TrackHists pfTrackHists_matched    = TrackHists("pfTracks_matched", fs);
  TrackHists offTrackHists_matched   = TrackHists("offTracks_matched", fs);
  TrackHists offTrackHists_unmatched = TrackHists("offTracks_unmatched", fs);
  TrackHists offTrackHists           = TrackHists("offTracks", fs);

  JetHists pfJetHistsPreOLap     = JetHists("pfJetsPreOLap",fs,true);
  JetHists pfJetHists            = JetHists("pfJets", fs);
  JetHists pfJetHists_matched    = JetHists("pfJets_matched" ,fs);
  JetHists pfJetHists_matchedB   = JetHists("pfJets_matchedB",fs);
  JetHists pfJetHists_matchedL   = JetHists("pfJets_matchedL",fs);

  JetHists caloJetHistsPreOLap     = JetHists("caloJetsPreOLap",fs, true);
  JetHists caloJetHists            = JetHists("caloJets",fs);
  JetHists caloJetHists_matched    = JetHists("caloJets_matched" ,fs);
  JetHists caloJetHists_matchedB   = JetHists("caloJets_matchedB",fs);
  JetHists caloJetHists_matchedL   = JetHists("caloJets_matchedL",fs);

  JetHists offJetHistsPreOLap = JetHists("offJetsPreOLap",fs, true);

  JetHists offJetHists   = JetHists("offJets",  fs);
  JetHists offJetHists_B = JetHists("offJets_B",fs);
  JetHists offJetHists_L = JetHists("offJets_L",fs);

  JetHists offJetHists_matched        = JetHists("offJets_matched",  fs);
  JetHists offJetHists_matchedJet     = JetHists("offJets_matchedJet",  fs);
  JetHists offJetHists_B_matched      = JetHists("offJets_B_matched",fs);
  JetHists offJetHists_L_matched      = JetHists("offJets_L_matched",fs);

  JetHists offJetHists_matched_online60   = JetHists("offJets_matched_online60",  fs);
  JetHists offJetHists_B_matched_online60 = JetHists("offJets_B_matched_online60",  fs);
  JetHists offJetHists_L_matched_online60 = JetHists("offJets_L_matched_online60",  fs);

  JetHists offJetHists_matched_online70   = JetHists("offJets_matched_online70",  fs);
  JetHists offJetHists_B_matched_online70 = JetHists("offJets_B_matched_online70",  fs);
  JetHists offJetHists_L_matched_online70 = JetHists("offJets_L_matched_online70",  fs);

  JetHists offJetHists_matched_online80   = JetHists("offJets_matched_online80",  fs);
  JetHists offJetHists_B_matched_online80 = JetHists("offJets_B_matched_online80",  fs);
  JetHists offJetHists_L_matched_online80 = JetHists("offJets_L_matched_online80",  fs);

  JetHists offJetHists_matched_online90   = JetHists("offJets_matched_online90",  fs);
  JetHists offJetHists_B_matched_online90 = JetHists("offJets_B_matched_online90",  fs);
  JetHists offJetHists_L_matched_online90 = JetHists("offJets_L_matched_online90",  fs);

  JetHists offJetHists_offline70   = JetHists("offJets_offline70",  fs);
  JetHists offJetHists_offline70_B = JetHists("offJets_offline70_B",fs);
  JetHists offJetHists_offline70_L = JetHists("offJets_offline70_L",fs);

  JetHists offJetHists_offline70_matched        = JetHists("offJets_offline70_matched",  fs);
  JetHists offJetHists_offline70_B_matched      = JetHists("offJets_offline70_B_matched",fs);
  JetHists offJetHists_offline70_L_matched      = JetHists("offJets_offline70_L_matched",fs);

  JetHists offJetHists_offline70_matched_online60   = JetHists("offJets_offline70_matched_online60",  fs);
  JetHists offJetHists_offline70_B_matched_online60 = JetHists("offJets_offline70_B_matched_online60",  fs);
  JetHists offJetHists_offline70_L_matched_online60 = JetHists("offJets_offline70_L_matched_online60",  fs);

  JetHists offJetHists_offline70_matched_online70   = JetHists("offJets_offline70_matched_online70",  fs);
  JetHists offJetHists_offline70_B_matched_online70 = JetHists("offJets_offline70_B_matched_online70",  fs);
  JetHists offJetHists_offline70_L_matched_online70 = JetHists("offJets_offline70_L_matched_online70",  fs);

  JetHists offJetHists_offline70_matched_online80   = JetHists("offJets_offline70_matched_online80",  fs);
  JetHists offJetHists_offline70_B_matched_online80 = JetHists("offJets_offline70_B_matched_online80",  fs);
  JetHists offJetHists_offline70_L_matched_online80 = JetHists("offJets_offline70_L_matched_online80",  fs);

  JetHists offJetHists_offline70_matched_online90   = JetHists("offJets_offline70_matched_online90",  fs);
  JetHists offJetHists_offline70_B_matched_online90 = JetHists("offJets_offline70_B_matched_online90",  fs);
  JetHists offJetHists_offline70_L_matched_online90 = JetHists("offJets_offline70_L_matched_online90",  fs);

  //
  //  Event Displays
  //
  EventDisplayData eventDisplay = EventDisplayData("offline");


  cout << " In procNtupleExample " << endl;

  int nEventThisFile = tree->GetEntries();
  int maxEvents = parser.integerValue("maxEvents");
  int outputEvery = parser.integerValue("outputEvery");
  cout <<  "Number of input events: " << nEventThisFile << endl;

  for(int entry = 0; entry<nEventThisFile; entry++){

    if(entry %outputEvery == 0 || debug)
      cout << "Processed .... "<<entry<<" Events"<<endl;
    if( (maxEvents > 0) && (entry > maxEvents))
      break;

    hCutFlowHist ->Fill( cf_bin_all, 1 );
    tree->GetEntry( entry );

    eventData.SetEvent();

    if(debug){
      cout << "RunNumber: "  << eventData.runNumber << endl;
      cout << "EventNumber: " << eventData.eventNumber << endl;
    }

    //
    // Fill All events
    //
    eventHists.Fill(eventData);
    hCutFlowHist ->Fill( cf_bin_trig, 1 );

    // Converting from "row-level" info to "column-level" info
    std::vector<NtupleAna::LeptonData> elecs  = elecDB.GetLeps(30);
    std::vector<NtupleAna::LeptonData> muons  = muonDB.GetLeps(20);
    std::vector<NtupleAna::JetData>    offJets = offJetsDB.GetJets();
    std::vector<NtupleAna::JetData>    pfJets   = pfJetsDB.GetJets();
    std::vector<NtupleAna::JetData>    caloJets = caloJetsDB.GetJets();

    if((elecs.size()+muons.size()) < 2)  continue;

    for(JetData& offJet : offJets){

      if(fabs(offJet.m_eta) > 2.4) continue;
      if(offJet.m_pt       < 35)   continue;

      offJetHistsPreOLap.Fill(offJet);

      if(NtupleAna::failOverlap(offJet,elecs)) continue;
      if(NtupleAna::failOverlap(offJet,muons)) continue;

      if(makeEventDisplays) eventDisplay.AddJet("offJetAll", offJet, true);

      // Match offline to online
      float dR = 1e6;
      JetData* matchedJet = nullptr;
      for(JetData& pfJet : pfJets){
	float this_dR = pfJet.m_vec.DeltaR(offJet.m_vec);
	if (this_dR < dR){
	  dR = this_dR;
	  matchedJet = &pfJet;
	}
      }

      if( dR < 0.4){
	matchedJet->m_matchedJet = &offJet;
	matchedJet->m_match_dR   = dR;
	offJet.m_matchedJet = matchedJet;
	offJet.m_match_dR   = dR;
      }

      // match tracks if we matched jets
      if(offJet.m_matchedJet){

	if(makeEventDisplays){
	  eventDisplay.AddJet("offJet", offJet);
	  eventDisplay.AddJet("offMatchJet", *offJet.m_matchedJet);
	}

	for(TrackData& offTrack: offJet.m_tracks){
	  //need to check that the track (with matching resolution cone r=0.01) is in region where R=0.3 circles inside the two jets overlap!
	  //if offTrack.dR > 0.29 - offJet.match_dR: continue
	  if(offTrack.m_dR                                > 0.29) continue; // offTrack is not in cone of offJet
	  if(offTrack.m_vec.DeltaR(offJet.m_matchedJet->m_vec) > 0.29) continue; // offTrack is not in cone of pfJet
	  offTrackHists.Fill(offTrack);

	  if(makeEventDisplays) eventDisplay.AddTrk("offJet_Trks", offTrack);

	  float dR = 1e6;
	  float dR2 = 1e6;
	  TrackData* matchedTrack  = nullptr;
	  TrackData* secondClosest = nullptr;

	  for(TrackData& pfTrack: offJet.m_matchedJet->m_tracks){
	    //need to check that the track (with matching resolution cone r=0.01) is in region where R=0.3 circles inside the two jets overlap!
	    if(pfTrack.m_dR                     > 0.29) continue; // pfTrack is not in cone of pfJet
	    if(pfTrack.m_vec.DeltaR(offJet.m_vec) > 0.29) continue; // pfTrack is not in cone of offJet

	    //this_dR = ((offTrack.eta - pfTrack.eta)**2 + (offTrack.dPhi(pfTrack))**2)**0.5
	    float this_dR = offTrack.m_vec.DeltaR(pfTrack.m_vec);
	    if(this_dR > dR && this_dR < dR2){
	      dR2 = this_dR;
	      secondClosest = &pfTrack;
	    }

	    if(this_dR < dR){
	      dR2 = dR;
	      secondClosest = matchedTrack;

	      dR  = this_dR;
	      matchedTrack = &pfTrack;
	    }
	  }

	  if( dR > 0.01){
	    //if dR > 1e5:
	    offTrackHists_unmatched.Fill(offTrack);
	    if(makeEventDisplays) eventDisplay.AddTrk("offJet_TrksNoMatch", offTrack);
	    continue;
	  }

	  if(makeEventDisplays) eventDisplay.AddTrk("offJet_TrksMatch", offTrack);
	  matchedTrack->m_matchedTrack = &offTrack;
	  offTrack.m_matchedTrack     = matchedTrack;
	  offTrack.m_secondClosest    = secondClosest;
	  offTrack.m_nMatches              += 1;
	  offTrack.m_matchedTrack->m_nMatches += 1;
	  offTrackHists_matched.Fill(offTrack);
	  pfTrackHists_matched .Fill(*offTrack.m_matchedTrack);
	}

	for(const TrackData& pfTrack: offJet.m_matchedJet->m_tracks){
	  //need to check that the track (with matching resolution cone r=0.01) is in region where R=0.3 circles inside the two jets overlap!
	  if(pfTrack.m_dR                     > 0.29) continue; // pfTrack is not in cone of pfJet
	  if(pfTrack.m_vec.DeltaR(offJet.m_vec) > 0.29) continue; // pfTrack is not in cone of offJet

	  pfTrackHists.Fill(pfTrack); //all pftracks in matched jets
	  pfTrackHists.FillMatchStats(pfTrack); //check how often we match pfTracks to more than one offTrack
	  if(!pfTrack.m_nMatches){
	    pfTrackHists_unmatched.Fill(pfTrack); //all unmatched pftracks
	    pfTrackHists_unmatched.FillMatchStats(pfTrack);
	  }else{
	    pfTrackHists_matched.FillMatchStats(pfTrack);
	  }
	}
      }

      // Fill offJetHists
      // offJets_ROC
      //  wp60 DeepCSV > 0.76 (actual efficiency = 0.604992842829)
      //  wp70 DeepCSV > 0.56 (actual efficiency = 0.717503578586)
      //  wp80 DeepCSV > 0.36 (actual efficiency = 0.808474091039)
      //  wp90 DeepCSV > 0.12 (actual efficiency = 0.912533638706)
      offJetHists.Fill(offJet);
      if(offJet.m_deepcsv > 0.56) offJetHists_offline70.Fill(offJet);
      if(offJet.m_hadronFlavour == 5){
	offJetHists_B.Fill(offJet);
	if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B.Fill(offJet);
      }else{
	offJetHists_L.Fill(offJet);
	if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L.Fill(offJet);
      }

      if(offJet.m_matchedJet){
	offJetHists_matched.Fill(offJet);
	offJetHists_matchedJet.Fill(*offJet.m_matchedJet);
	if( offJet.m_deepcsv > 0.56) offJetHists_offline70_matched.Fill(offJet);

	if(offJet.m_hadronFlavour == 5){
	  offJetHists_B_matched.Fill(offJet);
	  if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B_matched.Fill(offJet);
	}else{
	  offJetHists_L_matched.Fill(offJet);
	  if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L_matched.Fill(offJet);
	}

        // pfJets_matched_ROC
	//  wp60 DeepCSV > 0.64 (actual efficiency = 0.610276798066)
	//  wp70 DeepCSV > 0.48 (actual efficiency = 0.708732661175)
	//  wp80 DeepCSV > 0.28 (actual efficiency = 0.814155211306)
	//  wp90 DeepCSV > 0.08 (actual efficiency = 0.924525480128)
	if(offJet.m_matchedJet->m_deepcsv > 0.64){ //approximate 60% Online WP
	  offJetHists_matched_online60.Fill(offJet);
	  if( offJet.m_deepcsv > 0.56) offJetHists_offline70_matched_online60.Fill(offJet);
	  if(offJet.m_hadronFlavour == 5){
	    offJetHists_B_matched_online60.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B_matched_online60.Fill(offJet);
	  }else{
	    offJetHists_L_matched_online60.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L_matched_online60.Fill(offJet);
	  }
	}

	if(offJet.m_matchedJet->m_deepcsv > 0.48){ //approximate 70% Online WP
	  offJetHists_matched_online70.Fill(offJet);
	  if(offJet.m_deepcsv > 0.56) offJetHists_offline70_matched_online70.Fill(offJet);
	  if(offJet.m_hadronFlavour == 5){
	    offJetHists_B_matched_online70.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B_matched_online70.Fill(offJet);
	  }else{
	    offJetHists_L_matched_online70.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L_matched_online70.Fill(offJet);
	  }
	}

	if(offJet.m_matchedJet->m_deepcsv > 0.28){ //approximate 80% Online WP
	  offJetHists_matched_online80.Fill(offJet);
	  if(offJet.m_deepcsv > 0.56) offJetHists_offline70_matched_online80.Fill(offJet);
	  if(offJet.m_hadronFlavour == 5){
	    offJetHists_B_matched_online80.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B_matched_online80.Fill(offJet);
	  }else{
	    offJetHists_L_matched_online80.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L_matched_online80.Fill(offJet);
	  }
	}

	if(offJet.m_matchedJet->m_deepcsv > 0.08){ //approximate 90% Online WP
	  offJetHists_matched_online90.Fill(offJet);
	  if(offJet.m_deepcsv > 0.56) offJetHists_offline70_matched_online90.Fill(offJet);
	  if(offJet.m_hadronFlavour == 5){
	    offJetHists_B_matched_online90.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_B_matched_online90.Fill(offJet);
	  }else{
	    offJetHists_L_matched_online90.Fill(offJet);
	    if(offJet.m_deepcsv > 0.56) offJetHists_offline70_L_matched_online90.Fill(offJet);
	  }
	}
      }//m_matchedJet

      // Match offline to online
      for(JetData& caloJet : caloJets){
	float deltaR = caloJet.m_vec.DeltaR(offJet.m_vec);
	if( deltaR < 0.4){
	  caloJet.m_matchedJet = &offJet;
	  break;
	}
      }

    }//offJets

    //
    //  pf Jets
    //
    for(JetData& pfJet : pfJets){
      if(fabs(pfJet.m_eta) > 4.5) continue;
      if(pfJet.m_pt       < 30)   continue;


      pfJetHistsPreOLap.Fill(pfJet);

      if(NtupleAna::failOverlap(pfJet,elecs)) continue;
      if(NtupleAna::failOverlap(pfJet,muons)) continue;

      if(makeEventDisplays) eventDisplay.AddJet("pfJet", pfJet, true);

      pfJetHists.Fill(pfJet);

      if(pfJet.m_matchedJet){
	pfJetHists_matched .Fill(pfJet);

	if( pfJet.m_matchedJet->m_hadronFlavour == 5)
	  pfJetHists_matchedB.Fill(pfJet);
	else
	  pfJetHists_matchedL.Fill(pfJet);
      }
    }

    //
    //  calo Jets
    //
    for(JetData& caloJet : caloJets){
      if(fabs(caloJet.m_eta) > 2.5) continue;
      if(caloJet.m_pt       < 35)   continue;

      caloJetHistsPreOLap.Fill(caloJet);
      if(NtupleAna::failOverlap(caloJet,elecs)) continue;
      if(NtupleAna::failOverlap(caloJet,muons)) continue;

      caloJetHists.Fill(caloJet);

      if(caloJet.m_matchedJet){
	caloJetHists_matched .Fill(caloJet);

	if(caloJet.m_matchedJet->m_hadronFlavour == 5)
	  caloJetHists_matchedB.Fill(caloJet);
	else
	  caloJetHists_matchedL.Fill(caloJet);
      }
    }

    if(makeEventDisplays) eventDisplay.NewEvent();


  }


  if(makeEventDisplays)
    eventDisplay.Write();


  return 0;
}
