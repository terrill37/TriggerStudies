#include <iostream>
#include <iomanip>
#include <cstdio>
#include <TROOT.h>
#include <boost/bind.hpp>


#include "TriggerStudies/NtupleAna/interface/HH4bAnalysis.h"
#include "nTupleAnalysis/baseClasses/interface/helpers.h"

using std::cout; using std::endl;
using namespace TriggerStudies;
using std::vector;  using std::map; using std::string; using std::set;




HH4bAnalysis::HH4bAnalysis(TChain* _eventsRAW, TChain* _eventsAOD, fwlite::TFileService& fs, bool _debug){
  
  
  if(_debug) cout<<"In HH4bAnalysis constructor"<<endl;
  debug      = _debug;

  eventsRAW     = _eventsRAW;
  eventsRAW->SetBranchStatus("*", 0);

  eventsAOD     = _eventsAOD;
  eventsAOD->SetBranchStatus("*", 0);


  event      = new eventData(eventsRAW, eventsAOD, true, "2018", debug, "");
  treeEvents = eventsRAW->GetEntries();

  cutflow    = new nTupleAnalysis::cutflowHists("cutflow", fs);
  
  triggers     = new nTupleAnalysis::triggers("trigger", fs);
  mass_preCut  = new nTupleAnalysis::mass("L1_noDeepCSVCut", fs);
  mass_postCut = new nTupleAnalysis::mass("L1_wDeepCSVCut", fs);
  mass_trig1   = new nTupleAnalysis::mass("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);
  mass_trig2   = new nTupleAnalysis::mass("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1",fs);

  cutflow->AddCut("all");
  cutflow->AddCut("foundMatch");
  //cutflow->AddCut("passMuonCut");
  //cutflow->AddCut("passElecCut");
  //cutflow->AddCut("passLeptonCut");
  cutflow->AddCut("passNJetCut");
  cutflow->AddCut("passNBJetCut");



  //
  // hists
  //
  dir = fs.mkdir("HH4bAnalysis");

  hEvents                 = new nTupleAnalysis::eventHists("Events", fs);

  //h4b_all           = dir.make<TH1F>("mtt_off",            "BTagAnalysis/mtt_off;             mtt;   Entries", 100,-0.01, 2);

}


void HH4bAnalysis::monitor(long int e){
 
  
  //Monitor progress
  percent        = (e+1)*100/nEvents;
  duration       = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  eventRateAVE      = (e+1)/duration;
  timeRemaining  = (nEvents-e)/eventRateAVE;
  minutes = static_cast<int>(timeRemaining/60);
  seconds = static_cast<int>(timeRemaining - minutes*60);
  getrusage(who, &usage);
  usageMB = usage.ru_maxrss/1024;


  timeStep = (std::clock() - lastTime) / (double) CLOCKS_PER_SEC;
  eventStep = e - lastEvent;
  eventRate = eventStep/timeStep;
  //print status and flush stdout so that status bar only uses one line
  fprintf(stdout, "\rProcessed: %8li of %li ( %2li%% | %.0f events/s AVE | %.0f events/s  | done in %02i:%02i | memory usage: %li MB)       ",
	  e+1, nEvents, percent,   eventRateAVE,    eventRate, minutes, seconds,                usageMB);
  fflush(stdout);
  lastTime = std::clock();
  lastEvent = e+1;
}

int HH4bAnalysis::eventLoop(int maxEvents, int nSkipEvents){
  //Set Number of events to process. Take manual maxEvents if maxEvents is > 0 and less than the total number of events in the input files.
  nEvents = (maxEvents > 0 && maxEvents < treeEvents) ? maxEvents : treeEvents;

  cout << "\nProcess " << nEvents << " of " << treeEvents << " events.\n";

  start = std::clock();//2546000 //2546043
  lastTime = std::clock();
  lastEvent = 0;
  for(long int e = 0; e < nEvents; e++){

    if(e < nSkipEvents) continue;

    event->update(e);
    processEvent();
    if(debug) event->dump();

    //periodically update status
    //if( (e+1)%1 == 0 || e+1==nEvents || debug)
    if( (e+1)%10000 == 0 || e+1==nEvents || debug)
      monitor(e);

  }


  cout << endl;
  cout << "HH4bAnalysis::End of Event Loop" << endl;


  minutes = static_cast<int>(duration/60);
  seconds = static_cast<int>(duration - minutes*60);


  fprintf(stdout,"---------------------------\nProcessed in %02i:%02i", minutes, seconds);

  return 0;
}

int HH4bAnalysis::processEvent(){
  //constants
  float deepCSV_cut = 0.3;
  float eta_cut     = 4.0;
  float pt_cut      = 30 ;
  //float phi_cut = 3.14 ;
  
  int bsetList[2][32];
  for(int i =0; i<2 ; i++){
    for(int j =0; j<32; j++){
      bsetList[i][j] = std::bitset<32>(event->BitTrigger[i])[j];
    }
  }

  int triggerBit_L1[2]= {0,0};
  int triggerBit_1[2] = {1,19}; // {nBitTrigger, BitTrigger}
  int triggerBit_2[2] = {0,18}; // {nBitTrigger, BitTrigger}

  if(debug) cout << "processEvent start" << endl;

  cutflow->Fill("all", 1.0);

  if(event->run != event->runAOD)
    return 0;

  if(event->event != event->eventAOD)
    return 0;

  cutflow->Fill("foundMatch", 1.0);



  //
  //  Offline BTags
  //
  
  //check if L1 triggered  in event before continuing
  if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) return 0;

  if(debug) cout << "Count BTags " << endl;
  unsigned int nOffJetsForCut = 0;
  unsigned int nOffJetsTaggedForCut = 0;
 
  for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
    
    if(fabs(offJet->eta) > eta_cut) continue;
    if(offJet->pt       < pt_cut)       continue; // 40 ? 

    ++nOffJetsForCut;
    
    if(offJet->DeepCSV > deepCSV_cut) ++nOffJetsTaggedForCut; // FIXME OFFLINE BTAG CUT: BTagged pass

  }
  float eventWeight = 1.0;

  //  
  if(nOffJetsForCut < 3      ){
    if(debug) cout << "Fail NJet Cut" << endl;
    return 0;
  }

  else{
    for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
        mass_preCut -> Fill(offJet, eventWeight); //fill the deepCSV scores of the mass_preCut directory
        
        if(nOffJetsTaggedForCut >= 4){
            mass_postCut -> Fill(offJet, eventWeight); //fill the deepCSV scores of the mass_postCut directory
            }
        }
    }
  
  //cutflow->Fill("passNJetCut", eventWeight);
  if(debug) cout << "Pass NJet Cut " << endl;

  bool doOfflineBTagCut = true;
  if(doOfflineBTagCut){
  
    if(nOffJetsTaggedForCut < 3) {
      if(debug) cout << "Fail NBJet Cut" << endl;
      return 0;
    }
  //  cutflow->Fill("passNBJetCut", eventWeight);
  }

  //if(debug) cout << "Pass NBJet Cut " << endl;
  
  if(nOffJetsForCut >= 4){          // at least four jets
    if(nOffJetsTaggedForCut >=4){   // at least four btagged jets
        if(debug) cout<<"pass 4 tagged"<<endl;
        
        float mass_post=0;            //variable storage of mass
        TLorentzVector momentum_post; //four momentum of variable
        int index = 1;                //index to check only four pass
        
        float massTrig1=0;
        float massTrig2=0;
        TLorentzVector momentum_trig1;
        TLorentzVector momentum_trig2;

        for(const nTupleAnalysis::jetPtr& offJet : event->offJets){ //loop through jets in event
          if(offJet->DeepCSV <= deepCSV_cut) continue; //checks that passes DeepCSV cut
          
          momentum_post += offJet->p;
          if(debug) cout<<"index in mass loop"<<index;
          if(index >= 4) break; // break loop after fourth jet added

          if(debug) cout<<momentum_post(0)<<endl;
         
        }

        mass_post = momentum_post.M();
        mass_postCut->FillMass(mass_post);
        
      
        //int bitList[2][32] = {bset_0_array,bset_1_array};

        //for first trigger
        if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
          for(const nTupleAnalysis::jetPtr& offJet : event->offJets){ //loop through jets in event
            if(offJet->DeepCSV <= deepCSV_cut) continue; //checks that passes DeepCSV cut
            
            momentum_trig1 += offJet->p;
            if(debug) cout<<"index in mass loop"<<index;
            if(index >= 4) break; // break loop after fourth jet added

            if(debug) cout<<momentum_trig1(0)<<endl;
           
          }
          massTrig1 = momentum_trig1.M();
          mass_trig1->FillMass(massTrig1);
        }
        
        //for second trigger
        if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
          for(const nTupleAnalysis::jetPtr& offJet : event->offJets){ //loop through jets in event
            if(offJet->DeepCSV <= deepCSV_cut) continue; //checks that passes DeepCSV cut
            
            momentum_trig2 += offJet->p;
            if(debug) cout<<"index in mass loop"<<index;
            if(index >= 4) break; // break loop after fourth jet added

            if(debug) cout<<momentum_trig2(0)<<endl;
           
          }
          massTrig2 = momentum_trig2.M();
          mass_trig2->FillMass(massTrig2);
        }



    }
  }


  //
  // Fill All events
  //
  if(debug) cout << "Fill All Events " << endl;
  //hEvents->Fill(event->offPVs.size(),  0.0, eventWeight);

  // Calculate m4b
  //float m4b = 1.0;//xx
  // Fill m4b
  //h4b_all->Fill(m4b);
  if(debug){  
    std::bitset<32> bset(event->BitTrigger[0]);

    cout << " Processing event " << endl;
    cout << " \T BitTrigger[0] " << event->BitTrigger[0] << endl;
    cout << " bit value: "       << bset << endl;
    cout << " \T BitTrigger[1] " << event->BitTrigger[1] << endl;
  }
  //iterate through bit trigger list
  for(int i = 0; i < 2; i++){
    //convert decimal into binary
    std::bitset<32> bset(event->BitTrigger[i]);
    
    //iterate through list of binary values
    for(long unsigned int j = 0; j < 32; j++){
      if(bset[j]==1){
        triggers -> Fill(j + 32*i);
        //if(std::bitset<32> (event->BitTrigger[0])[0] == 1){
          //triggers_L1Cut -> Fill(j+32*i);
        //}
        }
    }
  }

  // L1 L1_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV_2p4_v1
  //if(passTriggerBit(0,4)){
  //   h4b_L1->Fill(m4b);    
  // HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1 # 1 19
  //  if(passTriggerBit(1,19)){
  //   h4b_HLT_noBTag->Fill(m4b);    
  // HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1 # 0 18
  //  if(passTriggerBit(1,19)){
  //   h4b_HLT_wBTag->Fill(m4b);    
  ///}
  //}


  return 0;
}

//passTrigBit(unsigned int nBit, unsigned int trigIndex 
//  
//  return bool(event->BitTrigger[nBit] & (1 << trigIndex));


HH4bAnalysis::~HH4bAnalysis(){
  cout << "HH4bAnalysis::destroyed" << endl;
}


