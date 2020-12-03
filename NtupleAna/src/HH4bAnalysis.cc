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
  
  mass_preCut  = new nTupleAnalysis::mass("noCuts", fs);
  
  L1_untagged = new nTupleAnalysis::mass("L1_untagged", fs);
  trig1   = new nTupleAnalysis::mass("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);
  trig2   = new nTupleAnalysis::mass("HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1",fs);
  trig3 = new nTupleAnalysis::mass("HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);

  deepCut_noL1 = new nTupleAnalysis::mass("tagged_noL1", fs);
  L1_deepCut_tagged = new nTupleAnalysis::mass("tagged_L1", fs);
  trig1_tagged   = new nTupleAnalysis::mass("tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);
  trig2_tagged   = new nTupleAnalysis::mass("tagged_HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV0p5_2p4_v1",fs);
    trig3_tagged = new nTupleAnalysis::mass("tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);


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
 // cout<<"event"<<endl;
  
  //constants
  float deepCSV_cut = 0.3;
  float eta_cut     = 4.0;
  float pt_cut      = 30 ;
  
  for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
    triggers -> Fillpt_all(offJet->pt);
    if(offJet->pt < pt_cut) continue;
    triggers -> Fillpt_cut(offJet->pt);
  }

  //initial pt for all events
  for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
    triggers -> Fillpt_initial(offJet->pt);
    break;
  }

  int bsetList[2][32];
  for(int i =0; i<2 ; i++){
    for(int j =0; j<32; j++){
      bsetList[i][j] = std::bitset<32>(event->BitTrigger[i])[j];
    }
  }

  int triggerBit_L1[2]= {0,4};
  int triggerBit_1[2] = {1,19}; // {nBitTrigger, BitTrigger}
  int triggerBit_2[2] = {0,18}; // {nBitTrigger, BitTrigger}
  int triggerBit_3[2] = {1,21};

  if(debug) cout << "processEvent start" << endl;

  cutflow->Fill("all", 1.0);

  if(event->run != event->runAOD)
    return 0;

  if(event->event != event->eventAOD)
    return 0;

  cutflow->Fill("foundMatch", 1.0);
  //cout<<"event matched"<<endl;


  //
  //  Offline BTags
  //
  
  //check if L1 triggered  in event before continuing
 // if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) return 0;

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
  
  //fill deepCSV plots
  else{
    for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
        if(fabs(offJet->eta) > eta_cut) continue;
        if(offJet->pt       < pt_cut)       continue; // 40 ? 

        mass_preCut -> Fill(offJet, eventWeight); //fill the deepCSV scores of the mass_preCut directory
        //four jet requirement
        if(nOffJetsTaggedForCut >= 4){
            // pass DeepCSV cut
            if(offJet->DeepCSV < deepCSV_cut) continue;
            deepCut_noL1 -> Fill(offJet,eventWeight);

            //Pass L1
            if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
            L1_deepCut_tagged->Fill(offJet, eventWeight);

            //first trigger
            if(bsetList[triggerBit_1[0]][triggerBit_1[1]] ==1){
              trig1_tagged -> Fill(offJet, eventWeight);
            }
            //second trigger
            if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
              trig2_tagged -> Fill(offJet, eventWeight);
            }
            //third trigger
            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
              trig3_tagged -> Fill(offJet, eventWeight);
            }
            
          }
      }
    }

    if(nOffJetsForCut >= 4){
        for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
            if(fabs(offJet->eta) > eta_cut) continue;
            if(offJet->pt       < pt_cut)       continue;

            if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
            L1_untagged -> Fill(offJet, eventWeight);

            if(bsetList[triggerBit_1[0]][triggerBit_1[1]] ==1){
                trig1 -> Fill(offJet, eventWeight);
            }

            if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
                trig2 -> Fill(offJet, eventWeight);
            }

            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
                trig3 -> Fill(offJet, eventWeight);
            }
        }
    }
  
  cutflow->Fill("passNJetCut", eventWeight);
  if(debug) cout << "Pass NJet Cut " << endl;

  bool doOfflineBTagCut = true;
  if(doOfflineBTagCut){
  
    if(nOffJetsTaggedForCut < 3) {
      if(debug) cout << "Fail NBJet Cut" << endl;
      return 0;
    }
    cutflow->Fill("passNBJetCut", eventWeight);
  }

  //no cuts on Trigger or DeepCSV
  float m_preCut=0;                 //variable storage of mass
  TLorentzVector momentum_preCut;   //four momentum of variable
  float Ht_preCut=0;                //Ht values
  std::vector<float> pt_preCut;

  //cuts only on DeepCSV
  float mass_deepCut_noL1=0;
  TLorentzVector momentum_deepCut_noL1;
  float Ht_deepCut_noL1=0;
  std::vector<float> pt_deepCut_noL1;
  
  //
  //tagged
  //

  //cuts on L1 and DeepCSV
  float mass_L1_deepCut_tagged=0;
  TLorentzVector momentum_L1_deepCut_tagged;
  float Ht_L1_deepCut_tagged=0;
  std::vector<float> pt_L1_deepCut_tagged;
  
  //cuts on L1, DeepCSV and trigger 1
  float mass_trig1_tagged = 0;
  TLorentzVector momentum_trig1_tagged;
  float Ht_trig1_tagged=0;
  std::vector<float> pt_trig1_tagged;
  
  //cuts on L1, DeepCSV, and trigger 2
  float mass_trig2_tagged = 0;
  TLorentzVector momentum_trig2_tagged;
  float Ht_trig2_tagged=0;
  std::vector<float> pt_trig2_tagged;
  
  //cuts on L1, DeepCSV, and trigger 3
  float mass_trig3_tagged = 0;
  TLorentzVector momentum_trig3_tagged;
  float Ht_trig3_tagged = 0;
  std::vector<float> pt_trig3_tagged;
  
  //
  //untagged
  //
  
  //cuts on L1
  float mass_L1_untagged = 0;
  TLorentzVector momentum_L1_untagged;
  float Ht_L1_untagged = 0;
  std::vector<float> pt_L1_untagged;
  
  //cuts on L1 and trigger 1
  float mass_trig1 = 0;
  TLorentzVector momentum_trig1;
  float Ht_trig1=0;
  std::vector<float> pt_trig1;
  
  //cuts on L1 and trigger 2
  float mass_trig2 = 0;
  TLorentzVector momentum_trig2;
  float Ht_trig2=0;
  std::vector<float> pt_trig2;
  
  //cuts on L1 and trigger 3
  float mass_trig3 = 0;
  TLorentzVector momentum_trig3;
  float Ht_trig3=0;
  std::vector<float> pt_trig3;


  if(nOffJetsForCut >= 4){ // at least four jets
  if(debug) cout<<"pass 4 jets"<<endl;
      
      int index = 1; //index to check only four pass; resets after each loop
      
      for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
          if(fabs(offJet->eta) > eta_cut) continue;
          if(offJet->pt       < pt_cut)       continue; // 40 ? 

          momentum_preCut += offJet->p;
          pt_preCut.push_back(offJet->pt);
          if(index >=4){
              index=1;
              m_preCut = momentum_preCut.M();
              mass_preCut ->FillMass(m_preCut);
              mass_preCut ->Fillpts(pt_preCut);
              break;
          }
          index++;
      } 

      //L1_noCut
      if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
          //L1_untagged
          for(const nTupleAnalysis::jetPtr& offJet : event->offJets){    
              if(fabs(offJet->eta) > eta_cut) continue;
              if(offJet->pt       < pt_cut)       continue; // 40 ? 

              momentum_L1_untagged += offJet->p;
              pt_L1_untagged.push_back(offJet->pt);
              if(index>=4){
                  index=1;
                  mass_L1_untagged = momentum_L1_untagged.M();
                  L1_untagged ->FillMass(mass_L1_untagged);
                  L1_untagged ->Fillpts(pt_L1_untagged);
                  break;
                }
                index++;
          }
          //trig1
          //cout<<"trigger mass bit: "<<bsetList[triggerBit_1[0]][triggerBit_1[1]]<<endl;
          if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  momentum_trig1 += offJet->p;
                  pt_trig1.push_back(offJet->pt);
                  if(index>=4){
                      index =1;
                      mass_trig1 = momentum_trig1.M();
                      trig1 -> FillMass(mass_trig1);
                      trig1 -> Fillpts(pt_trig1);
                      break;
                  }
                  index++;
              }
          }

          //trig2
          if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  momentum_trig2 += offJet->p;
                  pt_trig2.push_back(offJet->pt);
                  if(index>=4){
                      index=1;
                      mass_trig2 = momentum_trig2.M();
                      trig2 -> FillMass(mass_trig2);
                      trig2 -> Fillpts(pt_trig2);
                      break;
                  }
                  index++;
              }
          }
          
          //trig3
          if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  momentum_trig3 += offJet->p;
                  pt_trig3.push_back(offJet->pt);
                  if(index>=4){
                      index=1;
                      mass_trig3 = momentum_trig3.M();
                      trig3 -> FillMass(mass_trig3);
                      trig3 -> Fillpts(pt_trig3);
                      break;
                  }
                  index++;
              }
          }

      }
  
  }          // end of at least four jets
  
  if(nOffJetsTaggedForCut >=4){   // at least four btagged jets
      if(debug) cout<<"pass 4 tagged"<<endl;
      
      int index = 1; //index to check only four pass; resets after each loop 

      //deepCut_noL1
      for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
          if(fabs(offJet->eta) > eta_cut) continue;
          if(offJet->pt       < pt_cut)       continue; // 40 ? 

          if(offJet->DeepCSV < deepCSV_cut) continue;
          momentum_deepCut_noL1 += offJet->p;
          pt_deepCut_noL1.push_back(offJet->pt);
          if(index >=4){
              index=1;
              mass_deepCut_noL1 = momentum_deepCut_noL1.M();
              deepCut_noL1->FillMass(mass_deepCut_noL1);
              deepCut_noL1->Fillpts(pt_deepCut_noL1);
              break;
          }
          index++;
          
      }

      //L1_noCut
      if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
          //L1_deepCut
          for(const nTupleAnalysis::jetPtr& offJet : event->offJets){    
              if(fabs(offJet->eta) > eta_cut) continue;
              if(offJet->pt       < pt_cut)       continue; // 40 ? 

              if(offJet->DeepCSV < deepCSV_cut) continue;
              momentum_L1_deepCut_tagged += offJet->p;
              pt_L1_deepCut_tagged.push_back(offJet->pt);
              if(index>=4){
                  index=1;
                  mass_L1_deepCut_tagged = momentum_L1_deepCut_tagged.M();
                  L1_deepCut_tagged ->FillMass(mass_L1_deepCut_tagged);
                  L1_deepCut_tagged ->Fillpts(pt_L1_deepCut_tagged);
                  break;
                }
                index++;
          }
          
          //trig1
          //cout<<"trigger mass bit: "<<bsetList[triggerBit_1[0]][triggerBit_1[1]]<<endl;
          if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  if(offJet->DeepCSV < deepCSV_cut) continue;
                  momentum_trig1_tagged += offJet->p;
                  pt_trig1_tagged.push_back(offJet->pt);
                  if(index>=4){
                      index =1;
                      mass_trig1_tagged = momentum_trig1_tagged.M();
                      trig1_tagged -> FillMass(mass_trig1_tagged);
                      trig1_tagged -> Fillpts(pt_trig1_tagged);
                      break;
                  }
                  index++;
              }
          }

          //trig2
          if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  
                  if(offJet->DeepCSV < deepCSV_cut) continue;
                  momentum_trig2_tagged += offJet->p;
                  pt_trig2_tagged.push_back(offJet->pt);
                  if(index>=4){
                      index=1;
                      mass_trig2_tagged = momentum_trig2_tagged.M();
                      trig2_tagged -> FillMass(mass_trig2_tagged);
                      trig2_tagged -> Fillpts(pt_trig2_tagged);
                      break;
                  }
                  index++;
              }
          }
          
          //trig3
          if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
              for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
                  if(fabs(offJet->eta) > eta_cut) continue;
                  if(offJet->pt       < pt_cut)       continue; // 40 ? 

                  
                  if(offJet->DeepCSV < deepCSV_cut) continue;
                  momentum_trig3_tagged += offJet->p;
                  pt_trig3_tagged.push_back(offJet->pt);
                  if(index>=4){
                      index=1;
                      mass_trig3_tagged = momentum_trig3_tagged.M();
                      trig3_tagged -> FillMass(mass_trig3_tagged);
                      trig3_tagged -> Fillpts(pt_trig3_tagged);
                      break;
                  }
                  index++;
              }
          }

      }
  
  }
  
    
    //
    //
    //Ht filling section
    //doesn't require index
    //does require pt>ptCut
    //
    //
  
  //tagged filling
  if(nOffJetsTaggedForCut >= 4){
    for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
        if(fabs(offJet->eta) > eta_cut) continue;
        if(offJet->pt       < pt_cut)       continue; 
        
        //cut on DeepCSV
        if(offJet->DeepCSV < deepCSV_cut) continue;
        Ht_deepCut_noL1 += offJet -> pt;

        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
        Ht_L1_deepCut_tagged += offJet->pt;
        
        //trigger 1 cut
        if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
            Ht_trig1_tagged += offJet->pt;
        }

        //trigger 2 cut
        if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
            Ht_trig2_tagged += offJet->pt;
        }

        //trigger 3 cut
        if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
            Ht_trig3_tagged += offJet->pt;
        }

    }
  }
  
  //untagged filling
  if(nOffJetsForCut>=4){
    for(const nTupleAnalysis::jetPtr& offJet : event->offJets){
        //if(fabs(offJet->eta) > eta_cut) continue;
        //cut on pt
        if(offJet->pt       < pt_cut)       continue; 
        Ht_preCut += offJet->pt;                        //precuts should be moved to correspond with untagged events

        //cut on DeepCSV
        //if(offJet->DeepCSV < deepCSV_cut) continue;
        //Ht_deepCut_noL1 += offJet->pt;
        
        //Cut on L1
        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
        Ht_L1_untagged += offJet->pt;
        
        //trigger 1 cut
        if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
            Ht_trig1 += offJet->pt;
        }

        //trigger 2 cut
        if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
            Ht_trig2 += offJet->pt;
        }

        //trigger 3 cut
        if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
            Ht_trig3 += offJet->pt;
        }

    }
  }

    //Fill Ht's
    if(Ht_preCut != 0){
        mass_preCut -> FillHt(Ht_preCut);
        }
    if(Ht_deepCut_noL1 != 0){
        deepCut_noL1-> FillHt(Ht_deepCut_noL1);
        }
    if(Ht_L1_deepCut_tagged!=0){
        L1_deepCut_tagged  -> FillHt(Ht_L1_deepCut_tagged);
        }
    if(Ht_trig1_tagged!=0){
        trig1_tagged       -> FillHt(Ht_trig1_tagged);
    }
    if(Ht_trig2_tagged!=0){
        trig2_tagged       -> FillHt(Ht_trig2_tagged);
    }
    if(Ht_trig3_tagged!=0){
        trig3_tagged       -> FillHt(Ht_trig3_tagged);
    }
    if(Ht_L1_untagged!=0){
        L1_untagged  -> FillHt(Ht_L1_untagged);
        }
    if(Ht_trig1!=0){
        trig1       -> FillHt(Ht_trig1);
    }
    if(Ht_trig2_tagged!=0){
        trig2       -> FillHt(Ht_trig2);
    }
    if(Ht_trig3_tagged!=0){
        trig3       -> FillHt(Ht_trig3);
    }

   // }
  //}


  //
  // Fill All events
  //
  if(debug) cout << "Fill All Events " << endl;
  //hEvents->Fill(event->offPVs.size(),  0.0, eventWeight);

  if(debug){  
    std::bitset<32> bset(event->BitTrigger[0]);

    cout << " Processing event " << endl;
    cout << " BitTrigger[0] " << event->BitTrigger[0] << endl;
    cout << " bit value: "       << bset << endl;
    cout << " BitTrigger[1] " << event->BitTrigger[1] << endl;
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


