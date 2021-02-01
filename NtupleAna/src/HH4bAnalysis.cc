#include <iostream>
#include <iomanip>
#include <cstdio>
#include <TROOT.h>
#include <boost/bind.hpp>


#include "TriggerStudies/NtupleAna/interface/HH4bAnalysis.h"
#include "nTupleAnalysis/baseClasses/interface/helpers.h"
#include "TriggerStudies/NtupleAna/interface/HH4bAnalysis_struct.h"

using std::cout; using std::endl;
using namespace TriggerStudies;
using std::vector;  using std::map; using std::string; using std::set;




HH4bAnalysis::HH4bAnalysis(TChain* _eventsRAW, fwlite::TFileService& fs, bool _debug, std::string jetDetailString){
  
  
  if(_debug) cout<<"In HH4bAnalysis constructor"<<endl;
  debug      = _debug;

  eventsRAW     = _eventsRAW;
  eventsRAW->SetBranchStatus("*", 0);

  event      = new eventData(eventsRAW, nullptr, true, "2018", debug, jetDetailString);
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
  trig3_tagged   = new nTupleAnalysis::mass("tagged_HLT_QuadPFPuppiJet_75_60_45_40_2p4_v1",fs);

  jetCount = new nTupleAnalysis::events("Jet_Counts",fs);
  
  triggers->AddTrig("All");
  triggers->AddTrig("L1");
  triggers->AddTrig("Quad");
  triggers->AddTrig("Quad+HT");
  triggers->AddTrig("Quad+HT+TriplePuppi");

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
  //float deepCSV_cut = 0.3;
  int flavour_b = 5;
  int flavour_c = 4;
  float eta_cut     = 2.4;
  float pt_cut      = 30 ;
  
  unsigned int minJet = 4;
  unsigned int minBJet = 4;


  int bsetList[2][32];
  for(int i =0; i<2 ; i++){
    for(int j =0; j<32; j++){
      bsetList[i][j] = std::bitset<32>(event->BitTrigger[i])[j];
    }
  }
  //will need to change 
  int triggerBit_L1[2]= {0,4};
  int triggerBit_1[2] = {1,19}; // {nBitTrigger, BitTrigger}, Quad + Ht
  int triggerBit_2[2] = {0,18}; // {nBitTrigger, BitTrigger}, Quad + Ht + triplePuppi
  int triggerBit_3[2] = {1,21}; // Quad only


  //
  //For Truth particle testing
  //
    
  //cout<<"before GenJet Loop"<<endl;
  
  std::vector<int> jetCounter (9) ;
  for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
    jetCounter.at(0)++;
    
    if(puppiJet->genJet_p.Pt()       < pt_cut)       continue;
    jetCounter.at(1)++;

    if(puppiJet->genJet_p.Eta()      > eta_cut)      continue;
    jetCounter.at(2)++;

    if(puppiJet->flavour == flavour_b){
        jetCounter.at(3)++;

        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
          jetCounter.at(4)++;
        } 

        if(bsetList[triggerBit_3[0]][triggerBit_3[1]] == 1){
          jetCounter.at(5)++;
        }

        if(bsetList[triggerBit_1[0]][triggerBit_1[1]] == 1){
          jetCounter.at(6)++;
        }

        if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
          jetCounter.at(7)++;
        }
    }
  }
  jetCount -> FillAll(jetCounter);
  
  if(debug){
    cout<<"Jet Count: "<<jetCounter.at(0)<<endl;
  }

  if(debug){
    for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
      cout<<"genJet"<<puppiJet->genJet_p.Pt() << "/" <<puppiJet->genJet_p.Eta() << "/" <<puppiJet->phi<<endl;
      cout<<"in genJet loop"<<endl;
    }
  }

  for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
    triggers -> Fillpt_all(puppiJet->genJet_p.Pt());
    if(puppiJet->genJet_p.Pt() < pt_cut) continue;
    triggers -> Fillpt_cut(puppiJet->genJet_p.Pt());
  }

  //initial pt for all events
  for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
    triggers -> Fillpt_initial(puppiJet->genJet_p.Pt());
    break;
  }

  
  if(debug) cout << "processEvent start" << endl;

  triggers->Fill_trigCount("All");
  if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
    triggers->Fill_trigCount("L1");
  }
  if(bsetList[triggerBit_3[0]][triggerBit_3[1]] == 1){
    triggers->Fill_trigCount("Quad");
  }
  if(bsetList[triggerBit_1[0]][triggerBit_1[1]] == 1){
    triggers->Fill_trigCount("Quad+HT");
  }
  if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
    triggers->Fill_trigCount("Quad+HT+TriplePuppi");
  }


  cutflow->Fill("all", 1.0);
  //cout<<"event matched"<<endl;


  //
  //  BTags
  //

  if(debug) cout << "Count BTags " << endl;
  unsigned int nTruthForCut = 0;
  unsigned int nTruthTaggedForCut = 0;
 
  for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
    
    if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
    if(puppiJet->genJet_p.Pt()       < pt_cut)       continue; // 40 ? 

    ++nTruthForCut;
    
    if(puppiJet->flavour == flavour_b) ++nTruthTaggedForCut; // ONLINE BTAG CUT: BTagged pass

  }
  
  if (debug) cout << "nTruthTaggedForCut: " << nTruthTaggedForCut <<endl;
  
  float eventWeight = 1.0;
  
  
 
  
  //fill flavour plots
  if(nTruthForCut >= minJet){
    for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
        if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
        if(puppiJet->genJet_p.Pt()       < pt_cut)       continue; // 40 ? 

        mass_preCut -> Fill(puppiJet, eventWeight); //fill the flavour vals of the mass_preCut
        
        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
            L1_untagged -> Fill(puppiJet, eventWeight);

            if(bsetList[triggerBit_1[0]][triggerBit_1[1]] ==1){
                trig1 -> Fill(puppiJet, eventWeight);
            }

            if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
                trig2 -> Fill(puppiJet, eventWeight);
            }

            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
                trig3 -> Fill(puppiJet, eventWeight);
            }
        }

        //four jet requirement that pass flavour b
        if(nTruthTaggedForCut >= minBJet){
            // pass flavour cut
            if(puppiJet->flavour != flavour_b) continue;
            deepCut_noL1 -> Fill(puppiJet,eventWeight);

            //Pass L1
            if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
            L1_deepCut_tagged->Fill(puppiJet, eventWeight);

            //first trigger
            if(bsetList[triggerBit_1[0]][triggerBit_1[1]] ==1){
              trig1_tagged -> Fill(puppiJet, eventWeight);
            }
            //second trigger
            if(bsetList[triggerBit_2[0]][triggerBit_2[1]] == 1){
              trig2_tagged -> Fill(puppiJet, eventWeight);
            }
            //third trigger
            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
              trig3_tagged -> Fill(puppiJet, eventWeight);
            }
            
          }
      }
    }
  
  //cutflow->Fill("passNJetCut", eventWeight);
  //if(debug) cout << "Pass NJet Cut " << endl;

  //no cuts on Trigger or flavour
  HH4bStruct preCut; 
  //cuts only on flavour
  HH4bStruct noL1_deepCut;
    
  //
  //tagged
  //

  //cuts on L1 and flavour
  HH4bStruct tagged_L1_deepCut; 
  //cuts on L1, flavour and trigger 1
  HH4bStruct tagged_trig1;
  //cuts on L1, flavour, and trigger 2
  HH4bStruct tagged_trig2;
  //cuts on L1, flavour, and trigger 3
  HH4bStruct tagged_trig3;
  
  //
  //untagged
  //
  
  //cuts on L1
  HH4bStruct untagged_L1;
  //cuts on L1 and trigger 1
  HH4bStruct untagged_trig1;
  //cuts on L1 and trigger 2
  HH4bStruct untagged_trig2;
  //cuts on L1 and trigger 3
  HH4bStruct untagged_trig3;

  if(nTruthForCut >= minJet){ // at least four jets
      if(debug) cout<<"pass 4 jets"<<endl;
      
     
      //cout<<"before mass loops"<<endl;
      for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
          if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
          if(puppiJet->genJet_p.Pt()       < pt_cut)       continue;
          
          preCut.momentum += puppiJet->p;
          preCut.pt.push_back(puppiJet->genJet_p.Pt());
          if(preCut.index == (int)minJet){
              preCut.mass = preCut.momentum.M(); 
              mass_preCut ->FillMass(preCut.mass);
              mass_preCut ->Fillpts(preCut.pt);
          }
          preCut.index ++;
          
          //L1
          if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
            untagged_L1.momentum += puppiJet->p;
            untagged_L1.pt.push_back(puppiJet->genJet_p.Pt());
            if(untagged_L1.index == (int)minJet){
                untagged_L1.mass = untagged_L1.momentum.M();
                L1_untagged ->FillMass(untagged_L1.mass);
                L1_untagged ->Fillpts(untagged_L1.pt);
            }
            untagged_L1.index++;
            
            //trig1
            if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
                untagged_trig1.momentum += puppiJet->p;
                untagged_trig1.pt.push_back(puppiJet->genJet_p.Pt());
                if(untagged_trig1.index == (int)minJet){
                    untagged_trig1.mass = untagged_trig1.momentum.M();
                    trig1 -> FillMass(untagged_trig1.mass);
                    trig1 -> Fillpts(untagged_trig1.pt);
                }
                untagged_trig1.index++;
            }

            //trig2
            if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
                untagged_trig2.momentum += puppiJet->p;
                untagged_trig2.pt.push_back(puppiJet->genJet_p.Pt());
                if(untagged_trig2.index == (int)minJet){
                    untagged_trig2.mass = untagged_trig2.momentum.M();
                    trig2 -> FillMass(untagged_trig2.mass);
                    trig2 -> Fillpts(untagged_trig2.pt);
                }
                untagged_trig2.index++;
            }
            
            //trig3
            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
                untagged_trig3.momentum += puppiJet->p;
                untagged_trig3.pt.push_back(puppiJet->genJet_p.Pt());
                if(untagged_trig3.index == (int)minJet){
                    untagged_trig3.mass = untagged_trig3.momentum.M();
                    trig3 -> FillMass(untagged_trig3.mass);
                    trig3 -> Fillpts(untagged_trig3.pt);
                }
                untagged_trig3.index++;
            }


          }
      } 
  }// end of at least four jets
  
  if(nTruthTaggedForCut >= minBJet){   // at least four btagged jets
      if(debug) cout<<"pass 4 tagged"<<endl; 

      //deepCut_noL1
      for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
          if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
          if(puppiJet->genJet_p.Pt()       < pt_cut)       continue; // 40 ? 

          if(puppiJet->flavour != flavour_b) continue; //check if b
          noL1_deepCut.momentum += puppiJet->p;
          noL1_deepCut.pt.push_back(puppiJet->genJet_p.Pt());
          if(noL1_deepCut.index == (int)minBJet){
              noL1_deepCut.mass = noL1_deepCut.momentum.M();
              deepCut_noL1->FillMass(noL1_deepCut.mass);
              deepCut_noL1->Fillpts(noL1_deepCut.pt);
          }
          noL1_deepCut.index++;
          
          //L1
          if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] == 1){
            tagged_L1_deepCut.momentum += puppiJet->p;
            tagged_L1_deepCut.pt.push_back(puppiJet->genJet_p.Pt());
            if(tagged_L1_deepCut.index == (int)minBJet){
                tagged_L1_deepCut.mass = tagged_L1_deepCut.momentum.M();
                L1_deepCut_tagged ->FillMass(tagged_L1_deepCut.mass);
                L1_deepCut_tagged ->Fillpts(tagged_L1_deepCut.pt);
            }
            tagged_L1_deepCut.index++;
          
            //trig1
            if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
              tagged_trig1.momentum += puppiJet->p;
              tagged_trig1.pt.push_back(puppiJet->genJet_p.Pt());
              if(tagged_trig1.index == (int)minBJet){
                  tagged_trig1.mass = tagged_trig1.momentum.M();
                  trig1_tagged-> FillMass(tagged_trig1.mass);
                  trig1_tagged-> Fillpts(tagged_trig1.pt);
              }
              tagged_trig1.index++;
            }
            
            //trig2
            if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
              tagged_trig2.momentum += puppiJet->p;
              tagged_trig2.pt.push_back(puppiJet->genJet_p.Pt());
              if(tagged_trig2.index == (int)minBJet){
                  tagged_trig2.mass = tagged_trig2.momentum.M();
                  trig2_tagged -> FillMass(tagged_trig2.mass);
                  trig2_tagged -> Fillpts(tagged_trig2.pt);
              }
              tagged_trig2.index++;
            }

            //trig3
            if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
                tagged_trig3.momentum += puppiJet->p;
                tagged_trig3.pt.push_back(puppiJet->genJet_p.Pt());
                if(tagged_trig3.index == (int)minBJet){
                    tagged_trig3.mass = tagged_trig3.momentum.M();
                    trig3_tagged -> FillMass(tagged_trig3.mass);
                    trig3_tagged -> Fillpts(tagged_trig3.pt);
                }
                tagged_trig3.index++;
            }
          }
      } 
  } //end of at least 4 b's 
    
    //
    //
    //Ht filling section
    //doesn't require index
    //does require pt>ptCut
    //
    //
  
  //tagged filling
  if(nTruthTaggedForCut >= minBJet){
    for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
        if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
        if(puppiJet->genJet_p.Pt()       < pt_cut)       continue; 
        
        //cut on flavour
        if(puppiJet->flavour != flavour_b) continue;
        noL1_deepCut.Ht += puppiJet -> genJet_p.Pt();

        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
        tagged_L1_deepCut.Ht += puppiJet->genJet_p.Pt();
        
        //trigger 1 cut
        if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
            tagged_trig1.Ht += puppiJet->genJet_p.Pt();
        }

        //trigger 2 cut
        if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
            tagged_trig2.Ht += puppiJet->genJet_p.Pt();
        }

        //trigger 3 cut
        if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
            tagged_trig3.Ht += puppiJet->genJet_p.Pt();
        }

    }
  }
  
  //untagged filling
  if(nTruthForCut >= minJet){
    for(const nTupleAnalysis::jetPtr& puppiJet : event->puppiJets){
        if(fabs(puppiJet->genJet_p.Eta()) > eta_cut) continue;
        //cut on pt
        if(puppiJet->genJet_p.Pt()       < pt_cut)       continue; 
        preCut.Ht += puppiJet->genJet_p.Pt();                        //precuts should be moved to correspond with untagged events

        //Cut on L1
        if(bsetList[triggerBit_L1[0]][triggerBit_L1[1]] != 1) continue;
        untagged_L1.Ht += puppiJet->genJet_p.Pt();
        
        //trigger 1 cut
        if(bsetList[triggerBit_1[0]][triggerBit_1[1]]==1){
            untagged_trig1.Ht += puppiJet->genJet_p.Pt();
        }

        //trigger 2 cut
        if(bsetList[triggerBit_2[0]][triggerBit_2[1]]==1){
            untagged_trig2.Ht += puppiJet->genJet_p.Pt();
        }

        //trigger 3 cut
        if(bsetList[triggerBit_3[0]][triggerBit_3[1]]==1){
            untagged_trig3.Ht += puppiJet->genJet_p.Pt();
        }

    }
  }

    //Fill Ht's
    if(preCut.Ht != 0){
        mass_preCut -> FillHt(preCut.Ht);
        }
    if(noL1_deepCut.Ht != 0){
        deepCut_noL1-> FillHt(noL1_deepCut.Ht);
        }
    if(tagged_L1_deepCut.Ht!=0){
        L1_deepCut_tagged  -> FillHt(tagged_L1_deepCut.Ht);
        }
    if(tagged_trig1.Ht!=0){
        trig1_tagged       -> FillHt(tagged_trig1.Ht);
    }
    if(tagged_trig2.Ht!=0){
        trig2_tagged       -> FillHt(tagged_trig2.Ht);
    }
    if(tagged_trig3.Ht!=0){
        trig3_tagged       -> FillHt(tagged_trig3.Ht);
    }
    if(untagged_L1.Ht!=0){
        L1_untagged  -> FillHt(untagged_L1.Ht);
        }
    if(untagged_trig1.Ht!=0){
        trig1       -> FillHt(untagged_trig1.Ht);
    }
    if(untagged_trig2.Ht!=0){
        trig2       -> FillHt(untagged_trig2.Ht);
    }
    if(untagged_trig3.Ht!=0){
        trig3       -> FillHt(untagged_trig3.Ht);
    }

   // }
  //}
  
  // Make sure at least four jets exist in Event  
  if(nTruthForCut < minJet      ){
    if(debug) cout << "Fail NJet Cut" << endl;
    return 0;
  }

  bool doOnlineBTagCut = true;
  if(doOnlineBTagCut){
  
    if(nTruthTaggedForCut < minBJet) { 
      if(debug) cout << "Fail NBJet Cut" << endl;
      return 0;
    }
    cutflow->Fill("passNBJetCut", eventWeight);
  }

  //
  // Fill All events
  //
  if(debug) cout << "Fill All Events " << endl;
 

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
        }
    }
  }

  return 0;
}

HH4bAnalysis::~HH4bAnalysis(){
  cout << "HH4bAnalysis::destroyed" << endl;
}


