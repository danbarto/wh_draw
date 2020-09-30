#include "core/wh_functions.hpp"
#include <algorithm> //std::min
#include <math.h>

#include "TFile.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"


using namespace std;

namespace WH_Functions{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Miscellaneous
float medDeepCSV2016 = 0.6321;
float medDeepCSV2017 = 0.4941;
float medDeepCSV2018 = 0.4184;

float deepTag2016 = 0.8945;
float deepTag2017 = 0.8695;
float deepTag2018 = 0.8365;

  const NamedFunc pass_any_mt_variation("pass_any_mt_variation",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.mt_met_lep()>150. || b.mt_met_lep_jdown()>150. || b.mt_met_lep_jup()>150. || b.mt_met_lep_resdown()>150. || b.mt_met_lep_resup()>150.  ) return 1.;
  else return 0.;
  });

  const NamedFunc pass_any_met_variation("pass_any_met_variation",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.pfmet()>125. || b.pfmet_jdown()>125. || b.pfmet_jup()>125. || b.pfmet_resdown()>125. || b.pfmet_resup()>125.  ) return 1.;
  else return 0.;
  });

  const NamedFunc pass_any_njets_variation("pass_any_njets_variation",[](const Baby &b) -> NamedFunc::ScalarType{
  if( (b.ngoodjets()>=2 && b.ngoodjets()<=3)
      || (b.jup_ngoodjets()>=2 && b.jup_ngoodjets()<=3)
      || (b.jdown_ngoodjets()>=2 && b.jdown_ngoodjets()<=3)
      // || (b.resup_ngoodjets()>=2 && b.resup_ngoodjets()<=3)
      // || (b.resdown_ngoodjets()>=2 && b.resdown_ngoodjets()<=3)
    ) return 1.;
  else return 0.;
  });

  //Not really correct
      const NamedFunc pass_any_nb_variation("pass_any_nb_variation",[](const Baby &b) -> NamedFunc::ScalarType{
  if( (b.ngoodbtags()>=1)
      // || (b.jup_ngoodbtags()>=1)
      // || (b.jdown_ngoodbtags()>=1 )
      // || (b.resup_ngoodbtags()>=1 )
      // || (b.resdown_ngoodbtags()>=1)
    ) return 1.;
  else return 0.;
  });

  const NamedFunc nanoWeight("nanoWeight",[](const Baby &b) -> NamedFunc::ScalarType{
    //For 06_12 production
    float weight=1;
    //w_lumi_scale1fb is a weight unique for each MC sample (formed from n_events and xsec)
    //it can be negative, so must use absolute value.

    //TTJets_1lep_tbar_met80_a18v1 
    if (abs(0.000679649-abs(b.w_lumi_scale1fb()))<0.000001) weight = 1.020;
    //TTJets_1lep_tbar_f17v2
    else if (abs(abs(b.w_lumi_scale1fb())-0.00324188)<0.00001) weight = 1.014;
    // in 2017, some mass points have no nano
    else if (b.year()==2017 && b.mass_stop()==200  && b.mass_lsp()==25)  weight = 1/0.989598335734 ;
    else if (b.year()==2017 && b.mass_stop()==275  && b.mass_lsp()==75)  weight = 1/0.957338965153 ;
    else if (b.year()==2017 && b.mass_stop()==300  && b.mass_lsp()==50)  weight = 1/0.961527621195 ;
    else if (b.year()==2017 && b.mass_stop()==375  && b.mass_lsp()==25)  weight = 1/0.951265229616 ;
    else if (b.year()==2017 && b.mass_stop()==675  && b.mass_lsp()==1)   weight = 1/0.97484237468  ;
    else if (b.year()==2017 && b.mass_stop()==675  && b.mass_lsp()==225) weight = 1/0.940583363342 ;
    else if (b.year()==2017 && b.mass_stop()==725  && b.mass_lsp()==150) weight = 1/0.969029492202 ;
    else if (b.year()==2017 && b.mass_stop()==725  && b.mass_lsp()==275) weight = 1/0.962572303505 ;
    else if (b.year()==2017 && b.mass_stop()==725  && b.mass_lsp()==325) weight = 1/0.96882402384  ;
    else if (b.year()==2017 && b.mass_stop()==725  && b.mass_lsp()==550) weight = 1/0.991338321477 ;
    else if (b.year()==2017 && b.mass_stop()==725  && b.mass_lsp()==575) weight = 1/0.989381003202 ;
    else if (b.year()==2017 && b.mass_stop()==825  && b.mass_lsp()==175) weight = 1/0.959898615009 ;
    else if (b.year()==2017 && b.mass_stop()==875  && b.mass_lsp()==400) weight = 1/0.972290031009 ;
    else if (b.year()==2017 && b.mass_stop()==900  && b.mass_lsp()==600) weight = 1/0.967833570056 ;
    else if (b.year()==2017 && b.mass_stop()==925  && b.mass_lsp()==25)  weight = 1/0.972826712723 ;
    else if (b.year()==2017 && b.mass_stop()==950  && b.mass_lsp()==550) weight = 1/0.972117321249 ;
    else if (b.year()==2017 && b.mass_stop()==1000 && b.mass_lsp()==125) weight = 1/0.970782940802 ;
    else if (b.year()==2017 && b.mass_stop()==1000 && b.mass_lsp()==425) weight = 1/0.974963353855 ;

    return weight;
  });


  const NamedFunc failTauVetos("failTauVetos",[](const Baby &b) -> NamedFunc::ScalarType{
    if(!b.PassTauVeto() || !b.PassTauVeto()) return 1.;
    else return 0.;
    });
  const NamedFunc yearWeight("yearWeight",[](const Baby &b) -> NamedFunc::ScalarType{
    
    float weight=0;
    if(b.genmet()==-9999) weight =1.;//DATA
    else{
      float totalLumi = 35.9+41.6+59.7;
    
      if(b.year()==2016) weight = 35.9 / totalLumi;
      else if(b.year()==2017)  weight = 41.6 / totalLumi;
     else if(b.year()==2018)  weight = 59.7 / totalLumi;
    }
    return weight;
    });

  //Temporary hack
  const NamedFunc ST_up("ST_up",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int nTop = 0; 
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i)) == 6) nTop++;
      }
      if(nTop==1) return 1.5;
      else return 1.;
  });
  const NamedFunc ST_off("ST_off",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int nTop = 0; 
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i)) == 6) nTop++;
      }
      if(nTop==1) return 0.;
      else return 1.;
  });
  const NamedFunc ST_down("ST_down",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int nTop = 0; 
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i)) == 6) nTop++;
      }
      if(nTop==1) return 0.5;
      else return 1.;
  });
  const NamedFunc fake_up("fake_up",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int ndeepfakes=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        //if (b.ak4pfjets_deepCSV()->at(i) > 0.6324 && abs(b.ak4pfjets_parton_flavor()->at(i)) != 5) ndeepfakes++;
	if(abs(b.ak4pfjets_parton_flavor()->at(i))!=5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))) ndeepfakes++;
      }
      if(ndeepfakes>0) return 1.5;
      else return 1.;
  });

  const NamedFunc fake_down("fake_down",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int ndeepfakes=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
	if(abs(b.ak4pfjets_parton_flavor()->at(i))!=5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))) ndeepfakes++;
	//if (b.ak4pfjets_deepCSV()->at(i) > 0.6324 && abs(b.ak4pfjets_parton_flavor()->at(i)) != 5) ndeepfakes++;
      }
      if(ndeepfakes>0) return 0.5;
      else return 1.;
  });

  const NamedFunc VV_up("VV_up",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int nBoson = 0; 
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if( (abs(b.gen_id()->at(i)) == 23 || abs(b.gen_id()->at(i)) == 24) && (abs(b.gen_motherid()->at(i))!=6) ) nBoson++;
      }
      if(nBoson>=2) return 1.25;
      else return 1.;
  });
  const NamedFunc VV_down("VV_down",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      int nBoson = 0; 
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if( (abs(b.gen_id()->at(i)) == 23 || abs(b.gen_id()->at(i)) == 24) && (abs(b.gen_motherid()->at(i))!=6) ) nBoson++;
      }
      if(nBoson>=2) return 0.75;
      else return 1.;
  });

    const NamedFunc ttbar_genmet_fix("ttbar_genmet_fix",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
      //genmet 200=1
      //genmet 540=1.4
      if(b.year()==2016 || b.genmet() <125) return 1.;

      else{
        float ratio_to_2016 = ((1.4-1)/(540-200)) * (b.genmet()-200) + 1.;
        if(b.genmet() > 540) ratio_to_2016=1.4;
        return 1./ratio_to_2016;
      }
  });
      const NamedFunc ttbar_genmet_antifix("ttbar_genmet_antifix",[](const Baby &b) -> NamedFunc::ScalarType{
      if(b.genmet()==-9999) return 1.; //data
    //genmet 200=1
    //genmet 540=1.4
    if(b.year()==2016 || b.genmet() <125) return 1.;

    else{
      float ratio_to_2016 = ((1.4-1)/(540-200)) * (b.genmet()-200) + 1.;
      if(b.genmet() > 540) ratio_to_2016=1.4;
      return ratio_to_2016;
    }
});

  const NamedFunc genmct("genmct",[](const Baby &b) -> NamedFunc::ScalarType{
  float gen_mct=0;
    for (unsigned i(0); i<b.gen_pt()->size(); i++){
    if (abs(b.gen_id()->at(i)) == 5 && abs(b.gen_motherid()->at(i)) == 6 ){
      for (unsigned j(i+1); j<b.gen_pt()->size(); j++){
         if (abs(b.gen_id()->at(j)) == 5 && abs(b.gen_motherid()->at(j)) == 6 ){
          gen_mct = sqrt(2*b.gen_pt()->at(j)*b.gen_pt()->at(i)*(1+cos(deltaPhi(b.gen_phi()->at(i),b.gen_phi()->at(j)))));
          break;
         }
      }
      break;
    }
  }
  return gen_mct;
  });

  const NamedFunc mct_genpt("mct_genpt",[](const Baby &b) -> NamedFunc::ScalarType{
  float mctgenpt=0;
    for (unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
    if (abs(b.ak4pfjets_parton_flavor()->at(i)) == 5){
      for (unsigned j(i+1); j<b.ak4pfjets_pt()->size(); j++){
         if (abs(b.ak4pfjets_parton_flavor()->at(j)) == 5){
          mctgenpt = sqrt(2*b.ak4pfjets_genpt()->at(j)*b.ak4pfjets_genpt()->at(i)*(1+cos(deltaPhi(b.ak4pfjets_phi()->at(i),b.ak4pfjets_phi()->at(j)))));
          break;
         }
      }
      break;
    }
  }
  return mctgenpt;
  });

  const NamedFunc nTightb("nTightb",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.8953*(b.year()==2016) + 0.8001*(b.year()==2017) + 0.7527*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });
   const NamedFunc nTightb_jup("nTightb_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt_jup()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.8953*(b.year()==2016) + 0.8001*(b.year()==2017) + 0.7527*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });
   const NamedFunc nTightb_jdown("nTightb_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt_jdown()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.8953*(b.year()==2016) + 0.8001*(b.year()==2017) + 0.7527*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });


   const NamedFunc nMedb("nMedb",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

   const NamedFunc nMedb_jup("nMedb_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_pt_jup()->size(); i++){
        if (b.ak4pfjets_pt_jup()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

   const NamedFunc nMedb_jdown("nMedb_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_pt_jdown()->size(); i++){
        if (b.ak4pfjets_pt_jdown()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

   const NamedFunc nLooseb("nLooseb",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.2217*(b.year()==2016) + 0.1522*(b.year()==2017) + 0.1241*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

      const NamedFunc nLooseb_jup("nLooseb_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt_jup()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.2217*(b.year()==2016) + 0.1522*(b.year()==2017) + 0.1241*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

       const NamedFunc nLooseb_jdown("nLooseb_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      int nb=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt_jdown()->at(i)>30. && b.ak4pfjets_deepCSV()->at(i) > (0.2217*(b.year()==2016) + 0.1522*(b.year()==2017) + 0.1241*(b.year()==2018))){
           nb++;
        }
      }
      return nb;
    });

  const NamedFunc nEventsGluonSplit("nEventsGluonSplit",[](const Baby &b) -> NamedFunc::ScalarType{
    int nevent = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_motherid()->at(i)==21){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)&&b.gen_motheridx()->at(j)==b.gen_motheridx()->at(i)&&(b.gen_pt()->at(j)>30||b.gen_pt()->at(i)>30)){
            nevent++;
          }//Close if statement for opposite b, same mother index, pt cut
        }//Close for loop over second half of particles in event
      }//Close if statement that find b with gluon mother
    }//Close for loop over all particles in event
    return nevent;

    });

  const NamedFunc nModEventsGluonSplit("nModEventsGluonSplit",[](const Baby &b) -> NamedFunc::ScalarType{
    int nevent = 0;
    int ngenb = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5){
        ngenb++;
      }
      if(abs(b.gen_id()->at(i))==5&&b.gen_motherid()->at(i)==21&&b.gen_pt()->at(i)>20){
        nevent++;
      }//Close if statement that find b with gluon mother
    }//Close for loop over all particles in event

    if(ngenb==0){
      nevent = -1;
    }

    return nevent;

    });

  const NamedFunc NHighPtNu("NHighPtNu",[](const Baby &b) -> NamedFunc::ScalarType{
      int nnu=0;
        for (unsigned i(0); i<b.gen_pt()->size(); i++){
        if (abs(b.gen_motherid()->at(i))==24 && ( abs(b.gen_id()->at(i)) == 12 || abs(b.gen_id()->at(i)) == 14 || abs(b.gen_id()->at(i)) == 16) && b.gen_pt()->at(i) > 200 ) nnu++;
      }
      return nnu;
      });

  const NamedFunc HighNuPt("HighNuPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float max_nupt=0;
        for (unsigned i(0); i<b.gen_pt()->size(); i++){
         if (abs(b.gen_motherid()->at(i))==24 && ( abs(b.gen_id()->at(i)) == 12 || abs(b.gen_id()->at(i)) == 14 || abs(b.gen_id()->at(i)) == 16) &&  b.gen_pt()->at(i) > max_nupt) max_nupt = b.gen_pt()->at(i);
      }
      return max_nupt;
      });


  const NamedFunc zpt("zpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float z_pt=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 23) z_pt = b.gen_pt()->at(i);
    }
    return z_pt;
    });

  const NamedFunc wpt("wpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float w_pt=-1;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 24) w_pt = b.gen_pt()->at(i);
    }
    return w_pt;
    });

  
  const NamedFunc wpt_lnu("wpt_lnu",[](const Baby &b) -> NamedFunc::ScalarType{
    // this only works like this because the gen collection is so pruned. otherwise a more careful check of the lepton/neutrino history would be necessary
    float w_pt=-1;
    TLorentzVector n1;
    TLorentzVector l1;
    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 11 ||  abs(b.gen_id()->at(i)) == 13 || abs(b.gen_id()->at(i)) == 15 ){
        l1.SetPtEtaPhiM(b.gen_pt()->at(i),b.gen_eta()->at(i),b.gen_phi()->at(i),0);
      }
      if ( abs(b.gen_id()->at(i)) == 12 ||  abs(b.gen_id()->at(i)) == 14 || abs(b.gen_id()->at(i)) == 16 ){
        n1.SetPtEtaPhiM(b.gen_pt()->at(i),b.gen_eta()->at(i),b.gen_phi()->at(i),0);
      }
    }
    w_pt = (l1+n1).Pt(); 
    return w_pt;
    });

  const NamedFunc higgs_pt("higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float h_pt=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 25) h_pt = b.gen_pt()->at(i);
    }
    return h_pt;
    });

  // count the number of b's in the fat jet
  const NamedFunc nBInFatJet("nBInFatJet",[](const Baby &b) -> NamedFunc::ScalarType{
      int nBinFat=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_pt()->at(i) > 200){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if (abs(b.ak4pfjets_hadron_flavor()->at(j))==5){
                        nBinFat++;
                    }
                }
            }
        }
      }
      return nBinFat;
    });

  const NamedFunc higgsMistagSF("higgsMistagSF",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PData=1.;
      float PMC=1.;
      float eff=0.;
      float SF=1.;
      float SF_unc=0.;
      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
        int nBinFat=0;
      	if (b.FatJet_pt_nom()->at(i) > 250){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
                        nBinFat++;
                    }
                }
            }
        }
        if (nBinFat==0){
            eff=0.0;
            SF=1.;
            SF_unc=0.;
        }
        else if (nBinFat==1){
            eff=0.25;
            SF=1.13;
            SF_unc=0.28;
        }
        else if (nBinFat==2){
            eff=0.67;
            SF=0.95;
            SF_unc=0.08;
        }
        //else {eff=0.;} // this should really never happen

        int delta;
        delta=0;
        if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018) ) && eff>0){
            PMC *= eff;
            PData *= eff*(SF+delta*SF_unc);
        }
        else {
            PMC *= (1-eff);
            PData *= (1-eff*(SF+delta*SF_unc));
        }
      }
      float weight=PData/PMC;
      return weight;
    });


  const NamedFunc higgsMistagSFUp("higgsMistagSFUp",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PData=1.;
      float PMC=1.;
      float eff=0.;
      float SF=1.;
      float SF_unc=0.;
      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
        int nBinFat=0;
      	if (b.FatJet_pt_nom()->at(i) > 250){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
                        nBinFat++;
                    }
                }
            }
        }
        if (nBinFat==0){
            eff=0.0;
            SF=1.;
            SF_unc=0.;
        }
        else if (nBinFat==1){
            eff=0.24;
            SF=1.13;
            SF_unc=0.28;
        }
        else if (nBinFat==2){
            eff=0.67;
            SF=0.95;
            SF_unc=0.08;
        }
        //else {eff=0.;} // this should really never happen

        int delta;
        delta=1;
        if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018) ) && eff>0 ){
            PMC *= eff;
            PData *= eff*(SF+delta*SF_unc);
        }
        else {
            PMC *= (1-eff);
            PData *= (1-eff*(SF+delta*SF_unc));
        }
      }
      float weight=PData/PMC;
      return weight;
    });

  const NamedFunc higgsMistagSFDown("higgsMistagSFDown",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PData=1.;
      float PMC=1.;
      float eff=0.;
      float SF=1.;
      float SF_unc=0.;
      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
        int nBinFat=0;
      	if (b.FatJet_pt_nom()->at(i) > 250){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
                        nBinFat++;
                    }
                }
            }
        }
        if (nBinFat==0){
            eff=0.0;
            SF=1.;
            SF_unc=0.;
        }
        else if (nBinFat==1){
            eff=0.24;
            SF=1.13;
            SF_unc=0.28;
        }
        else if (nBinFat==2){
            eff=0.67;
            SF=0.95;
            SF_unc=0.08;
        }
        //else {eff=0.} // this should really never happen

        int delta;
        delta=-1;
        if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018) ) && eff>0. ){
            PMC *= eff;
            PData *= eff*(SF+delta*SF_unc);
        }
        else {
            PMC *= (1-eff);
            PData *= (1-eff*(SF+delta*SF_unc));
        }
      }
      float weight=PData/PMC;
      return weight;
    });


  // get the mistag probablility for higgs
  const NamedFunc higgsMistagProb("higgsMistagProb",[](const Baby &b) -> NamedFunc::ScalarType{
      float prob=1;
      float mistag=0;
      bool lepsInFatJet=false;
      // load the map. need to fix hardcoding the path
      TFile *mapFile = new TFile("/home/users/dspitzba/WH/CMSSW_10_2_9/src/WH_studies/Analysis/python/eff_pt_mass_allYears_QCD_combined.root", "READ");
      TString eff_2b = "eff_pt_mass_2b";
      TString eff_1b = "eff_pt_mass_1b";
      TString eff_0b = "eff_pt_mass_0b";
      if (b.year()==2016 || true){
        eff_2b += TString("_2016");
        eff_1b += TString("_2016");
        eff_0b += TString("_2016");
      }
      else if (b.year()==2017){
        eff_2b += TString("_2017");
        eff_1b += TString("_2017");
        eff_0b += TString("_2017");
      }
      else if (b.year()==2018){
        eff_2b += TString("_2018");
        eff_1b += TString("_2018");
        eff_0b += TString("_2018");
      }
      TH2* effMap_2b = 0;
      TH2* effMap_1b = 0;
      TH2* effMap_0b = 0;
      mapFile->GetObject(eff_2b, effMap_2b);
      mapFile->GetObject(eff_1b, effMap_1b);
      mapFile->GetObject(eff_0b, effMap_0b);
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
        lepsInFatJet=false;
        int nBInFat=0;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet){
            // skip fat jet if leptons are found inside
            continue;
        }
        //if (b.weight()>0.01){
        //    std::cout << "Found large weight" << std::endl;
        ////    continue;
        //}
      	if (b.ak8pfjets_pt()->at(i) > 200){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    // have to use b-tag status as true flavor is not stored for the sub jets used for the efficiency measurement (nanoAOD/nanoAOD-tools)
                    if (b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){
                        nBInFat++;
                    }
                }
            }
            // get the rate for the fat jet
            if (nBInFat==2){
                mistag = effMap_2b->GetBinContent(effMap_2b->GetXaxis()->FindBin(b.ak8pfjets_pt()->at(i)), effMap_2b->GetYaxis()->FindBin(b.ak8pfjets_m()->at(i)));
            }
            else if (nBInFat==1){
                mistag = effMap_1b->GetBinContent(effMap_1b->GetXaxis()->FindBin(b.ak8pfjets_pt()->at(i)), effMap_1b->GetYaxis()->FindBin(b.ak8pfjets_m()->at(i)));
                //mistag = 0;
            }
            else {
                mistag = effMap_0b->GetBinContent(effMap_0b->GetXaxis()->FindBin(b.ak8pfjets_pt()->at(i)), effMap_0b->GetYaxis()->FindBin(b.ak8pfjets_m()->at(i)));
                //mistag = 0;
            }
            prob = prob*(1-mistag);
        }
      }
      // clean up
      mapFile->Close();
      delete mapFile;
      //std::cout << "prob: " << (1-prob) << std::endl;
      return 1-prob;
    });

  // get the number of b-tagged jets in fat jets
  const NamedFunc nBTagInFat("nBTagInFat",[](const Baby &b) -> NamedFunc::ScalarType{
      int nBInFat=0;
      bool lepsInFatJet=false;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
        lepsInFatJet=false;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet){
            //std::cout << "Found leptons in my fat jet. Passing." << std::endl;
            continue;
        }
      	if (b.ak8pfjets_pt()->at(i) > 200){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if (b.ak4pfjets_deepCSV()->at(j)>(0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))){
                        nBInFat++;
                    }
                }
            }
        }
      }
      return nBInFat;
    });



  const NamedFunc OneDoubleBFat("OneDoubleBFat",[](const Baby &b) -> NamedFunc::ScalarType{
      int nBInFat=0;
      int nDoubleBFat=0;
      bool lepsInFatJet=false;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
        int nBinThisFat=0;
        lepsInFatJet=false;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet){
            //std::cout << "Found leptons in my fat jet. Passing." << std::endl;
            continue;
        }
      	if (b.ak8pfjets_pt()->at(i) > 200){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if (b.ak4pfjets_deepCSV()->at(j)>(0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))){
                        nBInFat++;
                        nBinThisFat++;
                    }
                }
            }
        }
        if (nBinThisFat==2) nDoubleBFat++;
      }
      return nDoubleBFat==1;
    });

  const NamedFunc OneDoubleBFatNanoV0("OneDoubleBFatNanoV0",[](const Baby &b) -> NamedFunc::ScalarType{
      int nDoubleBFat=0;
      bool lepsInFatJet=false;
      for (unsigned i(0); i<b.FatJet_pt()->size(); i++){
        int nBinThisFat=0;
        lepsInFatJet=false;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet){
            //std::cout << "Found leptons in my fat jet. Passing." << std::endl;
            continue;
        }
      	if (b.FatJet_pt()->at(i) > 200){
            for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
                if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
                    if (b.ak4pfjets_deepCSV()->at(j)>(0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))){
                        nBinThisFat++;
                    }
                }
            }
        }
        if (nBinThisFat==2) nDoubleBFat++;
      }
      return nDoubleBFat==1;
    });

  const NamedFunc OneDoubleBFatNano("OneDoubleBFatNano",[](const Baby &b) -> NamedFunc::ScalarType{
      int nDoubleBFat=0;
      bool lepsInFatJet=false;
      for (unsigned i(0); i<b.FatJet_pt()->size(); i++){
        int nBinThisFat=0;
        lepsInFatJet=false;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet){
            //std::cout << "Found leptons in my fat jet. Passing." << std::endl;
            continue;
        }
      	if (b.FatJet_pt()->at(i) > 200){
            //std::cout << "SubJet index 1 " << b.FatJet_subJetIdx1()->at(i) << std::endl;
            //std::cout << "SubJet index 2 " << b.FatJet_subJetIdx2()->at(i) << std::endl;
            //std::cout << "WTF" << std::endl;
            //std::cout << "SubJet length " << b.nSubJet() << std::endl; // I don't yet store the SubJets... only in the next round. MCT enforces this anyway?!
            if (b.FatJet_subJetIdx1()->at(i)>=0 && b.SubJet_btagDeepB()->size()>b.FatJet_subJetIdx1()->at(i)){
                if (b.SubJet_btagDeepB()->at( b.FatJet_subJetIdx1()->at(i) ) > (0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))){
                    nBinThisFat++;
                }
            }
            if (b.FatJet_subJetIdx2()->at(i)>=0 && b.SubJet_btagDeepB()->size()>b.FatJet_subJetIdx2()->at(i)){
                if (b.SubJet_btagDeepB()->at( b.FatJet_subJetIdx2()->at(i) ) > (0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))){
                    nBinThisFat++;
                }
            }
        }
        if (nBinThisFat==2) nDoubleBFat++;
      }
      return nDoubleBFat==1;
    });

  // nFatJets with pt>250
  const NamedFunc nFatJet250("nFatJet250",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      //      for (unsigned i(0); i<b.FatJet_pt()->size(); i++){
      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
	//if (b.FatJet_pt()->at(i) > 250) nloose++;
	if (b.FatJet_pt_nom()->at(i) > 250) nloose++;
      }
      return nloose;
    });

  std::vector<unsigned> binsSF = {200,300,400,500,600};
//  std::vector<float> const signalSF2016 = {0.99, 1.00, 0.97, 0.91, 0.95};
//  std::vector<float> const signalSF2016_unc = {0.04, 0.06, 0.03, 0.05, 0.04};
//  std::vector<float> const signalSF2017 = {1.05, 1.05, 1.01, 1.06, 1.13};
//  std::vector<float> const signalSF2017_unc = {0.04, 0.02, 0.03, 0.03, 0.05};
//  std::vector<float> const signalSF2018 = {1.32, 1.35, 1.22, 1.31, 1.30};
//  //  std::vector<float> const signalSF2018_unc = {0.05, 0.03, 0.03, 0.04, 0.04};
//  std::vector<float> const signalSF2018_unc = {0.3, 0.3, 0.3, 0.3, 0.3};

//  std::vector<float> const signalSF_LP_2016 = {0.99, 1.00, 0.97, 0.91, 0.95};     // old
//  std::vector<float> const signalSF_LP_2016_unc = {0.04, 0.06, 0.03, 0.05, 0.04}; // old
//  std::vector<float> const signalSF_LP_2017 = {1.05, 1.05, 1.01, 1.06, 1.13};	  // old
//  std::vector<float> const signalSF_LP_2017_unc = {0.04, 0.02, 0.03, 0.03, 0.05}; // old
//  std::vector<float> const signalSF_LP_2018 = {1.32, 1.35, 1.22, 1.31, 1.30};	  // old
//  std::vector<float> const signalSF_LP_2018_unc = {0.05, 0.03, 0.03, 0.04, 0.04}; // old
// 										  // old
//wrong  std::vector<float> const signalSF_HP_2016 = {1.06, 1.06, 1.00, 0.98, 1.01};	  // old
//wrong  std::vector<float> const signalSF_HP_2016_unc = {0.07, 0.06, 0.09, 0.08, 0.05}; // old
//wrong  std::vector<float> const signalSF_HP_2017 = {0.92, 0.96, 0.99, 1.05, 1.17};	  // old
//wrong  std::vector<float> const signalSF_HP_2017_unc = {0.08, 0.09, 0.07, 0.07, 0.11}; // old
//wrong  std::vector<float> const signalSF_HP_2018 = {0.97, 0.85, 1.04, 1.03, 1.03};	  // old
//wrong  std::vector<float> const signalSF_HP_2018_unc = {0.12, 0.07, 0.09, 0.08, 0.09}; // old

  std::vector<float> const signalSF_LP_2016 = {1.04, 0.97, 1.02, 0.95, 0.98};     // new
  std::vector<float> const signalSF_LP_2016_unc = {0.11, 0.09, 0.08, 0.10, 0.10}; // new
  std::vector<float> const signalSF_LP_2017 = {1.01, 1.10, 1.03, 0.95, 0.96};	  // new
  std::vector<float> const signalSF_LP_2017_unc = {0.11, 0.13, 0.10, 0.12, 0.10}; // new
  std::vector<float> const signalSF_LP_2018 = {0.91, 1.2, 1.05, 0.93, 1.05};	  // new
  std::vector<float> const signalSF_LP_2018_unc = {0.12, 0.15, 0.09, 0.10, 0.12}; // new
 										  // new
  std::vector<float> const signalSF_HP_2016 = {1.06, 1.06, 1.00, 0.98, 1.01};	  // new
  std::vector<float> const signalSF_HP_2016_unc = {0.07, 0.06, 0.09, 0.08, 0.05}; // new
  std::vector<float> const signalSF_HP_2017 = {0.92, 0.96, 0.99, 1.05, 1.17};	  // new
  std::vector<float> const signalSF_HP_2017_unc = {0.08, 0.09, 0.07, 0.07, 0.11}; // new
  std::vector<float> const signalSF_HP_2018 = {0.97, 0.85, 1.04, 1.03, 1.03};	  // new
  std::vector<float> const signalSF_HP_2018_unc = {0.12, 0.07, 0.09, 0.08, 0.09}; // new

  float const LP_cutoff = 0.8;
  float const HP_cutoff = 0.97;

  std::vector<unsigned> binsEff = {200,400};		     
  std::vector<float> const fastEff2016_750_1 = {0.702, 0.838};
  std::vector<float> const fastEff2017_750_1 = {0.74, 0.865};
  std::vector<float> const fastEff2018_750_1 = {0.816, 0.865};
  std::vector<float> const fullEff2016_750_1 = {0.742, 0.84};
  std::vector<float> const fullEff2018_750_1 = {0.824, 0.909};
  //  std::vector<float> const fullEff2017_750_1 = fullEff2018_750_1;//placeholder
  std::vector<float> const fullEff2017_750_1 = {0.7933, 0.88};//placeholder

  const NamedFunc signalHiggsMistagSF("signalHiggsMistagSF",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PFast = 1.;
      float PFull = 1.;
      float PData=1.;
      float PMC=1.;
      float effFast=0.;
      float effFull=0.;
      float SFFast = 1.;
      float SFFast_unc = 0.;
      float SFFull = 1.;
      float SFFull_unc = 0.;
      
      // Find Higgs idx
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }

      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
	// Gen Higgs
	if ((b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0) && (deltaR(b.FatJet_eta()->at(i), b.FatJet_phi()->at(i), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8)) { // consider only Fat Jets with gen H
	  float pt = 0;
	  if ((pt = b.FatJet_pt_nom()->at(i)) > 250) {
	    unsigned binSF = 0;
	    unsigned binEff = 0;
	    for(unsigned iBinSF(0); iBinSF<binsSF.size(); iBinSF++) { // Lookup pt bin for SF
	      if(binsSF.size() == (iBinSF+1)) {
		binSF = iBinSF;
	      } else if ((binsSF[iBinSF] < pt) && (binsSF[iBinSF+1] > pt)) {
		binSF = iBinSF;
		break;
	      }
	    }
	    for(unsigned iBinEff(0); iBinEff<binsEff.size(); iBinEff++) { // Lookup pt bin for eff
	      if(binsEff.size() == (iBinEff+1)) {
		binEff = iBinEff;
	      } else if ((binsEff[iBinEff] < pt) && (binsEff[iBinEff+1] > pt)) {
		binEff = iBinEff;
		break;
	      }
	    }
	    if(2016==b.year()) {
	      effFast = fastEff2016_750_1[binEff];
	      effFull = fullEff2016_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2016[binSF];
		SFFull_unc = signalSF_HP_2016_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2016[binSF];
		SFFull_unc = signalSF_LP_2016_unc[binSF];
	      }
//	      SFFull = signalSF2016[binSF];
//	      SFFull_unc = signalSF2016_unc[binSF];
	    } else if(2017==b.year()) {
	      effFast = fastEff2017_750_1[binEff];
	      effFull = fullEff2017_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2017[binSF];
		SFFull_unc = signalSF_HP_2017_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2017[binSF];
		SFFull_unc = signalSF_LP_2017_unc[binSF];
	      }
//	      SFFull = signalSF2017[binSF];
//	      SFFull_unc = signalSF2017_unc[binSF];
	    } else if(2018==b.year()) {
	      effFast = fastEff2018_750_1[binEff];
	      effFull = fullEff2018_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2018[binSF];
		SFFull_unc = signalSF_HP_2018_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2018[binSF];
		SFFull_unc = signalSF_LP_2018_unc[binSF];
	      }
//	      SFFull = signalSF2018[binSF];
//	      SFFull_unc = signalSF2018_unc[binSF];	    
//	    } else { // This should never happen
	    }
	    SFFast = effFull/effFast;
	    SFFast_unc = 0.05;
	    int delta;
	    delta=0;
	    if (!(effFull > 0)) {printf("effFull is %f\n", effFull);}
	    if (!(effFast > 0)) {printf("effFast is %f\n", effFast);}
	    if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018))){// && (effFull > 0) && (effFast > 0)){
	      PFast *= effFast;
	      PFull *= effFast*(SFFast+delta*SFFast_unc);
	      PMC *= effFull;
	      PData *= effFull*(SFFull+delta*SFFull_unc);
	    }
	    else {
	      PFast *= (1-effFast);
	      PFull *= (1-effFast*(SFFast+delta*SFFast_unc));
	      PMC *= (1-effFull);
	      float zero = 0.0;
	      PData *= std::max((1-effFull*(SFFull+delta*SFFull_unc)), zero);
	      //	      PData *= (1-effFull*(SFFull+delta*SFFull_unc));
	    }
	  }
	  // No Gen Higgs
	} else {
	  int nBinFat=0;
	  if (b.FatJet_pt_nom()->at(i) > 250){
	    for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
	      if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
		if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
		  nBinFat++;
		}
	      }
	    }
	  }
	  if (nBinFat==0){
	    effFull=0.0;
	    SFFull=1.;
	    SFFull_unc=0.;
	  }
	  else if (nBinFat==1){
	    effFull=0.25;
	    SFFull=1.13;
	    SFFull_unc=0.28;
	  }
	  else if (nBinFat==2){
	    effFull=0.67;
	    SFFull=0.95;
	    SFFull_unc=0.08;
	  }
	  //else {eff=0.;} // this should really never happen
	  
	  int delta;
	  delta=0;
	  if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018)) && (effFull > 0) && (effFast > 0)){
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= effFull;
	    PData *= effFull*(SFFull+delta*SFFull_unc);
	  }
	  else {
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= (1-effFull);
	    PData *= (1-effFull*(SFFull+delta*SFFull_unc));
	  }
	}
      }
      float weight=(PFull/PFast)*(PData/PMC);
      return weight;
    });

  const NamedFunc signalHiggsMistagSFUp("signalHiggsMistagSFUp",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PFast = 1.;
      float PFull = 1.;
      float PData=1.;
      float PMC=1.;
      float effFast=0.;
      float effFull=0.;
      float SFFast = 1.;
      float SFFast_unc = 0.;
      float SFFull = 1.;
      float SFFull_unc = 0.;
      
      // Find Higgs idx
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }

      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
	// Gen Higgs
	if ((b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0) && (deltaR(b.FatJet_eta()->at(i), b.FatJet_phi()->at(i), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8)) { // consider only Fat Jets with gen H
	  float pt = 0;
	  if ((pt = b.FatJet_pt_nom()->at(i)) > 250) {
	    unsigned binSF = 0;
	    unsigned binEff = 0;
	    for(unsigned iBinSF(0); iBinSF<binsSF.size(); iBinSF++) { // Lookup pt bin for SF
	      if(binsSF.size() == (iBinSF+1)) {
		binSF = iBinSF;
	      } else if ((binsSF[iBinSF] < pt) && (binsSF[iBinSF+1] > pt)) {
		binSF = iBinSF;
		break;
	      }
	    }
	    for(unsigned iBinEff(0); iBinEff<binsEff.size(); iBinEff++) { // Lookup pt bin for eff
	      if(binsEff.size() == (iBinEff+1)) {
		binEff = iBinEff;
	      } else if ((binsEff[iBinEff] < pt) && (binsEff[iBinEff+1] > pt)) {
		binEff = iBinEff;
		break;
	      }
	    }
	    if(2016==b.year()) {
	      effFast = fastEff2016_750_1[binEff];
	      effFull = fullEff2016_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2016[binSF];
		SFFull_unc = signalSF_HP_2016_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2016[binSF];
		SFFull_unc = signalSF_LP_2016_unc[binSF];
	      }
//	      SFFull = signalSF2016[binSF];
//	      SFFull_unc = signalSF2016_unc[binSF];
	    } else if(2017==b.year()) {
	      effFast = fastEff2017_750_1[binEff];
	      effFull = fullEff2017_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2017[binSF];
		SFFull_unc = signalSF_HP_2017_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2017[binSF];
		SFFull_unc = signalSF_LP_2017_unc[binSF];
	      }
//	      SFFull = signalSF2017[binSF];
//	      SFFull_unc = signalSF2017_unc[binSF];
	    } else if(2018==b.year()) {
	      effFast = fastEff2018_750_1[binEff];
	      effFull = fullEff2018_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2018[binSF];
		SFFull_unc = signalSF_HP_2018_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2018[binSF];
		SFFull_unc = signalSF_LP_2018_unc[binSF];
	      }
//	      SFFull = signalSF2018[binSF];
//	      SFFull_unc = signalSF2018_unc[binSF];	    
//	    } else { // This should never happen
	    }
//
//	    if(2016==b.year()) {
//	      effFast = fastEff2016_750_1[binEff];
//	      effFull = fullEff2016_750_1[binEff];
//	      SFFull = signalSF2016[binSF];
//	      SFFull_unc = signalSF2016_unc[binSF];
//	    } else if(2017==b.year()) {
//	      effFast = fastEff2017_750_1[binEff];
//	      effFull = fullEff2017_750_1[binEff];
//	      SFFull = signalSF2017[binSF];
//	      SFFull_unc = signalSF2017_unc[binSF];
//	    } else if(2018==b.year()) {
//	      effFast = fastEff2018_750_1[binEff];
//	      effFull = fullEff2018_750_1[binEff];
//	      SFFull = signalSF2018[binSF];
//	      SFFull_unc = signalSF2018_unc[binSF];	    
////	    } else { // This should never happen
//	    }
	    SFFast = effFull/effFast;
	    SFFast_unc = 0.05;
	    int delta;
	    delta = 1;
	    if (!(effFull > 0)) {printf("effFull is %f\n", effFull);}
	    if (!(effFast > 0)) {printf("effFast is %f\n", effFast);}
	    if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018))){// && (effFull > 0) && (effFast > 0)){
	      PFast *= effFast;
	      PFull *= effFast*(SFFast+delta*SFFast_unc);
	      PMC *= effFull;
	      PData *= effFull*(SFFull+delta*SFFull_unc);
	    }
	    else {
	      PFast *= (1-effFast);
	      PFull *= (1-effFast*(SFFast+delta*SFFast_unc));
	      PMC *= (1-effFull);
	      float zero = 0.0;
	      PData *= std::max((1-effFull*(SFFull+delta*SFFull_unc)), zero);
	    }
	  }
	  // No Gen Higgs
	} else {
	  int nBinFat=0;
	  if (b.FatJet_pt_nom()->at(i) > 250){
	    for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
	      if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
		if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
		  nBinFat++;
		}
	      }
	    }
	  }
	  if (nBinFat==0){
	    effFull=0.0;
	    SFFull=1.;
	    SFFull_unc=0.;
	  }
	  else if (nBinFat==1){
	    effFull=0.25;
	    SFFull=1.13;
	    SFFull_unc=0.28;
	  }
	  else if (nBinFat==2){
	    effFull=0.67;
	    SFFull=0.95;
	    SFFull_unc=0.08;
	  }
	  //else {eff=0.;} // this should really never happen
	  
	  int delta;
	  delta = 1;
	  if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018)) && (effFull > 0) && (effFast > 0)){
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= effFull;
	    PData *= effFull*(SFFull+delta*SFFull_unc);
	  }
	  else {
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= (1-effFull);
	    PData *= (1-effFull*(SFFull+delta*SFFull_unc));
	  }
	}
      }
      float weight=(PFull/PFast)*(PData/PMC);
      return weight;
    });

  const NamedFunc signalHiggsMistagSFDown("signalHiggsMistagSFDown",[](const Baby &b) -> NamedFunc::ScalarType{
      // based on method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
      float PFast = 1.;
      float PFull = 1.;
      float PData=1.;
      float PMC=1.;
      float effFast=0.;
      float effFull=0.;
      float SFFast = 1.;
      float SFFast_unc = 0.;
      float SFFull = 1.;
      float SFFull_unc = 0.;
      
      // Find Higgs idx
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }

      for (unsigned i(0); i<b.FatJet_pt_nom()->size(); i++){
	// Gen Higgs
	if ((b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0) && (deltaR(b.FatJet_eta()->at(i), b.FatJet_phi()->at(i), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8)) { // consider only Fat Jets with gen H
	  float pt = 0;
	  if ((pt = b.FatJet_pt_nom()->at(i)) > 250) {
	    unsigned binSF = 0;
	    unsigned binEff = 0;
	    for(unsigned iBinSF(0); iBinSF<binsSF.size(); iBinSF++) { // Lookup pt bin for SF
	      if(binsSF.size() == (iBinSF+1)) {
		binSF = iBinSF;
	      } else if ((binsSF[iBinSF] < pt) && (binsSF[iBinSF+1] > pt)) {
		binSF = iBinSF;
		break;
	      }
	    }
	    for(unsigned iBinEff(0); iBinEff<binsEff.size(); iBinEff++) { // Lookup pt bin for eff
	      if(binsEff.size() == (iBinEff+1)) {
		binEff = iBinEff;
	      } else if ((binsEff[iBinEff] < pt) && (binsEff[iBinEff+1] > pt)) {
		binEff = iBinEff;
		break;
	      }
	    }
	    if(2016==b.year()) {
	      effFast = fastEff2016_750_1[binEff];
	      effFull = fullEff2016_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2016[binSF];
		SFFull_unc = signalSF_HP_2016_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2016[binSF];
		SFFull_unc = signalSF_LP_2016_unc[binSF];
	      }
//	      SFFull = signalSF2016[binSF];
//	      SFFull_unc = signalSF2016_unc[binSF];
	    } else if(2017==b.year()) {
	      effFast = fastEff2017_750_1[binEff];
	      effFull = fullEff2017_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2017[binSF];
		SFFull_unc = signalSF_HP_2017_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2017[binSF];
		SFFull_unc = signalSF_LP_2017_unc[binSF];
	      }
//	      SFFull = signalSF2017[binSF];
//	      SFFull_unc = signalSF2017_unc[binSF];
	    } else if(2018==b.year()) {
	      effFast = fastEff2018_750_1[binEff];
	      effFull = fullEff2018_750_1[binEff];
	      if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >HP_cutoff) {
		SFFull = signalSF_HP_2018[binSF];
		SFFull_unc = signalSF_HP_2018_unc[binSF];
	      } else if (b.FatJet_deepTagMD_HbbvsQCD()->at(i) >LP_cutoff) {
		SFFull = signalSF_LP_2018[binSF];
		SFFull_unc = signalSF_LP_2018_unc[binSF];
	      }
//	      SFFull = signalSF2018[binSF];
//	      SFFull_unc = signalSF2018_unc[binSF];	    
//	    } else { // This should never happen
	    }
//
//	    if(2016==b.year()) {
//	      effFast = fastEff2016_750_1[binEff];
//	      effFull = fullEff2016_750_1[binEff];
//	      SFFull = signalSF2016[binSF];
//	      SFFull_unc = signalSF2016_unc[binSF];
//	    } else if(2017==b.year()) {
//	      effFast = fastEff2017_750_1[binEff];
//	      effFull = fullEff2017_750_1[binEff];
//	      SFFull = signalSF2017[binSF];
//	      SFFull_unc = signalSF2017_unc[binSF];
//	    } else if(2018==b.year()) {
//	      effFast = fastEff2018_750_1[binEff];
//	      effFull = fullEff2018_750_1[binEff];
//	      SFFull = signalSF2018[binSF];
//	      SFFull_unc = signalSF2018_unc[binSF];	    
////	    } else { // This should never happen
//	    }
	    SFFast = effFull/effFast;
	    SFFast_unc = 0.05;
	    int delta;
	    delta = -1;
	    if (!(effFull > 0)) {printf("effFull is %f\n", effFull);}
	    if (!(effFast > 0)) {printf("effFast is %f\n", effFast);}
	    if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018))){// && (effFull > 0) && (effFast > 0)){
	      PFast *= effFast;
	      PFull *= effFast*(SFFast+delta*SFFast_unc);
	      PMC *= effFull;
	      PData *= effFull*(SFFull+delta*SFFull_unc);
	    }
	    else {
	      PFast *= (1-effFast);
	      PFull *= (1-effFast*(SFFast+delta*SFFast_unc));
	      PMC *= (1-effFull);
	      float zero = 0.0;
	      PData *= std::max((1-effFull*(SFFull+delta*SFFull_unc)), zero);
	      //	      PData *= (1-effFull*(SFFull+delta*SFFull_unc));
	    }
	  }
	  // No Gen Higgs
	} else {
	  int nBinFat=0;
	  if (b.FatJet_pt_nom()->at(i) > 250){
	    for (unsigned j(0); j<b.ak4pfjets_eta()->size(); j++){
	      if (deltaR(b.FatJet_eta()->at(i),b.FatJet_phi()->at(i),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j))<0.8){
		if(b.ak4pfjets_deepCSV()->at(j)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
		  nBinFat++;
		}
	      }
	    }
	  }
	  if (nBinFat==0){
	    effFull=0.0;
	    SFFull=1.;
	    SFFull_unc=0.;
	  }
	  else if (nBinFat==1){
	    effFull=0.25;
	    SFFull=1.13;
	    SFFull_unc=0.28;
	  }
	  else if (nBinFat==2){
	    effFull=0.67;
	    SFFull=0.95;
	    SFFull_unc=0.08;
	  }
	  //else {eff=0.;} // this should really never happen
	  
	  int delta;
	  delta = -1;
	  if ( b.FatJet_deepTagMD_HbbvsQCD()->at(i)>(deepTag2017*(b.year()==2017) + deepTag2016*(b.year()==2016) + deepTag2018*(b.year()==2018)) && (effFull > 0) && (effFast > 0)){
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= effFull;
	    PData *= effFull*(SFFull+delta*SFFull_unc);
	  }
	  else {
	    PFast *= 1.;
	    PFull *= 1.;
	    PMC *= (1-effFull);
	    PData *= (1-effFull*(SFFull+delta*SFFull_unc));
	  }
	}
      }
      float weight=(PFull/PFast)*(PData/PMC);
      return weight;
    });

  //%%%%%%%%%%%%%%%%%%%%%%%%%
  // nFatJets with pt>250 and 0 b jets
  const NamedFunc nFatJet250b0("nFatJet250b0",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_0b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {
	if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	  //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	  int nBinFat = 0;					                
	  for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	    //	    if (deltaR(b.FatJet_eta()->at(i),			                
	    //		       b.FatJet_phi()->at(i),			                
	    if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
		       b.FatJet_phi()->at(highest_pt_idx),			                
		       b.ak4pfjets_eta()->at(j),			                
		       b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
	      if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						  + 0.6321*(b.year()==2016)
						  + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		nBinFat++;                                                      
	      }
	    }
	  }
	  if (0==nBinFat) {                                                     // Consider FatJets with 0 b-tags
	    nFatJet250_0b++;
	  }
	}
      }
      //    } 
      if(nFatJet250_0b > 1) printf("b0 fjerr: %d\n", nFatJet250_0b);
      return nFatJet250_0b;
    });

  // nFatJets with pt>250 and 1 b jet
  const NamedFunc nFatJet250b1("nFatJet250b1",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_1b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {
        old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {      
	if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	  //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	  int nBinFat = 0;					                
	  for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	    //	    if (deltaR(b.FatJet_eta()->at(i),			                
	    //		       b.FatJet_phi()->at(i),			                
	    if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
		       b.FatJet_phi()->at(highest_pt_idx),			                
		       b.ak4pfjets_eta()->at(j),			                
		       b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
	      if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						  + 0.6321*(b.year()==2016)
						  + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		nBinFat++;                                                      
	      }
	    }
	}
	  if (1==nBinFat) {                                                     // Consider FatJets with 1 b-tag
	    nFatJet250_1b++;
	  }
	}
      }
	//      } 
      if(nFatJet250_1b > 1) printf("b1 fjerr: %d\n", nFatJet250_1b);
      return nFatJet250_1b;
    });
  
  // nFatJets with pt>250 and 2 b jets
  const NamedFunc nFatJet250b2("nFatJet250b2",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_2b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {      
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {            
	if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	  //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	  int nBinFat = 0;					                
	  for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	    //	    if (deltaR(b.FatJet_eta()->at(i),			                
	    //		       b.FatJet_phi()->at(i),			                
	    if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
		       b.FatJet_phi()->at(highest_pt_idx),			                
		       b.ak4pfjets_eta()->at(j),			                
		       b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
	      if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						  + 0.6321*(b.year()==2016)
						  + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		nBinFat++;                                                      
	      }
	    }
	  }
	  if (2==nBinFat) {                                                     // Consider FatJets with 2 b-tags
	    nFatJet250_2b++;
	  }
	}
      }
      //    } 
      if(nFatJet250_2b > 1) printf("b2 fjerr: %d\n", nFatJet250_2b);
      return nFatJet250_2b;
    });

  // Higgs-tagged nFatJets with pt>250 and 0 b jets  
  const NamedFunc nHiggsFatJet250b0("nHiggsFatJet250b0",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_0b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {      
        old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {            
	if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								//      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								+ 0.8695*(b.year()==2017)
								+ 0.8365*(b.year()==2018))) {
	  if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	    //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	    int nBinFat = 0;					                
	    for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	      //	    if (deltaR(b.FatJet_eta()->at(i),			                
	      //		       b.FatJet_phi()->at(i),			                
	      if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			 b.FatJet_phi()->at(highest_pt_idx),			                
			 b.ak4pfjets_eta()->at(j),			                
			 b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						    + 0.6321*(b.year()==2016)
						    + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		  nBinFat++;                                                      
		}
	      }
	    }
	    if (0==nBinFat) {                                                     // Consider FatJets with 0 b-tags
	      nHiggsFatJet250_0b++;
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_0b > 1) printf("Higgs b0 fjerr: %d\n", nHiggsFatJet250_0b);
      return nHiggsFatJet250_0b;
    });
  
  // Higgs-tagged nFatJets with pt>250 and 1 b jet
  const NamedFunc nHiggsFatJet250b1("nHiggsFatJet250b1",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_1b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {        
	if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								//      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								+ 0.8695*(b.year()==2017)
								+ 0.8365*(b.year()==2018))) {
	  if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	    //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	    int nBinFat = 0;					                
	    for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	      //	    if (deltaR(b.FatJet_eta()->at(i),			                
	      //		       b.FatJet_phi()->at(i),			                
	      if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			 b.FatJet_phi()->at(highest_pt_idx),			                
			 b.ak4pfjets_eta()->at(j),			                
			 b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						    + 0.6321*(b.year()==2016)
						    + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		  nBinFat++;                                                      
		}
	      }
	    }
	    if (1==nBinFat) {                                                     // Consider FatJets with 1 b-tag
	      nHiggsFatJet250_1b++;
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_1b > 1) printf("Higgs b1 fjerr: %d\n", nHiggsFatJet250_1b);
      return nHiggsFatJet250_1b;
    });
    
  // Higgs-tagged nFatJets with pt>250 and 2 b jets
  const NamedFunc nHiggsFatJet250b2("nHiggsFatJet250b2",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_2b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {      
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
      if (b.FatJet_pt_nom()->size() > 0) {      
	if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								//      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								+ 0.8695*(b.year()==2017)
								+ 0.8365*(b.year()==2018))) {
	  if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	    //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	    int nBinFat = 0;					                
	    for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
	      //	    if (deltaR(b.FatJet_eta()->at(i),			                
	      //		       b.FatJet_phi()->at(i),			                
	      if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			 b.FatJet_phi()->at(highest_pt_idx),			                
			 b.ak4pfjets_eta()->at(j),			                
			 b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						    + 0.6321*(b.year()==2016)
						    + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		  nBinFat++;                                                      
		}
	      }
	    }
	    if (2==nBinFat) {                                                     // Consider FatJets with 2 b-tags
	      nHiggsFatJet250_2b++;
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_2b > 1) printf("Higgs b2 fjerr: %d\n", nHiggsFatJet250_2b);
      return nHiggsFatJet250_2b;
    });

  //%%%%%%%%%%%%%%% 
  // Considering only Fat Jets with gen H

  // nFatJets (with gen H) with pt>250
  const NamedFunc nGenHFatJet250("nGenHFatJet250",[](const Baby &b) -> NamedFunc::ScalarType{
      int nGenHFatJet250_ = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {      
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	  //	  printf("Ok I found a higgs, but is it at idx = 0?     %d\n", higgsIdx);
	}
      }
//      if (0 == higgsIdxCount) {
//	printf("No higgs in gen_id(). sample type: %d\n", b.SampleType());
//      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if (0 < higgsIdxCount) { // Only consider scenarios where higgsIdxCount was actually incremented. If you don't, and it wasn't incremented, you'd be finding deltaR for the first particle due to "unsigned higgsIdx = 0;"
	//	if (10 == b.SampleType()) {printf("Ok sample type is: %d and higgsIdx is: h%d\n", b.SampleType(), higgsIdx);}
	//	if (0 == higgsIdx) {printf("0 == higgs but not all the way in %d\n", higgsIdx);}
	if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	  //	  if (0 == higgsIdx) {printf("TOTAL I'm zero but I'm NOT all the way in %d %d\n", higgsIdx, higgsIdxCount);}				  
	  //  if (0 == higgsIdxCount) {printf("Ok so the sizes are non-zero but yet there is no higgsIdxCount. If I don't see this, then that is good\n");}
	  if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) { // consider only Fat Jets with gen H
	    if (b.FatJet_pt_nom()->size() > 0) {            
	      if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                                                                            // consider only > 250 GeV pt Fat Jets
		//		if (0 == higgsIdx) {printf("TOTAL I'm zero AND I'm all the way in %d %d %d\n", higgsIdx, higgsIdxCount, b.gen_id()->at(higgsIdx));}
//		printf("pfmet: %f\nmct: %f\nmbb: %f\nngoodbtags: %d\nnvetoleps: %d\nngoodjets: %d\nmt_met_lep: %f\n",
//		       b.pfmet(),
//		       b.mct(),
//		       b.mbb(),
//		       b.ngoodbtags(),
//		       b.nvetoleps(),
//		       b.ngoodjets(),
//		       b.mt_met_lep());
		nGenHFatJet250_++;
	      }
	    }
	  }
	}
      } //else {printf("No higgs\n");}
      return nGenHFatJet250_;
    });
  
  // Higgs-tagged nFatJets (with gen H) with pt>250 
  const NamedFunc nGenHHiggsFatJet250("nGenHHiggsFatJet250",[](const Baby &b) -> NamedFunc::ScalarType{
      int nGenHHiggsFatJet250_ = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {      
        old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	  //	  printf("Ok I found a higgs, but is it at idx = 0?     %d\n", higgsIdx);
	}	
      }
//      if (0 == higgsIdxCount) {
//	printf("No higgs in gen_id(). sample type: %d\n", b.SampleType());
//      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if (0 < higgsIdxCount) { // Only consider scenarios where higgsIdxCount was actually incremented. If you don't, and it wasn't incremented, you'd be finding deltaR for the first particle due to "unsigned higgsIdx = 0;"
	//	if (10 == b.SampleType()) {printf("Ok sample type is: %d and higgsIdx is: h%d\n", b.SampleType(), higgsIdx);}
	//	if (0 == higgsIdx) {printf("0 == higgs but not all the way in %d\n", higgsIdx);}
	if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	  //	  if (0 == higgsIdx) {printf("HIGGS I'm zero but I'm NOT all the way in %d %d\n", higgsIdx, higgsIdxCount);}				  
	  //	if (0 == higgsIdxCount) {printf("Ok so the sizes are non-zero but yet there is no higgsIdxCount. If I don't see this, then that is good\n");}	
	  if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) { // consider only Fat Jets with gen H
	    if (b.FatJet_pt_nom()->size() > 0) {            
	      if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								      + 0.8695*(b.year()==2017)
								      + 0.8365*(b.year()==2018))) {                                                         // considered only Higgs-Tagged Fat Jets
		if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                                                                          // consider only > 250 GeV pt Fat Jets
		  //		if (0 == higgsIdx) {printf("not zero %d\n", higgsIdx);}		 
		  //		  if (0 == higgsIdx) {printf("HIGGS I'm zero AND I'm all the way in %d %d %d\n", higgsIdx, higgsIdxCount, b.gen_id()->at(higgsIdx));}				  
//		  printf("pfmet: %f\nmct: %f\nmbb: %f\nngoodbtags: %d\nnvetoleps: %d\nngoodjets: %d\nmt_met_lep: %f\n",
//			 b.pfmet(),
//			 b.mct(),
//			 b.mbb(),
//			 b.ngoodbtags(),
//			 b.nvetoleps(),
//			 b.ngoodjets(),
//			 b.mt_met_lep());
		  nGenHHiggsFatJet250_++;
		}
	      }
	    }
	  }
	}
      } //else {printf("No higgs\n");}
      return nGenHHiggsFatJet250_;
    });

  // nFatJets (with gen H) with pt>250 and 0 b jets
  const NamedFunc nGenHFatJet250b0("nGenHFatJet250b0",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_0b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {
	    if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	      //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	      int nBinFat = 0;					                
	      for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		//	    if (deltaR(b.FatJet_eta()->at(i),			                
		//		       b.FatJet_phi()->at(i),			                
		if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			   b.FatJet_phi()->at(highest_pt_idx),			                
			   b.ak4pfjets_eta()->at(j),			                
			   b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		  if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						      + 0.6321*(b.year()==2016)
						      + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		    nBinFat++;                                                      
		  }
		}
	      }
	      if (0==nBinFat) {                                                     // Consider FatJets with 0 b-tags
		nFatJet250_0b++;
	      }
	    }
	  }
	}
      }
      //    } 
      if(nFatJet250_0b > 1) printf("b0 fjerr: %d\n", nFatJet250_0b);
      return nFatJet250_0b;
    });

  // nFatJets (with gen H) with pt>250 and 1 b jet
  const NamedFunc nGenHFatJet250b1("nGenHFatJet250b1",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_1b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {
        old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {      
	    if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	      //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	      int nBinFat = 0;					                
	      for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		//	    if (deltaR(b.FatJet_eta()->at(i),			                
		//		       b.FatJet_phi()->at(i),			                
		if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			   b.FatJet_phi()->at(highest_pt_idx),			                
			   b.ak4pfjets_eta()->at(j),			                
			   b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		  if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						      + 0.6321*(b.year()==2016)
						      + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		    nBinFat++;                                                      
		  }
		}
	}
	      if (1==nBinFat) {                                                     // Consider FatJets with 1 b-tag
		nFatJet250_1b++;
	      }
	    }
	  }
	}
      }
      //      } 
      if(nFatJet250_1b > 1) printf("b1 fjerr: %d\n", nFatJet250_1b);
      return nFatJet250_1b;
    });
  
  // nFatJets (with gen H) with pt>250 and 2 b jets
  const NamedFunc nGenHFatJet250b2("nGenHFatJet250b2",[](const Baby &b) -> NamedFunc::ScalarType{
      int nFatJet250_2b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {      
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {            
	    if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
	      //	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
	      int nBinFat = 0;					                
	      for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		//	    if (deltaR(b.FatJet_eta()->at(i),			                
		//		       b.FatJet_phi()->at(i),			                
		if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			   b.FatJet_phi()->at(highest_pt_idx),			                
			   b.ak4pfjets_eta()->at(j),			                
			   b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		  if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
						      + 0.6321*(b.year()==2016)
						      + 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		    nBinFat++;                                                      
		  }
		}
	      }
	      if (2==nBinFat) {                                                     // Consider FatJets with 2 b-tags
		nFatJet250_2b++;
	      }
	    }
	  }
	}
      }
      //    } 
      if(nFatJet250_2b > 1) printf("b2 fjerr: %d\n", nFatJet250_2b);
      return nFatJet250_2b;
    });
  
  // Higgs-tagged nFatJets (with gen H) with pt>250 and 0 b jets  
  const NamedFunc nGenHHiggsFatJet250b0("nGenHHiggsFatJet250b0",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_0b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;
      if (b.FatJet_pt_nom()->size() > 0) {      
        old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {            
	    if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								    //      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								    + 0.8695*(b.year()==2017)
								    + 0.8365*(b.year()==2018))) {
	      if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
		//	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
		int nBinFat = 0;					                
		for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		  //	    if (deltaR(b.FatJet_eta()->at(i),			                
		  //		       b.FatJet_phi()->at(i),			                
		  if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			     b.FatJet_phi()->at(highest_pt_idx),			                
			     b.ak4pfjets_eta()->at(j),			                
			     b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		    if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
							+ 0.6321*(b.year()==2016)
							+ 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		      nBinFat++;                                                      
		    }
		  }
		}
		if (0==nBinFat) {                                                     // Consider FatJets with 0 b-tags
		  nHiggsFatJet250_0b++;
		}
	      }
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_0b > 1) printf("Higgs b0 fjerr: %d\n", nHiggsFatJet250_0b);
      return nHiggsFatJet250_0b;
    });
  
  // Higgs-tagged nFatJets (with gen H) with pt>250 and 1 b jet
  const NamedFunc nGenHHiggsFatJet250b1("nGenHHiggsFatJet250b1",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_1b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {        
	    if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								    //      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								    + 0.8695*(b.year()==2017)
								    + 0.8365*(b.year()==2018))) {
	      if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
		//	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
		int nBinFat = 0;					                
		for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		  //	    if (deltaR(b.FatJet_eta()->at(i),			                
		  //		       b.FatJet_phi()->at(i),			                
		  if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			     b.FatJet_phi()->at(highest_pt_idx),			                
			     b.ak4pfjets_eta()->at(j),			                
			     b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		    if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
							+ 0.6321*(b.year()==2016)
							+ 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		      nBinFat++;                                                      
		    }
		  }
		}
		if (1==nBinFat) {                                                     // Consider FatJets with 1 b-tag
		  nHiggsFatJet250_1b++;
		}
	      }
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_1b > 1) printf("Higgs b1 fjerr: %d\n", nHiggsFatJet250_1b);
      return nHiggsFatJet250_1b;
    });
  
  // Higgs-tagged nFatJets (with Gen H) with pt>250 and 2 b jets
  const NamedFunc nGenHHiggsFatJet250b2("nGenHHiggsFatJet250b2",[](const Baby &b) -> NamedFunc::ScalarType{
      int nHiggsFatJet250_2b = 0;
      unsigned highest_pt_idx = 0;
      float old_pt = 0;      
      if (b.FatJet_pt_nom()->size() > 0) {      
	old_pt = b.FatJet_pt_nom()->at(0);
      }
      float new_pt = 0;
      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                // iterate all FatJets
	new_pt = b.FatJet_pt_nom()->at(i);
	if (new_pt > old_pt) {                                                  // find highest pt index
	  old_pt = new_pt;
	  highest_pt_idx = i;
	}
      }
      if (0 != highest_pt_idx) {
	printf("HIGHEST_PT_IDX NOT 0\n");
      }
      unsigned higgsIdx = 0;
      int higgsIdxCount = 0;
      for (unsigned i = 0; i < b.gen_id()->size(); i++) {
	if (b.gen_id()->at(i) == 25) {
	  higgsIdx = i;
	  higgsIdxCount++;
	}
      }
      if (1 < higgsIdxCount) {
	printf("More than one Higgs in baby: %d! Check!\n", higgsIdxCount);
      }
      if ((b.FatJet_eta()->size() > 0) && (b.FatJet_phi()->size() > 0) && (b.gen_eta()->size() > 0) && (b.gen_phi()->size() > 0)) {
	if (deltaR(b.FatJet_eta()->at(highest_pt_idx), b.FatJet_phi()->at(highest_pt_idx), b.gen_eta()->at(higgsIdx), b.gen_phi()->at(higgsIdx)) < 0.8) {
	  //      for (unsigned i = 0; i < b.FatJet_pt_nom()->size(); i++) {                  // iterate all FatJets
	  if (b.FatJet_pt_nom()->size() > 0) {      
	    if(b.FatJet_deepTagMD_HbbvsQCD()->at(highest_pt_idx) > (0.8945*(b.year()==2016)
								    //      if(b.FatJet_deepTagMD_HbbvsQCD()->at(i) > (0.8945*(b.year()==2016)
								    + 0.8695*(b.year()==2017)
								    + 0.8365*(b.year()==2018))) {
	      if (b.FatJet_pt_nom()->at(highest_pt_idx) > 250) {                                   // consider only > 250 GeV pt ones
		//	if (b.FatJet_pt_nom()->at(i) > 250) {                                   // consider only > 250 GeV pt ones
		int nBinFat = 0;					                
		for (unsigned j = 0; j < b.ak4pfjets_eta()->size(); j++) {            // iterate all Jets
		  //	    if (deltaR(b.FatJet_eta()->at(i),			                
		  //		       b.FatJet_phi()->at(i),			                
		  if (deltaR(b.FatJet_eta()->at(highest_pt_idx),			                
			     b.FatJet_phi()->at(highest_pt_idx),			                
			     b.ak4pfjets_eta()->at(j),			                
			     b.ak4pfjets_phi()->at(j)) < 0.8) {                       // match Jet and FatJet
		    if (b.ak4pfjets_deepCSV()->at(j) > (0.4941*(b.year()==2017)
							+ 0.6321*(b.year()==2016)
							+ 0.4184*(b.year()==2018))) { // check if Jet is b-tagged
		      nBinFat++;                                                      
		    }
		  }
		}
		if (2==nBinFat) {                                                     // Consider FatJets with 2 b-tags
		  nHiggsFatJet250_2b++;
		}
	      }
	    }
	  }
	}
      }
      //    } 
      if(nHiggsFatJet250_2b > 1) printf("Higgs b2 fjerr: %d\n", nHiggsFatJet250_2b);
      return nHiggsFatJet250_2b;
    });

  // returns highest fat jet pt
  const NamedFunc highest_FatJet_pt("highest_FatJet_pt",[](const Baby &b) -> NamedFunc::ScalarType{
      printf("in here\n");
      float highest_pt = 0;
      float new_highest_pt = 0;
      for (unsigned i = 0; b.FatJet_pt_nom()->size(); i++) {
	new_highest_pt = b.FatJet_pt_nom()->at(i);
	if (new_highest_pt > highest_pt) {
	  highest_pt = new_highest_pt;
	}
      }
      printf("highest pt is %f\n", highest_pt);
      return highest_pt;
    });

  //%%%%%%%%%%%%%%%%%%%%%%%%%

  // first attempt for boosted higgs part
  const NamedFunc nBoostedFatJet("nBoostedFatJet",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_pt()->at(i) > 200) nloose++;
      }
      return nloose;
    });

  // first attempt for boosted higgs part
  const NamedFunc HasBoostedHiggs("HasBoostedHiggs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      int nmedium=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.60 && b.ak8pfjets_pt()->at(i) > 200) nloose++;
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.80 && b.ak8pfjets_pt()->at(i) > 200) nmedium++;
      }
      if(nloose>=1 && nmedium>=1) return 1;
      else return 0;
    });

  // first attempt for boosted higgs part
  const NamedFunc HasOnHiggsJet("HasOnHiggsJet",[](const Baby &b) -> NamedFunc::ScalarType{
      int nmedium=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_m()->at(i) > 90 && b.ak8pfjets_m()->at(i) < 150 && b.ak8pfjets_pt()->at(i) > 200) nmedium++;
      }
      if(nmedium>=1) return 1;
      else return 0;
    });

  // boosted AK8 jet, higgs-tagged, with now lepton overlap
  const NamedFunc HasBoostedHiggsNoLepton("HasBoostedHiggsNoLepton",[](const Baby &b) -> NamedFunc::ScalarType{
      int nmedium=0;
      bool lepsInFatJet=false;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
        lepsInFatJet=false;
        for (unsigned j(0); j<b.leps_pt()->size(); j++){
            if (deltaR(b.ak8pfjets_eta()->at(i),b.ak8pfjets_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.8){
                lepsInFatJet=true;
            }
        }
        if (lepsInFatJet) { continue; }
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.80 && b.ak8pfjets_pt()->at(i) > 200) nmedium++;
      }
      if(nmedium>=1) return 1;
      else return 0;
    });

  const NamedFunc HasReallyBoostedHiggs("HasReallyBoostedHiggs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nmedium=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.90 && b.ak8pfjets_pt()->at(i) > 400) nmedium++;
      }
      if(nmedium>=1) return 1;
      else return 0;
    });

  // first attempt for boosted higgs part
  const NamedFunc nVLooseHiggsTags("nVLooseHiggsTags",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.50) nloose++;
      }
      return nloose;
    });

  // first attempt for boosted higgs part
  const NamedFunc nLooseHiggsTags("nLooseHiggsTags",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.80) nloose++;
      }
      return nloose;
    });

  const NamedFunc HasLooseBoostedHiggs("HasLooseBoostedHiggs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      int nmedium=0;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.80) nloose++;
      	if (b.ak8pfjets_deepdisc_hbb()->at(i) > 0.90) nmedium++;
      }
      if(nloose>=1 && nmedium>=0) return 1;
      else return 0;
    });

    const NamedFunc LeadingToppT("LeadingToppT",[](const Baby &b) -> NamedFunc::ScalarType{
    float top_pt=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 6 && b.gen_pt()->at(i) > top_pt ) top_pt = b.gen_pt()->at(i);
    }
    return top_pt;
    });


    const NamedFunc LeadingWpT("LeadingWpT",[](const Baby &b) -> NamedFunc::ScalarType{
    float w_pt=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 24 && b.gen_pt()->at(i) > w_pt ) w_pt = b.gen_pt()->at(i);
    }
    return w_pt;
    });

    const NamedFunc SubLeadingWpT("SubLeadingWpT",[](const Baby &b) -> NamedFunc::ScalarType{
    float w_pt0=0;
    float w_pt1=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
      if ( abs(b.gen_id()->at(i)) == 24 && b.gen_pt()->at(i) > w_pt0 ) {w_pt1= w_pt0;  w_pt0= b.gen_pt()->at(i);}
      else if ( abs(b.gen_id()->at(i)) == 24 && b.gen_pt()->at(i) > w_pt1 ) {w_pt1= b.gen_pt()->at(i);}
    }
    return w_pt1;
    });


    const NamedFunc HasHadronicTau("HasHadronicTau",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTau=0;
    int nLepFromTau=0;
      for (unsigned i(0); i<b.gen_pt()->size(); i++){
       if ( abs(b.gen_id()->at(i)) == 15 && abs(b.gen_motherid()->at(i)) == 24) nTau++;
       if ( (abs(b.gen_id()->at(i)) == 13 || abs(b.gen_id()->at(i)) == 11 ) && abs(b.gen_motherid()->at(i)) == 15) nLepFromTau++;
    }
    if(nTau - nLepFromTau > 0) return 1.;
    else return 0.;
    });


    const NamedFunc n_true_emu("n_true_emu",[](const Baby &b) -> NamedFunc::ScalarType{
    int nemu = 0;

     for (unsigned i(0); i<b.gen_id()->size(); i++){
        if ((abs(b.gen_id()->at(i))==11 || abs(b.gen_id()->at(i))==13) && (abs(b.gen_motherid()->at(i) == 24)|| abs(b.gen_motherid()->at(i))==15)) nemu++;
    }//Close if
    return nemu;
    });

      //Number of b-flavor jets in acceptance but failing tagger.
      const NamedFunc HasBFailedTag("HasBFailedTag",[](const Baby &b) -> NamedFunc::ScalarType{
      int nfail=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i) < (0.4941*(b.year()==2017) + 0.6321*(b.year()==2016) + 0.4184*(b.year()==2018))) nfail++;
      }
      return nfail;

    });


   /*
   * returns 1 if there are two b jets with at least one Medium and one Loose score
   * returns 0 otherwise
   */
  const NamedFunc HasMedLooseCSV("HasMedLooseCSV",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      int nmedium=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
      	if (b.ak4pfjets_deepCSV()->at(i) > 0.2217) nloose++;
      	if (b.ak4pfjets_deepCSV()->at(i) > 0.6321) nmedium++;
      }
      if(nloose>=2 && nmedium>=1) return 1;
      else return 0;
    });

  const NamedFunc HasMedMedDeepCSV("HasMedMedDeepCSV",[](const Baby &b) -> NamedFunc::ScalarType{
      int ndeepmedium=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.6324) ndeepmedium++;
      }
      if(ndeepmedium>=2) return 1;
      else return 0;
    });

  const NamedFunc HasLooseLooseDeepCSV("HasLooseLooseDeepCSV",[](const Baby &b) -> NamedFunc::ScalarType{
      int ndeeploose=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2219) ndeeploose++;
      }
      if(ndeeploose>=2) return 1;
      else return 0;
    });

  const NamedFunc HasMedLooseDeepCSV("HasMedLooseDeepCSV",[](const Baby &b) -> NamedFunc::ScalarType{
      int ndeeploose=0;
      int ndeepmedium=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2219) ndeeploose++;
        if (b.ak4pfjets_deepCSV()->at(i) > 0.6324) ndeepmedium++;
      }
      if(ndeeploose>=2 && ndeepmedium>=1) return 1;
      else return 0;
    });

  const NamedFunc HasLooseNoMedDeepCSV("HasLooseNoMedDeepCSV",[](const Baby &b) -> NamedFunc::ScalarType{
      int ndeeploose=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if ((b.ak4pfjets_deepCSV()->at(i) > 0.2219)&&(b.ak4pfjets_deepCSV()->at(i) < 0.6324)) ndeeploose++;
      }
      if(ndeeploose>=2) return 1;
      else return 0;
    });

  const NamedFunc nDeepMedBTagged("nDeepMedBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.6324) njets++;
      }
      return njets;
    });

  const NamedFunc nDeepLooseBTagged("nDeepLooseBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2219) njets++;
      }
      return njets;
    });

  const NamedFunc nDeepMedCTagged("nDeepMedCTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSVc()->size(); i++){
        if (b.ak4pfjets_deepCSVc()->at(i)/(b.ak4pfjets_deepCSVc()->at(i)+b.ak4pfjets_deepCSVl()->at(i)) > 0.155) njets++;
      }
      return njets;
    });
  const NamedFunc nDeepTightCTagged("nDeepTightCTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSVc()->size(); i++){
        if (b.ak4pfjets_deepCSVc()->at(i)/(b.ak4pfjets_deepCSVc()->at(i)+b.ak4pfjets_deepCSVl()->at(i)) > 0.59) njets++;
      }
      return njets;
    });

  const NamedFunc nDeepTightCvBTagged("nDeepTightCvBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSVc()->size(); i++){
        if (b.ak4pfjets_deepCSVc()->at(i)/(b.ak4pfjets_deepCSVc()->at(i)+b.ak4pfjets_deepCSV()->at(i)) > 0.19) njets++;
      }
      return njets;
    });

  const NamedFunc nDeepMedCvBTagged("nDeepMedCvBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSVc()->size(); i++){
        if (b.ak4pfjets_deepCSVc()->at(i)/(b.ak4pfjets_deepCSVc()->at(i)+b.ak4pfjets_deepCSV()->at(i)) > 0.14) njets++;
      }
      return njets;
    });
  const NamedFunc nDeepLooseCvBTagged("nDeepLooseCvBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak4pfjets_deepCSVc()->size(); i++){
        if (b.ak4pfjets_deepCSVc()->at(i)/(b.ak4pfjets_deepCSVc()->at(i)+b.ak4pfjets_deepCSV()->at(i)) > 0.05) njets++;
      }
      return njets;
    });

  const NamedFunc nAK8jets("nAK8jets",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak8pfjets_pt()->size(); i++){
        njets++;
      }
      return njets;
    });
  const NamedFunc nHiggsTag("nHiggsTag",[](const Baby &b) -> NamedFunc::ScalarType{
      int njets=0;
      for(unsigned i(0); i<b.ak8pfjets_pt()->size(); i++){
        if(b.ak8pfjets_deepdisc_hbb()->at(i)>0.8)njets++;
      }
      return njets;
    });


  const NamedFunc bJetPt("bJetPt",[](const Baby &b) -> NamedFunc::VectorType{
      vector<double> bjetpt;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5) bjetpt.push_back(b.ak4pfjets_pt()->at(i));
      }
      return bjetpt;
    });

  /*
  * Returns number of >loose bjets
  */
  const NamedFunc NBJets("NBJets",[](const Baby &b) -> NamedFunc::ScalarType{
      int nloose=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2217) nloose++;
      }
      return nloose;
    });

  /*
   * returns number of good, isolated leptons in event
   */
  const NamedFunc WHLeptons("WHLeptons",[](const Baby &b) -> NamedFunc::ScalarType{
      int nwhleptons=0;
      if (abs(b.lep1_pdgid())==11&&b.leps_pt()->at(0)>30&&b.lep1_relIso()*b.leps_pt()->at(0)<5) nwhleptons++;
      if (abs(b.lep1_pdgid())==13&&b.leps_pt()->at(0)>25&&b.lep1_relIso()*b.leps_pt()->at(0)<5&&abs(b.leps_eta()->at(0))<2.1) nwhleptons++;

      if (b.leps_pt()->size()==2){
	       if (abs(b.lep2_pdgid())==11&&b.leps_pt()->at(1)>30&&b.lep2_relIso()*b.leps_pt()->at(1)<5) nwhleptons++;
	       if (abs(b.lep2_pdgid())==13&&b.leps_pt()->at(1)>25&&b.lep2_relIso()*b.leps_pt()->at(1)<5&&abs(b.leps_eta()->at(1))<2.1) nwhleptons++;
      }
      return nwhleptons;
    });

  const NamedFunc WHMuonSigEff("WHMuonSigEff",[](const Baby &b) -> NamedFunc::VectorType{
      double nMuonsBarrel = 0;
      double nMuonsEndcap = 0;
      double nMuons = 0;

      vector< double > muonEtaCountsVector;
      
      if (abs(b.lep1_pdgid())==13&&b.leps_pt()->at(0)>25&&b.lep1_relIso()*b.leps_pt()->at(0)<5){
        nMuons++;
        if(abs(b.leps_eta()->at(0))<=1.442) nMuonsBarrel++;
        if(abs(b.leps_eta()->at(0))>1.442) nMuonsEndcap++;
      }

      if (b.leps_pt()->size()==2){
	       //if (abs(b.lep2_pdgid())==11&&b.leps_pt()->at(1)>30&&b.lep2_relIso()*b.leps_pt()->at(1)<5) nwhleptons++;
	       //if (abs(b.lep2_pdgid())==13&&b.leps_pt()->at(1)>25&&b.lep2_relIso()*b.leps_pt()->at(1)<5&&abs(b.leps_eta()->at(1))<2.1) nwhleptons++;
	       
         if (abs(b.lep2_pdgid())==13&&b.leps_pt()->at(1)>25&&b.lep2_relIso()*b.leps_pt()->at(1)<5){
           nMuons++;
           if(abs(b.leps_eta()->at(0))<=1.442) nMuonsBarrel++;
           if(abs(b.leps_eta()->at(0))>1.442) nMuonsEndcap++;
         }
      }

      muonEtaCountsVector.push_back(nMuons);
      muonEtaCountsVector.push_back(nMuonsBarrel);
      muonEtaCountsVector.push_back(nMuonsEndcap);
      
      return muonEtaCountsVector;
    });

  const NamedFunc WHElSigEff("WHElSigEff",[](const Baby &b) -> NamedFunc::VectorType{
      double nElesBarrel = 0;
      double nElesEndcap = 0;
      double nEles = 0;

      vector< double > EleEtaCountsVector;
      
      if (abs(b.lep1_pdgid())==11&&b.leps_pt()->at(0)>30&&b.lep1_relIso()*b.leps_pt()->at(0)<5){
        nEles++;
        if(abs(b.leps_eta()->at(0))<=1.442) nElesBarrel++;
        if(abs(b.leps_eta()->at(0))>1.442) nElesEndcap++;
      }

      if (b.leps_pt()->size()==2){
	       //if (abs(b.lep2_pdgid())==11&&b.leps_pt()->at(1)>30&&b.lep2_relIso()*b.leps_pt()->at(1)<5) nwhleptons++;
	       //if (abs(b.lep2_pdgid())==13&&b.leps_pt()->at(1)>25&&b.lep2_relIso()*b.leps_pt()->at(1)<5&&abs(b.leps_eta()->at(1))<2.1) nwhleptons++;
	       
         if (abs(b.lep2_pdgid())==11&&b.leps_pt()->at(1)>30&&b.lep2_relIso()*b.leps_pt()->at(1)<5){
           nEles++;
           if(abs(b.leps_eta()->at(0))<=1.442) nElesBarrel++;
           if(abs(b.leps_eta()->at(0))>1.442) nElesEndcap++;
         }
      }

      EleEtaCountsVector.push_back(nEles);
      EleEtaCountsVector.push_back(nElesBarrel);
      EleEtaCountsVector.push_back(nElesEndcap);
      
      return EleEtaCountsVector;
    });

  const NamedFunc nRealBs("nRealBs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5) nbquarks++;
      }
      return nbquarks;
    }); 

    const NamedFunc nRealBsfromTop("nRealBsfromTop",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
          float dr=999;
          int matchindex=-1;
          for(unsigned j(0); j<b.gen_id()->size(); j++){
            if(abs(b.gen_id()->at(j)) == 5){
              float this_dr = deltaR(b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
              if( this_dr < dr ){ dr = this_dr; matchindex=j;}
            }
        }
        if(dr < 0.4 && matchindex>=0 && abs(b.gen_motherid()->at(matchindex))==6) nbquarks++;
      }//found b-flavor jet
      }
      return nbquarks;
    });


    const NamedFunc nRealBtags("nRealBtags",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))) nbquarks++;
      }
      return nbquarks;
    });

    const NamedFunc max_genjet_bquark_pt_ratio("max_genjet_bquark_pt_ratio",[](const Baby &b) -> NamedFunc::ScalarType{
      float max_ratio=-1;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
          float dr=999;
          float this_ratio=-1;
          int matchindex=-1;
          for(unsigned j(0); j<b.gen_id()->size(); j++){
            if(abs(b.gen_id()->at(j)) == 5){
              float this_dr = deltaR(b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
              if( this_dr < dr ){ dr = this_dr; matchindex=j; this_ratio = (b.ak4pfjets_genpt()->at(i) / b.gen_pt()->at(j));}
            }
         }
        if(dr < 0.4 && matchindex>=0 && this_ratio > max_ratio) max_ratio = this_ratio;
      }//found b-flavor jet
      }
      return max_ratio;
    });

    const NamedFunc dR_jet_bquark_max_pt_ratio("dR_genjet_bquark_max_pt_ratio",[](const Baby &b) -> NamedFunc::ScalarType{
      float max_ratio=-1;
      float dr_max_ratio=-1;
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i)>(medDeepCSV2017*(b.year()==2017) + medDeepCSV2016*(b.year()==2016) + medDeepCSV2018*(b.year()==2018))){ //Then find closest gen b
          float dr=999;
          float this_ratio=-1;
          int matchindex=-1;
          for(unsigned j(0); j<b.gen_id()->size(); j++){
            if(abs(b.gen_id()->at(j)) == 5){
              float this_dr = deltaR(b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
              if( this_dr < dr ){ dr = this_dr; matchindex=j; this_ratio = (b.ak4pfjets_genpt()->at(i) / b.gen_pt()->at(j));}
            }
         }
        if(dr < 0.4 && matchindex>=0 && this_ratio > max_ratio){max_ratio = this_ratio; dr_max_ratio = dr;}
        }//found b-flavor jet
      }
      return dr_max_ratio;
    });




  const NamedFunc nGenLightLeps("nGenLightLeps",[](const Baby &b) -> NamedFunc::ScalarType{
      int nGenLeps=0;

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          nGenLeps++;
        }
      }

      return nGenLeps;
      //if function returns 2, 1 is lost lepton, assuming single reco lepton
    });

  const NamedFunc LostHadTaus("LostHadTaus",[](const Baby &b) -> NamedFunc::ScalarType{
      int nGenLepsfromW=0;
      int nGenLepsfromTau=0;

      int lostLepChecker = 0;

      int counter=0; 


      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&abs(b.gen_motherid()->at(i))==24){
          lostLepChecker++;
        }
      }

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&abs(b.gen_motherid()->at(i))==24){
          nGenLepsfromW++;
        }
      }

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&abs(b.gen_motherid()->at(i))==15){
          nGenLepsfromTau++;
        }
      }

      if(lostLepChecker==2){
        if((nGenLepsfromW+nGenLepsfromTau)<2){
          counter = 1;
          //if light leptons from W or tau is less than 2, then the lost lepton is a hadronic tau
        }else{
          counter = 0;
          //if there are two light leptons from either a W or a tau, then then lost lepton is just a light lepton and can be counted using the nGenLightLeps function
        }
      }else{
        counter = 0;
      }

      return counter;
    });

  const NamedFunc ptLostLeps("ptLostLeps",[](const Baby &b) -> NamedFunc::ScalarType{
      float ptLostLep=0.;

      int lostLepChecker = 0;

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&abs(b.gen_motherid()->at(i))==24){
          lostLepChecker++;
        }
      }

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        //check for lost light leptons first
        if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          //if we found an event with two light leptons, check for gen matching
          for(unsigned j(0); j<b.leps_eta()->size(); j++){
            if(deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.4){
              //this lepton is matched! so we want to continue the for loop
              continue;
            }else{
              //this lepton is the lost one! so we want to save the pt
              ptLostLep=b.gen_pt()->at(i);
            }//closes if/else statement
          }//closes for loop over reco leps
        }else if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==15)&&(abs(b.gen_motherid()->at(i))==24)){
          //this is a lost hadronic tau, so we just automatically save the pt of the tau
          ptLostLep=b.gen_pt()->at(i);
        }//closes if statement
      }//closes first for loop

      return ptLostLep;
      //returns pt of least energetic light lepton in gen-level event
    });

  const NamedFunc etaLostLeps("etaLostLeps",[](const Baby &b) -> NamedFunc::ScalarType{
      float etaLostLep=0.;

      int lostLepChecker = 0;

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&abs(b.gen_motherid()->at(i))==24){
          lostLepChecker++;
        }
      }

      for(unsigned i(0); i<b.gen_id()->size(); i++){
      //check for lost light leptons first
        if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          //if we found an event with two light leptons, check for gen matching
          for(unsigned j(0); j<b.leps_eta()->size(); j++){
            if(deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.4){
              //this lepton is matched! so we want to continue the for loop
              continue;
            }else{
              //this lepton is the lost one! so we want to save the eta
              etaLostLep=b.gen_eta()->at(i);
            }//closes if/else statement
          }//closes for loop over reco leps
        }else if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==15)&&(abs(b.gen_motherid()->at(i))==24)){
          //this is a lost hadronic tau, so we just automatically save the eta of the tau
          etaLostLep=b.gen_eta()->at(i);
        }//closes if statement
      }//closes first for loop

      return etaLostLep;
      //returns eta of least energetic light lepton in gen-level event
    });

  const NamedFunc causeLostLeps("causeLostLeps",[](const Baby &b) -> NamedFunc::VectorType{

      int causeVar = 0;
        //causeVar==1 for low pT
        //causeVar==2 for high |eta|
        //causeVar==0 for other
      int lostLepID = 0;
      int lostLepChecker = 0;

      vector<double> causeVec;

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&abs(b.gen_motherid()->at(i))==24){
          lostLepChecker++;
        }
      }

      //now we will check if lost due to being too forward
      for(unsigned i(0); i<b.gen_id()->size(); i++){
      //check for lost light leptons first
        if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          //if we found an event with two light leptons, check for gen matching
          for(unsigned j(0); j<b.leps_eta()->size(); j++){
            if(deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.4){
              //this lepton is matched! so we want to continue the for loop
              continue;
            }else{
              //this lepton is the lost one! so we want to check the eta threshold
              if((b.gen_id()->at(i)==11&&abs(b.gen_eta()->at(i))>2.4)||(b.gen_id()->at(i)==11&&(abs(b.gen_eta()->at(i))>1.4442&&abs(b.gen_eta()->at(i))<1.566))||(b.gen_id()->at(i)==13&&abs(b.gen_eta()->at(i))>2.4)){
                causeVar = 2;
                lostLepID = abs(b.gen_id()->at(i));
              }
            }//closes if/else statement
          }//closes for loop over reco leps
        }else if(lostLepChecker==2&&(abs(b.gen_id()->at(i))==15)&&(abs(b.gen_motherid()->at(i))==24)){
          //this is a lost hadronic tau, so we check the eta threshold
          if(abs(b.gen_eta()->at(i))>2.4){
                causeVar = 2;
                lostLepID = abs(b.gen_id()->at(i));
          }
        }//closes if statement
      }//closes first for loop

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        //check for lost light leptons first
        if(lostLepChecker==2&&causeVar!=2&&(abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          //if we found an event with two light leptons, check for gen matching
          for(unsigned j(0); j<b.leps_eta()->size(); j++){
            if(deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.4){
              //this lepton is matched! so we want to continue the for loop
              continue;
            }else{
              //this lepton is the lost one! so we need to enforce the pT threshold to check if that's why it's lost
              if((b.gen_id()->at(i)==11&&b.gen_pt()->at(i)<5)||(b.gen_id()->at(i)==13&&b.gen_pt()->at(i)<5)){
                causeVar = 1;
                lostLepID = abs(b.gen_id()->at(i));
              }
            }//closes if/else statement
          }//closes for loop over reco leps
        }else if(lostLepChecker==2&&causeVar!=2&&(abs(b.gen_id()->at(i))==15)&&(abs(b.gen_motherid()->at(i))==24)){
          //this is a lost hadronic tau, so we check the pT of the tau
          if(b.gen_pt()->at(i)<20){
            causeVar = 1;
            lostLepID = abs(b.gen_id()->at(i));
          }
        }//closes if statement
      }//closes first for loop

      for(unsigned i(0); i<b.gen_id()->size(); i++){
        //check for lost light leptons first
        if(lostLepChecker==2&&causeVar!=2&&causeVar!=1&&(abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13)&&(abs(b.gen_motherid()->at(i))==24||abs(b.gen_motherid()->at(i))==15)){
          //if we found an event with two light leptons, check for gen matching
          for(unsigned j(0); j<b.leps_eta()->size(); j++){
            if(deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j))<0.4){
              //this lepton is matched! so we want to continue the for loop
              continue;
            }else{
              //this lepton is the lost one!
              causeVar = 0;
              lostLepID = abs(b.gen_id()->at(i));
            }//closes if/else statement
          }//closes for loop over reco leps
        }else if(lostLepChecker==2&&causeVar!=2&&causeVar!=1&&(abs(b.gen_id()->at(i))==15)&&(abs(b.gen_motherid()->at(i))==24)){
          //this is a lost hadronic tau
          causeVar = 0;
          lostLepID = abs(b.gen_id()->at(i));
        }//closes if statement
      }//closes first for loop

      causeVec.push_back(lostLepID);
      causeVec.push_back(causeVar);

      return causeVec;
      //returns cause of losing lepton in gen-level event
    });


  const NamedFunc FatJet_HighestHScore("FatJet_HighestHScore",[](const Baby &b) -> NamedFunc::ScalarType{
    //vector<float>* v = b.FatJet_deepTag_H();
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200){
            v.push_back(b.FatJet_deepTag_H()->at(i));
        }
    }
    sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc FatJet_HighestMDHScore("FatJet_HighestMDHScore",[](const Baby &b) -> NamedFunc::ScalarType{
    //vector<float>* v = b.FatJet_deepTagMD_HbbvsQCD();
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200){
            v.push_back(b.FatJet_deepTagMD_HbbvsQCD()->at(i));
        }
    }
    sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc MediumMDHiggsTag_m("MediumMDHiggsTag_m",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_mass()->at(i));
        }
    }
    //sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc MediumMDHiggsTag_msoftdrop("MediumMDHiggsTag_msoftdrop",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_msoftdrop()->at(i));
        }
    }
    //sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc MediumMDHiggsTag_mass("MediumMDHiggsTag_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_mass()->at(i));
        }
    }
    //sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc LooseMDHiggsTag_msoftdrop("LooseMDHiggsTag_msoftdrop",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.5165){
            v.push_back(b.FatJet_msoftdrop()->at(i));
        }
    }
    //sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });



  const NamedFunc MediumMDHiggsTag_pt("MediumMDHiggsTag_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_pt()->at(i));
        }
    }
    sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc MediumMDHiggsTag_MDHscore("MediumMDHiggsTag_MDHscore",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_deepTagMD_HbbvsQCD()->at(i));
        }
    }
    sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc MediumMDHiggsTag_Hscore("MediumMDHiggsTag_Hscore",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            v.push_back(b.FatJet_deepTag_H()->at(i));
        }
    }
    sort(v.begin(), v.end(), greater<int>());

    if (v.size()>0){
        return v.at(0);
    }
    else {
        return 0;
    }
  });

  const NamedFunc nMediumMDHiggsTag("nMediumMDHiggsTag",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            nTags++;
        }
    }
    return nTags;
  });

  const NamedFunc HasMediumMDHiggsTag("HasMediumMDHiggsTag",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            nTags++;
        }
    }
    return nTags>0;
  });

  const NamedFunc HasMediumMDHiggsTag_medPt("HasMediumMDHiggsTag_medPt",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>250 && b.FatJet_deepTagMD_HbbvsQCD()->at(i)>0.8365){
            nTags++;
        }
    }
    return nTags>0;
  });

  const NamedFunc HasMediumMDBBTag_medPt("HasMediumMDBBTag_medPt",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>400 && b.FatJet_deepTagMD_bbvsLight()->at(i)>0.8365){
            nTags++;
        }
    }
    return nTags>0;
  });


  const NamedFunc HasHbbTag("HasHbbTag",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_btagHbb()->at(i)>0.7){
            nTags++;
        }
    }
    return nTags>0;
  });

  const NamedFunc nMediumHiggsTag("nMediumHiggsTag",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTag_H()->at(i)>0.8){
            nTags++;
        }
    }
    return nTags;
  });

  const NamedFunc HasMediumHiggsTag("HasMediumHiggsTag",[](const Baby &b) -> NamedFunc::ScalarType{
    int nTags = 0;
    vector<float> v;
    for(unsigned i(0); i<b.FatJet_pt()->size(); i++){
        if(b.FatJet_pt()->at(i)>200 && b.FatJet_deepTag_H()->at(i)>0.8){
            nTags++;
        }
    }
    return nTags>0;
  });


  const NamedFunc FatJet_ClosestMSD("FatJet_ClosestMSD",[](const Baby &b) -> NamedFunc::ScalarType{
    //vector<float>* v = b.FatJet_msoftdrop();
    vector<float> masses;

    for(unsigned i(0); i<b.FatJet_msoftdrop()->size(); i++){
        masses.push_back(abs(b.FatJet_msoftdrop()->at(i) - 125.));
        //if(abs(b.ak4pfjets_parton_flavor()->at(i))==5) nbquarks++;
    }

    sort(masses.begin(), masses.end(), less<int>());

    return (masses.at(0)+125.);
  });

  const NamedFunc sortedJetsPt_Leading("sortedJetsPt_Leading",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float>* v = b.ak4pfjets_pt();

    sort(v->begin(), v->end(), greater<int>());

    return v->at(0);
  });

  const NamedFunc sortedJetsPt_subLeading("sortedJetsPt_subLeading",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<float>* v = b.ak4pfjets_pt();

    sort(v->begin(), v->end(), greater<int>());

    return v->at(1);
  });

  const NamedFunc sortedJetsCSV_Leading("sortedJetsCSV_Leading",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    for(unsigned i(0);i<b.ak4pfjets_deepCSV()->size();i++){
      v.push_back(make_pair(b.ak4pfjets_deepCSV()->at(i),b.ak4pfjets_pt()->at(i)));
    }

    sort(v.begin(),v.end());

    int n = v.size();

    return v[n-1].second;
  });

  const NamedFunc sortedJetsCSV_subLeading("sortedJetsCSV_subLeading",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    for(unsigned i(0);i<b.ak4pfjets_deepCSV()->size();i++){
      v.push_back(make_pair(b.ak4pfjets_deepCSV()->at(i),b.ak4pfjets_pt()->at(i)));
    }

    sort(v.begin(),v.end());

    int n = v.size();

    return v[n-2].second;
  });

  const NamedFunc sortedJetsCSV_deltaR("sortedJetsCSV_deltaR",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<pair <float,float> > v_etaPhi;
    vector<pair <float,pair<float,float> > > v_ptEtaPhi;
    vector<pair <float, pair<float, pair <float,float> > > > v_csvPtEtaPhi;

    for(unsigned i(0);i<b.ak4pfjets_deepCSV()->size();i++){
      v_etaPhi.push_back(make_pair(b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i)));
      v_ptEtaPhi.push_back(make_pair(b.ak4pfjets_pt()->at(i),v_etaPhi[i]));
      v_csvPtEtaPhi.push_back(make_pair(b.ak4pfjets_deepCSV()->at(i),v_ptEtaPhi[i]));
    }
    sort(v_csvPtEtaPhi.begin(), v_csvPtEtaPhi.end());

    int n = v_csvPtEtaPhi.size();
    
    float deltaR_leading = deltaR(v_csvPtEtaPhi[n-1].second.second.first,v_csvPtEtaPhi[n-1].second.second.second,v_csvPtEtaPhi[n-2].second.second.first,v_csvPtEtaPhi[n-2].second.second.second);

    return deltaR_leading;
  });


  const NamedFunc deltaRLeadingJets("deltaRLeadingJets",[](const Baby &b) -> NamedFunc::ScalarType{
  
    float maxpt=0; int maxindex=-1;
    float secondmaxpt=0; int secondindex=-1;
    for (unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      if (b.ak4pfjets_pt()->at(i) > secondmaxpt && b.ak4pfjets_pt()->at(i)<=maxpt){
        secondmaxpt = b.ak4pfjets_pt()->at(i);
        secondindex=i;
      }
      else if (b.ak4pfjets_pt()->at(i) > maxpt){
        secondmaxpt=maxpt;
        secondindex=maxindex;
        maxpt = b.ak4pfjets_pt()->at(i);
        maxindex=i;
      }
    }
    return deltaR(b.ak4pfjets_eta()->at(maxindex),b.ak4pfjets_phi()->at(maxindex),b.ak4pfjets_eta()->at(secondindex),b.ak4pfjets_phi()->at(secondindex));
  });

    const NamedFunc deltaPhiLeadingJets("deltaPhiLeadingJets",[](const Baby &b) -> NamedFunc::ScalarType{

    float maxpt=0; int maxindex=-1;
    float secondmaxpt=0; int secondindex=-1;
    for (unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      if (b.ak4pfjets_pt()->at(i) > secondmaxpt && b.ak4pfjets_pt()->at(i)<=maxpt){
        secondmaxpt = b.ak4pfjets_pt()->at(i);
        secondindex=i;
      }
      else if (b.ak4pfjets_pt()->at(i) > maxpt){
        secondmaxpt=maxpt;
        secondindex=maxindex;
        maxpt = b.ak4pfjets_pt()->at(i);
        maxindex=i;
      }
    }
    return deltaPhi(b.ak4pfjets_phi()->at(maxindex),b.ak4pfjets_phi()->at(secondindex));
  });

  const NamedFunc bbmass("bbmass",[](const Baby &b) -> NamedFunc::ScalarType{
    float mass = -1;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)){
            //delR = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
            TLorentzVector v1,v2,sum;
            v1.SetPtEtaPhiM(b.gen_pt()->at(i),b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_m()->at(i));
            v2.SetPtEtaPhiM(b.gen_pt()->at(j),b.gen_eta()->at(j),b.gen_phi()->at(j),b.gen_m()->at(j));
            sum=v1+v2;
            mass = sum.M();
          }//Close if statement for opposite b
        }//Close for loop over second half of particles in event
        if (mass>-1) break;
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return mass;

    });

  const NamedFunc hasGenBs("hasGenBs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5 && (abs(b.gen_motherid()->at(i))!=23) && (abs(b.gen_motherid()->at(i))!=24) && (abs(b.gen_motherid()->at(i))!=6) ) nbquarks++; // only care about b's from gluons/protons.
      }
      if (nbquarks>0) return 1;
      return 0;
    });

  const NamedFunc nGenBs("nGenBs",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc nGenBsFromGluons("nGenBsFromGluons",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5&&abs(b.gen_motherid()->at(i))==21) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc nGenBs_ptG15("nGenBs_ptG15",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>15) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc nGenBsFromGluons_ptG15("nGenBsFromGluons_ptG15",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5&&abs(b.gen_motherid()->at(i))==21&&b.gen_pt()->at(i)>15) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc nGenBs_ptG30("nGenBs_ptG30",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>30) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc nGenBs_forward("nGenBs_forward",[](const Baby &b) -> NamedFunc::ScalarType{
    int nbquarks=0;
    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&abs(b.gen_eta()->at(i))>2.4) nbquarks++;
    }
    return nbquarks;
  });

  const NamedFunc nGenBs_central("nGenBs_central",[](const Baby &b) -> NamedFunc::ScalarType{
    int nbquarks=0;
    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&abs(b.gen_eta()->at(i))<2.4) nbquarks++;
    }
    return nbquarks;
  });

  const NamedFunc nGenBs_ptG30_central("nGenBs_ptG30_central",[](const Baby &b) -> NamedFunc::ScalarType{
    int nbquarks=0;
    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5 && abs(b.gen_eta()->at(i))<2.4 && b.gen_pt()->at(i)>30) nbquarks++;
    }
    return nbquarks;
  });



  const NamedFunc nGenBsFromGluons_ptG30("nGenBsFromGluons_ptG30",[](const Baby &b) -> NamedFunc::ScalarType{
      int nbquarks=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5&&abs(b.gen_motherid()->at(i))==21&&b.gen_pt()->at(i)>30) nbquarks++;
      }
      return nbquarks;
    });

  const NamedFunc genBpT("genBpT",[](const Baby &b) -> NamedFunc::VectorType{
      vector <double> bquarkpT;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5){
          bquarkpT.push_back (b.gen_pt()->at(i));
        }
      }
      return bquarkpT;
    });
  
  const NamedFunc genBeta("genBeta",[](const Baby &b) -> NamedFunc::VectorType{
      vector <double> bquarketa;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5){
          bquarketa.push_back (b.gen_eta()->at(i));
        }
      }
      return bquarketa;
    });

    const NamedFunc genBeta_mostForward("genBeta_mostForward",[](const Baby &b) -> NamedFunc::ScalarType{
      double bquarketa=0;
      for(unsigned i(0); i<b.gen_id()->size(); i++){
        if(abs(b.gen_id()->at(i))==5 && abs(b.gen_eta()->at(i))>abs(bquarketa)){
          bquarketa=b.gen_eta()->at(i);
        }
      }
      return bquarketa;
    });


  const NamedFunc genB_leadingpT("genB_leadingpT",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    for(unsigned i(0);i<b.gen_id()->size();i++){
      if(abs(b.gen_id()->at(i))==5){
        v.push_back(make_pair(b.gen_pt()->at(i),b.gen_id()->at(i)));
      }
    }

    sort(v.begin(),v.end());

    int n = v.size();

    float bquarkpT = v[n-1].first;

    return bquarkpT;
  });

  const NamedFunc genB_subleadingpT("genB_subleadingpT",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    for(unsigned i(0);i<b.gen_id()->size();i++){
      if(abs(b.gen_id()->at(i))==5){
        v.push_back(make_pair(b.gen_pt()->at(i),b.gen_id()->at(i)));
      }
    }

    sort(v.begin(),v.end());

    int n = v.size();

    float bquarkpT = v[n-2].first;

    return bquarkpT;
  });

  const NamedFunc bDeltaRGluonSplit("bDeltaRGluonSplit",[](const Baby &b) -> NamedFunc::ScalarType{
    float delR = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_motherid()->at(i)==21){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)&&b.gen_motheridx()->at(j)==b.gen_motheridx()->at(i)&&(b.gen_pt()->at(j)>30||b.gen_pt()->at(i)>30)){
            delR = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
          }//Close if statement for opposite b, same mother index, pt cut
        }//Close for loop over second half of particles in event
      }//Close if statement that find b with gluon mother
    }//Close for loop over all particles in event
    return delR;

    });

  const NamedFunc bDeltaRHiggs("bDeltaRHiggs",[](const Baby &b) -> NamedFunc::ScalarType{
    float delR = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_motherid()->at(i)==25){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)&&b.gen_motheridx()->at(j)==b.gen_motheridx()->at(i)&&(b.gen_pt()->at(j)>30||b.gen_pt()->at(i)>30)){
            delR = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
          }//Close if statement for opposite b, same mother index, pt cut
        }//Close for loop over second half of particles in event
      }//Close if statement that find b with gluon mother
    }//Close for loop over all particles in event
    return delR;

    });

  const NamedFunc bDeltaR("bDeltaR",[](const Baby &b) -> NamedFunc::ScalarType{
    float delR = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)){
            delR = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
          }//Close if statement for opposite b
        }//Close for loop over second half of particles in event
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return delR;

    });

  const NamedFunc bDeltaRfromtop("bDeltaRfromtop",[](const Baby &b) -> NamedFunc::ScalarType{
    float delR = 0;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5 && abs(b.gen_motherid()->at(i))==6){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i) && abs(b.gen_motherid()->at(j))==6){

            delR = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.gen_eta()->at(j),b.gen_phi()->at(j));
          }//Close if statement for opposite b
        }//Close for loop over second half of particles in event
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return delR;

    });

  const NamedFunc minDR_lep_bquark("minDR_lep_bquark",[](const Baby &b) -> NamedFunc::ScalarType{
    float delR = 999;
    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5 && abs(b.gen_motherid()->at(i))==6){
        for(unsigned j(0);j<b.leps_eta()->size();j++){
          float this_dr = deltaR(b.gen_eta()->at(i),b.gen_phi()->at(i),b.leps_eta()->at(j),b.leps_phi()->at(j));
          if (this_dr < delR) delR = this_dr; 
        }//Close for loop over leptons
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return delR;

    });

  const NamedFunc bMother("bMother",[](const Baby &b) -> NamedFunc::VectorType{
    vector <double> mother_id;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5){
        if(abs(b.gen_motherid()->at(i))==2212){
          mother_id.push_back (30);
        }else{
          mother_id.push_back (b.gen_motherid()->at(i));
        }//Close if/else statement for redefining proton
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return mother_id;

    });

  const NamedFunc bMother_pt15("bMother_pt15",[](const Baby &b) -> NamedFunc::VectorType{
    vector <double> mother_id;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>15){
        if(abs(b.gen_motherid()->at(i))==2212){
          mother_id.push_back (30);
        }else{
          mother_id.push_back (b.gen_motherid()->at(i));
        }//Close if/else statement for redefining proton
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return mother_id;

    });

  const NamedFunc bMother_pt30("bMother_pt30",[](const Baby &b) -> NamedFunc::VectorType{
    vector <double> mother_id;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>30){
        if(abs(b.gen_motherid()->at(i))==2212){
          mother_id.push_back (30);
        }else{
          mother_id.push_back (b.gen_motherid()->at(i));
        }//Close if/else statement for redefining proton
      }//Close if statement that find b
    }//Close for loop over all particles in event
    return mother_id;

    });

  const NamedFunc leadingBMother_pt20("leadingBMother_pt20",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    int ngenb = 0;

    for(unsigned i(0);i<b.gen_id()->size();i++){
      if(abs(b.gen_id()->at(i))==5){
        ngenb++;
      }
      if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>20){
        v.push_back(make_pair(b.gen_pt()->at(i),b.gen_motherid()->at(i)));
      }
    }

    sort(v.begin(),v.end());

    int n = v.size();

    float bmom = 0.;
    
    if(n!=0){
      bmom = v[n-1].second;
    }

    if(ngenb==0.){
      bmom = 19;
    }

    return bmom;
  });

  const NamedFunc subleadingBMother_pt20("subleadingBMother_pt20",[](const Baby &b) -> NamedFunc::ScalarType{
    vector< pair <float,float> > v;
    int ngenb = 0;

    for(unsigned i(0);i<b.gen_id()->size();i++){
      if(abs(b.gen_id()->at(i))==5){
        ngenb++;
      }
      if(abs(b.gen_id()->at(i))==5&&b.gen_pt()->at(i)>20){
        v.push_back(make_pair(b.gen_pt()->at(i),b.gen_motherid()->at(i)));
      }
    }

    sort(v.begin(),v.end());

    int n = v.size();

    float bmom = 0.;
    
    if(n>=2){
      bmom = v[n-2].second;
    }

    if(ngenb==0.){
      bmom = 19;
    }

    return bmom;
  });

  const NamedFunc bDeltaPhi("bDeltaPhi",[](const Baby &b) -> NamedFunc::ScalarType{
      float delphi=0;
      vector < pair < float,float > > v_csvPhi;
      
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5){
          v_csvPhi.push_back(make_pair(b.ak4pfjets_deepCSV()->at(i),b.ak4pfjets_phi()->at(i)));
        }
      }

      sort(v_csvPhi.begin(),v_csvPhi.end());

      int n = v_csvPhi.size();

      if(n>=2){
        delphi = deltaPhi(v_csvPhi[n-1].second,v_csvPhi[n-2].second);
        return delphi;
      }else return 0;

    });

  const NamedFunc bmetMinDeltaPhi("bmetMinDeltaPhi",[](const Baby &b) -> NamedFunc::ScalarType{
      float delphi=10;
      
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5&&delphi>b.dphi_ak4pfjet_met()->at(i)){
          delphi = b.dphi_ak4pfjet_met()->at(i);
        }
      }

      return delphi;

    });

  const NamedFunc nHeavy("nHeavy",[](const Baby &b) -> NamedFunc::ScalarType{
      float njets=0;
      
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==5||abs(b.ak4pfjets_parton_flavor()->at(i))==4){
          njets++;
        }
      }

      return njets;

    });

  const NamedFunc nLight("nLight",[](const Baby &b) -> NamedFunc::ScalarType{
      float njets=0;
      
      for(unsigned i(0); i<b.ak4pfjets_parton_flavor()->size(); i++){
        if(abs(b.ak4pfjets_parton_flavor()->at(i))==1||abs(b.ak4pfjets_parton_flavor()->at(i))==2||abs(b.ak4pfjets_parton_flavor()->at(i))==3){
          njets++;
        }
      }

      return njets;

    });

  const NamedFunc gluBTagged("gluBTagged",[](const Baby &b) -> NamedFunc::ScalarType{
    vector <pair <float, float> > gen_bs;

    for(unsigned i(0); i<b.gen_id()->size(); i++){
      if(abs(b.gen_id()->at(i))==5&&b.gen_motherid()->at(i)==21){
        for(unsigned j(i+1);j<b.gen_id()->size();j++){
          if(b.gen_id()->at(j)==-b.gen_id()->at(i)&&b.gen_motheridx()->at(j)==b.gen_motheridx()->at(i)&&(b.gen_pt()->at(j)>30||b.gen_pt()->at(i)>30)){
            gen_bs.push_back(make_pair(b.gen_eta()->at(i),b.gen_phi()->at(i)));
            gen_bs.push_back(make_pair(b.gen_eta()->at(j),b.gen_phi()->at(j)));
          }//Close if statement for opposite b
        }//Close for loop over second half of particles in event
      }//Close if statement that find b
    }//Close for loop over all particles in event

    float ntagged = 0;
    float delR = -1;
    float mindelR = 10;
    float dcsv = -1;

    int idx1 = 0;
    int idx2 = 0;

    if (gen_bs.size()==0){
      ntagged = -1;
    }else{
      for(unsigned i(0);i<2;i++){
        for(unsigned j(0);j<b.ak4pfjets_eta()->size();j++){
          delR = deltaR(gen_bs[i].first,gen_bs[i].second,b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j));
          if(delR < mindelR){
            mindelR = delR;
            dcsv = b.ak4pfjets_deepCSV()->at(j);
            if(i==0){
              idx1=j;
            }else{
              idx2=j;
            }
          }//Close if statement for smallest dR
        }//Close for loop over reco objects
        if(dcsv > 0.6324) ntagged++;
      }//Close for loop over gen objects

      if(idx1==idx2){
        ntagged = ntagged-1;
      }
    }

    return ntagged;

    });

  const NamedFunc outsideHiggsWindow("outsideHiggsWindow",[](const Baby &b) -> NamedFunc::ScalarType{

    int outsideWindow = 0.;

    if(b.mbb()<90||b.mbb()>150||b.weight()!=1.){
      outsideWindow=1.;
    }

    return outsideWindow;
  });

  const NamedFunc passTriggers("passTriggers",[](const Baby &b) -> NamedFunc::ScalarType{

    int pass = 0.;

    if(  b.HLT_SingleEl()==1||b.HLT_SingleMu()==1||b.HLT_MET_MHT()==1){
      pass=1.;
    }

    return pass;
  });

  const NamedFunc nJetsGood("nJetsGood",[](const Baby &b) -> NamedFunc::ScalarType{

    int pass = 0.;
    pass = b.ngoodjets();

    return pass;
  });

  const NamedFunc mht("mht",[](const Baby &b) -> NamedFunc::ScalarType{
    double mht_var = 0;
    double x = 0;
    double y = 0;

    for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      TLorentzVector v1;
      v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
      x -= v1.Px();
      y -= v1.Py();
    }//Close for loop over all jets in event

    for(unsigned i(0); i<b.leps_pt()->size(); i++){
      TLorentzVector v2;
      v2.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
      x -= v2.Px();
      y -= v2.Py();
    }//Close for loop over all leps in event

    mht_var = sqrt((x*x)+(y*y));

    return mht_var;

    });

  const NamedFunc mht_phi("mht_phi",[](const Baby &b) -> NamedFunc::ScalarType{
    double mht_var = 0;
    double x = 0;
    double y = 0;

    for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      TLorentzVector v1;
      v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
      x -= v1.Px();
      y -= v1.Py();
    }//Close for loop over all jets in event

    for(unsigned i(0); i<b.leps_pt()->size(); i++){
      TLorentzVector v2;
      v2.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
      x -= v2.Px();
      y -= v2.Py();
    }//Close for loop over all leps in event

    mht_var = atan2(x, y);
    return mht_var;

    });

  const NamedFunc W_pt_lep_met("W_pt_lep_met",[](const Baby &b) -> NamedFunc::ScalarType{
    double W_pt_var = 0;

    TLorentzVector lep;
    TLorentzVector met;
    TLorentzVector W_cand;

    lep.SetPtEtaPhiM(b.leps_pt()->at(0), b.leps_eta()->at(0), b.leps_phi()->at(0), 0.);
    met.SetPtEtaPhiM(b.pfmet(), 0., b.pfmet_phi(), 0.); //b.genmet_phi()

    W_cand = lep + met;

    W_pt_var = W_cand.Pt();
    return W_pt_var;

    });

  const NamedFunc W_pt_lep_mht("W_pt_lep_mht",[](const Baby &b) -> NamedFunc::ScalarType{
    double W_pt_var = 0;

    TLorentzVector lep;
    TLorentzVector met;
    TLorentzVector W_cand;

    double var_mht_pt = 0;
    double var_mht_phi = 0;
    double x = 0;
    double y = 0;

    for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      TLorentzVector v1;
      v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
      x -= v1.Px();
      y -= v1.Py();
    }//Close for loop over all jets in event

    for(unsigned i(0); i<b.leps_pt()->size(); i++){
      TLorentzVector v2;
      v2.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
      x -= v2.Px();
      y -= v2.Py();
    }//Close for loop over all leps in event

    var_mht_phi = atan2(x, y);
    var_mht_pt = sqrt((x*x)+(y*y));

    lep.SetPtEtaPhiM(b.leps_pt()->at(0), b.leps_eta()->at(0), b.leps_phi()->at(0), 0.);
    met.SetPtEtaPhiM(var_mht_pt, 0., var_mht_phi, 0.); //b.genmet_phi()

    W_cand = lep + met;

    W_pt_var = W_cand.Pt();
    return W_pt_var;

    });

  const NamedFunc mt_lep_mht("mt_lep_mht",[](const Baby &b) -> NamedFunc::ScalarType{
    double mt_var = 0;

    double var_mht_pt = 0;
    double var_mht_phi = 0;
    double x = 0;
    double y = 0;

    for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
      TLorentzVector v1;
      v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
      x -= v1.Px();
      y -= v1.Py();
    }//Close for loop over all jets in event

    for(unsigned i(0); i<b.leps_pt()->size(); i++){
      TLorentzVector v2;
      v2.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
      x -= v2.Px();
      y -= v2.Py();
    }//Close for loop over all leps in event

    var_mht_phi = atan2(x, y);
    var_mht_pt = sqrt((x*x)+(y*y));

    mt_var = sqrt(2*b.leps_pt()->at(0)*var_mht_pt * (1-(cos(b.leps_phi()->at(0)-var_mht_phi))));
    return mt_var;

    });

  const NamedFunc mt_lep_met_rec("mt_lep_met_rec",[](const Baby &b) -> NamedFunc::ScalarType{
    double mt_var = 0;

    mt_var = sqrt(2*b.leps_pt()->at(0)*b.pfmet() * (1-(cos(b.leps_phi()->at(0)-b.pfmet_phi()))));
    return mt_var;

    });

  const NamedFunc dijet_mass("dijet_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    float mass = -1;
    if (b.ak4pfjets_pt()->size()>1){
        TLorentzVector v1,v2,sum;
        v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(0),b.ak4pfjets_eta()->at(0),b.ak4pfjets_phi()->at(0),0);
        v2.SetPtEtaPhiM(b.ak4pfjets_pt()->at(1),b.ak4pfjets_eta()->at(1),b.ak4pfjets_phi()->at(1),0);
        sum=v1+v2;
        mass = sum.M();
        }
    return mass;
    });

  const NamedFunc dilepton_mass("dilepton_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    float mass = -1;
    float newmass = -1;
    for(unsigned i(0); i<b.leps_pt()->size(); i++){
        for(unsigned j(i+1);j<b.leps_pt()->size();j++){
          if(b.leps_pdgid()->at(j)==-b.leps_pdgid()->at(i)){
            TLorentzVector v1,v2,sum;
            v1.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
            v2.SetPtEtaPhiM(b.leps_pt()->at(j),b.leps_eta()->at(j),b.leps_phi()->at(j),0);
            sum=v1+v2;
            mass = sum.M();
            if (abs(mass-91.2) < abs(newmass-91.2)) newmass=mass;
          }//Close if statement for opposite pdgId
        }//Close for loop over second half of particles in event
    }//Close for loop over all particles in event
    return newmass;

    });

  const NamedFunc dilepton_pt("dilepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float mass = -1;
    float newmass = -1;
    float dl_pt = -1;
    for(unsigned i(0); i<b.leps_pt()->size(); i++){
        for(unsigned j(i+1);j<b.leps_pt()->size();j++){
          if(b.leps_pdgid()->at(j)==-b.leps_pdgid()->at(i)){
            TLorentzVector v1,v2,sum;
            v1.SetPtEtaPhiM(b.leps_pt()->at(i),b.leps_eta()->at(i),b.leps_phi()->at(i),0);
            v2.SetPtEtaPhiM(b.leps_pt()->at(j),b.leps_eta()->at(j),b.leps_phi()->at(j),0);
            sum=v1+v2;
            mass = sum.M();
            if (abs(mass-91.2) < abs(newmass-91.2)){
                newmass=mass;
                dl_pt = sum.Pt();
            }
          }//Close if statement for opposite pdgId
        }//Close for loop over second half of particles in event
    }//Close for loop over all particles in event
    return dl_pt;

    });

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Basic Jet Pt

  /*
   * returns pt of highest pt jet
   */
  const NamedFunc LeadingJetPt("LeadingJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt) maxpt =  b.ak4pfjets_pt()->at(i);
      }
      return maxpt;
    });
  

  /*
   * returns pt of second highest pt jet
   */
  const NamedFunc SubLeadingJetPt("SubLeadingJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      float secondmaxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > secondmaxpt && b.ak4pfjets_pt()->at(i)<=maxpt){
          secondmaxpt = b.ak4pfjets_pt()->at(i);
        }
        else if (b.ak4pfjets_pt()->at(i) > maxpt){
          secondmaxpt=maxpt;
          maxpt =  b.ak4pfjets_pt()->at(i);
        }
      }
      return secondmaxpt;
    });

  /*
   * returns pt of the third highest pt jet
   */
  const NamedFunc SubSubLeadingJetPt("SubSubLeadingJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      float secondmaxpt=0;
      float thirdmaxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > thirdmaxpt && b.ak4pfjets_pt()->at(i)<=secondmaxpt){
          thirdmaxpt = b.ak4pfjets_pt()->at(i);
        }
        else if (b.ak4pfjets_pt()->at(i) > maxpt){
          thirdmaxpt=secondmaxpt;
          secondmaxpt=maxpt;
          maxpt=b.ak4pfjets_pt()->at(i);
        }
        else if (b.ak4pfjets_pt()->at(i) <= maxpt && b.ak4pfjets_pt()->at(i) > secondmaxpt){
          thirdmaxpt = secondmaxpt;
          secondmaxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return thirdmaxpt;
    });


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Funcs to look at number of 3rd jets which are actually B's

  /*
   * Pt of the non-bjet with the highest pt
   */
  const NamedFunc LeadingNonBJetPt("LeadingNonBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });

  const NamedFunc LeadingFakeNonBJetPt("LeadingFakeNonBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && abs(b.ak4pfjets_hadron_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });

   const NamedFunc LeadingRealNonBJetPt("LeadingRealNonBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && abs(b.ak4pfjets_hadron_flavor()->at(i))!=5 && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });

   const NamedFunc LeadingNonBJetPt_med("LeadingNonBJetPt_med",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && b.ak4pfjets_deepCSV()->at(i) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });
   const NamedFunc LeadingNonBJetPt_med_jup("LeadingNonBJetPt_med_jup",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxpt=0;
    for (unsigned i(0); i<b.ak4pfjets_pt_jup()->size(); i++){
      if (b.ak4pfjets_pt_jup()->at(i) > maxpt && b.ak4pfjets_deepCSV()->at(i) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
         maxpt = b.ak4pfjets_pt_jup()->at(i);
      }
    }
    return maxpt;
  });
    const NamedFunc LeadingNonBJetPt_med_jdown("LeadingNonBJetPt_med_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxpt=0;
    for (unsigned i(0); i<b.ak4pfjets_pt_jdown()->size(); i++){
      if (b.ak4pfjets_pt_jdown()->at(i) > maxpt && b.ak4pfjets_deepCSV()->at(i) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
         maxpt = b.ak4pfjets_pt_jdown()->at(i);
      }
    }
    return maxpt;
  });

   const NamedFunc max_ak8pfjets_deepdisc_hbb("max_ak8pfjets_deepdisc_hbb",[](const Baby &b) -> NamedFunc::ScalarType{
      float max_disc= -0.05;
      for (unsigned i(0); i<b.ak8pfjets_deepdisc_hbb()->size(); i++){
        if (b.ak8pfjets_deepdisc_hbb()->at(i) > max_disc){
           max_disc = b.ak8pfjets_deepdisc_hbb()->at(i);
        }
      }
      return max_disc;
    });

   const NamedFunc PassThirdJetHighpTVeto("PassThirdJetHighpTVeto",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && b.ak4pfjets_deepCSV()->at(i) < 0.6321){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      if(b.ak4pfjets_deepCSV()->size() <3 || maxpt<100) return true;
      else return false;
    });

  const NamedFunc LeadingFakeNonBJetPt_med("LeadingFakeNonBJetPt_med",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && abs(b.ak4pfjets_hadron_flavor()->at(i))==5 && b.ak4pfjets_deepCSV()->at(i) < 0.6321){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });

   const NamedFunc LeadingRealNonBJetPt_med("LeadingRealNonBJetPt_med",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > maxpt && abs(b.ak4pfjets_hadron_flavor()->at(i))!=5 && b.ak4pfjets_deepCSV()->at(i) < 0.6321){
           maxpt = b.ak4pfjets_pt()->at(i);
        }
      }
      return maxpt;
    });

   const NamedFunc nNonBTagJets("nNonBTagJets",[](const Baby &b) -> NamedFunc::ScalarType{
      float njet=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           njet++;
        }
      }
      return njet;
    });
  
   const NamedFunc signature_pT("signature_pT",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector jet1;
      TLorentzVector jet2;
      TLorentzVector lep1;
      TLorentzVector lep2;
      TLorentzVector sum;

      unsigned idx(0);

      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           jet1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
           idx = i;
           break;
        }
      }

      for (unsigned j(0); j<b.ak4pfjets_deepCSV()->size(); j++){
        if (j>idx && b.ak4pfjets_deepCSV()->at(j) < (0.6321*(b.year()==2016) + 0.4941*(b.year()==2017) + 0.4184*(b.year()==2018))){
           jet2.SetPtEtaPhiM(b.ak4pfjets_pt()->at(j),b.ak4pfjets_eta()->at(j),b.ak4pfjets_phi()->at(j),b.ak4pfjets_m()->at(j));
           break;
        }
      }

      if (b.leps_pt()->size()>=2){
        lep1.SetPtEtaPhiM(b.leps_pt()->at(0),b.leps_eta()->at(0),b.leps_phi()->at(0),0.);
        lep2.SetPtEtaPhiM(b.leps_pt()->at(1),b.leps_eta()->at(1),b.leps_phi()->at(1),0.);
      }

      sum = (jet1+jet2+lep1+lep2);

      return sum.Pt();
    });


   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // Funcs to look at number of B-tags which are fakes

   const NamedFunc LeadingBJetPt("LeadingBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > maxCSV && b.ak4pfjets_deepCSV()->at(i) > 0.2217){
           maxCSV = b.ak4pfjets_deepCSV()->at(i);
           maxpt = b.ak4pfjets_pt()->at(i);

        }
      }
      return maxpt;
    });

  const NamedFunc LeadingRealBJetPt("LeadingRealBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      float output=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > maxCSV && b.ak4pfjets_deepCSV()->at(i) > 0.2217){
           maxCSV = b.ak4pfjets_deepCSV()->at(i);
           maxpt = b.ak4pfjets_pt()->at(i);
           if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) output=maxpt;
           else if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) output=0;
        }
      }
      return output;
    });

   const NamedFunc LeadingFakeBJetPt("LeadingFakeBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      float output=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > maxCSV && b.ak4pfjets_deepCSV()->at(i) > 0.2217){
           maxCSV = b.ak4pfjets_deepCSV()->at(i);
           maxpt = b.ak4pfjets_pt()->at(i);
           if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) output=maxpt;
           else if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) output=0;
        }
      }
      return output;
    });

    /*
   * returns sub-leading b jet pt
   */
  const NamedFunc SubLeadingBJetPt("SubLeadingBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      float secondmaxCSV=0;
      float secondmaxpt=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > secondmaxCSV && b.ak4pfjets_deepCSV()->at(i)<=maxCSV){
            secondmaxCSV = b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt = b.ak4pfjets_pt()->at(i);
        }
        else if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > maxCSV){
            secondmaxCSV=maxCSV;
            maxCSV =  b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt=maxpt;
            maxpt =  b.ak4pfjets_pt()->at(i);
        }
      }
      return secondmaxpt;
    });

  const NamedFunc SubLeadingRealBJetPt("SubLeadingRealBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      bool maxisb=false;
      float secondmaxCSV=0;
      float secondmaxpt=0;
      float output=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > secondmaxCSV && b.ak4pfjets_deepCSV()->at(i)<=maxCSV){
            secondmaxCSV = b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt = b.ak4pfjets_pt()->at(i);
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) output=secondmaxpt;
            else if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) output=0;
        }
        else if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > maxCSV){
            secondmaxCSV=maxCSV;
            maxCSV =  b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt=maxpt;
            maxpt =  b.ak4pfjets_pt()->at(i);
            if (maxisb) output=secondmaxpt;
            if (!maxisb) output=0;
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) maxisb=true;
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) maxisb=false;
        }
      }
      return output;
    });

  const NamedFunc SubLeadingFakeBJetPt("SubLeadingFakeBJetPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float maxCSV=0;
      float maxpt=0;
      bool maxisb=false;
      float secondmaxCSV=0;
      float secondmaxpt=0;
      float output=0;
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > secondmaxCSV && b.ak4pfjets_deepCSV()->at(i)<=maxCSV){
            secondmaxCSV = b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt = b.ak4pfjets_pt()->at(i);
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) output=secondmaxpt;
            else if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) output=0;
         }
        else if (b.ak4pfjets_deepCSV()->at(i) > 0.2217 && b.ak4pfjets_deepCSV()->at(i) > maxCSV){
            secondmaxCSV=maxCSV;
            maxCSV =  b.ak4pfjets_deepCSV()->at(i);
            secondmaxpt=maxpt;
            maxpt =  b.ak4pfjets_pt()->at(i);
            if (!maxisb) output=secondmaxpt;
            if (maxisb) output=0;
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))==5) maxisb=true;
            if (abs(b.ak4pfjets_hadron_flavor()->at(i))!=5) maxisb=false;
        }
      }
      return output;
    });


   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // Looking at gen level ISR

    /*
    * Pt ISR using truth info
    */
    const NamedFunc genISRPt("genISRPt",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector *obj1 = new TLorentzVector();
      TLorentzVector *obj2 = new TLorentzVector();
      TLorentzVector ISR;
      int obj1_idx=1000;
      int obj2_idx=1000;
      float ISRpt=-1;
      for (unsigned i(0); i<b.gen_id()->size(); i++){
        if (abs(b.gen_id()->at(i))==1000024) {obj1_idx=i;}
        if (abs(b.gen_id()->at(i))==1000023) {obj2_idx=i;}
        if (b.gen_id()->at(i)==6) {obj1_idx=i;}
        if (b.gen_id()->at(i)==-6) {obj2_idx=i;}
      }
      if (obj1_idx!=1000 && obj2_idx!=1000) {
        obj1->SetPtEtaPhiM(b.gen_pt()->at(obj1_idx), b.gen_eta()->at(obj1_idx), b.gen_phi()->at(obj1_idx), b.gen_m()->at(obj1_idx));
        obj2->SetPtEtaPhiM(b.gen_pt()->at(obj2_idx), b.gen_eta()->at(obj2_idx), b.gen_phi()->at(obj2_idx), b.gen_m()->at(obj2_idx));
        ISR = -(*obj1+*obj2);
        ISRpt = abs(ISR.Pt());
      }
      delete obj1;
      delete obj2;
      return ISRpt;
    });

    const NamedFunc genISRgenMETdPhi("genISRMETdPhi",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector *obj1 = new TLorentzVector();
      TLorentzVector *obj2 = new TLorentzVector();
      TLorentzVector ISR;
      TLorentzVector MET;
      int obj1_idx=1000;
      int obj2_idx=1000;
      int nmet_sources=0;
      float dR_ISR_MET=-1;
      for (unsigned i(0); i<b.gen_id()->size(); i++){
        if (abs(b.gen_id()->at(i))==1000024) {obj1_idx=i;}
        if (abs(b.gen_id()->at(i))==1000023) {obj2_idx=i;}
        if (b.gen_id()->at(i)==6) {obj1_idx=i;}
        if (b.gen_id()->at(i)==-6) {obj2_idx=i;}
        if (abs(b.gen_id()->at(i))==12||abs(b.gen_id()->at(i))==14||abs(b.gen_id()->at(i))==16
            ||abs(b.gen_id()->at(i))==1000022) {
          TLorentzVector *nu = new TLorentzVector();
          nu->SetPtEtaPhiM(b.gen_pt()->at(i), b.gen_eta()->at(i), b.gen_phi()->at(i), b.gen_m()->at(i));
          MET = MET + *nu;
          nmet_sources++;
          delete nu;
        }
      }
      if (obj1_idx!=1000 && obj2_idx!=1000) {
        obj1->SetPtEtaPhiM(b.gen_pt()->at(obj1_idx), b.gen_eta()->at(obj1_idx), b.gen_phi()->at(obj1_idx), b.gen_m()->at(obj1_idx));
        obj2->SetPtEtaPhiM(b.gen_pt()->at(obj2_idx), b.gen_eta()->at(obj2_idx), b.gen_phi()->at(obj2_idx), b.gen_m()->at(obj2_idx));
        ISR = -(*obj1+*obj2);
      }
      if (nmet_sources!=0) dR_ISR_MET = ISR.DeltaPhi(MET);
      delete obj1;
      delete obj2;
      return dR_ISR_MET;
    });

    const NamedFunc genISRrecoMETdPhi("genISRrecoMETdPhi",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector *obj1 = new TLorentzVector();
      TLorentzVector *obj2 = new TLorentzVector();
      TLorentzVector ISR;
      TLorentzVector MET;
      int obj1_idx=1000;
      int obj2_idx=1000;
      int nmet_sources=0;
      float dR_ISR_MET=-1;
      for (unsigned i(0); i<b.gen_id()->size(); i++){
        if (abs(b.gen_id()->at(i))==1000024) {obj1_idx=i;}
        if (abs(b.gen_id()->at(i))==1000023) {obj2_idx=i;}
        if (b.gen_id()->at(i)==6) {obj1_idx=i;}
        if (b.gen_id()->at(i)==-6) {obj2_idx=i;}
        if (abs(b.gen_id()->at(i))==12||abs(b.gen_id()->at(i))==14||abs(b.gen_id()->at(i))==16
            ||abs(b.gen_id()->at(i))==1000022) {
          TLorentzVector *nu = new TLorentzVector();
          nu->SetPtEtaPhiM(b.gen_pt()->at(i), b.gen_eta()->at(i), b.gen_phi()->at(i), b.gen_m()->at(i));
          MET = MET + *nu;
          nmet_sources++;
          delete nu;
        }
      }
      if (obj1_idx!=1000 && obj2_idx!=1000) {
        obj1->SetPtEtaPhiM(b.gen_pt()->at(obj1_idx), b.gen_eta()->at(obj1_idx), b.gen_phi()->at(obj1_idx), b.gen_m()->at(obj1_idx));
        obj2->SetPtEtaPhiM(b.gen_pt()->at(obj2_idx), b.gen_eta()->at(obj2_idx), b.gen_phi()->at(obj2_idx), b.gen_m()->at(obj2_idx));
        ISR = -(*obj1+*obj2);
      }
      if (nmet_sources!=0) dR_ISR_MET = ISR.DeltaPhi(MET);
      delete obj1;
      delete obj2;
      return dR_ISR_MET;
    });


    const NamedFunc genISRrecoISRdPhi("genISRrecoISRdPhi",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector *obj1 = new TLorentzVector();
      TLorentzVector *obj2 = new TLorentzVector();
      TLorentzVector gISR;
      TLorentzVector rISR;
      int obj1_idx=1000;
      int obj2_idx=1000;
      int rISR_idx=1000;
      float rISRpt=0;
      float gISR_rISR_dPhi=-1;
      for (unsigned i(0); i<b.gen_id()->size(); i++){
        if (abs(b.gen_id()->at(i))==1000024) {obj1_idx=i;}
        if (abs(b.gen_id()->at(i))==1000023) {obj2_idx=i;}
        if (b.gen_id()->at(i)==6) {obj1_idx=i;}
        if (b.gen_id()->at(i)==-6) {obj2_idx=i;}
      }
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > rISRpt && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
           rISRpt = b.ak4pfjets_pt()->at(i);
           rISR_idx = i;
        }
      }
      if (obj1_idx!=1000 && obj2_idx!=1000 && rISR_idx!=1000) {
        obj1->SetPtEtaPhiM(b.gen_pt()->at(obj1_idx), b.gen_eta()->at(obj1_idx), 
                            b.gen_phi()->at(obj1_idx), b.gen_m()->at(obj1_idx));
        obj2->SetPtEtaPhiM(b.gen_pt()->at(obj2_idx), b.gen_eta()->at(obj2_idx), 
                            b.gen_phi()->at(obj2_idx), b.gen_m()->at(obj2_idx));
        rISR.SetPtEtaPhiM(b.ak4pfjets_pt()->at(rISR_idx), b.ak4pfjets_eta()->at(rISR_idx), 
                            b.ak4pfjets_phi()->at(rISR_idx), b.ak4pfjets_m()->at(rISR_idx));
        gISR = -(*obj1+*obj2);
        gISR_rISR_dPhi = rISR.DeltaPhi(gISR);
      }
      delete obj1;
      delete obj2;
      return gISR_rISR_dPhi;
    });


    const NamedFunc genISRrecoISRDeltaPt("genISRrecoISRDeltaPt",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector *obj1 = new TLorentzVector();
      TLorentzVector *obj2 = new TLorentzVector();
      TLorentzVector gISR;
      TLorentzVector rISR;
      int obj1_idx=1000;
      int obj2_idx=1000;
      int rISR_idx=1000;
      float rISRpt=0;
      float gISR_rISR_DeltaPt=1000;
      for (unsigned i(0); i<b.gen_id()->size(); i++){
        if (abs(b.gen_id()->at(i))==1000024) {obj1_idx=i;}
        if (abs(b.gen_id()->at(i))==1000023) {obj2_idx=i;}
        if (b.gen_id()->at(i)==6) {obj1_idx=i;}
        if (b.gen_id()->at(i)==-6) {obj2_idx=i;}
      }
      for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++){
        if (b.ak4pfjets_pt()->at(i) > rISRpt && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
           rISRpt = b.ak4pfjets_pt()->at(i);
           rISR_idx = i;
        }
      }
      if (obj1_idx!=1000 && obj2_idx!=1000 && rISR_idx!=1000) {
        obj1->SetPtEtaPhiM(b.gen_pt()->at(obj1_idx), b.gen_eta()->at(obj1_idx), 
                            b.gen_phi()->at(obj1_idx), b.gen_m()->at(obj1_idx));
        obj2->SetPtEtaPhiM(b.gen_pt()->at(obj2_idx), b.gen_eta()->at(obj2_idx), 
                            b.gen_phi()->at(obj2_idx), b.gen_m()->at(obj2_idx));
        rISR.SetPtEtaPhiM(b.ak4pfjets_pt()->at(rISR_idx), b.ak4pfjets_eta()->at(rISR_idx), 
                            b.ak4pfjets_phi()->at(rISR_idx), b.ak4pfjets_m()->at(rISR_idx));
        gISR = -(*obj1+*obj2);
        gISR_rISR_DeltaPt = rISR.Pt()-gISR.Pt();
      }
      delete obj1;
      delete obj2;
      return gISR_rISR_DeltaPt;
    });

    const NamedFunc truthOrigin3rdJet("truthOrigin3rdJet",[](const Baby &b) -> NamedFunc::ScalarType{
        TLorentzVector *gOBJ = new TLorentzVector();
        TLorentzVector *rISR = new TLorentzVector();
        int jet_origin=0;
        float dR_rISR_gOBJ = 0.5;
        int gOBJ_idx=1000;
        int rISR_idx=1000;
        bool gOBJ_parent_of_leptons = false;
        float rISRpt=0;
        for (unsigned i(0); i<b.ak4pfjets_deepCSV()->size(); i++) {
          if (b.ak4pfjets_pt()->at(i) > rISRpt && b.ak4pfjets_deepCSV()->at(i) < 0.2217){
            rISRpt = b.ak4pfjets_pt()->at(i);
            rISR_idx = i;
          }
        }
        if (rISR_idx!=1000) { 
          rISR->SetPtEtaPhiM(b.ak4pfjets_pt()->at(rISR_idx), b.ak4pfjets_eta()->at(rISR_idx), 
                            b.ak4pfjets_phi()->at(rISR_idx), b.ak4pfjets_m()->at(rISR_idx));
        }
        for (unsigned i(0); i<b.gen_id()->size(); i++) {
          gOBJ->SetPtEtaPhiM(b.gen_pt()->at(i), b.gen_eta()->at(i), 
                            b.gen_phi()->at(i), b.gen_m()->at(i));
          if (gOBJ->DeltaR(*rISR)<dR_rISR_gOBJ && abs(b.gen_id()->at(i))!=24 && abs(b.gen_id()->at(i))!=6 
              && abs(b.gen_id()->at(i))!=12 && abs(b.gen_id()->at(i))!=14 && abs(b.gen_id()->at(i))!=16
              && abs(b.gen_id()->at(i))!=1000024 && abs(b.gen_id()->at(i))!=1000023 && abs(b.gen_id()->at(i))!=1000022
              && abs(b.gen_id()->at(i))!=22) {
                gOBJ_idx = i;
                dR_rISR_gOBJ = gOBJ->DeltaR(*rISR);
          }
        }
        for (unsigned i(0); i<b.gen_id()->size(); i++) {
          if ((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||
              abs(b.gen_id()->at(i))==15)&& b.gen_motheridx()->at(i)==gOBJ_idx) gOBJ_parent_of_leptons = true;
        }
        if (gOBJ_idx!=1000) {
          if (abs(b.gen_id()->at(gOBJ_idx))==5 && abs(b.gen_motherid()->at(gOBJ_idx))==6) jet_origin=1;
          else if (abs(b.gen_id()->at(gOBJ_idx))<=6 && abs(b.gen_motherid()->at(gOBJ_idx))==24 && abs(b.gen_gmotherid()->at(gOBJ_idx))==6) jet_origin=2;
          else if (abs(b.gen_id()->at(gOBJ_idx))==15 && abs(b.gen_motherid()->at(gOBJ_idx))==24 && abs(b.gen_gmotherid()->at(gOBJ_idx))==6 && gOBJ_parent_of_leptons) jet_origin=3;
          else if (abs(b.gen_id()->at(gOBJ_idx))==15 && abs(b.gen_motherid()->at(gOBJ_idx))==24 && abs(b.gen_gmotherid()->at(gOBJ_idx))==6) jet_origin=4;
          else jet_origin=5;
        }
        delete gOBJ;
        delete rISR;

        return jet_origin;
    });


//Functions for Data Quality Checks
      const NamedFunc pfmet_uncorr_func("pfmet_uncorr_func",[](const Baby &b) -> NamedFunc::ScalarType{
      float pfmet_uncorr_var = 0;
      if(b.pfmet_uncorr()!= -9999){
        pfmet_uncorr_var  = b.pfmet_uncorr();
      }else{
        pfmet_uncorr_var = b.pfmet();
      }
      return pfmet_uncorr_var;
    });

  const NamedFunc mt_met_lep_uncorr("mt_met_lep_uncorr",[](const Baby &b) -> NamedFunc::ScalarType{
      float mt_met_lep_uncorr_var = 0;
      if(b.pfmet_uncorr()!= -9999){
        float phi1 = b.leps_phi()->at(0);
        float phi2 = b.pfmet_uncorr_phi();
        float Et1  = b.leps_pt()->at(0);
        float Et2  = b.pfmet_uncorr();

        mt_met_lep_uncorr_var = sqrt(2*Et1*Et2*(1.0-cos(phi1-phi2)));
      }else{
        mt_met_lep_uncorr_var = b.mt_met_lep();
      }
      return mt_met_lep_uncorr_var;
    });

  const NamedFunc noPrefireWeight("noPrefireWeight",[](const Baby &b) -> NamedFunc::ScalarType{
      float noPrefireCuts = 0;

      noPrefireCuts = b.pass()/b.w_L1();

      return noPrefireCuts;
    });

  const NamedFunc HasHEMjet("HasHEMjet",[](const Baby &b) -> NamedFunc::ScalarType{
      int njet=0;

      for (unsigned i(0); i<b.ak4pfjets_phi()->size(); i++){
        if (b.run()>=319077 && b.ak4pfjets_phi()->at(i)>-1.6 && b.ak4pfjets_phi()->at(i) < -0.8 && b.ak4pfjets_eta()->at(i) > -2.4 && b.ak4pfjets_eta()->at(i) < -1.4){
          njet++;
        }
      }
      return njet;
    });

  const NamedFunc HasHEMevent("HasHEMevent",[](const Baby &b) -> NamedFunc::ScalarType{
      int nevent=0;

      for (unsigned i(0); i<b.vetoleps_pdgid()->size(); i++){
        if (b.run()>=319077 && abs(b.vetoleps_pdgid()->at(i))==11 && b.vetoleps_phi()->at(i)>-1.6 && b.vetoleps_phi()->at(i) < -0.8 && b.vetoleps_eta()->at(i) > -2.4 && b.vetoleps_eta()->at(i) < -1.4){
          nevent++;
        }
      }
      return nevent;
    });

  const NamedFunc mcHEMWeight("mcHEMWeight",[](const Baby &b) -> NamedFunc::ScalarType{
      float weight=0.;
      int nEls = 0;
      int nJets = 0;

      if(b.year()==2018&&b.genmet()!=-9999){

        for (unsigned i(0); i<b.vetoleps_pdgid()->size(); i++){
          if (abs(b.vetoleps_pdgid()->at(i))==11 && b.vetoleps_phi()->at(i)>-1.6 && b.vetoleps_phi()->at(i) < -0.8 && b.vetoleps_eta()->at(i) > -2.4 && b.vetoleps_eta()->at(i) < -1.4){
            nEls++;
          }
        }

        for (unsigned i(0); i<b.ak4pfjets_phi()->size(); i++){
          if (b.ak4pfjets_phi()->at(i)>-1.6 && b.ak4pfjets_phi()->at(i) < -0.8 && b.ak4pfjets_eta()->at(i) > -2.4 && b.ak4pfjets_eta()->at(i) < -1.4){
            nJets++;
          }
        }

      }

      if(nEls>0||nJets>0){
        weight=0.352;
        }else{
          weight = 1.;
        }
      return weight;
    });

  const NamedFunc muPt("muPt",[](const Baby &b) -> NamedFunc::ScalarType{
    float pt = 0;

    if(abs(b.lep1_pdgid())==13){
      pt = b.leps_pt()->at(0);
    }//Close if
    return pt;

    });

  const NamedFunc elPt("elPt",[](const Baby &b) -> NamedFunc::ScalarType{
    float pt = 0;

    if(abs(b.lep1_pdgid())==11){
      pt = b.leps_pt()->at(0);
    }//Close if
    return pt;

    });

  const NamedFunc pfmet_div_genmet("pfmet_div_genmet",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()/b.genmet();

      return met_var;

  });

  const NamedFunc pfmet_resup_div_genmet("pfmet_resup_div_genmet",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet_resup()/b.genmet();

      return met_var;

  });

  const NamedFunc pfmet_resdown_div_genmet("pfmet_resdown_div_genmet",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet_resdown()/b.genmet();

      return met_var;

  });

  const NamedFunc pfmet_subt_resup("pfmet_subt_resup",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()-b.pfmet_resup();

      return met_var;

  });

  const NamedFunc pfmet_subt_resdown("pfmet_subt_resdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()-b.pfmet_resdown();

      return met_var;

  });

  const NamedFunc pfmet_subt_genmet("pfmet_subt_genmet",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()-b.genmet();

      return met_var;

  });

  const NamedFunc pfmet_subt_genmet_subt_lepPt("pfmet_subt_genmet_subt_lepPt",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;
      float pt_var1=0.;
      float pt_var2=0.;

    for (unsigned i(0); i<b.gen_id()->size(); i++){
        if ((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&b.gen_pt()->at(i)>pt_var1){
          pt_var2 = pt_var1;
          pt_var1 = b.gen_pt()->at(i);
        }else if((abs(b.gen_id()->at(i))==11||abs(b.gen_id()->at(i))==13||abs(b.gen_id()->at(i))==15)&&b.gen_pt()->at(i)>pt_var2){
          pt_var2 = b.gen_pt()->at(i);
        }
      }

      met_var = b.pfmet()-b.genmet()-pt_var2;

      return met_var;

  });

  const NamedFunc pfmet_subt_jup("pfmet_subt_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()-b.pfmet_jup();

      return met_var;

  });

  const NamedFunc pfmet_subt_jdown("pfmet_subt_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0.;

      met_var = b.pfmet()-b.pfmet_jdown();

      return met_var;

  });

  const NamedFunc mt_subt_jup("mt_subt_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      float mt_var=0.;

      mt_var = b.mt_met_lep()-b.mt_met_lep_jup();

      return mt_var;

  });

  const NamedFunc mt_subt_jdown("mt_subt_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float mt_var=0.;

      mt_var = b.mt_met_lep()-b.mt_met_lep_jdown();

      return mt_var;

  });

  const NamedFunc mct_subt_jup("mct_subt_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      float mct_var=0.;

      mct_var = b.mct()-b.jup_mct();

      return mct_var;

  });

  const NamedFunc mct_subt_jdown("mct_subt_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float mct_var=0.;

      mct_var = b.mct()-b.jdown_mct();

      return mct_var;

  });

  const NamedFunc mbb_subt_jup("mbb_subt_jup",[](const Baby &b) -> NamedFunc::ScalarType{
      float mbb_var=0.;

      mbb_var = b.mbb()-b.jup_mbb();

      return mbb_var;

  });

  const NamedFunc mbb_subt_jdown("mbb_subt_jdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float mbb_var=0.;

      mbb_var = b.mbb()-b.jdown_mbb();

      return mbb_var;

  });

  const NamedFunc mct_subt_jerup("mct_subt_jerup",[](const Baby &b) -> NamedFunc::ScalarType{
      float mct_var=0.;

      mct_var = b.mct()-b.jerup_mct();

      return mct_var;

  });

  const NamedFunc mct_subt_jerdown("mct_subt_jerdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float mct_var=0.;

      mct_var = b.mct()-b.jerdown_mct();

      return mct_var;

  });

  const NamedFunc mbb_subt_jerup("mbb_subt_jerup",[](const Baby &b) -> NamedFunc::ScalarType{
      float mbb_var=0.;

      mbb_var = b.mbb()-b.jerup_mbb();

      return mbb_var;

  });

  const NamedFunc mbb_subt_jerdown("mbb_subt_jerdown",[](const Baby &b) -> NamedFunc::ScalarType{
      float mbb_var=0.;

      mbb_var = b.mbb()-b.jerdown_mbb();

      return mbb_var;

  });

   const NamedFunc genPt_ak8_higgsTagger("genPt_ak8_higgsTagger",[](const Baby &b) -> NamedFunc::ScalarType{
      float nHiggs=0;
      int year = b.year();
      float higgsWP = 0;

      if(year==2016){
        higgsWP = 0.8945;
      }else if(year==2017){
        higgsWP = 0.8695;
      }else if(year==2018){
        higgsWP = 0.8365;
      }

      for (unsigned i(0); i<b.FatJet_deepTagMD_HbbvsQCD()->size(); i++){
        if (b.FatJet_deepTagMD_HbbvsQCD()->at(i)>=higgsWP&&b.FatJet_genPt()->at(i)>250){
           nHiggs++;
        }else if(b.FatJet_genPt()->at(i)<0.&&b.FatJet_deepTagMD_HbbvsQCD()->at(i)>=higgsWP&&b.FatJet_pt_nom()->at(i)>250){
          nHiggs++;
        }
      }

      return nHiggs;
    });

   const NamedFunc reduced_MET("reduced_MET",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0;
      int year = b.year();

      if(year==2016){
        met_var = b.pfmet();
      }else{
        met_var = 0.97*b.pfmet();
      }

      return met_var;
    });

   const NamedFunc reducedMET_mT("reducedMET_mT",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0;
      float mt_var=0;
      int year = b.year();

      if(year==2016){
        met_var = b.pfmet();
      }else{
        met_var = 0.97*b.pfmet();
      }

      mt_var = sqrt(2*b.leps_pt()->at(0)*met_var * (1-(cos(b.leps_phi()->at(0)-b.pfmet_phi()))));
      return mt_var;

    });

   const NamedFunc fastsim_MET("fastsim_MET",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0;
      double x = 0;
      double y = 0;
      TLorentzVector met;
      TLorentzVector jets;
      TLorentzVector finalvar;
      int year = b.year();

      if(year==2016){
        met_var = b.pfmet();
      }else{
        for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
          TLorentzVector v1;
          if(b.ak4pfjets_hadron_flavor()->at(i)==5){
            v1.SetPtEtaPhiM((0.03*b.ak4pfjets_pt()->at(i)),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
          }/*else{
            v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
          }*/
          x += v1.Px();
          y += v1.Py();
        }//for loop over jet pt vector

        jets.SetPtEtaPhiM(sqrt((x*x)+(y*y)), 0., atan2(x,y), 0.);
        met.SetPtEtaPhiM(b.pfmet(),0.,b.pfmet_phi(),0.);
        finalvar = met+jets;
        met_var = finalvar.Pt();
      }

      return met_var;
    });

   const NamedFunc fastsim_MT("fastsim_MT",[](const Baby &b) -> NamedFunc::ScalarType{
      float met_var=0;
      float mt_var=0;
      double x = 0;
      double y = 0;
      TLorentzVector met;
      TLorentzVector jets;
      TLorentzVector finalvar;
      int year = b.year();

      if(year==2016){
        met_var = b.pfmet();
      }else{
        for(unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
          TLorentzVector v1;
          if(b.ak4pfjets_hadron_flavor()->at(i)==5){
            v1.SetPtEtaPhiM((0.03*b.ak4pfjets_pt()->at(i)),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
          }/*else{
            v1.SetPtEtaPhiM(b.ak4pfjets_pt()->at(i),b.ak4pfjets_eta()->at(i),b.ak4pfjets_phi()->at(i),b.ak4pfjets_m()->at(i));
          }*/
          x += v1.Px();
          y += v1.Py();
        }//for loop over jet pt vector

        jets.SetPtEtaPhiM(sqrt((x*x)+(y*y)), 0., atan2(x,y), 0.);
        met.SetPtEtaPhiM(b.pfmet(),0.,b.pfmet_phi(),0.);
        finalvar = met+jets;
        met_var = finalvar.Pt();
      }

      mt_var = sqrt(2*b.leps_pt()->at(0)*met_var * (1-(cos(b.leps_phi()->at(0)-b.pfmet_phi()))));

      return mt_var;
    });


    const NamedFunc fastsim_MCT("fastsim_MCT",[](const Baby &b) -> NamedFunc::ScalarType{
  
    float mct_var=0;
    int year = b.year();
    
    if(year==2016){
      mct_var = b.mct();
    }else{
      for (unsigned i(0); i<b.ak4pfjets_pt()->size(); i++){
        if (b.ak4pfjets_hadron_flavor()->at(i) == 5){
          for (unsigned j(i+1); j<b.ak4pfjets_pt()->size(); j++){
            if (b.ak4pfjets_hadron_flavor()->at(j) == 5){
              mct_var = sqrt(2*(0.97*b.ak4pfjets_pt()->at(j))*(0.97*b.ak4pfjets_pt()->at(i))*(1+cos(deltaPhi(b.ak4pfjets_phi()->at(i),b.ak4pfjets_phi()->at(j)))));
            }
          }
        }
      }
    }

    return mct_var;
    });

 
}
