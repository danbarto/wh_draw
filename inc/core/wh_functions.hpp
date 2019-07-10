#ifndef H_WH_FUNCTIONS
#define H_WH_FUNCTIONS

#include <cstddef>
#include <string>
#include "core/named_func.hpp"
 
namespace WH_Functions{

  // Miscilaneous
  extern const NamedFunc HasMedLooseCSV;
  extern const NamedFunc NHighPtNu;
  extern const NamedFunc HasMedLooseCSV;
  extern const NamedFunc HasMedMedDeepCSV;
  extern const NamedFunc HasExactMedMedDeepCSV;
  extern const NamedFunc HasLooseLooseDeepCSV;
  extern const NamedFunc HasMedLooseDeepCSV;
  extern const NamedFunc HasLooseNoMedDeepCSV;
  extern const NamedFunc nDeepMedBTagged;
  extern const NamedFunc nDeepLooseBTagged;
  extern const NamedFunc bJetPt;
  extern const NamedFunc HasNoBs;
  extern const NamedFunc WHLeptons;
  extern const NamedFunc NBJets;
  extern const NamedFunc NHighPtNu;
  extern const NamedFunc nRealBs;
  extern const NamedFunc sortedJetsPt_Leading;
  extern const NamedFunc sortedJetsPt_subLeading;
  extern const NamedFunc sortedJetsCSV_Leading;
  extern const NamedFunc sortedJetsCSV_subLeading;
  extern const NamedFunc sortedJetsCSV_deltaR;
  extern const NamedFunc nEventsGluonSplit;
  extern const NamedFunc nGenBs;
  extern const NamedFunc nGenBsFromGluons;
  extern const NamedFunc nGenBs_ptG15;
  extern const NamedFunc nGenBsFromGluons_ptG15;
  extern const NamedFunc nGenBs_ptG30;
  extern const NamedFunc nGenBsFromGluons_ptG30;
  extern const NamedFunc genBpT;
  extern const NamedFunc bDeltaRGluonSplit;
  extern const NamedFunc bDeltaR;
  extern const NamedFunc bMother;
  extern const NamedFunc bMother_pt15;
  extern const NamedFunc bMother_pt30;
  extern const NamedFunc genB_leadingpT;
  extern const NamedFunc genB_subleadingpT;
  extern const NamedFunc bDeltaPhi;
  extern const NamedFunc bmetMinDeltaPhi;
  extern const NamedFunc nHeavy;
  extern const NamedFunc nLight;
  extern const NamedFunc gluBTagged;
  extern const NamedFunc nModEventsGluonSplit;
  extern const NamedFunc leadingBMother_pt20;
  extern const NamedFunc subleadingBMother_pt20;

  // Basic Jet Pt
  extern const NamedFunc LeadingJetPt;
  extern const NamedFunc SubLeadingJetPt;
  extern const NamedFunc SubSubLeadingJetPt;

  // Funcs to look at number of 3rd jets which are actually B's
  extern const NamedFunc LeadingNonBJetPt;
  extern const NamedFunc LeadingFakeNonBJetPt;
  extern const NamedFunc LeadingRealNonBJetPt;
  extern const NamedFunc LeadingNonBJetPt_med;
  extern const NamedFunc LeadingFakeNonBJetPt_med;
  extern const NamedFunc LeadingRealNonBJetPt_med;

  // Funcs to look at number of B-tags which are fakes
  extern const NamedFunc LeadingBJetPt;
  extern const NamedFunc LeadingRealBJetPt;
  extern const NamedFunc LeadingFakeBJetPt;
  extern const NamedFunc LeadingBJetPt_med;
  extern const NamedFunc LeadingRealBJetPt_med;
  extern const NamedFunc LeadingFakeBJetPt_med;
  extern const NamedFunc SubLeadingBJetPt;
  extern const NamedFunc SubLeadingRealBJetPt;
  extern const NamedFunc SubLeadingFakeBJetPt;

  // Looking at gen level ISR
  extern const NamedFunc genISRPt;
  extern const NamedFunc genISRgenMETdPhi;
  extern const NamedFunc genISRrecoMETdPhi;
  extern const NamedFunc genISRrecoISRdPhi;
  extern const NamedFunc genISRrecoISRDeltaPt;
  extern const NamedFunc truthOrigin3rdJet;
  
}

#endif
