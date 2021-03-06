///// table_preds: Generates tables with MC/data yields, bkg predictions
/////              ABCDs are defined in abcd_method, with planecuts (typicaly MET bins),
////               bincuts (typically Nb/Njets bins), and abcdcuts (the cuts for the 4 regions)

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>
#include <string>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "RooStats/NumberCountingUtils.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TVector3.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/abcd_method.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"

using namespace std;

namespace{
  bool only_mc = false;
  bool only_kappa = false;
  bool split_bkg = false;
  bool only_dilepton = false;
  bool do_leptons = false;
  bool do_signal = true;
  bool low_abcd = false;
  bool ichep_nbm = false;
  bool unblind = false;
  bool debug = false;
  bool do_ht = false;
  bool do_correction = false;
  bool table_preview = false;
  TString skim = "standard";
  TString json = "full";
  TString only_method = "";
  TString mc_lumi = "";
  string sys_wgts_file = "txt/sys_weights.cfg";
  string mm_scen = "";
  float lumi=35.9;
  bool quick_test = false;
  pair<string, string> sig_nc = make_pair("2100","100");
  pair<string, string> sig_c = make_pair("1900","1250");
}

string GetScenario(const string &method);

TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
		   vector<vector<float> > yieldsPlane);
void plotKappaMCData(abcd_method &abcd, vector<vector<vector<float> > >  &kappas, 
		     vector<vector<vector<float> > >  &kappas_mm, vector<vector<vector<float> > >  &kmcdat);
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > > &kappas);
vector<vector<float> > findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
				 vector<vector<vector<float> > > &kappas, 
				 vector<vector<vector<float> > > &kappas_mm, 
				 vector<vector<vector<float> > > &kmcdat, 
				 vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
		vector<vector<vector<float> > > &preds);
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg);

void GetOptions(int argc, char *argv[]);

//// Defining st because older ntuples don't have it
const NamedFunc st("st", [](const Baby &b) -> NamedFunc::ScalarType{
    float stvar = b.ht();
    for (const auto &pt: *(b.leps_pt())) stvar += pt; 
    return stvar;
  });

//// Number of spurious muons
const NamedFunc nbadmu("nbadmu", [](const Baby &b) -> NamedFunc::ScalarType{
    int nbad=0;
    for (const auto &bad: *(b.mus_bad())) nbad += bad; 
    return nbad;
  });


int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string ntupletag="";
  // if(only_method!="" && !only_method.Contains("allmet") && !only_method.Contains("onemet")){
  //   if(only_method.Contains("met150")) ntupletag = "_metG150";
  //   else ntupletag = "_metG200";
  // }

  if(mm_scen == ""){
    cout << " ======== Doing all mis-measurement scenarios ======== \n" << endl;
    only_mc = true;
  }else if(mm_scen == "data"){
    cout << " ======== Comparing MC and actual DATA ======== \n" << endl;
  }else if(mm_scen == "no_mismeasurement" || mm_scen == "off" || mm_scen == "mc_as_data"){
    cout << " ======== No mismeasurement applied ======== \n" << endl;
    if(mm_scen == "mc_as_data") only_mc = true;
  }else{
    cout << " ======== Doing mis-measurement scenario " << mm_scen << " ======== \n" << endl;
    only_mc = true;
  }

  const NamedFunc top_jet_mass_cut("top_jet_mass_cut", [](const Baby &b) -> NamedFunc::ScalarType{
    int hadtop_id(0);
    double mass(-1);
    TVector3 hadtop(0,0,0), jet(0,0,0);
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find which top is hadronic (opposite sign as lepton)
      if(b.mc_id()->at(imc) == 11 || b.mc_id()->at(imc) == 13 || b.mc_id()->at(imc) == 15) hadtop_id = 6;
      else if(b.mc_id()->at(imc) == -11 || b.mc_id()->at(imc) == -13 || b.mc_id()->at(imc) == -15) hadtop_id = -6;
    }
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find hadronic top
      if(b.mc_id()->at(imc) == hadtop_id) hadtop.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
    }
    for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++){ // Truth-match jet to hadronic top
      jet.SetPtEtaPhi(b.ak8jets_pt()->at(ijet),b.ak8jets_eta()->at(ijet),b.ak8jets_phi()->at(ijet));
      if(jet.DeltaR(hadtop) < 0.6) mass = b.ak8jets_m()->at(ijet);
    }
    if(mass > 105 && mass < 210) return 1;
    else return 0;
    });
  
  const NamedFunc top_jet_nom_score("top_jet_nom_score", [](const Baby &b) -> NamedFunc::ScalarType{
    int hadtop_id(0);
    double score(-1);
    TVector3 hadtop(0,0,0), jet(0,0,0);
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find which top is hadronic (opposite sign as lepton)
      if(b.mc_id()->at(imc) == 11 || b.mc_id()->at(imc) == 13 || b.mc_id()->at(imc) == 15) hadtop_id = 6;
      else if(b.mc_id()->at(imc) == -11 || b.mc_id()->at(imc) == -13 || b.mc_id()->at(imc) == -15) hadtop_id = -6;
    }
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find hadronic top
      if(b.mc_id()->at(imc) == hadtop_id) hadtop.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
    }
    for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++){ // Truth-match jet to hadronic top
      jet.SetPtEtaPhi(b.ak8jets_pt()->at(ijet),b.ak8jets_eta()->at(ijet),b.ak8jets_phi()->at(ijet));
      if(jet.DeltaR(hadtop) < 0.6) score = b.ak8jets_nom_bin_top()->at(ijet);
    }
    return score;
    });
  
  const NamedFunc top_jet_decor_score("top_jet_decor_score", [](const Baby &b) -> NamedFunc::ScalarType{
    int hadtop_id(0);
    double score(-1);
    TVector3 hadtop(0,0,0), jet(0,0,0);
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find which top is hadronic (opposite sign as lepton)
      if(b.mc_id()->at(imc) == 11 || b.mc_id()->at(imc) == 13 || b.mc_id()->at(imc) == 15) hadtop_id = 6;
      else if(b.mc_id()->at(imc) == -11 || b.mc_id()->at(imc) == -13 || b.mc_id()->at(imc) == -15) hadtop_id = -6;
    }
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find hadronic top
      if(b.mc_id()->at(imc) == hadtop_id) hadtop.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
    }
    for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++){ // Truth-match jet to hadronic top
      jet.SetPtEtaPhi(b.ak8jets_pt()->at(ijet),b.ak8jets_eta()->at(ijet),b.ak8jets_phi()->at(ijet));
      if(jet.DeltaR(hadtop) < 0.6) score = b.ak8jets_decor_bin_top()->at(ijet);
    }
    return score;
    });
  
  const NamedFunc is_300ak8("is_300ak8", [](const Baby &b) -> NamedFunc::ScalarType{
  	bool is300(false);
  	if(b.ak8jets_pt()->size() > 0) {
  		for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
  			if(b.ak8jets_pt()->at(iak8) > 300) is300 = true;
      }	}
  	return is300;
  	});
  
  const NamedFunc ak8_nLoose_nom("ak8_nLoose_nom", [](const Baby &b) -> NamedFunc::ScalarType{
  	int pass(0);
  	if(b.ak8jets_pt()->size() > 0) {
  		for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
  			if(b.ak8jets_pt()->at(iak8) > 300 && b.ak8jets_nom_bin_top()->at(iak8) > 0.1883) pass++;
  	}	}
  	return pass;
  	});
  
  const NamedFunc ak8_nLoose_decor("ak8_nLoose_decor", [](const Baby &b) -> NamedFunc::ScalarType{
  	int pass(0);
  	if(b.ak8jets_pt()->size() > 0) {
  		for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
  			if(b.ak8jets_pt()->at(iak8) > 300 && b.ak8jets_decor_bin_top()->at(iak8) > 0.04738) pass++;
  	}	}
  	return pass;
  	});
  
  vector<string> scenarios = ConfigParser::GetOptSets(sys_wgts_file);
  NamedFunc w = "weight*eff_trig";
  map<string, NamedFunc> weights, corrections;
  auto central = Functions::Variation::central;
  weights.emplace("no_mismeasurement", w);
  corrections.emplace("no_mismeasurement", 1.);
  if(mm_scen == ""){
    for(const auto &scen: scenarios){
      weights.emplace(scen, w*Functions::MismeasurementWeight(sys_wgts_file, scen));
      corrections.emplace(scen, Functions::MismeasurementCorrection(sys_wgts_file, scen, central));
    }
  }else if(mm_scen == "data"){
    scenarios = vector<string>{mm_scen};
  }else if(mm_scen != "no_mismeasurement"){
    scenarios = vector<string>{mm_scen};
    weights.emplace(mm_scen, w*Functions::MismeasurementWeight(sys_wgts_file, mm_scen));
    corrections.emplace(mm_scen, Functions::MismeasurementCorrection(sys_wgts_file, mm_scen, central));
  }

  //// Capybara
//  string foldersig(bfolder+"/cms2r0/babymaker/babies/2016_08_10/T1tttt/merged_mcbase_standard/");
//  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");
//  string folderdata(bfolder+"/cms2r0/babymaker/babies/2017_01_27/data/merged_database_stdnj5/");
  //string folderdata(bfolder+"/cms2r0/babymaker/babies/2017_01_21/data/merged_database_stdnj5/");
  
  //// Bear 
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/");
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/");
  string folderdata("");//bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/");

  // DeepAK8 study
  string ak8_mc_standard_path(bfolder+"/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
  string ak8_tag("");// ("ak8");
  if(ak8_tag != "") {
		foldersig = ak8_mc_standard_path;
		foldermc  = ak8_mc_standard_path;
  	}


  // Old 2015 data
  if(skim.Contains("2015")){
    ntupletag="stdnj5";
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_mcbase_stdnj5/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_stdnj5/";
    if(only_method.Contains("old")) {
      ntupletag="1lht500met200";
      foldermc = bfolder+"/cms2r0/babymaker/babies/2015_11_28/mc/skim_1lht500met200/";
      folderdata = bfolder+"/cms2r0/babymaker/babies/2016_02_04/data/singlelep/combined/skim_1lht500met200/";
    }
  }

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline_s = "mj14>250 && nleps>=1 && met>100 && njets>=5 && st<10000 && pass_ra2_badmu && met/met_calo<5";
  //string baseline_s = "mj14>250 && nleps>=1 && met>100 && njets>=5 && st<10000";
  if(skim.Contains("mj12")) ReplaceAll(baseline_s, "mj14","mj");
  if(skim.Contains("met100")) ReplaceAll(baseline_s, "150","100");

  NamedFunc baseline=baseline_s;
  if(do_ht) baseline = baseline && "ht>500";
  else baseline = baseline && st>500;
  // AK8 Selections
  if     (ak8_tag == "ak8_is300ak8") baseline = baseline && is_300ak8; 
  else if(ak8_tag == "ak8_2Lnom"   ) baseline = baseline && ak8_nLoose_nom   >= 2;
  else if(ak8_tag == "ak8_2Ldecor" ) baseline = baseline && ak8_nLoose_decor >= 2;
  else if(ak8_tag == "ak8_1Lnom"   ) baseline = baseline && ak8_nLoose_nom   >= 1;
  else if(ak8_tag == "ak8_1Ldecor" ) baseline = baseline && ak8_nLoose_decor >= 1;

  auto proc_t1nc = Process::MakeShared<Baby_full>("("+sig_nc.first+","+sig_nc.second+")", Process::Type::signal,     colors("t1tttt"), 
    {foldersig+"*mGluino-"+sig_nc.first+"_mLSP-"+sig_nc.second+"_*.root"},  baseline && "stitch_met");
  auto proc_t1c  = Process::MakeShared<Baby_full>("("+sig_c.first+","+sig_c.second+")",  Process::Type::signal,     kMagenta, 
    {foldersig+"*mGluino-"+sig_c.first+"_mLSP-"+sig_c.second+"_*.root"}, baseline && "stitch_met");
  auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep",    Process::Type::background, colors("tt_1l"),  
    {foldermc+ "*_TTJets*SingleLept*.root"},      baseline && "stitch_met && pass");  
  auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep",    Process::Type::background, colors("tt_2l"),  
    {foldermc+ "*_TTJets*DiLept*.root"},          baseline && "stitch_met && pass");

  // Filling all other processes
  vector<string> vnames_other({
    "_WJetsToLNu",
     "_ST_",
     "_TTW",
     "_TTZ",
     "DYJetsToLL",
    "_ZJet",
     "_ttHTobb_M125_",
     "_TTGJets",
     "_TTTT",
    "_WH_HToBB",
     "_ZH_HToBB",
     "_WWTo",
     "_WZ",
     "_ZZ_"
  });
  // QCD changed name in Capybara (2016_08_10)
  if(skim.Contains("2015")) vnames_other.push_back("QCD_HT");
  else{
    // vnames_other.push_back("QCD_HT*0_Tune");
    // vnames_other.push_back("QCD_HT*Inf_Tune");
  }
  set<string> names_other;
  for(auto name : vnames_other)
    names_other.insert(name = foldermc + "*" + name + "*" + ntupletag + "*.root");
  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    names_other, baseline && "stitch_met && pass");

  //// All MC files to make pseudodata
  set<string> names_allmc = names_other;
  names_allmc.insert(foldermc + "*_TTJets*Lept*" + ntupletag + "*.root");

  string trigs = "trig_ra4";
  if(skim.Contains("2015")) trigs = "(trig[4]||trig[8]||trig[28]||trig[14])";

  // Setting luminosity
  string jsonCuts = "nonblind";
  if(skim.Contains("2015")) lumi = 2.3;
  else if(json=="0p869"){
    lumi = 0.869;
    jsonCuts = "nonblind";
  } else if(json=="2p8"){
    lumi = 2.8;
    jsonCuts = "json2p6";
  } else if(json=="1p5"){
    lumi = 1.5;
    jsonCuts = "json4p0&&!json2p6";
  } else if(json=="4p3"){
    lumi = 4.3;
    jsonCuts = "json4p0";
  } else if(json=="3p4"){
    lumi = 3.4;
    jsonCuts = "json7p65&&!json4p0";
  } else if(json=="7p7"){
    lumi = 7.7;
    jsonCuts = "json7p65";
  } else if(json=="12p9"){
    lumi = 12.9;
    jsonCuts = "json12p9";
  } else if(json=="full"){
    lumi = 35.9;
    jsonCuts = "1";
  }
  if(mc_lumi!="") lumi = mc_lumi.Atof();


  if(only_method.Contains("old")) trigs = "(trig[4]||trig[8])";
  if(!skim.Contains("2015")) trigs += " && "+jsonCuts;

  set<string> names_data({folderdata+"*"+ntupletag+"*.root"});
  if(only_mc){
    names_data = names_allmc;
    if(quick_test) names_data = set<string>({foldermc+"*_TTJets_Tune*"+ntupletag+"*.root"});
    trigs = quick_test ? "1" : "stitch_met"; 
  }
  auto proc_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    names_data,baseline && trigs && "pass");
  //No bad muons: "pass && n_mus_bad==0. && n_mus_bad_dupl==0. && n_mus_bad_trkmu==0."
  
  //// Use this process to make quick plots. Requires being run without split_bkg
  auto proc_bkg = Process::MakeShared<Baby_full>("All_bkg", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets_Tune*"+ntupletag+"*.root"}, baseline && " pass");

  vector<shared_ptr<Process> > all_procs;
  if(!quick_test) all_procs = vector<shared_ptr<Process> >{proc_tt1l, proc_tt2l, proc_other};
  else {
    all_procs = vector<shared_ptr<Process> >{proc_bkg};
    split_bkg = false;
  }
  if (do_signal){
    all_procs.push_back(proc_t1nc);
    all_procs.push_back(proc_t1c);
  }
//   if     (tag == "ak8_is300ak8") all_procs = {proc_t1c_is300ak8, proc_t1nc_is300ak8, proc_tt1l_is300ak8, proc_tt2l_is300ak8};
//   else if(tag == "ak8_2Lnom")    all_procs = {proc_t1c_2Lnom,    proc_t1nc_2Lnom,    proc_tt1l_2Lnom,    proc_tt2l_2Lnom   };
//   else if(tag == "ak8_2Ldecor")  all_procs = {proc_t1c_2Ldecor,  proc_t1nc_2Ldecor,  proc_tt1l_2Ldecor,  proc_tt2l_2Ldecor };

  all_procs.push_back(proc_data);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining basic cuts //////////////////////////////////////////
  // baseline defined above

  ////// MET cuts
  TString c_vvlowmet = "met>100 && met<=150";
  TString c_vlowmet  = "met>150 && met<=200";
  TString c_lowmet   = "met>200 && met<=350";
  TString c_midmet   = "met>350 && met<=500";
  TString c_higmet   = "met>500";

  ////// Nb cuts
  TString c_vlownb = "nbm==0";
  TString c_lownb  = "nbm==1";
  TString c_midnb  = "nbm==2";
  TString c_hignb  = "nbm>=3";

  ////// Njets cuts
  TString c_vlownj = "njets>=4 && njets<=5";
  TString c_lownj = "njets>=6 && njets<=7";
  TString c_hignj = "njets>=8";

  ////// ABCD cuts
  vector<TString> abcdcuts_std  = {"mt<=140 && mj14<=400 &&                                   nj_all_1l",
                                   "mt<=140 && mj14> 400 &&                                   nj_1l",
                                   "mt>140 && mj14<=400 &&                                   nj_all_1l",
                                   "mt>140 && mj14> 400 &&                                   nj_1l"};

  vector<TString> abcdcuts_mtcr  = {"mt<=100 && mj14<=400 &&                                   nj_all_1l",
                                   "mt<=100 && mj14> 400 &&                                   nj_1l",
                                   "mt>100 && mt<=140  && mj14<=400 &&                                   nj_all_1l",
                                   "mt>100 && mt<=140  && mj14> 400 &&                                   nj_1l"};

  vector<TString> abcdcuts_veto = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1 && nj_all_1l",
                                   "mt<=140 && mj14> 400 && nleps==1 && nveto==0 && nbm>=1 && nj_1l",
                                   "mt>140  && mj14<=400 && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  nj_all_1l",
                                   "mt>140  && mj14> 400 && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  nj_1l"};

  vector<TString> abcdcuts_2l   = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1 && nj_all_1l",
                                   "mt<=140 && mj14> 400 && nleps==1 && nveto==0 && nbm>=1 && nj_1l",
                                   "           mj14<=400 && nleps==2             && nbm<=2 && nj_all_2l",
                                   "           mj14> 400 && nleps==2             && nbm<=2 && nj_2l"};

  vector<TString> abcdcuts_2lveto;
  for(size_t ind=0; ind<2; ind++) abcdcuts_2lveto.push_back(abcdcuts_2l[ind]);
  for(size_t ind=2; ind<abcdcuts_2l.size(); ind++){
    abcdcuts_2lveto.push_back("(("+abcdcuts_2l[ind]+") || ("+abcdcuts_veto[ind]+"))");
    abcdcuts_2lveto.back().ReplaceAll("((  ","((");
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD methods //////////////////////////////////////////
  vector<abcd_method> abcds;
  vector<TString> abcdcuts, metcuts, bincuts;
  vector<vector<TString> > vbincuts;
  bool doVBincuts = false;
  PlotMaker pm;

  ///// Running over these methods
  vector<TString> methods_all = {"m2lveto", "m2lonly", "mvetoonly", "signal", "signal_nb1", "signal_nb2",
                                 "m2lvetomet150", "m2lonlymet150", "mvetoonlymet150", "m1lmet150",
				 "m5j", "agg_himet", "agg_mixed", "agg_himult", "agg_1b"};
 
  //vector<TString> methods_std = {"signalmet100onebin", "m5jmet100onebin", 
  //			 "m2lvetoonebin", "nb1l", "njets1lmet100x200", "njets1lmet200x500",
  //			 "njets2lveto", "inclvetoonly"};  
  vector<TString> methods_std = {"signal200", "signal350", "signal500"};

  vector<TString> methods = methods_std;

  if(only_method!="") methods = vector<TString>({only_method+"_highmj_abcd", only_method+"_lowmj_abcd"});
  if(do_leptons){
    for(auto name: methods){
      name += "_el";
      methods.push_back(name);
      name.ReplaceAll("_el", "_mu");
      methods.push_back(name);
      if(name.Contains("2l")){
        name.ReplaceAll("_mu", "_emu");
        methods.push_back(name);
      }
    }
  }

  vector<TString> methods_tmp = methods;
  methods.clear();
  for(const auto &scenario: scenarios){
    for(const auto &method: methods_tmp){
      methods.push_back(method+"_scen_"+scenario.c_str());
    }
  }

  TString njets = "N#lower[-0.1]{_{jets}}";
  TString nbs = "N#lower[-0.1]{_{b}}";
  int firstSigBin;
  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    mm_scen = GetScenario(method.Data());
    
    TString basecuts = "", caption = "", abcd_title;
    doVBincuts = false;
    firstSigBin = -1; //// First MET bin that is a signal bin

    //////// General assignments to all methods
    if(method.Contains("2l") || method.Contains("veto")) {
      metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet};
      if(only_mc) metcuts.push_back(c_higmet);
      if(method.Contains("twomet")) metcuts = vector<TString>{"met>100 && met<=200", "met>200 && met<=500"};
      bincuts = vector<TString>{c_lownj, c_hignj}; // 2l nj cuts automatically lowered in abcd_method
      if(method.Contains("onebin")) bincuts = vector<TString>{"njets>=6"};
      caption = "Dilepton validation regions. D3 and D4 have ";
      abcd_title = "Dilepton";
    } else {
      if(only_dilepton) continue;
      abcdcuts = abcdcuts_std;
      basecuts = "nleps==1 && nveto==0 && nbm>=1";
    }

    
    /////// Methods to check Nb
    if(method.Contains("nb1l")) {
      metcuts = vector<TString>{"met>200&&met<=500&&njets==5", "met>200&&met<=500&&njets>=6"};
      firstSigBin = 1;
      bincuts = vector<TString>{"nbm==1", "nbm==2", "nbm>=3"};
      caption = "Signal search regions + $\\njets=5$";
      abcd_title = "Signal + "+njets+"=5 ("+nbs+" bins)";
    }
    /////// Methods to check Njets
    if(method.Contains("njets1lmet100x200")) {
      metcuts = vector<TString>{"met>100&&met<=150", "met>100&&met<=150","met>150&&met<=200", "met>150&&met<=200"};
      vbincuts = vector<vector<TString> >{{"njets==5"}, {"njets>=6&&njets<=8", "njets>=9"},
					  {"njets==5"}, {"njets>=6&&njets<=8", "njets>=9"}}; 
      doVBincuts = true;
      caption = "Low MET";
      abcd_title = "Low MET ("+njets+" bins)";
    }
    if(method.Contains("njets1lmet200x500")) {
      metcuts = vector<TString>{"met>200&&met<=500", "met>200&&met<=500"};
      vbincuts = vector<vector<TString> >{{"njets==5"}, {"njets>=6&&njets<=8", "njets>=9"}}; 
      doVBincuts = true;
      firstSigBin = 1;
      caption = "Signal search regions + $\\njets=5$";
      abcd_title = "Signal + "+njets+"=5 ";
    }
    if(method.Contains("njets2l")) {
      metcuts = vector<TString>{"met>100&&met<=200", "met>200&&met<=500"};
      bincuts = vector<TString>{"njets>=6&&njets<=8", "njets>=9"}; 
    }
    if(method.Contains("inclvetoonly")) {
      metcuts = vector<TString>{"met>100&&met<=200", "met>200&&met<=500"};
      bincuts = vector<TString>{"njets>=6"}; 
    }
    //////// Dilepton methods
    if(method.Contains("2lonly")) {
      abcdcuts = abcdcuts_2l;
      caption += "two reconstructed leptons";
      abcd_title = "Dilepton (ll)";
    }
    if(method.Contains("2lveto")) {
      abcdcuts = abcdcuts_2lveto;
      caption += "either two reconstructed leptons, or one lepton and one track";
      abcd_title = "Dilepton (ll+lv)";
    }
    if(method.Contains("vetoonly")) {
      abcdcuts = abcdcuts_veto;
      caption += "one lepton and one track";
      abcd_title = "Dilepton (lv)";
    }
    if(method.Contains("2lcombined")) {
      metcuts = vector<TString>{"met>200&&met<=500"};
      bincuts = vector<TString>{"njets>=6"}; // 2l nj cuts automatically lowered in abcd_method
      abcdcuts = abcdcuts_2l;
      caption += "two reconstructed leptons";
    }
    if(method.Contains("2lold")) {
      metcuts = vector<TString>{"met>200&&met<=400"};
      abcdcuts = abcdcuts_2l;
      abcdcuts[0].ReplaceAll("&& nveto==0 ","");
      abcdcuts[1].ReplaceAll("&& nveto==0 ","");
      caption += "two reconstructed leptons";
    }

    if(method.Contains("2lvetocombined")) {
      metcuts = vector<TString>{"met>150&&met<=500"};
      bincuts = vector<TString>{"njets>=6"}; // 2l nj cuts automatically lowered in abcd_method
      abcdcuts = abcdcuts_2lveto;
      caption += "either two reconstructed leptons, or one lepton and one track";
    }

    //////// Single lepton methods, all use the standard ABCD plane and nleps==1&&nveto==0&&nbm>=1
    if(method.Contains("signal")) {
      metcuts = vector<TString>{c_lowmet, c_midmet, c_higmet};
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj,
                                c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj,
                                c_hignb+" && "+c_lownj, c_hignb+" && "+c_hignj};
      caption = "Signal search regions";
      abcd_title = "Signal + low MET";
      firstSigBin = 0;
      if(method.Contains("nb1")) {
        bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj};
        caption += " for $\\nb=1$";
      }
      if(method.Contains("nb2")) {
        bincuts = vector<TString>{c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj,
                                  c_hignb+" && "+c_lownj, c_hignb+" && "+c_hignj};
        caption += " for $\\nb\\geq2$";
      }
      if(method.Contains("allmet")) {
				metcuts = vector<TString>{c_vlowmet, c_lowmet, c_midmet, c_higmet};
				caption = "Signal search regions plus $150<\\met\\leq200$ GeV";
				firstSigBin = 1;
      } // allmetsignal
      if(method.Contains("met100")) {
				metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet, c_higmet};
				caption = "Signal search regions plus $100<\\met\\leq200$ GeV";
				firstSigBin = 2;
      } // allmetsignal
//      if(method.Contains("met350")) {
//				metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, "met>350"};
//				caption = "Signal search regions plus $100<\\met\\leq200$ GeV";
//				firstSigBin = 2;
//      } // allmetsignal
      if(method.Contains("200")) metcuts = vector<TString>{c_lowmet};
      if(method.Contains("350")) metcuts = vector<TString>{c_midmet};
      if(method.Contains("500")) metcuts = vector<TString>{c_higmet};
      if(method.Contains("nb0")) {
        metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet, c_higmet};
        bincuts = vector<TString>{"nbm==0&&njets>=6"};
        basecuts = "nleps==1 && nveto==0";
        caption = "Signal search regions plus $100<\\met\\leq200$ GeV for $\\Nb==0$";
        firstSigBin = 2;
      } // allmetsignal
      if(method.Contains("onemet")) {
        metcuts = vector<TString>{"met>200"};
        caption = "Signal search regions plus $150<\\met\\leq200$ GeV";
        firstSigBin = 0;
      } // allmetsignal
      if(method.Contains("onebin")) bincuts = vector<TString>{"njets>=6"};
    } // signal

    if(method.Contains("mtcr")) {
      abcdcuts = abcdcuts_mtcr;
      metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet, c_higmet};
      bincuts = vector<TString>{"njets>=6"};
      caption = "m_{T} 100-140 control regions";
      abcd_title = "m_{T} CR";
      firstSigBin = 0;
    }

    if(method.Contains("m5j")) {
      metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet};
      bincuts = vector<TString>{"njets==5"};
      caption = "Validation regions with $1\\ell, \\njets=5$";
      abcd_title = njets+" = 5";
      if(method.Contains("met100")) {
        metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet};
        caption = "Validation regions with $1\\ell, \\njets=5$";
      } // allmetsignal
      if(method.Contains("onebin")) bincuts = vector<TString>{"njets==5"};
      if(only_mc) metcuts.push_back(c_higmet);     
    }

    if(method.Contains("m6j")) {
      metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet};
      bincuts = vector<TString>{"njets==6"};
      caption = "Validation regions with $1\\ell, \\njets=6$";
      abcd_title = njets+" = 6";
      if(only_mc) metcuts.push_back(c_higmet);     
    }

    if(method.Contains("m1b")) {
      metcuts = vector<TString>{c_vvlowmet, c_vlowmet, c_lowmet, c_midmet};
      bincuts = vector<TString>{c_lownb+" && "+"njets>=6"};
      caption = "Validation regions with $1\\ell, \\nb=1$";
      abcd_title = "N_{b} = 1";
      if(only_mc) metcuts.push_back(c_higmet);     
    }

    ////// Aggregate regions (single lepton). The nbm, njets integration in R1/R3 is done in abcd_method
    if(method.Contains("agg_himet")) {
      metcuts = vector<TString>{"met>500"};
      bincuts = vector<TString>{"nbm>=3&&njets>=6"};
      caption = "High-\\met aggregate region with $1\\ell$, $\\met>500\\text{ GeV}$, $\\njets\\geq6$, $\\nb\\geq3$";
      firstSigBin = 0;
    }
    if(method.Contains("agg_mixed")) {
      metcuts = vector<TString>{"met>350"};
      bincuts = vector<TString>{"nbm>=2&&njets>=9"};
      caption = "Mixed aggregate region with $1\\ell$, $\\met>350\\text{ GeV}$, $\\njets\\geq9$, $\\nb\\geq2$";
      firstSigBin = 0;
    }
    if(method.Contains("agg_himult")) {
      metcuts = vector<TString>{"met>200"};
      bincuts = vector<TString>{"nbm>=3&&njets>=9"};
      caption = "High-multiplicity aggregate region with $1\\ell$, $\\met>200\\text{ GeV}$, $\\njets\\geq9$, $\\nb\\geq3$";
      firstSigBin = 0;
    }
    if(method.Contains("agg_1b")) {
      metcuts = vector<TString>{"met>500"};
      bincuts = vector<TString>{"nbm>=1&&njets>=9"};
      caption = "Single b-tag aggregate region with $1\\ell$, $\\met>500\\text{ GeV}$, $\\njets\\geq9$, $\\nb\\geq1$";
      firstSigBin = 0;
    }

    //////// MET150 methods
    if(method.Contains("met150")) {
      metcuts = vector<TString>{c_vlowmet};
      caption.ReplaceAll("regions", "region for very low \\met");
    }
    if(method.Contains("m1lmet150")) {
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj,
                                c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj};
      caption = "Single lepton validation region for very low \\met";
    }
    if(skim.Contains("2015")) {
      caption += ". Data taken in 2015";
    }

    //////// Pushing all cuts to then find the yields
    if(doVBincuts) abcds.push_back(abcd_method(method, metcuts, vbincuts, abcdcuts, caption, basecuts, abcd_title));
    else {
      abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts, abcd_title));
      abcds.back().printCuts();
    }
    abcds.back().setFirstSignalBin(firstSigBin);

    if(skim.Contains("mj12")) {
      abcds.back().setMj12();
      abcds.back().caption += ". Using $M_J^{1.2}$";
    }
    if(method.Contains("noint")) abcds.back().setIntNbNj(false);
    if(method.Contains("_el") || method.Contains("_mu") || method.Contains("_emu")) abcds.back().setLeptons();
    if(method.Contains("_el"))  abcds.back().caption += ". Only electrons";
    if(method.Contains("_mu"))  abcds.back().caption += ". Only muons";
    if(method.Contains("_emu")) abcds.back().caption += ". Only $e/\\mu$ pairs in D3 and D4";
    //if(method.Contains("agg_")) abcds.back().int_nbnj = false; // Only changes caption since there is only 1 bin


    vector<TableRow> table_cuts, table_cuts_mm;
    NamedFunc correction = do_correction ? corrections.at(mm_scen) : NamedFunc(1.);
    for(size_t icut=0; icut < abcds.back().allcuts.size(); icut++){
      // Changing b-tag working point
      string totcut = abcds.back().allcuts[icut].Data();
      if(!ichep_nbm) ReplaceAll(totcut, "nbm", "nbm_moriond");
      // NamedFunc totcut="1";
      // if(Contains(totcut_s, "nbm==1")){
      // 	ReplaceAll(totcut_s, "&&nbm==1", "");
      // 	totcut = totcut_s && Functions::nbm_moriond == 1.;
      // }
      // if(Contains(totcut_s, "nbm==2")){
      // 	ReplaceAll(totcut_s, "&&nbm==2", "");
      // 	totcut = totcut_s && Functions::nbm_moriond == 2.;
      // }
      // if(Contains(totcut_s, "nbm>=3")){
      // 	ReplaceAll(totcut_s, "&&nbm>=3", "");
      // 	totcut = totcut_s && Functions::nbm_moriond >= 3.;
      // }
      //// Adding cuts to table for yield calculation
      table_cuts.push_back(TableRow(abcds.back().allcuts[icut].Data(), totcut,
				    0,0,weights.at("no_mismeasurement")*correction));
      if(only_mc) table_cuts_mm.push_back(TableRow(abcds.back().allcuts[icut].Data(), totcut,
						   0,0,weights.at(mm_scen)));
    }
    TString tname = "preds"; tname += iabcd;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, true, false);
    tname += mm_scen;
    if(only_mc) pm.Push<Table>(tname.Data(),  table_cuts_mm, all_procs, true, false);
  } // Loop over ABCD methods

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  bool single_thread = false;
  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true; 
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    mm_scen = GetScenario(methods.at(imethod).Data());
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    Table * yield_table;
    if(only_mc){
      yield_table = static_cast<Table*>(pm.Figures()[imethod*2].get());
      Table * yield_table_mm = static_cast<Table*>(pm.Figures()[imethod*2+1].get());
      allyields.push_back(yield_table_mm->Yield(proc_data.get(), lumi));
    } else {
      yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
      allyields.push_back(yield_table->DataYield());
    }
    allyields.push_back(yield_table->BackgroundYield(lumi));
    if(do_signal){
      allyields.push_back(yield_table->Yield(proc_t1nc.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_t1c.get(), lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt1l.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt2l.get(), lumi));
    }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, kappas_mm, kmcdat, preds;
    vector<vector<float> > yieldsPlane = findPreds(abcds[imethod], allyields, kappas, kappas_mm, kmcdat, preds);

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], allyields, TString(baseline.Name()), kappas, kappas_mm, preds);

    //// Makes table MC/Data yields, kappas, preds, Zbi
    if(!only_kappa) {
      TString fullname = printTable(abcds[imethod], allyields, kappas, preds, yieldsPlane);
      tablenames.push_back(fullname);
    }

    //// Plotting kappa comparison between MC and data
    plotKappaMCData(abcds[imethod], kappas, kappas_mm, kmcdat);

    //// Plotting MC kappa
    plotKappa(abcds[imethod], kappas);

  } // Loop over ABCD methods

  if(!only_kappa){
    //// Printing names of ouput files
    cout<<endl<<"===== Tables to be moved to the AN/PAS/paper:"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++){
      TString name=tablenames[ind]; name.ReplaceAll("fulltable","table");
      cout<<" mv "<<name<<"  ${tables_folder}"<<endl;
    }
    cout<<endl<<"===== Tables that can be compiled"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++)
      cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
    cout<<endl;
  }

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

string GetScenario(const string &method){
  string key = "_scen_";
  return method.substr(method.rfind(key)+key.size());
}

//// Prints table with results
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
// if split_bkg: [2/4] Other, [3/5] tt1l, [4/6] tt2l
TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
		   vector<vector<float> > yieldsPlane){
  //cout<<endl<<"Printing table (significance estimation can take a bit)"<<endl;
  //// Table general parameters
  int digits = 2;
  TString ump = " & ";
  bool do_zbi = true;
  if(!unblind) do_zbi = false;
  size_t Ncol = 6;
  if(do_zbi) Ncol++;
  if(do_signal) Ncol += 2;
  if(split_bkg) Ncol += 3;
  if(only_mc) {
    if(do_signal) Ncol -= 2;
    else Ncol -= 3;
  }
  TString blind_s = "$\\spadesuit$";

  //// Setting output file name
  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  TString lumi_s = RoundNumber(lumi, digits_lumi);
  TString outname = "tables/table_pred_lumi"+lumi_s; outname.ReplaceAll(".","p");
  if(skim.Contains("2015")) outname += "_2015";
  if(skim.Contains("mj12")) outname += "_mj12";
  if(unblind) outname += "_unblind";
  else outname += "_blind";
  if(do_ht) outname += "_ht500";
  if(ichep_nbm) outname += "_ichepnbm";
  outname += "_"+abcd.method+".tex";
  ofstream out(outname);

  //// Printing main table preamble
  if(abcd.method.Contains("signal") && Ncol>7 && !table_preview) out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l ";
  if(do_signal) out << "|cc";
//   if(do_signal) out << "|c";
  if(split_bkg) out << "|ccc";
  out << "|cc|ccc";
  if(do_zbi) out<<"c";
  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ ";
  if(do_signal) out << " & ("+sig_nc.first+","+sig_nc.second+") & ("+sig_c.first+","+sig_c.second+") ";
  if(split_bkg) out << " & Other & $t\\bar{t}(1\\ell)$ & $t\\bar{t}(2\\ell)$ ";
  out << "& $\\kappa$ & MC bkg. & Pred.";
  if(!only_mc) out << "& Obs. & Obs./MC "<<(do_zbi?"& Signi.":"");
  else if(do_signal) out << "& Signi.(NC) & Signi.(C)";
//   else if(do_signal) out << "& S(NC)";
  out << " \\\\ \\hline\\hline\n";

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Printing results////////////////////////////////////////////////
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    out<<endl<< "\\multicolumn{"<<Ncol<<"}{c}{$"<<CodeToLatex(abcd.planecuts[iplane].Data())
       <<"$ (Obs/MC = $"<<RoundNumber(yieldsPlane[iplane][0],2)<<"\\pm"<<RoundNumber(yieldsPlane[iplane][1],2)
       <<"$)}  \\\\ \\hline\n";
    for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
      for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        if(abcd.int_nbnj && iabcd%2==0 && ibin>0) continue;
        if(iabcd==3 && ibin==0) out << "\\hline" << endl;
        //// Printing bin name
        out << (iabcd<=1?"R":abcd.rd_letter) << iabcd+1 << ": ";
        if(iabcd%2==0 && abcd.int_nbnj)
          out << "All "<<(abcd.bincuts[iplane][ibin].Contains("nbm")?"\\nb, ":"")<<"\\njets" ;
        else {
          if(abcd.method.Contains("2lonly") && iabcd>=2) out<<"$"<<CodeToLatex(abcd.lowerNjets(abcd.bincuts[iplane][ibin]).Data())<<"$";
          else if(abcd.method.Contains("2lveto") && iabcd>=2){
            if(abcd.bincuts[iplane][ibin].Contains("6")) out<<"Low \\njets";
            else out<<"High \\njets";
          } else out<<"$"<<CodeToLatex(abcd.bincuts[iplane][ibin].Data())<<"$";
        }
        //// Printing signal yields
        if(do_signal)
          out<<ump<<RoundNumber(allyields[2][index].Yield(), digits)<< ump <<RoundNumber(allyields[3][index].Yield(), digits);
//           out<<ump<<RoundNumber(allyields[2][index].Yield(), digits);
        //// Printing Other, tt1l, tt2l
        if(split_bkg){
          size_t offset = (do_signal?2:0);
          out << ump <<RoundNumber(allyields[offset+2][index].Yield(), digits)
              << ump <<RoundNumber(allyields[offset+3][index].Yield(), digits)
              << ump <<RoundNumber(allyields[offset+4][index].Yield(), digits);
        }
        //// Printing kappa
        out<<ump;
        if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits)
                          << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits)
                          << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits) <<"}$ ";
        //// Printing MC Bkg yields
        out << ump << RoundNumber(allyields[1][index].Yield(), digits)<< ump;
        //// Printing background predictions
        if(iabcd==3) out << "$"    << RoundNumber(preds[iplane][ibin][0], digits)
                         << "^{+"  << RoundNumber(preds[iplane][ibin][1], digits)
                         << "}_{-" << RoundNumber(preds[iplane][ibin][2], digits) <<"}$ ";
        if(!only_mc){
          //// Printing observed events in data and Obs/MC ratio
          if(!unblind && iabcd==3 && abcd.signalplanes[iplane]) out << ump << blind_s<< ump << blind_s;
          else {
            out << ump << RoundNumber(allyields[0][index].Yield(), 0);
            TString ratio_s = "-";
            double Nobs = allyields[0][index].Yield(), Nmc = allyields[1][index].Yield();
            double Eobs = sqrt(Nobs), Emc = allyields[1][index].Uncertainty();
            if(Nobs==0) Eobs=1;
            if(Emc>0){
              double ratio = Nobs/Nmc;
              double Eratio = sqrt(pow(Eobs/Nmc,2) + pow(Nobs*Emc/Nmc/Nmc,2));
              ratio_s = "$"+RoundNumber(ratio, 2)+"\\pm"+RoundNumber(Eratio,2)+"$";
            }
            out << ump << ratio_s;
          }
          //// Printing Zbi significance
          if(do_zbi && iabcd==3) out << ump << Zbi(allyields[0][index].Yield(), preds[iplane][ibin][0], 
						   preds[iplane][ibin][1], preds[iplane][ibin][2]);
        } else {// if not only_mc
          if(iabcd==3 && do_signal) {
            out<<ump<<Zbi(allyields[0][index].Yield()+allyields[2][index].Yield(),preds[iplane][ibin][0],
			  preds[iplane][ibin][1], preds[iplane][ibin][2]);
            out<<ump<<Zbi(allyields[0][index].Yield()+allyields[3][index].Yield(),preds[iplane][ibin][0],
			  preds[iplane][ibin][1], preds[iplane][ibin][2]);
          }
        }
        out << "\\\\ \n";
      } // Loop over bin cuts
    } // Loop over ABCD cuts
    out << "\\hline\\hline\n";
  } // Loop over plane cuts
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //// Printing footer and closing file
  out<< "\\end{tabular}"<<endl;
  if(abcd.method.Contains("signal") && Ncol>7 && !table_preview) out << "}\n"; // For resizebox
  out.close();

  //// Copying header and table to the compilable file
  TString fullname = outname; fullname.ReplaceAll("table_","fulltable_");
  ofstream full(fullname);
  ifstream header("txt/header.tex");
  full<<header.rdbuf();
  header.close();
  if(!abcd.method.Contains("signal")) full << "\\usepackage[landscape]{geometry}\n\n";
  full << "\\begin{document}\n\n";
  if(table_preview){
    full << "\\begin{preview}\n";
  } else {
    full << "\\begin{table}\n\\centering\n";
    full << "\\caption{" << abcd.caption <<".}\\vspace{0.1in}\n\\label{tab:"<<abcd.method<<"}\n";
  }
  ifstream outtab(outname);
  full << outtab.rdbuf();
  outtab.close();
  if(table_preview) full << "\\end{preview}\n\n";
  else full << "\\end{table}\n\n";
  full << "\\end{document}\n";
  full.close();

  return fullname;
} // printTable
//// Estimating significance
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg){
  TString zbi_s;
  if(false){ // Old, bad Zbi
    double Nsig = Nobs-Nbkg;
    double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Eup_bkg/Nbkg);
    if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Eup_bkg);
    if(zbi<0) zbi=0;
    zbi_s = RoundNumber(zbi,1);
    if(zbi_s!="-") zbi_s = "$"+zbi_s+"\\sigma$";
    if(Nsig<=0 || Eup_bkg<=0) zbi_s = "-";
  } else zbi_s = "$"+RoundNumber(Significance(Nobs, Nbkg, Eup_bkg, Edown_bkg),1)+"\\sigma$";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Eup_bkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}

//// Makes kappa plots
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > > &kappas){

  bool label_up = false; //// Putting the MET labels at the bottom
  double markerSize = 1.1;

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(label_up) opts.BottomMargin(0.11);
  if(kappas.size() >= 3) {
    opts.CanvasWidth(1300);
    markerSize = 1.5;
    opts.YTitleOffset(0.6);
    opts.LeftMargin(0.1);
    if(kappas.size()>=5){
      opts.RightMargin(0.03);
      }

    cout<<"kappas.size() is "<<kappas.size()<<endl;
    cout<<"kappas[0].size() is "<<kappas[0].size()<<endl;
    cout<<"kappas[0][0].size() is "<<kappas[0][0].size()<<endl;
  }
  
  setPlotStyle(opts);

  struct kmarker{
    TString cut;
    int color;
    int style;
    vector<float> kappa;
  };
  //// k_ordered has all the kappas grouped in sets of nb cuts (typically, in bins of njets)
  vector<vector<vector<kmarker> > > k_ordered;
  vector<kmarker> ind_bcuts; // nb cuts actually used in the plot
  vector<float> zz; // Zero length vector for the kmarker constructor
  vector<kmarker> bcuts({{"nbm==1",4,20,zz}, {"nbm==2",2,21,zz}, {"nbm>=3",kGreen+3,22,zz}, 
								   {"nbm==0",kMagenta+2,23,zz}, 
								   {"nbl==0",kMagenta+2,23,zz}});

  int nbins = 0; // Total number of njets bins (used in the base histo)
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      TString bincut = abcd.bincuts[iplane][ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);
      bool found=false;
      for(size_t ib=0; ib<bcuts.size(); ib++){
        if(bincut.Contains(bcuts[ib].cut)){
          //// Storing the number of different nb cuts in ind_bcuts
          bool cutfound=false;
          for(size_t indb=0; indb<ind_bcuts.size(); indb++)
            if(ind_bcuts[indb].color == bcuts[ib].color) cutfound = true;
          if(!cutfound) ind_bcuts.push_back(bcuts[ib]);

          //// Cleaning the nb cut from the bincut
          bincut.ReplaceAll(bcuts[ib].cut+"&&","");
          for(size_t ik=0; ik<k_ordered[iplane].size(); ik++){
            //// Adding point to a given njets cut
            if(bincut==k_ordered[iplane][ik][0].cut){
              k_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]});
              found = true;
              break;
            } // if same njets cut
          } // Loop over existing ordered kappas
          //// If it doesn't correspond to any njet cut yet, create a new bin
          if(!found) {
            k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]}}));
            found = true;
            nbins++;
          }
        } // if bincut.Contains(bcuts[ib].cut)
      } // Loop over nb cuts

      //// If it doesn't correspond to any nb cut, create a new bin with default (color in [0], blue)
      if(!found) {
        k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kappas[iplane][ibin]}}));
        nbins++;
        if(ind_bcuts.size()==0) ind_bcuts.push_back(bcuts[0]);
      }
    } // Loop over bin cuts
  } // Loop over plane cuts

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
  if(k_ordered.size()>3) label.SetTextSize(0.04);


  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = 4.;//2.4;
  if(label_up) maxy = 2.6;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.SetYTitle("#kappa");
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<vector<double> > vx(ind_bcuts.size()), vexh(ind_bcuts.size()), vexl(ind_bcuts.size());
  vector<vector<double> > vy(ind_bcuts.size()), veyh(ind_bcuts.size()), veyl(ind_bcuts.size());
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, CodeToRootTex(k_ordered[iplane][ibin][0].cut.Data()).c_str());
      // xval is the x position of the first marker in the group
      double xval = bin, nbs = k_ordered[iplane][ibin].size(), minxb = 0.15, binw = 0;
      // If there is more than one point in the group, it starts minxb to the left of the center of the bin
      // binw is the distance between points in the njets group
      if(nbs>1) {
        xval -= minxb;
        binw = 2*minxb/(nbs-1);
      }
      for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
        //// Finding which TGraph this point goes into by comparing the color of the TGraph and the point
        for(size_t indb=0; indb<ind_bcuts.size(); indb++){
          if(ind_bcuts[indb].color == k_ordered[iplane][ibin][ib].color){
            vx[indb].push_back(xval);
            xval += binw;
            vexl[indb].push_back(0);
            vexh[indb].push_back(0);
            vy[indb].push_back(k_ordered[iplane][ibin][ib].kappa[0]);
            veyh[indb].push_back(k_ordered[iplane][ibin][ib].kappa[1]);
            veyl[indb].push_back(k_ordered[iplane][ibin][ib].kappa[2]);
          }
        } // Loop over nb cuts in ordered TGraphs
      } // Loop over nb cuts in kappa plot
    } // Loop over bin cuts

    // Drawing line separating MET planes
    line.SetLineStyle(2); line.SetLineWidth(2);
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);
    // Drawing MET labels
    if(label_up) label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., maxy-0.1, CodeToRootTex(abcd.planecuts[iplane].Data()).c_str());
    else if(kappas.size()<5) label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.26, CodeToRootTex(abcd.planecuts[iplane].Data()).c_str());
    else label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.45, CodeToRootTex(abcd.planecuts[iplane].Data()).c_str());
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  double legX(opts.LeftMargin()+0.026), legY(1-opts.TopMargin()-0.04), legSingle = 0.05;
  if(label_up) legY = 0.8;
  double legW = 0.22, legH = legSingle*(ind_bcuts.size()+1)/2;
  if(ind_bcuts.size()>3) legH = legSingle*((ind_bcuts.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t indb=0; indb<ind_bcuts.size(); indb++){
    graph[indb] = TGraphAsymmErrors(vx[indb].size(), &(vx[indb][0]), &(vy[indb][0]),
                                    &(vexl[indb][0]), &(vexh[indb][0]), &(veyl[indb][0]), &(veyh[indb][0]));
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(markerSize);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");
    leg.AddEntry(&graph[indb], CodeToRootTex(ind_bcuts[indb].cut.Data()).c_str(), "p");
  } // Loop over TGraphs
  if(ind_bcuts.size()>1) leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}");
  cmslabel.SetTextAlign(31);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");
  cmslabel.DrawLatex(1-opts.RightMargin()-0.155, 1-opts.TopMargin()+0.015,"#scale[0.76]{#font[82]{arXiv:xxxx.xxxxx}}");
  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_" + abcd.method;
  if(do_ht) fname  += "_ht500";
  if(ichep_nbm) fname += "_ichepnbm";
  fname += ".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

}

//// Makes kappa plots comparing MC and data
void plotKappaMCData(abcd_method &abcd, vector<vector<vector<float> > > &kappas, 
		     vector<vector<vector<float> > > &kappas_mm, vector<vector<vector<float> > > &kmcdat){

  bool label_up = false; //// Putting the MET labels at the bottom
  double markerSize = 1.1;

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(label_up) opts.BottomMargin(0.11);
  if(kappas.size() >= 1) { // Used to be 4
    opts.CanvasWidth(1300);
    markerSize = 1.5;
  }
  setPlotStyle(opts);

  struct kmarker{
    TString cut;
    int color;
    int style;
    vector<float> kappa;
  };
  //// k_ordered has all the kappas grouped in sets of nb cuts (typically, in bins of njets)
  vector<vector<vector<kmarker> > > k_ordered, kmd_ordered, k_ordered_mm;
  vector<kmarker> ind_bcuts; // nb cuts actually used in the plot
  vector<float> zz; // Zero length vector for the kmarker constructor
  // vector<kmarker> bcuts({{"nbm==1",2,21,zz}, {"nbm==2",4,20,zz}, {"nbm>=3",kGreen+3,22,zz}, 
  // 								   {"nbm==0",kMagenta+2,23,zz}, 
  // 								   {"nbl==0",kMagenta+2,23,zz}});
  vector<kmarker> bcuts({{"none",2,21,zz}});
   
  int cSignal = kBlue;
  float maxy = 2.4, fYaxis = 1.3;
  int nbins = 0; // Total number of njets bins (used in the base histo)
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    kmd_ordered.push_back(vector<vector<kmarker> >());
    k_ordered_mm.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      if(maxy < fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1])) 
        maxy = fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1]);
      if(maxy < fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1])) 
        maxy = fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1]);
      if(maxy < fYaxis*(kappas_mm[iplane][ibin][0])) 
        maxy = fYaxis*(kappas_mm[iplane][ibin][0]);
      TString bincut = abcd.bincuts[iplane][ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);
      bool found=false;
      for(size_t ib=0; ib<bcuts.size(); ib++){
        if(bincut.Contains(bcuts[ib].cut)){
          //// Storing the number of different nb cuts in ind_bcuts
          bool cutfound=false;
          for(size_t indb=0; indb<ind_bcuts.size(); indb++)
            if(ind_bcuts[indb].color == bcuts[ib].color) cutfound = true;
          if(!cutfound) ind_bcuts.push_back(bcuts[ib]);

          //// Cleaning the nb cut from the bincut
          bincut.ReplaceAll(bcuts[ib].cut+"&&","");
          for(size_t ik=0; ik<k_ordered[iplane].size(); ik++){
            //// Adding point to a given njets cut
            if(bincut==k_ordered[iplane][ik][0].cut){
              k_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]});
              kmd_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kmcdat[iplane][ibin]});
              if(unblind || !abcd.signalplanes[iplane])
                k_ordered_mm[iplane][ik].push_back({bincut, 1, bcuts[ib].style, kappas_mm[iplane][ibin]});
              found = true;
              break;
            } // if same njets cut
          } // Loop over existing ordered kappas
          //// If it doesn't correspond to any njet cut yet, create a new bin
          if(!found) {
            k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]}}));
            kmd_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kmcdat[iplane][ibin]}}));
            if(unblind || !abcd.signalplanes[iplane])
              k_ordered_mm[iplane].push_back(vector<kmarker>({{bincut, 1, bcuts[ib].style, kappas_mm[iplane][ibin]}}));
            found = true;
            nbins++;
          }
        } // if bincut.Contains(bcuts[ib].cut)
      } // Loop over nb cuts

      //// If it doesn't correspond to any nb cut, create a new bin with default (color in [0], blue)
      if(!found) {
        k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kappas[iplane][ibin]}}));
        kmd_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kmcdat[iplane][ibin]}}));
        if(unblind || !abcd.signalplanes[iplane])
          k_ordered_mm[iplane].push_back(vector<kmarker>({{bincut, 1, bcuts[0].style, kappas_mm[iplane][ibin]}}));
        nbins++;
        if(ind_bcuts.size()==0) ind_bcuts.push_back(bcuts[0]);
      }
    } // Loop over bin cuts
  } // Loop over plane cuts

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.045); label.SetTextFont(42); label.SetTextAlign(23);
  if(k_ordered.size()>3) label.SetTextSize(0.04);
  if(k_ordered.size()>4) label.SetTextSize(0.037);
  TLatex klab; klab.SetTextFont(42); klab.SetTextAlign(23);


  float minx = 0.5, maxx = nbins+0.5, miny = 0;
  if(label_up) maxy = 2.6;
  if(maxy > 6) maxy = 6;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  TString ytitle = "#kappa";
  if(mm_scen!="data") ytitle += " (Scen. = "+mm_scen+")";
  histo.SetTitleOffset(0.7,"y");
  histo.SetTitleSize(0.07,"y");
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<vector<double> > vx(ind_bcuts.size()), vexh(ind_bcuts.size()), vexl(ind_bcuts.size());
  vector<vector<double> > vy(ind_bcuts.size()), veyh(ind_bcuts.size()), veyl(ind_bcuts.size());
  vector<vector<double> > vx_mm(ind_bcuts.size()), vexh_mm(ind_bcuts.size()), vexl_mm(ind_bcuts.size());
  vector<vector<double> > vy_mm(ind_bcuts.size()), veyh_mm(ind_bcuts.size()), veyl_mm(ind_bcuts.size());
  vector<vector<double> > vx_kmd(ind_bcuts.size()), vexh_kmd(ind_bcuts.size()), vexl_kmd(ind_bcuts.size());
  vector<vector<double> > vy_kmd(ind_bcuts.size()), veyh_kmd(ind_bcuts.size()), veyl_kmd(ind_bcuts.size());
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, CodeToRootTex(k_ordered[iplane][ibin][0].cut.Data()).c_str());
      // xval is the x position of the first marker in the group
      double xval = bin, nbs = k_ordered[iplane][ibin].size(), minxb = 0.15, binw = 0;
      // If there is more than one point in the group, it starts minxb to the left of the center of the bin
      // binw is the distance between points in the njets group
      if(nbs>1) {
        xval -= minxb;
        binw = 2*minxb/(nbs-1);
      }
      for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
        //// Finding which TGraph this point goes into by comparing the color of the TGraph and the point
        for(size_t indb=0; indb<ind_bcuts.size(); indb++){
          if(ind_bcuts[indb].color == k_ordered[iplane][ibin][ib].color){
            vx[indb].push_back(xval);
            vexl[indb].push_back(0);
            vexh[indb].push_back(0);
            vy[indb].push_back(k_ordered[iplane][ibin][ib].kappa[0]);
            veyh[indb].push_back(k_ordered[iplane][ibin][ib].kappa[1]);
            veyl[indb].push_back(k_ordered[iplane][ibin][ib].kappa[2]);
	    
	    //// MC kappas with data uncertainties
            vx_kmd[indb].push_back(xval);
            vexl_kmd[indb].push_back(0);
            vexh_kmd[indb].push_back(0);
            vy_kmd[indb].push_back(kmd_ordered[iplane][ibin][ib].kappa[0]);
	    float ekmdUp = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[1],2) +
				pow(kmd_ordered[iplane][ibin][ib].kappa[1],2));
	    float ekmdDown = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[2],2) +
				  pow(kmd_ordered[iplane][ibin][ib].kappa[2],2));
            veyh_kmd[indb].push_back(ekmdUp);
            veyl_kmd[indb].push_back(ekmdDown);
	    
	    if(unblind || !abcd.signalplanes[iplane]) {
	      //// Data/pseudodata kappas
	      vx_mm[indb].push_back(xval+0.1);
	      vexl_mm[indb].push_back(0);
	      vexh_mm[indb].push_back(0);
	      vy_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[0]);
	      // veyh_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[1]);
	      // veyl_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[2]);
	      veyh_mm[indb].push_back(0);
	      veyl_mm[indb].push_back(0);
	    }

	    if(unblind || !abcd.signalplanes[iplane]) {
	      //// Printing difference between kappa and kappa_mm
	      float kap = k_ordered[iplane][ibin][ib].kappa[0], kap_mm = k_ordered_mm[iplane][ibin][ib].kappa[0];
	      TString text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-kap)*100,0,kap)+"%";
	      if(abcd.signalplanes[iplane])
          klab.SetTextColor(cSignal);
	      else klab.SetTextColor(1);
	      klab.SetTextSize(0.04);
	      if (!only_mc) klab.DrawLatex(xval, 0.952*maxy, text);
	      //// Printing stat uncertainty of kappa_mm/kappa
	      float kapUp = k_ordered[iplane][ibin][ib].kappa[1], kapDown = k_ordered[iplane][ibin][ib].kappa[2];
	      float kap_mmUp = k_ordered_mm[iplane][ibin][ib].kappa[1];
	      float kap_mmDown = k_ordered_mm[iplane][ibin][ib].kappa[2];
	      float errStat = (kap>kap_mm?sqrt(pow(kapDown,2)+pow(kap_mmUp,2)):sqrt(pow(kapUp,2)+pow(kap_mmDown,2)));
	      text = "#sigma_{st} = "+RoundNumber(errStat*100,0, kap)+"%";
        text = "#sigma_{st} = ^{+"+RoundNumber(ekmdUp*100,0, kap)+"%}_{-"+RoundNumber(ekmdDown*100,0, kap)+"%}";
	      klab.SetTextSize(0.05);
	      klab.DrawLatex(xval, 0.888*maxy, text);
	    } // If unblind || not signal bin

           xval += binw;
          }
        } // Loop over nb cuts in ordered TGraphs
      } // Loop over nb cuts in kappa plot
    } // Loop over bin cuts

    // Drawing line separating MET planes
    line.SetLineStyle(2); line.SetLineWidth(2);
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);
    // Drawing MET labels
    TString metlabel = CodeToRootTex(abcd.planecuts[iplane].Data()) + " GeV";
    if(label_up) label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., maxy-0.1, metlabel);
    else label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.10*maxy, metlabel);
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  TString lumi_s = RoundNumber(lumi, digits_lumi);
  double legX(opts.LeftMargin()+0.005), legY(1-0.03), legSingle = 0.05;
  legX = 0.32;
  if(only_mc) legX = 0.4;
  if(label_up) legY = 0.8;
  double legW = 0.35, legH = legSingle*(ind_bcuts.size()+1)/2;
  if(ind_bcuts.size()>3) legH = legSingle*((ind_bcuts.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.15); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  TGraphAsymmErrors graph_kmd[20]; // There's problems with vectors of TGraphs, so using an array
  TGraphAsymmErrors graph_mm[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t indb=0; indb<ind_bcuts.size(); indb++){
    graph_kmd[indb] = TGraphAsymmErrors(vx_kmd[indb].size(), &(vx_kmd[indb][0]), &(vy_kmd[indb][0]),
					&(vexl_kmd[indb][0]), &(vexh_kmd[indb][0]), &(veyl_kmd[indb][0]), 
					&(veyh_kmd[indb][0]));
    graph_kmd[indb].SetMarkerStyle(ind_bcuts[indb].style); graph_kmd[indb].SetMarkerSize(markerSize);
    graph_kmd[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph_kmd[indb].SetLineColor(1); graph_kmd[indb].SetLineWidth(2);
    graph_kmd[indb].Draw("p0 same");

    graph[indb] = TGraphAsymmErrors(vx[indb].size(), &(vx[indb][0]), &(vy[indb][0]),
                                    &(vexl[indb][0]), &(vexh[indb][0]), &(veyl[indb][0]), &(veyh[indb][0]));
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(markerSize);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");

    graph_mm[indb] = TGraphAsymmErrors(vx_mm[indb].size(), &(vx_mm[indb][0]), &(vy_mm[indb][0]),
				       &(vexl_mm[indb][0]), &(vexh_mm[indb][0]), &(veyl_mm[indb][0]), 
				       &(veyh_mm[indb][0]));
    graph_mm[indb].SetMarkerStyle(20); graph_mm[indb].SetMarkerSize(markerSize*1.2);
    graph_mm[indb].SetMarkerColor(1);
    graph_mm[indb].SetLineColor(1); graph_mm[indb].SetLineWidth(2);
    graph_mm[indb].Draw("p0 same");

    leg.AddEntry(&graph[indb], "MC", "p");
    TString data_s = (mm_scen=="data"||mm_scen=="off"||mm_scen=="no_mismeasurement"?"Data":"Pseudodata");
    leg.AddEntry(&graph_mm[indb], data_s+" "+lumi_s+" fb^{-1} (13 TeV)", "p");
    //leg.AddEntry(&graph[indb], CodeToRootTex(ind_bcuts[indb].cut.Data()).c_str(), "p");

  } // Loop over TGraphs
  leg.Draw();
  //if(ind_bcuts.size()>1) leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}";
  TString cmsSim = "#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Preliminary}}";
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  if(only_mc) cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,cmsSim);
  else cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,cmsPrel);
  cmslabel.SetTextAlign(31);
  //cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");
  cmslabel.SetTextSize(0.053);
  TString title = "#font[42]{"+abcd.title+"}";
  TString newSignal = "#color["; newSignal += cSignal; newSignal += "]{Signal}";
  title.ReplaceAll("Signal", newSignal);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.03, title);

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/dataKappa_" + abcd.method;
  if(do_ht) fname  += "_ht500";
  if(ichep_nbm) fname += "_ichepnbm";
  lumi_s.ReplaceAll(".","p");
  fname += "_lumi"+lumi_s;
  fname += ".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

}

//// Calculating kappa and Total bkg prediction
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
vector<vector<float> > findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
	       vector<vector<vector<float> > > &kmcdat, vector<vector<vector<float> > > &preds){
  // Powers for kappa:   ({R1, R2, D3, R4})
  vector<float> pow_kappa({ 1, -1, -1,  1});
  vector<float> pow_kk({ 1, -1, -1,  1});
  // Powers for TotBkg pred:({R1, R2, D3,  R1, R2, D3, D4})
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});

  float val(1.), valup(1.), valdown(1.);
  vector<vector<float> > yieldsPlane;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    //// Counting yields in plane without double-counting R1/R3 yields when integrated
    GammaParams NdataPlane, NmcPlane;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < 4; iabcd++){
	if( ! (abcd.int_nbnj && ibin>0 && (iabcd==0||iabcd==2)) ){
	  size_t index = abcd.indexBin(iplane, ibin, iabcd);
	  NdataPlane += allyields[0][index];
	  NmcPlane += allyields[1][index];
	}
      } // Loop over ABCD cuts
    } // Loop over bin cuts
    //cout<<"Plane "<<iplane<<": MC is "<<NmcPlane<<", data is "<<NdataPlane<<endl;
    float Nobs = NdataPlane.Yield(), Nmc = NmcPlane.Yield();
    float dataMC = Nobs/Nmc;
    float edataMC = sqrt(pow(sqrt(Nobs)/Nmc,2) + pow(Nobs*NmcPlane.Uncertainty()/Nmc/Nmc,2));
    yieldsPlane.push_back({dataMC, edataMC});
    

    kappas.push_back(vector<vector<float> >());
    kmcdat.push_back(vector<vector<float> >());
    kappas_mm.push_back(vector<vector<float> >());
    preds.push_back(vector<vector<float> >());
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      //// Pushing data yields for predictions
      for(size_t iabcd=0; iabcd < 3; iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[0][index].Yield());
        weights.back().push_back(1.);
      } // Loop over ABCD cuts

      vector<vector<float> > kentries;
      vector<vector<float> > kweights;
      vector<vector<float> > kkentries;
      vector<vector<float> > kkweights;
      vector<vector<float> > kentries_mm;
      vector<vector<float> > kweights_mm;
      //// Pushing MC yields for predictions and kappas
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
	// Renormalizing MC to data
	allyields[1][index] *= dataMC;

        // Yields for predictions
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[1][index].NEffective());
        weights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas
        kentries.push_back(vector<float>());
        kweights.push_back(vector<float>());
        kentries.back().push_back(allyields[1][index].NEffective());
        kweights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas on pseudodata
        kentries_mm.push_back(vector<float>());
        kweights_mm.push_back(vector<float>());
        kentries_mm.back().push_back(allyields[0][index].Yield());
        kweights_mm.back().push_back(1.);
        // Yields for kappas_mc normalized to data
        kkentries.push_back(vector<float>());
        kkweights.push_back(vector<float>());
        kkentries.back().push_back(allyields[1][index].Yield());
        kkweights.back().push_back(1.);

      } // Loop over ABCD cuts

      // Throwing toys to find predictions and uncertainties
      val = calcKappa(entries, weights, pow_totpred, valdown, valup);
      if(valdown<0) valdown = 0;
      preds[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries, kweights, pow_kappa, valdown, valup, false, true);
      if(valdown<0) valdown = 0;
      kappas[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries_mm, kweights_mm, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas_mm[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kkentries, kkweights, pow_kk, valdown, valup);
      if(valdown<0) valdown = 0;
      kmcdat[iplane].push_back(vector<float>({val, valup, valdown}));
    } // Loop over bin cuts
  } // Loop over plane cuts

  return yieldsPlane;
} // findPreds

// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
		vector<vector<vector<float> > > &preds){

  int digits = 3;
  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcd.method<<" ==================="<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<abcd.planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        cout<<"MC: "<<setw(8)<<RoundNumber(allyields[1][index].Yield(),digits)
            <<"  Data: "<<setw(4)<<RoundNumber(allyields[0][index].Yield(), 0)
            <<"  - "<< abcd.allcuts[index]<<endl;
      } // Loop over ABCD cuts
      cout<<"Kappa MC = "<<RoundNumber(kappas[iplane][ibin][0],digits)<<"+"<<RoundNumber(kappas[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(kappas[iplane][ibin][2],digits)
	  <<", Kappa Data = "<<RoundNumber(kappas_mm[iplane][ibin][0],digits)
	  <<"+"<<RoundNumber(kappas_mm[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(kappas_mm[iplane][ibin][2],digits)<<", Prediction = "
          <<RoundNumber(preds[iplane][ibin][0],digits)<<"+"<<RoundNumber(preds[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(preds[iplane][ibin][2],digits)<<endl;
      cout<<endl;
    } // Loop over bin cuts
  } // Loop over plane cuts

} // printDebug

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'}, // Method to run on (if you just want one)
      {"correct",      no_argument, 0, 'c'}, // Apply correction
      {"lumi",   required_argument, 0, 'l'}, // Luminosity to normalize MC with (no data)
      {"skim",   required_argument, 0, 's'}, // Which skim to use: standard, 2015 data
      {"json",   required_argument, 0, 'j'}, // Which JSON to use: 0p815, 2p6, 4p0, 7p65, 12p9
      {"split_bkg",    no_argument, 0, 'b'}, // Prints Other, tt1l, tt2l contributions
      {"no_signal",    no_argument, 0, 'n'}, // Does not print signal columns
      {"do_leptons",   no_argument, 0, 'p'}, // Does tables for e/mu/emu as well
      {"unblind",      no_argument, 0, 'u'}, // Unblinds R4/D4
      {"only_mc",      no_argument, 0, 'o'}, // Uses MC as data for the predictions
      {"only_kappa",   no_argument, 0, 'k'}, // Only plots kappa (no table)
      {"debug",        no_argument, 0, 'd'}, // Debug: prints yields and cuts used
      {"only_dilepton",no_argument, 0, '2'}, // Makes tables only for dilepton tests
      {"MJcut",  required_argument, 0, 'M'}, // MJcut for ABCD plane
      {"ht",           no_argument, 0, 0},   // Cuts on ht>500 instead of st>500
      {"mm",     required_argument, 0, 0},   // Mismeasurment scenario, 0 for data
      {"quick",        no_argument, 0, 0},   // Used inclusive ttbar for quick testing
      {"low_abcd",     no_argument, 0, 0}, // MJcut for ABCD plane
      {"ichep_nbm",    no_argument, 0, 0},   // Use ICHEP b-tagging working point
      {"preview",      no_argument, 0, 0},   // Table preview, no caption
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:cs:j:udbnl:p2ok", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      only_method = optarg;
      break;
    case 'c':
      do_correction = true;
      break;
    case 'l':
      mc_lumi = optarg;
      only_mc = true;
      break;
    case 'k':
      only_kappa = true;
      only_mc = true;
      break;
    case 's':
      skim = optarg;
      break;
    case 'j':
      json = optarg;
      break;
    case 'b':
      split_bkg = true;
      break;
    case 'o':
      only_mc = true;
      break;
    case '2':
      only_dilepton = true;
      break;
    case 'p':
      do_leptons = true;
      break;
    case 'n':
      do_signal = false;
      break;
    case 'u':
      unblind = true;
      break;
    case 'd':
      debug = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "ht"){
        do_ht = true;
      } else if(optname == "mm"){
        mm_scen = optarg;
      }else if(optname == "ichep_nbm"){
        ichep_nbm = true;
      }else if(optname == "preview"){
        table_preview = true;
      }else if(optname == "quick"){
        quick_test = true;
      }else if(optname == "low_abcd"){
        low_abcd = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
	    exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
