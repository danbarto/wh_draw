#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/wh_functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace WH_Functions;

int main(){
  gErrorIgnoreLevel = 6000;

  // paths, new babies
  string mc2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/s16v3/";
  string mc2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/f17v2/";
  string mc2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/a18v1/";
 
  string data2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2016*.root
  string data2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2017*.root
  string data2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2018*.root
 
//  string signal2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root
//  string signal2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root
//  string signal2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root

  string signal2016_dir_full_350_100 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh350_mLSP100_s16v3/";
  string signal2018_dir_full_350_100 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh350_mLSP100_a18v1/";

  string signal2016_dir_full_750_1 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh750_mLSP1_s16v3/";
  string signal2018_dir_full_750_1 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh750_mLSP1_a18v1/";
 
  string signal2016_dir_fast = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/";
  string signal2017_dir_fast = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/";
  string signal2018_dir_fast = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/";
 
  Palette colors("txt/colors.txt", "default");

  // Columns
  // names | type | color | file 

  //Data
  auto data2016 = Process::MakeShared<Baby_full>("2016 Data", Process::Type::data, colors("data"),{data2016_dir+"slim_data_2016*.root"},"pass&&(HLT_SingleEl==1||HLT_SingleMu==1||HLT_MET_MHT==1)");
  auto data2017 = Process::MakeShared<Baby_full>("2017 Data", Process::Type::data, colors("data"),{data2017_dir+"slim_data_2017*.root"},"pass&&(HLT_SingleEl==1||HLT_SingleMu==1||HLT_MET_MHT==1)");
  auto data2018 = Process::MakeShared<Baby_full>("2018 Data", Process::Type::data, colors("data"),{data2018_dir+"slim_data_2018*.root"},"pass&&(HLT_SingleEl==1||HLT_SingleMu==1||HLT_MET_MHT==1)"&&HasHEMevent==0.&&HasHEMjet==0.);
  auto dataComb = Process::MakeShared<Baby_full>("Combined Data", Process::Type::data, colors("data"),{data2016_dir+"slim_data_2016*.root",data2017_dir+"slim_data_2017*.root",data2018_dir+"slim_data_2018*.root"},"pass&&(HLT_SingleEl==1||HLT_SingleMu==1||HLT_MET_MHT==1)"&&HasHEMevent==0.&&HasHEMjet==0.);
 
  // SM WH 
  auto SM_WH_2016 = Process::MakeShared<Baby_full>("SM WH 2016", Process::Type::background, colors("SMWH"),{mc2016_dir+"slim_WplusH*root",mc2016_dir+"slim_WminusH*root"});
  auto SM_WH_2017 = Process::MakeShared<Baby_full>("SM WH 2017", Process::Type::background, colors("SMWH"),{mc2017_dir+"slim_WplusH*root",mc2017_dir+"slim_WminusH*root"});
  auto SM_WH_2018 = Process::MakeShared<Baby_full>("SM WH 2018", Process::Type::background, colors("SMWH"),{mc2018_dir+"slim_WplusH*root",mc2018_dir+"slim_WminusH*root"});
  auto SM_WH_Comb = Process::MakeShared<Baby_full>("SM WH Combined", Process::Type::background, colors("SMWH"),{mc2016_dir+"slim_WplusH*root",mc2016_dir+"slim_WminusH*root",mc2017_dir+"slim_WplusH*root",mc2017_dir+"slim_WminusH*root",mc2018_dir+"slim_WplusH*root",mc2018_dir+"slim_WminusH*root"});
 
  // all ttbar
  auto tt_2016_all = Process::MakeShared<Baby_full>("t#bar{t} (1,2l) 2016", Process::Type::background, colors("t1tttt"),{mc2016_dir+"slim*TTJets*root"}, "stitch");
  auto tt_2017_all = Process::MakeShared<Baby_full>("t#bar{t} (1,2l) 2017", Process::Type::background, colors("t1tttt"),{mc2017_dir+"slim*TTJets*root"}, "stitch");
  auto tt_2018_all = Process::MakeShared<Baby_full>("t#bar{t} (1,2l) 2018", Process::Type::background, colors("t1tttt"),{mc2018_dir+"slim*TTJets*root"}, "stitch");
  auto tt_comb_all = Process::MakeShared<Baby_full>("t#bar{t} (1,2l) Combined", Process::Type::background, colors("t1tttt"),{mc2016_dir+"slim*TTJets*root",mc2017_dir+"slim*TTJets*root",mc2018_dir+"slim*TTJets*root"}, "stitch");
 
//  //ttbar // separated high met // no stitch
//  auto tt1l_2016 = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2016", Process::Type::background, colors("tt_1l"),{mc2016_dir+"slim*TTJets_1lep_top_s*.root",mc2016_dir+"slim*TTJets_1lep_tbar_s*"});//,"stitch");
//  auto tt2l_2016 = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2016", Process::Type::background, colors("tt_2l"),{mc2016_dir+"slim*TTJets_2lep_s*.root"});//,"stitch");
// 
//  auto tt1l_2017 = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2017", Process::Type::background, colors("tt_1l"),{mc2017_dir+"slim*TTJets_1lep_top_f*.root",mc2017_dir+"slim*TTJets_1lep_tbar_f*"});//,"stitch");
//  auto tt2l_2017 = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2017", Process::Type::background, colors("tt_2l"),{mc2017_dir+"slim*TTJets_2lep_f*.root"});//,"stitch");
// 
//  auto tt1l_2018 = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2018", Process::Type::background, colors("tt_1l"),{mc2018_dir+"slim*TTJets_1lep_top_a*.root",mc2018_dir+"slim*TTJets_1lep_tbar_a*"});//,"stitch");
//  auto tt2l_2018 = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2018", Process::Type::background, colors("tt_2l"),{mc2018_dir+"slim*TTJets_2lep_a*.root"});//,"stitch");
// 
//  auto tt1l_2016_met = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2016 met", Process::Type::background, colors("tt_1l"),{mc2016_dir+"slim*TTJets_1lep_top_*met*.root",mc2016_dir+"slim*TTJets_1lep_tbar_*met*"});//,"stitch");
//  auto tt2l_2016_met = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2016 met", Process::Type::background, colors("tt_2l"),{mc2016_dir+"slim*TTJets_2lep_*met*.root"});//,"stitch");
// 									 
//  auto tt1l_2017_met = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2017 met", Process::Type::background, colors("tt_1l"),{mc2017_dir+"slim*TTJets_1lep_top_*met*.root",mc2017_dir+"slim*TTJets_1lep_tbar_*met*"});//,"stitch");
//  auto tt2l_2017_met = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2017 met", Process::Type::background, colors("tt_2l"),{mc2017_dir+"slim*TTJets_2lep_*met*.root"});//,"stitch");
// 									 
//  auto tt1l_2018_met = Process::MakeShared<Baby_full>("t#bar{t} (1l) 2018 met", Process::Type::background, colors("tt_1l"),{mc2018_dir+"slim*TTJets_1lep_top_*met*.root",mc2018_dir+"slim*TTJets_1lep_tbar_*met*"});//,"stitch");
//  auto tt2l_2018_met = Process::MakeShared<Baby_full>("t#bar{t} (2l) 2018 met", Process::Type::background, colors("tt_2l"),{mc2018_dir+"slim*TTJets_2lep_*met*.root"});//,"stitch");
// 
//  auto tt1l_Comb = Process::MakeShared<Baby_full>("t#bar{t} (1l) Combined", Process::Type::background, colors("tt_1l"),{mc2016_dir+"slim*TTJets_1lep_top_s*.root",mc2016_dir+"slim*TTJets_1lep_tbar_s*",mc2017_dir+"slim*TTJets_1lep_top_f*.root",mc2017_dir+"slim*TTJets_1lep_tbar_f*",mc2018_dir+"slim*TTJets_1lep_top_a*.root",mc2018_dir+"slim*TTJets_1lep_tbar_a*"});//,"stitch");
//  auto tt2l_Comb = Process::MakeShared<Baby_full>("t#bar{t} (2l) Combined", Process::Type::background, colors("tt_2l"),{mc2016_dir+"slim*TTJets_2lep_s*.root",mc2017_dir+"slim*TTJets_2lep_f*.root",mc2018_dir+"slim*TTJets_2lep_a*.root"});//,"stitch");
//  auto tt1l_Comb_met = Process::MakeShared<Baby_full>("t#bar{t} (1l) Combined met", Process::Type::background, colors("tt_1l"),{mc2016_dir+"slim*TTJets_1lep_top_*met*.root",mc2016_dir+"slim*TTJets_1lep_tbar_*met*",mc2017_dir+"slim*TTJets_1lep_top_*met*.root",mc2017_dir+"slim*TTJets_1lep_tbar_*met*",mc2018_dir+"slim*TTJets_1lep_top_*met*.root",mc2018_dir+"slim*TTJets_1lep_tbar_*met*"});//,"stitch");
//  auto tt2l_Comb_met = Process::MakeShared<Baby_full>("t#bar{t} (2l) Combined met", Process::Type::background, colors("tt_2l"),{mc2016_dir+"slim*TTJets_2lep_*met*.root",mc2017_dir+"slim*TTJets_2lep_*met*.root",mc2018_dir+"slim*TTJets_2lep_*met*.root"});//,"stitch");
 
 //WJets
  auto wjets_2016 = Process::MakeShared<Baby_full>("W+jets 2016", Process::Type::background, colors("wjets"), {mc2016_dir+"slim*_W*JetsToLNu_s16v3*",mc2016_dir+"slim*W*Jets_NuPt200_s16v*.root"},"stitch");
  auto wjets_2017 = Process::MakeShared<Baby_full>("W+jets 2017", Process::Type::background, colors("wjets"), {mc2017_dir+"slim*_W*JetsToLNu_f17v2*",mc2017_dir+"slim*W*Jets_NuPt200_f17v2*.root"},"stitch");
  auto wjets_2018 = Process::MakeShared<Baby_full>("W+jets 2018", Process::Type::background, colors("wjets"), {mc2018_dir+"slim_W*JetsToLNu_a18v1*",mc2018_dir+"slim*W*Jets_NuPt200_a18v1*.root"},"stitch");
  auto wjets_Comb = Process::MakeShared<Baby_full>("W+jets Combined", Process::Type::background, colors("wjets"), {mc2016_dir+"slim*_W*JetsToLNu_s16v3*",mc2016_dir+"slim*W*Jets_NuPt200_s16v*.root",mc2017_dir+"slim*_W*JetsToLNu_f17v2*",mc2017_dir+"slim*W*Jets_NuPt200_f17v2*.root",mc2018_dir+"slim_W*JetsToLNu_a18v1*",mc2018_dir+"slim*W*Jets_NuPt200_a18v1*.root"},"stitch");
 
 //Single top
  auto single_t_2016 = Process::MakeShared<Baby_full>("Single t 2016", Process::Type::background, colors("single_t"), {mc2016_dir+"slim*_ST_*.root"});
  auto single_t_2017 = Process::MakeShared<Baby_full>("Single t 2017", Process::Type::background, colors("single_t"), {mc2017_dir+"slim*_ST_*.root"});
  auto single_t_2018 = Process::MakeShared<Baby_full>("Single t 2018", Process::Type::background, colors("single_t"), {mc2018_dir+"slim*_ST_*.root"});
  auto single_t_Comb = Process::MakeShared<Baby_full>("Single t Combined", Process::Type::background, colors("single_t"), {mc2016_dir+"slim*_ST_*.root",mc2017_dir+"slim*_ST_*.root",mc2018_dir+"slim*_ST_*.root"});
 
  //diboson // excluding 2017 TTWZ // excluding all ZZ // putting ZZ back in
  auto diboson_2016 = Process::MakeShared<Baby_full>("Diboson 2016", Process::Type::background, colors("other"),{mc2016_dir+"slim*WW*.root", mc2016_dir+"slim*WZ*.root" , mc2016_dir+"slim*ZZ*L*.root"});
  auto diboson_2017 = Process::MakeShared<Baby_full>("Diboson 2017", Process::Type::background, colors("other"),{mc2017_dir+"slim*WWTo2*.root", mc2017_dir+"slim*_WZ*.root", mc2017_dir+"slim*ZZ*.root"});
  auto diboson_2018 = Process::MakeShared<Baby_full>("Diboson 2018", Process::Type::background, colors("other"),{mc2018_dir+"slim*WW*.root", mc2018_dir+"slim*WZ*.root" , mc2018_dir+"slim*ZZ*.root"});
  auto diboson_Comb = Process::MakeShared<Baby_full>("Diboson Combined", Process::Type::background, colors("other"),{mc2016_dir+"slim*WW*.root", mc2016_dir+"slim*WZ*.root",mc2016_dir+"slim*ZZ*L*.root",mc2017_dir+"slim*WWTo2*.root", mc2017_dir+"slim*_WZ*.root",mc2017_dir+"slim*ZZ*.root",mc2018_dir+"slim*WW*.root", mc2018_dir+"slim*WZ*.root",mc2018_dir+"slim*ZZ*.root"});//, ZZ2016_dir+"slim*ZZTo*L*root",ZZ2017_dir+"slim*ZZ*root",ZZ2018_dir+"slim*ZZ*root"});
 
  // ttV // excluding TTZToLL_
  auto ttV_2016 = Process::MakeShared<Baby_full>("t#bar{t}V 2016", Process::Type::background, colors("ttv"),{mc2016_dir+"slim*_TTWJets*.root", mc2016_dir+"slim*_TTZ*.root"});
  auto ttV_2017 = Process::MakeShared<Baby_full>("t#bar{t}V 2017", Process::Type::background, colors("ttv"), {mc2017_dir+"slim*_TTWJets*.root", mc2017_dir+"slim*_TTZToLLNuNu*.root", mc2017_dir+"slim*_TTZToQQ*.root"});
  auto ttV_2018 = Process::MakeShared<Baby_full>("t#bar{t}V 2018", Process::Type::background, colors("ttv"), {mc2018_dir+"slim*_TTWJets*.root", mc2018_dir+"slim*_TTZ*.root"});
  auto ttV_Comb = Process::MakeShared<Baby_full>("t#bar{t}V Combined", Process::Type::background, colors("ttv"), {mc2016_dir+"slim*_TTWJets*.root", mc2016_dir+"slim*_TTZ*.root",mc2017_dir+"slim*_TTWJets*.root", mc2017_dir+"slim*_TTZToLLNuNu*.root",mc2017_dir+"slim*_TTZToQQ*.root",mc2018_dir+"slim*_TTWJets*.root", mc2018_dir+"slim*_TTZ*.root"});
 
  //signals
  auto signal_2016_350_100_full = Process::MakeShared<Baby_full>("2016 Signal (350,100) Full", Process::Type::signal, colors("t1tttt"),{signal2016_dir_full_350_100+"slim*"},"pass");
  auto signal_2018_350_100_full = Process::MakeShared<Baby_full>("2018 Signal (350,100) Full", Process::Type::signal, colors("t1tttt"),{signal2018_dir_full_350_100+"slim*"},"pass");

  auto signal_2016_750_1_full = Process::MakeShared<Baby_full>("2016 Signal (750,1) Full", Process::Type::signal, colors("t2tttt"),{signal2016_dir_full_750_1+"slim*"},"pass");
  auto signal_2018_750_1_full = Process::MakeShared<Baby_full>("2018 Signal (750,1) Full", Process::Type::signal, colors("t2tttt"),{signal2018_dir_full_750_1+"slim*"},"pass");

  auto signal_2016_350_100_fast = Process::MakeShared<Baby_full>("2016 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2016_dir_fast+"slim*s16v3*"},"pass&&mass_stop==350&&mass_lsp==100");
  auto signal_2017_350_100_fast = Process::MakeShared<Baby_full>("2017 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2017_dir_fast+"slim*f17v2*"},"pass&&mass_stop==350&&mass_lsp==100");
  auto signal_2018_350_100_fast = Process::MakeShared<Baby_full>("2018 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2018_dir_fast+"slim*a18v1*"},"pass&&mass_stop==350&&mass_lsp==100");

  auto signal_2016_750_1_fast = Process::MakeShared<Baby_full>("2016 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2016_dir_fast+"slim*s16v3*"},"pass&&mass_stop==750&&mass_lsp==1");
  auto signal_2017_750_1_fast = Process::MakeShared<Baby_full>("2017 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2017_dir_fast+"slim*f17v2*"},"pass&&mass_stop==750&&mass_lsp==1");
  auto signal_2018_750_1_fast = Process::MakeShared<Baby_full>("2018 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2018_dir_fast+"slim*a18v1*"},"pass&&mass_stop==750&&mass_lsp==1");

//  auto signal_2016_225_75 = Process::MakeShared<Baby_full>("2016 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==225&&mass_lsp==75");
//  auto signal_2017_225_75 = Process::MakeShared<Baby_full>("2017 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==225&&mass_lsp==75");
//  auto signal_2018_225_75 = Process::MakeShared<Baby_full>("2018 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==225&&mass_lsp==75");
//  auto signal_comb_225_75 = Process::MakeShared<Baby_full>("Combined Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==225&&mass_lsp==75");
//  //auto signal_comb_225_75 = Process::MakeShared<Baby_full>("Combined Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&(mass_lsp+250)<=mass_stop");//mass_stop==225&&mass_lsp==75");
// 
//  auto signal_2016_650_300 = Process::MakeShared<Baby_full>("2016 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==650&&mass_lsp==300");
//  auto signal_2017_650_300 = Process::MakeShared<Baby_full>("2017 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==650&&mass_lsp==300");
//  auto signal_2018_650_300 = Process::MakeShared<Baby_full>("2018 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==650&&mass_lsp==300");
//  auto signal_comb_650_300 = Process::MakeShared<Baby_full>("Combined Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==650&&mass_lsp==300");
//    
//  auto signal_2016_700_1 = Process::MakeShared<Baby_full>("2016 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==700&&mass_lsp==1");
//  auto signal_2017_700_1 = Process::MakeShared<Baby_full>("2017 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==700&&mass_lsp==1");
//  auto signal_2018_700_1 = Process::MakeShared<Baby_full>("2018 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==700&&mass_lsp==1");
//  auto signal_comb_700_1 = Process::MakeShared<Baby_full>("Combined Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==700&&mass_lsp==1");
// 
//  auto signal_2016_all = Process::MakeShared<Baby_full>("2016 Signal All", Process::Type::signal, colors("t4tttt"),{signal2016_dir+"slim_*s16v3*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
//  auto signal_2017_all = Process::MakeShared<Baby_full>("2017 Signal All", Process::Type::signal, colors("t4tttt"),{signal2017_dir+"slim_*f17v2*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
//  auto signal_2018_all = Process::MakeShared<Baby_full>("2018 Signal All", Process::Type::signal, colors("t4tttt"),{signal2018_dir+"slim_*a18v1*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
//  auto signal_comb_all = Process::MakeShared<Baby_full>("Combined Signal All", Process::Type::signal, colors("t4tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
// 
//  auto signal_2016_all_mass = Process::MakeShared<Baby_full>("2016 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass");
//  auto signal_2017_all_mass = Process::MakeShared<Baby_full>("2017 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass");
//  auto signal_2018_all_mass = Process::MakeShared<Baby_full>("2018 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass");
//  auto signal_comb_all_mass = Process::MakeShared<Baby_full>("Combined Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass");
//
//  auto signal_2016_all_mass_diff_l250 = Process::MakeShared<Baby_full>("2016 Signal 250 > mStop - mLsp", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&(mass_stop<(250+mass_lsp))");
//  auto signal_2017_all_mass_diff_l250 = Process::MakeShared<Baby_full>("2017 Signal 250 > mStop - mLsp", Process::Type::signal, colors("t1tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&(mass_stop<(250+mass_lsp))");
//  auto signal_2018_all_mass_diff_l250 = Process::MakeShared<Baby_full>("2018 Signal 250 > mStop - mLsp", Process::Type::signal, colors("t1tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&(mass_stop<(250+mass_lsp))");
//  auto signal_comb_all_mass_diff_l250 = Process::MakeShared<Baby_full>("Combined Signal 250 > mStop - mLsp", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&(mass_stop<(250+mass_lsp))");
// 
//  auto signal_2016_all_mass_diff_ge250_le500 = Process::MakeShared<Baby_full>("2016 Signal 500 >= mStop - mLsp >= 250", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&((mass_stop>=(250+mass_lsp))&&(mass_stop<=(500+mass_lsp)))");
//  auto signal_2017_all_mass_diff_ge250_le500 = Process::MakeShared<Baby_full>("2017 Signal 500 >= mStop - mLsp >= 250", Process::Type::signal, colors("t2tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&((mass_stop>=(250+mass_lsp))&&(mass_stop<=(500+mass_lsp)))");
//  auto signal_2018_all_mass_diff_ge250_le500 = Process::MakeShared<Baby_full>("2018 Signal 500 >= mStop - mLsp >= 250", Process::Type::signal, colors("t2tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&((mass_stop>=(250+mass_lsp))&&(mass_stop<=(500+mass_lsp)))");
//  auto signal_comb_all_mass_diff_ge250_le500 = Process::MakeShared<Baby_full>("Combined Signal 500 >= mStop - mLsp >= 250", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&((mass_stop>=(250+mass_lsp))&&(mass_stop<=(500+mass_lsp)))");
// 
//  auto signal_2016_all_mass_diff_g500 = Process::MakeShared<Baby_full>("2016 Signal mStop - mLsp > 500", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&(mass_stop>(500+mass_lsp))");
//  auto signal_2017_all_mass_diff_g500 = Process::MakeShared<Baby_full>("2017 Signal mStop - mLsp > 500", Process::Type::signal, colors("t3tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&(mass_stop>(500+mass_lsp))");
//  auto signal_2018_all_mass_diff_g500 = Process::MakeShared<Baby_full>("2018 Signal mStop - mLsp > 500", Process::Type::signal, colors("t3tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&(mass_stop>(500+mass_lsp))");
//  auto signal_comb_all_mass_diff_g500 = Process::MakeShared<Baby_full>("Combined Signal mStop - mLsp > 500", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&(mass_stop>(500+mass_lsp))");
  
  // Column combinations by year
  vector<shared_ptr<Process> > sample_list_2016        = {data2016,/*tt_2016_all,*//*tt2l_2016,tt1l_2016,tt2l_2016_met,tt1l_2016_met,*//*single_t_2016,wjets_2016,ttV_2016,*/diboson_2016,/*signal_2016_225_75,signal_2016_650_300,signal_2016_700_1,signal_2016_all,signal_2016_all_mass,signal_2016_all_mass_diff_l250,signal_2016_all_mass_diff_ge250_le500,signal_2016_all_mass_diff_g500,*/signal_2016_350_100_full,signal_2016_750_1_full,signal_2016_350_100_fast,signal_2016_750_1_fast};
  vector<shared_ptr<Process> > sample_list_2017        = {data2017,/*tt_2017_all,*//*tt2l_2017,tt1l_2017,tt2l_2017_met,tt1l_2017_met,*//*single_t_2017,wjets_2017,ttV_2017,*/diboson_2017,/*signal_2017_225_75,signal_2017_650_300,signal_2017_700_1,signal_2017_all,signal_2017_all_mass,signal_2017_all_mass_diff_l250,signal_2017_all_mass_diff_ge250_le500,signal_2017_all_mass_diff_g500,*//*signal_2017_350_100_full,signal_2017_750_1_full,*/signal_2017_350_100_fast,signal_2017_750_1_fast};
  vector<shared_ptr<Process> > sample_list_2018        = {data2018,/*tt_2018_all,*//*tt2l_2018,tt1l_2018,tt2l_2018_met,tt1l_2018_met,*//*single_t_2018,wjets_2018,ttV_2018,*/diboson_2018,/*signal_2018_225_75,signal_2018_650_300,signal_2018_700_1,signal_2018_all,signal_2018_all_mass,signal_2018_all_mass_diff_l250,signal_2018_all_mass_diff_ge250_le500,signal_2018_all_mass_diff_g500,*/signal_2018_350_100_full,signal_2018_750_1_full,signal_2018_350_100_fast,signal_2018_750_1_fast};
  vector<shared_ptr<Process> > sample_list_2016_nodata = {         /*tt_2016_all,*//*tt2l_2016,tt1l_2016,tt2l_2016_met,tt1l_2016_met,*//*single_t_2016,wjets_2016,ttV_2016,*/diboson_2016,/*signal_2016_225_75,signal_2016_650_300,signal_2016_700_1,signal_2016_all,signal_2016_all_mass,signal_2016_all_mass_diff_l250,signal_2016_all_mass_diff_ge250_le500,signal_2016_all_mass_diff_g500,*/signal_2016_350_100_full,signal_2016_750_1_full,signal_2016_350_100_fast,signal_2016_750_1_fast};
  vector<shared_ptr<Process> > sample_list_2017_nodata = {         /*tt_2017_all,*//*tt2l_2017,tt1l_2017,tt2l_2017_met,tt1l_2017_met,*//*single_t_2017,wjets_2017,ttV_2017,*/diboson_2017,/*signal_2017_225_75,signal_2017_650_300,signal_2017_700_1,signal_2017_all,signal_2017_all_mass,signal_2017_all_mass_diff_l250,signal_2017_all_mass_diff_ge250_le500,signal_2017_all_mass_diff_g500,*//*signal_2017_350_100_full,signal_2017_750_1_full,*/signal_2017_350_100_fast,signal_2017_750_1_fast};
  vector<shared_ptr<Process> > sample_list_2018_nodata = {         /*tt_2018_all,*//*tt2l_2018,tt1l_2018,tt2l_2018_met,tt1l_2018_met,*//*single_t_2018,wjets_2018,ttV_2018,*/diboson_2018,/*signal_2018_225_75,signal_2018_650_300,signal_2018_700_1,signal_2018_all,signal_2018_all_mass,signal_2018_all_mass_diff_l250,signal_2018_all_mass_diff_ge250_le500,signal_2018_all_mass_diff_g500,*/signal_2018_350_100_full,signal_2018_750_1_full,signal_2018_350_100_fast,signal_2018_750_1_fast};
  // Column combinations all years
//  vector<shared_ptr<Process> > sample_list_comb        = {dataComb, tt_comb_all, /*tt2l_Comb, tt1l_Comb, tt2l_Comb_met, tt1l_Comb_met,*/ single_t_Comb, wjets_Comb, ttV_Comb, diboson_Comb, /*signal_comb_225_75, signal_comb_650_300, signal_comb_700_1, signal_comb_all, signal_comb_all_mass, */signal_comb_all_mass_diff_l250, signal_comb_all_mass_diff_ge250_le500, signal_comb_all_mass_diff_g500};
//  vector<shared_ptr<Process> > sample_list_comb_nodata = {          tt_comb_all, /*tt2l_Comb, tt1l_Comb, tt2l_Comb_met, tt1l_Comb_met,*/ single_t_Comb, wjets_Comb, ttV_Comb, diboson_Comb, /*signal_comb_225_75, signal_comb_650_300, signal_comb_700_1, signal_comb_all, signal_comb_all_mass, */signal_comb_all_mass_diff_l250, signal_comb_all_mass_diff_ge250_le500, signal_comb_all_mass_diff_g500};
 
//  vector<shared_ptr<Process> > dummy = {signal_2016_all_mass};
//  vector<shared_ptr<Process> > signal_sample_list_2016 = {signal_2016_all_mass_diff_l250,signal_2016_all_mass_diff_ge250_le500,signal_2016_all_mass_diff_g500};
//  vector<shared_ptr<Process> > signal_sample_list_2017 = {signal_2017_all_mass_diff_l250,signal_2017_all_mass_diff_ge250_le500,signal_2017_all_mass_diff_g500};
//  vector<shared_ptr<Process> > signal_sample_list_2018 = {signal_2018_all_mass_diff_l250,signal_2018_all_mass_diff_ge250_le500,signal_2018_all_mass_diff_g500};
//  vector<shared_ptr<Process> > signal_sample_list_comb = {signal_comb_all_mass_diff_l250,signal_comb_all_mass_diff_ge250_le500,signal_comb_all_mass_diff_g500};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Overflow(OverflowType::both)
    //    .Stack(StackType::data_norm)
  //    .Stack(StackType::lumi_shapes)
    .Stack(StackType::signal_overlay)
    .RatioMinimum(0.5).RatioMaximum(1.5) //
    .FileExtensions({"pdf","root","png"}); //
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
//  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
//    .Bottom(BottomType::off) //
//    .ShowBackgroundError(false);
//  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
//  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
//  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info/*lin_shapes_info*/};
  PlotOpt log_lumi_info_SR = log_lumi_info().Bottom(BottomType::off);
  PlotOpt lin_lumi_info_SR = lin_lumi_info().Bottom(BottomType::off);
  vector<PlotOpt> all_plot_types_SR = {log_lumi_info_SR, lin_lumi_info_SR/*lin_shapes_info*/};

  // To make plots and tables
  PlotMaker * pm2016 = new PlotMaker();
  PlotMaker * pm2017 = new PlotMaker();
  PlotMaker * pm2018 = new PlotMaker();
  //  PlotMaker * pmComb = new PlotMaker();

  //Preselection
  NamedFunc basic_cut_SR = "pfmet>125&&mct>200&&mbb>90&&mbb<150&&ngoodbtags==2&&nvetoleps==1&&pass&&PassTrackVeto&&PassTauVeto&&hasNano>0"&&WHLeptons==1; 

  vector<NamedFunc> sels_SR =
    {
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nHiggsFatJet250b2>=1,
      
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nHiggsFatJet250b2>=1,

      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nHiggsFatJet250b2>=1,


//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHHiggsFatJet250b2>=1,
//      
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHHiggsFatJet250b2>=1,
//
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHHiggsFatJet250b2>=1,

      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep>150"&&nGenHHiggsFatJet250>=1,
      
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep>150"&&nGenHHiggsFatJet250>=1,

      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep>150"&&nGenHHiggsFatJet250>=1,
    };
  
  vector<NamedFunc> sels_CR = 
    {
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nHiggsFatJet250b2>=1, 
            
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nHiggsFatJet250b2>=1,

      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nFatJet250b2>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nHiggsFatJet250b0>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nHiggsFatJet250b1>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nHiggsFatJet250b2>=1, 


//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHHiggsFatJet250b2>=1, 
//            
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHHiggsFatJet250b2>=1,
//
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHFatJet250b2>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHHiggsFatJet250b0>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHHiggsFatJet250b1>=1,
//      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHHiggsFatJet250b2>=1, 

      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets==2&&mt_met_lep<150"&&nGenHHiggsFatJet250>=1,
            
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets==3&&mt_met_lep<150"&&nGenHHiggsFatJet250>=1,
 
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHFatJet250>=1,
      basic_cut_SR&&"ngoodjets>=2&&ngoodjets<=3&&mt_met_lep<150"&&nGenHHiggsFatJet250>=1,
    };

  vector<NamedFunc> looseSel = 
    {
      "pfmet>125&&nvetoleps==1&&pass&&PassTrackVeto&&PassTauVeto&&hasNano>0"&&WHLeptons==1&&nGenHFatJet250>=1,
      "pfmet>125&&nvetoleps==1&&pass&&PassTrackVeto&&PassTauVeto&&hasNano>0"&&WHLeptons==1&&nGenHHiggsFatJet250>=1
    };
  // Weights
  //  NamedFunc weight_2016  = "w_noBtagSF"*yearWeight;
  NamedFunc weight_other = "weight"*yearWeight;
  // Met bins
//  vector<string> metbins = {"1","pfmet>200"};
  // all plot makers
  vector<PlotMaker *> pms = {pm2016, pm2017, pm2018};//, pmComb};
  // all years
  vector<string> years = {"y2016", "y2017", "y2018"};//, "yComb"};
  //  vector<string> years_nodata = {"y2016", "y2017", "y2018", "yComb"};
  // all samples
  vector<vector<shared_ptr<Process> > > samples_Run2 = {sample_list_2016,
   							sample_list_2017,
   							sample_list_2018};//,
  //   							sample_list_comb};
  vector<vector<shared_ptr<Process> > > samples_Run2_nodata = {sample_list_2016_nodata,
   							       sample_list_2017_nodata,
   							       sample_list_2018_nodata};//,
  //   							       sample_list_comb_nodata};
//  vector<vector<shared_ptr<Process> > > signal_samples_Run2 = //{dummy};
//    {signal_sample_list_2016,
//     signal_sample_list_2017,
//     signal_sample_list_2018,
//     signal_sample_list_comb};
  
//  vector<double> binning = {200,400,600,800,1000};
//  for (uint iyear = 0; iyear < years.size(); iyear++) {
//    //    if (0 != iyear) {break;}
//    for (uint isel = 0; isel < sels_CR.size(); isel++) {
//      pms[iyear]->Push<Hist1D>(Axis(binning, "FatJet_pt[0]", "FatJet p_{T} [GeV]"), sels_SR[isel], samples_Run2_nodata[iyear], all_plot_types_SR).Weight(weight_other).Tag(years[iyear]);
//      pms[iyear]->Push<Hist1D>(Axis(binning, "FatJet_pt[0]", "FatJet p_{T} [GeV]"), sels_CR[isel], samples_Run2[iyear]       , all_plot_types).Weight(weight_other).Tag(years[iyear]);
////      pms[iyear]->Push<Hist1D>(Axis(binning, "FatJet_pt[0]", "FatJet p_{T} [GeV]"), sels_SR[isel], signal_samples_Run2[iyear], all_plot_types).Weight(weight_other).Tag(years[iyear]);
////      pms[iyear]->Push<Hist1D>(Axis(binning, "FatJet_pt[0]", "FatJet p_{T} [GeV]"), sels_CR[isel], signal_samples_Run2[iyear], all_plot_types).Weight(weight_other).Tag(years[iyear]);
//    }
//  }
  vector<double> binning_loose = {200,300,400,500,600,700,800,2000};
  for (uint iyear = 0; iyear < years.size(); iyear++) {
    for (uint isel = 0; isel < looseSel.size(); isel++) {
      pms[iyear]->Push<Hist1D>(Axis(binning_loose, "FatJet_pt[0]", "FatJet p_{T} [GeV]"), looseSel[isel], samples_Run2_nodata[iyear], all_plot_types_SR).Weight(weight_other).Tag(years[iyear]);
    }
  }


  
  string plotPath = "mistag/mT/ANsignal/FatJet_pT";
  double lumi2016 = 35.9;
  pm2016->MakePlots(lumi2016, plotPath);
  double lumi2017 = 41.6;
  pm2017->MakePlots(lumi2017, plotPath);
  double lumi2018 = 59.7;
  pm2018->MakePlots(lumi2018, plotPath);
//  double lumiComb = 137.2;
//  pmComb->MakePlots(lumiComb, plotPath);
}
