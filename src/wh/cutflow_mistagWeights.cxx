
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "core/wh_functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace WH_Functions;

namespace{
  bool single_thread = false;
}


int main(){
  gErrorIgnoreLevel = 6000;

  // paths, new babies
  string mc2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/s16v3/";
  string mc2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/f17v2/";
  string mc2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/a18v1/";

  string data2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2016*.root
  string data2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2017*.root
  string data2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim_data_2018*.root

  string signal2016_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root
  string signal2017_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root
  string signal2018_dir = "/home/users/dspitzba/wh_babies/babies_v33_4_2020_05_27/"; // slim*SMS*.root

  string signal2016_dir_full_350_100 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh350_mLSP100_s16v3/";
  string signal2017_dir_full_350_100 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh350_mLSP100_f17v2/";
  string signal2018_dir_full_350_100 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh350_mLSP100_a18v1/";

  string signal2016_dir_full_750_1 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh750_mLSP1_s16v3/";
  string signal2017_dir_full_750_1 = "/home/users/dspitzba/wh_babies/babies_v33_8_2020_07_16/SMS_TChiWH_mCh750_mLSP1_f17v2/";
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

  //ttbar // separated high met // no stitch
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
  auto signal_2016_225_75 = Process::MakeShared<Baby_full>("2016 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==225&&mass_lsp==75");
  auto signal_2017_225_75 = Process::MakeShared<Baby_full>("2017 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==225&&mass_lsp==75");
  auto signal_2018_225_75 = Process::MakeShared<Baby_full>("2018 Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==225&&mass_lsp==75");
  auto signal_comb_225_75 = Process::MakeShared<Baby_full>("Combined Signal (225,75)", Process::Type::signal, colors("t1tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==225&&mass_lsp==75");

  auto signal_2016_650_300 = Process::MakeShared<Baby_full>("2016 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==650&&mass_lsp==300");
  auto signal_2017_650_300 = Process::MakeShared<Baby_full>("2017 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==650&&mass_lsp==300");
  auto signal_2018_650_300 = Process::MakeShared<Baby_full>("2018 Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==650&&mass_lsp==300");
  auto signal_comb_650_300 = Process::MakeShared<Baby_full>("Combined Signal (650,300)", Process::Type::signal, colors("t2tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==650&&mass_lsp==300");
    
  auto signal_2016_700_1 = Process::MakeShared<Baby_full>("2016 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass&&mass_stop==700&&mass_lsp==1");
  auto signal_2017_700_1 = Process::MakeShared<Baby_full>("2017 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass&&mass_stop==700&&mass_lsp==1");
  auto signal_2018_700_1 = Process::MakeShared<Baby_full>("2018 Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==700&&mass_lsp==1");
  auto signal_comb_700_1 = Process::MakeShared<Baby_full>("Combined Signal (700,1)", Process::Type::signal, colors("t3tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass&&mass_stop==700&&mass_lsp==1");

  auto signal_2016_all = Process::MakeShared<Baby_full>("2016 Signal All", Process::Type::signal, colors("t4tttt"),{signal2016_dir+"slim_*s16v3*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
  auto signal_2017_all = Process::MakeShared<Baby_full>("2017 Signal All", Process::Type::signal, colors("t4tttt"),{signal2017_dir+"slim_*f17v2*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
  auto signal_2018_all = Process::MakeShared<Baby_full>("2018 Signal All", Process::Type::signal, colors("t4tttt"),{signal2018_dir+"slim_*a18v1*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");
  auto signal_comb_all = Process::MakeShared<Baby_full>("Combined Signal All", Process::Type::signal, colors("t4tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"(pass&&mass_stop==225&&mass_lsp==75)||(pass&&mass_stop==650&&mass_lsp==300)||(pass&&mass_stop==700&&mass_lsp==1)");

  auto signal_2016_all_mass = Process::MakeShared<Baby_full>("2016 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2016_dir+"slim_*s16v3*.root"},"pass");
  auto signal_2017_all_mass = Process::MakeShared<Baby_full>("2017 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2017_dir+"slim_*f17v2*.root"},"pass");
  auto signal_2018_all_mass = Process::MakeShared<Baby_full>("2018 Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2018_dir+"slim_*a18v1*.root"},"pass");
  auto signal_comb_all_mass = Process::MakeShared<Baby_full>("Combined Signal All Mass", Process::Type::signal, colors("t5tttt"),{signal2016_dir+"slim_*s16v3*.root",signal2017_dir+"slim_*f17v2*.root",signal2018_dir+"slim_*a18v1*.root"},"pass");

  auto signal_2016_350_100_full = Process::MakeShared<Baby_full>("2016 Signal (350,100) Full", Process::Type::signal, colors("t1tttt"),{signal2016_dir_full_350_100+"slim*"},"pass");
  auto signal_2017_350_100_full = Process::MakeShared<Baby_full>("2017 Signal (350,100) Full", Process::Type::signal, colors("t1tttt"),{signal2017_dir_full_350_100+"slim*"},"pass");
  auto signal_2018_350_100_full = Process::MakeShared<Baby_full>("2018 Signal (350,100) Full", Process::Type::signal, colors("t1tttt"),{signal2018_dir_full_350_100+"slim*"},"pass");

  auto signal_2016_750_1_full = Process::MakeShared<Baby_full>("2016 Signal (750,1) Full", Process::Type::signal, colors("t2tttt"),{signal2016_dir_full_750_1+"slim*"},"pass");
  auto signal_2017_750_1_full = Process::MakeShared<Baby_full>("2017 Signal (750,1) Full", Process::Type::signal, colors("t2tttt"),{signal2017_dir_full_750_1+"slim*"},"pass");
  auto signal_2018_750_1_full = Process::MakeShared<Baby_full>("2018 Signal (750,1) Full", Process::Type::signal, colors("t2tttt"),{signal2018_dir_full_750_1+"slim*"},"pass");

  auto signal_2016_350_100_fast = Process::MakeShared<Baby_full>("2016 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2016_dir_fast+"slim*s16v3*"},"pass&&mass_stop==350&&mass_lsp==100");
  auto signal_2017_350_100_fast = Process::MakeShared<Baby_full>("2017 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2017_dir_fast+"slim*f17v2*"},"pass&&mass_stop==350&&mass_lsp==100");
  auto signal_2018_350_100_fast = Process::MakeShared<Baby_full>("2018 Signal (350,100) Fast", Process::Type::signal, colors("t3tttt"),{signal2018_dir_fast+"slim*a18v1*"},"pass&&mass_stop==350&&mass_lsp==100");

  auto signal_2016_750_1_fast = Process::MakeShared<Baby_full>("2016 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2016_dir_fast+"slim*s16v3*"},"pass&&mass_stop==750&&mass_lsp==1");
  auto signal_2017_750_1_fast = Process::MakeShared<Baby_full>("2017 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2017_dir_fast+"slim*f17v2*"},"pass&&mass_stop==750&&mass_lsp==1");
  auto signal_2018_750_1_fast = Process::MakeShared<Baby_full>("2018 Signal (750,1) Fast", Process::Type::signal, colors("t4tttt"),{signal2018_dir_fast+"slim*a18v1*"},"pass&&mass_stop==750&&mass_lsp==1");
  
  // Column combinations by year
  vector<shared_ptr<Process> > sample_list_2016        = {data2016,tt_2016_all,/*tt2l_2016,tt1l_2016,tt2l_2016_met,tt1l_2016_met,*/single_t_2016,wjets_2016,ttV_2016,SM_WH_2016,diboson_2016,signal_2016_225_75,signal_2016_650_300,signal_2016_700_1};//,signal_2016_all,signal_2016_all_mass};
  vector<shared_ptr<Process> > sample_list_2017        = {data2017,tt_2017_all,/*tt2l_2017,tt1l_2017,tt2l_2017_met,tt1l_2017_met,*/single_t_2017,wjets_2017,ttV_2017,SM_WH_2017,diboson_2017,signal_2017_225_75,signal_2017_650_300,signal_2017_700_1};//,signal_2017_all,signal_2017_all_mass};
  vector<shared_ptr<Process> > sample_list_2018        = {data2018,tt_2018_all,/*tt2l_2018,tt1l_2018,tt2l_2018_met,tt1l_2018_met,*/single_t_2018,wjets_2018,ttV_2018,SM_WH_2018,diboson_2018,signal_2018_225_75,signal_2018_650_300,signal_2018_700_1};//,signal_2018_all,signal_2018_all_mass};
  vector<shared_ptr<Process> > sample_list_2016_nodata = {         tt_2016_all,/*tt2l_2016,tt1l_2016,tt2l_2016_met,tt1l_2016_met,*/single_t_2016,wjets_2016,ttV_2016,SM_WH_2016,diboson_2016};//,signal_2016_225_75,signal_2016_650_300,signal_2016_700_1};//,signal_2016_all,signal_2016_all_mass};
  vector<shared_ptr<Process> > sample_list_2017_nodata = {         tt_2017_all,/*tt2l_2017,tt1l_2017,tt2l_2017_met,tt1l_2017_met,*/single_t_2017,wjets_2017,ttV_2017,SM_WH_2017,diboson_2017};//,signal_2017_225_75,signal_2017_650_300,signal_2017_700_1};//,signal_2017_all,signal_2017_all_mass};
  vector<shared_ptr<Process> > sample_list_2018_nodata = {         tt_2018_all,/*tt2l_2018,tt1l_2018,tt2l_2018_met,tt1l_2018_met,*/single_t_2018,wjets_2018,ttV_2018,SM_WH_2018,diboson_2018};//,signal_2018_225_75,signal_2018_650_300,signal_2018_700_1};//,signal_2018_all,signal_2018_all_mass};
  // Column combinations all years
  vector<shared_ptr<Process> > sample_list_comb        = {dataComb, tt_comb_all, /*tt2l_Comb, tt1l_Comb, tt2l_Comb_met, tt1l_Comb_met,*/ single_t_Comb, wjets_Comb, ttV_Comb, SM_WH_Comb, diboson_Comb, signal_comb_225_75,signal_comb_650_300,signal_comb_700_1};//,signal_comb_all,signal_comb_all_mass};
  vector<shared_ptr<Process> > sample_list_comb_nodata = {          tt_comb_all, /*tt2l_Comb, tt1l_Comb, tt2l_Comb_met, tt1l_Comb_met,*/ single_t_Comb, wjets_Comb, ttV_Comb, SM_WH_Comb, diboson_Comb};//, signal_comb_225_75,signal_comb_650_300,signal_comb_700_1};//,signal_comb_all,signal_comb_all_mass};

  vector<shared_ptr<Process> > signal_sample_list_2016 = {         diboson_2016,signal_2016_350_100_full,signal_2016_750_1_full,signal_2016_350_100_fast,signal_2016_750_1_fast};
  vector<shared_ptr<Process> > signal_sample_list_2017 = {         diboson_2017,signal_2017_350_100_full,signal_2017_750_1_full,signal_2017_350_100_fast,signal_2017_750_1_fast};
  vector<shared_ptr<Process> > signal_sample_list_2018 = {         diboson_2018,signal_2018_350_100_full,signal_2018_750_1_full,signal_2018_350_100_fast,signal_2018_750_1_fast};


  // To make tables and plots
  PlotMaker pm;
  PlotMaker pm2016;
  PlotMaker pm2017;
  PlotMaker pm2018;

  NamedFunc basic_cut = "pass&&nvetoleps==1&&PassTrackVeto&&PassTauVeto&&hasNano>0"&&WHLeptons==1;
 
  // Making tables. Adding rows. 
//  Table & cutflow_SR_2016_MistagSF = pm2016.Push<Table>("cutflow_SR_2016_MistagSF", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*higgsMistagSF), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*higgsMistagSF),  
//	}, sample_list_2016_nodata, false);
//
//  Table & cutflow_SR_2017_MistagSF = pm2017.Push<Table>("cutflow_SR_2017_MistagSF", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*higgsMistagSF), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*higgsMistagSF),  
//	}, sample_list_2017_nodata, false);
//
//  Table & cutflow_SR_2018_MistagSF = pm2018.Push<Table>("cutflow_SR_2018_MistagSF", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*higgsMistagSF), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*higgsMistagSF),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*higgsMistagSF),  
//	}, sample_list_2018_nodata, false);
//
//  Table & cutflow_SR_comb_MistagSF = pm.Push<Table>("cutflow_SR_comb_MistagSF", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"*higgsMistagSF*yearWeight),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF*yearWeight),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*higgsMistagSF*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*higgsMistagSF*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*higgsMistagSF*yearWeight), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*higgsMistagSF*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*higgsMistagSF*yearWeight),  
//	}, sample_list_comb_nodata, false);
//
//  Table & cutflow_SR_2016 = pm2016.Push<Table>("cutflow_SR_2016", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
//	}, sample_list_2016_nodata, false);
//
//  Table & cutflow_SR_2017 = pm2017.Push<Table>("cutflow_SR_2017", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
//	}, sample_list_2017_nodata, false);
//
//  Table & cutflow_SR_2018 = pm2018.Push<Table>("cutflow_SR_2018", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
//	}, sample_list_2018_nodata, false);
//
//  Table & cutflow_SR_comb = pm.Push<Table>("cutflow_SR_comb", vector<TableRow>{
//      TableRow("", basic_cut, 0, 0, "weight"*yearWeight),
//	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*yearWeight),
//	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*yearWeight), 
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*yearWeight),
//	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*yearWeight),  
//	}, sample_list_comb_nodata, false);






  Table & signal_cutflow_SR_2016_MistagSF = pm2016.Push<Table>("signal_cutflow_SR_2016_MistagSF", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSF*1.0),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*1.0),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*1.0),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSF*1.0),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSF*1.0), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSF*1.0),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSF*1.0),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSF*1.0),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSF*1.0),  
	}, signal_sample_list_2016, false);

  Table & signal_cutflow_SR_2017_MistagSF = pm2017.Push<Table>("signal_cutflow_SR_2017_MistagSF", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSF*0.99),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*0.99),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*0.99),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSF*0.99),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSF*0.99), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSF*0.99),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSF*0.99),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSF*0.99),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSF*0.99),  
	}, signal_sample_list_2017, false);

  Table & signal_cutflow_SR_2018_MistagSF = pm2018.Push<Table>("signal_cutflow_SR_2018_MistagSF", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSF*0.99),//0.91),  
	}, signal_sample_list_2018, false);

  Table & signal_cutflow_SR_2016_MistagSFUp = pm2016.Push<Table>("signal_cutflow_SR_2016_MistagSFUp", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFUp*0.97),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFUp*0.97), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFUp*0.97),  
	}, signal_sample_list_2016, false);

  Table & signal_cutflow_SR_2017_MistagSFUp = pm2017.Push<Table>("signal_cutflow_SR_2017_MistagSFUp", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFUp*0.93),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFUp*0.93), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),  
	}, signal_sample_list_2017, false);

  Table & signal_cutflow_SR_2018_MistagSFUp = pm2018.Push<Table>("signal_cutflow_SR_2018_MistagSFUp", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFUp*0.93),//0.73),//*0.85),//),//*1.18),  
	}, signal_sample_list_2018, false);

  Table & signal_cutflow_SR_2016_MistagSFDown = pm2016.Push<Table>("signal_cutflow_SR_2016_MistagSFDown", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFDown*1.02), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),   
	}, signal_sample_list_2016, false);

  Table & signal_cutflow_SR_2017_MistagSFDown = pm2017.Push<Table>("signal_cutflow_SR_2017_MistagSFDown", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFDown*1.02), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFDown*1.02),  
	}, signal_sample_list_2017, false);

  Table & signal_cutflow_SR_2018_MistagSFDown = pm2018.Push<Table>("signal_cutflow_SR_2018_MistagSFDown", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"*signalHiggsMistagSFDown*1.04),//0.99),//),//*1.01),  
	}, signal_sample_list_2018, false);

  Table & signal_cutflow_SR_2016 = pm2016.Push<Table>("signal_cutflow_SR_2016", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"),  
	}, signal_sample_list_2016, false);

  Table & signal_cutflow_SR_2017 = pm2017.Push<Table>("signal_cutflow_SR_2017", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"),  
	}, signal_sample_list_2017, false);

  Table & signal_cutflow_SR_2018 = pm2018.Push<Table>("signal_cutflow_SR_2018", vector<TableRow>{
      TableRow("", basic_cut, 0, 0, "weight"),
	TableRow("2-3 jets", basic_cut&&"ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200", 0, 0, "weight"), 
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150", 0, 0, "weight"),
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, no Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs==0", 0, 0, "weight"),  
	TableRow("2-3 jets, 2 b tags, 125$<$MET, 200$<$mct, 150$<$mt, 90$<$mbb$<$150, at least one Higgs", basic_cut&&"ngoodbtags==2&&ngoodjets>=2&&ngoodjets<=3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&nHiggs>=1", 0, 0, "weight"),  
	}, signal_sample_list_2018, false);



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(single_thread) pm.multithreaded_ = false;
  // luminosity
  //  double lumi = 137.2;
  double lumi2016 = 35.9;
  double lumi2017 = 41.5;
  double lumi2018 = 60.0;  
  //  pm.MakePlots(lumi, "cutflows/higgsMistagSF");
  pm2016.MakePlots(lumi2016, "cutflows/higgsMistagSF");
  pm2017.MakePlots(lumi2017, "cutflows/higgsMistagSF");
  pm2018.MakePlots(lumi2018, "cutflows/higgsMistagSF");

//  vector<GammaParams> yields_2016 = cutflow_SR_2016.BackgroundYield(lumi2016);
//  for(const auto &yield: yields_2016) {
//    cout << "my yield 2016: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_2017 = cutflow_SR_2017.BackgroundYield(lumi2017);
//  for(const auto &yield: yields_2017) {
//    cout << "my yield 2017: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_2018 = cutflow_SR_2018.BackgroundYield(lumi2018);
//  for(const auto &yield: yields_2018) {
//    cout << "my yield 2018: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_comb = cutflow_SR_comb.BackgroundYield(lumi);
//  for(const auto &yield: yields_comb) {
//    cout << "my yield comb: " << yield << endl;
//  }
//  
//  vector<GammaParams> yields_SR_2016_MistagSF = cutflow_SR_2016_MistagSF.BackgroundYield(lumi2016);
//  for(const auto &yield: yields_SR_2016_MistagSF) {
//    cout << "my yield 2016: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_SR_2017_MistagSF = cutflow_SR_2017_MistagSF.BackgroundYield(lumi2017);
//  for(const auto &yield: yields_SR_2017_MistagSF) {
//    cout << "my yield 2017: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_SR_2018_MistagSF = cutflow_SR_2018_MistagSF.BackgroundYield(lumi2018);
//  for(const auto &yield: yields_SR_2018_MistagSF) {
//    cout << "my yield 2018: " << yield << endl;
//  }
//  cout << "---" << endl;
//  vector<GammaParams> yields_SR_comb_MistagSF = cutflow_SR_comb_MistagSF.BackgroundYield(lumi);
//  for(const auto &yield: yields_SR_comb_MistagSF) {
//    cout << "my yield comb: " << yield << endl;
//  }


  vector<GammaParams> signal_yields_2016 = signal_cutflow_SR_2016.BackgroundYield(lumi2016);
  for(const auto &yield: signal_yields_2016) {
    cout << "my yield 2016: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_2017 = signal_cutflow_SR_2017.BackgroundYield(lumi2017);
  for(const auto &yield: signal_yields_2017) {
    cout << "my yield 2017: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_2018 = signal_cutflow_SR_2018.BackgroundYield(lumi2018);
  for(const auto &yield: signal_yields_2018) {
    cout << "my yield 2018: " << yield << endl;
  }
  
  vector<GammaParams> signal_yields_SR_2016_MistagSF = signal_cutflow_SR_2016_MistagSF.BackgroundYield(lumi2016);
  for(const auto &yield: signal_yields_SR_2016_MistagSF) {
    cout << "my yield 2016: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2017_MistagSF = signal_cutflow_SR_2017_MistagSF.BackgroundYield(lumi2017);
  for(const auto &yield: signal_yields_SR_2017_MistagSF) {
    cout << "my yield 2017: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2018_MistagSF = signal_cutflow_SR_2018_MistagSF.BackgroundYield(lumi2018);
  for(const auto &yield: signal_yields_SR_2018_MistagSF) {
    cout << "my yield 2018: " << yield << endl;
  }

  vector<GammaParams> signal_yields_SR_2016_MistagSFUp = signal_cutflow_SR_2016_MistagSFUp.BackgroundYield(lumi2016);
  for(const auto &yield: signal_yields_SR_2016_MistagSFUp) {
    cout << "my yield 2016: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2017_MistagSFUp = signal_cutflow_SR_2017_MistagSFUp.BackgroundYield(lumi2017);
  for(const auto &yield: signal_yields_SR_2017_MistagSFUp) {
    cout << "my yield 2017: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2018_MistagSFUp = signal_cutflow_SR_2018_MistagSFUp.BackgroundYield(lumi2018);
  for(const auto &yield: signal_yields_SR_2018_MistagSFUp) {
    cout << "my yield 2018: " << yield << endl;
  }

  vector<GammaParams> signal_yields_SR_2016_MistagSFDown = signal_cutflow_SR_2016_MistagSFDown.BackgroundYield(lumi2016);
  for(const auto &yield: signal_yields_SR_2016_MistagSFDown) {
    cout << "my yield 2016: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2017_MistagSFDown = signal_cutflow_SR_2017_MistagSFDown.BackgroundYield(lumi2017);
  for(const auto &yield: signal_yields_SR_2017_MistagSFDown) {
    cout << "my yield 2017: " << yield << endl;
  }
  cout << "---" << endl;
  vector<GammaParams> signal_yields_SR_2018_MistagSFDown = signal_cutflow_SR_2018_MistagSFDown.BackgroundYield(lumi2018);
  for(const auto &yield: signal_yields_SR_2018_MistagSFDown) {
    cout << "my yield 2018: " << yield << endl;
  }
}
