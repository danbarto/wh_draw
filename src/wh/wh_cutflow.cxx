
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
  

  double lumi = 35.9;

  string data_dir = "/home/users/rheller/wh_babies/babies_v31_1_2019_04_03/";
  string mc_dir = "/home/users/rheller/wh_babies/babies_v30_9_2019_04_03/";
  string signal_dir = "/home/users/rheller/wh_babies/babies_signal_2019_04_03/";

  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    signal_dir = "/net/cms29";
  }
 
  Palette colors("txt/colors.txt", "default");

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, colors("data"),
    {data_dir+"*data_2016*.root"},"pass&&(HLT_SingleEl==1||HLT_SingleMu==1)");
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_dir+"*TTJets_1lep*ext1*.root"});
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_dir+"*TTJets_2lep*ext1*.root"});
  auto wjets_low_nu = Process::MakeShared<Baby_full>("W+jets, nu pT $<$ 200", Process::Type::background, colors("qcd"),
    {mc_dir+"*slim_W*JetsToLNu_s16v3*.root"},NHighPtNu==0.);
  auto wjets_high_nu = Process::MakeShared<Baby_full>("W+jets, nu pT $>=$ 200", Process::Type::background, colors("wjets"),
    {mc_dir+"*slim_W*Jets_NuPt200*.root"});
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_dir+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_dir+"*_TTWJets*.root", mc_dir+"*_TTZ*.root"});
  auto diboson = Process::MakeShared<Baby_full>("Diboson", Process::Type::background, colors("other"),
    {mc_dir+"*WW*.root", mc_dir+"*WZ*.root",mc_dir+"*ZZ*.root"});

   auto tchiwh_225_75 = Process::MakeShared<Baby_full>("TChiWH(225,75)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==225&&mass_lsp==75");

   auto tchiwh_250_1 = Process::MakeShared<Baby_full>("TChiWH(250,1)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==250&&mass_lsp==1");

   auto tchiwh_350_100 = Process::MakeShared<Baby_full>("TChiWH(350,100)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==350&&mass_lsp==100");

  auto tchiwh_nc = Process::MakeShared<Baby_full>("TChiWH(500,1)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==500&&mass_lsp==1");
  tchiwh_nc->SetMarkerStyle(21);
  tchiwh_nc->SetMarkerSize(0.9);
  auto tchiwh_c = Process::MakeShared<Baby_full>("TChiWH(500,125)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==500&&mass_lsp==125");
  auto tchiwh_700_1 = Process::MakeShared<Baby_full>("TChiWH(700,1)", Process::Type::signal, colors("t1tttt"),
    {signal_dir+"*SMS-TChiWH_WToLNu_HToBB_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*.root"},"mass_stop==700&&mass_lsp==1");


  vector<shared_ptr<Process> > sample_list = {data,ttv,single_t,diboson,wjets_low_nu, wjets_high_nu,tt1l,tt2l,tchiwh_225_75,tchiwh_250_1,tchiwh_350_100, tchiwh_nc, tchiwh_c};
  vector<shared_ptr<Process> > short_sample_list = {wjets_low_nu,wjets_high_nu,tt1l,tt2l,tchiwh_225_75,tchiwh_250_1,tchiwh_350_100, tchiwh_nc, tchiwh_c,tchiwh_700_1};
  vector<shared_ptr<Process> > new_sample_list = {single_t,wjets_low_nu, wjets_high_nu,tt1l,tt2l,tchiwh_225_75,tchiwh_250_1,tchiwh_350_100, tchiwh_nc, tchiwh_c,tchiwh_700_1};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
                                    log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style2D().Stack(StackType::data_norm).Title(TitleType::preliminary)};
  vector<PlotOpt> bkg_pts = {style2D().Stack(StackType::lumi_shapes).Title(TitleType::info)};
  
  PlotMaker pm;



  Table & cutflow = pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("1 good lep,",WHLeptons>=1),
      TableRow("2nd lep veto",WHLeptons==1),
      TableRow("Pass track veto", WHLeptons==1&&"PassTrackVeto"),
      TableRow("Pass tau veto", WHLeptons==1&&"PassTrackVeto&&PassTauVeto"),
      TableRow("==2jets", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2"),
      TableRow("one loose, one med", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasMedLooseCSV>0.),
      TableRow("==2btags", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2"&& HasMedLooseCSV>0.),
      TableRow("M$_{\\rm b\\bar{b}}$ window", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150"&& HasMedLooseCSV>0.),
      TableRow("M$_{\\rm CT}>170$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170"&& HasMedLooseCSV>0.),
      TableRow("E$_{\\rm T}^{\\rm miss}>125$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125"&& HasMedLooseCSV>0.),
      TableRow("M$_{\\rm T}>150$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&mt_met_lep>150"&& HasMedLooseCSV>0.),
      TableRow("E$_{\\rm T}^{\\rm miss}>200$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>200&&mt_met_lep>150"&& HasMedLooseCSV>0.)


        }, sample_list,false);

  Table & cutflow_no_weight = pm.Push<Table>("cutflow_raw", vector<TableRow>{
      TableRow("1 good lep,",WHLeptons>=1,0,0,"1"),
      TableRow("2nd lep veto",WHLeptons==1,0,0,"1"),
      TableRow("Pass track veto", WHLeptons==1&&"PassTrackVeto",0,0,"1"),
      TableRow("Pass tau veto", WHLeptons==1&&"PassTrackVeto&&PassTauVeto",0,0,"1"),
      TableRow("==2jets", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2",0,0,"1"),
      TableRow("one loose, one med", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasMedLooseCSV>0.,0,0,"1"),
      TableRow("==2btags", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2"&& HasMedLooseCSV>0.,0,0,"1"),
      TableRow("M$_{\\rm b\\bar{b}}$ window", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150"&& HasMedLooseCSV>0.,0,0,"1"),
      TableRow("M$_{\\rm CT}>$ 170 GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170"&& HasMedLooseCSV>0.,0,0,"1"),
      TableRow("E$_{\\rm T}^{\\rm miss}>125$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125"&& HasMedLooseCSV>0.,0,0,"1"),
      TableRow("M$_{\\rm T}>150$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&mt_met_lep>150"&& HasMedLooseCSV>0.,0,0,"1"),
      TableRow("E$_{\\rm T}^{\\rm miss}>200$ GeV", WHLeptons==1&&"PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>200&&mt_met_lep>150"&& HasMedLooseCSV>0.,0,0,"1")

        }, sample_list,false);

  Table & cutflow_no_btag_sf = pm.Push<Table>("cutflow_no_btag_sf", vector<TableRow>{
      TableRow("1 good lep",WHLeptons>=1,0,0,"w_noBtagSF"),
      TableRow("2nd lep veto",WHLeptons==1&&"nvetoleps==1",0,0,"w_noBtagSF"),
      TableRow("Pass track veto", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto",0,0,"w_noBtagSF"),
      TableRow("Pass tau veto", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto",0,0,"w_noBtagSF"),
      TableRow("==2jets", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2",0,0,"w_noBtagSF"),
      TableRow("one loose, one med btag", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("2 med btag (deepCSV)", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("2 loose btag (deepCSV)", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("one loose, one med btag (deepCSV)", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
     // TableRow("==2btags", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("M$_{\\rm b\\bar{b}}$ window", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("M$_{\\rm CT}>$ 170 GeV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("E$_{\\rm T}^{\\rm miss}>125$ GeV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("M$_{\\rm T}>150$ GeV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&mt_met_lep>150"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("$125<$E$_{\\rm T}^{\\rm miss}<200$ GeV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&pfmet<200&&mt_met_lep>150"&& HasMedLooseCSV>0.,1,0,"w_noBtagSF"),
      TableRow("E$_{\\rm T}^{\\rm miss}>200$ GeV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&mbb>90&&mbb<150&&mct>170&&pfmet>200&&mt_met_lep>150"&& HasMedLooseCSV>0.,0,0,"w_noBtagSF")

        }, sample_list,false);


  Table & cutflow_deepBTagger_preselection = pm.Push<Table>("cutflow_deepBTagger_preselection", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mt_met_lep>50" && HasMedLooseCSV>0.,0,0,"w_noBtagSF"),
      TableRow("2 loose btag", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mt_met_lep>50" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("2 loose btag, no med", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mt_met_lep>50" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("one loose, one med btag", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mt_met_lep>50" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("2 med btag", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mt_met_lep>50" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF")
	
        }, short_sample_list,true);

  Table & cutflow_deepBTagger_signalRegion = pm.Push<Table>("cutflow_deepBTagger_signalRegion", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150&&ngoodbtags==2",0,0,"w_noBtagSF*w_BtagSF_medmed")
	
        }, new_sample_list,true);

    Table & cutflow_deepBTagger_signalRegion_debug = pm.Push<Table>("cutflow_deepBTagger_signalRegion_debug", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH")
	
        }, new_sample_list,true);

Table & cutflow_deepBTagger_signalRegion_debug2 = pm.Push<Table>("cutflow_deepBTagger_signalRegion_debug2", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>170&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH")
	
        }, new_sample_list,true);

  Table & cutflow_deepBTagger_signalRegionLowMet = pm.Push<Table>("cutflow_deepBTagger_signalRegionLowMet", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_medmed")
	
      }, new_sample_list,true);

  Table & cutflow_deepBTagger_signalRegionMedMet = pm.Push<Table>("cutflow_deepBTagger_signalRegionMedMet", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_medmed")
	
      }, new_sample_list,true);

  Table & cutflow_deepBTagger_signalRegionHighMet = pm.Push<Table>("cutflow_deepBTagger_signalRegionHighMet", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_medmed")
	
      }, new_sample_list,true);

  Table & cutflow_deepBTagger_signalRegionHighMet_debug = pm.Push<Table>("cutflow_deepBTagger_signalRegionHighMet_debug", vector<TableRow>{
      TableRow("one loose, one med btag, CSV", WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag",                 WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 loose btag, no med",         WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasLooseNoMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("one loose, one med btag",      WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("2 med btag",                   WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150" && HasMedMedDeepCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH")
	
      }, new_sample_list,true);

  Table & cutflow_METmctOpt_2jet = pm.Push<Table>("cutflow_METmctOpt_2jet", vector<TableRow>{
      TableRow("MET$>$125, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$125, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$125, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$125, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$200, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$200, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$200, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$200, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$225, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>225&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF")
	
      }, short_sample_list,true);

  Table & cutflow_METmctOpt_3jet = pm.Push<Table>("cutflow_METmctOpt_3jet", vector<TableRow>{
      TableRow("MET$>$125, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$125, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$125, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$125, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$200, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$200, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$200, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$200, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH"),
      TableRow("MET$>$225, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>225&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseCSV>0.,0,0,"w_noBtagSF*w_BtagSF_WH")
	
      }, short_sample_list,true);

  Table & cutflow_METmctOpt_2jet_HighMET = pm.Push<Table>("cutflow_METmctOpt_2jet_HighMET", vector<TableRow>{
      TableRow("300$>$MET$>$200, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$300, mct$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$300, mct$>$225",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>225&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$300, mct$>$250",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>250&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$300, mct$>$275",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&mct>275&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF")
	
      }, short_sample_list,true);

  /*Table & cutflow_METscan_newBabies = pm.Push<Table>("cutflow_METscan_newBabies", vector<TableRow>{
      TableRow("200$>$MET$>$125",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("400$>$MET$>$300",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>300&&pfmet<400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$400",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&pfmet>400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,1,"w_noBtagSF"),
      TableRow("200$>$MET$>$125",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasExactMedMedDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasExactMedMedDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("400$>$MET$>$300",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>300&&pfmet<400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasExactMedMedDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$400",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasExactMedMedDeepCSV>0.,0,1,"w_noBtagSF"),
      TableRow("200$>$MET$>$125",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>125&&pfmet<200&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("300$>$MET$>$200",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>200&&pfmet<300&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("400$>$MET$>$300",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>300&&pfmet<400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF"),
      TableRow("MET$>$400",WHLeptons==1&&"nvetoleps==1&&PassTrackVeto&&PassTauVeto&&ngoodjets==3&&pfmet>400&&mct>200&&mt_met_lep>150&&mbb>90&&mbb<150"&&HasMedLooseDeepCSV>0.,0,0,"w_noBtagSF")
	
	}, short_sample_list,true);*/


//   Table & cutflow = pm.Push<Table>("cutflow", vector<TableRow>{
//       TableRow("1 good lep, 2 jets, $E_{\\rm T}^{\\rm miss}>100$ GeV",WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100"),
      
//       TableRow("Pass track veto", WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100&&PassTrackVeto"),
//       TableRow("Pass tau veto", WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto"),
//       TableRow("==2jets", WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2"),
//       TableRow("==2btags", WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2"),
//       TableRow("M$_{\\rm b#bar{b}}$ window", WHLeptons==1&&"nvetoleps==0&&ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150"),
//       TableRow("M$_{\\rm CT}$", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170"),
//       TableRow("$E_{\\rm T}^{\\rm miss}>125$ GeV", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125"),
//       TableRow("M$_{\\rm T}>150$ GeV", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&mt_met_lep>150")


//         }, sample_list);

// Table & cutflow_no_weight = pm.Push<Table>("cutflow_raw", vector<TableRow>{
//       TableRow("1 good lep, 2 jets, $E_{\\rm T}^{\\rm miss}>100$ GeV",WHLeptons==1&&"ngoodjets>=2&&pfmet>100",0,0,"1"),
      
//       TableRow("Pass track veto", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto",0,0,"1"),
//       TableRow("Pass tau veto", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto",0,0,"1"),
//       TableRow("==2jets", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2",0,0,"1"),
//       TableRow("==2btags", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2",0,0,"1"),
//       TableRow("M$_{\\rm b#bar{b}}$ window", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150",0,0,"1"),
//       TableRow("M$_{\\rm CT}$", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170",0,0,"1"),
//       TableRow("$E_{\\rm T}^{\\rm miss}>125$ GeV", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125",0,0,"1"),
//       TableRow("M$_{\\rm T}>150$ GeV", WHLeptons==1&&"ngoodjets>=2&&pfmet>100&&PassTrackVeto&&PassTauVeto&&ngoodjets==2&&ngoodbtags==2&&mbb>90&&mbb<150&&mct>170&&pfmet>125&&mt_met_lep>150",0,0,"1")


//         }, sample_list);


  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);

  vector<GammaParams> yields = cutflow.BackgroundYield(lumi);
  for(const auto &yield: yields){
    cout << yield << endl;
  }
  vector<GammaParams> yields_no_weight = cutflow_no_weight.BackgroundYield(lumi);
  for(const auto &yield: yields_no_weight){
    cout << yield << endl;
  }
  vector<GammaParams> yields_no_btag_sf = cutflow_no_btag_sf.BackgroundYield(lumi);
  for(const auto &yield: yields_no_btag_sf){
    cout << yield << endl;
  }
  vector<GammaParams> yields_deepBTagger_preselection = cutflow_deepBTagger_preselection.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_preselection){
    cout << yield << endl;
  }
  vector<GammaParams> yields_deepBTagger_signalRegion = cutflow_deepBTagger_signalRegion.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegion){
    cout << yield << endl;
  }

  vector<GammaParams> yields_deepBTagger_signalRegion_debug = cutflow_deepBTagger_signalRegion_debug.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegion_debug){
    cout << yield << endl;
  }
  
  vector<GammaParams> yields_deepBTagger_signalRegion_debug2 = cutflow_deepBTagger_signalRegion_debug2.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegion_debug2){
    cout << yield << endl;
  }

  vector<GammaParams> yields_deepBTagger_signalRegionLowMet = cutflow_deepBTagger_signalRegionLowMet.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegionLowMet){
    cout << yield << endl;
  }

  vector<GammaParams> yields_deepBTagger_signalRegionMedMet = cutflow_deepBTagger_signalRegionMedMet.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegionMedMet){
    cout << yield << endl;
  }

  vector<GammaParams> yields_deepBTagger_signalRegionHighMet = cutflow_deepBTagger_signalRegionHighMet.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegionHighMet){
    cout << yield << endl;
  }

vector<GammaParams> yields_deepBTagger_signalRegionHighMet_debug = cutflow_deepBTagger_signalRegionHighMet_debug.BackgroundYield(lumi);
  for(const auto &yield: yields_deepBTagger_signalRegionHighMet_debug){
    cout << yield << endl;
  }

  vector<GammaParams> yields_METmctOpt_2jet = cutflow_METmctOpt_2jet.BackgroundYield(lumi);
  for(const auto &yield: yields_METmctOpt_2jet){
    cout << yield << endl;
  }

  vector<GammaParams> yields_METmctOpt_3jet = cutflow_METmctOpt_3jet.BackgroundYield(lumi);
  for(const auto &yield: yields_METmctOpt_3jet){
    cout << yield << endl;
  }

  vector<GammaParams> yields_METmctOpt_2jet_HighMET = cutflow_METmctOpt_2jet_HighMET.BackgroundYield(lumi);
  for(const auto &yield: yields_METmctOpt_2jet){
    cout << yield << endl;
  }

  /*vector<GammaParams> yields_METscan_newBabies = cutflow_METscan_newBabies.BackgroundYield(lumi);
  for(const auto &yield: yields_METmctOpt_2jet){
    cout << yield << endl;
  }*/
  
}
