#include <set>
#include <string>
#include <iostream>
#include <memory>
#include <limits>

#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TError.h"
#include "TGraphAsymmErrors.h"

#include "timer.hpp"
#include "baby_full.hpp"

using namespace std;

namespace{
  float bignum = numeric_limits<float>::max();
}

TH2D RowNorm(const TH2D &hin){
  TH2D hout = hin;
  hout.SetName((hout.GetName()+string("_rownorm")).c_str());
  for(int irow = 1; irow <= hout.GetNbinsY(); ++irow){
    double y = 0.;
    for(int icol = 1; icol <= hout.GetNbinsX(); ++icol){
      y += hout.GetBinContent(icol, irow);
    }
    y = 100./y;
    for(int icol = 1; icol <= hout.GetNbinsY(); ++icol){
      hout.SetBinContent(icol, irow, hout.GetBinContent(icol, irow)*y);
      hout.SetBinError(icol, irow, hout.GetBinError(icol, irow)*y);
    }
  }
  hout.SetMinimum(0.);
  hout.SetMaximum(100.);
  return hout;
}

TH2D ColNorm(const TH2D &hin){
  TH2D hout = hin;
  hout.SetName((hout.GetName()+string("_colnorm")).c_str());
  for(int icol = 1; icol <= hout.GetNbinsY(); ++icol){
    double y = 0.;
    for(int irow = 1; irow <= hout.GetNbinsX(); ++irow){
      y += hout.GetBinContent(icol, irow);
    }
    y = 100./y;
    for(int irow = 1; irow <= hout.GetNbinsY(); ++irow){
      hout.SetBinContent(icol, irow, hout.GetBinContent(icol, irow)*y);
      hout.SetBinError(icol, irow, hout.GetBinError(icol, irow)*y);
    }
  }
  hout.SetMinimum(0.);
  hout.SetMaximum(100.);
  return hout;
}

int RegionIndex(const Baby &b, size_t i,
		double low_nj, double high_nj,
		double low_met, double high_met,
		bool cluster_leps){
  double low_nb = 0.5;
  if(b.mm_nleps()->at(i)==2){
    low_nb = -0.5;
  }
  if(cluster_leps){
    low_nj += 1.-b.mm_nleps()->at(i);
    high_nj += 1.-b.mm_nleps()->at(i);
  }
  bool pass_baseline = 
    b.mm_ht()->at(i)>500
    && b.mm_met()->at(i)>low_met && b.mm_met()->at(i)<=high_met
    && b.mm_njets()->at(i)>low_nj && b.mm_njets()->at(i)<=high_nj
    && b.mm_nbm()->at(i)>low_nb;

  if(!pass_baseline) return -1;
  double mj = cluster_leps ? b.mm_mj14_lep()->at(i) : b.mm_mj14_nolep()->at(i);
  if(mj<=250.) return 0;
  double mt = b.mm_mt()->at(i);
  switch(b.mm_nleps()->at(i)){
  case 0:
    if(mj<=400.) return 1;
    else return 2;
  case 1:
    if(mj<=400.){
      if(mt<=140.) return 3;
      else return 5;
    }else{
      if(mt<=140.) return 4;
      else return 6;
    }
  case 2:
    if(mj<=400.) return 7;
    else return 8;
  default: return -1;
  }
}

void Fill(bool pass, TH1D &h_pass, TH1D &h_total, double x, double w){
  if(pass) h_pass.Fill(x, w);
  h_total.Fill(x, w);
}

void Format(TH1D &h, bool no_stats = false){
  h.SetLineColor(kBlack);
  h.SetLineStyle(1);
  h.SetLineWidth(5);
  if(no_stats) h.SetStats(false);
}

void Format(TH2D &h, bool no_stats = false){
  //h.SetMarkerStyle(20);
  //h.SetMarkerSize(0.5);
  h.SetMarkerColor(kBlack);
  h.SetMinimum(0.);
  if(no_stats) h.SetStats(false);
}

void Print(TH1D &h, bool no_stats = false){
  Format(h, no_stats);
  TCanvas c;
  h.Draw("e1p");
  c.Print((string("plots/")+h.GetName()+".pdf").c_str());
}

void MakePositive(TH1D &h){
  for(int i=0; i < h.GetNbinsX(); ++i){
    double y = h.GetBinContent(i);
    if(y >= 0.) continue;
    double e = h.GetBinError(i);
    e=hypot(e,y);
    y=0.;
    h.SetBinContent(i, y);
    h.SetBinError(i, e);
  }
}

void SetAxisLabels(TAxis &a){
  a.SetBinLabel(1, "Very Low MJ");
  a.SetBinLabel(2, "0l, Low MJ");
  a.SetBinLabel(3, "0l, High MJ");
  a.SetBinLabel(4, "R1");
  a.SetBinLabel(5, "R2");
  a.SetBinLabel(6, "R3");
  a.SetBinLabel(7, "R4");
  a.SetBinLabel(8, "D3");
  a.SetBinLabel(9, "D4");
}

void PrintTransfer(TH2D &h, bool region_labels = false){
  if(region_labels){
    SetAxisLabels(*h.GetXaxis());
    SetAxisLabels(*h.GetYaxis());
  }
  Format(h, true);
  TCanvas c;
  h.Draw("col");
  h.Draw("text same");
  c.Print((string("plots/")+h.GetName()+".pdf").c_str());
  TH2D rn = RowNorm(h);
  rn.Draw("col");
  rn.Draw("text same");
  c.Print((string("plots/")+h.GetName()+"_rownorm.pdf").c_str());
  TH2D cn = ColNorm(h);
  cn.Draw("col");
  cn.Draw("text same");
  c.Print((string("plots/")+h.GetName()+"_colnorm.pdf").c_str());
}

void Print(TH2D &h, bool no_stats = false){
  Format(h, no_stats);
  h.SetTitle(("#rho="+to_string(h.GetCorrelationFactor())).c_str());
  TCanvas c;
  h.Draw("col");
  h.Draw("scat same");
  c.Print((string("plots/")+h.GetName()+".pdf").c_str());
}

void Print(TH1D &pass, TH1D &total){
  MakePositive(pass);
  MakePositive(total);
  Print(pass);
  Print(total);
  TCanvas c;
  TGraphAsymmErrors g(&pass, &total);
  g.SetTitle((string(";")+pass.GetXaxis()->GetTitle()+";Fraction Mismeasured").c_str());
  g.SetMinimum(0.);
  g.SetMaximum(1.);
  g.Draw("ap0");
  c.Print((string("plots/")+g.GetName()+".pdf").c_str());
}

int main(){
  gErrorIgnoreLevel = 6000;
  gStyle->SetOptStat("iouMRen");
  double stops[2] = {0., 1.};
  double red[2] = {1., 1.};
  double green[2] = {1., 0.};
  double blue[2] = {1., 0.};
  gStyle->SetNumberContours(256);
  gStyle->SetPaintTextFormat("4.2f");
  TColor::CreateGradientColorTable(2, stops, red, green, blue, 256);

  double lumi = 2.6;
  int igood = 0;
  int ibad = 2;

  string folder_mc="/net/cms26/cms26r0/babymaker/babies/mismeasured_v2/2016_06_14/mc/merged_mm_std_nj5mj250/";
  auto baby_nontt = make_shared<Baby_full>(set<string>{
      folder_mc+"*_DYJetsToLL*.root",
	folder_mc+"*_QCD_HT*.root",
	folder_mc+"*_ST*channel*.root",
	folder_mc+"*_TTGJets*.root",
	folder_mc+"*_TTTT*.root",
	folder_mc+"*_TTWJetsTo*.root",
	folder_mc+"*_TTZTo*.root",
	folder_mc+"*_WH_HToBB*.root",
	folder_mc+"*_WJetsToLNu_HT-*.root",
	folder_mc+"*_WWTo*.root",
	folder_mc+"*_WZTo*.root",
	folder_mc+"*_ZH*.root",
	folder_mc+"*_ZZ*.root",
	folder_mc+"*_ttHJetTobb*.root"
	});
  auto baby_tt = make_shared<Baby_full>(set<string>{
      folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"
	});

  TH1D h_1l_mm_lep_pt("h_1l_mm_lep_pt", ";Mismeasured Lepton p_{T} [GeV];Entries", 20, 0., 2000.);
  TH2D h_1l_mt_mj_lep("h_1l_mt_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mt_mj_nolep("h_1l_mt_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mm_mj_lep("h_1l_mm_mt_mm_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mm_mj_nolep("h_1l_mm_mt_mm_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mj_nolep_mj_lep("h_1l_mj_nolep_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mj_nolep_mm_mj_lep("h_1l_mm_mj_nolep_mm_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mt("h_1l_mm_mt_mt", ";Correct m_{T} [GeV];Mismeasured m_{T} [GeV]", 100, -1000., 1000., 100, 0., 2000.);
  TH2D h_1l_mm_mj_nolep_mj_nolep("h_1l_mm_mj_nolep_mj_nolep", ";Correct M_{J} (no lep) [GeV];Mismeasured M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mj_lep_mj_lep("h_1l_mm_mj_lep_mj_lep", ";Correct M_{J} (with lep) [GeV];Mismeasured M_{J} (with lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH1D h_2l_mm_lep_pt("h_2l_mm_lep_pt", ";Mismeasured Lepton p_{T} [GeV];Entries", 20, 0., 2000.);
  TH2D h_2l_mt_mj_lep("h_2l_mt_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mt_mj_nolep("h_2l_mt_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mm_mj_lep("h_2l_mm_mt_mm_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mm_mj_nolep("h_2l_mm_mt_mm_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mj_nolep_mj_lep("h_2l_mj_nolep_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mj_nolep_mm_mj_lep("h_2l_mm_mj_nolep_mm_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mt("h_2l_mm_mt_mt", ";Correct m_{T} [GeV];Mismeasured m_{T} [GeV]", 100, -1000., 1000., 100, 0., 2000.);
  TH2D h_2l_mm_mj_nolep_mj_nolep("h_2l_mm_mj_nolep_mj_nolep", ";Correct M_{J} (no lep) [GeV];Mismeasured M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mj_lep_mj_lep("h_2l_mm_mj_lep_mj_lep", ";Correct M_{J} (with lep) [GeV];Mismeasured M_{J} (with lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);

  TH2D h_1l_dmt_dmjlep("h_1l_dmt_dmjlep", ";#Delta M_{J} (with lep) [GeV];#Delta m_{T} [GeV]", 60, -100., 500., 110, -100., 1000.);
  TH2D h_1l_dmt_dmjnolep("h_1l_dmt_dmjnolep", ";#Delta M_{J} (no lep) [GeV];#Delta m_{T} [GeV]", 60, -300., 300., 110, -100., 1000.);
  TH2D h_1l_dmt_dmjlep_highmet("h_1l_dmt_dmjlep_highmet", ";#Delta M_{J} (with lep) [GeV];#Delta m_{T} [GeV]", 60, -100., 500., 110, -100., 1000.);
  TH2D h_1l_dmt_dmjnolep_highmet("h_1l_dmt_dmjnolep_highmet", ";#Delta M_{J} (no lep) [GeV];#Delta m_{T} [GeV]", 60, -300., 300., 110, -100., 1000.);
  TH2D h_1l_dmt_dmjlep_highnj("h_1l_dmt_dmjlep_highnj", ";#Delta M_{J} (with lep) [GeV];#Delta m_{T} [GeV]", 60, -100., 500., 110, -100., 1000.);
  TH2D h_1l_dmt_dmjnolep_highnj("h_1l_dmt_dmjnolep_highnj", ";#Delta M_{J} (no lep) [GeV];#Delta m_{T} [GeV]", 60, -300., 300., 110, -100., 1000.);

  TH2D h_transfer_mm_lep("h_transfer_mm_lep", "Baseline;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_all_lep("h_transfer_all_lep", "Baseline;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_mm_nolep("h_transfer_mm_nolep", "Baseline;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_all_nolep("h_transfer_all_nolep", "Baseline;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lowmet_mm_lep("h_transfer_lowmet_mm_lep", "Low MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lowmet_all_lep("h_transfer_lowmet_all_lep", "Low MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lowmet_mm_nolep("h_transfer_lowmet_mm_nolep", "Low MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lowmet_all_nolep("h_transfer_lowmet_all_nolep", "Low MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_medmet_mm_lep("h_transfer_medmet_mm_lep", "Med MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_medmet_all_lep("h_transfer_medmet_all_lep", "Med MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_medmet_mm_nolep("h_transfer_medmet_mm_nolep", "Med MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_medmet_all_nolep("h_transfer_medmet_all_nolep", "Med MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highmet_mm_lep("h_transfer_highmet_mm_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highmet_all_lep("h_transfer_highmet_all_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highmet_mm_nolep("h_transfer_highmet_mm_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highmet_all_nolep("h_transfer_highmet_all_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lownj_mm_lep("h_transfer_lownj_mm_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lownj_all_lep("h_transfer_lownj_all_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lownj_mm_nolep("h_transfer_lownj_mm_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_lownj_all_nolep("h_transfer_lownj_all_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highnj_mm_lep("h_transfer_highnj_mm_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highnj_all_lep("h_transfer_highnj_all_lep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highnj_mm_nolep("h_transfer_highnj_mm_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);
  TH2D h_transfer_highnj_all_nolep("h_transfer_highnj_all_nolep", "High MET;Correct Region;Mismeasured Region", 9, -0.5, 8.5, 9, -0.5, 8.5);

  TH1D h_num_1l_mt("h_num_1l_mt", "Numerator;m_{T} [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_1l_mt("h_den_1l_mt", "Denominator;m_{T} [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_1l_mj_lep("h_num_1l_mj_lep", "Numerator;M_{J} (with lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_1l_mj_lep("h_den_1l_mj_lep", "Denominator;M_{J} (with lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_1l_mj_nolep("h_num_1l_mj_nolep", "Numerator;M_{J} (no lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_1l_mj_nolep("h_den_1l_mj_nolep", "Denominator;M_{J} (no lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_1l_met("h_num_1l_met", "Numerator;MET [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_1l_met("h_den_1l_met", "Denominator;MET [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_1l_njets("h_num_1l_njets", "Numerator;N_{jets};Entries", 7, 4.5, 11.5);
  TH1D h_den_1l_njets("h_den_1l_njets", "Denominator;N_{jets};Entries", 7, 4.5, 11.5);
  TH1D h_num_1l_nbm("h_num_1l_nbm", "Numerator;N_{b};Entries", 5, -0.5, 4.5);
  TH1D h_den_1l_nbm("h_den_1l_nbm", "Denominator;N_{b};Entries", 5, -0.5, 4.5);
  TH1D h_num_2l_mt("h_num_2l_mt", "Numerator;m_{T} [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_2l_mt("h_den_2l_mt", "Denominator;m_{T} [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_2l_mj_lep("h_num_2l_mj_lep", "Numerator;M_{J} (with lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_2l_mj_lep("h_den_2l_mj_lep", "Denominator;M_{J} (with lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_2l_mj_nolep("h_num_2l_mj_nolep", "Numerator;M_{J} (no lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_2l_mj_nolep("h_den_2l_mj_nolep", "Denominator;M_{J} (no lep) [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_2l_met("h_num_2l_met", "Numerator;MET [GeV];Entries", 20, 0., 1000.);
  TH1D h_den_2l_met("h_den_2l_met", "Denominator;MET [GeV];Entries", 20, 0., 1000.);
  TH1D h_num_2l_njets("h_num_2l_njets", "Numerator;N_{jets};Entries", 7, 4.5, 11.5);
  TH1D h_den_2l_njets("h_den_2l_njets", "Denominator;N_{jets};Entries", 7, 4.5, 11.5);
  TH1D h_num_2l_nbm("h_num_2l_nbm", "Numerator;N_{b};Entries", 5, -0.5, 4.5);
  TH1D h_den_2l_nbm("h_den_2l_nbm", "Denominator;N_{b};Entries", 5, -0.5, 4.5);

  for(auto &b: {baby_nontt, baby_tt}){
    cout << "Getting entries..." << endl;
    long num_entries = b->GetEntries();
    cout << num_entries << " entries found." << endl;
    Timer timer(num_entries, 1.);
    for(long event = 0; event < num_entries; ++event){
      timer.Iterate();
      b->GetEntry(event);

      if(!b->stitch()) continue;
      double w = lumi*b->weight();

      bool mm = b->mm()->at(ibad);

      if(mm){
	h_transfer_mm_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., bignum, true),
			       RegionIndex(*b, ibad, 5.5, bignum, 200., bignum, true),
			       w);
	h_transfer_mm_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., bignum, false),
				 RegionIndex(*b, ibad, 5.5, bignum, 200., bignum, false),
				 w);
	h_transfer_lowmet_mm_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., 350., true),
				      RegionIndex(*b, ibad, 5.5, bignum, 200., 350., true),
				      w);
	h_transfer_lowmet_mm_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., 350., false),
					RegionIndex(*b, ibad, 5.5, bignum, 200., 350., false),
					w);
	h_transfer_medmet_mm_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 350., 500., true),
				      RegionIndex(*b, ibad, 5.5, bignum, 350., 500., true),
				      w);
	h_transfer_medmet_mm_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 350., 500., false),
					RegionIndex(*b, ibad, 5.5, bignum, 350., 500., false),
					w);
	h_transfer_highmet_mm_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 500., bignum, true),
				       RegionIndex(*b, ibad, 5.5, bignum, 500., bignum, true),
				       w);
	h_transfer_highmet_mm_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 500., bignum, false),
					 RegionIndex(*b, ibad, 5.5, bignum, 500., bignum, false),
					 w);
	h_transfer_lownj_mm_lep.Fill(RegionIndex(*b, igood, 5.5, 8.5, 200., bignum, true),
				     RegionIndex(*b, ibad, 5.5, 8.5, 200., bignum, true),
				     w);
	h_transfer_lownj_mm_nolep.Fill(RegionIndex(*b, igood, 5.5, 8.5, 200., bignum, false),
				       RegionIndex(*b, ibad, 5.5, 8.5, 200., bignum, false),
				       w);
	h_transfer_highnj_mm_lep.Fill(RegionIndex(*b, igood, 8.5, bignum, 200., bignum, true),
				      RegionIndex(*b, ibad, 8.5, bignum, 200., bignum, true),
				      w);
	h_transfer_highnj_mm_nolep.Fill(RegionIndex(*b, igood, 8.5, bignum, 200., bignum, false),
					RegionIndex(*b, ibad, 8.5, bignum, 200., bignum, false),
					w);
      }
      h_transfer_all_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., bignum, true),
			      RegionIndex(*b, ibad, 5.5, bignum, 200., bignum, true),
			      w);
      h_transfer_all_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., bignum, false),
				RegionIndex(*b, ibad, 5.5, bignum, 200., bignum, false),
				w);
      h_transfer_lowmet_all_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., 350., true),
				     RegionIndex(*b, ibad, 5.5, bignum, 200., 350., true),
				     w);
      h_transfer_lowmet_all_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 200., 350., false),
				       RegionIndex(*b, ibad, 5.5, bignum, 200., 350., false),
				       w);
      h_transfer_medmet_all_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 350., 500., true),
				     RegionIndex(*b, ibad, 5.5, bignum, 350., 500., true),
				     w);
      h_transfer_medmet_all_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 350., 500., false),
				       RegionIndex(*b, ibad, 5.5, bignum, 350., 500., false),
				       w);
      h_transfer_highmet_all_lep.Fill(RegionIndex(*b, igood, 5.5, bignum, 500., bignum, true),
				      RegionIndex(*b, ibad, 5.5, bignum, 500., bignum, true),
				      w);
      h_transfer_highmet_all_nolep.Fill(RegionIndex(*b, igood, 5.5, bignum, 500., bignum, false),
					RegionIndex(*b, ibad, 5.5, bignum, 500., bignum, false),
					w);
      h_transfer_lownj_all_lep.Fill(RegionIndex(*b, igood, 5.5, 8.5, 200., bignum, true),
				    RegionIndex(*b, ibad, 5.5, 8.5, 200., bignum, true),
				    w);
      h_transfer_lownj_all_nolep.Fill(RegionIndex(*b, igood, 5.5, 8.5, 200., bignum, false),
				      RegionIndex(*b, ibad, 5.5, 8.5, 200., bignum, false),
				      w);
      h_transfer_highnj_all_lep.Fill(RegionIndex(*b, igood, 8.5, bignum, 200., bignum, true),
				     RegionIndex(*b, ibad, 8.5, bignum, 200., bignum, true),
				     w);
      h_transfer_highnj_all_nolep.Fill(RegionIndex(*b, igood, 8.5, bignum, 200., bignum, false),
				       RegionIndex(*b, ibad, 8.5, bignum, 200., bignum, false),
				       w);

      if(mm
	 && (b->mm_nleps()->at(igood)>0 || b->mm_nleps()->at(ibad)>0)
	 && (b->mm_ht()->at(igood)>500 || b->mm_ht()->at(ibad)>500)
	 && (b->mm_met()->at(igood)>200 || b->mm_met()->at(ibad)>200)
	 && (b->mm_njets()->at(igood)>=5 || b->mm_njets()->at(ibad)>=5)){

	if(b->mm_nleps()->at(ibad)==1){
	  h_1l_dmt_dmjlep.Fill(b->mm_mj14_lep()->at(ibad)-b->mm_mj14_lep()->at(igood),
			       b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
			       10.*w);
	  h_1l_dmt_dmjnolep.Fill(b->mm_mj14_nolep()->at(ibad)-b->mm_mj14_nolep()->at(igood),
				 b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
				 10.*w);
	  if(b->mm_met()->at(ibad)>500.){
	    h_1l_dmt_dmjlep_highmet.Fill(b->mm_mj14_lep()->at(ibad)-b->mm_mj14_lep()->at(igood),
					 b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
					 10.*w);
	    h_1l_dmt_dmjnolep_highmet.Fill(b->mm_mj14_nolep()->at(ibad)-b->mm_mj14_nolep()->at(igood),
					   b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
					   10.*w);
	  }
	  if(b->mm_njets()->at(ibad)>8.5){
	    h_1l_dmt_dmjlep_highnj.Fill(b->mm_mj14_lep()->at(ibad)-b->mm_mj14_lep()->at(igood),
				 b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
				 10.*w);
	    h_1l_dmt_dmjnolep_highnj.Fill(b->mm_mj14_nolep()->at(ibad)-b->mm_mj14_nolep()->at(igood),
				   b->mm_mt()->at(ibad)-b->mm_mt()->at(igood),
				   10.*w);
	  }
	  h_1l_mm_lep_pt.Fill(b->mm_lep_pt()->at(ibad), w);
	  h_1l_mt_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mt()->at(igood), w);
	  h_1l_mt_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mt()->at(igood), w);
	  h_1l_mm_mt_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_1l_mm_mt_mm_mj_nolep.Fill(b->mm_mj14_nolep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_1l_mj_nolep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_nolep()->at(igood), w);
	  h_1l_mm_mj_nolep_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mj14_nolep()->at(ibad), w);
	  h_1l_mm_mt_mt.Fill(b->mm_mt()->at(igood), b->mm_mt()->at(ibad), w);
	  h_1l_mm_mj_nolep_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mj14_nolep()->at(ibad), w);
	  h_1l_mm_mj_lep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_lep()->at(ibad), w);
	}else if(b->mm_nleps()->at(ibad)==2){
	  h_2l_mm_lep_pt.Fill(b->mm_lep_pt()->at(ibad), w);
	  h_2l_mt_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mt()->at(igood), w);
	  h_2l_mt_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mt()->at(igood), w);
	  h_2l_mm_mt_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_2l_mm_mt_mm_mj_nolep.Fill(b->mm_mj14_nolep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_2l_mj_nolep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_nolep()->at(igood), w);
	  h_2l_mm_mj_nolep_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mj14_nolep()->at(ibad), w);
	  h_2l_mm_mt_mt.Fill(b->mm_mt()->at(igood), b->mm_mt()->at(ibad), w);
	  h_2l_mm_mj_nolep_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mj14_nolep()->at(ibad), w);
	  h_2l_mm_mj_lep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_lep()->at(ibad), w);
	}
      }

      if(b->mm_ht()->at(ibad)>500. && b->mm_met()->at(ibad)>200.
	 && b->mm_njets()->at(ibad)>=5 && b->mm_nbm()->at(ibad)>=1){
	if(b->mm_nleps()->at(ibad)==1){
	  Fill(mm, h_num_1l_mt, h_den_1l_mt, b->mm_mt()->at(ibad), w);
	  Fill(mm, h_num_1l_mj_lep, h_den_1l_mj_lep, b->mm_mj14_lep()->at(ibad), w);
	  Fill(mm, h_num_1l_mj_nolep, h_den_1l_mj_nolep, b->mm_mj14_nolep()->at(ibad), w);
	  Fill(mm, h_num_1l_met, h_den_1l_met, b->mm_met()->at(ibad), w);
	  Fill(mm, h_num_1l_njets, h_den_1l_njets, b->mm_njets()->at(ibad), w);
	  Fill(mm, h_num_1l_nbm, h_den_1l_nbm, b->mm_nbm()->at(ibad), w);
	}else if(b->mm_nleps()->at(ibad)>=2){
	  Fill(mm, h_num_2l_mt, h_den_2l_mt, b->mm_mt()->at(ibad), w);
	  Fill(mm, h_num_2l_mj_lep, h_den_2l_mj_lep, b->mm_mj14_lep()->at(ibad), w);
	  Fill(mm, h_num_2l_mj_nolep, h_den_2l_mj_nolep, b->mm_mj14_nolep()->at(ibad), w);
	  Fill(mm, h_num_2l_met, h_den_2l_met, b->mm_met()->at(ibad), w);
	  Fill(mm, h_num_2l_njets, h_den_2l_njets, b->mm_njets()->at(ibad), w);
	  Fill(mm, h_num_2l_nbm, h_den_2l_nbm, b->mm_nbm()->at(ibad), w);
	}
      }
    }
  }

  if(false){
    Print(h_num_1l_mt, h_den_1l_mt);
    Print(h_num_1l_mj_lep, h_den_1l_mj_lep);
    Print(h_num_1l_mj_nolep, h_den_1l_mj_nolep);
    Print(h_num_1l_met, h_den_1l_met);
    Print(h_num_1l_njets, h_den_1l_njets);
    Print(h_num_1l_nbm, h_den_1l_nbm);
    Print(h_num_2l_mt, h_den_2l_mt);
    Print(h_num_2l_mj_lep, h_den_2l_mj_lep);
    Print(h_num_2l_mj_nolep, h_den_2l_mj_nolep);
    Print(h_num_2l_met, h_den_2l_met);
    Print(h_num_2l_njets, h_den_2l_njets);
    Print(h_num_2l_nbm, h_den_2l_nbm);

    Print(h_1l_mm_lep_pt);
    Print(h_1l_mt_mj_lep);
    Print(h_1l_mt_mj_nolep);
    Print(h_1l_mm_mt_mm_mj_lep);
    Print(h_1l_mm_mt_mm_mj_nolep);
    Print(h_1l_mj_nolep_mj_lep);
    Print(h_1l_mm_mj_nolep_mm_mj_lep);
    Print(h_1l_mm_mt_mt);
    Print(h_1l_mm_mj_nolep_mj_nolep);
    Print(h_1l_mm_mj_lep_mj_lep);
    Print(h_2l_mm_lep_pt);
    Print(h_2l_mt_mj_lep);
    Print(h_2l_mt_mj_nolep);
    Print(h_2l_mm_mt_mm_mj_lep);
    Print(h_2l_mm_mt_mm_mj_nolep);
    Print(h_2l_mj_nolep_mj_lep);
    Print(h_2l_mm_mj_nolep_mm_mj_lep);
    Print(h_2l_mm_mt_mt);
    Print(h_2l_mm_mj_nolep_mj_nolep);
    Print(h_2l_mm_mj_lep_mj_lep);
  }

  Print(h_1l_dmt_dmjlep);
  Print(h_1l_dmt_dmjnolep);
  Print(h_1l_dmt_dmjlep_highmet);
  Print(h_1l_dmt_dmjnolep_highmet);
  Print(h_1l_dmt_dmjlep_highnj);
  Print(h_1l_dmt_dmjnolep_highnj);

  PrintTransfer(h_transfer_mm_lep, true);
  PrintTransfer(h_transfer_all_lep, true);
  PrintTransfer(h_transfer_mm_nolep, true);
  PrintTransfer(h_transfer_all_nolep, true);

  PrintTransfer(h_transfer_mm_lep, true);
  PrintTransfer(h_transfer_all_lep, true);
  PrintTransfer(h_transfer_mm_nolep, true);
  PrintTransfer(h_transfer_all_nolep, true);
  PrintTransfer(h_transfer_lowmet_mm_lep, true);
  PrintTransfer(h_transfer_lowmet_all_lep, true);
  PrintTransfer(h_transfer_lowmet_mm_nolep, true);
  PrintTransfer(h_transfer_lowmet_all_nolep, true);
  PrintTransfer(h_transfer_medmet_mm_lep, true);
  PrintTransfer(h_transfer_medmet_all_lep, true);
  PrintTransfer(h_transfer_medmet_mm_nolep, true);
  PrintTransfer(h_transfer_medmet_all_nolep, true);
  PrintTransfer(h_transfer_highmet_mm_lep, true);
  PrintTransfer(h_transfer_highmet_all_lep, true);
  PrintTransfer(h_transfer_highmet_mm_nolep, true);
  PrintTransfer(h_transfer_highmet_all_nolep, true);
  PrintTransfer(h_transfer_lownj_mm_lep, true);
  PrintTransfer(h_transfer_lownj_all_lep, true);
  PrintTransfer(h_transfer_lownj_mm_nolep, true);
  PrintTransfer(h_transfer_lownj_all_nolep, true);
  PrintTransfer(h_transfer_highnj_mm_lep, true);
  PrintTransfer(h_transfer_highnj_all_lep, true);
  PrintTransfer(h_transfer_highnj_mm_nolep, true);
  PrintTransfer(h_transfer_highnj_all_nolep, true);
}