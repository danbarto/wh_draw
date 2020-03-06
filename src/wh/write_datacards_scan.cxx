
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "core/wh_functions.hpp"

using namespace std;
using namespace WH_Functions;

namespace{
  // bool systematic_mode=true;
  // bool debug = true;
  float lumi=137.;
  bool discovery_mode=false;
  bool boosted =true;
  bool original_analysis=false;
  enum Bkgs {bkg, top, w, other};

}

void writeCard(vector<string> bin_names, vector<vector<GammaParams> > allyields,vector<vector<float> > mct_transfer_factors,vector<GammaParams>  sigyields, TString mass_tag,TString analysis_tag);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////


  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline         = "";
  NamedFunc baselinef     = "pass&&nvetoleps>=1&&ngoodjets>=2 && ngoodjets<=3 &&pfmet>125&&mt_met_lep>150 && ngoodbtags==2 && nWHLeptons>=1";
  if (original_analysis) baselinef     = "pass&&nvetoleps>=1&&ngoodjets>=2 && ngoodjets<=2 &&pfmet>125&&mt_met_lep>150 && ngoodbtags>=1 && nloosebtags==2 && nWHLeptons>=1";
  
  TString single_lep      = "nvetoleps==1&&PassTrackVeto&&PassTauVeto&&nWHLeptons==1";
  TString dilep           = "nvetoleps==2";

  string mc2016_dir       = "/home/users/rheller/wh_babies/babies_v31_2_2019_10_03/s16v3/";
  string mc2016_dir_ttmet = "/home/users/rheller/wh_babies/babies_v30_9_2019_10_05/s16v3/";

  string mc2017_dir       = "/home/users/rheller/wh_babies/babies_v31_2_2019_10_03/f17v2/";
  string mc2017_dir_ttmet = "/home/users/rheller/wh_babies/babies_v30_9_2019_10_05/f17v2/";

  string mc2018_dir       = "/home/users/rheller/wh_babies/babies_v31_2_2019_10_03/a18v1/";
  string mc2018_dir_ttmet = "/home/users/rheller/wh_babies/babies_v30_9_2019_10_07/a18v1/";

  string signal_dir2016   = "/home/users/rheller/wh_babies/babies_signal_s16v3_v32_2019_10_07/";
  string signal_dir2017   = "/home/users/rheller/wh_babies/babies_signal_f17v2_v32_0_2019_10_07/";
  string signal_dir2018   = "/home/users/rheller/wh_babies/babies_signal_a18v1_v32_0_2019_10_07/";



  auto all_top = {mc2016_dir+"*TTJets_1lep_top_s16v3*.root",mc2016_dir+"*TTJets_1lep_tbar_s16v3*",mc2016_dir_ttmet+"*TTJets_1lep_*met150*.root",mc2017_dir+"*TTJets_1lep_top_f17v2*.root",mc2017_dir+"*TTJets_1lep_tbar_f17v2*",mc2017_dir_ttmet+"*TTJets_1lep_*met150*.root",mc2018_dir+"*TTJets_1lep_top_a18v1*.root",mc2018_dir+"*TTJets_1lep_tbar_a18v1*",mc2018_dir_ttmet+"*TTJets_1lep_*met80*.root",mc2016_dir+"*_ST_*.root",mc2017_dir+"*_ST_*.root",mc2018_dir+"*_ST_*.root",mc2016_dir+"*TTJets_2lep_s16v3*.root", mc2016_dir_ttmet+"*TTJets_2lep_*met150*.root",mc2017_dir+"*TTJets_2lep_f17v2*.root", mc2017_dir_ttmet+"*TTJets_2lep_*met150*.root",mc2018_dir+"*TTJets_2lep_a18v1*.root",mc2018_dir_ttmet+"*TTJets_2lep_*met80*.root"};


  auto all_other = {mc2016_dir+"*WW*.root", mc2016_dir+"*WZ*.root",mc2016_dir+"*ZZ*.root",mc2017_dir+"*WW*.root", mc2017_dir+"*WZ*.root",mc2017_dir+"*ZZ*.root",mc2018_dir+"*WW*.root", mc2018_dir+"*WZ*.root",mc2018_dir+"*ZZ*.root",mc2016_dir+"*_TTWJets*.root", mc2016_dir+"*_TTZ*.root",mc2017_dir+"*_TTWJets*.root", mc2017_dir+"*_TTZ*.root",mc2018_dir+"*_TTWJets*.root", mc2018_dir+"*_TTZ*.root"};

  //// Contributions
  auto proc_wjets = Process::MakeShared<Baby_full>("W+jets 2016-2018", Process::Type::background, kCyan-3, {mc2016_dir+"*slim_W*JetsToLNu_s16v3*",mc2016_dir+"*W*Jets_NuPt200_s16v*.root",mc2017_dir+"*slim_W*JetsToLNu_f17v2*",mc2017_dir+"*W*Jets_NuPt200_f17v2*.root",mc2018_dir+"*slim_W*JetsToLNu_a18v1*",mc2018_dir+"*W*Jets_NuPt200_a18v1*.root"},"stitch&&evt!=74125994"&&baselinef);
  auto proc_top = Process::MakeShared<Baby_full>("top 2016-2018", Process::Type::background, kRed,all_top,"stitch"&&baselinef);

  auto proc_other = Process::MakeShared<Baby_full>("TTV and VV 2016-2018", Process::Type::background, kRed,all_other,baselinef);

  auto proc_sig = Process::MakeShared<Baby_full>("2016-2018 TChiWH(750,1)", Process::Type::signal, colors("t1tttt"),
  {signal_dir2016+"*TChiWH*.root",signal_dir2017+"*TChiWH*.root",signal_dir2018+"*TChiWH*.root"},"mass_stop==750&&mass_lsp==1"&&baselinef);

  auto proc_sig_compress = Process::MakeShared<Baby_full>("2016-2018 TChiWH(225,75)", Process::Type::signal, colors("t1tttt"),
  {signal_dir2016+"*TChiWH*.root",signal_dir2017+"*TChiWH*.root",signal_dir2018+"*TChiWH*.root"},"mass_stop==225&&mass_lsp==75"&&baselinef);

    auto proc_sig_medium = Process::MakeShared<Baby_full>("2016-2018 TChiWH(350,150)", Process::Type::signal, colors("t1tttt"),
  {signal_dir2016+"*TChiWH*.root",signal_dir2017+"*TChiWH*.root",signal_dir2018+"*TChiWH*.root"},"mass_stop==350&&mass_lsp==150"&&baselinef);

  auto proc_sig_all = Process::MakeShared<Baby_full>("2016-2018 TChiWH", Process::Type::signal, colors("t1tttt"),
  {signal_dir2016+"*TChiWH*.root",signal_dir2017+"*TChiWH*.root",signal_dir2018+"*TChiWH*.root"},baselinef);


  vector<shared_ptr<Process> > all_procs_bkg = {proc_top,proc_wjets,proc_other};
  vector<shared_ptr<Process> > sig_procs = {proc_sig_all};//{proc_sig,proc_sig_medium,proc_sig_compress};


  TChain sig_tree("t");
  sig_tree.Add("/home/users/rheller/wh_babies/babies_signal_s16v3_v32_2019_10_07/slim_SMS_TChiWH_s16v3_0.root");

  long nentries(sig_tree.GetEntries());
  cout<<"Got "<<nentries<<" entries."<<endl;
  // TString outfolder = outpath;
  // outfolder.Remove(outfolder.Last('/')+1, outfolder.Length());
  // if(outfolder!="") gSystem->mkdir(outfolder, kTRUE);

  //Project sig_tree into mass plane to find points for skim
  TH2F * mass_plane = new TH2F("mchi_vs_mlsp","mchi_vs_mlsp",1000,-0.5,999.5,1000,-0.5,999.5);
  sig_tree.Project("mchi_vs_mlsp","mass_stop:mass_lsp","","colz");

  int ini_x = mass_plane->FindFirstBinAbove(0,1);
  int last_x = mass_plane->FindLastBinAbove(0,1);
  int ini_y = mass_plane->FindFirstBinAbove(0,2);
  int last_y = mass_plane->FindLastBinAbove(0,2);

  //load all pairs as cuts into vector
  vector<TString> pair_cuts;
  vector<TString> mass_tag;
  int Npoints=0;
  for(int iy=ini_y; iy<=last_y; iy++){
    for(int ix=ini_x; ix<=last_x; ix++){
      if(mass_plane->GetBinContent(ix,iy) > 0){
        int mchi = static_cast<int>(mass_plane->GetYaxis()->GetBinCenter(iy));
        int mlsp = static_cast<int>(mass_plane->GetXaxis()->GetBinCenter(ix));
        if(mchi!=750) continue;
        pair_cuts.push_back(Form("mass_stop==%i&&mass_lsp==%i",mchi,mlsp));
        mass_tag.push_back(Form("mChi-%i_mLSP-%i_",mchi,mlsp));
        cout<<"Found mass point "<<mass_tag.back()<<endl;
  Npoints++;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Defining cuts ///////////////////////////////////////////////
  // baseline defined above

	  vector<NamedFunc> weights; //size_t nsels; 
	  vector<TString> numerators,denominators;
	  vector<TString> leglabels;


  TString analysis_tag = "nominal";
  if(boosted) analysis_tag = "boosted";
  if(original_analysis) analysis_tag="original";



  // vector<NamedFunc> metbins = {"pfmet>125&&pfmet<=200","pfmet>200&&pfmet<=300","pfmet>300"}; 
  vector<NamedFunc> metbins = {"pfmet>125&&pfmet<=200","pfmet>200&&pfmet<=300","pfmet>300&&pfmet<400","pfmet>400"}; analysis_tag+="_4metbins";
  vector<NamedFunc> boosted_metbins = {"pfmet>125&&pfmet<=200","pfmet>200&&pfmet<=300","pfmet>300&&pfmet<400","pfmet>400"}; //analysis_tag+="2bin_boostedtest";
  // vector<NamedFunc> metbins = {"pfmet>125&&pfmet<=200","pfmet>200&&pfmet<=300","pfmet>300&&pfmet<450","pfmet>450"}; analysis_tag+="_4metbins_450";
 // if(original_analysis) vector<NamedFunc> metbins = {"pfmet>125&&pfmet<=200","pfmet>200"}; 
	//vector<NamedFunc> metbins = {"pfmet>125&&pfmet<=200","pfmet>200&&pfmet<=300","pfmet>350"};	analysis_tag+="_met350";
  // if(metbins.size()>3) analysis_tag+="_4metbins_450";
	vector<NamedFunc> njetbins = {"ngoodjets==2","ngoodjets==3"&&LeadingNonBJetPt_med<100.};
  if (original_analysis) njetbins = {"ngoodjets==2"};
	vector<NamedFunc> deepAK8bins = {max_ak8pfjets_deepdisc_hbb<=0.8,max_ak8pfjets_deepdisc_hbb>0.8};

	vector<string> metnames={"lowmet","medmet","highmet","vhighmet"};
	vector<string> njetnames={"nj2","nj3"};
	vector<string> htagnames={"res","boos"};

    weights= { "weight * w_pu" * yearWeight};

    leglabels = {"Nominal"};


  gSystem->mkdir("statistics/"+analysis_tag);
  gSystem->mkdir("statistics/"+analysis_tag+"/datacards/");
    // nsels = 2*weights.size();

    TString signal_region = single_lep+"&&mct>200&&mbb>90&&mbb<150";
    TString mct_control_region = single_lep+"&&mct>150&&mct<=200&&mbb>90&&mbb<150";    
    if(original_analysis) { signal_region = single_lep+"&&mct>170&&mbb>90&&mbb<150";
                             mct_control_region = single_lep+"&&mct>120&&mct<=170&&mbb>90&&mbb<150";}

    TString signal_region_boost = single_lep+"&&mct>200&&mbb>90&&mbb<150";
    TString mct_control_region_boost = single_lep+"&&mct<=200&&mbb>90&&mbb<150";

    // for(uint ivar=0;ivar<weights.size();ivar++){
    //   numerators.push_back(numerator);
    //   denominators.push_back(denominator);
    // }
 

  vector<int> proc_options = {0,1};//allbkg, top 
  vector<TString> tags = {"allbkg","top"};

  vector<vector<vector<int> > > indices;

  for(int i=0;i<static_cast<int>(denominators.size());i++){
    //numerator is always top only; W predicted separately.
    indices.push_back(vector<vector<int> >({{proc_options[bkg],2*i,-1},{proc_options[top],2*i+1,1}}));
  }
  cout<<"Boosted met bin size "<<boosted_metbins.size()<<endl;
  PlotMaker pm;
  vector<NamedFunc> allcuts;
  vector<NamedFunc> allcuts_sig;
  vector<TableRow> table_rows;
  vector<TableRow> table_rows_sig;
  vector<string> bin_names;
  vector<string> bin_names_sig;
  for(uint inj=0;inj<njetbins.size();inj++){
		for(uint imet=0;imet<metbins.size();imet++){
			for(uint ideepAK8=0;ideepAK8<deepAK8bins.size();ideepAK8++){
			if(!boosted && ideepAK8>0) continue;
      NamedFunc totcut = "1";
			NamedFunc totcut_sig = "1";
			string bin_name="";

			if(boosted){
				if(ideepAK8==0) totcut = njetbins[inj] && metbins[imet] && deepAK8bins[ideepAK8] && signal_region;
				else if (ideepAK8>0 && imet<boosted_metbins.size()) totcut = njetbins[inj] && boosted_metbins[imet] && deepAK8bins[ideepAK8] && signal_region_boost;
        else if (ideepAK8>0 && imet>=boosted_metbins.size()) continue; // stop adding boosted case after all boosted_metbins have been added.
        bin_name = njetnames[inj]+"_"+metnames[imet]+"_"+htagnames[ideepAK8];
				// bin_names.push_back(njetnames[inj]+"_"+metnames[imet]+"_"+htagnames[ideepAK8]);
			}
			else{ totcut = njetbins[inj] && metbins[imet] && signal_region;
        bin_name = njetnames[inj]+"_"+metnames[imet];
				// bin_names.push_back(njetnames[inj]+"_"+metnames[imet]);
			}

			allcuts.push_back(totcut);
			table_rows.push_back(TableRow("", totcut,0,0,weights[0]));
      bin_names.push_back(bin_name);

      for(int isig=0;isig<Npoints;isig++){
        totcut_sig = pair_cuts[isig] && totcut; 
        allcuts_sig.push_back(totcut_sig);
        bin_names_sig.push_back(mass_tag[isig].Data()+bin_name);
        table_rows_sig.push_back(TableRow("", totcut_sig,0,0,weights[0]));
      }
			
			if(boosted){
				if(ideepAK8==0) totcut= njetbins[inj] && metbins[imet] && deepAK8bins[ideepAK8] && mct_control_region;
				else if (ideepAK8>0 && imet<boosted_metbins.size()) totcut = njetbins[inj] && boosted_metbins[imet] && deepAK8bins[ideepAK8] && mct_control_region_boost;
				else if (ideepAK8>0 && imet>=boosted_metbins.size()) continue;
        bin_names.push_back("CR_lowmct_"+njetnames[inj]+"_"+metnames[imet]+"_"+htagnames[ideepAK8]);
			}
			else{
				totcut= njetbins[inj] && metbins[imet] && mct_control_region;
				bin_names.push_back("CR_lowmct_"+njetnames[inj]+"_"+metnames[imet]);
			}
			allcuts.push_back(totcut);
			table_rows.push_back(TableRow("", totcut,0,0,weights[0]));
			

		}
  }
}
  TString tname = "all_bins";
  pm.Push<Table>(tname.Data(),  table_rows_sig, sig_procs, false, false);
  pm.Push<Table>(tname.Data(),  table_rows, all_procs_bkg, false, false);
 
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////
  cout<<"Analysis tag is "<<analysis_tag<<endl;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  size_t nrows= table_rows.size();
	size_t nrows_sig= table_rows_sig.size();
  int nsig = mass_tag.size();
  int nbins = nrows_sig/nsig;
    // allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
   
  vector<GammaParams>  sigyields;
  vector<vector<GammaParams> > sig_by_mass(nsig,vector<GammaParams> (2*nbins));
  Table * yield_table_sig = static_cast<Table*>(pm.Figures()[0].get());
  sigyields = yield_table_sig->Yield(proc_sig_all.get(), lumi);

  for(size_t irow=0; irow<nrows_sig; irow++){ 
    cout<<"signal, "<<bin_names_sig[irow]<<": "<<setw(7)<<RoundNumber(sigyields[irow].Yield(), 2)<<endl;
    int bin_num = 2*(irow/nsig); //hack to skip CR bins
    int sig_num = irow - bin_num*nsig/2;
    cout<<"bin_num, sig_num: "<<bin_num<<" "<<sig_num<<endl;
    sig_by_mass[sig_num][bin_num] = sigyields[irow];
  }

  for(int isig =0;isig<nsig;isig++){
    cout<<mass_tag[isig]<<endl;
    for(size_t ibin =0;ibin<nrows/2;ibin++){
    cout<<bin_names[2*ibin]<<": "<<sig_by_mass[isig][2*ibin]<<endl;
    }
  }




  cout<<"Finished dumping yields"<<endl;
   //Extract yields
    vector<vector<GammaParams> > allyields(4,vector<GammaParams> (nrows));
    Table * yield_table = static_cast<Table*>(pm.Figures()[1].get());
	allyields[bkg] = yield_table->BackgroundYield(lumi);
	allyields[top] = yield_table->Yield(proc_top.get(), lumi);
	allyields[w] = yield_table->Yield(proc_wjets.get(), lumi);
	allyields[other] = yield_table->Yield(proc_other.get(), lumi);


	// sigyields[1] = yield_table_sig->Yield(proc_sig_medium.get(), lumi);
	// sigyields[2] = yield_table_sig->Yield(proc_sig_compress.get(), lumi);


  // vector<vector<GammaParams> > sigyields(sig_procs.size(),vector<GammaParams> (nrows_sig));
  // Table * yield_table_sig = static_cast<Table*>(pm.Figures()[1].get());
  // sigyields[0] = yield_table_sig->Yield(proc_sig_all.get(), lumi);
  // // sigyields[1] = yield_table_sig->Yield(proc_sig_medium.get(), lumi);
  // // sigyields[2] = yield_table_sig->Yield(proc_sig_compress.get(), lumi);

	vector<vector<float> > mct_transfer_factors;
	float val(1.), valup(1.), valdown(1.);
    //calculate mct transfer factors
    for(size_t ibin=0; ibin<nrows; ibin+=2){ //NB iterating by 2
		vector<vector<float> > entries;
		vector<vector<float> > sumweights;
		vector<float> powers;

		//Numerator
		entries.push_back(vector<float>());
		sumweights.push_back(vector<float>());
		entries.back().push_back(allyields[top][ibin].NEffective()); //numerator includes top yield only
		sumweights.back().push_back(allyields[top][ibin].Weight());
		powers.push_back(1); //+1 indicates numerator

		//Denominator
		entries.push_back(vector<float>());
		sumweights.push_back(vector<float>());
		entries.back().push_back(allyields[bkg][ibin+1].NEffective()); //denominator includes all processes
		sumweights.back().push_back(allyields[bkg][ibin+1].Weight());
		powers.push_back(-1); //-1 indicates denominator

		 // Throwing toys to find ratios and uncertainties
		 //if(igraph!=0)
		 val = calcKappa(entries, sumweights, powers, valdown, valup);
		 if(valdown<0) valdown = 0;
		 mct_transfer_factors.push_back(vector<float>({val, valdown, valup}));

      } // Loop over indices


	  for(size_t irow=0; irow<nrows; irow++){ 
	    // allyields[bkg] = yield_table->BackgroundYield(lumi);
	    cout<<bin_names[irow]<<", MC: "    <<setw(7)<<RoundNumber(allyields[bkg][irow].Yield(),3)<<" pm "<<setw(7)<<RoundNumber(allyields[bkg][irow].Uncertainty(), 3);
	    //cout<<"allyields "<<allyields[bkg]  [idens]<<endl;
	    // allyields[top] = yield_table->Yield(proc_top.get(), lumi);
		  cout<<", top: "    <<setw(7)<<RoundNumber(allyields[top][irow].Yield(), 3);
		  // cout<<", signal NC: "    <<setw(7)<<RoundNumber(sigyields[0][irow].Yield(), 3);
		  // cout<<", signal medium: "    <<setw(7)<<RoundNumber(sigyields[1][irow].Yield(), 3);
		  // cout<<", signal compressed: "    <<setw(7)<<RoundNumber(sigyields[2][irow].Yield(), 3)<<endl;
		  if(irow%2==1) cout<<"mCT ratio: "<<setw(7)<<RoundNumber(mct_transfer_factors[(irow-1)/2][0],3)<<endl;	
	    // allyields[w] = yield_table->Yield(proc_wjets.get(), lumi);
	    // allyields[other] = yield_table->Yield(proc_other.get(), lumi);

	  } // Loop over rows

    for(int isig=0;isig<nsig;isig++){
       writeCard(bin_names,allyields,mct_transfer_factors,sig_by_mass[isig],mass_tag[isig],analysis_tag);
    }
	
  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding yields took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;

} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void writeCard(vector<string> bin_names, vector<vector<GammaParams> > allyields,vector<vector<float> > mct_transfer_factors,vector<GammaParams> sigyields,TString mass_tag,TString analysis_tag){

 	TString outpath = Form("statistics/%s/datacards/datacard_%s_.txt",analysis_tag.Data(),mass_tag.Data());
 	// if(boosted) outpath+="boosted_";
 	// if(discovery_mode) outpath+= "significance_";
  // outpath+=mass_tag;
 	// outpath+=".txt";

    cout<<"open "<<outpath<<endl;

    uint nbins = bin_names.size();

    unsigned wname(21), wdist(2), wbin(12);
    for (size_t ibin(0); ibin<nbins; ibin+=2) 
      if(bin_names[ibin].length() > wbin) wbin = bin_names[ibin].length();
    wbin+=1;
    unsigned digit = 2;
    // if (unblind) digit = 0;
    cout<<mct_transfer_factors[0][0]<<endl;
    // --------- write header
    ofstream fcard(outpath);
    int nbg=3;
    int nsyst = 3 + nbins/2 + 3*nbins/2; //one nuisance for every met/njet bin + 3 per bin
    fcard<<"imax "<<nbins/2<<"  number of channels\n";
    fcard<<"jmax "<<nbg<<"  number of backgrounds\n";
    fcard<<"kmax "<<nsyst<<"  number of nuisance parameters\n";
    fcard<<"shapes * * FAKE\n";
    fcard<<endl<<left<<setw(wname)<<"bin"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<bin_names[ibin];
    fcard<<endl<<left<<setw(wname)<<"observation"<<setw(wdist)<<" ";
    
    if(!discovery_mode) for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<RoundNumber(allyields[bkg][ibin].Yield(),digit);
 	  else for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<RoundNumber(allyields[bkg][ibin].Yield()+sigyields[ibin].Yield(),digit);

    fcard<<endl<<endl<<left<<setw(wname)<<"bin"<<setw(wdist)<<" ";
     for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<bin_names[ibin]<<left<<setw(wbin)<<bin_names[ibin]<<left<<setw(wbin)<<bin_names[ibin]<<left<<setw(wbin)<<bin_names[ibin];
    fcard<<endl<<left<<setw(wname)<<"process"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<"sig"<<left<<setw(wbin)<<"top"<<left<<setw(wbin)<<"wjets"<<left<<setw(wbin)<<"other";
    fcard<<endl<<left<<setw(wname)<<"process"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<setw(wbin)<<"0"<<setw(wbin)<<"1"<<setw(wbin)<<"2"<<setw(wbin)<<"3";
    fcard<<endl<<left<<setw(wname)<<"rate"<<setw(wdist)<<" ";
     for (size_t ibin(0); ibin<nbins; ibin+=2) 
        fcard<<setw(wbin)<<Form("%.2f",sigyields[ibin].Yield())<<setw(wbin)<<Form("%.2f",allyields[1][ibin].Yield())<<setw(wbin)<<Form("%.2f",allyields[2][ibin].Yield())<<setw(wbin)<<Form("%.2f",allyields[3][ibin].Yield())<<setw(wbin);
  
    fcard<<endl;
    float sys_lumi = 1.05;
    float sys_filler = 1.15;
    float sys_sig_filler = 1.15;
    unsigned wsyst(14); unsigned wsystype(wname-wsyst);
    fcard<<endl<<left<<setw(wsyst)<<"lumi"<<setw(wsystype)<<"lnN"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<Form("%.2f",sys_lumi)<<left<<setw(wbin)<<"-"<<left<<setw(wbin)<<Form("%.2f",sys_lumi)<<left<<setw(wbin)<<Form("%.2f",sys_lumi);
    fcard<<endl<<left<<setw(wsyst)<<"bkg_flat"<<setw(wsystype)<<"lnN"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<"-"<<left<<setw(wbin)<<Form("%.2f",sys_filler)<<left<<setw(wbin)<<Form("%.2f",sys_filler)<<left<<setw(wbin)<<Form("%.2f",sys_filler);
    fcard<<endl<<left<<setw(wsyst)<<"sig_flat"<<setw(wsystype)<<"lnN"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin+=2) fcard<<left<<setw(wbin)<<Form("%.2f",sys_sig_filler)<<left<<setw(wbin)<<"-"<<left<<setw(wbin)<<"-"<<left<<setw(wbin)<<"-";
    
    //Control region stats
    for (size_t ibin(0); ibin<nbins; ibin+=2){
    	fcard<<endl<<left<<setw(wname)<<Form("mCT_CR_stat%i gmN %.0f",static_cast<int>(ibin),allyields[bkg][ibin+1].Yield())<<setw(wdist)<<" ";
    	for(size_t j(0);j<1+(nbg+1)*(ibin/2);j++) fcard<<left<<setw(wbin)<<"-";
    	fcard<<left<<setw(wbin)<<Form("%.3f",mct_transfer_factors[ibin/2][0]);
    	for(size_t j(0);j<(nbg+1)*nbins/2 - (2+(nbg+1)*(ibin/2));j++) fcard<<left<<setw(wbin)<<"-";
    }
    //Uncorrelated BG systematics
    for (size_t ibin(0); ibin<nbins; ibin+=2){
      for(int k(1);k<=nbg;k++){
      fcard<<endl<<left<<setw(wname)<<Form("bkg_flat%i_proc%i lnN",static_cast<int>(ibin),k)<<setw(wdist)<<" ";
      
       for(size_t j(0);j<k+(nbg+1)*(ibin/2);j++) fcard<<left<<setw(wbin)<<"-";
       fcard<<left<<setw(wbin)<<Form("%.2f",sys_filler);
       for(size_t j(0);j<(nbg+1)*nbins/2 - (k+1+(nbg+1)*(ibin/2));j++) fcard<<left<<setw(wbin)<<"-";

      }
    }

    fcard<<endl<<endl;



    cout<<"Wrote headers"<<endl;
    fcard.close();

}
