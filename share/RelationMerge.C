#include <iostream>
#include "DmpLoadParameters.h"
#include <fstream>
#include <vector>
#include "TMath.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"


void Split(string fname)
{
   DmpParameterSteering  steering;
   DmpParameterHolder  rel_par;
   DAMPE::LoadParameters(fname,rel_par,steering);
   string gfName = "good_"+fname;
   ofstream goodF(gfName.c_str(),std::ios::out);
   string ffName = "failed_"+fname;
   ofstream failedF(ffName.c_str(),std::ios::out);

   for(DmpParameterHolder::iterator it = rel_par.begin();it != rel_par.end();++it){
     if(it->second.at(8) == -1){
       failedF<<it->first;
       for(size_t i=0;i < it->second.size();++i){
         failedF<<"\t\t"<<it->second.at(i);
       }
       failedF<<std::endl;
     }else{
       goodF<<it->first;
       for(size_t i=0;i < it->second.size();++i){
         goodF<<"\t\t"<<it->second.at(i);
       }
       goodF<<std::endl;
     }
   }
}

void FindBest(std::vector<std::string> filenames,string outPutFileName="allGood.txt")
{
  int nFiles = filenames.size();
  std::vector<DmpParameterHolder>  pars(nFiles);
  DmpParameterSteering  steering;
  for(int i=0;i<nFiles;++i){
    DAMPE::LoadParameters(filenames[i],pars.at(i),steering);
  }
  ofstream goodF(outPutFileName.c_str(),std::ios::out);
  for(DmpParameterHolder::iterator it = pars[0].begin();it!=pars[0].end();++it){
    for(int i=0;i<nFiles;++i){
      if(pars[i][it->first].at(8) == -1) continue;
      if(TMath::Abs(pars[i][it->first].at(4)) > 40) continue;
      if(pars[i][it->first].at(5) == -999) continue;
      if(pars[i][it->first].at(6) == -999) continue;
      if(it->second.at(8) == -1){
        it->second = pars[i][it->first];
        continue;
      }
      if(it->second.at(6) == -999 || it->second.at(6) == 999){
        it->second = pars[i][it->first];
        continue;
      }
      bool comp0 = TMath::Abs(pars[i][it->first].at(4)) < TMath::Abs(it->second.at(4)); // p0
      bool comp1 = pars[i][it->first].at(5) < it->second.at(5); // p1
      double deltchi = TMath::Abs(TMath::Abs(pars[i][it->first].at(6))-1) - TMath::Abs(TMath::Abs(it->second.at(6)) - 1); // chi2
      bool  comp2 = deltchi < 0? 1: 0;
      bool comp3 = pars[i][it->first].at(7) > 0.666 * it->second.at(7); // entries
      bool comp4 = (pars[i][it->first].at(7) > 3*it->second.at(7)) && (deltchi < 1.5);
      if((comp0 && comp1 && comp3) || (comp0 && comp2 && comp3) || (comp1 && comp2 && comp3) || comp4){
        it->second = pars[i][it->first];
      }
    }
    goodF<<it->first;
    for(size_t n=0;n<it->second.size();++n){
      goodF<<"\t\t"<<it->second.at(n);
    }
    goodF<<std::endl;
  }
}

void Static(string fname)
{
  DmpParameterSteering  steering;
  DmpParameterHolder  rel_par;
  DAMPE::LoadParameters((std::string)fname,rel_par,steering);
  TString fn = fname;
  fn.Remove(fn.Last('.'));
  fname = fn+" p0 VS p1";
  TH2F *h2_p0p1 = new TH2F(fname.c_str(),"p0 VS p1;p1;p0",400,0,0.04,1200,-60,60);
  fname = fn+" chi2";
  TH2F *map_chi2 = new TH2F(fname.c_str(),"chi2;bar;layer",22,0,22,14,0,14);
  fname = fn+" p0";
  TH2F *map_p0 = new TH2F(fname.c_str(),"p0;bar;layer",22,0,22,14,0,14);
  fname = fn+" p1";
  TH2F *map_p1 = new TH2F(fname.c_str(),"p1;bar;layer",22,0,22,14,0,14);
  TString isPSD = fn;
  isPSD.ToUpper();
  if(isPSD.Contains("PSD")){
    delete map_chi2;
    fname = fn+" chi2";
    map_chi2 = new TH2F(fname.c_str(),"chi2;bar;layer",42,0,42,2,0,2);
    delete map_p0;
    fname = fn+" p0";
    map_p0 = new TH2F(fname.c_str(),"p0;bar;layer",42,0,42,2,0,2);
    delete map_p1;
    fname = fn+" p1";
    map_p1 = new TH2F(fname.c_str(),"p1;bar;layer",42,0,42,2,0,2);
  }

  for(DmpParameterHolder::iterator it = rel_par.begin();it != rel_par.end();++it){
    h2_p0p1->Fill(it->second.at(5),it->second.at(4));
    map_p0->Fill(it->second.at(1),it->second.at(0),it->second.at(4));
    map_p1->Fill(it->second.at(1),it->second.at(0),it->second.at(5));
    map_chi2->Fill(it->second.at(1),it->second.at(0),it->second.at(6));
  }

  gStyle->SetOptStat(111111);
  TCanvas *c1 = new TCanvas(fn,fn,600,800);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetGrid();
  h2_p0p1->SetMarkerStyle(23);
  h2_p0p1->SetMarkerSize(0.8);
  h2_p0p1->Draw();
  c1->cd(2);
  map_chi2->Draw("colz");
  c1->cd(3);
  map_p0->Draw("colz");
  c1->cd(4);
  map_p1->Draw("colz");
}

