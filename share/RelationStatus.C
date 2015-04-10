/*
 *  $Id: RelationStatus.C, 2015-04-10 10:20:49 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 10/04/2015
*/

#include <iostream>

#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"

#include "DmpLoadParameters.h"

void RelationStatus(TString fname,int entries_cut = 100)
{
   DmpParameterSteering  steering;
   DmpParameterHolder  rel_par;
   if(fname.Contains("Bgo")) DAMPE::Bgo::LoadRelation((std::string)fname,rel_par,steering);
   if(fname.Contains("Psd")) DAMPE::Psd::LoadRelation((std::string)fname,rel_par,steering);

   TH2F*    h_p0p1[2][2];
   TH2F*    h_chi[2][2];
   TH2F*    h_pileupEvts[2][2];
   TH2F*    h_pileupEvtsRatio[2][2];
   TCanvas *c[4];

   TString name(fname);
   name.Remove(name.Last('-'));
   for(short s=0;s<2;++s){
     for(short dy =0;dy<2;++dy){
       TString xx= name+Form("Side%d_Dy%d-Dy%d",s,dy*3+2,(dy+1)*3+2);
       c[s*2+dy+1] =new TCanvas(xx,xx);
       h_p0p1[s][dy] = new TH2F(xx+" p0 VS P1","p0 VS p1;p0;p1",3000,0.015,0.045,200,-20,40);
       h_p0p1[s][dy]->SetMarkerStyle(23);
       h_p0p1[s][dy]->SetMarkerSize(0.8);
       if(fname.Contains("Bgo")){
         h_chi[s][dy]  = new TH2F(xx+" chi2/NDF","chi2/NDF;bar ID;layer ID",22,0,22,14,0,14);
         h_pileupEvts[s][dy] = new TH2F(xx+" pileup events","pileup events;bar ID;layer ID",22,0,22,14,0,14);
         h_pileupEvtsRatio[s][dy] = new TH2F(xx+" pileup events ratio","pileup ratio;bar ID;layer ID",22,0,22,14,0,14);
       }
       if(fname.Contains("Psd")){
         h_chi[s][dy]  = new TH2F(xx+" chi2/NDF","chi2/NDF;bar ID;layer ID",41,0,41,2,0,2);
         h_pileupEvts[s][dy] = new TH2F(xx+" pileup events","pileup events;bar ID;layer ID",41,0,41,2,0,2);
         h_pileupEvtsRatio[s][dy] = new TH2F(xx+" pileup events ratio","pileup ratio;bar ID;layer ID",41,0,41,2,0,2);
       }
     }
   }
   //gStyle->SetOptStat(000000);
   gStyle->SetOptStat("emru");

   int entries = 0, errorEvts = 0;
   short lid =0,bid = 0,sid = 0,dyid =0;
   for(DmpParameterHolder::iterator it = rel_par.begin(); it != rel_par.end();++it){
     entries = it->second.at(7);
     if(entries < entries_cut) continue;
     lid  = it->second.at(0);
     bid  = it->second.at(1);
     sid  = it->second.at(2);
     dyid = (it->second.at(3) - 2) /3;
     h_p0p1[sid][dyid]->Fill(it->second.at(5),it->second.at(4));
     h_chi[sid][dyid]->Fill(bid,lid,it->second.at(6));
     errorEvts = it->second.at(8);
     h_pileupEvts[sid][dyid]->Fill(bid,lid,errorEvts);
     h_pileupEvtsRatio[sid][dyid]->Fill(bid,lid,(double)errorEvts /  (double)entries * 100);
   }

   for(short s=0;s<2;++s){
     for(short d=0;d<2;++d){
       c[s*2+d+1]->Divide(2,2);
       c[s*2+d+1]->cd(1);
       h_p0p1[s][d]->Draw();
       c[s*2+d+1]->cd(2);
       h_pileupEvts[s][d]->Draw("colz");
       c[s*2+d+1]->cd(3);
       h_chi[s][d]->Draw("colz");
       c[s*2+d+1]->cd(4);
       h_pileupEvtsRatio[s][d]->Draw("colz");
     }
   }

}

