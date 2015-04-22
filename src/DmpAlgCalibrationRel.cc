/*
 *  $Id: DmpAlgCalibrationRel.cc, 2015-04-01 10:03:59 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

#include <algorithm>

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TMath.h"

#include "DmpEvtBgoRaw.h"
#include "DmpAlgCalibrationRel.h"
#include "DmpDataBuffer.h"
#include "DmpLoadParameters.h"
#include "DmpBgoBase.h"
#include "DmpPsdBase.h"
#include "DmpParameterBgo.h"
#include "DmpParameterPsd.h"
#include "DmpCore.h"
#include "DmpTimeConvertor.h"
#define  NBinX  3200
#define  NBinY  90
//#define  ChiSCut 5
#define  XLowCut 3000
#define  EntriesCut  80

//-------------------------------------------------------------------
DmpAlgCalibrationRel::DmpAlgCalibrationRel()
 :DmpVAlg("DmpAlgCalibrationRel"),
  fEvtBgo(0),
  fEvtPsd(0)
{
  gRootIOSvc->SetOutFileKey("CalRel");

}

//-------------------------------------------------------------------
DmpAlgCalibrationRel::~DmpAlgCalibrationRel(){
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationRel::Initialize(){
  // read input data
  std::string inpath = "Event/Rec0/";
  fEvtBgo = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject(inpath+"Bgo"));
  if(0 == fEvtBgo){
    fEvtBgo = new DmpEvtBgoRaw();
    gDataBuffer->LinkRootFile(inpath+"Bgo",fEvtBgo);
  }
  fEvtPsd = dynamic_cast<DmpEvtPsdRaw*>(gDataBuffer->ReadObject(inpath+"Psd"));
  if(0 == fEvtPsd){
    fEvtPsd = new DmpEvtPsdRaw();
    gDataBuffer->LinkRootFile(inpath+"Psd",fEvtPsd);
  }

  char name[200];
  // Bgo create Hist map
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        snprintf(name,200,"Bgo_L%02d_B%02d_S%d_Dy2-Dy5",l,b,s);
        fBgoRelHist[l][b][s][0] = new TH2D(name,name,NBinX,0,16000,NBinY,0,450);
        snprintf(name,200,"Bgo_L%02d_B%02d_S%d_Dy5-Dy8",l,b,s);
        fBgoRelHist[l][b][s][1] = new TH2D(name,name,NBinX,0,16000,NBinY,0,450);
      }
    }
  }
  // Psd
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s=0;s<DmpParameterPsd::kSideNo;++s){
        snprintf(name,200,"Psd_L%02d_B%02d_S%d_Dy5-Dy8",l,b,s);
        fPsdRelHist[l][b][s] = new TH2D(name,name,NBinX,0,16000,NBinY,0,450);
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationRel::ProcessThisEvent(){
  // bgo bar
  static short gid_dy5=0,gid_dy2 = 0,gid_dy8=0;
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b = 0;b<DmpParameterBgo::kBarNo;++b){
      for(short s = 0;s<2;++s){
        gid_dy5 = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,5);
        if(fEvtBgo->fADC.find(gid_dy5) == fEvtBgo->fADC.end()){
          continue;
        }
        gid_dy8 = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        if(fEvtBgo->fADC.find(gid_dy8) != fEvtBgo->fADC.end()){
          fBgoRelHist[l][b][s][1]->Fill(fEvtBgo->fADC[gid_dy8],fEvtBgo->fADC[gid_dy5]);
        }
        gid_dy2 = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,2);
        if(fEvtBgo->fADC.find(gid_dy2) != fEvtBgo->fADC.end()){
          fBgoRelHist[l][b][s][0]->Fill(fEvtBgo->fADC[gid_dy5],fEvtBgo->fADC[gid_dy2]);
        }
      }
    }
  }

  // Psd Rel bar
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b = 0;b<DmpParameterPsd::kStripNo;++b){
      for(short s = 0;s<2;++s){
        gid_dy5 = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,5);
        if(fEvtPsd->fADC.find(gid_dy5) == fEvtPsd->fADC.end()){
          continue;
        }
        gid_dy8 = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        if(fEvtPsd->fADC.find(gid_dy8) != fEvtPsd->fADC.end()){
          fPsdRelHist[l][b][s]->Fill(fEvtPsd->fADC[gid_dy8],fEvtPsd->fADC[gid_dy5]);
        }
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationRel::Finalize(){
  TF1 *lxg_f = new TF1("linear","pol1",XLowCut,13000);
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_DyRelationHist.root";
  TFile *histFile = gRootIOSvc->GetOutputRootFile();//new TFile(histFileName.c_str(),"RECREATE");
  std::string epsFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_DyRelationHist.eps";
  //TPostScript *epsF = new TPostScript(epsFileName.c_str(),113);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("e");
  TCanvas *c0 = new TCanvas("c0","c0",0,0,600,800);
  std::string epsFileName0 = epsFileName+"[";
  c0->Print(epsFileName0.c_str());

  // create output txtfile      BGO
  histFile->mkdir("Bgo");
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".rel";
  o_RelData_Bgo.open(name.c_str(),std::ios::out);
  o_RelData_Bgo<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_RelData_Bgo<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<std::endl;
  o_RelData_Bgo<<Mark_D<<std::endl;
  TCanvas *c[2][2];
  for(int i=0;i<4;++i){
    c[i/2][i%2] = new TCanvas(Form("c%d",i),Form("c%d",i),600,800);
    c[i/2][i%2]->Divide(4,6);
  }
  double p0,p1,chi2;
  int fitStat;
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s = 0;s<DmpParameterBgo::kSideNo;++s){
        for(short nd=0;nd<2;++nd){
          o_RelData_Bgo<<DmpBgoBase::ConstructGlobalDynodeID(l,b,s,nd*3+2)<<"\t\t"<<Form("%d\t\t%d\t\t%d\t\t%d",l,b,s,nd*3+2);
          c[s][nd]->cd(b+1);
          if(fBgoRelHist[l][b][s][nd]->GetEntries() > EntriesCut){
            fitStat = fBgoRelHist[l][b][s][nd]->Fit(lxg_f,"RQB");
            fBgoRelHist[l][b][s][nd]->Write();
            p0 = lxg_f->GetParameter(0);
            p1 = lxg_f->GetParameter(1);
            chi2 = lxg_f->GetChisquare() / lxg_f->GetNDF();
            if(TMath::Abs(chi2) > 5){
              p0=-999;
              p1=-999;
              //chi2 = -999;
              fitStat = -1;
            }
          }else{
            p0=-999;
            p1=-999;
            chi2 = -999;
            fitStat = -1;
          }
          fBgoRelHist[l][b][s][nd]->Draw("colz");
          o_RelData_Bgo<<"\t\t"<<p0<<"\t\t"<<p1<<"\t\t"<<chi2<<"\t\t"<<fBgoRelHist[l][b][s][nd]->GetEntries()<<"\t\t"<<fitStat<<std::endl;
          //delete fBgoRelHist[l][b][s][nd];
        }
      }
    }
    for(int i=0;i<4;++i){
      c[i/2][i%2]->Print(epsFileName.c_str(),"Portrait");//Update();
    }
  }
  o_RelData_Bgo<<Mark_N<<std::endl;
  o_RelData_Bgo.close();

  // create output txtfile      PSD
  histFile->mkdir("Psd");
  histFile->cd("Psd");
  name = "Psd_"+gRootIOSvc->GetInputStem()+".rel";
  o_RelData_Psd.open(name.c_str(),std::ios::out);
  o_RelData_Psd<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_RelData_Psd<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<std::endl;
  o_RelData_Psd<<Mark_D<<std::endl;
  TCanvas *cpsd[2];
  for(int i=0;i<2;++i){
    cpsd[i] = new TCanvas(Form("cPsd%d",i),"",600,800);
    cpsd[i]->Divide(6,7);
  }
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s = 0;s<DmpParameterPsd::kSideNo;++s){
        o_RelData_Psd<<DmpPsdBase::ConstructGlobalDynodeID(l,b,s,5)<<"\t\t"<<Form("%d\t\t%d\t\t%d\t\t5",l,b,s);
        cpsd[s]->cd(b+1);
        if(fPsdRelHist[l][b][s]->GetEntries() > EntriesCut){
          fPsdRelHist[l][b][s]->Fit(lxg_f,"RQB");
          fPsdRelHist[l][b][s]->Write();
          p0 = lxg_f->GetParameter(0);
          p1 = lxg_f->GetParameter(1);
          chi2 = lxg_f->GetChisquare() / lxg_f->GetNDF();
          if(TMath::Abs(chi2) > 5){
            p0=-999;
            p1=-999;
            //chi2 = -999;
            fitStat = -1;
          }
        }else{
          p0=-999;
          p1=-999;
          chi2 = -999;
          fitStat = -1;
        }
        fPsdRelHist[l][b][s]->Draw("colz");
        o_RelData_Psd<<"\t\t"<<p0<<"\t\t"<<p1<<"\t\t"<<chi2<<"\t\t"<<fPsdRelHist[l][b][s]->GetEntries()<<"\t\t"<<fitStat<<std::endl;
        //delete fPsdRelHist[l][b][s];
      }
    }
    cpsd[0]->Print(epsFileName.c_str(),"Portrait");//Update();
    cpsd[1]->Print(epsFileName.c_str(),"Portrait");//Update();
  }
  o_RelData_Psd<<Mark_N<<std::endl;
  o_RelData_Psd.close();
  TCanvas *c2 = new TCanvas("c2","c2");
  std::string epsFileName1 = epsFileName+"]";
  c2->Print(epsFileName1.c_str());

  return true;
}

