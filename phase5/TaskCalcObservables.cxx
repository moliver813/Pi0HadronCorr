#if defined(__CINT__)
#define _SYS_TYPES_H_
#endif
#define DTR TMath::DegToRad()
// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TSpectrum.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>
#include <TLatex.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TMath.h>
#include <TMinuit.h>
#include <THnSparse.h>
#include <TDatime.h>

// C,C++ Stuff
#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

#include "TaskCalcObservablesGraphicsTools.cxx"
#include "TaskCalcObservablesMathTools.cxx"
#include "TaskCalcObservables.h"

using namespace std;

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;
/// \cond CLASSIMP
ClassImp(TaskCalcObservables);


TaskCalcObservables::TaskCalcObservables() {
  gROOT->SetBatch(kTRUE);

  SetStyle();
}


void TaskCalcObservables::SetStyle() {
//  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(53);  //standard is 1
  gStyle->SetCanvasColor(10);
  //  TGaxis::SetMaxDigits(4);  //..ELI I don't remember why I wanted 4 changed to 2
  //gStyle->SetPadTopMargin(0.07);//0.05
  //gStyle->SetPadBottomMargin(0.18);//0.15
  //  gStyle->SetPadRightMargin(0.045);
  gStyle->SetPadRightMargin(0.08);
  //gStyle->SetPadLeftMargin(0.21);
  gStyle->SetFrameFillColor(10);
//  gStyle->SetLabelSize(0.05,"X");
//  gStyle->SetLabelSize(0.05,"Y");
//  gStyle->SetTitleSize(0.07,"X");
//  gStyle->SetTitleSize(5.0,"Y");
  TGaxis::SetMaxDigits(3); //2
  gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas
}

void TaskCalcObservables::LoadHistograms() {
  cout<<"Loading Histograms"<<endl;

  // Loading our observable settings:
  TH1D * VariableInfo = (TH1D *) fInputFileCentral->Get("VariableInfo");
  if (VariableInfo) {
    fObservable = VariableInfo->GetBinContent(1);
  } else {
    cout<<"No Variable Input TH1D found in Input File All."<<endl;
    cout<<"Using default values for Observable Info and name."<<endl;
    //fObservable=1; // the Z_t
    fObservable=2; // the pTa
  }
  if (fObservable==0)      fObservableName = "p_{T} (GeV/#it{c})";
  else if (fObservable==1) fObservableName = "z_{T}";
  else if (fObservable==2) fObservableName = "p_{T}^{a} (GeV/#it{c})";


//dPhi_RPFMethod0_Full_AllEP_RPFSub_ObsBin0
//dPhi_RPFMethod0_Full_EP0_RPFSub_ObsBin0

  

  nObsBins= 0;
  // Check the number of observable bins safely.
  for (int i = 0; i < 13; i++) {
    TH1D * fLocal = 0;
    TString fLocalName = Form("dPhi_RPFMethod%d_Full_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);

  // EP0,1,2 use the rescale histograms
// dPhi_RPFMethod1_NearEta_EP2_RPFSub_ObsBin1_Rescale

    fLocal  = (TH1D *) fInputFileCentral->Get(fLocalName);
    if (!fLocal) {
      printf("  Histo %s not found. End of search.\n",fLocalName.Data());
      break;
    }
    nObsBins++;
  }


  printf("  Observable Bins: %d\n",nObsBins);

  // Loading the nearside Delta Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
//      fLocalName = Form("Proj_PtBin%d_EP-1_NearSideDEta_ObsBin%d",iPtBin,i);

    }
  }




  // Loading the Full Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    vector<TH2F *> fLocalVectorVar = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_Full_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_Full_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);

      //TString fLocalNameVar = Form("%s_Variants",fLocalName.Data());
      TString fLocalNameVar = Form("dPhi_RPFMethod%d_Full_EP%d_RPFSub_ObsBin%d_Variants_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalNameVar = Form("dPhi_RPFMethod%d_Full_AllEP_RPFSub_ObsBin%d_Variants",iRPFMethod,i);
      TH2F * fLocalVariant = 0;
      fLocalVariant = (TH2F *) fInputFileCentral->Get(fLocalNameVar);
      if (!fLocalVariant) {
        fprintf(stderr,"Missing histogram %s\n",fLocalNameVar.Data());
      }
      fLocalVectorVar.push_back(fLocalVariant);

    }
    fFullDPhiProj_Sub.push_back(fLocalVector);
    fFullDPhiProj_Sub_RPFVar.push_back(fLocalVectorVar);
  }
  // Loading the Near Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    vector<TH2F *> fLocalVectorVar = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_NearEta_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_NearEta_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);

      //TString fLocalNameVar = Form("%s_Variants",fLocalName.Data());
      TString fLocalNameVar = Form("dPhi_RPFMethod%d_NearEta_EP%d_RPFSub_ObsBin%d_Variants_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalNameVar = Form("dPhi_RPFMethod%d_NearEta_AllEP_RPFSub_ObsBin%d_Variants",iRPFMethod,i);
      TH2F * fLocalVariant = 0;
      fLocalVariant = (TH2F *) fInputFileCentral->Get(fLocalNameVar);
      if (!fLocalVariant) {
        fprintf(stderr,"Missing histogram %s\n",fLocalNameVar.Data());
      }
      fLocalVectorVar.push_back(fLocalVariant);
    }
    fNearEtaDPhiProj_Sub.push_back(fLocalVector);
    fNearEtaDPhiProj_Sub_RPFVar.push_back(fLocalVectorVar);
  }
  // Loading the Far Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    vector<TH2F *> fLocalVectorVar = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_FarEta_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_FarEta_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);

      //TString fLocalNameVar = Form("%s_Variants",fLocalName.Data());
      TString fLocalNameVar = Form("dPhi_RPFMethod%d_FarEta_EP%d_RPFSub_ObsBin%d_Variants_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalNameVar = Form("dPhi_RPFMethod%d_FarEta_AllEP_RPFSub_ObsBin%d_Variants",iRPFMethod,i);
      TH2F * fLocalVariant = 0;
      fLocalVariant = (TH2F *) fInputFileCentral->Get(fLocalNameVar);
      if (!fLocalVariant) {
        fprintf(stderr,"Missing histogram %s\n",fLocalNameVar.Data());
      }
      fLocalVectorVar.push_back(fLocalVariant);
    }
    fFarEtaDPhiProj_Sub.push_back(fLocalVector);
    fFarEtaDPhiProj_Sub_RPFVar.push_back(fLocalVectorVar);
  }

  if (fFullDPhiProj_Sub_RPFVar[0][0] != 0) {
    nRPFVariants = fFullDPhiProj_Sub_RPFVar[0][0]->GetNbinsY();
  } else {
    printf("No Variants found\n");
  }

}

void TaskCalcObservables::InitArrays() {

  Double_t fZtStep = 1.0/(7 - 1.0);
  Double_t fXiStep = 2.5/(8 - 1.0);

  Double_t * ObsArray = 0;

  Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,20,23,30,60};
  Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0.03,0.08,0.16,0.29,0.5,0.84,1.39,2.};
  //Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
  Double_t array_XI_BinsValue[kXiNBINS+1]   ={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};
  Double_t array_HPT_BinsValue[kNoHPtBins+1]={0.2,0.4,0.8,1.5,2.5,4,7,11,17};


  if (fObservable == 0) ObsArray = array_G_BinsValue;
  else if (fObservable == 1) ObsArray = array_ZT_BinsValue;
  //else if (fObservable == 2) ObsArray = array_XI_BinsValue;
  else if (fObservable == 2) ObsArray = array_HPT_BinsValue;

  for (Int_t i = 0; i <= nObsBins; i++) {
    fObsBins.push_back(ObsArray[i]);
  }


  fNSRmsEP = {};
  fNSYieldsEP = {};

  fNSYieldsInc = new TGraphErrors(nObsBins);
  fNSYieldsInc->SetName("NSYieldsInc");
  fNSYieldsInc->SetTitle("Near-side Yields");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("NSYieldsEP%d",i));
    newGraph->SetTitle(Form("Near-side Yields (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fNSYieldsEP.push_back(newGraph);
  }
  fASYieldsInc = new TGraphErrors(nObsBins);
  fASYieldsInc->SetName("ASYieldsInc");
  fASYieldsInc->SetTitle("Away-side Yields");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("ASYieldsEP%d",i));
    newGraph->SetTitle(Form("Away-side Yields (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fASYieldsEP.push_back(newGraph);
  }

  fNSRmsInc = new TGraphErrors(nObsBins);
  fNSRmsInc->SetName("NSRmsInc");
  fNSRmsInc->SetTitle("Near-side width");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("NSRmsEP%d",i));
    newGraph->SetTitle(Form("Near-side RMS (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fNSRmsEP.push_back(newGraph);
  }
  fASRmsInc = new TGraphErrors(nObsBins);
  fASRmsInc->SetName("ASRmsInc");
  fASRmsInc->SetTitle("Away-side width");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("ASRmsEP%d",i));
    newGraph->SetTitle(Form("Away-side RMS (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fASRmsEP.push_back(newGraph);
  }

  // Initialize the variant TGraphErrors
  if (fFullDPhiProj_Sub_RPFVar[0][0] != 0) {
    nRPFVariants = fFullDPhiProj_Sub_RPFVar[0][0]->GetNbinsY();
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {

      vector<TGraphErrors *> fNSRmsEPVar = {};
      vector<TGraphErrors *> fNSYieldsEPVar = {};
      vector<TGraphErrors *> fASRmsEPVar = {};
      vector<TGraphErrors *> fASYieldsEPVar = {};

      TGraphErrors * fNSYieldsIncVariant = new TGraphErrors(nObsBins);
      fNSYieldsIncVariant->SetName(Form("NSYieldsInc_RPFVar%d",iVar));
      fNSYieldsIncVariant->SetTitle("Near-side Yields");
      fNSYieldsInc_RPFVariants.push_back(fNSYieldsIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("NSYieldsEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Near-side Yields (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fNSYieldsEPVar.push_back(newGraph);
      }
      fNSYieldsEP_RPFVariants.push_back(fNSYieldsEPVar);

      TGraphErrors * fASYieldsIncVariant = new TGraphErrors(nObsBins);
      fASYieldsIncVariant->SetName(Form("ASYieldsInc_RPFVar%d",iVar));
      fASYieldsIncVariant->SetTitle("Away-side Yields");
      fASYieldsInc_RPFVariants.push_back(fASYieldsIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("ASYieldsEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Away-side Yields (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fASYieldsEPVar.push_back(newGraph);
      }
      fASYieldsEP_RPFVariants.push_back(fASYieldsEPVar);

      TGraphErrors * fNSRmsIncVariant = new TGraphErrors(nObsBins);
      fNSRmsIncVariant->SetName(Form("NSRmsInc_RPFVar%d",iVar));
      fNSRmsIncVariant->SetTitle("Near-side width");
      fNSRmsInc_RPFVariants.push_back(fNSRmsIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("NSRmsEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Near-side RMS (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fNSRmsEPVar.push_back(newGraph);
      }
      fNSRmsEP_RPFVariants.push_back(fNSRmsEPVar);
      TGraphErrors * fASRmsIncVariant = new TGraphErrors(nObsBins);
      fASRmsIncVariant->SetName(Form("ASRmsInc_RPFVar%d",iVar));
      fASRmsIncVariant->SetTitle("Away-side width");
      fASRmsInc_RPFVariants.push_back(fASRmsIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("ASRmsEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Away-side RMS (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fASRmsEPVar.push_back(newGraph);
      }
      fASRmsEP_RPFVariants.push_back(fASRmsEPVar);

    }
  }

}



void TaskCalcObservables::DrawDirectComparisons() {
  TCanvas * cDirCompareCanvas = new TCanvas("DirCompareCanvas","DirCompareCanvas");
  TLegend * legDirCompareNearside = new TLegend(0.7,0.55,0.93,0.85);
  TLegend * legDirCompareAwayside = new TLegend(0.7,0.65,0.93,0.85);




  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    fFullDPhiProj_Sub[iObsBin][kNEPBins]->SetTitle(Form("%.1f #leq p_{T}^{assoc} < %.1f GeV/#it{c}",fObsBins[iObsBin],fObsBins[iObsBin+1]));
    fNearEtaDPhiProj_Sub[iObsBin][kNEPBins]->SetTitle(Form("%.1f #leq p_{T}^{assoc} < %.1f GeV/#it{c}",fObsBins[iObsBin],fObsBins[iObsBin+1]));
    fFarEtaDPhiProj_Sub[iObsBin][kNEPBins]->SetTitle(Form("%.1f #leq p_{T}^{assoc} < %.1f GeV/#it{c}",fObsBins[iObsBin],fObsBins[iObsBin+1]));
    for (int iEP = 0; iEP < kNEPBins+1; iEP++) {

      // Test
      // Would make sense to do this elsewhere, function description does not imply changes
      fFullDPhiProj_Sub[iObsBin][iEP]->Rebin(nRebinDPhi);
      fNearEtaDPhiProj_Sub[iObsBin][iEP]->Rebin(nRebinDPhi);
      fFarEtaDPhiProj_Sub[iObsBin][iEP]->Rebin(nRebinDPhi);

      if (nRebinDPhi != 1) {
        fFullDPhiProj_Sub[iObsBin][iEP]->Scale(1./nRebinDPhi);
        fNearEtaDPhiProj_Sub[iObsBin][iEP]->Scale(1./nRebinDPhi);
        fFarEtaDPhiProj_Sub[iObsBin][iEP]->Scale(1./nRebinDPhi);
      }


      fFullDPhiProj_Sub[iObsBin][iEP]->SetMarkerStyle(kEPMarkerList[iEP]);
      fNearEtaDPhiProj_Sub[iObsBin][iEP]->SetMarkerStyle(kEPMarkerList[iEP]);
      fFarEtaDPhiProj_Sub[iObsBin][iEP]->SetMarkerStyle(kEPMarkerList[iEP]);

      fFullDPhiProj_Sub[iObsBin][iEP]->SetMarkerColor(kEPColorList[iEP]);
      fNearEtaDPhiProj_Sub[iObsBin][iEP]->SetMarkerColor(kEPColorList[iEP]);
      fFarEtaDPhiProj_Sub[iObsBin][iEP]->SetMarkerColor(kEPColorList[iEP]);

      fFullDPhiProj_Sub[iObsBin][iEP]->SetLineColor(kEPColorList[iEP]);
      fNearEtaDPhiProj_Sub[iObsBin][iEP]->SetLineColor(kEPColorList[iEP]);
      fFarEtaDPhiProj_Sub[iObsBin][iEP]->SetLineColor(kEPColorList[iEP]);
    }
  }


  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {

    // Draw NearSide (nearEta, -pi/2 < DeltaPhi < pi/2) with different event plane bins


    std::vector<TH1D * > fGraphHists = {};

    legDirCompareNearside->Clear();
 //   legDirCompare->SetHeader(Form("%.1f - %.1f",fObsBins[iObsBin+1],fObsBins[iObsBin]));
    fNearEtaDPhiProj_Sub[iObsBin][kNEPBins]->Draw();
    fGraphHists.push_back(fNearEtaDPhiProj_Sub[iObsBin][kNEPBins]);

    legDirCompareNearside->AddEntry(fNearEtaDPhiProj_Sub[iObsBin][kNEPBins],"Inclusive","lp");
    // Set Limits

    fNearEtaDPhiProj_Sub[iObsBin][kNEPBins]->GetXaxis()->SetRangeUser(-TMath::Pi()/2,TMath::Pi()/2);

    for (int iEP = 0; iEP < kNEPBins; iEP++) {
      fNearEtaDPhiProj_Sub[iObsBin][iEP]->Draw("SAME");
      fGraphHists.push_back(fNearEtaDPhiProj_Sub[iObsBin][iEP]);
      legDirCompareNearside->AddEntry(fNearEtaDPhiProj_Sub[iObsBin][iEP],fEPBinTitles[iEP],"lp");
    }
    ZoomCommonMinMax(fGraphHists);
    legDirCompareNearside->Draw("SAME");
    cDirCompareCanvas->Print(Form("%s/CompareDPhi_ObsBin%d_NearSide.pdf",fOutputDir.Data(),iObsBin));
    cDirCompareCanvas->Print(Form("%s/CompareDPhi_ObsBin%d_NearSide.png",fOutputDir.Data(),iObsBin));
    cDirCompareCanvas->Print(Form("%s/CFiles/CompareDPhi_ObsBin%d_NearSide.C",fOutputDir.Data(),iObsBin));

    // Draw AwaySide (FullEta, pi/2 < DeltaPhi < 3pi/2) with different event plane bins
    fGraphHists = {};
    cDirCompareCanvas->Clear();
    legDirCompareAwayside->Clear();
    fFullDPhiProj_Sub[iObsBin][kNEPBins]->Draw();
    fGraphHists.push_back(fFullDPhiProj_Sub[iObsBin][kNEPBins]);
    legDirCompareAwayside->AddEntry(fFullDPhiProj_Sub[iObsBin][kNEPBins],"Inclusive","lp");
    fFullDPhiProj_Sub[iObsBin][kNEPBins]->GetXaxis()->SetRangeUser(TMath::Pi()/2,3*TMath::Pi()/2-1e-3);
    for (int iEP = 0; iEP < kNEPBins; iEP++) {
      fFullDPhiProj_Sub[iObsBin][iEP]->Draw("SAME");
      fGraphHists.push_back(fFullDPhiProj_Sub[iObsBin][iEP]);
      legDirCompareAwayside->AddEntry(fFullDPhiProj_Sub[iObsBin][iEP],fEPBinTitles[iEP],"lp");
    }
    ZoomCommonMinMax(fGraphHists);
    legDirCompareAwayside->Draw("SAME");
    cDirCompareCanvas->Print(Form("%s/CompareDPhi_ObsBin%d_AwaySide.pdf",fOutputDir.Data(),iObsBin));
    cDirCompareCanvas->Print(Form("%s/CompareDPhi_ObsBin%d_AwaySide.png",fOutputDir.Data(),iObsBin));
    cDirCompareCanvas->Print(Form("%s/CFiles/CompareDPhi_ObsBin%d_AwaySide.C",fOutputDir.Data(),iObsBin));





  }


}




void TaskCalcObservables::CalculateResults() {
  TCanvas * cCalcQACanvas = new TCanvas("CalcQACanvas");

  // Calculate the yields and RMSs
  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    for (int iEPBin = 0; iEPBin < kNEPBins+1; iEPBin++) {
      printf("Preparing to calculate results for observable bin %d, event plane bin %d\n",iObsBin,iEPBin);
      CalculateResultsObsBinEPBin(iObsBin,iEPBin,cCalcQACanvas);
    }
  }

  CalculateRPFErrorYieldsRms();

  CreateRatioAndDifferencesGraphs();

  // Calculate Ratios
  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    printf("Preparing to calculate ratios for obs bin %d.\n",iObsBin);
    CalculateResultsRatiosObsBin(iObsBin,cCalcQACanvas);
  }

  CalculateRPFErrorRatios();

  // Calculate other Axes Ratios?



  // Calculate Differences
  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    printf("Preparing to calculate differences for obs bin %d.\n",iObsBin);
    CalculateResultsDifferencesObsBin(iObsBin,cCalcQACanvas);
  }

}

/**
  * Remove the first nSkipPoints points
  */
void TaskCalcObservables::CleanResults() {

  vector<TGraphErrors * > gResultsGraphs = {
    OutOverIn_AS, OutOverIn_NS, MidOverIn_AS, MidOverIn_NS,
    OutOverIn_AS_RPFError, OutOverIn_NS_RPFError, MidOverIn_AS_RPFError, MidOverIn_NS_RPFError,
    RMSOutOverIn_AS, RMSOutOverIn_NS, RMSMidOverIn_AS, RMSMidOverIn_NS,
    OutMinusIn_AS, OutMinusIn_NS, MidMinusIn_AS, MidMinusIn_NS
  };

  for (TGraphErrors * graph : gResultsGraphs) {
    for (int i = 0; i < nSkipPoints; i++) {
      graph->RemovePoint(0); // After removing the first point, the second point will be the first point
    }
  }
}


void TaskCalcObservables::CalculateResultsObsBinEPBin(int iObsBin, int iEPBin, TCanvas * canv) {
  canv->Clear();

  double xValue = (fObsBins[iObsBin+1] + fObsBins[iObsBin]) / 2.;
  double xWidth = (fObsBins[iObsBin+1] - fObsBins[iObsBin]) / 2.;


  // For NS, use Near Eta
  // For AS, use Full
  // If more than the Far eta nearside region is used for the background, this should be revisited.

  TH1D * histNearSide = fNearEtaDPhiProj_Sub[iObsBin][iEPBin];
  TH1D * histAwaySide = fFullDPhiProj_Sub[iObsBin][iEPBin];

  // RPF Variations
  TH2F * histNearSideVar = fNearEtaDPhiProj_Sub_RPFVar[iObsBin][iEPBin];
  TH2F * histAwaySideVar = fFullDPhiProj_Sub_RPFVar[iObsBin][iEPBin];


  printf("Calculating Near Side Results with histogram %s\n",histNearSide->GetName());

  Double_t fYieldNS     = 0;
  Double_t fYieldNS_Err = 0;

  Double_t fRmsNS     = 0;
  Double_t fRmsNS_Err = 0;

//  histNearSide->GetXaxis()->SetRangeUser(-fYieldRangeNS,fYieldRangeNS);
 
  Int_t iMinRangeBin = histNearSide->FindBin(-fYieldRangeNS);
  Int_t iMaxRangeBin = histNearSide->FindBin(fYieldRangeNS);

  fYieldNS = histNearSide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldNS_Err,"width");

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fYieldNS_Var = 0;
    Double_t fYieldNS_Var_Err = 0;

    // Integrate over the Y-axis bin iVar
    fYieldNS_Var = histNearSideVar->IntegralAndError(iMinRangeBin,iMaxRangeBin,iVar+1,iVar+1,fYieldNS_Var_Err);

    if (iEPBin == kNEPBins) {
      fNSYieldsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fYieldNS_Var);
      fNSYieldsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fYieldNS_Var_Err);
    } else {
      fNSYieldsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fYieldNS_Var);
      fNSYieldsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fYieldNS_Var_Err);
    }
  }

  histNearSide->SetMarkerColor(kNonSelectColor);
  histNearSide->SetLineColor(kNonSelectColor);
  histNearSide->DrawCopy();

  histNearSide->GetXaxis()->SetRangeUser(-fRmsRangeNS,fRmsRangeNS);
  fRmsNS     = histNearSide->GetRMS();
  fRmsNS_Err = histNearSide->GetRMSError();

  Int_t iMinRangeRMS = histNearSide->GetXaxis()->FindBin(-fRmsRangeNS);
  Int_t iMaxRangeRMS = histNearSide->GetXaxis()->FindBin(fRmsRangeNS);

  TProfile * histNearSideRmsProfile = 0;
  if (nRPFVariants > 0) {
    histNearSideRmsProfile = histNearSideVar->ProfileY("_pfy",iMinRangeRMS,iMaxRangeRMS,"s");
  }

  // Draw a clone with a different color after setting the range?
  histNearSide->SetMarkerColor(kSelectColor);
  histNearSide->SetLineColor(kSelectColor);
  histNearSide->Draw("SAME");

  printf("  Found Near-side Yield = %e \\pm %e and RMS = %f \\pm %f\n",fYieldNS,fYieldNS_Err,fRmsNS,fRmsNS_Err);


  if (iEPBin == kNEPBins) {
    fNSYieldsInc->SetPoint(iObsBin,xValue,fYieldNS);
    fNSYieldsInc->SetPointError(iObsBin,xWidth,fYieldNS_Err);

    fNSRmsInc->SetPoint(iObsBin,xValue,fRmsNS);
    fNSRmsInc->SetPointError(iObsBin,xWidth,fRmsNS_Err);
  } else {
    fNSYieldsEP[iEPBin]->SetPoint(iObsBin,xValue,fYieldNS);
    fNSYieldsEP[iEPBin]->SetPointError(iObsBin,xWidth,fYieldNS_Err);

    fNSRmsEP[iEPBin]->SetPoint(iObsBin,xValue,fRmsNS);
    fNSRmsEP[iEPBin]->SetPointError(iObsBin,xWidth,fRmsNS_Err);
  }


  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fRmsNS_Var = 0;
    Double_t fRmsNS_Var_Err = 0;
    fRmsNS_Var = histNearSideRmsProfile->GetBinContent(iVar+1);
    fRmsNS_Var_Err = histNearSideRmsProfile->GetBinError(iVar+1);

    if (iEPBin == kNEPBins) {
      fNSRmsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fRmsNS_Var);
      fNSRmsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fRmsNS_Var_Err);
    } else {
      fNSRmsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fRmsNS_Var);
      fNSRmsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fRmsNS_Var_Err);
    }
  }

  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_NearSide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_NearSide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_NearSide.C",fOutputDir.Data(),iObsBin,iEPBin));


  printf("Calculating Away Side Results with histogram %s\n",histAwaySide->GetName());

  TProfile * histAwaySideRmsProfile = 0;
  if (nRPFVariants > 0) {
    histAwaySideRmsProfile = histAwaySideVar->ProfileY("_pfy",iMinRangeRMS,iMaxRangeRMS,"s");
  }

  Double_t fYieldAS     = 0;
  Double_t fYieldAS_Err = 0;

  Double_t fRmsAS     = 0;
  Double_t fRmsAS_Err = 0;

  iMinRangeBin = histAwaySide->FindBin(TMath::Pi()-fYieldRangeAS);
  iMaxRangeBin = histAwaySide->FindBin(TMath::Pi()+fYieldRangeAS);

  fYieldAS = histAwaySide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldAS_Err,"width");

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fYieldAS_Var = 0;
    Double_t fYieldAS_Var_Err = 0;
    Double_t fRmsAS_Var = 0;
    Double_t fRmsAS_Var_Err = 0;

    // Integrate over the Y-axis bin iVar
    fYieldAS_Var = histAwaySideVar->IntegralAndError(iMinRangeBin,iMaxRangeBin,iVar+1,iVar+1,fYieldAS_Var_Err);
    // probably have to project a histogram to get the RMS

    if (iEPBin == kNEPBins) {
      fASYieldsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fYieldAS_Var);
      fASYieldsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fYieldAS_Var_Err);
    } else {
      fASYieldsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fYieldAS_Var);
      fASYieldsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fYieldAS_Var_Err);
    }
  }


  histAwaySide->SetMarkerColor(kNonSelectColor);
  histAwaySide->SetLineColor(kNonSelectColor);
  histAwaySide->DrawCopy();

  histAwaySide->GetXaxis()->SetRangeUser(TMath::Pi()-fRmsRangeAS,TMath::Pi()+fRmsRangeAS);
  fRmsAS     = histAwaySide->GetRMS();
  fRmsAS_Err = histAwaySide->GetRMSError();

  histAwaySide->SetMarkerColor(kSelectColor);
  histAwaySide->SetFillColor(kSelectColor);
  histAwaySide->SetLineColor(kSelectColor);
  histAwaySide->Draw("SAME F");

  printf("  Found Away-side Yield = %e \\pm %e and RMS = %f \\pm %f\n",fYieldNS,fYieldNS_Err,fRmsNS,fRmsNS_Err);

  if (iEPBin == kNEPBins) {
    fASYieldsInc->SetPoint(iObsBin,xValue,fYieldAS);
    fASYieldsInc->SetPointError(iObsBin,xWidth,fYieldAS_Err);

    fASRmsInc->SetPoint(iObsBin,xValue,fRmsAS);
    fASRmsInc->SetPointError(iObsBin,xWidth,fRmsAS_Err);
  } else {
    fASYieldsEP[iEPBin]->SetPoint(iObsBin,xValue,fYieldAS);
    fASYieldsEP[iEPBin]->SetPointError(iObsBin,xWidth,fYieldAS_Err);

    fASRmsEP[iEPBin]->SetPoint(iObsBin,xValue,fRmsAS);
    fASRmsEP[iEPBin]->SetPointError(iObsBin,xWidth,fRmsAS_Err);
  }




  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fRmsAS_Var = 0;
    Double_t fRmsAS_Var_Err = 0;
    fRmsAS_Var = histAwaySideRmsProfile->GetBinContent(iVar+1);
    fRmsAS_Var_Err = histAwaySideRmsProfile->GetBinError(iVar+1);

    if (iEPBin == kNEPBins) {
      fASRmsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fRmsAS_Var);
      fASRmsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fRmsAS_Var_Err);
    } else {
      fASRmsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fRmsAS_Var);
      fASRmsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fRmsAS_Var_Err);
    }
  }

  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_AwaySide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_AwaySide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_AwaySide.C",fOutputDir.Data(),iObsBin,iEPBin));


  // Now study delta eta projection of nearside peak




}

void TaskCalcObservables::CreateRatioAndDifferencesGraphs() {

  OutOverIn_AS = new TGraphErrors(nObsBins);
  OutOverIn_NS = new TGraphErrors(nObsBins);
  MidOverIn_AS = new TGraphErrors(nObsBins);
  MidOverIn_NS = new TGraphErrors(nObsBins);

  TString sSideName[2]  = {"AS","NS"};
  TString sSideTitle[2] = {"Away-Side","Near-Side"};
  TString sNumerator[2] = {"Out","Mid"};

  int iNumerator = 0; // 0 for Out
  int iSide = 0; // 0 for AwaySide

  OutOverIn_AS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutOverIn_AS->SetLineColor(kOutInColor);
  OutOverIn_AS->SetMarkerColor(kOutInColor);
  OutOverIn_AS->SetMarkerStyle(kOutOverInMarker);
  OutOverIn_AS_RPFError = (TGraphErrors *) OutOverIn_AS->Clone(Form("%s_RPFError",OutOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) OutOverIn_AS->Clone(Form("%s_RPFVar%d",OutOverIn_AS->GetName(),iVar));
    OutOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  MidOverIn_AS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  MidOverIn_AS->SetLineColor(kMidInColor);
  MidOverIn_AS->SetMarkerColor(kMidInColor);
  MidOverIn_AS->SetMarkerStyle(kMidOverInMarker);
  MidOverIn_AS_RPFError = (TGraphErrors *) MidOverIn_AS->Clone(Form("%s_RPFError",MidOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) MidOverIn_AS->Clone(Form("%s_RPFVar%d",MidOverIn_AS->GetName(),iVar));
    MidOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  OutOverIn_NS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutOverIn_NS->SetLineColor(kOutInColor);
  OutOverIn_NS->SetMarkerColor(kOutInColor);
  OutOverIn_NS->SetMarkerStyle(kOutOverInMarker);
  OutOverIn_NS_RPFError = (TGraphErrors *) OutOverIn_NS->Clone(Form("%s_RPFError",OutOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) OutOverIn_NS->Clone(Form("%s_RPFVar%d",OutOverIn_NS->GetName(),iVar));
    OutOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  MidOverIn_NS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  MidOverIn_NS->SetLineColor(kMidInColor);
  MidOverIn_NS->SetMarkerColor(kMidInColor);
  MidOverIn_NS->SetMarkerStyle(kMidOverInMarker);
  MidOverIn_NS_RPFError = (TGraphErrors *) MidOverIn_NS->Clone(Form("%s_RPFError",MidOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) MidOverIn_NS->Clone(Form("%s_RPFVar%d",MidOverIn_NS->GetName(),iVar));
    MidOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  RMSOutOverIn_AS = new TGraphErrors(nObsBins);
  RMSOutOverIn_NS = new TGraphErrors(nObsBins);
  RMSMidOverIn_AS = new TGraphErrors(nObsBins);
  RMSMidOverIn_NS = new TGraphErrors(nObsBins);


  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  RMSOutOverIn_AS->SetLineColor(kOutInColor);
  RMSOutOverIn_AS->SetMarkerColor(kOutInColor);
  RMSOutOverIn_AS->SetMarkerStyle(kOutOverInMarker);
  RMSOutOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSOutOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsOutOverIn_AS_RPFError = (TGraphErrors *) RMSOutOverIn_AS->Clone(Form("%s_RPFError",RMSOutOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RMSOutOverIn_AS->Clone(Form("%s_RPFVar%d",RMSOutOverIn_AS->GetName(),iVar));
    RmsOutOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  RMSMidOverIn_AS->SetLineColor(kMidInColor);
  RMSMidOverIn_AS->SetMarkerColor(kMidInColor);
  RMSMidOverIn_AS->SetMarkerStyle(kMidOverInMarker);
  RMSMidOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSMidOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsMidOverIn_AS_RPFError = (TGraphErrors *) RMSMidOverIn_AS->Clone(Form("%s_RPFError",RMSMidOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RMSMidOverIn_AS->Clone(Form("%s_RPFVar%d",RMSMidOverIn_AS->GetName(),iVar));
    RmsMidOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  RMSOutOverIn_NS->SetLineColor(kOutInColor);
  RMSOutOverIn_NS->SetMarkerColor(kOutInColor);
  RMSOutOverIn_NS->SetMarkerStyle(kOutOverInMarker);
  RMSOutOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSOutOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsOutOverIn_NS_RPFError = (TGraphErrors *) RMSOutOverIn_NS->Clone(Form("%s_RPFError",RMSOutOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RMSOutOverIn_NS->Clone(Form("%s_RPFVar%d",RMSOutOverIn_NS->GetName(),iVar));
    RmsOutOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  RMSMidOverIn_NS->SetLineColor(kMidInColor);
  RMSMidOverIn_NS->SetMarkerColor(kMidInColor);
  RMSMidOverIn_NS->SetMarkerStyle(kMidOverInMarker);
  RMSMidOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSMidOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsMidOverIn_NS_RPFError = (TGraphErrors *) RMSMidOverIn_NS->Clone(Form("%s_RPFError",RMSMidOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RMSMidOverIn_NS->Clone(Form("%s_RPFVar%d",RMSMidOverIn_NS->GetName(),iVar));
    RmsMidOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  OutMinusIn_AS = new TGraphErrors(nObsBins);
  OutMinusIn_NS = new TGraphErrors(nObsBins);
  MidMinusIn_AS = new TGraphErrors(nObsBins);
  MidMinusIn_NS = new TGraphErrors(nObsBins);


  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  OutMinusIn_AS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutMinusIn_AS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iNumerator = 1; // 1 for Mid
  MidMinusIn_AS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidMinusIn_AS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  OutMinusIn_NS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutMinusIn_NS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iNumerator = 1; // 1 for Mid
  MidMinusIn_NS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidMinusIn_NS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

}



// Calculate Ratios
void TaskCalcObservables::CalculateResultsRatiosObsBin(int iObsBin, TCanvas * canv) {


  double xValue = fASYieldsEP[0]->GetX()[iObsBin];
  double xErr   = fASYieldsEP[0]->GetEX()[iObsBin];

  double YieldInPlane      = fASYieldsEP[0]->GetY()[iObsBin];
  double YieldInPlaneError = fASYieldsEP[0]->GetEY()[iObsBin];

  double YieldMidPlane      = fASYieldsEP[1]->GetY()[iObsBin];
  double YieldMidPlaneError = fASYieldsEP[1]->GetEY()[iObsBin];

  double YieldOutPlane      = fASYieldsEP[2]->GetY()[iObsBin];
  double YieldOutPlaneError = fASYieldsEP[2]->GetEY()[iObsBin];

  if (YieldInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 yield in-plane\n");
    return;  
  }
  double YieldOutOverIn = YieldOutPlane / YieldInPlane;
  double YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  double YieldMidOverIn = YieldMidPlane / YieldInPlane;
  double YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  OutOverIn_AS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  OutOverIn_AS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  MidOverIn_AS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  MidOverIn_AS->SetPointError(iObsBin,xErr,YieldMidOverInError);

  
  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {

    YieldInPlane      = fASYieldsEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    YieldInPlaneError = fASYieldsEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    YieldMidPlane      = fASYieldsEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    YieldMidPlaneError = fASYieldsEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    YieldOutPlane      = fASYieldsEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    YieldOutPlaneError = fASYieldsEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (YieldInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 yield in-plane\n");
      continue;  
    }

    YieldOutOverIn = YieldOutPlane / YieldInPlane;
    YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    OutOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldOutOverIn);
    OutOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldOutOverInError);

    MidOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldMidOverIn);
    MidOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldMidOverInError);

  }

  // Repeat for Near-Side

  YieldInPlane      = fNSYieldsEP[0]->GetY()[iObsBin];
  YieldInPlaneError = fNSYieldsEP[0]->GetEY()[iObsBin];

  YieldMidPlane      = fNSYieldsEP[1]->GetY()[iObsBin];
  YieldMidPlaneError = fNSYieldsEP[1]->GetEY()[iObsBin];

  YieldOutPlane      = fNSYieldsEP[2]->GetY()[iObsBin];
  YieldOutPlaneError = fNSYieldsEP[2]->GetEY()[iObsBin];

  if (YieldInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 yield in-plane\n");
    return;  
  }

  YieldOutOverIn = YieldOutPlane / YieldInPlane;
  YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  YieldMidOverIn = YieldMidPlane / YieldInPlane;
  YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  OutOverIn_NS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  OutOverIn_NS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  MidOverIn_NS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  MidOverIn_NS->SetPointError(iObsBin,xErr,YieldMidOverInError);



  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {

    YieldInPlane      = fNSYieldsEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    YieldInPlaneError = fNSYieldsEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    YieldMidPlane      = fNSYieldsEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    YieldMidPlaneError = fNSYieldsEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    YieldOutPlane      = fNSYieldsEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    YieldOutPlaneError = fNSYieldsEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (YieldInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 yield in-plane\n");
      continue;  
    }

    YieldOutOverIn = YieldOutPlane / YieldInPlane;
    YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    OutOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldOutOverIn);
    OutOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldOutOverInError);

    MidOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldMidOverIn);
    MidOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldMidOverInError);

  }



  // RMS 
  // Away Side

  YieldInPlane      = fASRmsEP[0]->GetY()[iObsBin];
  YieldInPlaneError = fASRmsEP[0]->GetEY()[iObsBin];

  YieldMidPlane      = fASRmsEP[1]->GetY()[iObsBin];
  YieldMidPlaneError = fASRmsEP[1]->GetEY()[iObsBin];

  YieldOutPlane      = fASRmsEP[2]->GetY()[iObsBin];
  YieldOutPlaneError = fASRmsEP[2]->GetEY()[iObsBin];

  if (YieldInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 width in-plane\n");
    return;  
  }

  YieldOutOverIn = YieldOutPlane / YieldInPlane;
  YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  YieldMidOverIn = YieldMidPlane / YieldInPlane;
  YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  RMSOutOverIn_AS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  RMSOutOverIn_AS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  RMSMidOverIn_AS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  RMSMidOverIn_AS->SetPointError(iObsBin,xErr,YieldMidOverInError);

  // Calculation with RPF variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    YieldInPlane      = fASRmsEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    YieldInPlaneError = fASRmsEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    YieldMidPlane      = fASRmsEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    YieldMidPlaneError = fASRmsEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    YieldOutPlane      = fASRmsEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    YieldOutPlaneError = fASRmsEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (YieldInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 AS Rms in-plane\n");
      continue;  
    }
    YieldOutOverIn = YieldOutPlane / YieldInPlane;
    YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));
    
    RmsOutOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldOutOverIn);
    RmsOutOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldOutOverInError);

    RmsMidOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldMidOverIn);
    RmsMidOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldMidOverInError);

  }

  // Near Side

  YieldInPlane      = fNSRmsEP[0]->GetY()[iObsBin];
  YieldInPlaneError = fNSRmsEP[0]->GetEY()[iObsBin];

  YieldMidPlane      = fNSRmsEP[1]->GetY()[iObsBin];
  YieldMidPlaneError = fNSRmsEP[1]->GetEY()[iObsBin];

  YieldOutPlane      = fNSRmsEP[2]->GetY()[iObsBin];
  YieldOutPlaneError = fNSRmsEP[2]->GetEY()[iObsBin];

  if (YieldInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 width in-plane\n");
    return;  
  }

  YieldOutOverIn = YieldOutPlane / YieldInPlane;
  YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  YieldMidOverIn = YieldMidPlane / YieldInPlane;
  YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));


  RMSOutOverIn_NS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  RMSOutOverIn_NS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  RMSMidOverIn_NS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  RMSMidOverIn_NS->SetPointError(iObsBin,xErr,YieldMidOverInError);

  // Calculation with RPF variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    YieldInPlane      = fNSRmsEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    YieldInPlaneError = fNSRmsEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    YieldMidPlane      = fNSRmsEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    YieldMidPlaneError = fNSRmsEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    YieldOutPlane      = fNSRmsEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    YieldOutPlaneError = fNSRmsEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (YieldInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 NS Rms in-plane\n");
      continue;  
    }
    YieldOutOverIn = YieldOutPlane / YieldInPlane;
    YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));
    
    RmsOutOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldOutOverIn);
    RmsOutOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldOutOverInError);

    RmsMidOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldMidOverIn);
    RmsMidOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldMidOverInError);

  }


//  TGraphErrors * RMSOutOverIn_AS;
//  TGraphErrors * RMSOutOverIn_NS;
//  TGraphErrors * RMSMidOverIn_AS;
//  TGraphErrors * RMSMidOverIn_NS;


}

/**
  * Calculate the errors in the yields and Rms
  */
void TaskCalcObservables::CalculateRPFErrorYieldsRms() {

  //vector<vector<TH2F *>> fFullDPhiProj_Sub_RPFVar

// TGraphErrors * fNSYieldsInc; ///< Near-Side Yields in all EP
//  vector<TGraphErrors *> fNSYieldsEP; ///< Near-Side Yields in the EP bins

  fNSYieldsInc_RPFError = (TGraphErrors *) fNSYieldsInc->Clone(Form("%s_RPF_Error",fNSYieldsInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSYieldsEP_RPFError.push_back((TGraphErrors *) fNSYieldsEP[iEPBin]->Clone(Form("%s_RPF_Error",fNSYieldsEP[iEPBin]->GetName())));
  }

  fASYieldsInc_RPFError = (TGraphErrors *) fASYieldsInc->Clone(Form("%s_RPF_Error",fASYieldsInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASYieldsEP_RPFError.push_back((TGraphErrors *) fASYieldsEP[iEPBin]->Clone(Form("%s_RPF_Error",fASYieldsEP[iEPBin]->GetName())));
  }

  fNSRmsInc_RPFError = (TGraphErrors *) fNSRmsInc->Clone(Form("%s_RPF_Error",fNSRmsInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSRmsEP_RPFError.push_back((TGraphErrors *) fNSRmsEP[iEPBin]->Clone(Form("%s_RPF_Error",fNSRmsEP[iEPBin]->GetName())));
  }

  fASRmsInc_RPFError = (TGraphErrors *) fASRmsInc->Clone(Form("%s_RPF_Error",fASRmsInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASRmsEP_RPFError.push_back((TGraphErrors *) fASRmsEP[iEPBin]->Clone(Form("%s_RPF_Error",fASRmsEP[iEPBin]->GetName())));
  }


  TGraphErrors * fNSYieldsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fNSYieldsEP_RPFErr_Mean = {};
  TGraphErrors * fASYieldsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fASYieldsEP_RPFErr_Mean = {};

  TGraphErrors * fNSRmsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fNSRmsEP_RPFErr_Mean = {};
  TGraphErrors * fASRmsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fASRmsEP_RPFErr_Mean = {};


  printf("Calculating RPF Systematic for nearside yields\n");

  // Use these:
  // fNSYieldsInc_RPFVariants
  // fNSYieldsEP_RPFVariants
  // fASYieldsInc_RPFVariants
  // fASYieldsEP_RPFVariants

  fNSYieldsInc_RPFErr_Mean = ProduceSystematicFromGraphs(fNSYieldsInc_RPFVariants,fNSYieldsInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSYieldsEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fNSYieldsEP_RPFVariants[iEPBin],fNSYieldsEP_RPFError[iEPBin]));
  }
  fASYieldsInc_RPFErr_Mean = ProduceSystematicFromGraphs(fASYieldsInc_RPFVariants,fASYieldsInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASYieldsEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fASYieldsEP_RPFVariants[iEPBin],fASYieldsEP_RPFError[iEPBin]));
  }

  fNSRmsInc_RPFErr_Mean = ProduceSystematicFromGraphs(fNSRmsInc_RPFVariants,fNSRmsInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSRmsEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fNSRmsEP_RPFVariants[iEPBin],fNSRmsEP_RPFError[iEPBin]));
  }
  fASRmsInc_RPFErr_Mean = ProduceSystematicFromGraphs(fASRmsInc_RPFVariants,fASRmsInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASRmsEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fASRmsEP_RPFVariants[iEPBin],fASRmsEP_RPFError[iEPBin]));
  }

  printf("Finished calculating RPF fit error for first calculated variables.\n");

  fNSYieldsInc_RPFErr_Mean->SetLineColor(kBlack);
  


}

/**
  * Calculate the errors in the ratios from the RPF variations
  */
void TaskCalcObservables::CalculateRPFErrorRatios() {

// goal: use 
//   vector<TGraphErrors *> OutOverIn_AS_RPFVariants;
// to fill
//   TGraphErrors * OutOverIn_AS_RPFError

  // FIXME these will return the average value as the points. If we so choose, we could use this as the "final" result, or just overwrite the points with the central value already calculated
  // or maybe pass in the central result


  TGraphErrors * OutOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * MidOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * OutOverIn_NS_RPFErr_Mean = 0;
  TGraphErrors * MidOverIn_NS_RPFErr_Mean = 0;

  for (int i = 0; i < OutOverIn_AS->GetN(); i++) {
    OutOverIn_AS_RPFError->SetPoint(i,OutOverIn_AS->GetX()[i],OutOverIn_AS->GetY()[i]);
    MidOverIn_AS_RPFError->SetPoint(i,MidOverIn_AS->GetX()[i],MidOverIn_AS->GetY()[i]);
  }
  for (int i = 0; i < OutOverIn_NS->GetN(); i++) {
    OutOverIn_NS_RPFError->SetPoint(i,OutOverIn_NS->GetX()[i],OutOverIn_NS->GetY()[i]);
    MidOverIn_NS_RPFError->SetPoint(i,MidOverIn_NS->GetX()[i],MidOverIn_NS->GetY()[i]);
  }

  printf("Calculating RPF Systematic for OutOverIn_AS\n");
  OutOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(OutOverIn_AS_RPFVariants,OutOverIn_AS_RPFError);
  printf("Calculating RPF Systematic for MidOverIn_AS\n");
  MidOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(MidOverIn_AS_RPFVariants,MidOverIn_AS_RPFError);

  printf("Calculating RPF Systematic for OutOverIn_NS\n");
  OutOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(OutOverIn_NS_RPFVariants,OutOverIn_NS_RPFError);
  printf("Calculating RPF Systematic for MidOverIn_NS\n");
  MidOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(MidOverIn_NS_RPFVariants,MidOverIn_NS_RPFError);


  // Red / Blue, Green / Blue
 /* int kOutInErrorColor = kOutInColor;
  int kMidInErrorColor = kMidInColor;

  int kOutInFillStyle = 3245;
  int kMidInFillStyle = 3254;
  int kOutInNSFillStyle = 3245;
  int kMidInNSFillStyle = 3254;
*/

  // Settings for systematic uncertainty. E4
  OutOverIn_AS_RPFError->SetFillColor(kOutInErrorColor);
  OutOverIn_AS_RPFError->SetFillStyle(kOutInFillStyle);
  OutOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
//  OutOverIn_AS_RPFError->SetDrawOption("E4");
  
  MidOverIn_AS_RPFError->SetFillColor(kMidInErrorColor);
  MidOverIn_AS_RPFError->SetFillStyle(kMidInFillStyle);
  MidOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
//  MidOverIn_AS_RPFError->SetDrawOption("E4");


  OutOverIn_NS_RPFError->SetFillColor(kOutInErrorColor);
  OutOverIn_NS_RPFError->SetFillStyle(kOutInNSFillStyle);
  OutOverIn_NS_RPFError->SetMarkerColor(kBlue-9);
//  OutOverIn_AS_RPFError->SetDrawOption("E4");
  
  MidOverIn_NS_RPFError->SetFillColor(kMidInErrorColor);
  MidOverIn_NS_RPFError->SetFillStyle(kMidInNSFillStyle);
  MidOverIn_NS_RPFError->SetMarkerColor(kBlue-9);



  TGraphErrors * RmsOutOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * RmsMidOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * RmsOutOverIn_NS_RPFErr_Mean = 0;
  TGraphErrors * RmsMidOverIn_NS_RPFErr_Mean = 0;

  for (int i = 0; i < RMSOutOverIn_AS->GetN(); i++) {
    RmsOutOverIn_AS_RPFError->SetPoint(i,RMSOutOverIn_AS->GetX()[i],RMSOutOverIn_AS->GetY()[i]);
    RmsMidOverIn_AS_RPFError->SetPoint(i,RMSMidOverIn_AS->GetX()[i],RMSMidOverIn_AS->GetY()[i]);
  }
  for (int i = 0; i < RMSOutOverIn_NS->GetN(); i++) {
    RmsOutOverIn_NS_RPFError->SetPoint(i,RMSOutOverIn_NS->GetX()[i],RMSOutOverIn_NS->GetY()[i]);
    RmsMidOverIn_NS_RPFError->SetPoint(i,RMSMidOverIn_NS->GetX()[i],RMSMidOverIn_NS->GetY()[i]);
  }

  printf("Calculating RPF Systematic for RMS OutOverIn_AS\n");
  RmsOutOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsOutOverIn_AS_RPFVariants,RmsOutOverIn_AS_RPFError);
  printf("Calculating RPF Systematic for RMS MidOverIn_AS\n");
  RmsMidOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsMidOverIn_AS_RPFVariants,RmsMidOverIn_AS_RPFError);

  printf("Calculating RPF Systematic for RMS OutOverIn_NS\n");
  RmsOutOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsOutOverIn_NS_RPFVariants,RmsOutOverIn_NS_RPFError);
  printf("Calculating RPF Systematic for RMS MidOverIn_NS\n");
  RmsMidOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsMidOverIn_NS_RPFVariants,RmsMidOverIn_NS_RPFError);


  RmsOutOverIn_AS_RPFError->SetFillColor(kOutInErrorColor);
  RmsOutOverIn_AS_RPFError->SetFillStyle(kOutInFillStyle);
  RmsOutOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
  
  RmsMidOverIn_AS_RPFError->SetFillColor(kMidInErrorColor);
  RmsMidOverIn_AS_RPFError->SetFillStyle(kMidInFillStyle);
  RmsMidOverIn_AS_RPFError->SetMarkerColor(kBlue-9);

  RmsOutOverIn_NS_RPFError->SetFillColor(kOutInErrorColor);
  RmsOutOverIn_NS_RPFError->SetFillStyle(kOutInNSFillStyle);
  RmsOutOverIn_NS_RPFError->SetMarkerColor(kBlue-9);
  
  RmsMidOverIn_NS_RPFError->SetFillColor(kMidInErrorColor);
  RmsMidOverIn_NS_RPFError->SetFillStyle(kMidInNSFillStyle);
  RmsMidOverIn_NS_RPFError->SetMarkerColor(kBlue-9);


}


// Calculate Differences
void TaskCalcObservables::CalculateResultsDifferencesObsBin(int iObsBin, TCanvas * canv) {

//  TGraphErrors * OutMinusIn_AS;
//  TGraphErrors * OutMinusIn_NS;
//  TGraphErrors * MidMinusIn_AS;
//  TGraphErrors * MidMinusIn_NS;


}


void TaskCalcObservables::SetGraphStyles() {

  // All EP together
  fNSYieldsInc->SetLineColor(kEPColorList[3]);
  fNSYieldsInc->SetMarkerColor(kEPColorList[3]);
  fNSYieldsInc->SetMarkerStyle(kEPMarkerList[3]);
  fNSRmsInc->SetLineColor(kEPColorList[3]);
  fNSRmsInc->SetMarkerColor(kEPColorList[3]);
  fNSRmsInc->SetMarkerStyle(kEPMarkerList[3]);
  fASYieldsInc->SetLineColor(kEPColorList[3]);
  fASYieldsInc->SetMarkerColor(kEPColorList[3]);
  fASYieldsInc->SetMarkerStyle(kEPMarkerList[3]);
  fASRmsInc->SetLineColor(kEPColorList[3]);
  fASRmsInc->SetMarkerColor(kEPColorList[3]);
  fASRmsInc->SetMarkerStyle(kEPMarkerList[3]);
  // By EP Bin
  for (int i = 0; i < kNEPBins; i++) {
    fNSYieldsEP[i]->SetLineColor(kEPColorList[i]);
    fNSYieldsEP[i]->SetMarkerColor(kEPColorList[i]);
    fNSYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fNSRmsEP[i]->SetLineColor(kEPColorList[i]);
    fNSRmsEP[i]->SetMarkerColor(kEPColorList[i]);
    fNSRmsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fASYieldsEP[i]->SetLineColor(kEPColorList[i]);
    fASYieldsEP[i]->SetMarkerColor(kEPColorList[i]);
    fASYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fASRmsEP[i]->SetLineColor(kEPColorList[i]);
    fASRmsEP[i]->SetMarkerColor(kEPColorList[i]);
    fASRmsEP[i]->SetMarkerStyle(kEPMarkerList[i]);

  }

}


// take array of {inc,in,mid,out} ?
// and array of RPFuncert {inc,in,mid,out} 
// array of systematic errors.
// or should it be array of arrays of systematic errors( one per inc,in,etc and one per systematic error type (v3,purity,etc)
  // Make it all{in,mid,out,incl} to match up with arrays elsewhere
void TaskCalcObservables::DrawObservable(vector<TGraphErrors *> fObsGraphs, vector<TGraphErrors *> fObsRPFErrors, vector<TGraphErrors *> fObsSysErrors) {

  bool bHasSysErrors = (fObsSysErrors.size() > 0);

  TCanvas * cObservable = new TCanvas("cObservable","cObservable",fYieldCanvasWidth,fYieldCanvasHeight);

  if (fObsGraphs.size() != 4) {
    printf("DrawObservable Code not being used right\n");
    return;
  }

  int ObsType = 0; // 0 for yields, 1 for Rms
  TString sName = fObsGraphs[0]->GetName();
  if (sName.Contains("Rms")) ObsType = 1;

  if (ObsType == 0) cObservable->SetLogy(1);
  else cObservable->SetLogy(0);

  for (int i = 0; i < (int) fObsGraphs.size(); i++) {

    TString sTitleShort = fEPBinTitlesShort[i];

    TGraphErrors * fObsGraph = fObsGraphs[i];
    fObsGraph->SetMarkerStyle(kEPMarkerList[i]);
    fObsGraph->SetMarkerColor(kEPColorList[i]);
    fObsGraph->SetLineColor(kEPColorList[i]);

    TGraphErrors * fObsRPFError = fObsRPFErrors[i];
    fObsRPFError->SetMarkerStyle(kEPMarkerList[i]);
    fObsRPFError->SetMarkerColor(kEPColorList[i]);
    fObsRPFError->SetLineColor(kEPColorList[i]);
    fObsRPFError->SetFillColor(kEPRPFFillColorList[i]);
    fObsRPFError->SetFillStyle(kEPRPFFillStyleList[i]);

    TGraphErrors * fObsSysError = 0;
    if (bHasSysErrors) {
      printf("Intend to code something here to plot the other systematic errors\n");

    }


    fObsGraph->Draw("ALP");
    fObsRPFError->Draw("E3 SAME");





    cObservable->Print(Form("%s/%sWithErrors.pdf",fOutputDir.Data(),fObsGraph->GetName()));
    cObservable->Print(Form("%s/%sWithErrors.png",fOutputDir.Data(),fObsGraph->GetName()));

  }





}


void TaskCalcObservables::DrawResults() {
  printf("Saving some results plots\n");

  SetGraphStyles();

  float fRatioMin = 0.4;
  float fRatioMax = 1.6;


  // Drawing Individual Observable Plots
  vector<TGraphErrors *> fNSYieldsGraphs = {};
  vector<TGraphErrors *> fNSYieldsRPFErrors = {};
  vector<TGraphErrors *> fNSRmsGraphs = {};
  vector<TGraphErrors *> fNSRmsRPFErrors = {};
  vector<TGraphErrors *> fASYieldsGraphs = {};
  vector<TGraphErrors *> fASYieldsRPFErrors = {};
  vector<TGraphErrors *> fASRmsGraphs = {};
  vector<TGraphErrors *> fASRmsRPFErrors = {};

  for (int i = 0; i < 4; i++) {
    if (i == 3) {
      fNSYieldsGraphs.push_back(fNSYieldsInc);
      fNSYieldsRPFErrors.push_back(fNSYieldsInc_RPFError);   
      fNSRmsGraphs.push_back(fNSRmsInc);
      fNSRmsRPFErrors.push_back(fNSRmsInc_RPFError);   
      fASYieldsGraphs.push_back(fASYieldsInc);
      fASYieldsRPFErrors.push_back(fASYieldsInc_RPFError);   
      fASRmsGraphs.push_back(fASRmsInc);
      fASRmsRPFErrors.push_back(fASRmsInc_RPFError);   
    } else {
      fNSYieldsGraphs.push_back(fNSYieldsEP[i]);
      fNSYieldsRPFErrors.push_back(fNSYieldsEP_RPFError[i]);
      fNSRmsGraphs.push_back(fNSRmsEP[i]);
      fNSRmsRPFErrors.push_back(fNSRmsEP_RPFError[i]);
      fASYieldsGraphs.push_back(fASYieldsEP[i]);
      fASYieldsRPFErrors.push_back(fASYieldsEP_RPFError[i]);
      fASRmsGraphs.push_back(fASRmsEP[i]);
      fASRmsRPFErrors.push_back(fASRmsEP_RPFError[i]);
    }
  }

  DrawObservable(fNSYieldsGraphs,fNSYieldsRPFErrors);
  DrawObservable(fNSRmsGraphs,fNSRmsRPFErrors);
  DrawObservable(fASYieldsGraphs,fASYieldsRPFErrors);
  DrawObservable(fASRmsGraphs,fASRmsRPFErrors);



  // Undecided about whether to have the same canvas sizes for yields and widths apparently
  TCanvas * cResults = new TCanvas("cResults","cResults",fYieldCanvasWidth,fYieldCanvasHeight);

  TMultiGraph * mgNSY = new TMultiGraph();
  TMultiGraph * mgNSRms = new TMultiGraph();
  TMultiGraph * mgASY = new TMultiGraph();
  TMultiGraph * mgASRms = new TMultiGraph();

  TString sDrawStyle = "AP";
  TLegend * legYields = new TLegend(0.55,0.45,0.85,0.85);
  TLegend * legRms = new TLegend(0.55,0.65,0.85,0.85);




  cResults->cd();
  // Near-side Yields
  mgNSY->Add(fNSYieldsInc);
  legYields->AddEntry(fNSYieldsInc,"Inclusive","LP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legYields->AddEntry(fNSYieldsEP[iEPBin],fEPBinTitles[iEPBin],"LP");
    mgNSY->Add(fNSYieldsEP[iEPBin]);
  }
  mgNSY->Draw(sDrawStyle);
  //cResults->SetLogy(1);
  mgNSY->GetXaxis()->SetTitle(fNSYieldsInc->GetXaxis()->GetTitle());
  mgNSY->GetYaxis()->SetTitle(fNSYieldsInc->GetTitle());
  legYields->Draw("SAME");

  cResults->Print(Form("%s/NearsideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYields.C",fOutputDir.Data()));
  //cResults->SetLogy(0);

  // Near-side Widths
  mgNSRms->Add(fNSRmsInc);
  legRms->AddEntry(fNSRmsInc,"Inclusive","LP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legRms->AddEntry(fNSRmsEP[iEPBin],fEPBinTitles[iEPBin],"LP");
    mgNSRms->Add(fNSRmsEP[iEPBin]);
  }
  mgNSRms->Draw(sDrawStyle);
  mgNSRms->GetXaxis()->SetTitle(fNSRmsInc->GetXaxis()->GetTitle());
  mgNSRms->GetYaxis()->SetTitle(fNSRmsInc->GetTitle());
  legRms->Draw("SAME");
  cResults->Print(Form("%s/NearsideRms.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRms.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRms.C",fOutputDir.Data()));

  // Away-side Yields
  mgASY->Add(fASYieldsInc);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    mgASY->Add(fASYieldsEP[iEPBin]);
  }
  mgASY->Draw(sDrawStyle);
  mgASY->GetXaxis()->SetTitle(fASYieldsInc->GetXaxis()->GetTitle());
  mgASY->GetYaxis()->SetTitle(fASYieldsInc->GetTitle());
  legYields->Draw("SAME");

  //cResults->SetLogy(1);

  cResults->Print(Form("%s/AwaysideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYields.C",fOutputDir.Data()));

  //cResults->SetLogy(0);
  // Away-side Widths
  mgASRms->Add(fASRmsInc);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    mgASRms->Add(fASRmsEP[iEPBin]);
  }
  mgASRms->Draw(sDrawStyle);
  mgASRms->GetXaxis()->SetTitle(fASRmsInc->GetXaxis()->GetTitle());
  mgASRms->GetYaxis()->SetTitle(fASRmsInc->GetTitle());
  legRms->Draw("SAME");
  cResults->Print(Form("%s/AwaysideRms.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRms.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRms.C",fOutputDir.Data()));

  // =============================================================
  // Ratios
  // =============================================================

  // OutOverIn and MidOverIn AS
  TMultiGraph * mgASYieldRatios = new TMultiGraph();

  TLegend * legASYieldRatios = new TLegend(0.55,0.65,0.85,0.85);

  if (nRPFVariants > 0) {
    mgASYieldRatios->Add(OutOverIn_AS_RPFError);
    mgASYieldRatios->Add(MidOverIn_AS_RPFError);
  }

  mgASYieldRatios->Add(OutOverIn_AS);
  mgASYieldRatios->Add(MidOverIn_AS);

  legASYieldRatios->SetHeader("Awayside Yield Ratios","C");
  legASYieldRatios->SetBorderSize(0);
  legASYieldRatios->AddEntry(OutOverIn_AS,"Out/In (Stat. Error)","LP");
  legASYieldRatios->AddEntry(MidOverIn_AS,"Mid/In (Stat. Error)","LP");

  mgASYieldRatios->GetXaxis()->SetTitle(OutOverIn_AS->GetXaxis()->GetTitle());
  mgASYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
  mgASYieldRatios->Draw(sDrawStyle);

  if (nRPFVariants > 0) {
    printf("Drawing the error from RPF Variants\n");
    OutOverIn_AS_RPFError->Draw("E3 SAME");
    MidOverIn_AS_RPFError->Draw("E3 SAME");
 //   mgASYieldRatios->Add(OutOverIn_AS_RPFError);
 //   mgASYieldRatios->Add(MidOverIn_AS_RPFError);
    //mgASYieldRatios->Draw(Form("%s SAME",sDrawStyle.Data()));
   // mgASYieldRatios->Draw("SAME");

    legASYieldRatios->AddEntry(OutOverIn_AS_RPFError,"RPF Error (Out/In)","F");
    legASYieldRatios->AddEntry(MidOverIn_AS_RPFError,"RPF Error (Mid/In)","F");

  } else {
    printf("  Not adding RPF variation error because I don't have it\n");
  }
  //mgASYieldRatios->Draw("ALP");
  //mgASYieldRatios->Draw("P SAME");

  legASYieldRatios->Draw("SAME");
  mgASYieldRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/AwaysideYieldRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldRatios.C",fOutputDir.Data()));






  // OutOverIn and MidOverIn NS
  TMultiGraph * mgNSYieldRatios = new TMultiGraph();

  TLegend * legYieldRatios = new TLegend(0.55,0.65,0.85,0.85);

  mgNSYieldRatios->Add(OutOverIn_NS);
  mgNSYieldRatios->Add(MidOverIn_NS);

  legYieldRatios->SetBorderSize(0);
  legYieldRatios->SetHeader("Nearside Yield Ratios","C");
  legYieldRatios->AddEntry(OutOverIn_NS,"Out/In (Stat. Error)","LP");
  legYieldRatios->AddEntry(MidOverIn_NS,"Mid/In (Stat. Error)","LP");

  mgNSYieldRatios->GetXaxis()->SetTitle(OutOverIn_NS->GetXaxis()->GetTitle());
  mgNSYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
  mgNSYieldRatios->Draw(sDrawStyle);

  if (nRPFVariants > 0) {
    OutOverIn_NS_RPFError->Draw("E3 SAME");
    MidOverIn_NS_RPFError->Draw("E3 SAME");

    legYieldRatios->AddEntry(OutOverIn_NS_RPFError,"RPF Error (Out/In)","F");
    legYieldRatios->AddEntry(MidOverIn_NS_RPFError,"RPF Error (Mid/In)","F");
  }




  legYieldRatios->Draw("SAME");
  mgNSYieldRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/NearsideYieldRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYieldRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYieldRatios.C",fOutputDir.Data()));


  // Draw RMS Ratios
  TMultiGraph * mgNSRmsRatios = new TMultiGraph();
  TLegend * legNSRmsRatios = new TLegend(0.55,0.65,0.85,0.85);

  mgNSRmsRatios->Add(RMSOutOverIn_NS);
  mgNSRmsRatios->Add(RMSMidOverIn_NS);
  legNSRmsRatios->SetBorderSize(0);
  legNSRmsRatios->SetHeader("Nearside RMS Ratios","C");
  legNSRmsRatios->AddEntry(RMSOutOverIn_NS,"Out/In (Stat. Error)","LP");
  legNSRmsRatios->AddEntry(RMSMidOverIn_NS,"Mid/In (Stat. Error)","LP");

  mgNSRmsRatios->GetXaxis()->SetTitle(RMSOutOverIn_NS->GetXaxis()->GetTitle());
  mgNSRmsRatios->GetYaxis()->SetTitle("RMS Ratio");


  mgNSRmsRatios->Draw(sDrawStyle);

  if (nRPFVariants > 0) {
    RmsOutOverIn_NS_RPFError->Draw("E3 SAME");
    RmsMidOverIn_NS_RPFError->Draw("E3 SAME");

    legNSRmsRatios->AddEntry(RmsOutOverIn_NS_RPFError,"RPF Error (Out/In)","F");
    legNSRmsRatios->AddEntry(RmsMidOverIn_NS_RPFError,"RPF Error (Mid/In)","F");
  }


  legNSRmsRatios->Draw("SAME");
  mgNSRmsRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/NearsideRmsRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsRatios.C",fOutputDir.Data()));



  // Awayside
  TMultiGraph * mgASRmsRatios = new TMultiGraph();
  TLegend * legASRmsRatios = new TLegend(0.55,0.65,0.85,0.85);

  mgASRmsRatios->Add(RMSOutOverIn_AS);
  mgASRmsRatios->Add(RMSMidOverIn_AS);
  legASRmsRatios->SetBorderSize(0);
  legASRmsRatios->SetHeader("Awayside RMS Ratios","C");
  legASRmsRatios->AddEntry(RMSOutOverIn_AS,"Out/In (Stat. Error)","LP");
  legASRmsRatios->AddEntry(RMSMidOverIn_AS,"Mid/In (Stat. Error)","LP");

  mgASRmsRatios->GetXaxis()->SetTitle(RMSOutOverIn_AS->GetXaxis()->GetTitle());
  mgASRmsRatios->GetYaxis()->SetTitle("RMS Ratio");

  mgASRmsRatios->Draw(sDrawStyle);

  if (nRPFVariants > 0) {
    RmsOutOverIn_AS_RPFError->Draw("E3 SAME");
    RmsMidOverIn_AS_RPFError->Draw("E3 SAME");

    legASRmsRatios->AddEntry(RmsOutOverIn_AS_RPFError,"RPF Error (Out/In)","F");
    legASRmsRatios->AddEntry(RmsMidOverIn_AS_RPFError,"RPF Error (Mid/In)","F");
  }


  legASRmsRatios->Draw("SAME");
  mgASRmsRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/AwaysideRmsRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRmsRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRmsRatios.C",fOutputDir.Data()));



}

void TaskCalcObservables::SaveOutput() {

  if (!fOutputFile) {
    fprintf(stderr,"Error: No output file set!\n");
    return;
  }

  if (fASYieldsInc) fOutputFile->Add(fASYieldsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fASYieldsEP[i]) fOutputFile->Add(fASYieldsEP[i]);
  }
  if (fASRmsInc) fOutputFile->Add(fASRmsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fASRmsEP[i]) fOutputFile->Add(fASRmsEP[i]);
  }


  if (fNSYieldsInc) fOutputFile->Add(fNSYieldsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fNSYieldsEP[i]) fOutputFile->Add(fNSYieldsEP[i]);
  }
  if (fNSRmsInc) fOutputFile->Add(fNSRmsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fNSRmsEP[i]) fOutputFile->Add(fNSRmsEP[i]);
  }

  if (OutOverIn_AS) fOutputFile->Add(OutOverIn_AS);
  if (MidOverIn_AS) fOutputFile->Add(MidOverIn_AS);
  if (OutOverIn_NS) fOutputFile->Add(OutOverIn_NS);
  if (MidOverIn_NS) fOutputFile->Add(MidOverIn_NS);

  if (OutOverIn_AS_RPFError) fOutputFile->Add(OutOverIn_AS_RPFError);
  if (MidOverIn_AS_RPFError) fOutputFile->Add(MidOverIn_AS_RPFError);
  if (OutOverIn_NS_RPFError) fOutputFile->Add(OutOverIn_NS_RPFError);
  if (MidOverIn_NS_RPFError) fOutputFile->Add(MidOverIn_NS_RPFError);
  //if (OutOverIn_NS) fOutputFile->Add(OutOverIn_NS);
  //if (MidOverIn_NS) fOutputFile->Add(MidOverIn_NS);


  if (RMSOutOverIn_AS) fOutputFile->Add(RMSOutOverIn_AS);
  if (RMSMidOverIn_AS) fOutputFile->Add(RMSMidOverIn_AS);
  if (RMSOutOverIn_NS) fOutputFile->Add(RMSOutOverIn_NS);
  if (RMSMidOverIn_NS) fOutputFile->Add(RMSMidOverIn_NS);

  if (OutMinusIn_AS) fOutputFile->Add(OutMinusIn_AS);
  if (MidMinusIn_AS) fOutputFile->Add(MidMinusIn_AS);
  if (OutMinusIn_NS) fOutputFile->Add(OutMinusIn_NS);
  if (MidMinusIn_NS) fOutputFile->Add(MidMinusIn_NS);


  fOutputFile->Write();

}


void TaskCalcObservables::Run() {

  LoadHistograms();

  InitArrays();

  DrawDirectComparisons();

  CalculateResults();

  CleanResults();

  DrawResults();

  SaveOutput();

  cout<<"Done!"<<endl;
}

