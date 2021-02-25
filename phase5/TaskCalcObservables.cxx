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

  // Loading the Full Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_Full_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_Full_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);
    }
    fFullDPhiProj_Sub.push_back(fLocalVector);
  }
  // Loading the Near Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_NearEta_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_NearEta_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);
    }
    fNearEtaDPhiProj_Sub.push_back(fLocalVector);
  }
  // Loading the Far Eta Projections
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    for (Int_t j = 0; j < kNEPBins+1; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("dPhi_RPFMethod%d_FarEta_EP%d_RPFSub_ObsBin%d_Rescale",iRPFMethod,j,i);
      if (j == kNEPBins) fLocalName = Form("dPhi_RPFMethod%d_FarEta_AllEP_RPFSub_ObsBin%d",iRPFMethod,i);
      fLocal = (TH1D *) fInputFileCentral->Get(fLocalName);
      if (!fLocal) {
        fprintf(stderr,"Missing histogram %s\n",fLocalName.Data());
      }
      fLocalVector.push_back(fLocal);
    }
    fFarEtaDPhiProj_Sub.push_back(fLocalVector);
  }


}

void TaskCalcObservables::InitArrays() {

  Double_t fZtStep = 1.0/(7 - 1.0);
  Double_t fXiStep = 2.5/(8 - 1.0);

  Double_t * ObsArray = 0;

  Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,20,23,30,60};
  Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
  Double_t array_XI_BinsValue[kXiNBINS+1]   ={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};
  Double_t array_HPT_BinsValue[kNoHPtBins+1]={0.15,0.4,0.8,1.45,2.5,4.2,6.95,11.4,18.6};


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

  CreateRatioAndDifferencesGraphs();

  // Calculate Ratios
  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    printf("Preparing to calculate ratios for obs bin %d.\n",iObsBin);
    CalculateResultsRatiosObsBin(iObsBin,cCalcQACanvas);
  }

  // Calculate Differences
  for (int iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
    printf("Preparing to calculate differences for obs bin %d.\n",iObsBin);
    CalculateResultsDifferencesObsBin(iObsBin,cCalcQACanvas);
  }

}

void TaskCalcObservables::CalculateResultsObsBinEPBin(int iObsBin, int iEPBin, TCanvas * canv) {
  canv->Clear();

  // For NS, use Near Eta
  // For AS, use Full
  // If more than the Far eta nearside region is used for the background, this should be revisited.

  TH1D * histNearSide = fNearEtaDPhiProj_Sub[iObsBin][iEPBin];
  TH1D * histAwaySide = fFullDPhiProj_Sub[iObsBin][iEPBin];

  printf("Calculating Near Side Results with histogram %s\n",histNearSide->GetName());

  Double_t fYieldNS     = 0;
  Double_t fYieldNS_Err = 0;

  Double_t fRmsNS     = 0;
  Double_t fRmsNS_Err = 0;

//  histNearSide->GetXaxis()->SetRangeUser(-fYieldRangeNS,fYieldRangeNS);
 
  Int_t iMinRangeBin = histNearSide->FindBin(-fYieldRangeNS);
  Int_t iMaxRangeBin = histNearSide->FindBin(fYieldRangeNS);

  fYieldNS = histNearSide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldNS_Err,"width");

  histNearSide->SetMarkerColor(kNonSelectColor);
  histNearSide->SetLineColor(kNonSelectColor);
  histNearSide->DrawCopy();

  histNearSide->GetXaxis()->SetRangeUser(-fRmsRangeNS,fRmsRangeNS);
  fRmsNS     = histNearSide->GetRMS();
  fRmsNS_Err = histNearSide->GetRMSError();

  // Draw a clone with a different color after setting the range?
  histNearSide->SetMarkerColor(kSelectColor);
  histNearSide->SetLineColor(kSelectColor);
  histNearSide->Draw("SAME");

  printf("  Found Near-side Yield = %e \\pm %e and RMS = %f \\pm %f\n",fYieldNS,fYieldNS_Err,fRmsNS,fRmsNS_Err);

  double xValue = (fObsBins[iObsBin+1] + fObsBins[iObsBin]) / 2.;
  double xWidth = (fObsBins[iObsBin+1] - fObsBins[iObsBin]) / 2.;

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


  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_NearSide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_NearSide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_NearSide.C",fOutputDir.Data(),iObsBin,iEPBin));

  printf("Calculating Away Side Results with histogram %s\n",histAwaySide->GetName());

  Double_t fYieldAS     = 0;
  Double_t fYieldAS_Err = 0;

  Double_t fRmsAS     = 0;
  Double_t fRmsAS_Err = 0;

  iMinRangeBin = histAwaySide->FindBin(TMath::Pi()-fYieldRangeAS);
  iMaxRangeBin = histAwaySide->FindBin(TMath::Pi()+fYieldRangeAS);

  fYieldAS = histAwaySide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldAS_Err,"width");

  histAwaySide->SetMarkerColor(kNonSelectColor);
  histAwaySide->SetLineColor(kNonSelectColor);
  histAwaySide->DrawCopy();

  histAwaySide->GetXaxis()->SetRangeUser(TMath::Pi()-fRmsRangeAS,TMath::Pi()+fRmsRangeAS);
  fRmsAS     = histAwaySide->GetRMS();
  fRmsAS_Err = histAwaySide->GetRMSError();

  histAwaySide->SetMarkerColor(kSelectColor);
  histAwaySide->SetLineColor(kSelectColor);
  histAwaySide->Draw("SAME");

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

  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_AwaySide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CalcQA_ObsBin%d_EPBin%d_AwaySide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_AwaySide.C",fOutputDir.Data(),iObsBin,iEPBin));

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
  OutOverIn_AS->SetLineColor(kBlack);
  OutOverIn_AS->SetMarkerColor(kBlack);
  OutOverIn_AS->SetMarkerStyle(kFullSquare);

  iNumerator = 1; // 1 for Mid
  MidOverIn_AS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  MidOverIn_AS->SetLineColor(kGray);
  MidOverIn_AS->SetMarkerColor(kGray);
  MidOverIn_AS->SetMarkerStyle(kOpenSquare);

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  OutOverIn_NS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutOverIn_NS->SetLineColor(kBlack);
  OutOverIn_NS->SetMarkerColor(kBlack);
  OutOverIn_NS->SetMarkerStyle(kFullSquare);

  iNumerator = 1; // 1 for Mid
  MidOverIn_NS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  MidOverIn_NS->SetLineColor(kGray);
  MidOverIn_NS->SetMarkerColor(kGray);
  MidOverIn_NS->SetMarkerStyle(kOpenSquare);


  RMSOutOverIn_AS = new TGraphErrors(nObsBins);
  RMSOutOverIn_NS = new TGraphErrors(nObsBins);
  RMSMidOverIn_AS = new TGraphErrors(nObsBins);
  RMSMidOverIn_NS = new TGraphErrors(nObsBins);



  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  RMSOutOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSOutOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iNumerator = 1; // 1 for Mid
  RMSMidOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSMidOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  RMSOutOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSOutOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iNumerator = 1; // 1 for Mid
  RMSMidOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RMSMidOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));


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

//fNSYieldsInc
//fNSYieldsEP

  double xValue = fASYieldsEP[0]->GetX()[iObsBin];
  double xErr   = fASYieldsEP[0]->GetEX()[iObsBin];

  double YieldInPlane      = fASYieldsEP[0]->GetY()[iObsBin];
  double YieldInPlaneError = fASYieldsEP[0]->GetEY()[iObsBin];

  double YieldMidPlane      = fASYieldsEP[1]->GetY()[iObsBin];
  double YieldMidPlaneError = fASYieldsEP[1]->GetEY()[iObsBin];

  double YieldOutPlane      = fASYieldsEP[2]->GetY()[iObsBin];
  double YieldOutPlaneError = fASYieldsEP[2]->GetEY()[iObsBin];

  // OutOverIn_NS
  // MidOverIn_NS

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



//  TGraphErrors * RMSOutOverIn_AS;
//  TGraphErrors * RMSOutOverIn_NS;
//  TGraphErrors * RMSMidOverIn_AS;
//  TGraphErrors * RMSMidOverIn_NS;


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
  fNSYieldsInc->SetLineColor(kEPColorList[0]);
  fNSYieldsInc->SetMarkerColor(kEPColorList[0]);
  fNSYieldsInc->SetMarkerStyle(kEPMarkerList[0]);
  fNSRmsInc->SetLineColor(kEPColorList[0]);
  fNSRmsInc->SetMarkerColor(kEPColorList[0]);
  fNSRmsInc->SetMarkerStyle(kEPMarkerList[0]);
  fASYieldsInc->SetLineColor(kEPColorList[0]);
  fASYieldsInc->SetMarkerColor(kEPColorList[0]);
  fASYieldsInc->SetMarkerStyle(kEPMarkerList[0]);
  fASRmsInc->SetLineColor(kEPColorList[0]);
  fASRmsInc->SetMarkerColor(kEPColorList[0]);
  fASRmsInc->SetMarkerStyle(kEPMarkerList[0]);
  // By EP Bin
  for (int i = 0; i < kNEPBins; i++) {
    fNSYieldsEP[i]->SetLineColor(kEPColorList[i+1]);
    fNSYieldsEP[i]->SetMarkerColor(kEPColorList[i+1]);
    fNSYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i+1]);
    fNSRmsEP[i]->SetLineColor(kEPColorList[i+1]);
    fNSRmsEP[i]->SetMarkerColor(kEPColorList[i+1]);
    fNSRmsEP[i]->SetMarkerStyle(kEPMarkerList[i+1]);
    fASYieldsEP[i]->SetLineColor(kEPColorList[i+1]);
    fASYieldsEP[i]->SetMarkerColor(kEPColorList[i+1]);
    fASYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i+1]);
    fASRmsEP[i]->SetLineColor(kEPColorList[i+1]);
    fASRmsEP[i]->SetMarkerColor(kEPColorList[i+1]);
    fASRmsEP[i]->SetMarkerStyle(kEPMarkerList[i+1]);

  }

}


void TaskCalcObservables::DrawResults() {
  printf("Saving some results plots\n");


  SetGraphStyles();

  // Undecided about whether to have the same canvas sizes for yields and widths apparently
  TCanvas * cResults = new TCanvas("cResults","cResults",fYieldCanvasWidth,fYieldCanvasHeight);

  TMultiGraph * mgNSY = new TMultiGraph();
  TMultiGraph * mgNSRms = new TMultiGraph();
  TMultiGraph * mgASY = new TMultiGraph();
  TMultiGraph * mgASRms = new TMultiGraph();

  TString sDrawStyle = "AP";
  TLegend * legYields = new TLegend(0.55,0.45,0.85,0.85);
  TLegend * legRms = new TLegend(0.55,0.65,0.85,0.85);

  // Near-side Yields
  mgNSY->Add(fNSYieldsInc);
  legYields->AddEntry(fNSYieldsInc,"Inclusive","LP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legYields->AddEntry(fNSYieldsEP[iEPBin],fEPBinTitles[iEPBin],"LP");
    mgNSY->Add(fNSYieldsEP[iEPBin]);
  }
  mgNSY->Draw(sDrawStyle);
  mgNSY->GetXaxis()->SetTitle(fNSYieldsInc->GetXaxis()->GetTitle());
  mgNSY->GetYaxis()->SetTitle(fNSYieldsInc->GetTitle());
  legYields->Draw("SAME");

  cResults->Print(Form("%s/NearsideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYields.C",fOutputDir.Data()));

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

  cResults->Print(Form("%s/AwaysideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYields.C",fOutputDir.Data()));

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

  mgASYieldRatios->Add(OutOverIn_AS);
  mgASYieldRatios->Add(MidOverIn_AS);

  legASYieldRatios->AddEntry(OutOverIn_AS,"Out/In","LP");
  legASYieldRatios->AddEntry(MidOverIn_AS,"Mid/In","LP");

  mgASYieldRatios->GetXaxis()->SetTitle(OutOverIn_AS->GetXaxis()->GetTitle());
  mgASYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
  mgASYieldRatios->Draw(sDrawStyle);
  legASYieldRatios->Draw("SAME");

  cResults->Print(Form("%s/AwaysideYieldRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldRatios.C",fOutputDir.Data()));






  // OutOverIn and MidOverIn NS
  TMultiGraph * mgNSYieldRatios = new TMultiGraph();

  TLegend * legYieldRatios = new TLegend(0.55,0.65,0.85,0.85);

  mgNSYieldRatios->Add(OutOverIn_NS);
  mgNSYieldRatios->Add(MidOverIn_NS);

  legYieldRatios->AddEntry(OutOverIn_NS,"Out/In","LP");
  legYieldRatios->AddEntry(MidOverIn_NS,"Mid/In","LP");

  mgNSYieldRatios->GetXaxis()->SetTitle(OutOverIn_NS->GetXaxis()->GetTitle());
  mgNSYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
  mgNSYieldRatios->Draw(sDrawStyle);
  legYieldRatios->Draw("SAME");

  cResults->Print(Form("%s/NearsideYieldRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYieldRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYieldRatios.C",fOutputDir.Data()));

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

  CalculateResults();

  DrawResults();

  SaveOutput();

  cout<<"Done!"<<endl;
}

