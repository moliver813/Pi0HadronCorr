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
#include <TString.h>

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



void RemoveXErrorBars(TGraphErrors * graph) {
  int n = graph->GetN();
  for (int i = 0; i < n; i++) {
    graph->SetPointError(i,0.0,graph->GetEY()[i]);
  }
}

double GetMaximumYFromGraphErrors(TGraphErrors * graph) {
  double yMax = -1;
  for (int i = 0; i < (int) graph->GetN(); i++) {
    double localY = graph->GetY()[i] + graph->GetEY()[i];
    yMax = max(yMax,localY);
  }
  return yMax;
}


/// \cond CLASSIMP
ClassImp(TaskCalcObservables);


TaskCalcObservables::TaskCalcObservables() {
  gROOT->SetBatch(kTRUE);

  SetStyle();
}


void TaskCalcObservables::SetStyle() {
  gStyle->SetOptTitle(0);
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

  // Test
  gStyle->SetColorModelPS(0); // RGB
  //gStyle->SetColorModelPS(1); // CMYK


  //kOutInColor = TColor::GetFreeColorIndex();
  //TColor * OutInColor = new TColor(kOutInColor,172, 151, 212);
  //OutInColor = new TColor(kOutInColor,172, 151, 212);

  kModel0Color = 1+TColor::GetFreeColorIndex();
  cout<<"Color Settings: Using index "<<kModel0Color<<" for starting custom colors"<<endl;
  Model0Color = new TColor(kModel0Color,253./255, 127./255, 127./255);
// FIXME test
  //Model0Color = new TColor(kModel0Color,151, 212, 151);

  kOutInErrorColor = 2+TColor::GetFreeColorIndex();
  //TColor * OutInColor = new TColor(kOutInColor,172, 151, 212);
  OutInErrorColor = new TColor(kOutInErrorColor,172./255, 151./255, 212./255);
  
  // TEST
  kOutInErrorColorNS = 3+TColor::GetFreeColorIndex();
  OutInErrorColorNS = new TColor(kOutInErrorColorNS,151./255, 212./255, 151./255);


  //Int_t ci = TColor::GetFreeColorIndex();
  //TColor *color = new TColor(ci, 1., 1., 1.);



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
  else if (fObservable==2) fObservableName = "#it{p}_{T}^{assoc} (GeV/#it{c})";
  //else if (fObservable==2) fObservableName = "p_{T}^{a} (GeV/#it{c})";

  
  // Load weight histogram for debug info
  if (fIsMCGenMode) {
    hWeight = (TH1I *) fInputFileCentral->Get("hWeight");
    if (hWeight) {
      printf("MCGen EvGen weight histogram has %f entries.\n",hWeight->GetEntries());
    } else {
      printf("Did not find weight histogram in file.\n");
    }
  }




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

/** 
  * Load the systematic uncertainty objects from the output of SysCompare tasks
  */
void TaskCalcObservables::LoadSystematics() {
  cout<<"Loading Systematics"<<endl;

  if (fSystematicsNames.size() == 0) {
    printf("No systematic uncertainty files enterred. Skipping.\n");
    return;
  }

  // Use the array of gResultsGraph as defined in CleanResults.
  // use the name of the graph to readout from file?
  // Or just code for each observable. Yield ratios are the most important

  int nSysSources = fSystematicsNames.size();
  for (int i = 0; i < nSysSources; i++) {
    TGraphErrors * fSysErrorByType = 0;
    // =======================================================
    // Ratios
    // =======================================================
    // AS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("OutOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    OutOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // AS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("MidOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    MidOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("OutOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    OutOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("MidOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    MidOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);

    // RMS
    // AS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("RMSOutOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    RmsOutOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // AS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("RMSMidOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    RmsMidOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("RMSOutOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    RmsOutOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("RMSMidOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    RmsMidOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);


    // Sigmas
    // AS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("SigmasOutOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    SigmasOutOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // AS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("SigmasMidOverIn_AS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    SigmasMidOverIn_AS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Out/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("SigmasOutOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    SigmasOutOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);
    // NS Mid/In
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("SigmasMidOverIn_NS_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    SigmasMidOverIn_NS_SysErrorBySource.push_back(fSysErrorByType);


    // =======================================================
    // Direct Observables (yields, widths)
    // =======================================================
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("ASYieldsInc_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fASYieldsInc_SysErrorBySource.push_back(fSysErrorByType);

    vector<TGraphErrors *> fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("ASYieldsEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fASYieldsEP_SysErrorBySource.push_back(fLocalVector);

    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("NSYieldsInc_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fNSYieldsInc_SysErrorBySource.push_back(fSysErrorByType);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("NSYieldsEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fNSYieldsEP_SysErrorBySource.push_back(fLocalVector);




    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("ASRmsInc_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fASRmsInc_SysErrorBySource.push_back(fSysErrorByType);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("ASRmsEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fASRmsEP_SysErrorBySource.push_back(fLocalVector);

    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("NSRmsInc_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fNSRmsInc_SysErrorBySource.push_back(fSysErrorByType);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("NSRmsEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fNSRmsEP_SysErrorBySource.push_back(fLocalVector);

    // Sigmas
    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("ASSigmasInc_SysErr");
    if (!fSysErrorByType) {
      fprintf(stderr,"Error: could not find ASSigmasInc_SysErr in systematics file %d\n",i);
    }

    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fASSigmasInc_SysErrorBySource.push_back(fSysErrorByType);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("ASSigmasEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fASSigmasEP_SysErrorBySource.push_back(fLocalVector);

    fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get("NSSigmasInc_SysErr");
    fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
    fNSSigmasInc_SysErrorBySource.push_back(fSysErrorByType);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fSysErrorByType = (TGraphErrors *) fSystematicsFiles[i]->Get(Form("NSSigmasEP%d_SysErr",j));
      fSysErrorByType->SetName(Form("%s_%s",fSysErrorByType->GetName(),fSystematicsNames[i].Data()));
      fLocalVector.push_back(fSysErrorByType);
    }
    fNSSigmasEP_SysErrorBySource.push_back(fLocalVector);


  }

  // May want to transpose all arrays with the event plane.
  // fASYieldsEP_SysErrorBySource, fNSYieldsEP_SysErrorBySource
  // fASRmsEP_SysErrorBySource, fNSRmsEP_SysErrorBySource
  fASYieldsEP_SysErrorBySource = Transpose2DTGraphErrArray(fASYieldsEP_SysErrorBySource);
  fNSYieldsEP_SysErrorBySource = Transpose2DTGraphErrArray(fNSYieldsEP_SysErrorBySource);
  fASRmsEP_SysErrorBySource = Transpose2DTGraphErrArray(fASRmsEP_SysErrorBySource);
  fNSRmsEP_SysErrorBySource = Transpose2DTGraphErrArray(fNSRmsEP_SysErrorBySource);
  fASSigmasEP_SysErrorBySource = Transpose2DTGraphErrArray(fASSigmasEP_SysErrorBySource);
  fNSSigmasEP_SysErrorBySource = Transpose2DTGraphErrArray(fNSSigmasEP_SysErrorBySource);

  cout<<"Finished loading systematics"<<endl;
}

TGraphErrors * TaskCalcObservables::BuildConstantGraph(TGraphErrors * graph, double value) {
  TGraphErrors * newGraph = (TGraphErrors *) graph->Clone(Form("%s_TrackingUncertainty",graph->GetName()));
  for (int i = 0; i < newGraph->GetN(); i++) {
    double oldYValue = newGraph->GetY()[i];
    newGraph->SetPointError(i,0,value * oldYValue);
  }
  return newGraph;
}


//void TaskCalcObservables::AddTrackingUncertaintiesIndiv(TGraphErrors * graph, vector<TGraphErrors *> fErrorsArray, double value) {
//  TGraphErrors * fTrackingErrorGraph = BuildConstantGraph(graph,value);
//  fErrorsArray.push_back(fTrackingErrorGraph);
//}
void TaskCalcObservables::AddTrackingUncertaintiesIndiv(TGraphErrors * graph, vector<TGraphErrors *> * fErrorsArray, double value) {
  TGraphErrors * fTrackingErrorGraph = BuildConstantGraph(graph,value);
  fErrorsArray->push_back(fTrackingErrorGraph);
}


/**
  * Adds systematic uncertainties for the tracking
  * Global uncertainty added to yields
  * Event Plane Uncertainty applied to yield event plane ratios
  */
void TaskCalcObservables::AddTrackingUncertainties() {

  if (fIsMCGenMode) {
    printf("Not adding tracking uncertainty, as this is MCGen\n");
    return;
  }

  cout<<"Adding tracking uncertainties"<<endl;
  //if (fSystematicsNames.size() == 0) {
  if (fSystematicsFiles.size() == 0) {
    printf("Skipping adding tracking uncertainties because no other systematics have been added. This avoids an error.\n");
    return;
  }

  cout<<"Adding tracking uncertainties to ratios"<<endl;
//  AddTrackingUncertaintiesIndiv(OutOverIn_AS,OutOverIn_AS_SysErrorBySource,fTrackingEventPlaneUncertainty);
  AddTrackingUncertaintiesIndiv(OutOverIn_AS,&OutOverIn_AS_SysErrorBySource,fTrackingEventPlaneUncertainty);
  AddTrackingUncertaintiesIndiv(OutOverIn_NS,&OutOverIn_NS_SysErrorBySource,fTrackingEventPlaneUncertainty);
  
  
  cout<<"Adding tracking uncertainties to absolute yields"<<endl;
  AddTrackingUncertaintiesIndiv(fASYieldsInc,&fASYieldsInc_SysErrorBySource,fGlobalTrackingUncertainty);
  AddTrackingUncertaintiesIndiv(fNSYieldsInc,&fNSYieldsInc_SysErrorBySource,fGlobalTrackingUncertainty);
  for (int j = 0; j < kNEPBins; j++) {
    //AddTrackingUncertaintiesIndiv(fASYieldsEP[j],fASYieldsEP_SysErrorBySource[j],fGlobalTrackingUncertainty);
    AddTrackingUncertaintiesIndiv(fASYieldsEP[j],&fASYieldsEP_SysErrorBySource[j],fGlobalTrackingUncertainty);
    AddTrackingUncertaintiesIndiv(fNSYieldsEP[j],&fNSYieldsEP_SysErrorBySource[j],fGlobalTrackingUncertainty);
  }

 /* for (int j = 0; j < kNEPBins; j++) {
    fASYieldsEP_SysErrorBySource[j]
    fNSYieldsEP_SysErrorBySource[j]
    fASRmsEP_SysErrorBySource[j]
    fNSRmsEP_SysErrorBySource[j]
    fASSigmasEP_SysErrorBySource[j]
    fNSSigmasEP_SysErrorBySource[j]
  }*/
  
  cout<<"Done addint tracking efficiency uncertainties"<<endl;
  fSystematicsNames.push_back("Tracking Efficiency");
}

/** 
  * Load the systematic uncertainty objects from the output of SysCompare tasks
  */
void TaskCalcObservables::LoadModels() {
  cout<<"Loading models"<<endl;

  if (fModelNames.size() == 0) {
    printf("No Model files entered.\n");
    return;
  }
  int nModels = fModelNames.size();
  for (int i = 0; i < nModels; i++) {
    TGraphErrors * fModelGraph = 0;

    // =======================================================
    // Ratios
    // =======================================================

    // Yields

    // AS Out/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("OutOverIn_AS");
    if (fModelGraph == 0) {
      printf("Warning: could not find OutOverIn_AS in Model file %s\n",fModelNames[i].Data());
      return;
    }
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    OutOverIn_AS_Models.push_back(fModelGraph);
    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("MidOverIn_AS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    MidOverIn_AS_Models.push_back(fModelGraph);
    // NS Out/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("OutOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    OutOverIn_NS_Models.push_back(fModelGraph);
    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("MidOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    MidOverIn_NS_Models.push_back(fModelGraph);

    // RMS 
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("RMSOutOverIn_AS");
    if (fModelGraph == 0) {
      printf("Warning: could not find OutOverIn_AS in Model file %s\n",fModelNames[i].Data());
    }
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);

    //fModelGraph->SetFillStyle(0);
    //fModelGraph->SetFillColor(kModelColors[i]);
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);


    RmsOutOverIn_AS_Models.push_back(fModelGraph);
    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("RMSMidOverIn_AS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    RmsMidOverIn_AS_Models.push_back(fModelGraph);
    // NS Out/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("RMSOutOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    RmsOutOverIn_NS_Models.push_back(fModelGraph);
    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("RMSMidOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    RmsMidOverIn_NS_Models.push_back(fModelGraph);



    // Sigma Ratios
    // AS Out/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("SigmasOutOverIn_AS");
    if (fModelGraph == 0) {
      printf("Warning: could not find OutOverIn_AS in Model file %s\n",fModelNames[i].Data());
    }
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    SigmasOutOverIn_AS_Models.push_back(fModelGraph);

    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("SigmasMidOverIn_AS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    SigmasMidOverIn_AS_Models.push_back(fModelGraph);
    // NS Out/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("SigmasOutOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    SigmasOutOverIn_NS_Models.push_back(fModelGraph);
    // AS Mid/In
    fModelGraph = (TGraphErrors * ) fModelFiles[i]->Get("SigmasMidOverIn_NS");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fModelGraph->SetLineColor(kModelColors[i]);
    fModelGraph->SetFillColor(kModelColors[i]);
    //fModelGraph->SetFillColorAlpha(kModelColors[i],0.5);
    fModelGraph->SetMarkerColor(kModelColors[i]);
    fModelGraph->SetMarkerStyle(kModelStyles[i]);
    SigmasMidOverIn_NS_Models.push_back(fModelGraph);



    // =======================================================
    // Direct Observables (yields, widths)
    // =======================================================
    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("ASYieldsInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fASYieldsInc_Models.push_back(fModelGraph);

    vector<TGraphErrors *> fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("ASYieldsEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fASYieldsEP_Models.push_back(fLocalVector);

    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("NSYieldsInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fNSYieldsInc_Models.push_back(fModelGraph);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("NSYieldsEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fNSYieldsEP_Models.push_back(fLocalVector);

    /*fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("NSRmsEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fNSRmsEP_Models.push_back(fLocalVector);
*/
    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("ASRmsInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fASRmsInc_Models.push_back(fModelGraph);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("ASRmsEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fASRmsEP_Models.push_back(fLocalVector);

    // AS Sigma
    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("ASSigmasInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fASSigmasInc_Models.push_back(fModelGraph);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("ASSigmasEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fASSigmasEP_Models.push_back(fLocalVector);

    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("NSRmsInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fNSRmsInc_Models.push_back(fModelGraph);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("NSRmsEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fNSRmsEP_Models.push_back(fLocalVector);

    // NS Sigma
    fModelGraph = (TGraphErrors *) fModelFiles[i]->Get("NSSigmasInc");
    fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
    fNSSigmasInc_Models.push_back(fModelGraph);

    fLocalVector = {};
    for (int j = 0; j < kNEPBins; j++) {
      fModelGraph = (TGraphErrors *) fModelFiles[i]->Get(Form("NSSigmasEP%d",j));
      fModelGraph->SetName(Form("%s_%s",fModelGraph->GetName(),fModelNames[i].Data()));
      fLocalVector.push_back(fModelGraph);
    }
    fNSSigmasEP_Models.push_back(fLocalVector);

  }

  fASYieldsEP_Models = Transpose2DTGraphErrArray(fASYieldsEP_Models);
  fNSYieldsEP_Models = Transpose2DTGraphErrArray(fNSYieldsEP_Models);
  fASRmsEP_Models = Transpose2DTGraphErrArray(fASRmsEP_Models);
  fNSRmsEP_Models = Transpose2DTGraphErrArray(fNSRmsEP_Models);
  fASSigmasEP_Models = Transpose2DTGraphErrArray(fASSigmasEP_Models);
  fNSSigmasEP_Models = Transpose2DTGraphErrArray(fNSSigmasEP_Models);

  cout<<"Finished loading models"<<endl;
}

TGraphErrors * TaskCalcObservables::CalculateSystematicIndiv(TGraphErrors * graph, vector<TGraphErrors *> graph_sys_errors) {
//TGraphErrors * TaskCalcObservables::CalculateSystematicIndiv(TGraphErrors * graph, vector<TGraphErrors *> graph_sys_errors, TF1 * UncertFit) {

  int nSysSources = graph_sys_errors.size();

  if (graph == 0) {
    printf("No graph given\n");
    return 0;
  }

  TGraphErrors * graph_sys = (TGraphErrors *) graph->Clone(Form("%s_SysError",graph->GetName()));
  // Second graph that has uncertainty as value
  TGraphErrors * graph_sys_only = (TGraphErrors *) graph->Clone(Form("%s_SysErrorOnly",graph->GetName()));

  if (nSysSources == 0) {
    printf("No systematic uncertainty files entered. Skipping Calculation for graph %s.\n",graph->GetName());
    return graph_sys;
    //return 0;
  }

  printf("Calculating Systematics for %s, which has %d points\n",graph->GetName(),graph->GetN());
  //printf("\tGiven %d sources.\n",(int) graph_sys_errors.size());

  float minXDefault = 0;
  float maxXDefault = 17;

  TString sFitName = Form("%s_UncertFit",graph->GetName());
  TString sFormula = "[0]";

  TString sBernPol3 = Form("[0]*TMath::Power(%f-x,3) + 3*[1]*x*TMath::Power(%f-x,2)\
  + 3*[2]*x*x*(%f-x) + [3]*TMath::Power(x,3)",maxXDefault,maxXDefault,maxXDefault);

  sFormula = sBernPol3;


//  TGraphErrors * graph_sys = (TGraphErrors *) graph->Clone(Form("%s_SysError",graph->GetName()));
  // Second graph that has uncertainty as value
//  TGraphErrors * graph_sys_only = (TGraphErrors *) graph->Clone(Form("%s_SysErrorOnly",graph->GetName()));
  for (int i = 0; i < nObsBins-nSkipPoints; i++) {
    // FIXME need to deal with deleted points
    double initX    = graph->GetX()[i];
    double initXErr = graph->GetEX()[i];
    double initY    = graph->GetY()[i];
    //double initYErr = graph->GetEY()[i];

    double finalYErr = 0;
    
    for (int j = 0; j < nSysSources; j++) {
      double localYErr = graph_sys_errors[j]->GetEY()[i];
      printf("\tAdding localYErr = %f\n",localYErr);
      finalYErr += localYErr * localYErr;
    }
    printf("Systematics: sum of squared = %f\n",finalYErr);
    if (finalYErr > 0) finalYErr = sqrt(finalYErr);
    graph_sys->SetPoint(i,initX,initY);
    graph_sys->SetPointError(i,initXErr,finalYErr);

  }
  //graph_sys_only->Fit(fit,"+");
  graph_sys->SetLineColor(kViolet);
  return graph_sys;
}

/**
  * Add the errors from the SysErrorBySource arrays in quadrature to observables like the
  * OutOverIn_AS and create the OutOverIn_AS_SysError tgraphs
  */
void TaskCalcObservables::CalculateSystematics() {

  int nSysSources = fSystematicsNames.size();

  if (nSysSources == 0) {
    printf("No systematic uncertainty files entered. Skipping Calculation.\n");
    return;
  }


  // Start with OutOverIn_AS
  TGraphErrors * graph = OutOverIn_AS;
  vector<TGraphErrors *> graph_sys_errors = OutOverIn_AS_SysErrorBySource;
  OutOverIn_AS_SysError = CalculateSystematicIndiv(OutOverIn_AS,OutOverIn_AS_SysErrorBySource);
  //OutOverIn_AS_SysError->SetFillColorAlpha(kOutInErrorColor,kOutInErrorAlpha);

  OutOverIn_AS_SysError->SetFillColor(kOutInErrorColor);
  //OutOverIn_AS_SysError->SetFillColorAlpha(kOutInErrorColor,0.5); // FIXME test

  OutOverIn_AS_SysError->SetFillStyle(kOutInSysFillStyle);
  OutOverIn_AS_SysError->SetMarkerColor(kOutInErrorColor);
  OutOverIn_AS_SysError->SetLineColor(kOutInErrorColor);
  //OutOverIn_AS_SysError->SetMarkerColor(kBlue-9);
  //OutOverIn_AS_SysError->SetLineColor(kBlue-9);

  // MidOverIn_AS
  graph = MidOverIn_AS;
  graph_sys_errors = MidOverIn_AS_SysErrorBySource;
  MidOverIn_AS_SysError = CalculateSystematicIndiv(MidOverIn_AS, MidOverIn_AS_SysErrorBySource);
  //MidOverIn_AS_SysError->SetFillColorAlpha(kOutInErrorColor,kOutInErrorAlpha);
  MidOverIn_AS_SysError->SetFillColor(kMidInErrorColor);
  MidOverIn_AS_SysError->SetFillStyle(kOutInSysFillStyle);
  MidOverIn_AS_SysError->SetMarkerColor(kBlue-9);
  
  // OutOverIn_NS
  graph = OutOverIn_NS;
  graph_sys_errors = OutOverIn_NS_SysErrorBySource;
  OutOverIn_NS_SysError = CalculateSystematicIndiv(OutOverIn_NS,OutOverIn_NS_SysErrorBySource);
  //OutOverIn_NS_SysError->SetFillColorAlpha(kOutInErrorColorNS,kOutInErrorAlpha);
  OutOverIn_NS_SysError->SetFillColor(kOutInErrorColorNS);
  OutOverIn_NS_SysError->SetFillStyle(kOutInSysFillStyle);
  OutOverIn_NS_SysError->SetMarkerColor(kOutInErrorColorNS);
  OutOverIn_NS_SysError->SetLineColor(kOutInErrorColorNS);

  // MidOverIn_NS
  graph = MidOverIn_NS;
  graph_sys_errors = MidOverIn_NS_SysErrorBySource;
  MidOverIn_NS_SysError = CalculateSystematicIndiv(MidOverIn_NS,MidOverIn_NS_SysErrorBySource);
  //MidOverIn_NS_SysError->SetFillColorAlpha(kMidInErrorColorNS,kMidInErrorAlpha);
  MidOverIn_NS_SysError->SetFillColor(kMidInErrorColorNS);
  //MidOverIn_NS_SysError->SetFillStyle(kOutInSysFillStyle);
  MidOverIn_NS_SysError->SetMarkerColor(kMidInErrorColorNS);

  // =======================================================================================================
  // RMS
  // =======================================================================================================
  // RmsOutOverIn_AS
  graph = RmsOutOverIn_AS;
  graph_sys_errors = RmsOutOverIn_AS_SysErrorBySource;
  RmsOutOverIn_AS_SysError = CalculateSystematicIndiv(RmsOutOverIn_AS,RmsOutOverIn_AS_SysErrorBySource);
  RmsOutOverIn_AS_SysError->SetFillColor(kOutInErrorColor);
  //RmsOutOverIn_AS_SysError->SetFillColorAlpha(kOutInErrorColor,kOutInErrorAlpha);
 // RmsOutOverIn_AS_SysError->SetFillColor(kOutInSysFillStyle);
  RmsOutOverIn_AS_SysError->SetMarkerColor(kOutInErrorColor);

  // RmsMidOverIn_AS
  graph = RmsMidOverIn_AS;
  graph_sys_errors = RmsMidOverIn_AS_SysErrorBySource;
  RmsMidOverIn_AS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  RmsMidOverIn_AS_SysError->SetFillColor(kMidInErrorColor);
  //RmsMidOverIn_AS_SysError->SetFillColorAlpha(kMidInErrorColor,kMidInErrorAlpha);
 // RmsMidOverIn_AS_SysError->SetFillColor(kOutInSysFillStyle);
  RmsMidOverIn_AS_SysError->SetMarkerColor(kMidInErrorColor);

  // RmsOutOverIn_NS
  graph = RmsOutOverIn_NS;
  graph_sys_errors = RmsOutOverIn_NS_SysErrorBySource;
  RmsOutOverIn_NS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
 // RmsOutOverIn_NS_SysError->SetFillColor(kOutInErrorColor);
 // RmsOutOverIn_NS_SysError->SetFillColor(kOutInSysFillStyle);
  RmsOutOverIn_NS_SysError->SetMarkerColor(kBlue-9);

  // RmsMidOverIn_NS
  graph = RmsMidOverIn_NS;
  graph_sys_errors = RmsMidOverIn_NS_SysErrorBySource;
  RmsMidOverIn_NS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
 // RmsMidOverIn_NS_SysError->SetFillColor(kOutInErrorColor);
 // RmsMidOverIn_NS_SysError->SetFillColor(kOutInSysFillStyle);
  RmsMidOverIn_NS_SysError->SetMarkerColor(kBlue-9);

  // =======================================================================================================
  // Sigmas
  // =======================================================================================================
  // SigmasOutOverIn_AS
  graph = SigmasOutOverIn_AS;
  graph_sys_errors = SigmasOutOverIn_AS_SysErrorBySource;
  SigmasOutOverIn_AS_SysError = CalculateSystematicIndiv(SigmasOutOverIn_AS,SigmasOutOverIn_AS_SysErrorBySource);
  SigmasOutOverIn_AS_SysError->SetFillColor(kOutInErrorColor);
  //SigmasOutOverIn_AS_SysError->SetFillColorAlpha(kOutInErrorColor,kOutInErrorAlpha);
 // SigmasOutOverIn_AS_SysError->SetFillColor(kOutInSysFillStyle);
  SigmasOutOverIn_AS_SysError->SetMarkerColor(kOutInErrorColor);

  // SigmasMidOverIn_AS
  graph = SigmasMidOverIn_AS;
  graph_sys_errors = SigmasMidOverIn_AS_SysErrorBySource;
  SigmasMidOverIn_AS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  SigmasMidOverIn_AS_SysError->SetFillColor(kMidInErrorColor);
  //SigmasMidOverIn_AS_SysError->SetFillColorAlpha(kMidInErrorColor,kMidInErrorAlpha);
 // SigmasMidOverIn_AS_SysError->SetFillColor(kOutInSysFillStyle);
  SigmasMidOverIn_AS_SysError->SetMarkerColor(kMidInErrorColor);

  // SigmasOutOverIn_NS
  graph = SigmasOutOverIn_NS;
  graph_sys_errors = SigmasOutOverIn_NS_SysErrorBySource;
  SigmasOutOverIn_NS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
 // SigmasOutOverIn_NS_SysError->SetFillColor(kOutInErrorColor);
 // SigmasOutOverIn_NS_SysError->SetFillColor(kOutInSysFillStyle);
  SigmasOutOverIn_NS_SysError->SetMarkerColor(kBlue-9);

  // SigmasMidOverIn_NS
  graph = SigmasMidOverIn_NS;
  graph_sys_errors = SigmasMidOverIn_NS_SysErrorBySource;
  SigmasMidOverIn_NS_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
 // SigmasMidOverIn_NS_SysError->SetFillColor(kOutInErrorColor);
 // SigmasMidOverIn_NS_SysError->SetFillColor(kOutInSysFillStyle);
  SigmasMidOverIn_NS_SysError->SetMarkerColor(kBlue-9);








  cout<<"Starting calculating systematics for the direct observables"<<endl;

  // =======================================================================================================
  // Direct Observables
  // =======================================================================================================
  graph = fNSYieldsInc;
  graph_sys_errors = fNSYieldsInc_SysErrorBySource;
  fNSYieldsInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
//  fNSYieldsInc_SysError->SetFillColorAlpha(,);
//  fNSYieldsInc_SysError->SetMarkerColor(,);
  printf("HELO\n");
  vector<TGraphErrors *> fLocalVector = {};

  vector<vector<TGraphErrors * >> fErrorsEPVector = {};

  for (int i = 0; i < kNEPBins; i++) {
    graph = fNSYieldsEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fNSYieldsEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fNSYieldsEP_SysError.push_back(fLocalSysError);
  }

  printf("HELO 2\n");
  graph = fASYieldsInc;
  graph_sys_errors = fASYieldsInc_SysErrorBySource;
  fASYieldsInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  fLocalVector = {};
  for (int i = 0; i < kNEPBins; i++) {
    graph = fASYieldsEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fASYieldsEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fASYieldsEP_SysError.push_back(fLocalSysError);
  }

  graph = fNSRmsInc;
  graph_sys_errors = fNSRmsInc_SysErrorBySource;
  fNSRmsInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  fLocalVector = {};
  for (int i = 0; i < kNEPBins; i++) {
    graph = fNSRmsEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fNSRmsEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fNSRmsEP_SysError.push_back(fLocalSysError);
  }

  graph = fASRmsInc;
  graph_sys_errors = fASRmsInc_SysErrorBySource;
  fASRmsInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  fLocalVector = {};
  for (int i = 0; i < kNEPBins; i++) {
    graph = fASRmsEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fASRmsEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fASRmsEP_SysError.push_back(fLocalSysError);
  }
  printf("Finished Yield and RMS Systematics\n");

  graph = fNSSigmasInc;
  graph_sys_errors = fNSSigmasInc_SysErrorBySource;
  fNSSigmasInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  fLocalVector = {};
  for (int i = 0; i < kNEPBins; i++) {
    graph = fNSSigmasEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fNSSigmasEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fNSSigmasEP_SysError.push_back(fLocalSysError);
  }
  printf("Finished NS Sigma Systematics\n");
  graph = fASSigmasInc;
  graph_sys_errors = fASSigmasInc_SysErrorBySource;
  fASSigmasInc_SysError = CalculateSystematicIndiv(graph,graph_sys_errors);
  fLocalVector = {};
  for (int i = 0; i < kNEPBins; i++) {
    graph = fASSigmasEP[i];
    graph_sys_errors = {};
    fErrorsEPVector = fASSigmasEP_SysErrorBySource;
    if (fErrorsEPVector.size() > 0) {
      graph_sys_errors = fErrorsEPVector[i];
    }
    TGraphErrors * fLocalSysError = CalculateSystematicIndiv(graph,graph_sys_errors);
    fASSigmasEP_SysError.push_back(fLocalSysError);
  }

  printf("Finished AS Sigma Systematics\n");




}

void TaskCalcObservables::InitArrays() {

  //Double_t fZtStep = 1.0/(7 - 1.0);
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


  // Load EPR (code copied from phase4)
   // EPR Set EP2, from QnVector MB (Train unknown, prior to T59)
  //  Double_t fEPRes_Set_0[4][6]  = {
   //                   {  0.765960, 0.619163,  0.509267, 0.348666, 0.318429, 0.187868},  
   //                   {  0.858157, 0.822691, 0.692985, 0.580624,  0.502229,  0.375755},
   //                   {  0.832549,  0.771133,  0.639423,  0.507014,  0.439729,  0.305388},
   //                   {  0.704550,  0.445893,  0.380824,  0.196809,  0.211605,  0.084895}};
    // EPR Set EP2, from QnVector MB (T59)
    Double_t fEPRes_Set_0[4][6]  = {
    //{      R_1,      R_2,      R_3,      R_4,      R_5,      R_6}
      { 0.765625, 0.619035, 0.509148, 0.348548, 0.318585, 0.187958},
      { 0.858052, 0.822638, 0.693049, 0.580605, 0.502134, 0.375592},
      { 0.832519, 0.771184, 0.639522, 0.507076, 0.439909, 0.305442},
      { 0.704443, 0.445964, 0.380898, 0.196861, 0.211537, 0.084192}
    };
    Double_t fEPRes_Set_0_Err[4][6]  = {
      { 0.000434, 0.000548, 0.000690, 0.000979, 0.001097, 0.001653},
      { 0.000262, 0.000233, 0.000343, 0.000435, 0.000517, 0.000681},
      { 0.000275, 0.000272, 0.000380, 0.000495, 0.000581, 0.000791},
      { 0.000271, 0.000424, 0.000516, 0.000847, 0.000876, 0.001653}
    };


    // EPR Set EP2, from QnVector EGA (T59)
    Double_t fEPRes_Set_1[4][6]  = {
                      {  0.767842, 0.622406, 0.525535, 0.390981, 0.352418, 0.236852},
                      {  0.781089, 0.651664, 0.549970, 0.423973, 0.376066, 0.259773},
                      {  0.779045, 0.646402, 0.544922, 0.414330, 0.369320, 0.252558},
                      {  0.762544, 0.610570, 0.515899, 0.376056, 0.343955, 0.226055}};
    Double_t fEPRes_Set_1_Err[4][6]  = {
      { 0.000986, 0.001280, 0.001574, 0.002121, 0.002375, 0.003367},
      { 0.000894, 0.001117, 0.001381, 0.001821, 0.002049, 0.002843},
      { 0.001326, 0.001662, 0.002064, 0.002738, 0.003081, 0.004337},
      { 0.001927, 0.002536, 0.003097, 0.004162, 0.004653, 0.006587}
    };

    // Read off from some plot
    Double_t fEPRes_Set_2[4][6]  = {
                      {  0,  0.73,  0.62, 0.275,  0.0,  0.0},
                      {  0, 0.885, 0.605, 0.245,  0.0,  0.0},
                      {  0,  0.85,  0.49,  0.21,  0.0,  0.0},
                      {  0,  0.58,  0.22,  0.08,  0.0,  0.0}  };

    // MCGen (EPR = 1)
    Double_t fEPRes_Set_3[4][6]  = {
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0}  };
    Double_t fEPRes_Set_3_Err[4][6]  = {
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0}  };


    // EPR Set EP3 from QnVector MB (T59)
    Double_t fEP3Res_Set_0[4][6]  = {
                      {  0.845080, 0.525658, 0.393520, 0.391907, 0.259402, 0.139690},
                      {  0.842722, 0.513905, 0.372819, 0.375413, 0.246513, 0.125048},
                      {  0.833544, 0.461511, 0.265899, 0.299722, 0.200492, 0.065516},
                      {  0.825777, 0.407792, 0.087980, 0.204440, 0.162711, 0.014950}};

    // EPR Set EP3 from QnVector EGA (T59)
    Double_t fEP3Res_Set_1[4][6]  = {
                      {  0.835607, 0.469750, 0.275954, 0.307249, 0.209391, 0.087959},
                      {  0.834803, 0.466532, 0.271198, 0.304266, 0.209051, 0.092139},
                      {  0.834278, 0.461766, 0.257766, 0.294599, 0.200117, 0.075660},
                      {  0.832765, 0.452335, 0.234578, 0.280758, 0.201374, 0.092235}};


    // EPR Set EP4 from QnVector MB (T59)
    Double_t fEP4Res_Set_0[4][6]  = {
                      {  0.926312, 0.712967, 0.384993, 0.141703, 0.308486, 0.569209},
                      {  0.925901, 0.711299, 0.380526, 0.133776, 0.308111, 0.574361},
                      {  0.924596, 0.706943, 0.369987, 0.095097, 0.306386, 0.600262},
                      {  0.923124, 0.703032, 0.365919, 0.029369, 0.330021, 0.652974}};

    // EPR Set EP4 from QnVector EGA (T59)
    Double_t fEP4Res_Set_1[4][6]  = {
                      {  0.924481, 0.707366, 0.374302, 0.101455, 0.312070, 0.598036},
                      {  0.924401, 0.706881, 0.372533, 0.097955, 0.312920, 0.605455},
                      {  0.924679, 0.708054, 0.375002, 0.096597, 0.312343, 0.609336},
                      {  0.924669, 0.708298, 0.375728, 0.090440, 0.311984, 0.603048}};


    // MC (1.0)
    Double_t fEP3Res_Set_3[4][6]  = {
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
      Double_t fEP3Res_Set_3_Err[4][6]  = {
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
    Double_t fEP4Res_Set_3[4][6]  = {
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
      Double_t fEP4Res_Set_3_Err[4][6]  = {
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};


    // Era 2
    // 18qr MB JEQVector
    // Only have Cent2 for now
    Double_t fEPRes_Set_4[4][6]  = {
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.821406, 0.749172, 0.622964, 0.491571 , 0.428120, 0.296947},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0}  };
    Double_t fEPRes_Set_4_Err[4][6]  = {
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.000138, 0.000143, 0.000194, 0.000250, 0.000291, 0.000392},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0}  };

    // Era 2 
    // 18qr EGA JEQVector
    Double_t fEPRes_Set_5[4][6]  = {
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {  0.791194, 0.685201, 0.572346, 0.433487, 0.387502, 0.269563},
                      {  0, 1.0, 1.0, 1.0, 1.0, 1.0}  };
    Double_t fEPRes_Set_5_Err[4][6]  = {
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.000525, 0.000604, 0.000779, 0.001023, 0.001165, 0.001567},
                      {  0, 0.0, 0.0, 0.0, 0.0, 0.0}  };




    if (fIsMCGenMode) iEPRSet = 3;



    switch (iEPRSet) {
      default:
      case 0: // From QnVector MB
        printf("Using EPR Set 0 (MB)\n");
        memcpy(fEPRes,fEPRes_Set_0[iCentBin], sizeof(fEPRes));
        memcpy(fEPRes_Err,fEPRes_Set_0_Err[iCentBin], sizeof(fEPRes_Err));
        memcpy(fEP3Res,fEP3Res_Set_0[iCentBin], sizeof(fEP3Res));
        memcpy(fEP4Res,fEP4Res_Set_0[iCentBin], sizeof(fEP4Res));


       /* fEPRes[0] = {  0,  0.73,  0.62, 0.275,  0.0,  0.0};  
        fEPRes[1] = {  0, 0.885, 0.605, 0.245,  0.0,  0.0};
        fEPRes[2] = {  0,  0.85,  0.49,  0.21,  0.0,  0.0};
        fEPRes[3] = {  0,  0.58,  0.22,  0.08,  0.0,  0.0}; */
        break;
      case 1: // From QnVector EGA 
        printf("Using EPR Set 1 (EGA)\n");
        memcpy(fEPRes,fEPRes_Set_1[iCentBin], sizeof(fEPRes));
        memcpy(fEPRes_Err,fEPRes_Set_1_Err[iCentBin], sizeof(fEPRes_Err));
        memcpy(fEP3Res,fEP3Res_Set_1[iCentBin], sizeof(fEP3Res));
        memcpy(fEP4Res,fEP4Res_Set_1[iCentBin], sizeof(fEP4Res));

        //old// From Raymond's analysis. Cent0M with QnVector correction
        break;
      case 2: // read off of that one graph. TPC R_n values
        memcpy(fEPRes,fEPRes_Set_2[iCentBin], sizeof(fEPRes));
        break;
      case 3: // Full resolution. Ideal for MC
        printf("Using EPR Set 3 (MCGen, All EPRs set to 1)\n");
        memcpy(fEPRes,fEPRes_Set_3[iCentBin], sizeof(fEPRes));
        memcpy(fEPRes_Err,fEPRes_Set_3_Err[iCentBin], sizeof(fEPRes_Err));
        memcpy(fEP3Res,fEPRes_Set_3[iCentBin], sizeof(fEP3Res));
        memcpy(fEP4Res,fEPRes_Set_3[iCentBin], sizeof(fEP4Res));
        break;
      case 4: // 18qr MB JEQVector
        printf("Using EPR calculated from 18qr Pass3 MB+SC V0M\n");
        memcpy(fEPRes,fEPRes_Set_4[iCentBin], sizeof(fEPRes));
        memcpy(fEPRes_Err,fEPRes_Set_4_Err[iCentBin], sizeof(fEPRes_Err));
        memcpy(fEP3Res,fEPRes_Set_4[iCentBin], sizeof(fEP3Res));
        memcpy(fEP4Res,fEPRes_Set_4[iCentBin], sizeof(fEP4Res));
        break;
      case 5: // 18qr EGA JEQVector
        printf("Using EPR calculated from 18qr Pass3 EGA-trigger V0M\n");
        memcpy(fEPRes,fEPRes_Set_5[iCentBin], sizeof(fEPRes));
        memcpy(fEPRes_Err,fEPRes_Set_5_Err[iCentBin], sizeof(fEPRes_Err));
        memcpy(fEP3Res,fEPRes_Set_5[iCentBin], sizeof(fEP3Res));
        memcpy(fEP4Res,fEPRes_Set_5[iCentBin], sizeof(fEP4Res));
        break;
    }

    for (Int_t i = 0; i < 6; i++) {
      printf("Loading resolution R_%d = %f\n",i+1,fEPRes[i]);
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
  fASYieldsInc->GetXaxis()->SetTitle(fObservableName.Data());
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

  fNSSigmasInc = new TGraphErrors(nObsBins);
  fNSSigmasInc->SetName("NSSigmasInc");
  fNSSigmasInc->SetTitle("Near-side width");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("NSSigmasEP%d",i));
    newGraph->SetTitle(Form("Near-side Sigmas (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fNSSigmasEP.push_back(newGraph);
  }
  fASSigmasInc = new TGraphErrors(nObsBins);
  fASSigmasInc->SetName("ASSigmasInc");
  fASSigmasInc->SetTitle("Away-side width");
  for (int i = 0; i < kNEPBins; i++) {
    TGraphErrors * newGraph = new TGraphErrors(nObsBins);
    newGraph->SetName(Form("ASSigmasEP%d",i));
    newGraph->SetTitle(Form("Away-side Sigmas (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
    fASSigmasEP.push_back(newGraph);
  }

  // Initialize the variant TGraphErrors
  if (fFullDPhiProj_Sub_RPFVar[0][0] != 0) {
    nRPFVariants = fFullDPhiProj_Sub_RPFVar[0][0]->GetNbinsY();
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {

      vector<TGraphErrors *> fNSRmsEPVar = {};
      vector<TGraphErrors *> fNSSigmasEPVar = {};
      vector<TGraphErrors *> fNSYieldsEPVar = {};
      vector<TGraphErrors *> fASRmsEPVar = {};
      vector<TGraphErrors *> fASYieldsEPVar = {};
      vector<TGraphErrors *> fASSigmasEPVar = {};

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


      TGraphErrors * fNSSigmasIncVariant = new TGraphErrors(nObsBins);
      fNSSigmasIncVariant->SetName(Form("NSSigmasInc_RPFVar%d",iVar));
      fNSSigmasIncVariant->SetTitle("Near-side width");
      fNSSigmasInc_RPFVariants.push_back(fNSSigmasIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("NSSigmasEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Near-side Sigmas (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fNSSigmasEPVar.push_back(newGraph);
      }
      fNSSigmasEP_RPFVariants.push_back(fNSSigmasEPVar);
      TGraphErrors * fASSigmasIncVariant = new TGraphErrors(nObsBins);
      fASSigmasIncVariant->SetName(Form("ASSigmasInc_RPFVar%d",iVar));
      fASSigmasIncVariant->SetTitle("Away-side width");
      fASSigmasInc_RPFVariants.push_back(fASSigmasIncVariant);
      for (int i = 0; i < kNEPBins; i++) {
        TGraphErrors * newGraph = new TGraphErrors(nObsBins);
        newGraph->SetName(Form("ASSigmasEP%d_RPFVar%d",i,iVar));
        newGraph->SetTitle(Form("Away-side Sigmas (%s);%s",fEPBinTitles[i].Data(),fObservableName.Data()));
        fASSigmasEPVar.push_back(newGraph);
      }
      fASSigmasEP_RPFVariants.push_back(fASSigmasEPVar);



    }
  }

}



void TaskCalcObservables::DrawDirectComparisons() {
  cout<<"Drawing Direct Comparisons"<<endl;

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


  cout << "Done drawing direct comparisons" << endl;
}

/*

 * iSide = 0 -> nearside
 * iSide = 1 -> awayside
 */
TF1 * TaskCalcObservables::FitAndCalculateSigmaFromHisto(TH1D * hist, double * Sigma, double * SigmaErr, int iSide) {

  //printf("Calculating sigma for histogram %s\n",hist->GetName());

  // Basic gaussian
  TString sForm1Gaus = "[0]*TMath::Gaus(x,0,[1]) + [2]";
  if (iSide == 1) {
    sForm1Gaus = "[0]*TMath::Gaus(x,TMath::Pi(),[1]) + [2]";
  }

  TString sForm1GausGaus = "[0]*TMath::Gaus(x,0,[1]) + [2]*TMath::Gaus(x,TMath::Pi(),[3])+[4]";
  // Add the wrap around peaks:
  sForm1GausGaus += "[0]*TMath::Gaus(x,2*TMath::Pi(),[1]) + [2]*TMath::Gaus(x,-TMath::Pi(),[3])";


  TString sGenGausAS = "[2] + [0] * (2.*[3] * pow(TMath::Gamma(3./[3]),0.5)) / (2.*[1]*pow(TMath::Gamma(1./[3]),1.5)) * TMath::Exp(-pow(abs(x-TMath::Pi())/([1]*pow(TMath::Gamma(1./[3])/TMath::Gamma(3./[3]),0.5)),[3]))";



  TString sFitFunction = sForm1Gaus;

  double fFuncMin = -TMath::Pi()/2;
  double fFuncMax = TMath::Pi()/2;

  if (iSide == 1) {
    fFuncMin = TMath::Pi()/2;
    fFuncMax = 3.*TMath::Pi()/2;
  }



  if (iSide == 1) {
    switch (iAwaysideFitOption) {
      case 1: // Gen Gaus
//[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.*[6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) +TMath::Exp(-pow(abs(x-3*TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6]))  ) + [7]"
       // sFitFunction = "[0]*TMath::Gaus(x,TMath::Pi(),[1]) + [2]";
        //sFitFunction = "[2] + [0] * (2.*[3] * pow(TMath::Gamma(3./[3]),0.5)) / (2.*[1]*pow(TMath::Gamma(1./[3]),1.5)) * TMath::Exp(-pow(abs(x-TMath::Pi())/([1]*pow(TMath::Gamma(1./[3])/TMath::Gamma(3./[6]),0.5)),[3]))";
        sFitFunction = sGenGausAS;
        printf("hello there\n");
        break;
      default:
        break;
    }


  }


  //printf("Debug sigma fit: hist has range %f - %f\n",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());




  int iSigmaPar = 1;
  int iPedestalPar = 2;

  TF1 * fit = new TF1(Form("%s_Fit",hist->GetName()),sFitFunction,fFuncMin,fFuncMax);

  fit->SetParName(iSigmaPar,"sigma");

  // Initialization
  double fMax = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParameter(0,fMax);
  fit->SetParameter(iSigmaPar,0.05);
  fit->SetParameter(iPedestalPar,0);

  fit->SetParLimits(0,0.,3*fMax);
  fit->SetParLimits(iSigmaPar,0.001,2.0);
  fit->SetParLimits(iPedestalPar,0.,fMax);

  // this should be fixed to 0, but temporarily allow pedestal, mostly for JEWEL KR.
  //if (!fIsMCGenMode) fit->FixParameter(iPedestalPar,0.);
  // Now KR DEtascaling7 removes pedestal correctly
  fit->FixParameter(iPedestalPar,0.);

  // Additional settings

  if (iSide == 1) {
    switch (iAwaysideFitOption) {
      case 1: // Gen Gaus 
        fit->SetParName(3,"beta");
        fit->SetParameter(3,2.);
        fit->SetParLimits(3,1.,4.);
        printf("hello there again\n");
      break;
      default:
      break;
    }


  }



  hist->Fit(fit,"");

  *Sigma = fit->GetParameter(iSigmaPar);
  *SigmaErr = fit->GetParError(iSigmaPar);

  TCanvas * cQA = new TCanvas();
  hist->GetXaxis()->SetRangeUser(fFuncMin,fFuncMax);
  hist->Draw();
  fit->Draw("SAME");
  TLegend * leg = new TLegend(0.5,0.3,0.5,0.3);
  leg->SetHeader(hist->GetTitle(),"c");
  leg->AddEntry(hist,Form("%s",hist->GetName()));
  leg->Draw("SAME");
  
  if (iSide == 0) {
    cQA->Print(Form("%s/QA/SigmaFitQA_%s_Nearside.png",fOutputDir.Data(),hist->GetName()));
  } else {
    cQA->Print(Form("%s/QA/SigmaFitQA_%s_Awayside.png",fOutputDir.Data(),hist->GetName()));
  }
  
  return fit;
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

  if (bApplyEPRCorrection && !fIsMCGenMode) {
    fSystematicsNames.push_back("EPR Correction");
    OutOverIn_AS_EPRError = new TGraphErrors(nObsBins);
    OutOverIn_NS_EPRError = new TGraphErrors(nObsBins);
    OutOverIn_AS_EPRError->SetName("OutOverIn_AS_EPRError");
    OutOverIn_NS_EPRError->SetName("OutOverIn_NS_EPRError");


    OutOverIn_AS_SysErrorBySource.push_back(OutOverIn_AS_EPRError);
    OutOverIn_NS_SysErrorBySource.push_back(OutOverIn_NS_EPRError);
  }

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

  // FIXME this does not include all the RPFError graphs.

  vector<TGraphErrors * > gResultsGraphs = {
    fNSYieldsInc, fASYieldsInc, fNSRmsInc, fASRmsInc, fNSSigmasInc,fASSigmasInc,
    fNSYieldsInc_RPFError, fASYieldsInc_RPFError, fNSRmsInc_RPFError, fASRmsInc_RPFError,
    OutOverIn_AS, OutOverIn_NS, MidOverIn_AS, MidOverIn_NS,
    OutOverIn_AS_RPFError, OutOverIn_NS_RPFError, MidOverIn_AS_RPFError, MidOverIn_NS_RPFError,
    RmsOutOverIn_AS, RmsOutOverIn_NS, RmsMidOverIn_AS, RmsMidOverIn_NS,
    RmsOutOverIn_AS_RPFError, RmsOutOverIn_NS_RPFError, RmsMidOverIn_AS_RPFError, RmsMidOverIn_NS_RPFError,
    SigmasOutOverIn_AS, SigmasOutOverIn_NS, SigmasMidOverIn_AS, SigmasMidOverIn_NS,
    SigmasOutOverIn_AS_RPFError, SigmasOutOverIn_NS_RPFError, SigmasMidOverIn_AS_RPFError, SigmasMidOverIn_NS_RPFError,
    OutMinusIn_AS, OutMinusIn_NS, MidMinusIn_AS, MidMinusIn_NS
  };

  // Adding the EP yield, width arrays to this array.
  for (int iEP = 0; iEP < kNEPBins; iEP++) {
    gResultsGraphs.push_back(fNSYieldsEP[iEP]);
    gResultsGraphs.push_back(fASYieldsEP[iEP]);
    gResultsGraphs.push_back(fNSRmsEP[iEP]);
    gResultsGraphs.push_back(fASRmsEP[iEP]);
    gResultsGraphs.push_back(fNSSigmasEP[iEP]);
    gResultsGraphs.push_back(fASSigmasEP[iEP]);
    gResultsGraphs.push_back(fNSYieldsEP_RPFError[iEP]);
    gResultsGraphs.push_back(fASYieldsEP_RPFError[iEP]);
    gResultsGraphs.push_back(fNSRmsEP_RPFError[iEP]);
    gResultsGraphs.push_back(fASRmsEP_RPFError[iEP]);
  }

  gResultsGraphs.push_back(OutOverIn_NS_EPRError);
  gResultsGraphs.push_back(OutOverIn_AS_EPRError);

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

  if (iSigmaNSRangeBinChange != 0) {
    printf("Shifting bin of nearside yield range shift by sigma shift %d\n",iSigmaNSRangeBinChange);
    iMinRangeBin -= iSigmaNSRangeBinChange;
    iMaxRangeBin += iSigmaNSRangeBinChange;
  }

  fYieldNS = histNearSide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldNS_Err,"width");

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fYieldNS_Var = 0;
    Double_t fYieldNS_Var_Err = 0;

    // Integrate over the Y-axis bin iVar
    fYieldNS_Var = histNearSideVar->IntegralAndError(iMinRangeBin,iMaxRangeBin,iVar+1,iVar+1,fYieldNS_Var_Err,"width");


    if (bNormalizeByEtaRange) {
      fYieldNS_Var = fYieldNS_Var * 2. * fEtaRangeNearside;
      fYieldNS_Var_Err = fYieldNS_Var_Err * 2. * fEtaRangeNearside;
    }

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

  histNearSide->GetXaxis()->SetTitleOffset(0.5);

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


  printf("  Found Near-side Yield = %e \\pm %e and RMS = %f \\pm %f\n",fYieldNS,fYieldNS_Err,fRmsNS,fRmsNS_Err);


  double fNearSideSigma = 0;
  double fNearSideSigma_Err = 0;

  int iMinRangeSigmaNS = histNearSide->GetXaxis()->FindBin(-fSigmaRangeNS);
  int iMaxRangeSigmaNS = histNearSide->GetXaxis()->FindBin(fSigmaRangeNS);
  
  if (iSigmaNSRangeBinChange != 0) {
    printf("Shifting bin of nearside sigma shift by %d\n",iSigmaNSRangeBinChange);
    iMinRangeSigmaNS -= iSigmaNSRangeBinChange;
    iMaxRangeSigmaNS += iSigmaNSRangeBinChange;
  }


  double fActualMinRangeSigmaNS = histNearSide->GetXaxis()->GetBinLowEdge(iMinRangeSigmaNS);
  double fActualMaxRangeSigmaNS = histNearSide->GetXaxis()->GetBinUpEdge(iMaxRangeSigmaNS);

  printf("Debug: default sigma NS fit range: %f - %f\n",-fSigmaRangeNS,fSigmaRangeNS);
  printf("Debug: printing nearside sigma over bins %d - %d (%f %f)\n",iMinRangeSigmaNS,iMaxRangeSigmaNS,fActualMinRangeSigmaNS,fActualMaxRangeSigmaNS);

  //histNearSide->GetXaxis()->SetRangeUser(-fSigmaRangeNS,fSigmaRangeNS);
  histNearSide->GetXaxis()->SetRange(iMinRangeSigmaNS,iMaxRangeSigmaNS);

  TF1 * fitNearSidePeak = FitAndCalculateSigmaFromHisto(histNearSide,&fNearSideSigma,&fNearSideSigma_Err,0);

  // Draw a clone with a different color after setting the range
  histNearSide->SetMarkerColor(kSelectColor);
  histNearSide->SetLineColor(kSelectColor);
  histNearSide->Draw("SAME");
  printf("Sigma: Calculated nearside sigma = %f \\pm %f\n",fNearSideSigma,fNearSideSigma_Err);

//  bool bNormalizeByEtaRange = true;

//  float fEtaRangeNearside = 0.8;
//  float fEtaRangeAwayside =  1.35;
  if (bNormalizeByEtaRange) {
    fYieldNS = fYieldNS * 2. * fEtaRangeNearside;
    fYieldNS_Err = fYieldNS_Err * 2. * fEtaRangeNearside;
  }

  if (iEPBin == kNEPBins) {
    fNSYieldsInc->SetPoint(iObsBin,xValue,fYieldNS);
    fNSYieldsInc->SetPointError(iObsBin,xWidth,fYieldNS_Err);

    fNSRmsInc->SetPoint(iObsBin,xValue,fRmsNS);
    fNSRmsInc->SetPointError(iObsBin,xWidth,fRmsNS_Err);

    fNSSigmasInc->SetPoint(iObsBin,xValue,fNearSideSigma);
    fNSSigmasInc->SetPointError(iObsBin,xWidth,fNearSideSigma_Err);

  } else {
    fNSYieldsEP[iEPBin]->SetPoint(iObsBin,xValue,fYieldNS);
    fNSYieldsEP[iEPBin]->SetPointError(iObsBin,xWidth,fYieldNS_Err);

    fNSRmsEP[iEPBin]->SetPoint(iObsBin,xValue,fRmsNS);
    fNSRmsEP[iEPBin]->SetPointError(iObsBin,xWidth,fRmsNS_Err);

    fNSSigmasEP[iEPBin]->SetPoint(iObsBin,xValue,fNearSideSigma);
    fNSSigmasEP[iEPBin]->SetPointError(iObsBin,xWidth,fNearSideSigma_Err);
  }


  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fRmsNS_Var = 0;
    Double_t fRmsNS_Var_Err = 0;
    fRmsNS_Var = histNearSideRmsProfile->GetBinError(iVar+1);
    fRmsNS_Var_Err = 0;
    //fRmsNS_Var = histNearSideRmsProfile->GetBinContent(iVar+1);
    //fRmsNS_Var_Err = histNearSideRmsProfile->GetBinError(iVar+1);

    if (iEPBin == kNEPBins) {
      fNSRmsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fRmsNS_Var);
      fNSRmsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fRmsNS_Var_Err);
    } else {
      fNSRmsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fRmsNS_Var);
      fNSRmsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fRmsNS_Var_Err);
    }
  }
  canv->cd();
  TLegend * legNS = new TLegend(0.55,0.45,0.9,0.85);
  legNS->SetHeader(Form("%.1f #leq p_{T}^{a} < %.1f",fObsBins[iObsBin],fObsBins[iObsBin+1]),"c");
  legNS->AddEntry(histNearSide,"Near-side","pe");
  legNS->AddEntry((TObject *) 0,Form("Yield = %.2e #pm %.2e",fYieldNS,fYieldNS_Err),"");
  legNS->AddEntry((TObject *) 0,Form("RMS   = %.2e #pm %.2e",fRmsNS,fRmsNS_Err),"");
  legNS->AddEntry((TObject *) 0,Form("Sigma = %.2e #pm %.2e",fNearSideSigma,fNearSideSigma_Err),"");

  legNS->Draw("SAME");
  canv->SetGridx();
  canv->SetGridy();

  canv->Print(Form("%s/QA/CalcQA_ObsBin%d_EPBin%d_NearSide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/QA/CalcQA_ObsBin%d_EPBin%d_NearSide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_NearSide.C",fOutputDir.Data(),iObsBin,iEPBin));

  // Nearside Sigma variants
  printf("Calculating sigma NS for variants\n");

  if (bFitSigmaSlices) {
    TString sForm1Gaus = "[0]*TMath::Gaus(x,0,[1]) + [2]";
    double fFuncMin = -fSigmaRangeNS;
    double fFuncMax = fSigmaRangeNS;

    TF1 * fitNS = new TF1(Form("%s_Fit",histNearSide->GetName()),sForm1Gaus,fFuncMin,fFuncMax);

    double fMax = histNearSide->GetBinContent(histNearSide->GetMaximumBin());

    // Setting up the thing
    int iSigmaPar = 1;
    int iPedestalPar = 2;

    fitNS->SetParName(iSigmaPar,"sigma");

    fitNS->SetParameter(0,fMax);
    fitNS->SetParameter(iSigmaPar,0.1);
    fitNS->SetParameter(iPedestalPar,0);
    fitNS->SetParLimits(0,0.,1.5*fMax);
    fitNS->SetParLimits(iSigmaPar,0.001,1.0);
    fitNS->SetParLimits(iPedestalPar,0.,fMax);

    //if (!fIsMCGenMode) fitNS->FixParameter(iPedestalPar,0.);
    fitNS->FixParameter(iPedestalPar,0.);

    TObjArray SlicesFitReturn;
    histNearSideVar->GetXaxis()->SetRangeUser(-fSigmaRangeNS,fSigmaRangeNS);
    histNearSideVar->FitSlicesX(fitNS,0,-1,0,"QNR",&SlicesFitReturn);
    // scale, sigma, pedestal, chi2
    TH1F * hSlicesSigma = (TH1F *) SlicesFitReturn[1];
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {
      if (iEPBin == kNEPBins) {
        fNSSigmasInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,hSlicesSigma->GetBinContent(iVar));
        fNSSigmasInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,hSlicesSigma->GetBinError(iVar));
      } else {
        fNSSigmasEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,hSlicesSigma->GetBinContent(iVar));
        fNSSigmasEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,hSlicesSigma->GetBinError(iVar));
      }
    }

  } else {
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {
      double fSigmasNS_Var = 0;
      double fSigmasNS_Var_Err = 0;

      TH1D * histNearSideVarProj = (TH1D *) histNearSideVar->ProjectionX(Form("%s_Var%d",histNearSideVar->GetName(),iVar),iVar+1,iVar+1,"e");
      histNearSideVarProj->GetXaxis()->SetRangeUser(-fSigmaRangeNS,fSigmaRangeNS);
      TF1 * fitNearSidePeakVar = FitAndCalculateSigmaFromHisto(histNearSideVarProj,&fSigmasNS_Var,&fSigmasNS_Var_Err,0);

      if (iEPBin == kNEPBins) {
        fNSSigmasInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fSigmasNS_Var);
        fNSSigmasInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fSigmasNS_Var_Err);
      } else {
        fNSSigmasEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fSigmasNS_Var);
        fNSSigmasEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fSigmasNS_Var_Err);
      }


      delete histNearSideVarProj;
    }
  }


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


  if (iSigmaASRangeBinChange != 0) {
    printf("Shifting bin of awayside yield range shift by sigma shift %d\n",iSigmaNSRangeBinChange);
    iMinRangeBin -= iSigmaASRangeBinChange;
    iMaxRangeBin += iSigmaASRangeBinChange;
  }


  fYieldAS = histAwaySide->IntegralAndError(iMinRangeBin,iMaxRangeBin,fYieldAS_Err,"width");

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fYieldAS_Var = 0;
    Double_t fYieldAS_Var_Err = 0;
    //Double_t fRmsAS_Var = 0;
    //Double_t fRmsAS_Var_Err = 0;

    // Integrate over the Y-axis bin iVar
    fYieldAS_Var = histAwaySideVar->IntegralAndError(iMinRangeBin,iMaxRangeBin,iVar+1,iVar+1,fYieldAS_Var_Err,"width");
    // probably have to project a histogram to get the RMS

    if (bNormalizeByEtaRange) {
      fYieldAS_Var = fYieldAS_Var * 2. * fEtaRangeAwayside;
      fYieldAS_Var_Err = fYieldAS_Var_Err * 2. * fEtaRangeAwayside;
    }
    printf("Debug AS Yield (EP=%d,iObsBin=%d,iVar=%d) = %.2e #pm %.2e\n",iEPBin,iObsBin,iVar,fYieldAS_Var,fYieldAS_Var_Err);
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

  // FIXME new
  histAwaySide->GetXaxis()->SetTitleOffset(0.5);

  histAwaySide->DrawCopy();

  histAwaySide->GetXaxis()->SetRangeUser(TMath::Pi()-fRmsRangeAS,TMath::Pi()+fRmsRangeAS);
  fRmsAS     = histAwaySide->GetRMS();
  fRmsAS_Err = histAwaySide->GetRMSError();



  printf("  Found Away-side Yield = %e \\pm %e and RMS = %f \\pm %f\n",fYieldNS,fYieldNS_Err,fRmsNS,fRmsNS_Err);


  double fAwaySideSigma = 0;
  double fAwaySideSigma_Err = 0;



  int iMinRangeSigmaAS = histAwaySide->GetXaxis()->FindBin(-fSigmaRangeAS+TMath::Pi());
  int iMaxRangeSigmaAS = histAwaySide->GetXaxis()->FindBin(fSigmaRangeAS+TMath::Pi());


  if (iSigmaASRangeBinChange != 0) {
    printf("Shifting bin of awayside sigma shift by %d\n",iSigmaASRangeBinChange);
    iMinRangeSigmaAS -= iSigmaASRangeBinChange;
    iMaxRangeSigmaAS += iSigmaASRangeBinChange;
  }

  double fActualMinRangeSigmaAS = histAwaySide->GetXaxis()->GetBinLowEdge(iMinRangeSigmaAS);
  double fActualMaxRangeSigmaAS = histAwaySide->GetXaxis()->GetBinUpEdge(iMaxRangeSigmaAS);

  printf("Debug: default sigma AS fit range: %f - %f\n",-fSigmaRangeAS+TMath::Pi(),fSigmaRangeAS+TMath::Pi());
  printf("Debug: printing awayside sigma over bins %d - %d (%f %f)\n",iMinRangeSigmaAS,iMaxRangeSigmaAS,fActualMinRangeSigmaAS,fActualMaxRangeSigmaAS);



  //histAwaySide->GetXaxis()->SetRangeUser(TMath::Pi()-fSigmaRangeAS,TMath::Pi()+fSigmaRangeAS);
  histAwaySide->GetXaxis()->SetRange(iMinRangeSigmaAS,iMaxRangeSigmaAS);

  TF1 * fitAwaySidePeak = FitAndCalculateSigmaFromHisto(histAwaySide,&fAwaySideSigma,&fAwaySideSigma_Err,1);


  histAwaySide->SetMarkerColor(kSelectColor);
  histAwaySide->SetFillColor(kSelectColor);
  histAwaySide->SetLineColor(kSelectColor);
  histAwaySide->Draw("SAME F");

  printf("Sigma: Calculated Awayside sigma = %f \\pm %f\n",fAwaySideSigma,fAwaySideSigma_Err);


  // Scale by the eta range to get 1/n dn/dphi


  if (bNormalizeByEtaRange) {
    fYieldAS = fYieldAS * 2. * fEtaRangeAwayside;
    fYieldAS_Err = fYieldAS_Err * 2. * fEtaRangeAwayside;
  }


  if (iEPBin == kNEPBins) {
    fASYieldsInc->SetPoint(iObsBin,xValue,fYieldAS);
    fASYieldsInc->SetPointError(iObsBin,xWidth,fYieldAS_Err);

    fASRmsInc->SetPoint(iObsBin,xValue,fRmsAS);
    fASRmsInc->SetPointError(iObsBin,xWidth,fRmsAS_Err);

    fASSigmasInc->SetPoint(iObsBin,xValue,fAwaySideSigma);
    fASSigmasInc->SetPointError(iObsBin,xWidth,fAwaySideSigma_Err);
  } else {
    fASYieldsEP[iEPBin]->SetPoint(iObsBin,xValue,fYieldAS);
    fASYieldsEP[iEPBin]->SetPointError(iObsBin,xWidth,fYieldAS_Err);

    fASRmsEP[iEPBin]->SetPoint(iObsBin,xValue,fRmsAS);
    fASRmsEP[iEPBin]->SetPointError(iObsBin,xWidth,fRmsAS_Err);

    fASSigmasEP[iEPBin]->SetPoint(iObsBin,xValue,fAwaySideSigma);
    fASSigmasEP[iEPBin]->SetPointError(iObsBin,xWidth,fAwaySideSigma_Err);
  }



  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    Double_t fRmsAS_Var = 0;
    Double_t fRmsAS_Var_Err = 0;
    fRmsAS_Var = histAwaySideRmsProfile->GetBinError(iVar+1);

    fRmsAS_Var_Err = 0;
    //fRmsAS_Var = histAwaySideRmsProfile->GetBinContent(iVar+1);
    //fRmsAS_Var_Err = histAwaySideRmsProfile->GetBinError(iVar+1);

    if (iEPBin == kNEPBins) {
      fASRmsInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fRmsAS_Var);
      fASRmsInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fRmsAS_Var_Err);
    } else {
      fASRmsEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fRmsAS_Var);
      fASRmsEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fRmsAS_Var_Err);
    }
  }


  canv->cd();
  TLegend * legAS = new TLegend(0.55,0.45,0.9,0.85);
  legAS->SetHeader(Form("%.1f #leq p_{T}^{a} < %.1f",fObsBins[iObsBin],fObsBins[iObsBin+1]),"c");
  legAS->AddEntry(histNearSide,"Away-side","pe");
  legAS->AddEntry((TObject *) 0,Form("Yield = %.2e #pm %.2e",fYieldAS,fYieldAS_Err),"");
  legAS->AddEntry((TObject *) 0,Form("RMS   = %.2e #pm %.2e",fRmsAS,fRmsAS_Err),"");
  legAS->AddEntry((TObject *) 0,Form("Sigma = %.2e #pm %.2e",fAwaySideSigma,fAwaySideSigma_Err),"");

  legAS->Draw("SAME");

  canv->Print(Form("%s/QA/CalcQA_ObsBin%d_EPBin%d_AwaySide.pdf",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/QA/CalcQA_ObsBin%d_EPBin%d_AwaySide.png",fOutputDir.Data(),iObsBin,iEPBin));
  canv->Print(Form("%s/CFiles/CalcQA_ObsBin%d_EPBin%d_AwaySide.C",fOutputDir.Data(),iObsBin,iEPBin));



  printf("Calculating sigma AS for variants\n");

  if (bFitSigmaSlices) {

    TString sForm1Gaus = "[0]*TMath::Gaus(x,TMath::Pi(),[1]) + [2]";
    //double fFuncMin = TMath::Pi()/2;
    //double fFuncMax = 3.*TMath::Pi()/2;
    //double fFuncMin = TMath::Pi() - fSigmaRangeAS;
    //double fFuncMax = TMath::Pi() + fSigmaRangeAS;

    double fFuncMin = fActualMinRangeSigmaAS;
    double fFuncMax = fActualMaxRangeSigmaAS;

    TF1 * fitAS = new TF1(Form("%s_Fit",histAwaySide->GetName()),sForm1Gaus,fFuncMin,fFuncMax);

    double fMax = histAwaySide->GetBinContent(histAwaySide->GetMaximumBin());

    // Setting up the thing
    int iSigmaPar = 1;
    int iPedestalPar = 2;

    fitAS->SetParName(iSigmaPar,"sigma");

    fitAS->SetParameter(0,fMax);
    fitAS->SetParameter(iSigmaPar,0.1);
    fitAS->SetParameter(iPedestalPar,0);

    fitAS->SetParLimits(0,0.,1.5*fMax);
    fitAS->SetParLimits(iSigmaPar,0.001,1.0);
    fitAS->SetParLimits(iPedestalPar,0.,fMax);

    //if (!fIsMCGenMode) fitAS->FixParameter(iPedestalPar,0.);
    fitAS->FixParameter(iPedestalPar,0.);

    TObjArray SlicesFitReturn;
    histAwaySideVar->GetXaxis()->SetRangeUser(TMath::Pi()-fSigmaRangeAS,TMath::Pi()+fSigmaRangeAS);
    histAwaySideVar->FitSlicesX(fitAS,0,-1,0,"NR",&SlicesFitReturn);
    // scale, sigma, pedestal, chi2
    TH1F * hSlicesSigma = (TH1F *) SlicesFitReturn[1];
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {
      if (iEPBin == kNEPBins) {
        fASSigmasInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,hSlicesSigma->GetBinContent(iVar));
        fASSigmasInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,hSlicesSigma->GetBinError(iVar));
      } else {
        fASSigmasEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,hSlicesSigma->GetBinContent(iVar));
        fASSigmasEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,hSlicesSigma->GetBinError(iVar));
      }
    }


  } else {
    for (int iVar = 0; iVar < nRPFVariants; iVar++) {
      double fSigmasAS_Var = 0;
      double fSigmasAS_Var_Err = 0;

      TH1D * histAwaySideVarProj = (TH1D *) histAwaySideVar->ProjectionX(Form("%s_Var%d",histAwaySideVar->GetName(),iVar),iVar+1,iVar+1,"e");
      histAwaySideVarProj->GetXaxis()->SetRangeUser(TMath::Pi()-fSigmaRangeAS,TMath::Pi()+fSigmaRangeAS);
      TF1 * fitAwaySidePeakVar = FitAndCalculateSigmaFromHisto(histAwaySideVarProj,&fSigmasAS_Var,&fSigmasAS_Var_Err,1);
      

      if (iEPBin == kNEPBins) {
        fASSigmasInc_RPFVariants[iVar]->SetPoint(iObsBin,xValue,fSigmasAS_Var);
        fASSigmasInc_RPFVariants[iVar]->SetPointError(iObsBin,xValue,fSigmasAS_Var_Err);
      } else {
        fASSigmasEP_RPFVariants[iVar][iEPBin]->SetPoint(iObsBin,xValue,fSigmasAS_Var);
        fASSigmasEP_RPFVariants[iVar][iEPBin]->SetPointError(iObsBin,xValue,fSigmasAS_Var_Err);
      }


      delete histAwaySideVarProj;
    }
  }



  delete histNearSideRmsProfile;
  delete histAwaySideRmsProfile;

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
  //OutOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;%s / In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;Yield Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
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
  MidOverIn_AS->SetTitle(Form("%s/In Yield Ratio (%s);%s;Yield Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
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
  OutOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;Yield Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutOverIn_NS->SetLineColor(kOutInColorNS);
  OutOverIn_NS->SetMarkerColor(kOutInColorNS);
  OutOverIn_NS->SetMarkerStyle(kOutOverInMarker);
  OutOverIn_NS->SetFillColor(kOutInColor);
  //OutOverIn_NS->SetFillColorAlpha(kOutInColor,kOutInErrorAlpha);
  OutOverIn_NS_RPFError = (TGraphErrors *) OutOverIn_NS->Clone(Form("%s_RPFError",OutOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) OutOverIn_NS->Clone(Form("%s_RPFVar%d",OutOverIn_NS->GetName(),iVar));
    OutOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  MidOverIn_NS->SetName(Form("%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidOverIn_NS->SetTitle(Form("%s/In Yield Ratio (%s);%s;Yield Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  MidOverIn_NS->SetLineColor(kMidInColor);
  MidOverIn_NS->SetMarkerColor(kMidInColor);
  MidOverIn_NS->SetMarkerStyle(kMidOverInMarker);
  MidOverIn_NS->SetFillColor(kMidInColor);
  //MidOverIn_NS->SetFillColorAlpha(kMidInColor,kMidInErrorAlpha);
  MidOverIn_NS_RPFError = (TGraphErrors *) MidOverIn_NS->Clone(Form("%s_RPFError",MidOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) MidOverIn_NS->Clone(Form("%s_RPFVar%d",MidOverIn_NS->GetName(),iVar));
    MidOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  RmsOutOverIn_AS = new TGraphErrors(nObsBins);
  RmsOutOverIn_NS = new TGraphErrors(nObsBins);
  RmsMidOverIn_AS = new TGraphErrors(nObsBins);
  RmsMidOverIn_NS = new TGraphErrors(nObsBins);


  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  RmsOutOverIn_AS->SetLineColor(kOutInColor);
  RmsOutOverIn_AS->SetMarkerColor(kOutInColor);
  RmsOutOverIn_AS->SetMarkerStyle(kOutOverInMarker);
  RmsOutOverIn_AS->SetFillColor(kOutInColor);
  //RmsOutOverIn_AS->SetFillColorAlpha(kOutInColor,kOutInErrorAlpha);
  RmsOutOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RmsOutOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;RMS Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsOutOverIn_AS_RPFError = (TGraphErrors *) RmsOutOverIn_AS->Clone(Form("%s_RPFError",RmsOutOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RmsOutOverIn_AS->Clone(Form("%s_RPFVar%d",RmsOutOverIn_AS->GetName(),iVar));
    RmsOutOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  RmsMidOverIn_AS->SetLineColor(kMidInColor);
  RmsMidOverIn_AS->SetMarkerColor(kMidInColor);
  RmsMidOverIn_AS->SetMarkerStyle(kMidOverInMarker);
  RmsMidOverIn_NS->SetFillColor(kMidInColor);
  //RmsMidOverIn_NS->SetFillColorAlpha(kMidInColor,kMidInErrorAlpha);
  RmsMidOverIn_AS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RmsMidOverIn_AS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;RMS Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsMidOverIn_AS_RPFError = (TGraphErrors *) RmsMidOverIn_AS->Clone(Form("%s_RPFError",RmsMidOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RmsMidOverIn_AS->Clone(Form("%s_RPFVar%d",RmsMidOverIn_AS->GetName(),iVar));
    RmsMidOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  RmsOutOverIn_NS->SetLineColor(kOutInColorNS);
  RmsOutOverIn_NS->SetMarkerColor(kOutInColorNS);
  RmsOutOverIn_NS->SetMarkerStyle(kOutOverInMarker);
  RmsOutOverIn_NS->SetFillColor(kOutInColorNS);
  //RmsOutOverIn_NS->SetFillColorAlpha(kOutInColorNS,kOutInErrorAlpha);
  RmsOutOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RmsOutOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;RMS Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsOutOverIn_NS_RPFError = (TGraphErrors *) RmsOutOverIn_NS->Clone(Form("%s_RPFError",RmsOutOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RmsOutOverIn_NS->Clone(Form("%s_RPFVar%d",RmsOutOverIn_NS->GetName(),iVar));
    RmsOutOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  RmsMidOverIn_NS->SetLineColor(kMidInColor);
  RmsMidOverIn_NS->SetMarkerColor(kMidInColor);
  RmsMidOverIn_NS->SetMarkerStyle(kMidOverInMarker);
  RmsMidOverIn_NS->SetFillColor(kMidInColor);
  //RmsMidOverIn_NS->SetFillColorAlpha(kMidInColor,kMidInErrorAlpha);
  RmsMidOverIn_NS->SetName(Form("RMS%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  RmsMidOverIn_NS->SetTitle(Form("RMS%s/In Width Ratio (%s);%s;RMS Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  RmsMidOverIn_NS_RPFError = (TGraphErrors *) RmsMidOverIn_NS->Clone(Form("%s_RPFError",RmsMidOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) RmsMidOverIn_NS->Clone(Form("%s_RPFVar%d",RmsMidOverIn_NS->GetName(),iVar));
    RmsMidOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  // Initialzing Sigma Ratios

  SigmasOutOverIn_AS = new TGraphErrors(nObsBins);
  SigmasOutOverIn_NS = new TGraphErrors(nObsBins);
  SigmasMidOverIn_AS = new TGraphErrors(nObsBins);
  SigmasMidOverIn_NS = new TGraphErrors(nObsBins);



  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  SigmasOutOverIn_AS->SetLineColor(kOutInColor);
  SigmasOutOverIn_AS->SetMarkerColor(kOutInColor);
  SigmasOutOverIn_AS->SetMarkerStyle(kOutOverInMarker);
  //SigmasOutOverIn_AS->SetFillColorAlpha(kOutInColor,kOutInErrorAlpha);
  SigmasOutOverIn_AS->SetFillColor(kOutInColor);
  SigmasOutOverIn_AS->SetName(Form("Sigmas%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  SigmasOutOverIn_AS->SetTitle(Form("Sigmas%s/In Width Ratio (%s);%s;Width Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  SigmasOutOverIn_AS_RPFError = (TGraphErrors *) SigmasOutOverIn_AS->Clone(Form("%s_RPFError",SigmasOutOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) SigmasOutOverIn_AS->Clone(Form("%s_RPFVar%d",SigmasOutOverIn_AS->GetName(),iVar));
    SigmasOutOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  SigmasMidOverIn_AS->SetLineColor(kMidInColor);
  SigmasMidOverIn_AS->SetMarkerColor(kMidInColor);
  SigmasMidOverIn_AS->SetMarkerStyle(kMidOverInMarker);
  SigmasMidOverIn_NS->SetFillColor(kMidInColor);
  //SigmasMidOverIn_NS->SetFillColorAlpha(kMidInColor,kMidInErrorAlpha);
  SigmasMidOverIn_AS->SetName(Form("Sigmas%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  SigmasMidOverIn_AS->SetTitle(Form("Sigmas%s/In Width Ratio (%s);%s;Width Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  SigmasMidOverIn_AS_RPFError = (TGraphErrors *) SigmasMidOverIn_AS->Clone(Form("%s_RPFError",SigmasMidOverIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) SigmasMidOverIn_AS->Clone(Form("%s_RPFVar%d",SigmasMidOverIn_AS->GetName(),iVar));
    SigmasMidOverIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  SigmasOutOverIn_NS->SetLineColor(kOutInColorNS);
  SigmasOutOverIn_NS->SetMarkerColor(kOutInColorNS);
  SigmasOutOverIn_NS->SetMarkerStyle(kOutOverInMarker);
  SigmasOutOverIn_NS->SetFillColor(kOutInColorNS);
  //SigmasOutOverIn_NS->SetFillColorAlpha(kOutInColorNS,kOutInErrorAlpha);
  SigmasOutOverIn_NS->SetName(Form("Sigmas%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  SigmasOutOverIn_NS->SetTitle(Form("Sigmas%s/In Width Ratio (%s);%s;Width Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  SigmasOutOverIn_NS_RPFError = (TGraphErrors *) SigmasOutOverIn_NS->Clone(Form("%s_RPFError",SigmasOutOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) SigmasOutOverIn_NS->Clone(Form("%s_RPFVar%d",SigmasOutOverIn_NS->GetName(),iVar));
    SigmasOutOverIn_NS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  SigmasMidOverIn_NS->SetLineColor(kMidInColor);
  SigmasMidOverIn_NS->SetMarkerColor(kMidInColor);
  SigmasMidOverIn_NS->SetMarkerStyle(kMidOverInMarker);
  SigmasMidOverIn_NS->SetFillColorAlpha(kMidInColor,kMidInErrorAlpha);
  SigmasMidOverIn_NS->SetName(Form("Sigmas%sOverIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  SigmasMidOverIn_NS->SetTitle(Form("Sigmas%s/In Width Ratio (%s);%s;Width Ratio (%s/In)",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  SigmasMidOverIn_NS_RPFError = (TGraphErrors *) SigmasMidOverIn_NS->Clone(Form("%s_RPFError",SigmasMidOverIn_NS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) SigmasMidOverIn_NS->Clone(Form("%s_RPFVar%d",SigmasMidOverIn_NS->GetName(),iVar));
    SigmasMidOverIn_NS_RPFVariants.push_back(VariationGraph);
  }



  // Initializing Yield Differences


  OutMinusIn_AS = new TGraphErrors(nObsBins);
  OutMinusIn_NS = new TGraphErrors(nObsBins);
  MidMinusIn_AS = new TGraphErrors(nObsBins);
  MidMinusIn_NS = new TGraphErrors(nObsBins);


  iSide = 0; // 1 for AS
  iNumerator = 0; // 0 for Mid
  OutMinusIn_AS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutMinusIn_AS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s - In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));
  OutMinusIn_AS_RPFError = (TGraphErrors *) OutMinusIn_AS->Clone(Form("%s_RPFError",OutMinusIn_AS->GetName()));
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    TGraphErrors * VariationGraph = 0;
    VariationGraph = (TGraphErrors *) OutMinusIn_AS->Clone(Form("%s_RPFVar%d",OutMinusIn_AS->GetName(),iVar));
    OutMinusIn_AS_RPFVariants.push_back(VariationGraph);
  }

  iNumerator = 1; // 1 for Mid
  MidMinusIn_AS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidMinusIn_AS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s - In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iSide = 1; // 1 for NS
  iNumerator = 0; // 0 for Mid
  OutMinusIn_NS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  OutMinusIn_NS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s - In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

  iNumerator = 1; // 1 for Mid
  MidMinusIn_NS->SetName(Form("%sMinusIn_%s",sNumerator[iNumerator].Data(),sSideName[iSide].Data()));
  MidMinusIn_NS->SetTitle(Form("%s - In Yield Difference (%s);%s;%s - In",sNumerator[iNumerator].Data(),sSideTitle[iSide].Data(),fObservableName.Data(),sNumerator[iNumerator].Data()));

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

  double fASEPRYieldOutOverInErr = 0; 

  if (YieldInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 yield in-plane\n");
    return;  
  }
  double YieldOutOverIn = YieldOutPlane / YieldInPlane;

  //double EPRRes2 = 0.77; // FIXME
  double EPRRes2 = fEPRes[1];
  double EPRDefaultError = 0.04;

  if (fIsMCGenMode) EPRDefaultError = 0;


  if (bApplyEPRCorrection && !fIsMCGenMode) {
    //printf("Appling EPR Correction with EPR = %f for away-side\n",EPRRes2);

    double EPRRes2Low = EPRRes2 - EPRDefaultError;
    double EPRRes2High = EPRRes2 + EPRDefaultError;
    double YieldOutOverInLow = ((1-1./EPRRes2Low) + (1+1./EPRRes2Low)*YieldOutOverIn) / ( (1+1./EPRRes2Low) + (1-1./EPRRes2Low)*YieldOutOverIn);
    double YieldOutOverInHigh = ((1-1./EPRRes2High) + (1+1./EPRRes2High)*YieldOutOverIn) / ( (1+1./EPRRes2High) + (1-1./EPRRes2High)*YieldOutOverIn);

    YieldOutOverIn = ((1-1./EPRRes2) + (1+1./EPRRes2)*YieldOutOverIn) / ( (1+1./EPRRes2) + (1-1./EPRRes2)*YieldOutOverIn);

    printf("Calculated AS YieldOutOverIn = %f with (low,high) = (%f,%f), AbsHalfDiff = %f\n",YieldOutOverIn,YieldOutOverInLow,YieldOutOverInHigh, 0.5*TMath::Abs(YieldOutOverInHigh-YieldOutOverInLow));

    fASEPRYieldOutOverInErr = 0.5*TMath::Abs(YieldOutOverInHigh-YieldOutOverInLow);
  }



  double YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));



  double YieldMidOverIn = YieldMidPlane / YieldInPlane;
  double YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  OutOverIn_AS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  OutOverIn_AS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  if (bApplyEPRCorrection && !fIsMCGenMode) {
    OutOverIn_AS_EPRError->SetPoint(iObsBin,xValue,YieldOutOverIn);
    OutOverIn_AS_EPRError->SetPointError(iObsBin,xErr,fASEPRYieldOutOverInErr);
  }

  MidOverIn_AS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  MidOverIn_AS->SetPointError(iObsBin,xErr,YieldMidOverInError);

  printf("Debug calculating yield and yield ratio RPF variations (iObsBin = %d)\n",iObsBin);

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {

    int nVarYields = fASYieldsEP_RPFVariants.size();
    if (nRPFVariants != nVarYields) {
      printf("Warning  nRPFVariants != nVarYields\n");
    }
    int nVarEventPlanes = fASYieldsEP_RPFVariants[iVar].size();
    if (nVarEventPlanes != 3) {
      printf("Warning nVarEventPlanes != 3 \n");
    }

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

    if (bApplyEPRCorrection && !fIsMCGenMode) {
      YieldOutOverIn = ((1-1./EPRRes2) + (1+1./EPRRes2)*YieldOutOverIn) / ( (1+1./EPRRes2) + (1-1./EPRRes2)*YieldOutOverIn);

    }


    YieldOutOverInError = TMath::Abs(YieldOutOverIn) * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = TMath::Abs(YieldMidOverIn) * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    if (YieldOutOverIn < 0 || YieldOutOverIn > 2) {
      printf("Outlier:   ");
      printf("%.2f #pm %.2f\n",YieldOutOverIn,YieldOutOverInError);

      printf("      YieldOut = %.2e #pm %.2e,  YieldIn = %.2e #pm %.2e \n",YieldOutPlane,YieldOutPlaneError,YieldInPlane,YieldInPlaneError);
    } else {
      printf("%.2f #pm %.2f\n",YieldOutOverIn,YieldOutOverInError);
    }

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


  double fNSEPRYieldOutOverInErr = 0; 

  if (bApplyEPRCorrection && !fIsMCGenMode) {


    //printf("Appling EPR Correction with EPR = %f for near-side\n",EPRRes2);

    double EPRRes2Low = EPRRes2 - EPRDefaultError;
    double EPRRes2High = EPRRes2 + EPRDefaultError;
    double YieldOutOverInLow = ((1-1./EPRRes2Low) + (1+1./EPRRes2Low)*YieldOutOverIn) / ( (1+1./EPRRes2Low) + (1-1./EPRRes2Low)*YieldOutOverIn);
    double YieldOutOverInHigh = ((1-1./EPRRes2High) + (1+1./EPRRes2High)*YieldOutOverIn) / ( (1+1./EPRRes2High) + (1-1./EPRRes2High)*YieldOutOverIn);


    YieldOutOverIn = ((1-1./EPRRes2) + (1+1./EPRRes2)*YieldOutOverIn) / ( (1+1./EPRRes2) + (1-1./EPRRes2)*YieldOutOverIn);


    printf("Calculated NS YieldOutOverIn = %f with (low,high) = (%f,%f), AbsHalfDiff = %f\n",YieldOutOverIn,YieldOutOverInLow,YieldOutOverInHigh, 0.5*TMath::Abs(YieldOutOverInHigh-YieldOutOverInLow));

    fNSEPRYieldOutOverInErr = 0.5*TMath::Abs(YieldOutOverInHigh-YieldOutOverInLow);

  }


  YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  YieldMidOverIn = YieldMidPlane / YieldInPlane;
  YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

  OutOverIn_NS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  OutOverIn_NS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  if (bApplyEPRCorrection && !fIsMCGenMode) {
    OutOverIn_NS_EPRError->SetPoint(iObsBin,xValue,YieldOutOverIn);
    OutOverIn_NS_EPRError->SetPointError(iObsBin,xErr,fNSEPRYieldOutOverInErr);
  }

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

    if (bApplyEPRCorrection && !fIsMCGenMode) {
      YieldOutOverIn = ((1-1./EPRRes2) + (1+1./EPRRes2)*YieldOutOverIn) / ( (1+1./EPRRes2) + (1-1./EPRRes2)*YieldOutOverIn);
    }

    YieldOutOverInError = YieldOutOverIn * TMath::Sqrt(TMath::Power(YieldOutPlaneError/YieldOutPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    YieldMidOverIn = YieldMidPlane / YieldInPlane;
    YieldMidOverInError = YieldMidOverIn * TMath::Sqrt(TMath::Power(YieldMidPlaneError/YieldMidPlane,2) + TMath::Power(YieldInPlaneError/YieldInPlane,2));

    OutOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldOutOverIn);
    OutOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldOutOverInError);

    MidOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,YieldMidOverIn);
    MidOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,YieldMidOverInError);

  }



  // RMS 
  // ========================================================
  // Away Side
  // ========================================================

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

  RmsOutOverIn_AS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  RmsOutOverIn_AS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  RmsMidOverIn_AS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  RmsMidOverIn_AS->SetPointError(iObsBin,xErr,YieldMidOverInError);

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





  // RMS 
  // ========================================================
  // Near Side
  // ========================================================
  // Too lazy to change the variable name from yield to width

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


  RmsOutOverIn_NS->SetPoint(iObsBin,xValue,YieldOutOverIn);
  RmsOutOverIn_NS->SetPointError(iObsBin,xErr,YieldOutOverInError);

  RmsMidOverIn_NS->SetPoint(iObsBin,xValue,YieldMidOverIn);
  RmsMidOverIn_NS->SetPointError(iObsBin,xErr,YieldMidOverInError);

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

  printf("Starting ratio calculation with sigmas\n");

  // Sigma 
  // ========================================================
  // Away Side
  // ========================================================

  double SigmaInPlane      = fASSigmasEP[0]->GetY()[iObsBin];
  double SigmaInPlaneError = fASSigmasEP[0]->GetEY()[iObsBin];

  double SigmaMidPlane      = fASSigmasEP[1]->GetY()[iObsBin];
  double SigmaMidPlaneError = fASSigmasEP[1]->GetEY()[iObsBin];

  double SigmaOutPlane      = fASSigmasEP[2]->GetY()[iObsBin];
  double SigmaOutPlaneError = fASSigmasEP[2]->GetEY()[iObsBin];

  //printf("  Debug: sigmas %f %f %f\n",SigmaInPlane,SigmaMidPlane,SigmaOutPlane);

  if (SigmaInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 width in-plane\n");
    return;  
  }

  double SigmaOutOverIn = SigmaOutPlane / SigmaInPlane;
  double SigmaOutOverInError = SigmaOutOverIn * TMath::Sqrt(TMath::Power(SigmaOutPlaneError/SigmaOutPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));

  double SigmaMidOverIn = SigmaMidPlane / SigmaInPlane;
  double SigmaMidOverInError = SigmaMidOverIn * TMath::Sqrt(TMath::Power(SigmaMidPlaneError/SigmaMidPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));

  SigmasOutOverIn_AS->SetPoint(iObsBin,xValue,SigmaOutOverIn);
  SigmasOutOverIn_AS->SetPointError(iObsBin,xErr,SigmaOutOverInError);

  SigmasMidOverIn_AS->SetPoint(iObsBin,xValue,SigmaMidOverIn);
  SigmasMidOverIn_AS->SetPointError(iObsBin,xErr,SigmaMidOverInError);


  // Calculation with RPF variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    SigmaInPlane      = fASSigmasEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    SigmaInPlaneError = fASSigmasEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    SigmaMidPlane      = fASSigmasEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    SigmaMidPlaneError = fASSigmasEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    SigmaOutPlane      = fASSigmasEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    SigmaOutPlaneError = fASSigmasEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (SigmaInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 AS Sigma in-plane (iObsBin = %d\n",iObsBin);
      continue;  
    }
    SigmaOutOverIn = SigmaOutPlane / SigmaInPlane;
    SigmaOutOverInError = SigmaOutOverIn * TMath::Sqrt(TMath::Power(SigmaOutPlaneError/SigmaOutPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));

    SigmaMidOverIn = SigmaMidPlane / SigmaInPlane;
    SigmaMidOverInError = SigmaMidOverIn * TMath::Sqrt(TMath::Power(SigmaMidPlaneError/SigmaMidPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));
    
    SigmasOutOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,SigmaOutOverIn);
    SigmasOutOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,SigmaOutOverInError);

    SigmasMidOverIn_AS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,SigmaMidOverIn);
    SigmasMidOverIn_AS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,SigmaMidOverInError);
  }





  // Sigma 
  // ========================================================
  // Near Side
  // ========================================================


  SigmaInPlane      = fNSSigmasEP[0]->GetY()[iObsBin];
  SigmaInPlaneError = fNSSigmasEP[0]->GetEY()[iObsBin];

  SigmaMidPlane      = fNSSigmasEP[1]->GetY()[iObsBin];
  SigmaMidPlaneError = fNSSigmasEP[1]->GetEY()[iObsBin];

  SigmaOutPlane      = fNSSigmasEP[2]->GetY()[iObsBin];
  SigmaOutPlaneError = fNSSigmasEP[2]->GetEY()[iObsBin];

  if (SigmaInPlane==0) {
    fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 sigma width in-plane (iObsBin = %d\n",iObsBin);
    return;  
  }

  SigmaOutOverIn = SigmaOutPlane / SigmaInPlane;
  SigmaOutOverInError = SigmaOutOverIn * TMath::Sqrt(TMath::Power(SigmaOutPlaneError/SigmaOutPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));

  SigmaMidOverIn = SigmaMidPlane / SigmaInPlane;
  SigmaMidOverInError = SigmaMidOverIn * TMath::Sqrt(TMath::Power(SigmaMidPlaneError/SigmaMidPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));


  SigmasOutOverIn_NS->SetPoint(iObsBin,xValue,SigmaOutOverIn);
  SigmasOutOverIn_NS->SetPointError(iObsBin,xErr,SigmaOutOverInError);

  SigmasMidOverIn_NS->SetPoint(iObsBin,xValue,SigmaMidOverIn);
  SigmasMidOverIn_NS->SetPointError(iObsBin,xErr,SigmaMidOverInError);

  // Calculation with RPF variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    SigmaInPlane      = fNSSigmasEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    SigmaInPlaneError = fNSSigmasEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    SigmaMidPlane      = fNSSigmasEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    SigmaMidPlaneError = fNSSigmasEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    SigmaOutPlane      = fNSSigmasEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    SigmaOutPlaneError = fNSSigmasEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    if (SigmaInPlane==0) {
      fprintf(stderr,"Error: Narrowly avoided division by zero with a 0 NS Sigma in-plane (iObsBin = %d)\n",iObsBin);
      continue;  
    }
    SigmaOutOverIn = SigmaOutPlane / SigmaInPlane;
    SigmaOutOverInError = SigmaOutOverIn * TMath::Sqrt(TMath::Power(SigmaOutPlaneError/SigmaOutPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));

    SigmaMidOverIn = SigmaMidPlane / SigmaInPlane;
    SigmaMidOverInError = SigmaMidOverIn * TMath::Sqrt(TMath::Power(SigmaMidPlaneError/SigmaMidPlane,2) + TMath::Power(SigmaInPlaneError/SigmaInPlane,2));
    
    SigmasOutOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,SigmaOutOverIn);
    SigmasOutOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,SigmaOutOverInError);

    SigmasMidOverIn_NS_RPFVariants[iVar]->SetPoint(iObsBin,xValue,SigmaMidOverIn);
    SigmasMidOverIn_NS_RPFVariants[iVar]->SetPointError(iObsBin,xErr,SigmaMidOverInError);

  }


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

  fNSSigmasInc_RPFError = (TGraphErrors *) fNSSigmasInc->Clone(Form("%s_RPF_Error",fNSSigmasInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSSigmasEP_RPFError.push_back((TGraphErrors *) fNSSigmasEP[iEPBin]->Clone(Form("%s_RPF_Error",fNSSigmasEP[iEPBin]->GetName())));
  }

  fASSigmasInc_RPFError = (TGraphErrors *) fASSigmasInc->Clone(Form("%s_RPF_Error",fASSigmasInc->GetName()));
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASSigmasEP_RPFError.push_back((TGraphErrors *) fASSigmasEP[iEPBin]->Clone(Form("%s_RPF_Error",fASSigmasEP[iEPBin]->GetName())));
  }

  TGraphErrors * fNSYieldsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fNSYieldsEP_RPFErr_Mean = {};
  TGraphErrors * fASYieldsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fASYieldsEP_RPFErr_Mean = {};

  TGraphErrors * fNSRmsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fNSRmsEP_RPFErr_Mean = {};
  TGraphErrors * fASRmsInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fASRmsEP_RPFErr_Mean = {};

  TGraphErrors * fNSSigmasInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fNSSigmasEP_RPFErr_Mean = {};
  TGraphErrors * fASSigmasInc_RPFErr_Mean = 0;
  vector<TGraphErrors *> fASSigmasEP_RPFErr_Mean = {};

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

  fNSSigmasInc_RPFErr_Mean = ProduceSystematicFromGraphs(fNSSigmasInc_RPFVariants,fNSSigmasInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fNSSigmasEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fNSSigmasEP_RPFVariants[iEPBin],fNSSigmasEP_RPFError[iEPBin]));
  }
  fASSigmasInc_RPFErr_Mean = ProduceSystematicFromGraphs(fASSigmasInc_RPFVariants,fASSigmasInc_RPFError);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    fASSigmasEP_RPFErr_Mean.push_back(ProduceSystematicFromGraphs(fASSigmasEP_RPFVariants[iEPBin],fASSigmasEP_RPFError[iEPBin]));
  }


  printf("Finished calculating RPF fit error for first calculated variables.\n");

  fNSYieldsInc_RPFErr_Mean->SetLineColor(kBlack);
  


}

/**
  * Calculate the errors in the ratios from the RPF variations
  */
void TaskCalcObservables::CalculateRPFErrorRatios() {

  cout<<"Calculating RPF Error for Ratios"<<endl;

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
    OutOverIn_AS_RPFError->SetPointError(i,OutOverIn_AS->GetEX()[i],0.0);
    MidOverIn_AS_RPFError->SetPointError(i,MidOverIn_AS->GetEX()[i],0.0);
  }
  for (int i = 0; i < OutOverIn_NS->GetN(); i++) {
    OutOverIn_NS_RPFError->SetPoint(i,OutOverIn_NS->GetX()[i],OutOverIn_NS->GetY()[i]);
    MidOverIn_NS_RPFError->SetPoint(i,MidOverIn_NS->GetX()[i],MidOverIn_NS->GetY()[i]);
    OutOverIn_NS_RPFError->SetPointError(i,OutOverIn_NS->GetEX()[i],0.0);
    MidOverIn_NS_RPFError->SetPointError(i,MidOverIn_NS->GetEX()[i],0.0);
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

  // Settings for systematic uncertainty. E4
  //OutOverIn_AS_RPFError->SetFillColor(kOutInErrorColor);
  //OutOverIn_AS_RPFError->SetFillColorAlpha(kOutInRPFErrorColor,kOutInRPFErrorAlpha);
  OutOverIn_AS_RPFError->SetFillColor(kOutInRPFErrorColor);
  //OutOverIn_AS_RPFError->SetFillStyle(kOutInFillStyle);
  OutOverIn_AS_RPFError->SetMarkerColor(kBlue-9);

  OutOverIn_AS_RPFError->SetMarkerStyle(kOpenSquare);
  OutOverIn_AS_RPFError->SetLineColor(kBlack);
 // OutOverIn_AS_RPFError->SetFillColor(kBlack);
  OutOverIn_AS_RPFError->SetMarkerColor(kBlack);
//  OutOverIn_AS_RPFError->SetDrawOption("E4");
  
  //MidOverIn_AS_RPFError->SetFillColorAlpha(kMidInRPFErrorColor,kMidInRPFErrorAlpha);
  MidOverIn_AS_RPFError->SetFillColor(kMidInRPFErrorColor);
  //MidOverIn_AS_RPFError->SetFillStyle(kMidInFillStyle);
  MidOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
  MidOverIn_AS_RPFError->SetMarkerStyle(kOpenSquare);
  MidOverIn_AS_RPFError->SetLineColor(kBlack);
  //MidOverIn_AS_RPFError->SetFillColor(kGray);
//  MidOverIn_AS_RPFError->SetDrawOption("E4");


  //OutOverIn_NS_RPFError->SetFillColorAlpha(kOutInRPFErrorColor,kOutInRPFErrorAlpha);
  OutOverIn_NS_RPFError->SetFillColor(kOutInRPFErrorColor);
  OutOverIn_NS_RPFError->SetFillStyle(kOutInNSFillStyle);
  OutOverIn_NS_RPFError->SetMarkerColor(kBlue-9);
//  OutOverIn_AS_RPFError->SetDrawOption("E4");
  
  MidOverIn_NS_RPFError->SetFillColor(kMidInRPFErrorColor);
  MidOverIn_NS_RPFError->SetFillStyle(kMidInNSFillStyle);
  MidOverIn_NS_RPFError->SetMarkerColor(kBlue-9);



  TGraphErrors * RmsOutOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * RmsMidOverIn_AS_RPFErr_Mean = 0;
  TGraphErrors * RmsOutOverIn_NS_RPFErr_Mean = 0;
  TGraphErrors * RmsMidOverIn_NS_RPFErr_Mean = 0;

  for (int i = 0; i < RmsOutOverIn_AS->GetN(); i++) {
    RmsOutOverIn_AS_RPFError->SetPoint(i,RmsOutOverIn_AS->GetX()[i],RmsOutOverIn_AS->GetY()[i]);
    RmsMidOverIn_AS_RPFError->SetPoint(i,RmsMidOverIn_AS->GetX()[i],RmsMidOverIn_AS->GetY()[i]);

    RmsOutOverIn_AS_RPFError->SetPointError(i,RmsOutOverIn_AS->GetEX()[i],0.0);
    RmsMidOverIn_AS_RPFError->SetPointError(i,RmsMidOverIn_AS->GetEX()[i],0.0);
  }
  for (int i = 0; i < RmsOutOverIn_NS->GetN(); i++) {
    RmsOutOverIn_NS_RPFError->SetPoint(i,RmsOutOverIn_NS->GetX()[i],RmsOutOverIn_NS->GetY()[i]);
    RmsMidOverIn_NS_RPFError->SetPoint(i,RmsMidOverIn_NS->GetX()[i],RmsMidOverIn_NS->GetY()[i]);

    RmsOutOverIn_NS_RPFError->SetPointError(i,RmsOutOverIn_NS->GetEX()[i],0.0);
    RmsMidOverIn_NS_RPFError->SetPointError(i,RmsMidOverIn_NS->GetEX()[i],0.0);
  }

  printf("Calculating RPF Systematic for RMS OutOverIn_AS\n");
  RmsOutOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsOutOverIn_AS_RPFVariants,RmsOutOverIn_AS_RPFError);
  printf("Calculating RPF Systematic for RMS MidOverIn_AS\n");
  RmsMidOverIn_AS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsMidOverIn_AS_RPFVariants,RmsMidOverIn_AS_RPFError);

  printf("Calculating RPF Systematic for RMS OutOverIn_NS\n");
  RmsOutOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsOutOverIn_NS_RPFVariants,RmsOutOverIn_NS_RPFError);
  printf("Calculating RPF Systematic for RMS MidOverIn_NS\n");
  RmsMidOverIn_NS_RPFErr_Mean = ProduceSystematicFromGraphs(RmsMidOverIn_NS_RPFVariants,RmsMidOverIn_NS_RPFError);


  RmsOutOverIn_AS_RPFError->SetFillColorAlpha(kOutInRPFErrorColor,kOutInRPFErrorAlpha);
  RmsOutOverIn_AS_RPFError->SetFillStyle(kOutInFillStyle);
  RmsOutOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
  RmsOutOverIn_AS_RPFError->SetLineColor(kBlack);
  
  RmsMidOverIn_AS_RPFError->SetFillColorAlpha(kMidInRPFErrorColor,kMidInRPFErrorAlpha);
  RmsMidOverIn_AS_RPFError->SetFillStyle(kMidInFillStyle);
  RmsMidOverIn_AS_RPFError->SetMarkerColor(kBlue-9);
  RmsMidOverIn_AS_RPFError->SetLineColor(kBlack);

  RmsOutOverIn_NS_RPFError->SetFillColorAlpha(kOutInRPFErrorColor,kOutInRPFErrorAlpha);
  RmsOutOverIn_NS_RPFError->SetFillStyle(kOutInNSFillStyle);
  RmsOutOverIn_NS_RPFError->SetMarkerColor(kBlue-9);
  RmsOutOverIn_NS_RPFError->SetLineColor(kBlack);
  
  RmsMidOverIn_NS_RPFError->SetFillColorAlpha(kMidInRPFErrorColor,kMidInRPFErrorAlpha);
  RmsMidOverIn_NS_RPFError->SetFillStyle(kMidInNSFillStyle);
  RmsMidOverIn_NS_RPFError->SetMarkerColor(kBlue-9);
  RmsMidOverIn_NS_RPFError->SetLineColor(kBlack);


}


// Calculate Differences
void TaskCalcObservables::CalculateResultsDifferencesObsBin(int iObsBin, TCanvas * canv) {

//  TGraphErrors * OutMinusIn_AS;
//  TGraphErrors * OutMinusIn_NS;
//  TGraphErrors * MidMinusIn_AS;
//  TGraphErrors * MidMinusIn_NS;

  double xValue = fASYieldsEP[0]->GetX()[iObsBin];
  double xErr   = fASYieldsEP[0]->GetEX()[iObsBin];

  double YieldInPlane      = fASYieldsEP[0]->GetY()[iObsBin];
  double YieldInPlaneError = fASYieldsEP[0]->GetEY()[iObsBin];

  double YieldMidPlane      = fASYieldsEP[1]->GetY()[iObsBin];
  double YieldMidPlaneError = fASYieldsEP[1]->GetEY()[iObsBin];

  double YieldOutPlane      = fASYieldsEP[2]->GetY()[iObsBin];
  double YieldOutPlaneError = fASYieldsEP[2]->GetEY()[iObsBin];


  double YieldOutMinusIn = YieldOutPlane - YieldInPlane;
  double YieldOutMinusInError = TMath::Sqrt(TMath::Power(YieldOutPlaneError,2) + TMath::Power(YieldInPlaneError,2));
  double YieldMidMinusIn = YieldMidPlane - YieldInPlane;
  double YieldMidMinusInError = TMath::Sqrt(TMath::Power(YieldMidPlaneError,2) + TMath::Power(YieldInPlaneError,2));


  OutMinusIn_AS->SetPoint(iObsBin,xValue,YieldOutMinusIn);
  OutMinusIn_AS->SetPointError(iObsBin,xErr,YieldOutMinusInError);
  MidMinusIn_AS->SetPoint(iObsBin,xValue,YieldMidMinusIn);
  MidMinusIn_AS->SetPointError(iObsBin,xErr,YieldMidMinusInError);

  // FIXME calculate this for OutMinusIn
  //if (bApplyEPRCorrection && !fIsMCGenMode) {
  //  OutOverIn_AS_EPRError->SetPoint(iObsBin,xValue,YieldOutOverIn);
  //  OutOverIn_AS_EPRError->SetPointError(iObsBin,xErr,fASEPRYieldOutOverInErr);
 // }

  // Calculation with the RPF Variants
  for (int iVar = 0; iVar < nRPFVariants; iVar++) {
    int nVarYields = fASYieldsEP_RPFVariants.size();
    if (nRPFVariants != nVarYields) {
      printf("Warning  nRPFVariants != nVarYields\n");
    }
    int nVarEventPlanes = fASYieldsEP_RPFVariants[iVar].size();
    if (nVarEventPlanes != 3) {
      printf("Warning nVarEventPlanes != 3 \n");
    }

    YieldInPlane      = fASYieldsEP_RPFVariants[iVar][0]->GetY()[iObsBin];
    YieldInPlaneError = fASYieldsEP_RPFVariants[iVar][0]->GetEY()[iObsBin];

    YieldMidPlane      = fASYieldsEP_RPFVariants[iVar][1]->GetY()[iObsBin];
    YieldMidPlaneError = fASYieldsEP_RPFVariants[iVar][1]->GetEY()[iObsBin];

    YieldOutPlane      = fASYieldsEP_RPFVariants[iVar][2]->GetY()[iObsBin];
    YieldOutPlaneError = fASYieldsEP_RPFVariants[iVar][2]->GetEY()[iObsBin];

    YieldOutMinusIn = YieldOutPlane - YieldInPlane;
    YieldOutMinusInError = TMath::Sqrt(TMath::Power(YieldOutPlaneError,2) + TMath::Power(YieldInPlaneError,2));
    YieldMidMinusIn = YieldMidPlane - YieldInPlane;
    YieldMidMinusInError = TMath::Sqrt(TMath::Power(YieldMidPlaneError,2) + TMath::Power(YieldInPlaneError,2));


  }

}

// Set the graph styles
void TaskCalcObservables::SetGraphStyles() {

  TString sYieldTitle = "Yields";
  sYieldTitle = "#frac{1}{N_{#pi^{0}}#frac{dN^{assoc}}{d#Delta#varphi}";

  // All EP together
  fNSYieldsInc->SetLineColor(kEPColorList[3]);
  fNSYieldsInc->SetMarkerColor(kEPColorList[3]);
  fNSYieldsInc->SetMarkerStyle(kEPMarkerList[3]);
  fNSYieldsInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fNSYieldsInc->GetYaxis()->SetTitle("#frac{1}{N_{#pi^{0}}}#frac{dN^{assoc}}{d#Delta#varphi}");

  fNSRmsInc->SetLineColor(kEPColorList[3]);
  fNSRmsInc->SetMarkerColor(kEPColorList[3]);
  fNSRmsInc->SetMarkerStyle(kEPMarkerList[3]);
  fNSRmsInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fNSRmsInc->GetYaxis()->SetTitle("RMS");

  fNSSigmasInc->SetLineColor(kEPColorList[3]);
  fNSSigmasInc->SetMarkerColor(kEPColorList[3]);
  fNSSigmasInc->SetMarkerStyle(kEPMarkerList[3]);
  fNSSigmasInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fNSSigmasInc->GetYaxis()->SetTitle("#sigma_{NS}");

  fASYieldsInc->SetLineColor(kEPColorList[3]);
  fASYieldsInc->SetMarkerColor(kEPColorList[3]);
  fASYieldsInc->SetMarkerStyle(kEPMarkerList[3]);
  fASYieldsInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fASYieldsInc->GetYaxis()->SetTitle("#frac{1}{N_{#pi^{0}}}#frac{dN^{assoc}}{d#Delta#varphi}");

  fASRmsInc->SetLineColor(kEPColorList[3]);
  fASRmsInc->SetMarkerColor(kEPColorList[3]);
  fASRmsInc->SetMarkerStyle(kEPMarkerList[3]);
  fASRmsInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fASRmsInc->GetYaxis()->SetTitle("RMS");

  fASSigmasInc->SetLineColor(kEPColorList[3]);
  fASSigmasInc->SetMarkerColor(kEPColorList[3]);
  fASSigmasInc->SetMarkerStyle(kEPMarkerList[3]);
  fASSigmasInc->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
  fASSigmasInc->GetYaxis()->SetTitle("#sigma_{AS}");
  // By EP Bin
  for (int i = 0; i < kNEPBins; i++) {
    fNSYieldsEP[i]->SetLineColor(kEPColorList[i]);
    fNSYieldsEP[i]->SetMarkerColor(kEPColorList[i]);
    fNSYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fNSYieldsEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fNSYieldsEP[i]->GetYaxis()->SetTitle("#frac{1}{N_{#pi^{0}}}#frac{dN^{assoc}}{d#Delta#varphi}");

    fNSRmsEP[i]->SetLineColor(kEPColorList[i]);
    fNSRmsEP[i]->SetMarkerColor(kEPColorList[i]);
    fNSRmsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fNSRmsEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fNSRmsEP[i]->GetYaxis()->SetTitle("RMS");

    fNSSigmasEP[i]->SetLineColor(kEPColorList[i]);
    fNSSigmasEP[i]->SetMarkerColor(kEPColorList[i]);
    fNSSigmasEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fNSSigmasEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fNSSigmasEP[i]->GetYaxis()->SetTitle("#sigma_{NS}");

    fASYieldsEP[i]->SetLineColor(kEPColorList[i]);
    fASYieldsEP[i]->SetMarkerColor(kEPColorList[i]);
    fASYieldsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fASYieldsEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fASYieldsEP[i]->GetYaxis()->SetTitle("#frac{1}{N_{#pi^{0}}}#frac{dN^{assoc}}{d#Delta#varphi}");

    fASRmsEP[i]->SetLineColor(kEPColorList[i]);
    fASRmsEP[i]->SetMarkerColor(kEPColorList[i]);
    fASRmsEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fASRmsEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fASRmsEP[i]->GetYaxis()->SetTitle("RMS");

    fASSigmasEP[i]->SetLineColor(kEPColorList[i]);
    fASSigmasEP[i]->SetMarkerColor(kEPColorList[i]);
    fASSigmasEP[i]->SetMarkerStyle(kEPMarkerList[i]);
    fASSigmasEP[i]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/#it{c})");
    fASSigmasEP[i]->GetYaxis()->SetTitle("#sigma_{AS}");

  }

}


// take array of {inc,in,mid,out} ?
// and array of RPFuncert {inc,in,mid,out} 
// array of systematic errors.
// or should it be array of arrays of systematic errors( one per inc,in,etc and one per systematic error type (v3,purity,etc)
  // Make it all{in,mid,out,incl} to match up with arrays elsewhere
  // For fModels: {in,mid,out,incl}, then {model1,model2,...}
void TaskCalcObservables::DrawObservable(vector<TGraphErrors *> fObsGraphs, vector<TGraphErrors *> fObsRPFErrors, vector<TGraphErrors *> fObsSysErrors, vector<vector<TGraphErrors *>> fModels) {

  printf("Drawing observable %s (%s)\n",fObsGraphs[0]->GetName(),fObsGraphs[0]->GetTitle());

  bool bHasSysErrors = (fObsSysErrors.size() > 0);

  bool bHasModels = (fModels.size() > 0);

  TCanvas * cObservable = new TCanvas("cObservable","cObservable",fYieldCanvasWidth,fYieldCanvasHeight);

  TLegend * legend = new TLegend(0.55,0.65,0.85,0.85);

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
    legend->Clear();

    TString sTitleShort = fEPBinTitlesShort[i];

    TMultiGraph * mg = new TMultiGraph();


    TGraphErrors * fObsGraph = fObsGraphs[i];
    fObsGraph->SetMarkerStyle(kEPMarkerList[i]);
    fObsGraph->SetMarkerColor(kEPColorList[i]);
    fObsGraph->SetLineColor(kEPColorList[i]);

    RemoveXErrorBars(fObsGraph);

    TGraphErrors * fObsRPFError = fObsRPFErrors[i];
    fObsRPFError->SetMarkerStyle(kEPMarkerList[i]);
    fObsRPFError->SetMarkerColor(kEPColorList[i]);
    //fObsRPFError->SetLineColor(kEPColorList[i]);
    fObsRPFError->SetFillColor(kEPRPFFillColorList[i]);
    fObsRPFError->SetFillStyle(kEPRPFFillStyleList[i]);

    TGraphErrors * fObsSysError = 0;


    fObsGraph->Draw("AP");
    //fObsRPFError->Draw("E3 SAME");
    //fObsRPFError->Draw("[]5 SAME");
    fObsRPFError->Draw("5 SAME");

    legend->AddEntry(fObsGraph,Form("%s (Stat. Error)",fObsGraph->GetTitle()),"LEP");
    legend->AddEntry(fObsRPFError,"RPF Error","F");

    TLegend * legAlice = DrawAliceLegend(fObsGraph,0,0,kAliceLegendWidth,kAliceLegendHeight);
    //TLegend * legAlice = DrawAliceLegend(fObsGraph,0.55,0.35,0.2,0.15);


    if (bHasSysErrors) {
      printf("Intend to code something here to plot the other systematic errors\n");

      printf("  i = %d, fObsSysErrors.size() = %d\n",i,(int) fObsSysErrors.size());

      fObsSysError = fObsSysErrors[i];
      fObsSysError->SetLineColor(kEPColorList[i]);
      fObsSysError->SetMarkerStyle(kEPMarkerList[i]);
      fObsSysError->SetFillColor(kEPSysFillColorList[i]);
      fObsSysError->SetFillStyle(kEPSysFillStyleList[i]);
 


      //fObsSysError->Draw("E3 SAME");
      //fObsSysError->Draw("5 SAME");
      fObsSysError->Draw("3 SAME");
      //fObsSysError->Draw("[]5 SAME");

      legend->AddEntry(fObsSysError,"Sys. Error","F");

    }

    legend->Draw("SAME");

    cObservable->Print(Form("%s/%sWithErrors.pdf",fOutputDir.Data(),fObsGraph->GetName()));
    cObservable->Print(Form("%s/%sWithErrors.png",fOutputDir.Data(),fObsGraph->GetName()));

    if (bHasModels) {

      for (int j = 0; j < (int) fModels[i].size(); j++) {
        TGraphErrors * fModel = fModels[i][j];
        printf("Drawing model comparison %d %d with object %s\n",i,j,fModel->GetName());
        fModel->SetMarkerStyle(kModelStyles[j]);
        //fModel->SetLineColor(kModelColorsEP[i][j]);
        //fModel->SetFillColorAlpha(kModelColorsEP[i][j],0.5);
        //fModel->SetMarkerColor(kModelColorsEP[i][j]);
        fModel->SetLineColor(kModelColors[j]);
        fModel->SetFillColor(kModelColors[j]);
        //fModel->SetFillStyle(kModelFillStyles[j]);
        fModel->SetFillStyle(1);
        //fModel->SetFillColorAlpha(kModelColors[j],0.5);
        fModel->SetMarkerColor(kModelColors[j]);

        fModel->Draw("SAME 3 L");
        legend->AddEntry(fModel,fModelTitles[j].Data(),"CFL");
        //legend->AddEntry(fModel,fModelTitles[j].Data(),"CLP4");
      }

      legAlice->Draw("SAME");
      cObservable->Print(Form("%s/%sWithModels.pdf",fOutputDir.Data(),fObsGraph->GetName()));
      cObservable->Print(Form("%s/%sWithModels.png",fOutputDir.Data(),fObsGraph->GetName()));
    }

  }

}

///
/// Draw an ALICE performance legend entry
//
//________________________________________________________________________
TLegend * TaskCalcObservables::DrawAliceLegend(TObject *obj, Float_t x, Float_t y, Float_t x_size, Float_t y_size)
{
  const char *kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
  const char *kCentList[5] = {"0-90%","0-10%","10-30%","30-50%","50-90%"}; // index=fCent+1
  TLegend * leg = 0;
  if (x == 0 && y == 0)
  // Automatic placement:
  leg =  new TLegend(x_size,y_size);
  else 
  // Manual placement:
  leg  = new TLegend(x,y,x+x_size,y+y_size);

  TDatime * time = new TDatime();
  const char * month = kMonthList[time->GetMonth()-1];

  //leg->SetHeader(Form("ALICE Performance - %d %s %d",time->GetDay(),month,time->GetYear()));
  //if (fPerformance) leg->AddEntry(Histo,Form("ALICE Performance %d %s %d",time->GetDay()-1,month,time->GetYear()),"");  
  //else leg->AddEntry(Histo,Form("Work in Progress %d %s %d",time->GetDay()-1,month,time->GetYear()),""); 

  // would do case switch, but have seen wierd stuff from root before
  if (iPlotStatus == 0) leg->AddEntry(obj,Form("Work in Progress %d %s %d",time->GetDay()-1,month,time->GetYear()),"");
  else if (iPlotStatus == 1) leg->AddEntry(obj,"ALICE Preliminary","");
  else if (iPlotStatus == 2) leg->AddEntry(obj,"This Note","");
  else if (iPlotStatus == 3) leg->AddEntry(obj,"This Thesis","");

  //if (fPreliminary) leg->AddEntry(obj,"ALICE Preliminary","");
  //if (fPreliminary) leg->AddEntry(obj,"This Thesis","");

  //else leg->AddEntry(obj,"This thesis",""); 

  else leg->AddEntry(obj,Form("Work in Progress %d %s %d",time->GetDay()-1,month,time->GetYear()),"");
 // leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 0-90%","");
  //leg->AddEntry(obj,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, %s",kCentList[iCentBin+1]),"");
  //leg->AddEntry(obj,"#pi^{0}-hadron Correlations","");
  leg->AddEntry(obj,Form("Pb-Pb #rightarrow #pi^{0} + h + X"),"");
  leg->AddEntry(obj,Form(" #sqrt{#it{s}_{NN}} = 5.02 TeV, %s",kCentList[iCentBin+1]),"");

  if (!fIsMCGenMode) {
    if (sTitle != "") {
      leg->AddEntry(obj,Form("%s %.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}",sTitle.Data(),PtBins[iPtBin-1],PtBins[iPtBin]),"");
    }
    else {
      leg->AddEntry(obj,Form("%.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}",PtBins[iPtBin-1],PtBins[iPtBin]),"");
    }
  }

  leg->SetTextSize(0.041); // 0.045
  leg->SetBorderSize(0);
  //leg->SetFillColorAlpha(10,0);
  leg->SetFillStyle(0);

  leg->Draw("SAME");

  return leg;
}




void TaskCalcObservables::DrawResults() {
  printf("Saving some results plots\n");

  SetGraphStyles();

  float fRatioMin = 0.5;
  float fRatioMax = 1.5;


  // Drawing Individual Observable Plots
  vector<TGraphErrors *> fNSYieldsGraphs = {};
  vector<TGraphErrors *> fNSYieldsRPFErrors = {};

  vector<TGraphErrors *> fNSRmsGraphs = {};
  vector<TGraphErrors *> fNSRmsRPFErrors = {};


  vector<TGraphErrors *> fASYieldsGraphs = {};
  vector<TGraphErrors *> fASYieldsRPFErrors = {};
  vector<TGraphErrors *> fASRmsGraphs = {};
  vector<TGraphErrors *> fASRmsRPFErrors = {};

  vector<TGraphErrors *> fNSSigmasGraphs = {};
  vector<TGraphErrors *> fNSSigmasRPFErrors = {};
  vector<TGraphErrors *> fASSigmasGraphs = {};
  vector<TGraphErrors *> fASSigmasRPFErrors = {};

  printf("hello\n");

  for (int i = 0; i < 4; i++) {
    if (i == 3) {
      fNSYieldsGraphs.push_back(fNSYieldsInc);
      fNSYieldsRPFErrors.push_back(fNSYieldsInc_RPFError);

      fNSRmsGraphs.push_back(fNSRmsInc);
      fNSRmsRPFErrors.push_back(fNSRmsInc_RPFError);

      fNSSigmasGraphs.push_back(fNSSigmasInc);
      fNSSigmasRPFErrors.push_back(fNSSigmasInc_RPFError);

      fASYieldsGraphs.push_back(fASYieldsInc);
      fASYieldsRPFErrors.push_back(fASYieldsInc_RPFError);

      fASRmsGraphs.push_back(fASRmsInc);
      fASRmsRPFErrors.push_back(fASRmsInc_RPFError);   

      fASSigmasGraphs.push_back(fASSigmasInc);
      fASSigmasRPFErrors.push_back(fASSigmasInc_RPFError);   
    } else {
      fNSYieldsGraphs.push_back(fNSYieldsEP[i]);
      fNSYieldsRPFErrors.push_back(fNSYieldsEP_RPFError[i]);

      fNSRmsGraphs.push_back(fNSRmsEP[i]);
      fNSRmsRPFErrors.push_back(fNSRmsEP_RPFError[i]);

      fNSSigmasGraphs.push_back(fNSSigmasEP[i]);
      fNSSigmasRPFErrors.push_back(fNSSigmasEP_RPFError[i]);

      fASYieldsGraphs.push_back(fASYieldsEP[i]);
      fASYieldsRPFErrors.push_back(fASYieldsEP_RPFError[i]);

      fASRmsGraphs.push_back(fASRmsEP[i]);
      fASRmsRPFErrors.push_back(fASRmsEP_RPFError[i]);

      fASSigmasGraphs.push_back(fASSigmasEP[i]);
      fASSigmasRPFErrors.push_back(fASSigmasEP_RPFError[i]);
    }
  }

  printf("hello 2\n");

  // Systematics
  vector<TGraphErrors *> fNSYieldsSysErrors = {};
  vector<TGraphErrors *> fASYieldsSysErrors = {};
  vector<TGraphErrors *> fNSRmsSysErrors = {};
  vector<TGraphErrors *> fASRmsSysErrors = {};
  vector<TGraphErrors *> fNSSigmasSysErrors = {};
  vector<TGraphErrors *> fASSigmasSysErrors = {};


  if (fNSYieldsInc_SysError != 0) {
    for (int i = 0; i < 4; i++) {
      if (i == 3) {
        fNSYieldsSysErrors.push_back(fNSYieldsInc_SysError);
        fASYieldsSysErrors.push_back(fASYieldsInc_SysError);
        fNSRmsSysErrors.push_back(fNSRmsInc_SysError);
        fASRmsSysErrors.push_back(fASRmsInc_SysError);
        fNSSigmasSysErrors.push_back(fNSSigmasInc_SysError);
        fASSigmasSysErrors.push_back(fASSigmasInc_SysError);

      } else {
        fNSYieldsSysErrors.push_back(fNSYieldsEP_SysError[i]);
        fASYieldsSysErrors.push_back(fASYieldsEP_SysError[i]);
        fNSRmsSysErrors.push_back(fNSRmsEP_SysError[i]);
        fASRmsSysErrors.push_back(fASRmsEP_SysError[i]);
        fNSSigmasSysErrors.push_back(fNSSigmasEP_SysError[i]);
        fASSigmasSysErrors.push_back(fASSigmasEP_SysError[i]);

      }

    }
  }

  printf("hello 3\n");

  // Models
  // First index is event plane (in,mid,out,incl), second is model
  vector<vector<TGraphErrors *>> fNSYieldsModels = {};
  vector<vector<TGraphErrors *>> fNSRmsModels = {};
  vector<vector<TGraphErrors *>> fNSSigmasModels = {};
  vector<vector<TGraphErrors *>> fASYieldsModels = {};
  vector<vector<TGraphErrors *>> fASRmsModels = {};
  vector<vector<TGraphErrors *>> fASSigmasModels = {};

  for (int i = 0; i < kNEPBins+1; i++) {
    vector<TGraphErrors *> fNSYieldsModelsWithinEPBin = {};
    vector<TGraphErrors *> fNSRmsModelsWithinEPBin = {};
    vector<TGraphErrors *> fNSSigmasModelsWithinEPBin = {};
    vector<TGraphErrors *> fASYieldsModelsWithinEPBin = {};
    vector<TGraphErrors *> fASRmsModelsWithinEPBin = {};
    vector<TGraphErrors *> fASSigmasModelsWithinEPBin = {};

    if (fNSYieldsInc_Models.size() > 0) { // Check if models are loaded

      if (i < kNEPBins) {
        for (int j = 0; j < (int) fNSYieldsEP_Models[i].size(); j++) {
          fNSYieldsModelsWithinEPBin.push_back(fNSYieldsEP_Models[i][j]);
        }
        for (int j = 0; j < (int) fNSRmsEP_Models[i].size(); j++) {
          fNSRmsModelsWithinEPBin.push_back(fNSRmsEP_Models[i][j]);
        }
        for (int j = 0; j < (int) fNSSigmasEP_Models[i].size(); j++) {
          fNSSigmasModelsWithinEPBin.push_back(fNSSigmasEP_Models[i][j]);
        }
        for (int j = 0; j < (int) fASYieldsEP_Models[i].size(); j++) {
          fASYieldsModelsWithinEPBin.push_back(fASYieldsEP_Models[i][j]);
        }
        for (int j = 0; j < (int) fASRmsEP_Models[i].size(); j++) {
          fASRmsModelsWithinEPBin.push_back(fASRmsEP_Models[i][j]);
        }
        for (int j = 0; j < (int) fASSigmasEP_Models[i].size(); j++) {
          fASSigmasModelsWithinEPBin.push_back(fASSigmasEP_Models[i][j]);
        }
      } else { // inclusive
        for (int j = 0; j < (int) fNSYieldsInc_Models.size(); j++) {
          fNSYieldsModelsWithinEPBin.push_back(fNSYieldsInc_Models[j]);
        }
        for (int j = 0; j < (int) fNSRmsInc_Models.size(); j++) {
          fNSRmsModelsWithinEPBin.push_back(fNSRmsInc_Models[j]);
        }
        for (int j = 0; j < (int) fNSSigmasInc_Models.size(); j++) {
          fNSSigmasModelsWithinEPBin.push_back(fNSSigmasInc_Models[j]);
        }
        for (int j = 0; j < (int) fASYieldsInc_Models.size(); j++) {
          fASYieldsModelsWithinEPBin.push_back(fASYieldsInc_Models[j]);
        }
        for (int j = 0; j < (int) fASRmsInc_Models.size(); j++) {
          fASRmsModelsWithinEPBin.push_back(fASRmsInc_Models[j]);
        }
        for (int j = 0; j < (int) fASSigmasInc_Models.size(); j++) {
          fASSigmasModelsWithinEPBin.push_back(fASSigmasInc_Models[j]);
        }
      }
    }

    fNSYieldsModels.push_back(fNSYieldsModelsWithinEPBin);
    fNSRmsModels.push_back(fNSRmsModelsWithinEPBin);
    fNSSigmasModels.push_back(fNSSigmasModelsWithinEPBin);
    fASYieldsModels.push_back(fASYieldsModelsWithinEPBin);
    fASRmsModels.push_back(fASRmsModelsWithinEPBin);
    fASSigmasModels.push_back(fASSigmasModelsWithinEPBin);
  }

  printf("hello finished organizing the things\n");


//  DrawObservable(fNSYieldsGraphs,fNSYieldsRPFErrors,fNSYieldsSysErrors);

  //DrawObservable(fNSRmsGraphs,fNSRmsRPFErrors,fNSRmsSysErrors);
  //DrawObservable(fASYieldsGraphs,fASYieldsRPFErrors,fASYieldsSysErrors);
  //DrawObservable(fASRmsGraphs,fASRmsRPFErrors,fASRmsSysErrors);
  //DrawObservable(fNSYieldsGraphs,fNSYieldsRPFErrors);
  //DrawObservable(fNSRmsGraphs,fNSRmsRPFErrors);
  //DrawObservable(fASYieldsGraphs,fASYieldsRPFErrors);
  //DrawObservable(fASRmsGraphs,fASRmsRPFErrors);

  DrawObservable(fNSYieldsGraphs,fNSYieldsRPFErrors,fNSYieldsSysErrors,fNSYieldsModels);
  DrawObservable(fNSRmsGraphs,fNSRmsRPFErrors,fNSRmsSysErrors,fNSRmsModels);
  DrawObservable(fASYieldsGraphs,fASYieldsRPFErrors,fASYieldsSysErrors,fASYieldsModels);
  DrawObservable(fASRmsGraphs,fASRmsRPFErrors,fASRmsSysErrors,fASRmsModels);

  // FIXME crash is here
  DrawObservable(fNSSigmasGraphs,fNSSigmasRPFErrors,fNSSigmasSysErrors,fNSSigmasModels);
  DrawObservable(fASSigmasGraphs,fASSigmasRPFErrors,fASSigmasSysErrors,fASSigmasModels);


  // Undecided about whether to have the same canvas sizes for yields and widths apparently
  TCanvas * cResults = new TCanvas("cResults","cResults",fYieldCanvasWidth,fYieldCanvasHeight);

  TMultiGraph * mgNSY = new TMultiGraph();
  TMultiGraph * mgNSRms = new TMultiGraph();
  TMultiGraph * mgNSSigmas = new TMultiGraph();
  TMultiGraph * mgASY = new TMultiGraph();
  TMultiGraph * mgASRms = new TMultiGraph();
  TMultiGraph * mgASSigmas = new TMultiGraph();

  TString sDrawStyle = "AP";
  TLegend * legYields = new TLegend(0.60,0.55,0.85,0.85);
  TLegend * legRms = new TLegend(0.55,0.65,0.85,0.85);
  TLegend * legSigmas = new TLegend(0.55,0.65,0.85,0.85);

  TString sDrawLabel = Form("%.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}",PtBins[iPtBin-1],PtBins[iPtBin]);

  legYields->SetHeader(sDrawLabel,"c");


  cResults->cd();
  // Near-side Yields
  mgNSY->Add(fNSYieldsInc);
  legYields->AddEntry(fNSYieldsInc,"Inclusive","LEP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legYields->AddEntry(fNSYieldsEP[iEPBin],fEPBinTitles[iEPBin],"LEP");
    mgNSY->Add(fNSYieldsEP[iEPBin]);
  }
  cResults->SetLogy(1);
  mgNSY->Draw(sDrawStyle);
  mgNSY->GetXaxis()->SetTitle(fNSYieldsInc->GetXaxis()->GetTitle());
  mgNSY->GetYaxis()->SetTitle(fNSYieldsInc->GetTitle());
  legYields->Draw("SAME");

  cResults->Print(Form("%s/NearsideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYields.C",fOutputDir.Data()));
  cResults->SetLogy(0);

  // Near-side Widths
  mgNSRms->Add(fNSRmsInc);
  legRms->AddEntry(fNSRmsInc,"Inclusive","LEP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legRms->AddEntry(fNSRmsEP[iEPBin],fEPBinTitles[iEPBin],"LEP");
    mgNSRms->Add(fNSRmsEP[iEPBin]);
  }
  mgNSRms->Draw(sDrawStyle);
  mgNSRms->GetXaxis()->SetTitle(fNSRmsInc->GetXaxis()->GetTitle());
  mgNSRms->GetYaxis()->SetTitle(fNSRmsInc->GetTitle());
  legRms->Draw("SAME");
  cResults->Print(Form("%s/NearsideRms.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRms.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRms.C",fOutputDir.Data()));

  mgNSSigmas->Add(fNSSigmasInc);
  legSigmas->AddEntry(fNSSigmasInc,"Inclusive","LEP");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    legSigmas->AddEntry(fNSSigmasEP[iEPBin],fEPBinTitles[iEPBin],"LEP");
    mgNSSigmas->Add(fNSSigmasEP[iEPBin]);
  }
  mgNSSigmas->Draw(sDrawStyle);
  mgNSSigmas->GetXaxis()->SetTitle(fNSSigmasInc->GetXaxis()->GetTitle());
  mgNSSigmas->GetYaxis()->SetTitle(fNSSigmasInc->GetTitle());
  legSigmas->Draw("SAME");
  cResults->Print(Form("%s/NearsideSigmas.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideSigmas.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideSigmas.C",fOutputDir.Data()));



  // Away-side Yields
  mgASY->Add(fASYieldsInc);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    mgASY->Add(fASYieldsEP[iEPBin]);
  }
  cResults->SetLogy(1);
  mgASY->Draw(sDrawStyle);
  mgASY->GetXaxis()->SetTitle(fASYieldsInc->GetXaxis()->GetTitle());
  mgASY->GetYaxis()->SetTitle(fASYieldsInc->GetTitle());
  legYields->Draw("SAME");


  cResults->Print(Form("%s/AwaysideYields.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYields.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYields.C",fOutputDir.Data()));

  cResults->SetLogy(0);
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

  mgASSigmas->Add(fASSigmasInc);
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    mgASSigmas->Add(fASSigmasEP[iEPBin]);
  }
  mgASSigmas->Draw(sDrawStyle);
  mgASSigmas->GetXaxis()->SetTitle(fASSigmasInc->GetXaxis()->GetTitle());
  mgASSigmas->GetYaxis()->SetTitle(fASSigmasInc->GetTitle());
  legSigmas->Draw("SAME");
  cResults->Print(Form("%s/AwaysideSigmas.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideSigmas.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideSigmas.C",fOutputDir.Data()));

  // =============================================================
  // Ratios
  // =============================================================

  // OutOverIn and MidOverIn AS
  TMultiGraph * mgASYieldRatios = new TMultiGraph();

  //TLegend * legASYieldRatios = new TLegend(0.25,0.65,0.45,0.85);
  //TLegend * legASYieldMidRatios = new TLegend(0.55,0.65,0.85,0.85);
  TLegend * legASYieldRatios = new TLegend(0.40,0.30);
  TLegend * legASYieldMidRatios = new TLegend(0.40,0.30);

  if (nRPFVariants > 0) {
    mgASYieldRatios->Add(OutOverIn_AS_RPFError);
   // mgASYieldRatios->Add(MidOverIn_AS_RPFError);
  }

  RemoveXErrorBars(OutOverIn_AS);

  mgASYieldRatios->Add(OutOverIn_AS);
  //mgASYieldRatios->Add(MidOverIn_AS);

  legASYieldRatios->SetHeader("Awayside Yield Ratios","C");
  legASYieldRatios->SetBorderSize(0);
  legASYieldRatios->AddEntry(OutOverIn_AS,"Out/In (Stat. Error)","LEP");
 // legASYieldRatios->AddEntry(MidOverIn_AS,"Mid/In (Stat. Error)","LP");

  legASYieldMidRatios->SetHeader("Awayside Yield Ratios","C");
  legASYieldMidRatios->SetBorderSize(0);
  legASYieldMidRatios->AddEntry(MidOverIn_AS,"Mid/In (Stat. Error)","LEP");

  mgASYieldRatios->GetXaxis()->SetTitle(OutOverIn_AS->GetXaxis()->GetTitle());
  mgASYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
//  mgASYieldRatios->Draw(sDrawStyle);

  //OutOverIn_AS->Draw(sDrawStyle);
  //MidOverIn_AS->Draw("SAME "+sDrawStyle);

  OutOverIn_AS->Draw("AP");
  //MidOverIn_AS->Draw("SAME LP");


  if (nRPFVariants > 0) {
    printf("Drawing the error from RPF Variants\n");
    OutOverIn_AS_RPFError->Draw("[]5 SAME");

    printf("Debug: (%f,%f) Err (%f,%f)\n",OutOverIn_AS_RPFError->GetX()[3],OutOverIn_AS_RPFError->GetY()[3],OutOverIn_AS_RPFError->GetEX()[3],OutOverIn_AS_RPFError->GetEY()[3]);

    //OutOverIn_AS_RPFError->Draw("E2 SAME");
    //OutOverIn_AS_RPFError->Draw("L[]5 SAME");
    //MidOverIn_AS_RPFError->Draw("E[]2 SAME");
    //OutOverIn_AS_RPFError->Draw("E2 SAME");
   // MidOverIn_AS_RPFError->Draw("E2 SAME");

//    mgASYieldRatios->Add(OutOverIn_AS_RPFError);
//    mgASYieldRatios->Add(MidOverIn_AS_RPFError);
    //mgASYieldRatios->Draw(Form("%s SAME",sDrawStyle.Data()));
 //   mgASYieldRatios->Draw("SAME");
    OutOverIn_AS->Draw("P SAME");

    legASYieldRatios->AddEntry(OutOverIn_AS_RPFError,"RPF Error (Out/In)","F");
//    legASYieldRatios->AddEntry(MidOverIn_AS_RPFError,"RPF Error (Mid/In)","F");


    //MidOverIn_AS->Draw("SAME LP");
    //OutOverIn_AS->Draw("SAME LP");
    //MidOverIn_AS->Draw("SAME "+sDrawStyle);
    //OutOverIn_AS->Draw("SAME "+sDrawStyle);

  } else {
    printf("  Not adding RPF variation error because I don't have it\n");
  }
  //mgASYieldRatios->Draw("ALP");
  //mgASYieldRatios->Draw("P SAME");

//  OutOverIn_AS_RPFError->Draw("[] SAME");
  mgASYieldRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
  //mgASYieldRatios->GetYaxis()->SetRangeUser(0.2,2.0);
  OutOverIn_AS->GetYaxis()->SetRangeUser(0.2,2.0);

  legASYieldRatios->Draw("SAME");
  TLegend * legendAlice = DrawAliceLegend(OutOverIn_AS,0.4,0.6,kAliceLegendWidth,kAliceLegendHeight);

  cResults->Print(Form("%s/AwaysideYieldOutOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldOutOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldOutOverInRatioRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    //OutOverIn_AS_SysError->Draw("LP[]5 SAME");
    OutOverIn_AS_SysError->Draw("[]5 SAME");
    legASYieldRatios->AddEntry(OutOverIn_AS_SysError,"Out/In (Sys. Error)","F");
    OutOverIn_AS->Draw("P SAME");
  }


  cResults->Print(Form("%s/AwaysideYieldOutOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldOutOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldOutOverInRatio.C",fOutputDir.Data()));

  if (fModelNames.size() > 0 && OutOverIn_AS_Models.size() > 0) {
    printf("Printing with model comparisons from %d models\n", (unsigned int) OutOverIn_AS_Models.size());

    OutOverIn_AS->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
    OutOverIn_AS->Draw("AP");
    for (unsigned int i = 0; i < fModelNames.size(); i++) {
      legASYieldRatios->AddEntry(OutOverIn_AS_Models[i],fModelTitles[i],"LF");
      //legASYieldRatios->AddEntry(OutOverIn_AS_Models[i],fModelNames[i],"LP");
      OutOverIn_AS_Models[i]->Draw("3 L SAME");
    }
    OutOverIn_AS_RPFError->Draw("[]5 SAME");
    OutOverIn_AS_SysError->Draw("[]5 SAME");
    OutOverIn_AS->Draw("P SAME");
    legASYieldRatios->Draw("SAME");
    legendAlice->Draw("SAME");

    cResults->Print(Form("%s/AwaysideYieldOutOverInRatioWithModels.pdf",fOutputDir.Data()));
    cResults->Print(Form("%s/AwaysideYieldOutOverInRatioWithModels.png",fOutputDir.Data()));
    cResults->Print(Form("%s/CFiles/AwaysideYieldOutOverInRatioWithModels.C",fOutputDir.Data()));
  }

  MidOverIn_AS->Draw("AP");
  if (nRPFVariants > 0) {
    printf("Drawing the error from RPF Variants\n");
    MidOverIn_AS_RPFError->Draw("[]5 SAME");
    legASYieldMidRatios->AddEntry(MidOverIn_AS_RPFError,"RPF Error (Mid/In)","F");
  }

  legASYieldMidRatios->Draw("SAME");

  cResults->Print(Form("%s/AwaysideYieldMidOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldMidOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldMidOverInRatioRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    MidOverIn_AS_SysError->Draw("L[]5 SAME");
    legASYieldMidRatios->AddEntry(MidOverIn_AS_SysError,"Mid/In (Sys. Error)","F");
    MidOverIn_AS->Draw("P SAME");
  }

  cResults->Print(Form("%s/AwaysideYieldMidOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideYieldMidOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideYieldMidOverInRatio.C",fOutputDir.Data()));

  if (fModelNames.size() > 0 && MidOverIn_AS_Models.size() > 0) {
    printf("Printing with model comparisons from %d models\n", (unsigned int) MidOverIn_AS_Models.size());

    for (unsigned int i = 0; i < fModelNames.size(); i++) {
      legASYieldMidRatios->AddEntry(MidOverIn_AS_Models[i],fModelNames[i],"LP");
      MidOverIn_AS_Models[i]->Draw("LP SAME");
    }

    cResults->Print(Form("%s/AwaysideYieldMidOverInRatioWithModels.pdf",fOutputDir.Data()));
    cResults->Print(Form("%s/AwaysideYieldMidOverInRatioWithModels.png",fOutputDir.Data()));
    cResults->Print(Form("%s/CFiles/AwaysideYieldMidOverInRatioWithModels.C",fOutputDir.Data()));
  }



  // OutOverIn and MidOverIn NS
  TMultiGraph * mgNSYieldRatios = new TMultiGraph();

  TLegend * legYieldRatios = new TLegend(0.55,0.65,0.85,0.85);

  RemoveXErrorBars(OutOverIn_NS);

  mgNSYieldRatios->Add(OutOverIn_NS);
  //mgNSYieldRatios->Add(MidOverIn_NS);

  legYieldRatios->SetBorderSize(0);
  legYieldRatios->SetHeader("Nearside Yield Ratios","C");
  legYieldRatios->AddEntry(OutOverIn_NS,"Out/In (Stat. Error)","LEP");
  //legYieldRatios->AddEntry(MidOverIn_NS,"Mid/In (Stat. Error)","LP");

  mgNSYieldRatios->GetXaxis()->SetTitle(OutOverIn_NS->GetXaxis()->GetTitle());
  mgNSYieldRatios->GetYaxis()->SetTitle("Yield Ratio");
  mgNSYieldRatios->Draw(sDrawStyle);

  if (nRPFVariants > 0) {
    OutOverIn_NS_RPFError->Draw("E3 SAME");
    //MidOverIn_NS_RPFError->Draw("E3 SAME");

    legYieldRatios->AddEntry(OutOverIn_NS_RPFError,"RPF Error (Out/In)","F");
    //legYieldRatios->AddEntry(MidOverIn_NS_RPFError,"RPF Error (Mid/In)","F");
  }

  legYieldRatios->Draw("SAME");
  mgNSYieldRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
  TLegend * legendAliceNearsideYieldRatio = DrawAliceLegend(OutOverIn_AS,0.4,0.3,kAliceLegendWidth,kAliceLegendHeight);

  cResults->Print(Form("%s/NearsideYieldRatiosRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYieldRatiosRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYieldRatiosRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    OutOverIn_NS_SysError->Draw("L[]5 SAME");
    legYieldRatios->AddEntry(OutOverIn_NS_SysError,"Out/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/NearsideYieldRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideYieldRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideYieldRatios.C",fOutputDir.Data()));






  if (fModelNames.size() > 0 && OutOverIn_NS_Models.size() > 0) {
    printf("Printing with model comparisons from %d models\n", (unsigned int) OutOverIn_NS_Models.size());

    for (unsigned int i = 0; i < fModelNames.size(); i++) {
      legYieldRatios->AddEntry(OutOverIn_NS_Models[i],fModelTitles[i],"LF");
      //legYieldRatios->AddEntry(OutOverIn_NS_Models[i],fModelNames[i],"LP");
      OutOverIn_NS_Models[i]->Draw("3 L SAME");
    }

    cResults->Print(Form("%s/NearsideYieldOutOverInRatioWithModels.pdf",fOutputDir.Data()));
    cResults->Print(Form("%s/NearsideYieldOutOverInRatioWithModels.png",fOutputDir.Data()));
    cResults->Print(Form("%s/CFiles/NearsideYieldOutOverInRatioWithModels.C",fOutputDir.Data()));
  }



  // Draw RMS Ratios together (out/in and mid/in)
  TMultiGraph * mgNSRmsRatiosBoth = new TMultiGraph();
  TLegend * legNSRmsRatiosBoth = new TLegend(0.55,0.65,0.85,0.85);

  RemoveXErrorBars(RmsOutOverIn_NS);
  RemoveXErrorBars(RmsMidOverIn_NS);

  mgNSRmsRatiosBoth->Add(RmsOutOverIn_NS);
  mgNSRmsRatiosBoth->Add(RmsMidOverIn_NS);
  legNSRmsRatiosBoth->SetBorderSize(0);
  legNSRmsRatiosBoth->SetHeader("Nearside RMS Ratios","C");
  legNSRmsRatiosBoth->AddEntry(RmsOutOverIn_NS,"Out/In (Stat. Error)","LEP");
  legNSRmsRatiosBoth->AddEntry(RmsMidOverIn_NS,"Mid/In (Stat. Error)","LEP");

  mgNSRmsRatiosBoth->GetXaxis()->SetTitle(RmsOutOverIn_NS->GetXaxis()->GetTitle());
  mgNSRmsRatiosBoth->GetYaxis()->SetTitle("RMS Ratio");

  mgNSRmsRatiosBoth->Draw(sDrawStyle);


  if (nRPFVariants > 0) {
    RmsOutOverIn_NS_RPFError->Draw("E3 SAME");
    RmsMidOverIn_NS_RPFError->Draw("E3 SAME");

    legNSRmsRatiosBoth->AddEntry(RmsOutOverIn_NS_RPFError,"RPF Error (Out/In)","F");
    legNSRmsRatiosBoth->AddEntry(RmsMidOverIn_NS_RPFError,"RPF Error (Mid/In)","F");
  }


  legNSRmsRatiosBoth->Draw("SAME");
  mgNSRmsRatiosBoth->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/NearsideRmsRatiosRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsRatiosRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsRatiosRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    RmsOutOverIn_NS_SysError->Draw("L[]5 SAME");
    legNSRmsRatiosBoth->AddEntry(RmsOutOverIn_NS_SysError,"Out/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/NearsideRmsRatios.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsRatios.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsRatios.C",fOutputDir.Data()));




  // Awayside
  TMultiGraph * mgASRmsRatios = new TMultiGraph();
  TLegend * legASRmsRatios = new TLegend(0.55,0.65,0.85,0.85);
  TLegend * legASRmsMidRatios = new TLegend(0.55,0.65,0.85,0.85);

  RemoveXErrorBars(RmsOutOverIn_AS);
  RemoveXErrorBars(RmsMidOverIn_AS);

  mgASRmsRatios->Add(RmsOutOverIn_AS);
  mgASRmsRatios->Add(RmsMidOverIn_AS);

  legASRmsRatios->SetBorderSize(0);
  legASRmsRatios->SetHeader("Awayside RMS Ratios","C");

  legASRmsMidRatios->SetBorderSize(0);
  legASRmsMidRatios->SetHeader("Awayside RMS Ratios","C");



  legASRmsRatios->AddEntry(RmsOutOverIn_AS,"Out/In (Stat. Error)","LEP");
  //legASRmsRatios->AddEntry(RmsMidOverIn_AS,"Mid/In (Stat. Error)","LP");

  //mgASRmsRatios->GetXaxis()->SetTitle(RmsOutOverIn_AS->GetXaxis()->GetTitle());
  //mgASRmsRatios->GetYaxis()->SetTitle("RMS Ratio");

  //mgASRmsRatios->Draw(sDrawStyle);

  RmsOutOverIn_AS->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
  RmsOutOverIn_AS->Draw(sDrawStyle);
  RemoveXErrorBars(RmsOutOverIn_AS);


  if (nRPFVariants > 0) {
    RmsOutOverIn_AS_RPFError->Draw("[]5 SAME");
    //RmsMidOverIn_AS_RPFError->Draw("E3 SAME");

    legASRmsRatios->AddEntry(RmsOutOverIn_AS_RPFError,"RPF Error (Out/In)","F");
  }


  legASRmsRatios->Draw("SAME");
  mgASRmsRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  TLegend * legendAliceAwaysideRmsRatio = DrawAliceLegend(RmsOutOverIn_AS,0,0,kAliceLegendWidth,kAliceLegendHeight);

  cResults->Print(Form("%s/AwaysideRmsOutOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRmsOutOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRmsOutOverInRatioRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    RmsOutOverIn_AS_SysError->Draw("[]5 SAME");
    legASRmsRatios->AddEntry(RmsOutOverIn_AS_SysError,"Out/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/AwaysideRmsOutOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRmsOutOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRmsOutOverInRatio.C",fOutputDir.Data()));

  
  if (fModelNames.size() > 0 && RmsOutOverIn_AS_Models.size() > 0) {
    printf("Printing with model comparisons from %d models\n", (unsigned int) RmsOutOverIn_AS_Models.size());

    for (unsigned int i = 0; i < fModelNames.size(); i++) {
      legASRmsRatios->AddEntry(RmsOutOverIn_AS_Models[i],fModelTitles[i],"LEP");
      //legYieldRatios->AddEntry(OutOverIn_NS_Models[i],fModelNames[i],"LP");
      RmsOutOverIn_AS_Models[i]->Draw("3 L SAME");
      //RmsOutOverIn_AS_Models[i]->Draw("LP SAME");
    }
    legendAliceAwaysideRmsRatio->Draw("SAME");

    cResults->Print(Form("%s/AwaysideRmsOutOverInRatioWithModels.pdf",fOutputDir.Data()));
    cResults->Print(Form("%s/AwaysideRmsOutOverInRatioWithModels.png",fOutputDir.Data()));
    cResults->Print(Form("%s/CFiles/AwaysideRmsOutOverInRatioWithModels.C",fOutputDir.Data()));
  }

  RmsMidOverIn_AS->Draw(sDrawStyle);
  RemoveXErrorBars(RmsMidOverIn_AS);

  if (nRPFVariants > 0) {
    RmsMidOverIn_AS_RPFError->Draw("L[]5 SAME");
    legASRmsMidRatios->AddEntry(RmsMidOverIn_AS_RPFError,"RPF Error (Mid/In)","F");
  }

  legASRmsMidRatios->Draw("SAME");

  cResults->Print(Form("%s/AwaysideRmsMidOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRmsMidOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRmsMidOverInRatioRPFError.C",fOutputDir.Data()));


  if (fSystematicsNames.size() > 0) {
    // FIXME make this exist
    RmsMidOverIn_AS_SysError->Draw("L[]5 SAME");
    legASRmsMidRatios->AddEntry(RmsMidOverIn_AS_SysError,"Mid/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/AwaysideRmsMidOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/AwaysideRmsMidOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/AwaysideRmsMidOverInRatio.C",fOutputDir.Data()));



  // Nearside RMS Ratios
  TMultiGraph * mgNSRmsRatios = new TMultiGraph();
  TLegend * legNSRmsRatios = new TLegend(0.55,0.65,0.85,0.85);
  TLegend * legNSRmsMidRatios = new TLegend(0.55,0.65,0.85,0.85);

  
  RemoveXErrorBars(RmsOutOverIn_NS);
  RemoveXErrorBars(RmsMidOverIn_NS);

  mgNSRmsRatios->Add(RmsOutOverIn_NS);
  mgNSRmsRatios->Add(RmsMidOverIn_NS);

  legNSRmsRatios->SetBorderSize(0);
  legNSRmsRatios->SetHeader("Nearside RMS Ratios","C");

  legNSRmsMidRatios->SetBorderSize(0);
  legNSRmsMidRatios->SetHeader("Nearside RMS Ratios","C");



  legNSRmsRatios->AddEntry(RmsOutOverIn_NS,"Out/In (Stat. Error)","LEP");
  //legNSRmsRatios->AddEntry(RmsMidOverIn_NS,"Mid/In (Stat. Error)","LP");

  //mgNSRmsRatios->GetXaxis()->SetTitle(RmsOutOverIn_AS->GetXaxis()->GetTitle());
  //mgNSRmsRatios->GetYaxis()->SetTitle("RMS Ratio");

  //mgNSRmsRatios->Draw(sDrawStyle);

  RmsOutOverIn_NS->Draw(sDrawStyle);
  RemoveXErrorBars(RmsOutOverIn_NS);

  if (nRPFVariants > 0) {
    RmsOutOverIn_NS_RPFError->Draw("[]5 SAME");
    //RmsMidOverIn_NS_RPFError->Draw("E3 SAME");

    legNSRmsRatios->AddEntry(RmsOutOverIn_NS_RPFError,"RPF Error (Out/In)","F");
  }


  legNSRmsRatios->Draw("SAME");
  mgNSRmsRatios->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);

  cResults->Print(Form("%s/NearsideRmsOutOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsOutOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsOutOverInRatioRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    RmsOutOverIn_NS_SysError->Draw("L[]5 SAME");
    legNSRmsRatios->AddEntry(RmsOutOverIn_NS_SysError,"Out/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/NearsideRmsOutOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsOutOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsOutOverInRatio.C",fOutputDir.Data()));

  
  if (fModelNames.size() > 0 && RmsOutOverIn_NS_Models.size() > 0) {
    printf("Printing with model comparisons from %d models\n", (unsigned int) RmsOutOverIn_NS_Models.size());

    for (unsigned int i = 0; i < fModelNames.size(); i++) {
      legNSRmsRatios->AddEntry(RmsOutOverIn_NS_Models[i],fModelTitles[i],"LF");
      //legYieldRatios->AddEntry(OutOverIn_NS_Models[i],fModelNames[i],"LP");
      RmsOutOverIn_NS_Models[i]->Draw("3 L SAME");
      //RmsOutOverIn_AS_Models[i]->Draw("LP SAME");
    }

    cResults->Print(Form("%s/NearsideRmsOutOverInRatioWithModels.pdf",fOutputDir.Data()));
    cResults->Print(Form("%s/NearsideRmsOutOverInRatioWithModels.png",fOutputDir.Data()));
    cResults->Print(Form("%s/CFiles/NearsideRmsOutOverInRatioWithModels.C",fOutputDir.Data()));
  }

  RmsMidOverIn_NS->Draw(sDrawStyle);
  RemoveXErrorBars(RmsMidOverIn_NS);

  if (nRPFVariants > 0) {
    RmsMidOverIn_NS_RPFError->Draw("L[]5 SAME");
    legNSRmsMidRatios->AddEntry(RmsMidOverIn_NS_RPFError,"RPF Error (Mid/In)","F");
  }

  legASRmsMidRatios->Draw("SAME");

  cResults->Print(Form("%s/NearsideRmsMidOverInRatioRPFError.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsMidOverInRatioRPFError.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsMidOverInnatioRPFError.C",fOutputDir.Data()));

  if (fSystematicsNames.size() > 0) {
    // FIXME make this exist
    RmsMidOverIn_NS_SysError->Draw("L[]5 SAME");
    legNSRmsMidRatios->AddEntry(RmsMidOverIn_NS_SysError,"Mid/In (Sys. Error)","F");
  }

  cResults->Print(Form("%s/NearsideRmsMidOverInRatio.pdf",fOutputDir.Data()));
  cResults->Print(Form("%s/NearsideRmsMidOverInRatio.png",fOutputDir.Data()));
  cResults->Print(Form("%s/CFiles/NearsideRmsMidOverInRatio.C",fOutputDir.Data()));



}


void TaskCalcObservables::ProduceSystematicUncertPlots() {

  if (fSystematicsNames.size() == 0) {
    printf("Not producing systematic uncertainty plots, as no systematics have been given\n");
    return;
  }


  TGraphErrors * first = OutOverIn_AS;
  printf("Trying to plot systematics for object %s\n",OutOverIn_AS->GetName());

  vector<TGraphErrors *> firstArray = OutOverIn_AS_SysErrorBySource;

  TGraphErrors * firstTotalSysUncert = OutOverIn_AS_SysError;

  ProduceSystematicUncertPlotIndiv(first,firstArray,firstTotalSysUncert,0.75,0.45);

  //ProduceSystematicUncertPlotIndiv(,fNSYieldsInc_SysErrorBySource);


  ProduceSystematicUncertPlotIndiv(OutOverIn_NS,OutOverIn_NS_SysErrorBySource,OutOverIn_NS_SysError,0.75,0.45);

  ProduceSystematicUncertPlotIndiv(RmsOutOverIn_AS,RmsOutOverIn_AS_SysErrorBySource,RmsOutOverIn_AS_SysError,0.75,0.5);
  ProduceSystematicUncertPlotIndiv(RmsOutOverIn_NS,RmsOutOverIn_NS_SysErrorBySource,RmsOutOverIn_NS_SysError,0.75,0.5);

  ProduceSystematicUncertPlotIndiv(SigmasOutOverIn_AS,SigmasOutOverIn_AS_SysErrorBySource,SigmasOutOverIn_AS_SysError,0.75,0.5);
  ProduceSystematicUncertPlotIndiv(SigmasOutOverIn_NS,SigmasOutOverIn_NS_SysErrorBySource,SigmasOutOverIn_NS_SysError,0.75,0.5);
  
  ProduceSystematicUncertPlotIndiv(fASYieldsInc,fASYieldsInc_SysErrorBySource,fASYieldsInc_SysError,0.75,0.5);
  ProduceSystematicUncertPlotIndiv(fNSYieldsInc,fNSYieldsInc_SysErrorBySource,fNSYieldsInc_SysError,0.75,0.5);


  if (fASYieldsEP_SysErrorBySource.size() > 0 && fASYieldsEP_SysError.size()) {
    for (int i = 0; i < kNEPBins; i++) {
      ProduceSystematicUncertPlotIndiv(fASYieldsEP[i],fASYieldsEP_SysErrorBySource[i],fASYieldsEP_SysError[i],0.75,0.5);
      ProduceSystematicUncertPlotIndiv(fNSYieldsEP[i],fNSYieldsEP_SysErrorBySource[i],fNSYieldsEP_SysError[i],0.75,0.5);
    }

    ProduceSystematicUncertPlotIndiv(fASRmsInc,fASRmsInc_SysErrorBySource,fASRmsInc_SysError,0.75,0.5);
    ProduceSystematicUncertPlotIndiv(fNSRmsInc,fNSRmsInc_SysErrorBySource,fNSRmsInc_SysError,0.75,0.5);

    ProduceSystematicUncertPlotIndiv(fASSigmasInc,fASSigmasInc_SysErrorBySource,fASSigmasInc_SysError,0.7,0.5);
    ProduceSystematicUncertPlotIndiv(fNSSigmasInc,fNSSigmasInc_SysErrorBySource,fNSSigmasInc_SysError,0.7,0.5);
    for (int i = 0; i < kNEPBins; i++) {
      ProduceSystematicUncertPlotIndiv(fASSigmasEP[i],fASSigmasEP_SysErrorBySource[i],fASSigmasEP_SysError[i],0.75,0.5);
      ProduceSystematicUncertPlotIndiv(fNSSigmasEP[i],fNSSigmasEP_SysErrorBySource[i],fNSSigmasEP_SysError[i],0.75,0.5);
    }
  }
}



double fBernsteinPoly3(double * x, double * p) {
  double value = 0;
  double xMin = p[0];
  double xMax = p[1];

  double y = (x[0] - xMin)/(xMax - xMin); // y in 0..1

  value += p[2] * TMath::Power(1-y,3);
  value += 3*p[3] * y * TMath::Power(1-y,2);
  value += 3*p[4] * TMath::Power(y,2) * (1-y);
  value += p[5] * TMath::Power(y,3);

  return value;
}

double fBernsteinPoly4(double * x, double * p) {
  double value = 0;
  double xMin = p[0];
  double xMax = p[1];

  double y = (x[0] - xMin)/(xMax - xMin); // y in 0..1

  value += p[2] * TMath::Power(1-y,4);
  value += 4*p[3] * y * TMath::Power(1-y,3);
  value += 6*p[4] * TMath::Power(y,2) * TMath::Power(1-y,2);
  value += 4*p[5] * TMath::Power(y,3) * (1-y);
  value += p[6] * TMath::Power(y,4);

  return value;
}


// Add input for the total systematic uncertainty
//void TaskCalcObservables::ProduceSystematicUncertPlotIndiv(TGraphErrors * fCentral, vector<TGraphErrors *> fSysErrors) {
void TaskCalcObservables::ProduceSystematicUncertPlotIndiv(TGraphErrors * fCentral, vector<TGraphErrors *> fSysErrors, TGraphErrors * fTotalSysErrors, float legendX, float legendY) {
  if (fSysErrors.size() == 0) {
    printf("Systematic uncertainties not loaded for %s, skipping.\n",fCentral->GetName());
    return;
  }
  if (fSysErrors.size() > fSystematicsNames.size()) {
    fprintf(stderr,"Error mismatch in number for systematic error types given and number named\n");
    return;
  }

  printf("Drawing systematic errors by source for object %s\n",fCentral->GetName());

  //TGraphErrors * fFullErrorPlot = 

  vector<TGraphErrors *> fErrorPlots = {};
  vector<TF1 * > fErrorPlotFits = {};
  // Create Graphs with 0 for the y-axis points
  // Or
  // Create Graphs or Histos where the 

  printf("About to fit the systematics\n");

  for (int i = 0; i < (int) fSysErrors.size(); i++) {
    TGraphErrors * fErrorPlot = (TGraphErrors *) fSysErrors[i]->Clone(Form("%s_ErrorOnly",fSysErrors[i]->GetName()));
    for (int j = 0; j < (int) fErrorPlot->GetN(); j++) {
      //double denom = fErrorPlot->GetY()[j];
      double denom = fCentral->GetY()[j];
      if (denom == 0) continue;
      fErrorPlot->SetPoint(j,fErrorPlot->GetX()[j],100.*abs(fErrorPlot->GetEY()[j])/denom);
      fErrorPlot->SetPointError(j,fErrorPlot->GetEX()[j],0.0);
    }
    //TH1F * fErrorHist = new TH1F(Form(),Form(),


    //float minXDefault = 1.5;
    //float maxXDefault = 17;
    float minXDefault = 0;
    float maxXDefault = 20;

    //TString sFitName = Form("%s_UncertFit",graph->GetName());
    TString sFormula = "[0]";


    // 1 2 1
    // 1 3 3 1
    // 1 4 6 4 1

    TString sBernPol3 = Form("[0]*TMath::Power(%f-x,3) + 3*[1]*x*TMath::Power(%f-x,2)\
    + 3*[2]*x*x*(%f-x) + [3]*TMath::Power(x,3)",maxXDefault,maxXDefault,maxXDefault);
    TString sBernPol4 = Form("[0]*TMath::Power(%f-x,4) + 4*[1]*(x-%f)*TMath::Power(%f-x,3)\
    + 6*[2]*TMath::Power(x-%f,2)*TMath::Power(%f-x,2) + 4*[3]*TMath::Power(x-%f,3)*(x-%f)+[4]*TMath::Power(x-%f,4)",maxXDefault,
      minXDefault,maxXDefault,
      minXDefault,maxXDefault,
      minXDefault,maxXDefault,
      minXDefault);
    //TString sBernPol4 = Form("[0]*TMath::Power(%f-x,4) + 4*[1]*x*TMath::Power(%f-x,3)\
    //+ 6*[2]*x*x*TMath::Power(%f-x,2) + 4*[3]*(%f-x)*TMath::Power(x,3)+[4]*TMath::Power(x,4)",maxXDefault,
    //  maxXDefault,maxXDefault,maxXDefault);

    //sFormula = sBernPol3;
    sFormula = sBernPol4;

    //TF1 * fErrorPlotFit = new TF1(Form("%s_Fit",fErrorPlot->GetName()),sFormula.Data(),0,17);
    //TF1 * fErrorPlotFit = new TF1(Form("%s_Fit",fErrorPlot->GetName()),"[0]+[1]*x+[2]*x*x",0,17);

    int nBaseParams=2;

    int nPar = 4 + nBaseParams; // 2 for min and max
    TF1 * fErrorPlotFit = new TF1(Form("%s_Fit",fErrorPlot->GetName()),fBernsteinPoly3,0,17,nPar);
 
    //int nPar = 5 + nBaseParams; // 2 for min and max
    //TF1 * fErrorPlotFit = new TF1(Form("%s_Fit",fErrorPlot->GetName()),fBernsteinPoly4,0,17,nPar);


    fErrorPlotFit->FixParameter(0,minXDefault);
    fErrorPlotFit->FixParameter(1,maxXDefault);


    float fMaxErrorRange = 50; // note: this is in percent

    //double fGraphMax = fErrorPlot->GetMaximum();
    double fGraphMax = GetMaximumYFromGraphErrors(fErrorPlot);

    printf("DebugError: found error graph maximum = %f\n",fGraphMax);
    fMaxErrorRange = 2 * fGraphMax;

    for (int j = nBaseParams; j < fErrorPlotFit->GetNpar(); j++) {
      fErrorPlotFit->SetParameter(j,2);
      fErrorPlotFit->SetParLimits(j,0.0,fMaxErrorRange);
    }

    TF1 * fPol0 = new TF1(Form("%s_Pol0Fit",fErrorPlot->GetName()),"[0]",0,17);
    fErrorPlot->Fit(fPol0,"W");

    //fErrorPlot->Fit(fErrorPlotFit,"W");

    fErrorPlots.push_back(fErrorPlot);

    //fErrorPlotFits.push_back(fErrorPlotFit);
    fErrorPlotFits.push_back(fPol0);
  }
  printf("\tBuilt error graphs\n");
  // Maybe save these for comparison?


  TCanvas * cSysError = new TCanvas("cSysError","cSysError");

  float fAliceLegendXMin = 0.25;
  float fAliceLegendYMin = 0.7;

  float legendWidth = 0.2;
  float legendHeight = 0.35;
  TLegend * lLegError = new TLegend(legendX, legendY, legendX + legendWidth, legendY + legendHeight);
  //lLegError->SetHeader(Form("%s Systematic Uncertainty by Source",fCentral->GetTitle()));
  TMultiGraph * mg = new TMultiGraph();

  TF1 * fTotalOfSysErrorFits = (TF1 * ) fErrorPlotFits[0]->Clone(Form("%s_TotalOfFits",fTotalSysErrors->GetName()));
  TGraphErrors * fTotalFromErrorFits = (TGraphErrors *) fTotalSysErrors->Clone(Form("%s_TotalFromFits",fTotalSysErrors->GetName()));

  TGraphErrors * fTotalSysErrorsClone = (TGraphErrors *) fTotalSysErrors->Clone(Form("%s_Clone",fTotalSysErrors->GetName()));
  fTotalSysErrorsClone->SetLineColor(kBlack);
  fTotalSysErrorsClone->SetMarkerColor(kBlack);
  fTotalSysErrorsClone->SetMarkerStyle(kOpenSquare);

  for (int j = 0; j < (int) fTotalSysErrorsClone->GetN(); j++) {
    double denom = fTotalSysErrorsClone->GetY()[j];
    if (denom == 0) continue;
    fTotalSysErrorsClone->SetPoint(j,fTotalSysErrorsClone->GetX()[j],100.*abs(fTotalSysErrorsClone->GetEY()[j])/denom);
    fTotalSysErrorsClone->SetPointError(j,fTotalSysErrorsClone->GetEX()[j],0.0);
  }


  for (int i = 0; i < (int) fTotalSysErrorsClone->GetN();i++) {
    double fTotalFromErrors = 0;
    double fXValue = fTotalSysErrorsClone->GetX()[i]; 
    for (int j = 0; j < (int) fErrorPlots.size(); j++) {
      fTotalFromErrors += TMath::Power(fErrorPlotFits[j]->Eval(fXValue),2);
    }
    double fYValue = TMath::Sqrt(fTotalFromErrors);
    fTotalFromErrorFits->SetPoint(i,fXValue,fYValue);
    
    // Setting the final uncertainty to what comes from summing the smoothed fit
    if (bSmoothFitUncertainty) {
      double fOldXErr = fTotalSysErrorsClone->GetEX()[i];
      // need to convert out of percentage and from ratio
      double fYAbsoluteValue = fTotalSysErrors->GetY()[i];
      double fNewYErr = 0.01 * fYValue * fYAbsoluteValue;

      double fFinalErr = max(fNewYErr,fOldXErr);

      fTotalSysErrors->SetPointError(i,fOldXErr,fNewYErr);

    }


  }


  fTotalFromErrorFits->SetLineColor(kBlack);
  fTotalFromErrorFits->SetMarkerColor(kBlack);

  lLegError->AddEntry(fTotalFromErrorFits,"Total From Fits","LP");
  lLegError->AddEntry(fTotalSysErrorsClone,"Total","LP");

  vector<TF1 * > fitFuncs = {};

  for (int i = 0; i < (int) fErrorPlots.size(); i++) {
    TF1 * fit = 0;

    //fit = (TF1 *) fErrorPlots[i]->FindObject(Form("%s_UncertFit",fErrorPlots[i]->GetName()));
    //if (fit == 0) {
    //  printf("Could not find fit for uncertainty graph %s\n",fErrorPlots[i]->GetName());
    //} else {
    fit = fErrorPlotFits[i];
    //} 
    
  
    fErrorPlots[i]->SetLineColor(kSysErrorColor[i]);
    fErrorPlots[i]->SetMarkerColor(kSysErrorColor[i]);
    fErrorPlots[i]->SetMarkerStyle(kSysErrorStyle[i]);
    mg->Add(fErrorPlots[i]);
    lLegError->AddEntry(fErrorPlots[i],fSystematicsNames[i],"PL");
  }
  mg->Add(fTotalSysErrorsClone);
  mg->Add(fTotalFromErrorFits);

  mg->GetXaxis()->SetTitle(fErrorPlots[0]->GetXaxis()->GetTitle());
  mg->GetYaxis()->SetTitle(Form("%s Sys, Uncert. by Source (%%)",fCentral->GetTitle()));
  mg->Draw("AP");
  for (int i = 0; i < (int) fErrorPlotFits.size(); i++) {
    fErrorPlotFits[i]->SetLineColor(kSysErrorColor[i]);
    fErrorPlotFits[i]->Draw("SAME");
  }
  printf("  the error plot fits had %d entries\n",(int) fErrorPlotFits.size());
  lLegError->Draw("SAME");

  TLegend * legendAliceNearsideYieldRatio = DrawAliceLegend(OutOverIn_AS,fAliceLegendXMin,fAliceLegendYMin,kAliceLegendWidth,kAliceLegendHeight);

  cSysError->Print(Form("%s/%s_SysErrByType.pdf",fOutputDir.Data(),fCentral->GetName()));
  cSysError->Print(Form("%s/%s_SysErrByType.png",fOutputDir.Data(),fCentral->GetName()));
}

// Combine raw statistical uncertainty with RPF to create total stat errors:
// e.g: fNSYieldsInc+/fNSYieldsInc_RPFError->fNSYieldsInc_StatError
void TaskCalcObservables::CombineStatisticalErrors() {
  // Start with Nearside Yields

  // Incl
  CombineStatisticalErrorsIndiv(fNSYieldsInc,fNSYieldsInc_RPFError); 
  CombineStatisticalErrorsIndiv(fASYieldsInc,fASYieldsInc_RPFError); 
  CombineStatisticalErrorsIndiv(fNSRmsInc,fNSRmsInc_RPFError); 
  CombineStatisticalErrorsIndiv(fASRmsInc,fASRmsInc_RPFError); 

  // By EPBin
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    CombineStatisticalErrorsIndiv(fNSYieldsEP[iEPBin],fNSYieldsEP_RPFError[iEPBin]); 
    CombineStatisticalErrorsIndiv(fASYieldsEP[iEPBin],fASYieldsEP_RPFError[iEPBin]); 
    CombineStatisticalErrorsIndiv(fNSRmsEP[iEPBin],fNSRmsEP_RPFError[iEPBin]); 
    CombineStatisticalErrorsIndiv(fASRmsEP[iEPBin],fASRmsEP_RPFError[iEPBin]); 
    CombineStatisticalErrorsIndiv(fNSSigmasEP[iEPBin],fNSSigmasEP_RPFError[iEPBin]); 
    CombineStatisticalErrorsIndiv(fASSigmasEP[iEPBin],fASSigmasEP_RPFError[iEPBin]); 
  }

  // Ratios
  CombineStatisticalErrorsIndiv(OutOverIn_AS,OutOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(OutOverIn_NS,OutOverIn_NS_RPFError);
  CombineStatisticalErrorsIndiv(MidOverIn_AS,MidOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(MidOverIn_NS,MidOverIn_NS_RPFError);

  CombineStatisticalErrorsIndiv(RmsOutOverIn_AS,RmsOutOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(RmsOutOverIn_NS,RmsOutOverIn_NS_RPFError);
  CombineStatisticalErrorsIndiv(RmsMidOverIn_AS,RmsMidOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(RmsMidOverIn_NS,RmsMidOverIn_NS_RPFError);

  CombineStatisticalErrorsIndiv(SigmasOutOverIn_AS,SigmasOutOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(SigmasOutOverIn_NS,SigmasOutOverIn_NS_RPFError);
  CombineStatisticalErrorsIndiv(SigmasMidOverIn_AS,SigmasMidOverIn_AS_RPFError);
  CombineStatisticalErrorsIndiv(SigmasMidOverIn_NS,SigmasMidOverIn_NS_RPFError);
}


void TaskCalcObservables::CombineStatisticalErrorsIndiv(TGraphErrors * fRawStatErrors, TGraphErrors * fRPFErrors) {

  int n = fRawStatErrors->GetN();
  for (int i = 0; i < n; i++) {
    double fInitStatError = fRawStatErrors->GetEY()[i];
    double fRPFStatError = fRPFErrors->GetEY()[i];

    double fFinalStatError = TMath::Sqrt(fInitStatError*fInitStatError + fRPFStatError*fRPFStatError);

    fRawStatErrors->SetPointError(i,fRawStatErrors->GetEX()[i],fFinalStatError);
  }

}

void TaskCalcObservables::DrawFinalResults() {
  TCanvas * cFinal = new TCanvas("FinalObs","FinalObs");

  //DrawFinalObservable(TGraphErrors * fObsGraph, TGraphErrors * fObsGraphSysErrors, vector<TGraphErrors*> fModels)
  // Start with Nearside Yields


  DrawFinalObservable(fNSYieldsInc,fNSYieldsInc_SysError,fNSYieldsInc_Models,cFinal,"NSYieldsInc","Nearside Yields (Incl.)");
  DrawFinalObservable(fASYieldsInc,fASYieldsInc_SysError,fASYieldsInc_Models,cFinal,"ASYieldsInc","Awayside Yields (Incl.)");
  DrawFinalObservable(fNSRmsInc,fNSRmsInc_SysError,fNSRmsInc_Models,cFinal,"NSRmsInc","Nearside RMS (Incl.)");
  DrawFinalObservable(fASRmsInc,fASRmsInc_SysError,fASRmsInc_Models,cFinal,"ASRmsInc","Awayside RMS (Incl.)");
  DrawFinalObservable(fNSSigmasInc,fNSSigmasInc_SysError,fNSSigmasInc_Models,cFinal,"NSSigmasInc","Nearside Sigmas (Incl.)");
  DrawFinalObservable(fASSigmasInc,fASSigmasInc_SysError,fASSigmasInc_Models,cFinal,"ASSigmasInc","Awayside Sigmas (Incl.)");


  
  printf("Now to draw observables by EPBin\n");
  for (int iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {

    if (fModelNames.size() > 0) {
      DrawFinalObservable(fNSYieldsEP[iEPBin],fNSYieldsEP_SysError[iEPBin],fNSYieldsEP_Models[iEPBin],cFinal,Form("NSYieldsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      DrawFinalObservable(fASYieldsEP[iEPBin],fASYieldsEP_SysError[iEPBin],fASYieldsEP_Models[iEPBin],cFinal,Form("ASYieldsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      DrawFinalObservable(fNSRmsEP[iEPBin],fNSRmsEP_SysError[iEPBin],fNSRmsEP_Models[iEPBin],cFinal,Form("NSRmsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      DrawFinalObservable(fASRmsEP[iEPBin],fASRmsEP_SysError[iEPBin],fASRmsEP_Models[iEPBin],cFinal,Form("ASRmsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data())); 

      DrawFinalObservable(fNSSigmasEP[iEPBin],fNSSigmasEP_SysError[iEPBin],fNSSigmasEP_Models[iEPBin],cFinal,Form("NSSigmasEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      DrawFinalObservable(fASSigmasEP[iEPBin],fASSigmasEP_SysError[iEPBin],fASSigmasEP_Models[iEPBin],cFinal,Form("ASSigmasEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data())); 
    } else {
      // This may need another check in case sys_errors are missing, and then the fNSYieldsEP_SysError array is empty, and can't be dereferenced
      //DrawFinalObservable(fNSYieldsEP[iEPBin],fNSYieldsEP_SysError[iEPBin],(vector<TGraphErrors*>) 0,cFinal,Form("NSYieldsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      //DrawFinalObservable(fASYieldsEP[iEPBin],fASYieldsEP_SysError[iEPBin],(vector<TGraphErrors*>) 0,cFinal,Form("ASYieldsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      //DrawFinalObservable(fNSRmsEP[iEPBin],fNSRmsEP_SysError[iEPBin],(vector<TGraphErrors*>)0,cFinal,Form("NSRmsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data()));
      //DrawFinalObservable(fASRmsEP[iEPBin],fASRmsEP_SysError[iEPBin],(vector<TGraphErrors*>) 0,cFinal,Form("ASRmsEP%d",iEPBin),Form("%s",fEPBinTitles[iEPBin].Data())); 
    }
  }

  printf("Now to begin drawing the ratios\n");

  TString sAwaysideDesc = "Away-side (|#Delta#varphi-#pi|<#pi/3)";
  TString sNearsideDesc = "Near-side (|#Delta#varphi|<#pi/3)";


  DrawFinalObservable(OutOverIn_AS,OutOverIn_AS_SysError,OutOverIn_AS_Models,cFinal,"ASYieldOutOverIn",Form("%s Yield Out/In",sAwaysideDesc.Data()));
  DrawFinalObservable(OutOverIn_NS,OutOverIn_NS_SysError,OutOverIn_NS_Models,cFinal,"NSYieldOutOverIn",Form("%s Yield Out/In",sNearsideDesc.Data()));
  DrawFinalObservable(MidOverIn_AS,MidOverIn_AS_SysError,MidOverIn_AS_Models,cFinal,"ASYieldMidOverIn",Form("%s Yield Mid/In",sAwaysideDesc.Data()));
  DrawFinalObservable(MidOverIn_NS,MidOverIn_NS_SysError,MidOverIn_NS_Models,cFinal,"NSYieldMidOverIn",Form("%s Yield Mid/In",sNearsideDesc.Data()));

  DrawFinalObservable(RmsOutOverIn_AS,RmsOutOverIn_AS_SysError,RmsOutOverIn_AS_Models,cFinal,"ASRmsOutOverIn",Form("%s RMS Out/In",sAwaysideDesc.Data()));
  DrawFinalObservable(RmsOutOverIn_NS,RmsOutOverIn_NS_SysError,RmsOutOverIn_NS_Models,cFinal,"NSRmsOutOverIn",Form("%s RMS Out/In",sNearsideDesc.Data()));
  DrawFinalObservable(RmsMidOverIn_AS,RmsMidOverIn_AS_SysError,RmsMidOverIn_AS_Models,cFinal,"ASRmsMidOverIn",Form("%s RMS Mid/In",sAwaysideDesc.Data()));
  DrawFinalObservable(RmsMidOverIn_NS,RmsMidOverIn_NS_SysError,RmsMidOverIn_NS_Models,cFinal,"NSRmsMidOverIn",Form("%s RMS Mid/In",sNearsideDesc.Data()));

  DrawFinalObservable(SigmasOutOverIn_AS,SigmasOutOverIn_AS_SysError,SigmasOutOverIn_AS_Models,cFinal,"ASSigmasOutOverIn",Form("%s Sigmas Out/In",sAwaysideDesc.Data()));
  DrawFinalObservable(SigmasOutOverIn_NS,SigmasOutOverIn_NS_SysError,SigmasOutOverIn_NS_Models,cFinal,"NSSigmasOutOverIn",Form("%s Sigmas Out/In",sNearsideDesc.Data()));
  DrawFinalObservable(SigmasMidOverIn_AS,SigmasMidOverIn_AS_SysError,SigmasMidOverIn_AS_Models,cFinal,"ASSigmasMidOverIn",Form("%s Sigmas Mid/In",sAwaysideDesc.Data()));
  DrawFinalObservable(SigmasMidOverIn_NS,SigmasMidOverIn_NS_SysError,SigmasMidOverIn_NS_Models,cFinal,"NSSigmasMidOverIn",Form("%s Sigmas Mid/In",sNearsideDesc.Data()));

  printf("Done drawing the final results\n");
}

void TaskCalcObservables::DrawFinalObservable(TGraphErrors * fObsGraph, TGraphErrors * fObsGraphSysErrors, vector<TGraphErrors*> fModels, TCanvas * cFinal, TString sFinalName = "", TString sFinalTitle = "") {
  printf("Drawing Final observable plot for %s\n",fObsGraph->GetName());
  if (fSystematicsNames.size() == 0) {
    printf("\tNo systematics provided, so none will be included.\n");
  }
  if (fModelNames.size() == 0) {
    printf("\tNo Models provided, so none will be included.\n");
  }



  double fXPlotMinRange = 1.0;
  double fXPlotMaxRange = 18;

  float fRatioMin = 0.45;
  float fRatioMax = 1.55;

  float fUnZoomYMin = 0;
  float fUnZoomYMax = 2.5;

  TString sName = fObsGraph->GetName();
  bool bIsRatio = sName.Contains("Ratio") || sName.Contains("Over");
  bool bIsYield = sName.Contains("Yield") || sName == "OutOverIn_AS" || sName == "OutOverIn_NS"; // for 

  bool bIsSigma = sName.Contains("Sigma");
  bool bIsRms   = sName.Contains("Rms",TString::kIgnoreCase);

  bool bIsAS = sName.Contains("AS");
  bool bIsNS = sName.Contains("NS");

  if (bIsRms && bIsRatio) {
    fRatioMin = 0.2;
    fRatioMax = 1.8;
  }
  if (bIsSigma && bIsRatio) {
    fRatioMin = 0.4;
    fRatioMax = 1.6;
  }

  if (bIsYield && !bIsRatio) cFinal->SetLogy(1);
  else cFinal->SetLogy(0);

  fObsGraph->SetMarkerColor(kBlack);
  fObsGraph->SetLineColor(kBlack);

  TH1F * hEmptyHist = new TH1F("hEmpty","",20,0,fXPlotMaxRange);
  hEmptyHist->SetLineColor(0);
  hEmptyHist->GetXaxis()->SetTitle(fObsGraph->GetXaxis()->GetTitle());
  hEmptyHist->GetYaxis()->SetTitle(fObsGraph->GetYaxis()->GetTitle());

  hEmptyHist->Draw();
  hEmptyHist->GetXaxis()->SetRangeUser(fXPlotMinRange,fXPlotMaxRange);

  TF1 * fConstantLine = new TF1("constant","[0]",0,18);
  fConstantLine->SetLineColor(kGray+3);
  fConstantLine->SetLineWidth(5);
  fConstantLine->SetLineStyle(kDashed);
  if (bIsRatio) {
    fConstantLine->SetParameter(0,1.);
  }
  if (bIsRatio) {
    fConstantLine->Draw("SAME");
  }


  printf("Drew the empty hist\n");
  TLegend * legTest = new TLegend(0.4,0.7,0.4,0.7);
  legTest->SetHeader("Test","c");

  if (fIsMCGenMode) {
    fObsGraph->SetFillColor(kAzure+1);
    //fObsGraph->SetFillColorAlpha(kAzure+1,0.5);
    fObsGraph->SetMarkerColor(kAzure+1);
    fObsGraph->SetMarkerStyle(kFullSquare);
    fObsGraph->SetFillStyle(1001);
    fObsGraph->SetLineColor(kAzure+1);
    fObsGraph->Draw("p 0 SAME");
    fObsGraph->Draw("[]5 0 SAME");
    //fObsGraph->Draw("p SAME");
    //fObsGraph->Draw("[]5 SAME");
  } else {
    fObsGraph->SetLineColor(kBlack);
    fObsGraph->SetMarkerColor(kBlack);
    fObsGraph->Draw("P SAME");
  }
    fObsGraph->SetLineColor(kBlack);
    fObsGraph->SetMarkerColor(kBlack);
  //fObsGraph->Draw("AP SAME");
  fObsGraph->GetXaxis()->SetRangeUser(fXPlotMinRange,fXPlotMaxRange);
  float legendX1 = 0.45;
  float legendX2 = 0.8;
  float legendY1 = 0.15;
  float legendY2 = 0.3;

  float fAliceLegendX = 0.45;
  float fAliceLegendY = 0.65;

  if (bIsSigma) {
    legendX1 = 0.75;
    legendX2 = 0.875;
    legendY1 = 0.40;
    legendY2 = 0.55;
  }
  if (bIsRatio) {
    if (bIsRms) {
      legendX1 = 0.60;
      legendX2 = 0.85;
      legendY1 = 0.15;
      legendY2 = 0.35;

    } else {
      legendX1 = 0.60;
      legendX2 = 0.85;
      legendY1 = 0.15;
      legendY2 = 0.35;
    }
    if (bIsSigma && bIsNS) {
      legendX1 = 0.60;
      legendX2 = 0.85;
      legendY1 = 0.55;
      legendY2 = 0.75;
    }
  }
  if (bIsYield) {
    if (bIsRatio) {
      // Default and AS values
      legendX1 = 0.4;
      legendX2 = 0.65;
      legendY1 = 0.12;
      legendY2 = 0.3;
      fAliceLegendX = 0.15;
      fAliceLegendY = 0.65;

      if (bIsNS) {
        legendX1 = 0.40;
        legendX2 = 0.40+0.275;
        legendY1 = 0.12;
        legendY2 = 0.3;
        fAliceLegendX = 0.2;
        //fAliceLegendX = 0.15;
        fAliceLegendY = 0.65;
      }

      printf("Setting legend settings for Yield Ratio for graph %s\n",sName.Data());
    } else {
      legendX1 = 0.6;
      legendX2 = 0.875;
      legendY1 = 0.40;
      legendY2 = 0.60;
    }
  }


  TLegend * legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
  if (sFinalTitle != "") legend->SetHeader(sFinalTitle);
  legend->SetTextSize(0.041); // 0.045
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  TLegend * legAlice = DrawAliceLegend(fObsGraph,fAliceLegendX,fAliceLegendY,kAliceLegendWidth,kAliceLegendHeight);

  if (bIsRatio) {
    hEmptyHist->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
  } else {
    
  }

  if (fSystematicsNames.size()>0) {
    if (fIsMCGenMode) {
      fObsGraphSysErrors->Draw("L SAME");
    } else {
      fObsGraphSysErrors->Draw("L[]5 0 SAME");
      //fObsGraphSysErrors->Draw("L[]5 SAME");
    }
  }

  //if (bIsRatio) {
  //  fConstantLine->Draw("SAME");
  //}

  // FIXME for MCGen, draw error bands

  if (fIsMCGenMode) {
    legend->AddEntry(fObsGraph,Form("%s %.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}",sTitle.Data(),PtBins[iPtBin-1],PtBins[iPtBin]),"LEP");
  } else {
    legend->AddEntry(fObsGraph,"Data","LEP");
    legend->AddEntry(fObsGraphSysErrors,"Sys. Uncert.","FE");
  }
  fObsGraph->Draw("SAME PE");
  legend->Draw("SAME");

  cFinal->SetTicks(1,1);

  cFinal->Print(Form("%s/Final%s.pdf",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/Final%s.eps",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/Final%s.png",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/CFiles/Final%s.C",fOutputDir.Data(),fObsGraph->GetName()));

  float fOldYMin = hEmptyHist->GetYaxis()->GetXmin();
  float fOldYMax = hEmptyHist->GetYaxis()->GetXmax();

  if (bIsRatio) {
    hEmptyHist->GetYaxis()->SetRangeUser(fUnZoomYMin,fUnZoomYMax);
    legend->SetX1(legend->GetX1()-0.15);
    legend->SetX2(legend->GetX2()+0.0);
    legend->SetY2(legend->GetY2()+0.05);
    legend->Draw("SAME");
  }

  cFinal->Print(Form("%s/Final%s_UnZoom.pdf",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/Final%s_UnZoom.png",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/Final%s_UnZoom.eps",fOutputDir.Data(),fObsGraph->GetName()));
  cFinal->Print(Form("%s/Final%s_UnZoom.C",fOutputDir.Data(),fObsGraph->GetName()));

  if (bIsRatio) {
    hEmptyHist->GetYaxis()->SetRangeUser(fOldYMin,fOldYMax);
  }

  if (fModelNames.size() > 0) {
    printf("Now to add the Models\n");
    if (bIsRatio) {
      hEmptyHist->GetYaxis()->SetRangeUser(fRatioMin,fRatioMax);
    } else {
      
    }
    hEmptyHist->Draw();
    if (bIsRatio) {
      fConstantLine->Draw("SAME");
    }
    fObsGraph->SetLineColor(kBlack);
    fObsGraph->SetMarkerColor(kBlack);
    fObsGraph->Draw("P SAME");
    //fObsGraph->Draw("AP");
    fObsGraph->GetXaxis()->SetRangeUser(fXPlotMinRange,fXPlotMaxRange);


    if (fSystematicsNames.size()>0) {
     // fObsGraphSysErrors->SetFillStyle(0);
     // fObsGraphSysErrors->SetLineColor(kBlack);
      fObsGraphSysErrors->Draw("L[]5 SAME");
    }


    TLegend * legAlice2 = DrawAliceLegend(fObsGraph,fAliceLegendX,fAliceLegendY,kAliceLegendWidth,kAliceLegendHeight);


    for (int i = 0; i < (int) fModels.size(); i++) {
      fModels[i]->Draw("3 L SAME");
      legend->AddEntry(fModels[i],fModelTitles[i],"LF");

    }

    //if (bIsRatio) {
    //  fConstantLine->Draw("SAME");
    //}

    fObsGraph->Draw("P SAME");

    legend->Draw("SAME");

    if (sName == "") {
      cFinal->Print(Form("%s/Final%sWithModels.pdf",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels.eps",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels.png",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/CFiles/Final%sWithModels.C",fOutputDir.Data(),fObsGraph->GetName()));
      if (bIsRatio) {
        hEmptyHist->GetYaxis()->SetRangeUser(fUnZoomYMin,fUnZoomYMax);
      }

      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.pdf",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.eps",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.png",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.C",fOutputDir.Data(),fObsGraph->GetName()));

    } else {
      cFinal->Print(Form("%s/Final%sWithModels.pdf",fOutputDir.Data(),sFinalName.Data()));
      cFinal->Print(Form("%s/Final%sWithModels.eps",fOutputDir.Data(),sFinalName.Data()));
      cFinal->Print(Form("%s/Final%sWithModels.png",fOutputDir.Data(),sFinalName.Data()));
      cFinal->Print(Form("%s/CFiles/Final%sWithModels.C",fOutputDir.Data(),sFinalName.Data()));
      if (bIsRatio) {
        hEmptyHist->GetYaxis()->SetRangeUser(fUnZoomYMin,fUnZoomYMax);
      }
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.pdf",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.eps",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.png",fOutputDir.Data(),fObsGraph->GetName()));
      cFinal->Print(Form("%s/Final%sWithModels_UnZoom.C",fOutputDir.Data(),fObsGraph->GetName()));

    }
  } else {
    printf("Not adding the models because they don't exist.\n");
  }

  printf("Done drawing the final plot for this observable\n");
  delete hEmptyHist;
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
  if (fASSigmasInc) fOutputFile->Add(fASSigmasInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fASSigmasEP[i]) fOutputFile->Add(fASSigmasEP[i]);
  }


  if (fNSYieldsInc) fOutputFile->Add(fNSYieldsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fNSYieldsEP[i]) fOutputFile->Add(fNSYieldsEP[i]);
  }
  if (fNSRmsInc) fOutputFile->Add(fNSRmsInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fNSRmsEP[i]) fOutputFile->Add(fNSRmsEP[i]);
  }
  if (fNSSigmasInc) fOutputFile->Add(fNSSigmasInc);
  for (int i = 0; i < kNEPBins; i++) {
    if (fNSSigmasEP[i]) fOutputFile->Add(fNSSigmasEP[i]);
  }

  if (OutOverIn_AS) fOutputFile->Add(OutOverIn_AS);
  if (MidOverIn_AS) fOutputFile->Add(MidOverIn_AS);
  if (OutOverIn_NS) fOutputFile->Add(OutOverIn_NS);
  if (MidOverIn_NS) fOutputFile->Add(MidOverIn_NS);

  if (OutOverIn_AS_RPFError) fOutputFile->Add(OutOverIn_AS_RPFError);
  if (MidOverIn_AS_RPFError) fOutputFile->Add(MidOverIn_AS_RPFError);
  if (OutOverIn_NS_RPFError) fOutputFile->Add(OutOverIn_NS_RPFError);
  if (MidOverIn_NS_RPFError) fOutputFile->Add(MidOverIn_NS_RPFError);

/*
  if (RmsOutOverIn_AS_RPFError) fOutputFile->Add(RmsOutOverIn_AS_RPFError);
  if (RmsMidOverIn_AS_RPFError) fOutputFile->Add(RmsMidOverIn_AS_RPFError);
  if (RmsOutOverIn_NS_RPFError) fOutputFile->Add(RmsOutOverIn_NS_RPFError);
  if (RmsMidOverIn_NS_RPFError) fOutputFile->Add(RmsMidOverIn_NS_RPFError);
  if (SigmasOutOverIn_AS_RPFError) fOutputFile->Add(SigmasOutOverIn_AS_RPFError);
  if (SigmasMidOverIn_AS_RPFError) fOutputFile->Add(SigmasMidOverIn_AS_RPFError);
  if (SigmasOutOverIn_NS_RPFError) fOutputFile->Add(SigmasOutOverIn_NS_RPFError);
  if (SigmasMidOverIn_NS_RPFError) fOutputFile->Add(SigmasMidOverIn_NS_RPFError);
*/

  //if (OutOverIn_NS) fOutputFile->Add(OutOverIn_NS);
  //if (MidOverIn_NS) fOutputFile->Add(MidOverIn_NS);

  if (RmsOutOverIn_AS) fOutputFile->Add(RmsOutOverIn_AS);
  if (RmsMidOverIn_AS) fOutputFile->Add(RmsMidOverIn_AS);
  if (RmsOutOverIn_NS) fOutputFile->Add(RmsOutOverIn_NS);
  if (RmsMidOverIn_NS) fOutputFile->Add(RmsMidOverIn_NS);

  if (SigmasOutOverIn_AS) fOutputFile->Add(SigmasOutOverIn_AS);
  if (SigmasMidOverIn_AS) fOutputFile->Add(SigmasMidOverIn_AS);
  if (SigmasOutOverIn_NS) fOutputFile->Add(SigmasOutOverIn_NS);
  if (SigmasMidOverIn_NS) fOutputFile->Add(SigmasMidOverIn_NS);

  if (OutMinusIn_AS) fOutputFile->Add(OutMinusIn_AS);
  if (MidMinusIn_AS) fOutputFile->Add(MidMinusIn_AS);
  if (OutMinusIn_NS) fOutputFile->Add(OutMinusIn_NS);
  if (MidMinusIn_NS) fOutputFile->Add(MidMinusIn_NS);


  bool SaveVariants = false;
    // Temporary
    if (SaveVariants) {
    for (auto graph : OutOverIn_AS_RPFVariants) {
      fOutputFile->Add(graph);
    }
    for (auto graph : OutOverIn_NS_RPFVariants) {
      fOutputFile->Add(graph);
    }
    for (auto list : fASYieldsEP_RPFVariants) {
      for (auto graph : list) {
        fOutputFile->Add(graph);
      }
    }
  }



  fOutputFile->Write();

}


void TaskCalcObservables::Run() {

  LoadHistograms();

  LoadSystematics();

  LoadModels();

  InitArrays();

  DrawDirectComparisons();

  CalculateResults();

  CleanResults();

  AddTrackingUncertainties();

  CalculateSystematics();

  ProduceSystematicUncertPlots();
  printf("Done with producing systematic uncertainty plots\n");

  DrawResults();

  CombineStatisticalErrors();

  DrawFinalResults();

  SaveOutput();

  cout<<"Done!"<<endl;
}

