
  #include <iostream>

#include "PionID.h"

#include "UserUtilities.cxx"

using std::vector;

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;
/// \cond CLASSIMP
ClassImp(PionID);


PionID::PionID():TObject(),sOutputDir("output")

//fObservable(-1),fObservableName(),nObsBins(0),fObsBins(),
//fFullDPhiPi0(), fFullDPhiSB(), fMassPtBinPi0(), fMassPtBinAll(), fMassPtBinSB(),
//fFullDPhiFinal(),
//fBackgroundSelection(0)
{
  fDebugLevel = 0;

  SetStyle();
}




//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PionID::SetStyle() {
//  gStyle->SetCanvasColor(kBlack);
//  gStyle->SetAxisColor(0);
//  TGaxis::SetMaxDigits(2);
 // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);

  gStyle->SetOptTitle(0); // New
  gStyle->SetOptStat(0);

 // gStyle->SetTitleOffset(0.6,"X"); // 0.7 
  gStyle->SetTitleOffset(0.6,"Y"); // 0.7 
//  gStyle->SetLabelSize(0.04,"X");
//  gStyle->SetTitleSize(0.01,"all"); //0.055
//  gStyle->SetTitleSize(0.045,"X"); // 0.05
//  gStyle->SetPadTopMargin(0.07);//0.05 //..ELI
//  gStyle->SetPadBottomMargin(0.18);//0.15
//  gStyle->SetPadRightMargin(0.04);
//  gStyle->SetPadLeftMargin(0.21);
  TGaxis::SetMaxDigits(3); //

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.10); //0.15

}

void PionID::DrawWIP(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size) {

  TLegend * leg = new TLegend(x,y,x+x_size,y+y_size);
//  leg->SetHeader("ALICE");
//  leg->AddEntry(Histo,"ALICE","");
  leg->AddEntry(Histo,"Work in Progress","");
  leg->AddEntry(Histo,"Pb-Pb #sqrt{s_{NN}} = 5.02 TeV","");
  if (sLabel.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel.Data()),"");
  if (sLabel2.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel2.Data()),"");
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(10,0);

  leg->Draw("SAME");
}

void PionID::DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size) {

  TLegend * leg = new TLegend(x,y,x+x_size,y+y_size);
  TDatime * time = new TDatime();
  string month = kMonthList[time->GetMonth()-1];
//  const char * month = kMonthList[time->GetMonth()-1];

  //leg->SetHeader(Form("ALICE Performance - %d %s %d",2,month,time->GetYear()));
  leg->SetHeader(Form("ALICE Work in Progress - %d %s %d",time->GetDay(),month.c_str(),time->GetYear()));
//  leg->AddEntry(Histo,"ALICE Performance","");
 // leg->AddEntry(Histo,Form("ALICE Performance %d %s %d",time->GetDay(),month.c_str(),time->GetYear()),"");
  leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
  if (sLabel.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel.Data()),"");
  if (sLabel2.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel2.Data()),"");
//  leg->AddEntry(Histo,"PbPb #sqrt{s_{NN}} = 5.02 TeV, EGA","");
  leg->SetTextSize(0.035); // 0.045 0.04
  leg->SetBorderSize(0);
//  leg->SetFillColorAlpha(10,0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
}

//________________________________________________________________________
void PionID::SetTH1Histo(TH1 * hist,TString xTitle,TString yTitle,Bool_t big)
{
  hist->SetStats(0);
  hist->SetTitle("");

  if(big==1)
  {
    hist->GetYaxis()->SetTitleOffset(0.76); //1.0
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(0.76); //0.82
    hist->GetXaxis()->SetLabelSize(0.04); //0.07
    hist->GetXaxis()->SetTitleSize(0.05);
  }
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  //..make nice font
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  if(xTitle!="")hist->GetXaxis()->SetTitle(xTitle);
  if(xTitle!="")hist->GetYaxis()->SetTitle(yTitle);

  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
}

//________________________________________________________________________
void PionID::AddFunctionToHist(TH1 * Hist, TF1 * func) {
  int nBins = Hist->GetNbinsX();

  for (int i = 1; i < nBins+1; i++) {
    double localX = Hist->GetXaxis()->GetBinCenter(i);
    if (localX >= func->GetXmin() && localX <= func->GetXmax()) {
      double localFuncValue = func->Eval(localX);
      Hist->SetBinContent(i,localFuncValue + Hist->GetBinContent(i));
    }
  }
}

//_______________________________
void PionID::CopyTF1Details(TF1 * target, TF1 * source) {
  int nParams = source->GetNpar();
  if (nParams > target->GetNpar()) {
    fprintf(stderr,"Error: trying to copy TF1 details to incorrectly initialized TF1.");
    return;
  }
  for (int i = 0; i < nParams; i++) {
    target->SetParameter(i,source->GetParameter(i));
    target->SetParError(i,source->GetParError(i));
    target->SetParName(i,source->GetParName(i));
  }
}

/*
void PionID::CopyTH2FEntries(TH2F * target, TH2F * source) {
  int nBinsX = target->GetXaxis()->GetNbins();
  int nBinsY = target->GetYaxis()->GetNbins();
  for (int i = 0; i <= nBinsX; i++) {
    for (int j = 0; j <= nBinsY; j++) {
 
    }
  }
}*/

TH1D * PionID::IntegralHist(TH1D * hist) {
  TH1D * intHist = (TH1D *) hist->Clone(TString::Format("%sInt",hist->GetName()));
  intHist->SetTitle(TString::Format("%s Integral",hist->GetName()));
  int nBins = hist->GetNbinsX();
  double integral = 0.0;
  for (int i = 1; i <= nBins; i++) {
    integral += TMath::Abs(hist->GetBinContent(i)); // should work even for negative vals
    intHist->SetBinContent(i,integral);
  }

  return intHist;
}

std::vector <double> PionID::FindEqualBins(TH1D * intHist, int nBins = 10) {
//  std::vector <double> bins = {intHist->GetXaxis()->GetBinLowEdge(1)};
  std::vector <double> bins = {};
  equalBinHeights = {};
  int nHistBins = intHist->GetNbinsX();
  double integral = intHist->GetBinContent(nHistBins); //last bin
  int j = 1;
  double last = 0;

  printf("Finding Equalish Bins ...\n");

  // Finding the first bin edge
  while(intHist->GetBinContent(j) == 0) j++;
  printf(" { %d,",j);
  bins.push_back(intHist->GetXaxis()->GetBinLowEdge(j));
  last = intHist->GetBinContent(j);
  equalBinHeights.push_back(last);
  j++;
  for (int i = 1; i < nBins; i++) {
    double target = last + (integral - last) / (nBins - i + 1); 
    while (intHist->GetBinContent(j) < target) {
      j++;
    }
    printf(" %d,",j);
    bins.push_back(intHist->GetXaxis()->GetBinUpEdge(j));
    equalBinHeights.push_back(intHist->GetBinContent(j));
    last = intHist->GetBinContent(j);
    j++;
  }
  printf(" %d}\n",nHistBins);
  bins.push_back(intHist->GetXaxis()->GetBinUpEdge(nHistBins));
  equalBinHeights.push_back(integral);
  return bins;
}

//  0 to choose underflow bin, 1 to choose first bin
// -1 to choose max bin, -2 to choose overflow
//// bkgChoice = 0 for Rotational background, 1 for Mixed Event background
// bkgChoice = 0 No background
//             1 Mixed Event background (Free Scale)
//             2 Rotational background (Free Scale)
//             3 Mixed Event background (Fixed Scale) 
//             4 Rotational background (Fixed Scale)
//             5 Pos. Swap background (Free Scale)
//             6 Pos. Swap background (Fixed Scale)



void PionID::PrintSettings() {
  cout<<"Printing Settings"<<endl;
  printf("\tOutputFile: %s\n",sOutputFileName.Data());
  printf("\tLambda Range Bins: %d %d\n",LambdaBinLow,LambdaBinHigh);
  printf("\tEnergy Range Bins: %d %d\n",EnergyBinLow,EnergyBinHigh);
  printf("\tAsym   Range Bins: %d %d\n",AsymBinLow,AsymBinHigh);
  printf("\tTheta  Range Bins: %d %d\n",OpeningAngleBinLow,OpeningAngleBinHigh);


}

/// Return 0 on failure
Bool_t PionID::LoadHistograms() {
  cout<<"Loading Histograms"<<endl;
  if (!fInputFile) {
    fprintf(stderr,"Missing Input File\n");
    return 1;
  }

  // Directory would be made here, but now handled by python

  // Finding the Pi0H ME list
  TString type = "Pi0";

  if (fSelectTriggerType == 1) {
    type = "MB";
  } else if (fSelectTriggerType == 2) {
    type = "GA";
  }

  if (bkgChoice == 0)                   { bkgType = 3; }  // Default to ME, but need to fix scale to 0 later.
  if (bkgChoice == 1 || bkgChoice == 3) { bkgType = 3; } // Mixed Event
  if (bkgChoice == 2 || bkgChoice == 4) { bkgType = 2; } // Rotational Background
  if (bkgChoice == 5 || bkgChoice == 6) { bkgType = 4; } // Pos Swapped

  if (bkgChoice == 3 || bkgChoice == 4 || bkgChoice == 6) { fitBkgSub = true; }

  TList * keys = fInputFile->GetListOfKeys();
  TObject * obj;
  TIter next(keys);
  if (sListName.Length() > 0) {
    type = sListName;
  }

  TString listName = "AliAnalysisTask_Pi0H_ME_tracks_caloClusters_histos"; // The Default list name
  while ((obj = next())) { // extra parentheses to keep root 6 from complaining 
    TString name = obj->GetName();
    printf("Found object: %s\n",name.Data());
    if (name.Contains(type.Data())) {
      listName = name;
    }
  }

  printf("Using List %s\n",listName.Data());

  TList * HistoList = (TList * ) fInputFile->Get(listName);
  if (!HistoList) {
    fprintf(stderr,"List %s not found!\n",listName.Data());
    return 1;
  }
  printf("Found List %s\n",HistoList->GetName());

  // Gleaning some information from the list name
  if (listName.Contains("TaskMB")) iThetaModelTrigger = 0;
  if (listName.Contains("TaskGA")) iThetaModelTrigger = 1;
  if (listName.Contains("TaskMC")) iThetaModelTrigger = 2;

  for (int i = 0; i < 4; i++) {
    if (listName.Contains(Form("Cent%d",i))) {
      iThetaModelCent = i;
      iMCPreAnalysis_Cent = i;
      break;
    }
  }

  fHistCentrality = (TH1F *) HistoList->FindObject("fHistCentrality");
  fHistEventHash = (TH1F *) HistoList->FindObject("HistEventHash");

  fClusEnergyMatchedTracks = (TH2F*) HistoList->FindObject("ClusterEnergyMatchedTracks");

  nOpeningAngleBinHigh = OpeningAngleBinHigh;
  nLambdaBinHigh = LambdaBinHigh;
  nEnergyBinHigh = EnergyBinHigh;

  // Extracting the THnSparse
  Pi0Cands = (THnSparse* ) HistoList->FindObject("Pi0Cands");
  if (Pi0Cands) {
    // Histograms of Interest
    // FIXME check if it is the new version.
    haveRotBkg = Pi0Cands->GetNdimensions() >= 7;
    printf("Checking for MC Info ...\n");
    // Iterate over axes
    for (Int_t i = 0; i < Pi0Cands->GetNdimensions(); i++) {
      if (strstr(Pi0Cands->GetAxis(i)->GetTitle(),"MC") != NULL) {
        haveMCStatus = true;
        iMCAxis = i;
      }
      if (strstr(Pi0Cands->GetAxis(i)->GetTitle(),"Patch") != NULL) {
        havePatchStatus = true;
        iPatchStatusAxis = i;
      }
    }
    printf("Patch Candidate Info Availability: %d\n",havePatchStatus);

    if (haveRotBkg) {
      haveMEBkg = Pi0Cands->GetAxis(6)->GetNbins() >= 3;
      if (!haveMEBkg) {
        printf("Pi0Cands THnSparse does not have Mixed Event cluster pairs.\n");
        if (bkgType == 3) {
          fprintf(stderr,"Error: Mixed Event Clusters requested from Root file without them!!.  Quitting ...\n");
          return 1;
        }
      }
      havePSBkg = Pi0Cands->GetAxis(6)->GetNbins() >= 4;
      if (!havePSBkg && bkgType == 4) {
          fprintf(stderr,"Error: Position Swapped Clusters requested from Root file without them!!.  Quitting ...\n");
          return 1;
      }
      // Getting Rotational Background
      TH2F * hMaxLambdaVsMass = (TH2F *) Pi0Cands->Projection(3,1);

      //FIXME check if this should be done sooner
      Pi0Cands->GetAxis(1)->SetRange(1.,Pi0Cands->GetAxis(1)->GetNbins());

      TCanvas * cCompareOpeningAngle = new TCanvas("cCompareOpeningAngle");

      Pi0Cands->GetAxis(6)->SetRange(bkgType,bkgType);
      TH1F * RotOpeningAngle = (TH1F *) Pi0Cands->Projection(2,"e");
      RotOpeningAngle->SetName("RotOpeningAngle");

      Pi0Cands->GetAxis(6)->SetRange(1,1);
      TH1F * SignalOpeningAngle = (TH1F *) Pi0Cands->Projection(2,"e");
      SignalOpeningAngle->SetName("SignalOpeningAngle");
      RotOpeningAngle->SetLineColor(kRed);
      RotOpeningAngle->SetMarkerColor(kRed);
      RotOpeningAngle->SetMarkerStyle(kOpenSquare);
      SignalOpeningAngle->SetMarkerStyle(kFullSquare);
      SignalOpeningAngle->Draw();

      if (RotOpeningAngle->GetEntries() > 0) {
        RotOpeningAngle->Scale(SignalOpeningAngle->GetEntries()/RotOpeningAngle->GetEntries());
      }
      RotOpeningAngle->Draw("SAME");

      fInvarMasspTNoCuts = (TH2F *) Pi0Cands->Projection(0,1);
      fInvarMasspTNoCuts->SetName("fInvarMasspTNoCuts");
      // Cuts
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

      //  0 to choose underflow bin, 1 to choose first bin
      // -1 to choose max bin, -2 to choose overflow

      // Apply Lambda Cut
      nLambdaBinHigh = LambdaBinHigh;
      if (nLambdaBinHigh == -1) nLambdaBinHigh = Pi0Cands->GetAxis(3)->GetNbins();
      if (nLambdaBinHigh == -2) nLambdaBinHigh = Pi0Cands->GetAxis(3)->GetNbins() + 1;
      printf("Lambda Range: [%.2f - %.2f] \t\tBins: %d %d\n",Pi0Cands->GetAxis(3)->GetBinLowEdge(LambdaBinLow),Pi0Cands->GetAxis(3)->GetBinUpEdge(nLambdaBinHigh),LambdaBinLow,nLambdaBinHigh);
      Pi0Cands->GetAxis(3)->SetRange(LambdaBinLow,nLambdaBinHigh);

      // Apply Min Energy Cut
      nEnergyBinHigh = EnergyBinHigh;
      if (nEnergyBinHigh == -1) nEnergyBinHigh = Pi0Cands->GetAxis(4)->GetNbins();
      if (nEnergyBinHigh == -2) nEnergyBinHigh = Pi0Cands->GetAxis(4)->GetNbins() + 1;
      printf("Energy Range: [%.2f - %.2f] \t\tBins: %d %d\n",Pi0Cands->GetAxis(4)->GetBinLowEdge(EnergyBinLow),Pi0Cands->GetAxis(4)->GetBinUpEdge(nEnergyBinHigh),EnergyBinLow,nEnergyBinHigh);
      Pi0Cands->GetAxis(4)->SetRange(EnergyBinLow,nEnergyBinHigh);

      // Apply Asymmetry Cut
      int nAsymBinHigh = AsymBinHigh;
      if (nAsymBinHigh == -1) nAsymBinHigh = Pi0Cands->GetAxis(5)->GetNbins();
      if (nAsymBinHigh == -2) nAsymBinHigh = Pi0Cands->GetAxis(5)->GetNbins() + 1;
      printf("Asym Range:   [%.2f - %.2f] \t\tBins: %d %d\n",Pi0Cands->GetAxis(5)->GetBinLowEdge(AsymBinLow),Pi0Cands->GetAxis(5)->GetBinUpEdge(nAsymBinHigh),AsymBinLow,nAsymBinHigh);
      Pi0Cands->GetAxis(5)->SetRange(AsymBinLow,nAsymBinHigh);

      // Apply Event Plane Bin Cut
      
      printf("EvtPl Range:  %d - %d\n",EventPlaneBinLow,EventPlaneBinHigh);
      Pi0Cands->GetAxis(EventPlaneAxis)->SetRange(EventPlaneBinLow,EventPlaneBinHigh);

      // Apply Opening Angle Cut
      nOpeningAngleBinHigh = OpeningAngleBinHigh;
      fOpeningAngleCut = Pi0Cands->GetAxis(2)->GetBinLowEdge(OpeningAngleBinLow);
      //  int nOpeningAngleBinHigh = OpeningAngleBinHigh;
      if (nOpeningAngleBinHigh == -1) nOpeningAngleBinHigh = Pi0Cands->GetAxis(2)->GetNbins();
      if (nOpeningAngleBinHigh == -2) nOpeningAngleBinHigh = Pi0Cands->GetAxis(2)->GetNbins() + 1;
      printf("OpeningAngle Range:   [%.3f - %.3f] \t\tBins: %d %d\n",Pi0Cands->GetAxis(2)->GetBinLowEdge(OpeningAngleBinLow),Pi0Cands->GetAxis(2)->GetBinUpEdge(nOpeningAngleBinHigh),OpeningAngleBinLow,nOpeningAngleBinHigh);
      Pi0Cands->GetAxis(2)->SetRange(OpeningAngleBinLow,nOpeningAngleBinHigh);

      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");


      // Setting for mixed background
      Pi0Cands->GetAxis(6)->SetRange(bkgType,bkgType);

      if (restrictToPatch && havePatchStatus) {
        Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(2,2);
        fInvarMasspTRotBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");
        printf("DEBUG RotBkg entries = %f\n",fInvarMasspTRotBkg->GetEntries());
        Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(1,Pi0Cands->GetAxis(iPatchStatusAxis)->GetNbins());
      } else {
        fInvarMasspTRotBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");

      }

      fInvarMasspTRotBkg->SetName("hInvarMasspTRotBkg");
      fInvarMasspTRotBkg->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background");

      // Adding MC RotBkg plots
      if (haveMCStatus) {
        // TH2F for real pi0 origin (valid for pos swap)
        // FIXME make this work for the new labelling method
        int nPi0Bin = 2;
        Pi0Cands->GetAxis(iMCAxis)->SetRange(nPi0Bin,nPi0Bin);
        fInvarMasspTRotBkgMCPi0 = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fInvarMasspTRotBkgMCPi0->SetName("hInvarMasspTRotBkgMCPi0");
        fInvarMasspTRotBkgMCPi0->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background, Pi0 Sources*");

        if (Pi0Cands->GetAxis(iMCAxis)->GetNbins()>11) {
//          Pi0Cands->GetAxis(iMCAxis)->SetRange(10,10);
          Pi0Cands->GetAxis(iMCAxis)->SetRange(11,11);
          fInvarMasspTRotBkgMCPi0EnergyPair = (TH2F *) Pi0Cands->Projection(0,1,"e");
          fInvarMasspTRotBkgMCPi0EnergyPair->SetName("hInvarMasspTRotBkgPi0EnergyPair");
          fInvarMasspTRotBkgMCPi0EnergyPair->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background, Same Energy Pi0 Source");
//          Pi0Cands->GetAxis(iMCAxis)->SetRange(11,11);
          Pi0Cands->GetAxis(iMCAxis)->SetRange(12,12);
          fInvarMasspTRotBkgMCPi0PosPair = (TH2F *) Pi0Cands->Projection(0,1,"e");
          fInvarMasspTRotBkgMCPi0PosPair->SetName("hInvarMasspTRotBkgPi0PosPair");
          fInvarMasspTRotBkgMCPi0PosPair->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background, Same Position Pi0 Source");
//          Pi0Cands->GetAxis(iMCAxis)->SetRange(12,12);
          Pi0Cands->GetAxis(iMCAxis)->SetRange(13,13);
          fInvarMasspTRotBkgMCEtaEnergyPair = (TH2F *) Pi0Cands->Projection(0,1,"e");
          fInvarMasspTRotBkgMCEtaEnergyPair->SetName("hInvarMasspTRotBkgEtaEnergyPair");
          fInvarMasspTRotBkgMCEtaEnergyPair->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background, Same Energy Eta Source");
//          Pi0Cands->GetAxis(iMCAxis)->SetRange(13,13);
          Pi0Cands->GetAxis(iMCAxis)->SetRange(14,14);
          fInvarMasspTRotBkgMCEtaPosPair = (TH2F *) Pi0Cands->Projection(0,1,"e");
          fInvarMasspTRotBkgMCEtaPosPair->SetName("hInvarMasspTRotBkgEtaPosPair");
          fInvarMasspTRotBkgMCEtaPosPair->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background, Same Position Eta Source");

          // Subtraction of MC Pi0 Same-Position Signal
          fInvarMasspTRotBkgMinusMCPi0Pos = (TH2F *) fInvarMasspTRotBkg->Clone("fInvarMasspTRotBkgMinusMCPi0Pos");
          fInvarMasspTRotBkgMinusMCPi0Pos->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background Minus True Same-Position Pi0 PS Peak");
          fInvarMasspTRotBkgMinusMCPi0Pos->Add(fInvarMasspTRotBkgMCPi0PosPair,-2); // the MC histogram is half-weighted

          // Subtraction of MC PS All
          fInvarMasspTRotBkgMinusMCAll = (TH2F *) fInvarMasspTRotBkg->Clone("fInvarMasspTRotBkgMinusMCAll");
          fInvarMasspTRotBkgMinusMCAll->SetTitle("#it{p}_{T} vs M_{#gamma#gamma} Rotated Background Minus All PosSwapped #pi^{0},#eta");
          fInvarMasspTRotBkgMinusMCAll->Add(fInvarMasspTRotBkgMCPi0PosPair,-2); // the MC histogram is half-weighted
          fInvarMasspTRotBkgMinusMCAll->Add(fInvarMasspTRotBkgMCPi0EnergyPair,-2); // the MC histogram is half-weighted
          fInvarMasspTRotBkgMinusMCAll->Add(fInvarMasspTRotBkgMCEtaPosPair,-2); // the MC histogram is half-weighted
          fInvarMasspTRotBkgMinusMCAll->Add(fInvarMasspTRotBkgMCEtaEnergyPair,-2); // the MC histogram is half-weighted
        }

        //TH2F for the real eta origin
        int nEtaLowBin = 4;
        int nEtaHighBin = 7;
        Pi0Cands->GetAxis(iMCAxis)->SetRange(nEtaLowBin,nEtaHighBin);
        fInvarMasspTRotBkgMCEta = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fInvarMasspTRotBkgMCEta->SetName("hInvarMasspTRotBkgMCEta");

        if (RemoveMCPi0PS) { // Trying to remove Pos Swapped pi0 peak
          fInvarMasspTRotBkg->Add(fInvarMasspTRotBkgMCPi0,-2);
        }
        if (RemoveMCEta) { // Removing eta contribution
          fInvarMasspTRotBkg->Add(fInvarMasspTRotBkgMCEta,-2); // -2 scale to account for half-weighting
//          fInvarMasspTRotBkg->Add(fInvarMasspTRotBkgMCEta,-1);
        }

        // Return to default:
        Pi0Cands->GetAxis(iMCAxis)->SetRange(1,Pi0Cands->GetAxis(iMCAxis)->GetNbins());
      }

      if (havePatchStatus) {
        Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(2,2);
        fInvarMassPtRotBkgGAPatch = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fInvarMassPtRotBkgGAPatch->SetName("hInvarMassPtRotBkgGAPatch");
        Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(1,1);
        fInvarMassPtRotBkgNoPatch = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fInvarMassPtRotBkgNoPatch->SetName("hInvarMassPtRotBkgNoPatch");
        Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(1,2);
      }

      hPairPtRotBkg = (TH1D *) Pi0Cands->Projection(0,"e");
      hPairPtRotBkg->SetName("hPairPtRotBkg");

      hPairOpeningAngleRotBkg = (TH1D *) Pi0Cands->Projection(2,"e");
      hPairOpeningAngleRotBkg->SetName("hPairOpeningAngleRotBkg");

      hPairPtOpAngleBkg = (TH2F *) Pi0Cands->Projection(2,0,"e");
      hPairPtOpAngleBkg->SetName("hPairPtOpAngleBkg");

      hPtAngleDetailRotBkg = (TH2F *) Pi0Cands->Projection(2,0,"e");
      hPtAngleDetailRotBkg->SetName("hPtAngleDetailRotBkg");

      hAngleEnergyCutRotBkg = (TH2F *) Pi0Cands->Projection(2,4,"e");
      hAngleEnergyCutRotBkg->SetName("hAngleEnergyCutRotBkg");


      // Returning to SE
      Pi0Cands->GetAxis(6)->SetRange(1,1); // not rotated.

    }

    // temporary switch FIXME
    if (restrictToPatch && havePatchStatus) {
      if (!restrictPatchOnlyInME)
      Pi0Cands->GetAxis(iPatchStatusAxis)->SetRange(2,2);
    }


    fInvarMasspT = (TH2F *) Pi0Cands->Projection(0,1,"e");
    fInvarMasspT->SetName("fInvarMasspT");
    fInvarMasspT->GetXaxis()->SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})");

    hPairPt = (TH1D *) Pi0Cands->Projection(0);
    hPairPt->SetName("hPairPt");

    if (haveMCStatus) {

      //Building histogram without pi->2gamma, eta->2gamma peaks
      Pi0Cands->GetAxis(iMCAxis)->SetRange(1,1);
      // Start with NoMatch, add other bkg later
      fInvarMassPtMCNoPeak = (TH2F *) Pi0Cands->Projection(0,1,"e");
      fInvarMassPtMCNoPeak->SetName("fInvarMassPtMCNoPeak");
      fInvarMassPtMCNoPeak->GetXaxis()->SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})");

      fInvarMassPtMCNoEta = (TH2F *) Pi0Cands->Projection(0,1,"e");
      fInvarMassPtMCNoEta->SetName("fInvarMassPtMCNoEta");
      fInvarMassPtMCNoEta->GetXaxis()->SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})");
      fInvarMassPtMCNoEtaBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");
      fInvarMassPtMCNoEtaBkg->SetName("fInvarMassPtMCNoEtaBkg");
      fInvarMassPtMCNoEtaBkg->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/#it{c}^{2})");

      //fInvarMasspTMCId
      for (Int_t k = 0; k < nMCId; k++) {
        Pi0Cands->GetAxis(iMCAxis)->SetRange(k+1,k+1);
        fInvarMassPtMCId.push_back((TH2F *) Pi0Cands->Projection(0,1,"e"));
        //fInvarMasspTMCId[k]->SetName(TString::Format("fInvarMasspTMC%s",sMCIdNames[k]));
        fInvarMassPtMCId[k]->SetName(Form("fInvarMasspTMC%s",sMCIdNames[k].Data()));
        fInvarMassPtMCId[k]->SetTitle(Form("#it{p}_{T} vs M_{#gamma#gamma} %s",sMCIdTitles[k].Data()));
        fInvarMassPtMCId[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/#it{c}^{2})");

        // Adding to NoPeak Histogram
        if ( k != 0 && k != 2 && k!= 4) { //No Peak already has NoMatch.  Don't add 2gamma modes
          fInvarMassPtMCNoPeak->Add(fInvarMassPtMCId[k]);
        }
        // Adding to No Eta Histogram
        if ( k != 0 && k != 4 && k != 5 && k != 6 && k != 7) {
          fInvarMassPtMCNoEta->Add(fInvarMassPtMCId[k]);
        }
        // Adding to No Eta Bkg Histogram (keep eta->2g peak)
        if ( k != 0 && k != 5 && k != 6 && k != 7) {
          fInvarMassPtMCNoEtaBkg->Add(fInvarMassPtMCId[k]);
        }
      }

      if (RemoveMCEta) { // Replace fInvarMasspT entries
        delete fInvarMasspT;
        fInvarMasspT = (TH2F *) fInvarMassPtMCNoEtaBkg->Clone("fInvarMasspT");
      }

      Pi0Cands->GetAxis(iMCAxis)->SetRange(1,Pi0Cands->GetAxis(iMCAxis)->GetNbins());
    }

    hPairOpeningAngle = (TH1D *) Pi0Cands->Projection(2);
    hPairOpeningAngle->SetName("hPairOpeningAngle");
    hPairOpeningAngle->Sumw2();

    if (useLambdaBkgOpenAngle){
      Pi0Cands->GetAxis(3)->SetRange(nHighLambdaBin,nHighLambdaBinMax);
      hPairPtOpAngle = (TH2F *) Pi0Cands->Projection(2,0,"e");
      Pi0Cands->GetAxis(3)->SetRange(nHighLambdaBin,nHighLambdaBinMax);
      Pi0Cands->GetAxis(3)->SetRange(LambdaBinLow,nLambdaBinHigh);
    } else {
      hPairPtOpAngle = (TH2F *) Pi0Cands->Projection(2,0,"e");
    }
    hPairPtOpAngle->SetName("hPairPtOpAngle");

    hPtAngleDetail = (TH2F *) Pi0Cands->Projection(2,0,"e");
    hPtAngleDetail->SetName("hPtAngleDetail");

    hAngleEnergyCut = (TH2F *) Pi0Cands->Projection(2,4,"e");
    hAngleEnergyCut->SetName("hAngleEnergyCut");

    // Look at large lambda, hope to avoid signal
    printf("Setting High Lambda range: %d %d\n",nHighLambdaBin,nHighLambdaBinMax);
    Pi0Cands->GetAxis(3)->SetRange(nHighLambdaBin,nHighLambdaBinMax);
    hPairLambdaBkgOpeningAngle = (TH1D *) Pi0Cands->Projection(2,"e");
    hPairLambdaBkgOpeningAngle->SetName("hPairLambdaBkgOpeningAngle");
    printf("Returning to signal Lambda range: %d %d\n",LambdaBinLow,nLambdaBinHigh);
    Pi0Cands->GetAxis(3)->SetRange(LambdaBinLow,nLambdaBinHigh);

    printf("Setting Low Pt Range %d %d\n",nLowPtBinMin,nLowPtBinMax);
    Pi0Cands->GetAxis(3)->SetRange(nLowPtBinMin,nLowPtBinMax);
    hPairPtBkgOpeningAngle = (TH1D *) Pi0Cands->Projection(2,"e");
    hPairPtBkgOpeningAngle->SetName("hPairPtBkgOpeningAngle");
    printf("Returning to full Pt Range %d %d\n",1,Pi0Cands->GetAxis(0)->GetNbins());
    Pi0Cands->GetAxis(0)->SetRange(1,Pi0Cands->GetAxis(0)->GetNbins());
  } else {
    fprintf(stdout,"Pi0Cands THnSparse not found!  Looking for fHistClusPairInvarMassPt\n");
    fInvarMasspT = (TH2F * ) HistoList->FindObject("fHistClusPairInvarMasspT");
    if (fInvarMasspT) {
      hPairPt = (TH1D *) fInvarMasspT->ProjectionY("hPairPt",1,fInvarMasspT->GetNbinsY());
    } else {
      fprintf(stderr,"fHistClusPairInvarMasspT not found!\n");
      return 1;
    }
  }

  // Create AngleScaled AngleCorrection here
  nThetaBins = Pi0Cands->GetAxis(2)->GetNbins();

  fHistEvsPt = (TH2F * ) HistoList->FindObject("fHistEvsPt");
  if (!fHistEvsPt) {
    fprintf(stderr,"Missing fHistEvsPt\n");
  }

  fClusEnergy = (TH1F *) HistoList->FindObject("ClusEnergy");
  if (!fClusEnergy) {
    fprintf(stderr,"Missing ClusEnergy\n");
  }
  //fClusEnergy->SetDirectory(0);

  ClusterProp = (THnSparse* ) HistoList->FindObject("ClusterProp");
  if (!ClusterProp) {
    fprintf(stderr,"Missing Cluster Prop ThnSparse.\n");
    if (!fHistEvsPt) {
      fprintf(stderr,"Missing both sources of cluster E energy info!  Exitting.\n");
      return 1;
    }
  }

  TH2F * fMAngle = (TH2F *) HistoList->FindObject("fMAngle");
  TH2F * fPtAngle = (TH2F *) HistoList->FindObject("fPtAngle");
  if (fMAngle) {
    hPairTheta = (TH1D *) fMAngle->ProjectionY("hPairTheta",1,fMAngle->GetNbinsX());
    hPairTheta->SetTitle("Opening Angle");
  } else if (fPtAngle) {
    hPairTheta = (TH1D *) fPtAngle->ProjectionY("hPairTheta",1,fPtAngle->GetNbinsX());
    hPairTheta->SetTitle("Opening Angle");
  } else {
    printf("Creating default theta distribution (flat)");
    hPairTheta = new TH1D("hPairTheta","Opening Angle",100,0,TMath::Pi());
    for (int i = 1; i <= hPairTheta->GetNbinsX(); i++) {
      hPairTheta->SetBinContent(i,1.);
    }
  }
  hPairTheta->Scale(1./hPairTheta->Integral("width"));
/*
  // Checking for the Position Swap Correction Histograms
  fUDist = (THnSparse* ) HistoList->FindObject("UDist");
  fUTildeDist = (THnSparse* ) HistoList->FindObject("UTildeDist");
  fVDist = (THnSparse* ) HistoList->FindObject("VDist");
  fVTildeDist = (THnSparse* ) HistoList->FindObject("VTildeDist");
  
  if (!fUDist || !fUTildeDist || !fVDist || !fVTildeDist) {
    printf("Position Swap Correction histograms not present.");
  }
*/
  fUScaleMatrix = (THnSparse* ) HistoList->FindObject("UScaleMatrix");
  fVScaleMatrix = (THnSparse* ) HistoList->FindObject("VScaleMatrix");

  if (!fUScaleMatrix || !fVScaleMatrix) {
    printf("Position Swap Correction matrices not found\n");
  }

  //fPSMassPtMap = (THnSparse *) HistoList->FindObject("fPSMassPtMap");
  //FIXME temporarily replaceing the PS map with the ES map
  fPSMassPtMap = (THnSparse *) HistoList->FindObject("fESMassPtMap");
  if (!fPSMassPtMap) { // use the normal one
    fPSMassPtMap = (THnSparse *) HistoList->FindObject("fPSMassPtMap");
  }


  if (!fPSMassPtMap) {
    printf("Direct PS Map not found\n");
  }

  fMatchDeltaPhiTrackPt     = (TH2F *) HistoList->FindObject("fMatchDeltaPhiTrackPt");
  fMatchDeltaEtaTrackPt     = (TH2F *) HistoList->FindObject("fMatchDeltaEtaTrackPt");
  fMatchCondDeltaPhiTrackPt = (TH2F *) HistoList->FindObject("fMatchCondDeltaPhiTrackPt");
  fMatchCondDeltaEtaTrackPt = (TH2F *) HistoList->FindObject("fMatchCondDeltaEtaTrackPt");

  // Loading Delta Psi histograms
  hPtEPAnglePionAcc = (TH2F *) HistoList->FindObject("PtEPAnglePionAcc");
  if (hPtEPAnglePionAcc) hPtEPAnglePionAcc->Rebin2D(nRebinDeltaPsi,1);
  hPtEPAngleMCPion = (TH2F *) HistoList->FindObject("PtEPAngleMCPion");
  if (hPtEPAngleMCPion) hPtEPAngleMCPion->Rebin2D(nRebinDeltaPsi,1);
  hPtEPAngleTrueRecMCPion = (TH2F *) HistoList->FindObject("PtEPAngleTrueRecMCPion");
  if (hPtEPAngleTrueRecMCPion) hPtEPAngleTrueRecMCPion->Rebin2D(nRebinDeltaPsi,1);

  // 3rd and 4th order
  hPtEP3AnglePionAcc = (TH2F *) HistoList->FindObject("PtEP3AnglePionAcc");
  if (hPtEP3AnglePionAcc) hPtEP3AnglePionAcc->Rebin2D(nRebinDeltaPsi,1);
  hPtEP4AnglePionAcc = (TH2F *) HistoList->FindObject("PtEP4AnglePionAcc");
  if (hPtEP4AnglePionAcc) hPtEP4AnglePionAcc->Rebin2D(nRebinDeltaPsi,1);


  hHistTrackPsiEPPtCent = (TH3F *) HistoList->FindObject("fHistTrackPsiEPPtCent");
  hHistTrackPsiEP3PtCent = (TH3F *) HistoList->FindObject("fHistTrackPsiEP3PtCent");
  hHistTrackPsiEP4PtCent = (TH3F *) HistoList->FindObject("fHistTrackPsiEP4PtCent");

  if (hHistTrackPsiEP3PtCent) {
    printf("Successfully found EP3 track histogram.\n");
  }


  hPtRPAnglePionAcc = (TH2F *) HistoList->FindObject("PtRPAnglePionAcc");
  if (hPtRPAnglePionAcc) hPtRPAnglePionAcc->Rebin2D(nRebinDeltaPsi,1);
  hPtRPAngleMCPion = (TH2F *) HistoList->FindObject("PtRPAngleMCPion");
  if (hPtRPAngleMCPion) hPtRPAngleMCPion->Rebin2D(nRebinDeltaPsi,1);
  hPtRPAngleTrueRecMCPion = (TH2F *) HistoList->FindObject("PtRPAngleTrueRecMCPion");
  if (hPtRPAngleTrueRecMCPion) hPtRPAngleTrueRecMCPion->Rebin2D(nRebinDeltaPsi,1);

  hHistTrackPsiRPPtCent = (TH3F *) HistoList->FindObject("fHistTrackPsiRPPtCent");


  return 0;
}

void PionID::LoadMCPreAnalysis() {

  for (Int_t i = 0; i < kNCentBins; i++) {


    TString sMCAnalysisFileName = sMCPreAnalysisFile;
    sMCAnalysisFileName.ReplaceAll("CentN",Form("Cent%d",i));
    printf("   PreAnalyzed MC Result loaded from file %s\n",sMCAnalysisFileName.Data());

    TFile * fMCAnalysisFile = TFile::Open(sMCAnalysisFileName);
    if (!fMCAnalysisFile) {
      printf("Could not find file %s. Disabling loading of parameters from MC.\n",sMCAnalysisFileName.Data());
      bUseMCPreAnalysis = false;
      return;
    }

    vector<TF1 *> fMCPreAnalysisPi0Fits_CentBin = {};
    for (int j = 0; j < 7 ; j++) { //load 4-5,5-7,7-9,9-11,11-14,14-17,17-20 FIXME magic number
      TString sObjectName = Form("MCPi0Mass_Fit_%d",j);
      TF1 * fLocalPi0Fit = (TF1 *) fMCAnalysisFile->Get(sObjectName);
      if (!fLocalPi0Fit) {
        fprintf(stderr,"Could not find MC Pi0 Fit for bin %d. Disabling loading of parameters from MC.\n",j);
        bUseMCPreAnalysis = false;
        return;
      }
      fMCPreAnalysisPi0Fits_CentBin.push_back(fLocalPi0Fit);
    }
    fMCPreAnalysisPi0Fits.push_back(fMCPreAnalysisPi0Fits_CentBin);
  }
}


void PionID::LoadThetaModelParameters() {
  // FIXME load tgraphs from file, make sure they can be browsed easily
  // string sThetaModelRootFile = "/home/moliver/cern/gammaHadron/wrk/phase1/AngleAnalysis/T40/GA/CentN/AngleAnalysis.root"
  // Replace CentN with Cent0,Cent1,Cent2,Cent3

  for (Int_t i = 0; i < kNCentBins; i++) {

    TString sAngleAnalysisFileName = sThetaModelRootFile;
    sAngleAnalysisFileName.ReplaceAll("CentN",Form("Cent%d",i));
    printf("  Theta Model loaded from file %s\n",sAngleAnalysisFileName.Data());

    TFile * fAngleAnalysisFile = TFile::Open(sAngleAnalysisFileName);
    if (!fAngleAnalysisFile) {
      printf("Could not find file %s\n",sAngleAnalysisFileName.Data());

      // Fill in null vectors
      vector<TGraphErrors *> lLocalMassPrimeArray;
      vector<TGraphErrors *> lLocalLambdaArray;
      for (Int_t j = 0; j < iNPtBinsThetaModel; j++) {
        lLocalMassPrimeArray.push_back(0);
        lLocalLambdaArray.push_back(0);
      }
      fThetaMassPrimeGraphs.push_back(lLocalMassPrimeArray);
      fThetaLambdaGraphs.push_back(lLocalLambdaArray);
      fMassPrimePar1Graphs.push_back(0);
      fLambdaPar1Graphs.push_back(0);
      continue;
    }

    vector<TGraphErrors *> lLocalMassPrimeArray;
    vector<TGraphErrors *> lLocalLambdaArray;

    for (Int_t j = 0; j < iNPtBinsThetaModel; j++) {
      TString sMassPrimeGraphName = Form("Param_0_Pt_%d",j);
      TGraphErrors * fMassPrimeGraph = (TGraphErrors *) fAngleAnalysisFile->Get(sMassPrimeGraphName);
      if (!fMassPrimeGraph) {
        printf("Could not find mass prime TGraphErrors %s.\n",sMassPrimeGraphName.Data());
        continue;
      }
      lLocalMassPrimeArray.push_back(fMassPrimeGraph);

      TString sLambdaGraphName = Form("Param_1_Pt_%d",j);
      TGraphErrors * fLambdaGraph = (TGraphErrors *) fAngleAnalysisFile->Get(sLambdaGraphName);
      if (!fLambdaGraph) {
        printf("Could not find  lambda TGraphErrors %s.\n",sLambdaGraphName.Data());
        continue;
      }
      lLocalLambdaArray.push_back(fLambdaGraph);
    }

    fThetaMassPrimeGraphs.push_back(lLocalMassPrimeArray);
    fThetaLambdaGraphs.push_back(lLocalLambdaArray);

    // The Pt Array
    TGraphErrors * fMassPrimePar1Graph = (TGraphErrors *) fAngleAnalysisFile->Get("MassPrimePar1Graph");
    TGraphErrors * fLambdaPar1Graph    = (TGraphErrors *) fAngleAnalysisFile->Get("LambdaPar1Graph");
    fMassPrimePar1Graphs.push_back(fMassPrimePar1Graph);
    fLambdaPar1Graphs.push_back(fLambdaPar1Graph);
  }
}

void PionID::OpeningAngleAnalysis() {
  cout<<"Opening Angle Analysis"<<endl;

  int kOpeningAngleMarkerStyle = kFullSquare;

  TCanvas * cOpeningAngle = new TCanvas("cOpeningAngle");
  cOpeningAngle->cd();
  NormalizeHistByBinWidth(hPairOpeningAngle);
  NormalizeHist(hPairOpeningAngle);

  hPairOpeningAngle->SetMarkerStyle(kOpeningAngleMarkerStyle);
  hPairOpeningAngle->Draw();
  hPairOpeningAngle->GetXaxis()->SetRangeUser(0.009,3.14);
  if (hPairOpeningAngleRotBkg) {

    // normalize these to dN/dTheta
    NormalizeHistByBinWidth(hPairOpeningAngleRotBkg);
    NormalizeHist(hPairOpeningAngleRotBkg);
    NormalizeHistByBinWidth(hPairLambdaBkgOpeningAngle);
    NormalizeHist(hPairLambdaBkgOpeningAngle);
    NormalizeHistByBinWidth(hPairPtBkgOpeningAngle);
    NormalizeHist(hPairPtBkgOpeningAngle);

    hPairOpeningAngleRotBkg->SetLineColor(kRed);
    hPairOpeningAngleRotBkg->SetMarkerColor(kRed);
    hPairOpeningAngleRotBkg->SetMarkerStyle(kOpeningAngleMarkerStyle);
    hPairOpeningAngleRotBkg->Draw("SAME");

    hPairLambdaBkgOpeningAngle->SetLineColor(kViolet);
    hPairLambdaBkgOpeningAngle->SetMarkerColor(kViolet);
    hPairLambdaBkgOpeningAngle->SetMarkerStyle(kOpeningAngleMarkerStyle);
    hPairLambdaBkgOpeningAngle->Draw("SAME");

    hPairPtBkgOpeningAngle->SetLineColor(kGreen);
    hPairPtBkgOpeningAngle->SetMarkerColor(kGreen);
    hPairPtBkgOpeningAngle->SetMarkerStyle(kOpeningAngleMarkerStyle);
//    hPairPtBkgOpeningAngle->Draw("SAME");

    TLegend * lOpeningAngle = new TLegend(0.65,0.55,0.90,0.8);
    lOpeningAngle->AddEntry(hPairOpeningAngle,"Same Event Pairs","lp");
    lOpeningAngle->AddEntry(hPairOpeningAngleRotBkg,"Rotated pairs","lp");
    lOpeningAngle->AddEntry(hPairLambdaBkgOpeningAngle,"#lambda Bkg Opening Angle","lp");
//    lOpeningAngle->AddEntry(hPairPtBkgOpeningAngle,"Low #it{p}_{T} Bkg Opening Angle","lp");

    lOpeningAngle->Draw("SAME");
  }
  cOpeningAngle->SetLogy();
  cOpeningAngle->SetLogx();
  cOpeningAngle->Print(Form("%s/OpeningAngle.pdf",sOutputDir.Data()));
  cOpeningAngle->Print(Form("%s/OpeningAngle.C",sOutputDir.Data()));

  TCanvas * cOpeningAngleCorrection = new TCanvas("cOpeningAngleCorrection");

  const int nPtOpAngleRebin_Pt = 1;
  const int nPtOpAngleRebin_Angle = 1;

  //  hPairPtOpAngle = (TH2F *) Pi0Cands->Projection(2,0,"e");
  if (hPairPtOpAngle) {

    hPairPtOpAngle->Rebin2D(nPtOpAngleRebin_Pt,nPtOpAngleRebin_Angle);
    hPairPtOpAngle->Rebin2D(nPtOpAngleRebin_Pt,nPtOpAngleRebin_Angle);

    h2DOpeningAngleCorr = (TH2F *) hPairPtOpAngle->Clone("h2DOpeningAngleCorr");
    // rebin???
    h2DOpeningAngleCorr->Divide(hPairPtOpAngleBkg);
  //  h2DOpeningAngleCorr->Divide(hPairPtOpAngle);

  }
  if (hPairOpeningAngleRotBkg) {

//      hPairOpeningAngleRotBkg->Sumw2();
//      hPairLambdaBkgOpeningAngle->Sumw2();
//      hPairPtBkgOpeningAngle->Sumw2();

    if (useLambdaBkgOpenAngle) {
      //Using high Lambda region to estimate proper theta dist
      hOpeningAngleCorrection = (TH1D *) hPairLambdaBkgOpeningAngle->Clone("hOpeningAngleCorrection");
    } else if (usePtBkgOpenAngle) {
      //Using low pt region to estimate proper theta dist
      hOpeningAngleCorrection = (TH1D *) hPairPtBkgOpeningAngle->Clone("hOpeningAngleCorrection");
    } else {
      //Using signal region to estimate proper theta dist
      hOpeningAngleCorrection = (TH1D *) hPairOpeningAngle->Clone("hOpeningAngleCorrection");
    }

   // hOpeningAngleCorrection = (TH1D *) hPairOpeningAngle->Clone("hOpeningAngleCorrection");

    hOpeningAngleCorrection->Sumw2();
    hOpeningAngleCorrection->SetTitle("Opening Angle Correction Factor");

    hOpeningAngleCorrection->Divide(hPairOpeningAngleRotBkg);

  }

  hOpeningAngleCorrection->Draw();
  cOpeningAngleCorrection->SetLogy(0);
  cOpeningAngleCorrection->Print(Form("%s/OpeningAngleRatio.C",sOutputDir.Data()));
  cOpeningAngleCorrection->Print(Form("%s/OpeningAngleRatio.pdf",sOutputDir.Data()));

  if (bApplyOpeningAngleCorrection) {
    // Only applying to RotBkg
    Pi0Cands->GetAxis(6)->SetRange(bkgType,bkgType);
   // for (int i = 1; i < nThetaBins; i++ ) { 
    for (int i = 0; i < 1 + OpeningAngleBinHigh - OpeningAngleBinLow; i++ ) {
    //for (int i = 0; i < 1 + nOpeningAngleBinHigh - OpeningAngleBinLow; i++ ) {
   // for (int i = OpeningAngleBinLow; i < nThetaBins; i++ ) { 
//    for (int i = OpeningAngleBinLow; i < nOpeningAngleBinHigh; i++ ) { 
//      fInvarMasspTRotBkgAngleScaled = 
 // OpeningAngleBinLow

      //  fInvarMasspTRotBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");
        //fInvarMasspTRotBkg->SetName("hInvarMasspTRotBkg");
      //Pi0Cands->GetAxis(2)->SetRange(i,i);

      if (bPtOpAngleCorr) {

        printf("Using 2-Dim Opening Angle Correction.  Theta bin %d\n",i);
        // instead of using a single scale, project a opening angle bin of the correction onto pT
        // Except root doesn't have a 1D hist * 2D hist function :| 
        int nTotalPtBins = Pi0Cands->GetAxis(0)->GetNbins();

        Pi0Cands->GetAxis(2)->SetRange(i+OpeningAngleBinLow,i+OpeningAngleBinLow);

        TH2F * fLocalInvarMasspTRotBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fLocalInvarMasspTRotBkg->SetName("local");

        //for (int j = 0; j < nPtBins; j++) { // or use finer pt bins???
        for (int j = 0; j < nTotalPtBins; j++) {

          double scale = 1;
          if (bApplyOpeningAngleCorrection) {
            scale = h2DOpeningAngleCorr->GetBinContent(j+1,i+1);  // i+1 is theta bin, j+1 is pt bin 
          }
//          printf("     OpeningAngle Scaling %s by separate pT bins %d ,scale = %f\n",fLocalInvarMasspTRotBkg->GetName(),j,scale);

          int nMassBins = Pi0Cands->GetAxis(1)->GetNbins();
          for (int k = 0; k < nMassBins; k++ ) {
            fLocalInvarMasspTRotBkg->SetBinContent(k+1,j+1,scale*fLocalInvarMasspTRotBkg->GetBinContent(k+1,j+1));
            // FIXME what about error propagation??
            double localError = 1; //??
  //            fLocalInvarMasspTRotBkg->SetBinError(fLocalInvarMasspTRotBkg->FindBin(k+1,j+1),localError);
          }
        }
        printf("     OpeningAngle Scaling %s by separate pT bins\n",fLocalInvarMasspTRotBkg->GetName());
        printf(" i = %d\n",i);
        // Store the local histogram as the final scaled histogram if it doesn't exist yet
        // Otherwise, add it
        if (!fInvarMasspTRotBkgAngleScaled) {
          fInvarMasspTRotBkgAngleScaled = fLocalInvarMasspTRotBkg;
          fInvarMasspTRotBkgAngleScaled->SetName("hInvarMasspTRotBkgAngleScaled");
        } else {
          fInvarMasspTRotBkgAngleScaled->Add(fLocalInvarMasspTRotBkg);
          delete fLocalInvarMasspTRotBkg;
        }
      } else {

        Pi0Cands->GetAxis(2)->SetRange(i+OpeningAngleBinLow,i+OpeningAngleBinLow);
        TH2F * fLocalInvarMasspTRotBkg = (TH2F *) Pi0Cands->Projection(0,1,"e");
        fLocalInvarMasspTRotBkg->SetName("local");

        double scale = 1;
        if (bApplyOpeningAngleCorrection) {
          scale = hOpeningAngleCorrection->GetBinContent(i+1);
        }
        fLocalInvarMasspTRotBkg->Scale(scale);
        printf("     OpeningAngle Scaling %s by %f\n",fLocalInvarMasspTRotBkg->GetName(),scale);
        printf(" i = %d\n",i);
        if (!fInvarMasspTRotBkgAngleScaled) {
          fInvarMasspTRotBkgAngleScaled = fLocalInvarMasspTRotBkg;
          fInvarMasspTRotBkgAngleScaled->SetName("hInvarMasspTRotBkgAngleScaled");
        } else {
          fInvarMasspTRotBkgAngleScaled->Add(fLocalInvarMasspTRotBkg);
          delete fLocalInvarMasspTRotBkg;
        }
      }
    }
  } else {
    printf("Not applying opening angle correction.\n");
  }
  Pi0Cands->GetAxis(3)->SetRange(1,nThetaBins);


}

void PionID::SetPtBins() {
  cout<<"Setting up Pt Bins"<<endl;

  TCanvas * c2 = new TCanvas("c2");
  int initBin  = fInvarMasspT->GetXaxis()->FindBin(0.10);
  int finalBin = fInvarMasspT->GetXaxis()->FindBin(0.17);

  TH1D * hpT;
  if (fInvarMasspTNoCuts) hpT = fInvarMasspTNoCuts->ProjectionY("hpT",initBin,finalBin);
  else hpT = fInvarMasspT->ProjectionY("hpT",initBin,finalBin);

  TH1D * hpTInt = IntegralHist(hpT);


  switch (PtBinChoice) {
  case 0: 
    Pi0PtBins = FindEqualBins(hpTInt,nPtBins);
    break;
  case 1:
    Pi0PtBins = {5,7,9,11,14,17}; // final bins
    nPtBins = 5;
    break;
  default:
  case 2:
    Pi0PtBins = {4,6,8,10,12,14,16,18,20};
    break;
  case 3:
    Pi0PtBins = {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    nPtBins = 16;
    break;
  case 4:
    Pi0PtBins = {4,5,7,9,11,14,17,20,22,30}; // const request the good one
    nPtBins = 7;
    break;
  }

  //Pi0PtBins = FindEqualBins(hpTInt,nPtBins);
  //Pi0PtBins = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20, 21, 22, 23, 24, 25, 27, 29, 31, 35};
  //Pi0PtBins = {5,7,9,11,14,17,22,30,60};
  // Final + extra:
  //Pi0PtBins = {3,4,5,7,9,11,14,17,22,30,60};
  //Pi0PtBins = {3,4,5,7,9,11,14,17,22,30};
  //Pi0PtBins = {1,2,3,4,5}; // const request
  //Pi0PtBins = {3,4,5,7,9,11,14,17,20,22,30}; // const request the good one
//  Pi0PtBins = {3,4,5,7,9,11,14,17,20,22,30}; // const request the good one
  //Pi0PtBins = {5,7,9,11,14,17,20,22,30}; // 
  //Pi0PtBins = {1,2,3,4,5,7,9,11,14,17,20,22,30}; // low pt included
  //Pi0PtBins = {3,4,5,7,9,11,14,17,18,19,20,22,30}; // const request

  hpTInt->Draw("HIST");

  printf("Using bins: ");
  for (int i = 0; i <= nPtBins; i++) {
    printf("%.3f ",Pi0PtBins[i]);
  }
  printf("\n");


}

void PionID::MakeBasicPlots() {
  cout<<"Making some basic plots"<<endl;


  // Getting Cluster Energy Spectra
  if (fClusEnergy) {
    fClusEnergy->Sumw2();
    hGammaE = (TH1D *) fClusEnergy->Clone("hGammaE");
  } else if (ClusterProp) {
    hGammaE = (TH1D * ) ClusterProp->Projection(0);
  } else {
    hGammaE = (TH1D * ) fHistEvsPt->ProjectionY("hGammaE");
  }
//  hGammaE->Sumw2();
//  hGammaE->Rebin(2);
  //hGammaE->Scale(1./hGammaE->Integral("width"));

  TCanvas * cGammaE = new TCanvas("cGammaE","cGammaE",fCanvasWidth,fCanvasHeight);
  cGammaE->cd();
  hGammaE->Draw();
  cGammaE->SetLogy();
  cGammaE->Print(Form("%s/ClusterEnergy.pdf",sOutputDir.Data()));
  cGammaE->Print(Form("%s/ClusterEnergy.C",sOutputDir.Data()));


  // Drawing plots for the track cluster matching

  // From my task: 
  TF1 * fFuncMyTaskPtDepEta = new TF1("MyEtaCutFunc","0.010 + TMath::Power((x + 4.07), -2.5)",0,30);
  TF1 * fFuncMyTaskPtDepPhi = new TF1("MyPhiCutFunc","0.015 + TMath::Power((x + 3.65), -2.)" ,0,30);
  TF1 * fFuncMyTaskPtDepEta2 = new TF1("MyEtaCutFunc2","-(0.010 + TMath::Power((x + 4.07), -2.5))",0,30);
  TF1 * fFuncMyTaskPtDepPhi2 = new TF1("MyPhiCutFunc2","-(0.015 + TMath::Power((x + 3.65), -2.))" ,0,30);
//    if (etaCut <= 0) etaCut = 0.010 + TMath::Power((trackPt + 4.07), -2.5);
//    if (phiCut <= 0) phiCut = 0.015 + TMath::Power((trackPt + 3.65), -2.);

  int kCutLineColor = kRed+1;
  int kCutLineStyle = 2;
  fFuncMyTaskPtDepEta->SetLineColor(kCutLineColor);
  fFuncMyTaskPtDepEta2->SetLineColor(kCutLineColor);
  fFuncMyTaskPtDepPhi->SetLineColor(kCutLineColor);
  fFuncMyTaskPtDepPhi2->SetLineColor(kCutLineColor);

  fFuncMyTaskPtDepEta->SetLineStyle(kCutLineStyle);
  fFuncMyTaskPtDepEta2->SetLineStyle(kCutLineStyle);
  fFuncMyTaskPtDepPhi->SetLineStyle(kCutLineStyle);
  fFuncMyTaskPtDepPhi2->SetLineStyle(kCutLineStyle);

  // The function in the EMCal correction framework:
  //https://alice-notes.web.cern.ch/node/411  Fig. 46   for mom<3GeV these cuts are wider than the standard cuts
  TF1* fFuncPtDepEta = new TF1("funcEta", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])",0,30);
  fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
  TF1* fFuncPtDepPhi = new TF1("funcPhi", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])",0,30);
  fFuncPtDepPhi->SetParameters(0.09, 0.015, 2. );

  fFuncPtDepEta->SetLineColor(kRed+1);
  fFuncPtDepPhi->SetLineColor(kRed+1);
  fFuncPtDepEta->SetLineStyle(2);
  fFuncPtDepPhi->SetLineStyle(2);

  // just tests to see what happens with the function
  double fTestMomentum = 4;
  double  phiCutlo = -fFuncPtDepPhi->Eval(fTestMomentum);
  double  phiCuthi = +fFuncPtDepPhi->Eval(fTestMomentum);
  double  etaCut   = fFuncPtDepEta->Eval(fTestMomentum);
  printf("Function tests returned %f %f %f\n",phiCutlo,phiCuthi,etaCut);     

  TCanvas * cTrackClusterCuts = new TCanvas("cTrackClusterCuts","cTrackClusterCuts",fCanvasWidth,fCanvasHeight);
  cTrackClusterCuts->Divide(2,2);
  bool bDrawTrackClusterPlots = true;
  if (bDrawTrackClusterPlots && fMatchDeltaEtaTrackPt) {
    cTrackClusterCuts->cd(1);
    fMatchDeltaEtaTrackPt->Draw("COLZ");
    fFuncMyTaskPtDepEta->Draw("SAME");
    fFuncMyTaskPtDepEta2->Draw("SAME");
//    fFuncPtDepEta->Draw("SAME");
    gPad->SetLogz();
    cTrackClusterCuts->cd(2);
    fMatchDeltaPhiTrackPt->Draw("COLZ");
    fFuncMyTaskPtDepPhi->Draw("SAME");
    fFuncMyTaskPtDepPhi2->Draw("SAME");
//    fFuncPtDepPhi->Draw("SAME");
    gPad->SetLogz();
    cTrackClusterCuts->cd(3);
    if (fMatchCondDeltaEtaTrackPt) fMatchCondDeltaEtaTrackPt->Draw("COLZ"); 
    fFuncMyTaskPtDepEta->Draw("SAME");
    fFuncMyTaskPtDepEta2->Draw("SAME");
    gPad->SetLogz();
    cTrackClusterCuts->cd(4);
    if (fMatchCondDeltaPhiTrackPt) fMatchCondDeltaPhiTrackPt->Draw("COLZ"); 
    fFuncMyTaskPtDepPhi->Draw("SAME");
    fFuncMyTaskPtDepPhi2->Draw("SAME");
    gPad->SetLogz();

    cTrackClusterCuts->Print(Form("%s/TrackClusterCuts.pdf",sOutputDir.Data()));
    cTrackClusterCuts->Print(Form("%s/TrackClusterCuts.C",sOutputDir.Data()));
  }
//  delete HistoList;// FIXME trying debug

//  TCanvas * cPairTheta = new TCanvas("cPairTheta");
//  cPairTheta->cd();
//  hPairTheta->Draw();

}

void PionID::InitializeWSignal() {
  cout<<"Initalizing W Signal TF2."<<endl;

  TString Formuoli = "";
  Bool_t bUseParametrizedModel = true; // FIXME

  // Ravioli, Ravioli, give me the ...
  if (bUseParametrizedModel) {
//    Formuoli = "[0]*TMath::Gaus(x,[3],[4])*(y > [1])*TMath::Exp(-[2]*y)";

    TString sParametrizedMass = "(y < [3]) * ([5]*y+[4]-[5]*[3]) + (y >= [3]) * ([6]*y+[4]-[6]*[3])";
    TString sParametrizedSigma = "(y < [7]) * ([9]*y+[8]-[9]*[7]) + (y >= [7]) * ([10]*y+[8]-[10]*[3])";

    Formuoli = TString::Format("[0]*TMath::Gaus(x,%s,%s)*(y > [1])*TMath::Exp(-[2]*y)",sParametrizedMass.Data(),sParametrizedSigma.Data());
//    Formuoli = "[0]*TMath::Gaus(x,(y < [3]) * ([5]*y+[4]-[5]*[3]) + (y >= [3]) * ([6]*y+[4]-[6]*[3]),(y < [7]) * ([9]*y+[8]-[9]*[7]) + (y >= [7]) * ([10]*y+[8]-[10]*[3]))*(y > [1])*TMath::Exp(-[2]*y)";

    // Don't need to have these as ROOT-recognized parameters

    // new parameters d,e,m1,m2 for mass, sigma
    // [3] = d_mass
    // [4] = e_mass
    // [5] = m1_mass
    // [6] = m2_mass
// mass = (y < [3]) * ([5]*y+[4]-[5]*[3]) + (y >= [3]) * ([6]*y+[4]-[6]*[3])
    // [7] = d_sigma
    // [8] = e_sigma
    // [9] = m1_sigma
    // [10] = m2_sigma
// sigma = (y < [7]) * ([9]*y+[8]-[9]*[7]) + (y >= [7]) * ([10]*y+[8]-[10]*[3])

    // The original mass function
  //TString formuoli ="(x < [0]) * ([2]*x+[1]-[2]*[0]) + (x >= [0]) * ([3]*x+[1]-[3]*[0])";
    // The original sigma function
  //TString formuoli ="(x < [0]) * ([2]*x+[1]-[2]*[0]) + (x >= [0]) * ([3]*x+[1]-[3]*[0])";

  } else {
    Formuoli = "[0]*TMath::Gaus(x,[3],[4])*(y > [1])*TMath::Exp(-[2]*y)";
  }
  
  // With Parametrized mass and sigma

  // scaling pt with pt0
//  TString Formuoli = "[0]*TMath::Gaus(x,[3],[4])*(y > [1])*TMath::Exp(-[2]*y/[1])";

  // Simple Gauss
 // TString Formuoli = "[0]*TMath::Gaus(x,[3],[4])*(y > [1])*TMath::Exp(-[2]*y/[1])";
  
  // Power Law Model
//  TString Formuoli = "[0]*TMath::Gaus(x,[3],[4])*(y > [1])*TMath::Power(y/[1],-[2])";

  W2DSignalModel = new TF2("WSignalTF2",Formuoli,0.,1.,0.,30.);

  // Pt Dependence Parameters
  W2DSignalModel->SetParName(0,"Yield");
  W2DSignalModel->SetParName(1,"pT0");
  W2DSignalModel->SetParName(2,"Power");
  W2DSignalModel->SetParameter(0,fWYield);
  W2DSignalModel->FixParameter(1,fWPt0);
  W2DSignalModel->SetParameter(2,fWPower);

  // Mass Dependence Parameters
  if (bUseParametrizedModel) {
    W2DSignalModel->SetParName(3,"d_mass");
    W2DSignalModel->SetParName(4,"e_mass");
    W2DSignalModel->SetParName(5,"m1_mass");
    W2DSignalModel->SetParName(6,"m2_mass");

    W2DSignalModel->SetParName(7,"d_sigma");
    W2DSignalModel->SetParName(8,"e_sigma");
    W2DSignalModel->SetParName(9,"m1_sigma");
    W2DSignalModel->SetParName(10,"m2_sigma");

    W2DSignalModel->FixParameter(3,fWD_Mass);
    W2DSignalModel->FixParameter(4,fWE_Mass);
    W2DSignalModel->FixParameter(5,fWM1_Mass);
    W2DSignalModel->FixParameter(6,fWM2_Mass);

    W2DSignalModel->FixParameter(7,fWD_Sigma);
    W2DSignalModel->FixParameter(8,fWE_Sigma);
    W2DSignalModel->FixParameter(9,fWM1_Sigma);
    W2DSignalModel->FixParameter(10,fWM2_Sigma);
  
  } else {
    W2DSignalModel->SetParName(3,"Mass");
    W2DSignalModel->SetParName(4,"Sigma");

    W2DSignalModel->SetParameter(3,fWMass);
    W2DSignalModel->SetParameter(4,fWSigma);
  }

  W2DSignalHist = new TH2F("W2DSignalHist","W Signal histogram;m_{#gamma#gamma} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})",100,0,1,100,0,30);
  // FIXME fill with 
  // W2DSignalModel->Integral(fMassEdgeLow,fMassEdgeHigh,fPtEdgeLow,fPtEdgeHigh) / (fMassBinSize * fPtBinSize);


  TCanvas * cW2D = new TCanvas("cW2D","cW2D");
  //W2DSignalHist->Draw("COLZ");
  //W2DSignalModel->Draw("SAME");
  W2DSignalModel->SetNpx(500);
  W2DSignalModel->SetNpy(500);
  W2DSignalModel->Draw("");

  cW2D->Print(Form("%s/PSWSignal.pdf",sOutputDir.Data()));
  cW2D->Print(Form("%s/PSWSignal.C",sOutputDir.Data()));

}

void PionID::DrawPSScaleCorrectionPlots() {
  cout<<"Position scaling Swap Correction section."<<endl;

  TCanvas * cPSC = new TCanvas("cPSC","cPSC",fCanvasWidth,fCanvasHeight);

/* I delete this histograms
  // 1 Dim Versions

  // Basic histograms to plot things
  TH1F * fUDisplay = (TH1F *) fUDist->Projection(0);
  fUDisplay->SetName("UDisplay_Raw");
  TH1F * fUTildeDisplay = (TH1F *) fUTildeDist->Projection(0);
  fUTildeDisplay->SetName("UTildeDisplay_Raw");
  TH1F * fVDisplay = (TH1F *) fVDist->Projection(0);
  fVDisplay->SetName("VDisplay_Raw");
  TH1F * fVTildeDisplay = (TH1F *) fVTildeDist->Projection(0);
  fVTildeDisplay->SetName("VTileDisplay_Raw");

  // Now apply cuts
  // Could shorten
  fUDist->GetAxis(1)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fUDist->GetAxis(2)->SetRange(EnergyBinLow,nEnergyBinHigh);
  fUTildeDist->GetAxis(1)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fUTildeDist->GetAxis(2)->SetRange(EnergyBinLow,nEnergyBinHigh);
  fVDist->GetAxis(1)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fVDist->GetAxis(2)->SetRange(EnergyBinLow,nEnergyBinHigh);
  fVTildeDist->GetAxis(1)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fVTildeDist->GetAxis(2)->SetRange(EnergyBinLow,nEnergyBinHigh);

  TH1F * fUCut = (TH1F *) fUDist->Projection(0);
  TH1F * fUTildeCut = (TH1F *) fUTildeDist->Projection(0);
  TH1F * fVCut = (TH1F *) fVDist->Projection(0);
  TH1F * fVTildeCut = (TH1F *) fVTildeDist->Projection(0);

  if (fUDisplay->Integral() > 0) fUDisplay->Scale(1.0/fUDisplay->Integral()); 
  if (fUTildeDisplay->Integral() > 0) fUTildeDisplay->Scale(1.0/fUTildeDisplay->Integral()); 
  if (fVDisplay->Integral() > 0) fVDisplay->Scale(1.0/fVDisplay->Integral()); 
  if (fVTildeDisplay->Integral() > 0) fVTildeDisplay->Scale(1.0/fVTildeDisplay->Integral()); 

  if (fUCut->Integral() > 0) fUCut->Scale(1.0/fUCut->Integral()); 
  if (fUTildeCut->Integral() > 0) fUTildeCut->Scale(1.0/fUTildeCut->Integral()); 
  if (fVCut->Integral() > 0) fVCut->Scale(1.0/fVCut->Integral()); 
  if (fVTildeCut->Integral() > 0) fVTildeCut->Scale(1.0/fVTildeCut->Integral()); 

  TLegend * legDisplay = new TLegend(0.5,0.4,0.9,0.9);
  TLegend * legCut = new TLegend(0.5,0.4,0.9,0.9);
  
  string fPSCName[4] = {"Same E Mass Scaling","Same E #it{p}_{T} Scaling","Same #vec{x} Mass Scaling","Same #vec{x} #it{p}_{T} Scaling"};
  TH1F * fPSCDisplayHist[4] = {fUDisplay,fUTildeDisplay,fVDisplay,fVTildeDisplay};
  TH1F * fPSCCutHist[4] = {fUCut,fUTildeCut,fVCut,fVTildeCut};

  Int_t kPSCColor[4] = {kRed+1,kOrange+8,kBlue+1,kTeal-1};
  Int_t kPSCStyle[4] = {21,21,21,21};
  Int_t kPSCCutStyle[4] = {25,25,25,25};

  for (int i = 0; i < 4; i++) {
    SetTH1Histo(fPSCDisplayHist[i],"Scale","Arb. Units",true);
    fPSCDisplayHist[i]->SetLineColor(kPSCColor[i]);
    fPSCDisplayHist[i]->SetMarkerColor(kPSCColor[i]);
    fPSCDisplayHist[i]->SetMarkerStyle(kPSCStyle[i]);
    legDisplay->AddEntry(fPSCDisplayHist[i],fPSCName[i].c_str(),"lp");

    SetTH1Histo(fPSCCutHist[i],"Scale (After Cuts)","Arb. Units",true);
    fPSCCutHist[i]->SetLineColor(kPSCColor[i]);
    fPSCCutHist[i]->SetMarkerColor(kPSCColor[i]);
    fPSCCutHist[i]->SetMarkerStyle(kPSCCutStyle[i]);
    legCut->AddEntry(fPSCCutHist[i],fPSCName[i].c_str(),"lp");
  }

//  fUDisplay->SetLineColor(kPSCColor[0]);
 // fUTildeDisplay->SetLineColor(kPSCStyle[1]);
//  fVDisplay->SetLineColor(kBlue+1);
//  fVTildeDisplay->SetLineColor(kTeal-1);

  fUDisplay->Draw();
  fUTildeDisplay->Draw("SAME");
  fVDisplay->Draw("SAME");
  fVTildeDisplay->Draw("SAME");

  legDisplay->Draw("SAME");  

  cPSC->Print(Form("%s/PosSwapCorr_Inc.pdf",sOutputDir.Data()));
  cPSC->Print(Form("%s/PosSwapCorr_Inc.C",sOutputDir.Data()));

  cPSC->Clear();
  fUCut->Draw();
  fUTildeCut->Draw("SAME");
  fVCut->Draw("SAME");
  fVTildeCut->Draw("SAME");

  legCut->Draw("SAME");  

  cPSC->Print(Form("%s/PosSwapCorr_Cut.pdf",sOutputDir.Data()));
  cPSC->Print(Form("%s/PosSwapCorr_Cut.C",sOutputDir.Data()));
  */

  // Draw 2D version
  fUScaleMatrix->GetAxis(2)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fUScaleMatrix->GetAxis(3)->SetRange(EnergyBinLow,nEnergyBinHigh);
  fVScaleMatrix->GetAxis(2)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fVScaleMatrix->GetAxis(3)->SetRange(EnergyBinLow,nEnergyBinHigh);

  fU2DCut = (TH2F *) fUScaleMatrix->Projection(1,0);
  fV2DCut = (TH2F *) fVScaleMatrix->Projection(1,0);

  fU2DCut->Scale(1.0/fU2DCut->Integral("width"));
  fV2DCut->Scale(1.0/fV2DCut->Integral("width"));

  fU2DCut->Draw("COLZ");

  cPSC->SetLogz();

  cPSC->Print(Form("%s/PosSwapUScaleMatrix_Cut.pdf",sOutputDir.Data()));
  cPSC->Print(Form("%s/PosSwapUScaleMatrix_Cut.C",sOutputDir.Data()));

  fV2DCut->Draw("COLZ");

  cPSC->Print(Form("%s/PosSwapVScaleMatrix_Cut.pdf",sOutputDir.Data()));
  cPSC->Print(Form("%s/PosSwapVScaleMatrix_Cut.C",sOutputDir.Data()));

  BuildPSScaleDistributions();

}

void PionID::BuildPSScaleDistributions() {
  cout<<"Building Energy-pair and Position-pair PosSwapped distributions"<<endl;

  if (!W2DSignalModel) {
    fprintf(stderr,"Error: Missing TF2 of signal.\n");
    return;
  }

  // Building the E Distribution
  fPSEnergyPairDist = (TH2F *) fInvarMasspT->Clone("PSEnergyPairDist");
  fPSEnergyPairDist->Reset();

//  fPSEnergyPairDist->Rebin2D(4,4); //FIXME

  Int_t nBinsM = fPSEnergyPairDist->GetXaxis()->GetNbins();
  Int_t nBinsP = fPSEnergyPairDist->GetYaxis()->GetNbins();

  printf("Integrating ... \n");
  // Loop over all mass and pt bins
  for (Int_t i = 1; i <= nBinsM; i++) {
    printf("%03d/%03d\n",i,nBinsM);
    // Local Mass Value
    Double_t lMass = fPSEnergyPairDist->GetXaxis()->GetBinCenter(i);
    for (Int_t j = 1; j <= nBinsP; j++) {
      Double_t Value = 1;
      // Local Pt Value
      // could do modification here for falling pt spectrum
      Double_t lPt   = fPSEnergyPairDist->GetYaxis()->GetBinCenter(j);
      Value = IntegrateScaleModification(fU2DCut,lMass,lPt); //error?
      fPSEnergyPairDist->SetBinContent(i,j,Value);
      // Set Error?
     // printf("%d %d Value = %f \n",i,j,Value);
    }  
  }

  //outFile->Add(fPSEnergyPairDist);
  //outFile->Add(fPSPosPairDist);

}
/**
  * Convolves the scaling modification distribution against the input signal distribution
  */
double PionID::IntegrateScaleModification(TH2 * fDist, Double_t lMass, Double_t lPt) {
  if (!fDist) return -1;
  Double_t Value = 0;

  Int_t nMassIntBins = 150;
  Double_t fMinIntMass = 0.040;
  Double_t fMaxIntMass = 0.70;
  Double_t fMassIntStep = (fMaxIntMass - fMinIntMass) / nMassIntBins;

  Int_t nPtIntBins = 45;
  Double_t fMinIntPt = 5.0;
  Double_t fMaxIntPt = 20.0;
  Double_t fPtIntStep = (fMaxIntPt - fMinIntPt) / nPtIntBins;



  // FIXME this also makes no sense 
  // set my own mass and pt bins, range
  Int_t nBinsM = fDist->GetXaxis()->GetNbins();
  Int_t nBinsP = fDist->GetYaxis()->GetNbins();

  Double_t fDistMaxX = fDist->GetXaxis()->GetXmax();
  Double_t fDistMinX = fDist->GetXaxis()->GetXmin();
  Double_t fDistMaxY = fDist->GetYaxis()->GetXmax();
  Double_t fDistMinY = fDist->GetYaxis()->GetXmin();
/*
  if (bPSCorrLogMod) {
    fDistMaxX = TMath::Exp(fDistMaxX);
    fDistMinX = TMath::Exp(fDistMinX);
    fDistMaxY = TMath::Exp(fDistMaxY);
    fDistMinY = TMath::Exp(fDistMinY);
  } 
*/

  // Loop over all mass and pt bins
  //for (Int_t i = 1; i <= nBinsM; i++) {
  for (Int_t i = 0; i < nMassIntBins; i++) {

    // Local Mass Value 
    // FIXME this makes no sense
//    Double_t mMass = fDist->GetXaxis()->GetBinCenter(i);
//    Double_t mMassWidth = fDist->GetXaxis()->GetBinWidth(i);
    Double_t mMass = fMinIntMass + i * fMassIntStep;
    Double_t mMassWidth = fMassIntStep;

    // Mass for Dist histogram (ratio of Source Mass to current mass)
    Double_t RatioMass = lMass / mMass;
    if (bPSCorrLogMod) {
      RatioMass = TMath::Log(lMass / mMass);
    }
//    for (Int_t j = 1; j <= nBinsP; j++) {
    for (Int_t j = 0; j < nPtIntBins; j++) {
      // Local Pt Value
      //Double_t mPt = fDist->GetYaxis()->GetBinCenter(j);
      //Double_t mPtWidth = fDist->GetYaxis()->GetBinWidth(j);
      Double_t mPt = fMinIntPt + j * fPtIntStep;
      Double_t mPtWidth = fPtIntStep;

      // Pt for Dist histogram (Source)
      Double_t RatioPt = lPt / mPt;
      if (bPSCorrLogMod) {
        RatioPt = TMath::Log(RatioPt);
      }
      // FIXME add in log mode
      // Check if sMass in range
    
      //printf("%d %d L(RMass) = %f L(RPt) = %f\n",i,j,RatioMass,RatioPt);

      if (RatioMass >= fDistMaxX || RatioMass < fDistMinX || RatioPt >= fDistMaxY || RatioPt < fDistMinY) continue;
      if (bPSCorrLogMod) {
        Value += mMassWidth * mPtWidth * W2DSignalModel->Eval(mMass,mPt) * fDist->Interpolate(RatioMass,RatioPt) / (mPt * mMass); 
      } else { 
        Value += mMassWidth * mPtWidth * W2DSignalModel->Eval(mPt,mMass) * fDist->Interpolate(RatioMass,RatioPt) / (mPt * mMass); 
      }
//      Value += mMassWidth * mPtWidth * W2DSignalModel->Eval(RatioMass,RatioPt) * fDist->GetBinContent(fDist->FindFixBin(RatioMass,RatioPt)) / (mPt * mMass); 
      // Division comes from dirac delta normalization (or change of variables)
    }
  }
  return Value;
}

void PionID::DrawPSDirectCorrectionPlots() {
  cout<<"Position Swap Direct Map Correction section"<<endl;

  TCanvas * cPSD = new TCanvas("cPSD","cPSD",2*fCanvasWidth,fCanvasHeight);

  fPSMassPtMap->GetAxis(4)->SetRange(LambdaBinLow,nLambdaBinHigh);
  fPSMassPtMap->GetAxis(5)->SetRange(EnergyBinLow,nEnergyBinHigh);

  fPSInitialMap = (TH2F *) fPSMassPtMap->Projection(1,0);
  fPSInitialMap->SetName("PSInitialMap");
  fPSFinalMap   = (TH2F *) fPSMassPtMap->Projection(3,2);
  fPSFinalMap->SetName("PSFinalMap");
 
  cPSD->Divide(2,1);
  cPSD->cd(1);
  fPSInitialMap->Draw("COLZ");
  cPSD->cd(2);
  fPSFinalMap->Draw("COLZ");

  cPSD->Print(Form("%s/PosSwapDirectMapInitial.pdf",sOutputDir.Data()));
  cPSD->Print(Form("%s/PosSwapDirectMapInitial.C",sOutputDir.Data()));

  BuildPSDirectDistributions();
}

void PionID::BuildPSDirectDistributions() {
  cout<<"Building Direct Mass-Pt map PosSwapped Distributions"<<endl;

  if (!W2DSignalModel) {
    fprintf(stderr,"Error: Missing TF2 of signal.\n");
    return;
  }

  TCanvas * cDebugPS = new TCanvas("cDebugPS","cDebugPS",fCanvasWidth,fCanvasHeight);
  cDebugPS->cd();

  // Histogram for the final result
  //fPSFinalPeak = (TH2F *) fInvarMasspT->Clone("PSFinalPeak");
  fPSFinalPeak = (TH2F *) fPSMassPtMap->Projection(3,2);
  fPSFinalPeak->SetName("PSFinalPeak");
  fPSFinalPeak->Reset();
  // FIXME make this less hard-coded
  //fPSFinalPeak->Rebin2D(2,1); // just to get from 

  //Int_t nRebinMass_Int = (Int_t) fPSFinalPeak->GetXaxis()->GetNbins()/(fInvarMasspT->GetXaxis()->GetNbins()/nRebinMass);
  Int_t nRebinMass_Int = 2;

  fPSFinalPeak->Rebin2D(nRebinMass_Int,1);

  Int_t nBinsM = fPSFinalPeak->GetXaxis()->GetNbins();
  Int_t nBinsP = fPSFinalPeak->GetYaxis()->GetNbins();

  // Code will be faster if we can make assumptions about the 
  // binning of fInvarMasspT and fPSMassPtMap

 // printf("Debug: fInvarMasspT has NBinsX = %03d NBinsY = %03d\n",fInvarMasspT->GetXaxis()->GetNbins(),fInvarMasspT->GetYaxis()->GetNbins());
  printf("       fPSFinalPeak has NBinsX = %03d NBinsY = %03d\n",nBinsM,nBinsP);
  printf("       fPSMassPtMap has NBinsX = %03d NBinsY = %03d\n",fPSMassPtMap->GetAxis(2)->GetNbins(),fPSMassPtMap->GetAxis(3)->GetNbins());

//  if (fInvarMasspT->GetXaxis()->GetNbins() != nBinsM || fInvarMasspT->GetYaxis()->GetNbins() != nBinsP) {
//    fprintf(stderr,"fInvarMassPt and fPSFinalPeak don't have the same binning. Check rebin settings.\n");
//    return;
//  }


  // Integral Ranges (make these setable)

  Double_t fMinIntMass = 0.040;
  Double_t fMaxIntMass = 0.20;
  //Double_t fMinIntPt = 5.0;
  Double_t fMinIntPt = fWPt0;
  Double_t fMaxIntPt = 20.0;

  Int_t fMinIntMassBin = fPSFinalPeak->GetXaxis()->FindFixBin(fMinIntMass);
  Int_t fMaxIntMassBin = fPSFinalPeak->GetXaxis()->FindFixBin(fMaxIntMass);
  Int_t fMinIntPtBin = fPSFinalPeak->GetYaxis()->FindFixBin(fMinIntPt);
  Int_t fMaxIntPtBin = fPSFinalPeak->GetYaxis()->FindFixBin(fMaxIntPt);

  printf("Integrating over direct map.\n");
  printf("  Integral Input Range: Mass %d - %d  (%2.2f - %2.2f GeV/#it{c}^2)\n",fMinIntMassBin,fMaxIntMassBin,fMinIntMass,fMaxIntMass);
  printf("                        pT   %d - %d  (%2.2f - %2.2f GeV/#it{c})\n",fMinIntPtBin,fMaxIntPtBin,fMinIntPt,fMaxIntPt);
  Double_t Value = 0;
//  for (Int_t i = 1; i <= nBinsM; i++) {
  for (Int_t i = fMinIntMassBin; i <= fMaxIntMassBin; i++) {
    printf("%03d/%03d\n",i,nBinsM);
    fPSMassPtMap->GetAxis(0)->SetRange(i,i);
    Double_t lMass = fPSFinalPeak->GetXaxis()->GetBinCenter(i);
    Double_t fMassBinSize = fPSFinalPeak->GetXaxis()->GetBinWidth(i);
//    for (Int_t j = 1; j <= nBinsP; j++) {
    for (Int_t j = fMinIntPtBin; j <= fMaxIntPtBin; j++) {
      // For each bin, we make a new final histogram
      // normalized to unity. Then we convolve with the 
      // signal model
      TH2F * fLocalHist = 0;
      Double_t lPt = fPSFinalPeak->GetYaxis()->GetBinCenter(j);
      Double_t fPtBinSize = fPSFinalPeak->GetYaxis()->GetBinWidth(i);
      fPSMassPtMap->GetAxis(1)->SetRange(j,j);
      fLocalHist = (TH2F *) fPSMassPtMap->Projection(3,2);

      fLocalHist->Rebin2D(nRebinMass_Int,1);

      // FIXME testing what happens with a delta-ish distribution
      //fLocalHist = (TH2F *) fPSMassPtMap->Projection(1,0);
      //fLocalHist->Reset(); //FIXME 
      //fLocalHist->SetBinContent(i,j,1.); //FIXME

      Double_t fIntegral = fLocalHist->Integral(); // Should this be "width"?
      Double_t fLocalWDistValue = 0;
      if (fIntegral != 0) { 

        // High performance value: FIXME
        //  should build this as a histogram first, access W2DSignalHist
        fLocalWDistValue = W2DSignalModel->Integral(lMass-0.5*fMassBinSize,lMass+0.5*fMassBinSize,lPt-0.5*fPtBinSize,lPt+0.5*fPtBinSize) / (fMassBinSize * fPtBinSize);
        // Low performance value:
        //fLocalWDistValue = W2DSignalModel->Eval(lMass,lPt);

        fLocalHist->Scale(fLocalWDistValue/fIntegral);
        // Debug
 /*       fLocalHist->SetTitle(Form("M_{init} = %2.2f #it{p}_{T,init} = %2.1f",lMass,lPt));
        fLocalHist->Draw();
        cDebugPS->Print(Form("debug/PSBuild_%d_%d.png",i,j));*/
      }
      fPSFinalPeak->Add(fLocalHist);     

//      Value = IntegrateDirectModification(,lMass,lPt);

      delete fLocalHist;
    }
  }
}

void PionID::GetThetaModelParameters(double fPt, double fThetaC, double &ThetaModelLambda, double &ThetaModelMPrime) {
  printf("GetThetaModel Called with Pt = %f, ThetaC = %f\n",fPt,fThetaC);
 
  // FIXME and also include an option with look up table
 
  // bUseThetaLookUpTable
  if (bUseThetaLookUpTable) {
  
    Double_t fLambdaValue = 0;
    Double_t fMassPrimeValue = 0;

    // Need to pick the appropriate pT and Theta bin


    // FIXME
    Int_t iThetaPtBin = 0; 
    Int_t iThetaThetaBin = 0;

    // Find Pt Bin in Pi0PtBins 
    //  get the pt bins from the MassPrimePar1Graph->GetXaxis().
    // then we don't have to worry too much about mismatching axes between here
    // and the angle analysis, though we should still worry
    for (int i = 0; i < iNPtBinsThetaModel; i++) {
      // This requires the fPt and bins be exact
      if (fPt == fMassPrimePar1Graphs[iThetaModelCent]->GetX()[i]) {
        iThetaPtBin = i;
        break;
      }
    }  

    if (iThetaPtBin == -1) {
      printf("Underflow in pt bin for theta\n");
      return;
    }

    TGraphErrors * fThetaLambdaGraph    = fThetaLambdaGraphs[iThetaModelCent][iThetaPtBin];
    TGraphErrors * fThetaMassPrimeGraph = fThetaMassPrimeGraphs[iThetaModelCent][iThetaPtBin];

    printf("  Theta Model for pT = %f, is using bin %d\n",fPt,iThetaPtBin);

    for (int i = 0; i < fThetaMassPrimeGraph->GetN(); i++) {
      if (fThetaC == fThetaMassPrimeGraph->GetX()[i]) {
        iThetaThetaBin = i;
        break;
      }
    }

    if (iThetaThetaBin == -1) {
      printf("Underflow in pt bin for theta\n");
      return;
    }
    printf("  Theta Model for ThetaC = %f, is using bin %d\n",fThetaC,iThetaThetaBin);
 
    fLambdaValue = fThetaLambdaGraph->GetY()[iThetaThetaBin];
    fMassPrimeValue = fThetaMassPrimeGraph->GetY()[iThetaThetaBin];
   
    ThetaModelLambda = fLambdaValue;
    ThetaModelMPrime = fMassPrimeValue;
    printf("theta model look up found lambda = %f and m' = %f\n",fLambdaValue,fMassPrimeValue);
    return;
  }

  // Pt Dependence parameters

  // Default (iThetaModelParamChoice == 0)
  // From T39 MB E_Min = 2 Same Event (0 < A < 1.0) 
  double pLambda_10 = -1.852351;
  double pLambda_11 = 9.795795;
  double pLambda_12 = -4.243938; 
  double pMPrime_10 = 2.; // Cluster Energy Min
  double pMPrime_11 = 1.213979; 
  double pMPrime_12 = 0.012710;

  if (iThetaModelParamChoice == 1) { // Latest and Greatest
    // From T40 HadCorr Cent0 (0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)
    if (iThetaModelTrigger == 0) { // MB
      switch (iThetaModelCent) {
      default:
      case 0:
				pLambda_10 = -1.343230; // \pm 0.086210
				pLambda_11 = 7.458780; // \pm 0.450105
				pLambda_12 = -2.594091; // \pm 0.465460
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 1.355498; // \pm 0.004236
				pMPrime_12 = 0.014616; // \pm 0.000100
        break;
      case 1:
        pLambda_10 = -1.136534; // \pm 0.064784
        pLambda_11 = 6.462728; // \pm 0.323681
        pLambda_12 = -1.679556; // \pm 0.353448
        pMPrime_10 = 2.000000; // \pm 0.000000
        pMPrime_11 = 1.225305; // \pm 0.005354
        pMPrime_12 = 0.014693; // \pm 0.000088
        break;
      case 2:
        pLambda_10 = -1.734734; // \pm 0.139330
        pLambda_11 = 9.624293; // \pm 0.796926
        pLambda_12 = -4.941174; // \pm 0.781095
        pMPrime_10 = 2.000000; // \pm 0.000000
        pMPrime_11 = 1.074171; // \pm 0.008286
        pMPrime_12 = 0.013348; // \pm 0.000126
        break;
      case 3:
        pLambda_10 = -2.128955; // \pm 0.288084
        pLambda_11 = 11.397052; // \pm 1.704273
        pLambda_12 = -5.990209; // \pm 1.517861
        pMPrime_10 = 2.000000; // \pm 0.000000
        pMPrime_11 = 0.971533; // \pm 0.014413
        pMPrime_12 = 0.011857; // \pm 0.000201
      }
    } else if (iThetaModelTrigger == 1) { // GA
      switch (iThetaModelCent) {
      default:
      case 0:
				pLambda_10 = -0.655594; // \pm 0.026680
				pLambda_11 = 4.388322; // \pm 0.129806
				pLambda_12 = -0.252177; // \pm 0.196335
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 1.058043; // \pm 0.005560
				pMPrime_12 = 0.014414; // \pm 0.000052
        break;
      case 1:
				pLambda_10 = -0.776916; // \pm 0.026300
				pLambda_11 = 5.047596; // \pm 0.134992
				pLambda_12 = -1.330238; // \pm 0.201064
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 0.904939; // \pm 0.004933
				pMPrime_12 = 0.014353; // \pm 0.000044
        break;
      case 2:
				pLambda_10 = -0.667893; // \pm 0.032602
				pLambda_11 = 4.410818; // \pm 0.158334
				pLambda_12 = -0.098910; // \pm 0.241189
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 0.643827; // \pm 0.008649
				pMPrime_12 = 0.013375; // \pm 0.000068
        break;
      case 3:
				pLambda_10 = -0.720727; // \pm 0.065700
				pLambda_11 = 4.670870; // \pm 0.329648
				pLambda_12 = -0.572612; // \pm 0.508737
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 0.533486; // \pm 0.015193
				pMPrime_12 = 0.012742; // \pm 0.000111
      }
    } else if (iThetaModelTrigger == 2) { // MC (not technically a trigger)
      // T38 (MC, Same Event, 0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)
      switch (iThetaModelCent) {
      default:
      case 0:
				pLambda_10 = -0.654782; // \pm 0.013737
				pLambda_11 = 4.433793; // \pm 0.064561
				pLambda_12 = 0.158259; // \pm 0.086945
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 1.087587; // \pm 0.002430
				pMPrime_12 = 0.014296; // \pm 0.000034

        break;
      case 1:
				pLambda_10 = -0.694652; // \pm 0.023331
				pLambda_11 = 4.697209; // \pm 0.114343
				pLambda_12 = -0.545263; // \pm 0.161476
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 1.082325; // \pm 0.003979
				pMPrime_12 = 0.014540; // \pm 0.000051

        break;
      case 2:
				pLambda_10 = -0.850918; // \pm 0.028144
				pLambda_11 = 5.538599; // \pm 0.148603
				pLambda_12 = -1.811707; // \pm 0.209824
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 1.037800; // \pm 0.004181
				pMPrime_12 = 0.014653; // \pm 0.000052

        break;
      case 3:
				pLambda_10 = -0.867271; // \pm 0.023312
				pLambda_11 = 5.564853; // \pm 0.121148
				pLambda_12 = -1.474318; // \pm 0.164174
				pMPrime_10 = 2.000000; // \pm 0.000000
				pMPrime_11 = 0.956292; // \pm 0.003741
				pMPrime_12 = 0.014507; // \pm 0.000042
      }
    }
  } 


  // From T40 MB HadCorr Cent0 (0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)
  if (iThetaModelParamChoice == 4) {
    pLambda_10 = -1.343230;
    pLambda_11 =  7.458780;
    pLambda_12 = -2.594091; 
    pMPrime_10 =  2.;
    pMPrime_11 =  1.355498; 
    pMPrime_12 =  0.014616;
  }

  // From T40 MB HadCorr Cent2 (0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)

  // From T40 GA HadCorr Cent0 (0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)

  // From T40 GA HadCorr Cent2 (0.1 < lambda < 0.7, 0 < alpha < 0.7, E_0 = 2)


  // From Same Event T33 MB E_min = 2.5
/*  if (iThetaModelParamChoice == 3) {
    pLambda_10 = -2.934598;
    pLambda_11 = 18.559780;
    pLambda_12 = -13.359569; 
    pMPrime_10 = 2.5; // Cluster Energy Min
    pMPrime_11 = 0.830040; 
    pMPrime_12 = 0.008016;
  }*/

  // From Same Event T33 MB E_min = 2
  // Older Parametrization
/*  double pLambda_10 = -2.574;
  double pLambda_11 = 4.314; 
  double pLambda_12 = 0.02344; 
  double pMPrime_10 = 9.072;
  double pMPrime_11 = 20;  */

  // From PosSwap T33 MB E_min = 2
/*  double pLambda_10 = -0.666;
  double pLambda_11 = 2.691; 
  double pLambda_12 = 0.0537; 
  double pMPrime_10 = 6.65;
  double pMPrime_11 = 14.12; 
*/
  // E cut = 1.5
/*  double pLambda_10 = 0.6293;
  double pLambda_11 = 4.389; 
  double pLambda_12 = 0.2827; 
  double pMPrime_10 = 4.824;
  double pMPrime_11 = 10.9; */

  // E cut = 2.0
  //double pMPrime_11 = 0.5; // Theoretical value

  // Theta Dependent Parameters
  double pLambda_0 = 0;
  double pLambda_1 = pLambda_10 + pLambda_11 / TMath::Sqrt(fPt - pLambda_12);
  //double pLambda_1 = pLambda_10 + pLambda_11 * TMath::Exp(-pLambda_12 * fPt);
//  double pLambda_1 = 0.7;
  double pMPrime_0 = 0;
  double pMPrime_1 = (1 + pMPrime_12 * fPt) *  TMath::Sqrt(pMPrime_10*(fPt - pMPrime_11 - pMPrime_10));
//  double pMPrime_1 = pMPrime_10 * TMath::ATan(fPt/pMPrime_11);
 // double pMPrime_1 = pMPrime_10 + pMPrime_11 * fPt;

  ThetaModelLambda = (pLambda_0 + pLambda_1 / fThetaC) / 1000.;
  ThetaModelMPrime = (pMPrime_0 + pMPrime_1 * fThetaC) * 1000.;
}


void PionID::ProduceDeltaPsiPlots() {
  /*
    TH2F * hPtEPAnglePionAcc=0;       // Pi0 Candidate
    TH2F * hPtEPAngleMCPion=0;        // MC Pi0
    TH2F * hPtEPAngleTrueRecMCPion=0; // MC Pi0 reconstructed as Pi0 Candidate
    TH3F * hHistTrackPsiEPPtCent=0;   // Tracks
    // Reaction Plane (Angle 0 in Data)
    TH2F * hPtRPAnglePionAcc=0;
    TH2F * hPtRPAngleMCPion=0;
    TH2F * hPtRPAngleTrueRecMCPion=0;
    TH3F * hHistTrackPsiRPPtCent=0;
*/
  //double NEvents = fHistCentrality->GetEntries();
  double NEvents = fHistEventHash->GetEntries();
  if (NEvents == 0) NEvents = 1;
  double OneOverNEvents = 1./NEvents;

  // Normalize to Number of Events
  if (hPtEPAnglePionAcc) hPtEPAnglePionAcc->Scale(OneOverNEvents);
  if (hPtEPAngleMCPion) hPtEPAngleMCPion->Scale(OneOverNEvents);
  if (hPtEPAngleTrueRecMCPion) hPtEPAngleTrueRecMCPion->Scale(OneOverNEvents);
  if (hPtRPAnglePionAcc) hPtRPAnglePionAcc->Scale(OneOverNEvents);
  if (hPtRPAngleMCPion) hPtRPAngleMCPion->Scale(OneOverNEvents);
  if (hPtRPAngleTrueRecMCPion) hPtRPAngleTrueRecMCPion->Scale(OneOverNEvents);

  if (hHistTrackPsiEPPtCent) hHistTrackPsiEPPtCent->Scale(OneOverNEvents);
  if (hHistTrackPsiEP3PtCent) hHistTrackPsiEP3PtCent->Scale(OneOverNEvents);

  // For each of the above, project onto DeltaPsi 6 pt bins
  for (int i = 0; i < kUsedPi0TriggerPtBins; i++) {
    // skipping the 4-5 bin is hardcoded here
    double fMinPt = Pi0PtBins[i+1];
    double fMaxPt = Pi0PtBins[i+2];

    int iMinBin = hPtEPAnglePionAcc->GetYaxis()->FindBin(fMinPt);
    int iMaxBin = hPtEPAnglePionAcc->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

    printf("Projecting Event Plane Histograms in bin from %.1f to %.1f\n",hPtEPAnglePionAcc->GetYaxis()->GetBinLowEdge(iMinBin),hPtEPAnglePionAcc->GetYaxis()->GetBinUpEdge(iMaxBin));

    TString sFormat = "%s_Proj_%d";
    TString sPtRange = Form("%.0f #leq #it{p}_{T} < %.0f GeV/#it{c}",fMinPt,fMaxPt);

    // Event Plane

    TH1F * hLocalPtEPAnglePionAcc_Proj = (TH1F *) hPtEPAnglePionAcc->ProjectionX(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
    hLocalPtEPAnglePionAcc_Proj->Sumw2();
    hLocalPtEPAnglePionAcc_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{EP} (%s)",sPtRange.Data()));
    hLocalPtEPAnglePionAcc_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
    hPtEPAnglePionAcc_Proj.push_back(hLocalPtEPAnglePionAcc_Proj);

    if (haveMCStatus && hPtEPAngleMCPion) {
      TH1F * hLocalPtEPAngleMCPion_Proj = (TH1F *) hPtEPAngleMCPion->ProjectionX(Form(sFormat.Data(),hPtEPAngleMCPion->GetName(),i),iMinBin,iMaxBin);
      hLocalPtEPAngleMCPion_Proj->SetTitle(Form("#pi_{0}^{MC True} #Delta#Psi_{EP} (%s)",sPtRange.Data()));
      hLocalPtEPAngleMCPion_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
      hPtEPAngleMCPion_Proj.push_back(hLocalPtEPAngleMCPion_Proj);

      if (hPtEPAngleTrueRecMCPion) {
        TH1F * hLocalPtEPAngleTrueRecMCPion_Proj = (TH1F *) hPtEPAngleTrueRecMCPion->ProjectionX(Form(sFormat.Data(),hPtEPAngleTrueRecMCPion->GetName(),i),iMinBin,iMaxBin);
        hLocalPtEPAngleTrueRecMCPion_Proj->SetTitle(Form("#pi_{0}^{MC True Rec} #Delta#Psi_{EP} (%s)",sPtRange.Data()));
        hLocalPtEPAngleTrueRecMCPion_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
        hPtEPAngleTrueRecMCPion_Proj.push_back(hLocalPtEPAngleTrueRecMCPion_Proj);
      }
    }

    // Reaction Plane
    if (hPtRPAnglePionAcc) {
      TH1F * hLocalPtRPAnglePionAcc_Proj = (TH1F *) hPtRPAnglePionAcc->ProjectionX(Form(sFormat.Data(),hPtRPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
      hLocalPtRPAnglePionAcc_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{RP} (%s)",sPtRange.Data()));
      hLocalPtRPAnglePionAcc_Proj->GetXaxis()->SetTitle("#Delta#Psi_{RP}"); // fixing typo
      hLocalPtRPAnglePionAcc_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
      hPtRPAnglePionAcc_Proj.push_back(hLocalPtRPAnglePionAcc_Proj);

      if (haveMCStatus && hPtRPAngleMCPion) {
        TH1F * hLocalPtRPAngleMCPion_Proj = (TH1F *) hPtRPAngleMCPion->ProjectionX(Form(sFormat.Data(),hPtRPAngleMCPion->GetName(),i),iMinBin,iMaxBin);
        hLocalPtRPAngleMCPion_Proj->SetTitle(Form("#pi_{0}^{MC True} #Delta#Psi_{RP} (%s)",sPtRange.Data()));
        hLocalPtRPAngleMCPion_Proj->GetXaxis()->SetTitle("#Delta#Psi_{RP}"); // fixing typo
        hLocalPtRPAngleMCPion_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
        hPtRPAngleMCPion_Proj.push_back(hLocalPtRPAngleMCPion_Proj);

        if (hPtRPAngleTrueRecMCPion) {
          TH1F * hLocalPtRPAngleTrueRecMCPion_Proj = (TH1F *) hPtRPAngleTrueRecMCPion->ProjectionX(Form(sFormat.Data(),hPtRPAngleTrueRecMCPion->GetName(),i),iMinBin,iMaxBin);
          hLocalPtRPAngleTrueRecMCPion_Proj->SetTitle(Form("#pi_{0}^{MC True Rec} #Delta#Psi_{RP} (%s)",sPtRange.Data()));
          hLocalPtRPAngleTrueRecMCPion_Proj->GetXaxis()->SetTitle("#Delta#Psi_{RP}"); // fixing typo
          hLocalPtRPAngleTrueRecMCPion_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
          hPtRPAngleTrueRecMCPion_Proj.push_back(hLocalPtRPAngleTrueRecMCPion_Proj);
        }
      }
    } 
  }

  // Tracks
  printf("Starting Track vs EP projections.\n");
//  const int kNTrackPtBins = 11;
//  std::vector <double> fTrackPtBins = {0.15, 0.25, 0.5,   1, 1.5, 2,   3, 4, 5,   6, 8, 10};

  if (hHistTrackPsiEPPtCent) {
    hHistTrackPsiEPPtCent->GetZaxis()->SetRange(iThetaModelCent+1,iThetaModelCent+1);
    hHistTrackPsiEPPt = (TH2F *) hHistTrackPsiEPPtCent->Project3D("yx");
    hHistTrackPsiEPPt->SetTitle("Track #it{p}_{T} vs #Delta#Psi_{EP}");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      printf("Projecting Tracks Event Plane Histograms in bin from %.1f to %.1f\n",hPtEPAnglePionAcc->GetYaxis()->GetBinLowEdge(iMinBin),hPtEPAnglePionAcc->GetYaxis()->GetBinUpEdge(iMaxBin));


      TString sFormat = "%s_Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

//(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
      hHistTrackPsiEPPtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEPAngleTrack_Proj = (TH1F *) hHistTrackPsiEPPtCent->Project3D("xe"); 
      hLocalPtEPAngleTrack_Proj->SetName(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i));
      //hLocalPtEPAnglePionAcc_Proj->Sumw2();
      hLocalPtEPAngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP} (%s)",sPtRange.Data()));
      hLocalPtEPAngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEPAngleTrack_Proj.push_back(hLocalPtEPAngleTrack_Proj);
    }
  }

  if (hHistTrackPsiEP3PtCent) {
    printf("Found the TH3 for tracks relative to 3rd order event plane\n");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      TString sFormat = "%s_EP3Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

      hHistTrackPsiEP3PtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEP3AngleTrack_Proj = (TH1F *) hHistTrackPsiEP3PtCent->Project3D("xe");
      hLocalPtEP3AngleTrack_Proj->SetName(Form(sFormat.Data(),hHistTrackPsiEP3PtCent->GetName(),i));
      //hLocalPtEP3AnglePionAcc_Proj->Sumw2();
      hLocalPtEP3AngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP,3} (%s)",sPtRange.Data()));
      hLocalPtEP3AngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEP3AngleTrack_Proj.push_back(hLocalPtEP3AngleTrack_Proj);

    }
  }
  if (hHistTrackPsiRPPtCent) {
    hHistTrackPsiRPPtCent->GetZaxis()->SetRange(iThetaModelCent+1,iThetaModelCent+1);
    hHistTrackPsiRPPt = (TH2F *) hHistTrackPsiRPPtCent->Project3D("yx");
    hHistTrackPsiRPPt->SetTitle("Track #it{p}_{T} vs #Delta#Psi_{RP}");
  }
}

void PionID::MeasureVn() {
  TCanvas * cVn = new TCanvas("cVn","cVn",fCanvasWidth,fCanvasHeight);

  // Event Plane Resolution
  // Copypasta'd from phase 4
  Double_t fEPRes_Set_0[4][6]  = {
    {  0.765960, 0.619163,  0.509267, 0.348666, 0.318429, 0.187868},
    {  0.858157, 0.822691, 0.692985, 0.580624,  0.502229,  0.375755},
    {  0.832549,  0.771133,  0.639423,  0.507014,  0.439729,  0.305388},
    {  0.704550,  0.445893,  0.380824,  0.196809,  0.211605,  0.084895}};
  Double_t fEPRes_R2 = fEPRes_Set_0[iThetaModelCent][1];
  Double_t fEPRes_R3 = fEPRes_Set_0[iThetaModelCent][2]; //
  Double_t fEPRes_R4 = fEPRes_Set_0[iThetaModelCent][3];
  Double_t fEPRes_R6 = fEPRes_Set_0[iThetaModelCent][5];

  gTrigger_Bv = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_Bv->SetName("Trigger_Bv");
  gTrigger_Bv->SetTitle("B Value from V_{n} fit");
  gTrigger_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_Bv->GetYaxis()->SetTitle("B");

  gTrigger_V2 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V2->SetName("Trigger_V2");
  gTrigger_V2->SetTitle("Calculated #tilde{v}_{2}^{trigger} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrigger_V2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V2->GetYaxis()->SetTitle("#tilde{v}_{2}");

  gTrigger_V4 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V4->SetName("Trigger_V4");
  gTrigger_V4->SetTitle("Calculated #tilde{v}_{4}^{trigger} (Event Plane method)");
  gTrigger_V4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V4->GetYaxis()->SetTitle("#tilde{v}_{4}");

  gTrigger_V6 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V6->SetName("Trigger_V6");
  gTrigger_V6->SetTitle("Calculated #tilde{v}_{6}^{trigger} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrigger_V6->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V6->GetYaxis()->SetTitle("#tilde{v}_{6}");


  gTrack_Bv = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrack_Bv->SetName("Track_Bv");
  gTrack_Bv->SetTitle("B Value from V_{n} fit");
  gTrack_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_Bv->GetYaxis()->SetTitle("B");

  gTrack_V2 = new TGraphErrors(kNTrackPtBins);
  gTrack_V2->SetName("Track_V2");
  gTrack_V2->SetTitle("Calculated #tilde{v}_{2}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrack_V2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V2->GetYaxis()->SetTitle("#tilde{v}_{2}");

  gTrack_V3 = new TGraphErrors(kNTrackPtBins);
  gTrack_V3->SetName("Track_V3");
  gTrack_V3->SetTitle("Calculated #tilde{v}_{3}^{Track} (Event Plane method)");
  gTrack_V3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V3->GetYaxis()->SetTitle("#tilde{v}_{3}");


  gTrack_V4 = new TGraphErrors(kNTrackPtBins);
  gTrack_V4->SetName("Track_V4");
  gTrack_V4->SetTitle("Calculated #tilde{v}_{4}^{Track} (Event Plane method)");
  gTrack_V4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V4->GetYaxis()->SetTitle("#tilde{v}_{4}");

  gTrack_V6 = new TGraphErrors(kNTrackPtBins);
  gTrack_V6->SetName("Track_V6");
  gTrack_V6->SetTitle("Calculated #tilde{v}_{6}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrack_V6->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V6->GetYaxis()->SetTitle("#tilde{v}_{6}");

  for (int i = 0; i < kUsedPi0TriggerPtBins; i++) {
    // skipping the 4-5 bin is hardcoded here
    double fMinPt = Pi0PtBins[i+1];
    double fMaxPt = Pi0PtBins[i+2];

    TH1F * hPionEP = hPtEPAnglePionAcc_Proj[i];
 //   TF1 * fitPionEP = new TF1(Form("Pion_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x))",hPionEP->GetXaxis()->GetXmin(),hPionEP->GetXaxis()->GetXmax());
   // TF1 * fitPionEP = new TF1(Form("Pion_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[3]*TMath::Cos(3*x)+2*[2]*TMath::Cos(4*x))",hPionEP->GetXaxis()->GetXmin(),hPionEP->GetXaxis()->GetXmax());
    TF1 * fitPionEP = new TF1(Form("Pion_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",hPionEP->GetXaxis()->GetXmin(),hPionEP->GetXaxis()->GetXmax());
    fitPionEP->SetParName(0,"B");
    fitPionEP->SetParName(1,"v_2");
    fitPionEP->SetParName(2,"v_4");
    fitPionEP->SetParName(3,"v_6");

    hPionEP->Fit(fitPionEP);
    fitPionEP->SetLineColor(kCyan);

    hPionEP->Draw();
    fitPionEP->Draw("SAME");

    gTrigger_Bv->SetPoint(i,(fMinPt+fMaxPt)/2.,fitPionEP->GetParameter(0));
    gTrigger_V2->SetPoint(i,(fMinPt+fMaxPt)/2.,fitPionEP->GetParameter(1)/fEPRes_R2);
    gTrigger_V4->SetPoint(i,(fMinPt+fMaxPt)/2.,fitPionEP->GetParameter(2)/fEPRes_R4);
    gTrigger_V6->SetPoint(i,(fMinPt+fMaxPt)/2.,fitPionEP->GetParameter(3)/fEPRes_R6);

    gTrigger_Bv->SetPointError(i,(fMaxPt-fMinPt)/2.,fitPionEP->GetParError(0));
    gTrigger_V2->SetPointError(i,(fMaxPt-fMinPt)/2.,fitPionEP->GetParError(1)/fEPRes_R2);
    gTrigger_V4->SetPointError(i,(fMaxPt-fMinPt)/2.,fitPionEP->GetParError(2)/fEPRes_R4);
    gTrigger_V6->SetPointError(i,(fMaxPt-fMinPt)/2.,fitPionEP->GetParError(3)/fEPRes_R6);

    cVn->Print(Form("%s/EPStudy_Pt_%.0f_%.0f.pdf",sOutputDir.Data(),fMinPt,fMaxPt));
    cVn->Print(Form("%s/CFiles/EPStudy_Pt_%.0f_%.0f.C",sOutputDir.Data(),fMinPt,fMaxPt));
  }

  cVn->Clear();

  gTrigger_V2->SetLineColor(kViolet+10);
  gTrigger_V2->SetMarkerColor(kViolet+10);
  gTrigger_V2->SetMarkerStyle(kFullDiamond);

  gTrigger_V4->SetLineColor(kCyan-3);
  gTrigger_V4->SetMarkerColor(kCyan-3);
  gTrigger_V4->SetMarkerStyle(kFullDoubleDiamond);

  gTrigger_V6->SetLineColor(kOrange+10);
  gTrigger_V6->SetMarkerColor(kOrange+10);
  gTrigger_V6->SetMarkerStyle(kFullTriangleUp);

  gTrigger_V2->Draw();
  gTrigger_V4->Draw("SAME");
  gTrigger_V6->Draw("SAME");

  TLegend * legVn = new TLegend(0.6,0.17,0.85,0.27);
  legVn->AddEntry(gTrigger_V2,"#tilde{v}_{2} (#pi^{0})","lp");
  legVn->AddEntry(gTrigger_V4,"#tilde{v}_{4} (#pi^{0})","lp");
  legVn->AddEntry(gTrigger_V6,"#tilde{v}_{6} (#pi^{0})","lp");

  float maxY = 0.05;
  float minY = -0.005;

  Double_t fMinX,fMinY,fMaxX,fMaxY;

  gTrigger_V2->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  gTrigger_V4->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  maxY += 0.05 * (maxY - minY);
  minY -= 0.1 * (maxY - minY);

  // GetEYHigh()[] ?
//  maxY = fmax(maxY,gTrigger_V2->GetMaximum());
///  maxY = fmax(maxY,gTrigger_V4->GetMaximum());
//  maxY = 1.3 * maxY;
//  minY = fmin(minY,-fabs(maxY));

  gTrigger_V2->GetYaxis()->SetRangeUser(minY,maxY);

  legVn->Draw("SAME");

  cVn->Print(Form("%s/EPStudy_V2_V4.pdf",sOutputDir.Data()));
  cVn->Print(Form("%s/CFiles/EPStudy_V2_V4.C",sOutputDir.Data()));

  TH1F * hTrackEP = 0;
  for (int i = 0; i < kNTrackPtBins; i++) {
    double fMinPt = fTrackPtBins[i];
    double fMaxPt = fTrackPtBins[i+1];

    hTrackEP = hPtEPAngleTrack_Proj[i];
    if (!hTrackEP) return;
   
    //TF1 * fitTrackEP = new TF1(Form("Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
    //TF1 * fitTrackEP = new TF1(Form("Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[3]*TMath::Cos(3*x)+2*[2]*TMath::Cos(4*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
    TF1 * fitTrackEP = new TF1(Form("Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
    fitTrackEP->SetParameter(0,hTrackEP->Integral("width") / (TMath::Pi() / 2));
    fitTrackEP->SetParameter(1,0.01);
    fitTrackEP->SetParameter(2,0.001);
    fitTrackEP->SetParameter(3,0.);

    fitTrackEP->SetParName(0,"B");
    fitTrackEP->SetParName(1,"v_2");
    fitTrackEP->SetParName(2,"v_4");
    fitTrackEP->SetParName(3,"v_6");

    hTrackEP->Fit(fitTrackEP);
    fitTrackEP->SetLineColor(kCyan);

    hTrackEP->Draw();
    fitTrackEP->Draw("SAME");

    gTrack_Bv->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(0));
    gTrack_V2->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEPRes_R2);
    gTrack_V4->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(2)/fEPRes_R4);
    gTrack_V6->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(3)/fEPRes_R6);

    gTrack_Bv->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(0));
    gTrack_V2->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEPRes_R2);
    gTrack_V4->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(2)/fEPRes_R4);
    gTrack_V6->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(3)/fEPRes_R6);

    cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f.pdf",sOutputDir.Data(),fMinPt,fMaxPt));
    cVn->Print(Form("%s/CFiles/EPStudy_Track_Pt_%.2f_%.2f.C",sOutputDir.Data(),fMinPt,fMaxPt));
  }
  cVn->Clear();
  legVn->Clear();

  // Calculate V3

  if (hHistTrackPsiEP3PtCent) {
    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      hTrackEP = hPtEP3AngleTrack_Proj[i];
      if (!hTrackEP) return;
      printf("Doing the v3 thing with histogram %s (%s)\n",hTrackEP->GetName(),hTrackEP->GetTitle());
      TF1 * fitTrackEP = new TF1(Form("Track_V3Fit_%d",i),"[0]*(1+2*[1]*TMath::Cos(3*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
      fitTrackEP->SetParameter(0,hTrackEP->Integral("width") / (TMath::Pi() / 2));
      fitTrackEP->SetParameter(1,0.01);

      fitTrackEP->SetParName(0,"B");
      fitTrackEP->SetParName(1,"v_3");

      hTrackEP->Fit(fitTrackEP);
      fitTrackEP->SetLineColor(kCyan);

      hTrackEP->Draw();
      fitTrackEP->Draw("SAME");

      gTrack_V3->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEPRes_R3);
      gTrack_V3->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEPRes_R3);

      cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f_v3.pdf",sOutputDir.Data(),fMinPt,fMaxPt));
      cVn->Print(Form("%s/CFiles/EPStudy_Track_Pt_%.2f_%.2f_v3.C",sOutputDir.Data(),fMinPt,fMaxPt));
    }

  }



  cVn->Clear();
  legVn->Clear();



  gTrack_V2->SetLineColor(kAzure-2);
  gTrack_V2->SetMarkerColor(kAzure-2);
  gTrack_V2->SetMarkerStyle(kOpenDiamond);

  gTrack_V4->SetLineColor(kCyan-6);
  gTrack_V4->SetMarkerColor(kCyan-6);
  gTrack_V4->SetMarkerStyle(kOpenDoubleDiamond);

  gTrack_V6->SetLineColor(kOrange+4);
  gTrack_V6->SetMarkerColor(kOrange+4);
  gTrack_V6->SetMarkerStyle(kOpenTriangleUp);

  legVn->AddEntry(gTrack_V2,"#tilde{v}_{2} (Track)","lp");
  legVn->AddEntry(gTrack_V4,"#tilde{v}_{4} (Track)","lp");
  legVn->AddEntry(gTrack_V6,"#tilde{v}_{6} (Track)","lp");

  gTrack_V2->Draw();
  gTrack_V4->Draw("SAME");
  gTrack_V6->Draw("SAME");
  legVn->Draw("SAME");

  maxY = 0.05;
  minY = -0.005;

  gTrack_V2->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  gTrack_V4->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  maxY += 0.05 * (maxY - minY);
  minY -= 0.1 * (maxY - minY);
  //maxY = 0.05;
  //maxY = fmax(maxY,gTrack_V2->GetMaximum());
  //maxY = fmax(maxY,gTrack_V4->GetMaximum());
  //maxY = 1.3 * maxY;
  //minY = fmin(-0.005,-fabs(maxY));

  gTrack_V2->GetYaxis()->SetRangeUser(minY,maxY);

  cVn->Print(Form("%s/EPStudy_Track_V2_V4.pdf",sOutputDir.Data()));
  cVn->Print(Form("%s/CFiles/EPStudy_Track_V2_V4.C",sOutputDir.Data()));


  // Draw just the V3 of tracks
  cVn->Clear();
  legVn->Clear();



  gTrack_V3->SetLineColor(kOrange+4);
  gTrack_V3->SetMarkerColor(kOrange+4);
  gTrack_V3->SetMarkerStyle(kOpenTriangleUp);
  gTrack_V3->Draw("ALP");

  cVn->Print(Form("%s/EPStudy_Track_V3.pdf",sOutputDir.Data()));
  cVn->Print(Form("%s/CFiles/EPStudy_Track_V3.C",sOutputDir.Data()));


  // Draw a combined plot to see if they match
  cVn->Clear();
  legVn->Clear();

  // FIXME use a TMultiGraph
  TMultiGraph * mg = new TMultiGraph();
  mg->Add(gTrigger_V2,"lp");
  mg->Add(gTrigger_V4,"lp");
  mg->Add(gTrigger_V6,"lp");

  mg->Add(gTrack_V2,"lp");
  mg->Add(gTrack_V4,"lp");
  mg->Add(gTrack_V6,"lp");

  legVn->AddEntry(gTrigger_V2,"#tilde{v}_{2} (#pi^{0})","lp");
  legVn->AddEntry(gTrack_V2,"#tilde{v}_{2} (Track)","lp");
  legVn->AddEntry(gTrigger_V4,"#tilde{v}_{4} (#pi^{0})","lp");
  legVn->AddEntry(gTrack_V4,"#tilde{v}_{4} (Track)","lp");
  legVn->AddEntry(gTrigger_V6,"#tilde{v}_{6} (#pi^{0})","lp");
  legVn->AddEntry(gTrack_V6,"#tilde{v}_{6} (Track)","lp");

/*
  gTrigger_V2->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  gTrigger_V4->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  gTrack_V2->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  gTrack_V4->ComputeRange(fMinX,fMinY,fMaxX,fMaxY);
  maxY = fmax(maxY,fMaxY);
  minY = fmin(minY,fMinY);
  maxY += 0.05 * (maxY - minY);
  minY -= 0.1 * (maxY - minY);
  */
//  gTrack_V2->GetXaxis()->SetRangeUser(0.0,20.0);
//  gTrack_V2->GetYaxis()->SetRangeUser(minY,maxY);

  mg->SetTitle("#tilde{v}_{n} (#pi^{0}_{Cand} Triggers and Charged Tracks)");
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  mg->GetYaxis()->SetTitle("#tilde{v}_{n}");
//  gTrack_V2->Draw();
//  gTrack_V4->Draw("SAME LP");
//  gTrigger_V2->Draw("SAME LP");
//  gTrigger_V4->Draw("SAME LP");

  minY = mg->GetYaxis()->GetXmin();
  maxY = mg->GetYaxis()->GetXmax();
  maxY += 0.05 * (maxY - minY);
  minY -= 0.10 * (maxY - minY);

  mg->GetYaxis()->SetRangeUser(minY,maxY);

  legVn->SetY1NDC(0.65);
  legVn->SetY2NDC(0.875);
  legVn->Draw("SAME");

  cVn->Print(Form("%s/EPStudy_Both_V2_V4.pdf",sOutputDir.Data()));
  cVn->Print(Form("%s/CFiles/EPStudy_Both_V2_V4.C",sOutputDir.Data()));
  printf("Done with measuring Vn\n");
}

void PionID::DoProjections() {
  cout<<"Projecting THnSparse"<<endl;
//void PionID::Pi0MassAnalysis() {
//  cout<<"Pi0 Mass Analysis"<<endl;

  for (int i = 0; i < nPtBins; i++) {
    double low_pt = Pi0PtBins[i];
    double high_pt = Pi0PtBins[i+1];

    double mean_pt = 0.5 * (Pi0PtBins[i] + Pi0PtBins[i+1]);

    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("++         Projecting %05.2f < pT < %05.2f GeV/#it{c}           ++\n",low_pt,high_pt);
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    

    // Any reweighting of MC background components could be done here
    //if (haveMCStatus && fMCReScaling != 1) {

    //}

    hInvarMasspTBin.push_back(fInvarMasspT->ProjectionX(TString::Format("MassPtBin_%d",i),fInvarMasspT->GetYaxis()->FindBin(low_pt),fInvarMasspT->GetYaxis()->FindBin(high_pt)));
    hInvarMasspTBin[i]->SetTitle(TString::Format("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));

    TH1D * localPairPt = (TH1D * ) hPairPt->Clone(TString::Format("PairPt_%d",i));
    localPairPt->GetXaxis()->SetRangeUser(low_pt,high_pt);

    if (haveMCStatus) {
      std::vector<TH1D *> localMCIdArray = {};
      for (Int_t k = 0; k < nMCId; k++) {
        localMCIdArray.push_back(fInvarMassPtMCId[k]->ProjectionX(Form("MassPtBin_MC%s_%d",sMCIdNames[k].Data(),i),fInvarMassPtMCId[k]->GetYaxis()->FindBin(low_pt),fInvarMassPtMCId[k]->GetYaxis()->FindBin(high_pt)));
        localMCIdArray[k]->Rebin(nRebinMass);
        localMCIdArray[k]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        localMCIdArray[k]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*localMCIdArray[k]->GetXaxis()->GetBinWidth(2)));
      }
      hInvarMassPtBinMCId.push_back(localMCIdArray);
      hInvarMassPtBinMCNoPeak.push_back(fInvarMassPtMCNoPeak->ProjectionX(Form("MassPtBin_MCNoPeak_%d",i),fInvarMassPtMCNoPeak->GetYaxis()->FindBin(low_pt),fInvarMassPtMCNoPeak->GetYaxis()->FindBin(high_pt)));
      hInvarMassPtBinMCNoPeak[i]->Rebin(nRebinMass);
      hInvarMassPtBinMCNoPeak[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMassPtBinMCNoPeak[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMassPtBinMCNoPeak[i]->GetXaxis()->GetBinWidth(2)));

      hInvarMassPtBinMCNoEta.push_back(fInvarMassPtMCNoEta->ProjectionX(Form("MassPtBin_MCNoEta_%d",i),fInvarMassPtMCNoEta->GetYaxis()->FindBin(low_pt),fInvarMassPtMCNoEta->GetYaxis()->FindBin(high_pt)));
      hInvarMassPtBinMCNoEta[i]->Rebin(nRebinMass);
      hInvarMassPtBinMCNoEta[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMassPtBinMCNoEta[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMassPtBinMCNoEta[i]->GetXaxis()->GetBinWidth(2)));


      hInvarMasspTBinRotBkgMCPi0.push_back(fInvarMasspTRotBkgMCPi0->ProjectionX(Form("MassPtBin_RotBkgPi0_%d",i),fInvarMasspTRotBkgMCPi0->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCPi0->GetYaxis()->FindBin(high_pt)));
      hInvarMasspTBinRotBkgMCPi0[i]->Rebin(nRebinMass);
      hInvarMasspTBinRotBkgMCPi0[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMasspTBinRotBkgMCPi0[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCPi0[i]->GetXaxis()->GetBinWidth(2)));

      hInvarMasspTBinRotBkgMCEta.push_back(fInvarMasspTRotBkgMCEta->ProjectionX(Form("MassPtBin_RotBkgEta_%d",i),fInvarMasspTRotBkgMCEta->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCEta->GetYaxis()->FindBin(high_pt)));
      hInvarMasspTBinRotBkgMCEta[i]->Rebin(nRebinMass);
      hInvarMasspTBinRotBkgMCEta[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMasspTBinRotBkgMCEta[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCEta[i]->GetXaxis()->GetBinWidth(2)));

      
      

      if (Pi0Cands->GetAxis(iMCAxis)->GetNbins()>11) {
        // New Ones
        hInvarMasspTBinRotBkgMCPi0EnergyPair.push_back(fInvarMasspTRotBkgMCPi0EnergyPair->ProjectionX(Form("MassPtBin_RotBkgPi0EnergyPair_%d",i),fInvarMasspTRotBkgMCPi0EnergyPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCPi0EnergyPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->SetLineColor(kPSPi0EnergyColor);
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->SetMarkerColor(kPSPi0EnergyColor);
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->GetXaxis()->GetBinWidth(2)));

        hInvarMasspTBinRotBkgMCEtaEnergyPair.push_back(fInvarMasspTRotBkgMCEtaEnergyPair->ProjectionX(Form("MassPtBin_RotBkgEtaEnergyPair_%d",i),fInvarMasspTRotBkgMCEtaEnergyPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCEtaEnergyPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->SetLineColor(kPSEtaEnergyColor);
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->SetMarkerColor(kPSEtaEnergyColor);
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->GetXaxis()->GetBinWidth(2)));

        hInvarMasspTBinRotBkgMCPi0PosPair.push_back(fInvarMasspTRotBkgMCPi0PosPair->ProjectionX(Form("MassPtBin_RotBkgPi0PosPair_%d",i),fInvarMasspTRotBkgMCPi0PosPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCPi0PosPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->SetLineColor(kPSPi0PosColor);
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->SetMarkerColor(kPSPi0PosColor);
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCPi0PosPair[i]->GetXaxis()->GetBinWidth(2)));

        hInvarMasspTBinRotBkgMCEtaPosPair.push_back(fInvarMasspTRotBkgMCEtaPosPair->ProjectionX(Form("MassPtBin_RotBkgEtaPosPair_%d",i),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->SetLineColor(kPSEtaPosColor);
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->SetMarkerColor(kPSEtaPosColor);
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMCEtaPosPair[i]->GetXaxis()->GetBinWidth(2)));

        hInvarMasspTBinRotBkgMinusMCPi0PosPair.push_back(fInvarMasspTRotBkgMinusMCPi0Pos->ProjectionX(Form("MassPtBin_RotBkgMinusMCPi0Pos_%d",i),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->SetLineColor(kRotBkgMinusMCPi0PosColor);
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->SetMarkerColor(kRotBkgMinusMCPi0PosColor);
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->SetMarkerStyle(kRotBkgMinusMCPi0PosStyle);
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->GetXaxis()->GetBinWidth(2)));

        hInvarMasspTBinRotBkgMinusMCAll.push_back(fInvarMasspTRotBkgMinusMCAll->ProjectionX(Form("MassPtBin_RotBkgMinusMCAll_%d",i),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgMCEtaPosPair->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinRotBkgMinusMCAll[i]->Rebin(nRebinMass);
        hInvarMasspTBinRotBkgMinusMCAll[i]->SetLineColor(kRotBkgMinusMCAllColor);
        hInvarMasspTBinRotBkgMinusMCAll[i]->SetMarkerColor(kRotBkgMinusMCAllColor);
        hInvarMasspTBinRotBkgMinusMCAll[i]->SetMarkerStyle(kRotBkgMinusMCAllStyle);
        hInvarMasspTBinRotBkgMinusMCAll[i]->SetTitle(Form("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        hInvarMasspTBinRotBkgMinusMCAll[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBinRotBkgMinusMCAll[i]->GetXaxis()->GetBinWidth(2)));
      }

    }

    if (haveRotBkg) {
      if (bApplyOpeningAngleCorrection) {
        hInvarMasspTBinRotBkg.push_back(fInvarMasspTRotBkgAngleScaled->ProjectionX(TString::Format("%s_%d",fInvarMasspTRotBkg->GetName(),i),fInvarMasspTRotBkgAngleScaled->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgAngleScaled->GetYaxis()->FindBin(high_pt)));
       // hInvarMasspTBinRotBkg.push_back(fInvarMasspTRotBkgAngleScaled->ProjectionX(TString::Format("%s_%d",fInvarMasspTRotBkgAngleScaled->GetName(),i),fInvarMasspTRotBkgAngleScaled->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkgAngleScaled->GetYaxis()->FindBin(high_pt)));
      } else {

        hInvarMasspTBinRotBkg.push_back(fInvarMasspTRotBkg->ProjectionX(TString::Format("%s_%d",fInvarMasspTRotBkg->GetName(),i),fInvarMasspTRotBkg->GetYaxis()->FindBin(low_pt),fInvarMasspTRotBkg->GetYaxis()->FindBin(high_pt)));
      }

      hInvarMasspTBinRotBkg[i]->SetTitle(TString::Format("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMasspTBinRotBkg[i]->Rebin(nRebinMass);

      if (fPSFinalPeak != 0) {
        hInvarMasspTBinPSCorr.push_back(fPSFinalPeak->ProjectionX(TString::Format("%s_%d",fPSFinalPeak->GetName(),i),fPSFinalPeak->GetYaxis()->FindBin(low_pt),fPSFinalPeak->GetYaxis()->FindBin(high_pt)));
        hInvarMasspTBinPSCorr[i]->SetTitle(TString::Format("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
        // Normalizing:
        Double_t lIntegral = hInvarMasspTBinPSCorr[i]->Integral("width");
        if (lIntegral != 0) hInvarMasspTBinPSCorr[i]->Scale(1.0/lIntegral);
      }

    }

    hInvarMasspTBin[i]->Rebin(nRebinMass);
    hInvarMasspTBin[i]->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/c^{2}",1000.*hInvarMasspTBin[i]->GetXaxis()->GetBinWidth(2)));

    ////
    if (scaleBkg) {

      //TH1D * hMassBkg = hInvarMasspTBinRotBkg[i];
      TH1D * hMassBkg = hInvarMasspTBinRotBkg[i];
      TH1D * hMassSE  = hInvarMasspTBin[i];
      double scale = hMassBkg->Integral(hMassBkg->GetXaxis()->FindBin(bkgScaleMin),hMassBkg->GetXaxis()->FindBin(bkgScaleMax));
      if (scale >= 1e-9) {
        scale = hMassSE->Integral(hMassSE->GetXaxis()->FindBin(bkgScaleMin),hMassSE->GetXaxis()->FindBin(bkgScaleMax)) / scale;
      }
      hMassBkg->Scale(scale);
      // Also Scaling the RotBkg MC Pi0 (if applicable)
      if (haveMCStatus && hInvarMasspTBinRotBkgMCPi0[i]) {

        hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->Scale(scale);
        hInvarMasspTBinRotBkgMinusMCAll[i]->Scale(scale);

        // FIXME use twice scale bc contribution is counted twice (once for each pair tracking choice)
        hInvarMasspTBinRotBkgMCPi0[i]->Scale(2.*scale);
        hInvarMasspTBinRotBkgMCEta[i]->Scale(2.*scale);
  
        // Scaling the Pos Swapped Pi0 and Eta components
        // FIXME does this also need the two factor?
        hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->Scale(2*scale);
        hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->Scale(2*scale);
        hInvarMasspTBinRotBkgMCPi0PosPair[i]->Scale(2*scale);
        hInvarMasspTBinRotBkgMCEtaPosPair[i]->Scale(2*scale);
      }


      Pi0BkgScaleArr.push_back(scale); //  Be careful about not double filling this.
      Pi0BkgScaleArrUn.push_back(0.);  // FIXME could actually calculate this uncertainty.  Just saying.
    }
    ////

    if (haveRotBkg) {  // check skip points??
      hInvarMassPtBinRotSub.push_back((TH1D *) hInvarMasspTBin[i]->Clone(Form("%s_RotSub",hInvarMasspTBin[i]->GetName())));
      //hInvarMassPtBinRotSub[i] = (TH1F *) hInvarMasspTBin[i]->Clone(Form("%s_RotSub",hInvarMasspTBin[i]->GetName()));
      hInvarMassPtBinRotSub[i]->SetTitle(Form("#splitline{Mixed/Rot Background Subtracted}{%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}}",Pi0PtBins[i],Pi0PtBins[i+1]));
//    hInvarMasspTBin[i]->SetTitle(TString::Format("",));
//    local_hBkgSub->SetTitle(Form("#splitline{%s}{Background Fit Subtracted}",local_hBkgSub->GetTitle()));
      if (bkgChoice != 0) hInvarMassPtBinRotSub[i]->Add(hInvarMasspTBinRotBkg[i],-1);
    }

    printf("Bin Integral = %.3f\n",hInvarMasspTBin[i]->Integral());
    if (0 == hInvarMasspTBin[i]->Integral()) {
      printf("Empty Bin, skipping\n");
      nSkipPoints++;
      continue;
    } else if (iFirstRealBin == -1) {
      iFirstRealBin = i;
      printf("Setting first real bin to %d\n.",iFirstRealBin);
    }
    // Removing the high pt bins
    if (high_pt > kHighPtCut) {
      nSkipPoints++;
      continue;
    }

    //=================================================================
    // After this, objects have indices i = 0 ... nPtBins - nSkipBins
    //=================================================================

    ptPointsForTGraph.push_back(localPairPt->GetMean());
    ptErrorsForTGraph.push_back(0.5*(high_pt-low_pt));
  }
}


void PionID::AnalyzeMatchedTracks() {

  if (!fClusEnergyMatchedTracks) {
    fprintf(stderr,"Missing Cluster Energy vs Number of Matched Tracks histogram.\n");
    return;
  }
//  TH2F * fClusEnergyMatchedTracksNorm = (TH2F *) fClusEnergyMatchedTracks->Clone("ClusEnergyMatchedTracksNorm");


  TH2F * fClusEnergyMatchedTracksNorm = Normalize2DHistByRow(fClusEnergyMatchedTracks);
  //TH2F * fClusEnergyMatchedTracksNorm = Normalize2DHistByCol(fClusEnergyMatchedTracks);

  TCanvas * cClusterMatched = new TCanvas("cClusterMatched","cClusterMatched",1200,900);


  fClusEnergyMatchedTracks->Draw("COLZ");

  cClusterMatched->SetLogz();

  cClusterMatched->Print(Form("%s/MatchedClusterEnergy.pdf",sOutputDir.Data()));
  cClusterMatched->Print(Form("%s/MatchedClusterEnergy.png",sOutputDir.Data()));

  fClusEnergyMatchedTracksNorm->Draw("COLZ");
  fClusEnergyMatchedTracksNorm->GetXaxis()->SetRangeUser(0,10);


  cClusterMatched->Print(Form("%s/MatchedClusterEnergyNorm.pdf",sOutputDir.Data()));
  cClusterMatched->Print(Form("%s/MatchedClusterEnergyNorm.png",sOutputDir.Data()));
}

void PionID::Pi0MassAnalysis() {
  cout<<"Pi0 Mass Analysis"<<endl;

  for (int i = 0; i < nPtBins - nSkipPoints; i++) {

    double low_pt = Pi0PtBins[i];
    double high_pt = Pi0PtBins[i+1];

    double mean_pt = 0.5 * (Pi0PtBins[i] + Pi0PtBins[i+1]);

    TF1 * localPi0Fit = 0;
    // FIXME FIXME normalize BKG
      
    TH1D * fitTarget = 0;
  
    if (fitBkgSub) {
      printf("Fitting Background-Subtracted Histogram ... \n");
      fitTarget = hInvarMassPtBinRotSub[i];
      if (!fitTarget) {
        fprintf(stderr,"Error, missing bkg Sub histogram\n");
        return;
      }
    } else {
      printf("Fitting Histogram without Background-Subtraction ...\n");
      fitTarget = hInvarMasspTBin[i];
    }
  
    bool useEta = false;
    if (0.5*(high_pt+low_pt) > kEtaThreshold) {
      useEta = true;
      printf("======================\nIncluding Eta in Fit\n======================\n");
    }

    TH1D * localBackgroundHist = 0;  // Background histogram for hist.
    // If the background has already been subtracted, leave this 0: 
    if (bkgChoice != 0 && !fitBkgSub)  {// Leave localBackgroundHist as 0 if we want no bkg.
      // if fitBkgSub==1, then we have already subtracted the bkg
      localBackgroundHist = hInvarMasspTBinRotBkg[i];
      printf("Debug: passing a bkg histogram to use with scale.\n"); //FIXME
    }
    // For Temporary Guidance
    //  int fPeakChoice;
    // 0  --  Gaus
    // 1  --  ExpDecay-Gaus 
    // 2  --  Breit-Wigner
    // 3  --  ? Crystal Ball?
    //    int fBkgChoice; // Choice of Residual Fit Background
    // 0  --  Poly(1)
    // 1  --  Poly(2)
    // 2  --  Poly(3)
    // 3  --  Poly(4)
    // 4  --  Exp * Poly(2)
    // 5  --  Exp * Poly(3)


    double fLowPt = Pi0PtBins[i];
    double fMeanPt = (Pi0PtBins[i+1] + Pi0PtBins[i]) / 2.;
    double fThetaModelLambda = 0;
    double fThetaModelMPrime = 0; 
    if (bEnableThetaModel) {
      GetThetaModelParameters(fLowPt, fOpeningAngleCut, fThetaModelLambda, fThetaModelMPrime);
//      GetThetaModelParameters(fMeanPt, fOpeningAngleCut, fThetaModelLambda, fThetaModelMPrime);
      printf("    ptbin %d Got Model Parameters Lambda = %f (1/MeV) and M' = %f (MeV)\n",i,fThetaModelLambda,fThetaModelMPrime);
    } else {
      printf("Theta Model Not Enabled\n");
    } 

    ParamBkg_Functor * fFunctor = 0;
    TF1 * localPi0Peak = 0;
    TF1 * localUnmodPi0Fit = 0; // Fit function with theta model disabled
    TF1 * localPi0Bkg  = 0;
    TF1 * localEtaPeak = 0;

    ParamBkg_Functor * localUnmodFunctor = 0;   

    TF1 * fMCPreAnalysis_Pi0Fit = 0;

    if (bUseMCPreAnalysis) {
      if (haveMCStatus) { // If this is MC, don't need to worry about matching centrality
      // This is already the same centrality range
        fMCPreAnalysis_Pi0Fit = fMCPreAnalysisPi0Fits[0][i];
      } else {
        fMCPreAnalysis_Pi0Fit = fMCPreAnalysisPi0Fits[iMCPreAnalysis_Cent][i];
      }
    }

    switch (fitMethod) {
      case 0: // temporary designation for new method
        localPi0Fit = fitPi0Peak(fitTarget,Form("Pi0Mass_Fit_%d",i),&fFunctor,useEta, fitPeakMethod, fitBkgMethod, localBackgroundHist, fThetaModelLambda, fThetaModelMPrime, true, fMCPreAnalysis_Pi0Fit); //
        //localPi0Fit = fitPi0Peak(fitTarget,Form("Pi0Mass_Fit_%d",i),&fFunctor,useEta, fitPeakMethod, fitBkgMethod, localBackgroundHist); // Without Theta Model
      //  localPi0Fit = fitPi0Peak(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta, fitPeakMethod, fitBkgMethod); //
        // FIXME could add ability to get fit AND Functor back.  Struct for fit,functor pairs??
        
        if (bEnableThetaModel) {
         // localUnmodFunctor; 
          localUnmodPi0Fit = fitPi0Peak(fitTarget,Form("Unmod_Pi0Mass_Fit_%d",i),&localUnmodFunctor,useEta, fitPeakMethod, fitBkgMethod, localBackgroundHist, 0., 0., false); // run without fitting
          // Copy parameters
          for (int z = 0; z < localPi0Fit->GetNpar(); z++) {
            localUnmodPi0Fit->SetParameter(z,localPi0Fit->GetParameter(z));
            localUnmodPi0Fit->SetParError(z,localPi0Fit->GetParError(z));
          }
          localUnmodPi0Fit->SetRange(0,localUnmodPi0Fit->GetXmax());
        }
        printf("Debug:  Functor for %d\n",i);
        printf("Debug:  Functor has %d params.\n",fFunctor->GetNParams());

        break;
      // old Code
      case 1:
        localPi0Fit = fitPi0Peak_1(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 2:
        localPi0Fit = fitPi0Peak_2(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 3:
        localPi0Fit = fitPi0Peak_3(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 4:
        localPi0Fit = fitPi0Peak_4(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 6:
        localPi0Fit = fitPi0Peak_6(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 7:
        localPi0Fit = fitPi0Peak_7(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 8:
        localPi0Fit = fitPi0Peak_8(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
      case 5:
      default:
        localPi0Fit = fitPi0Peak_5(fitTarget,Form("Pi0Mass_Fit_%d",i),useEta);
        break;
    }

    if (!localPi0Fit) {
      fprintf(stderr,"Error: something bad happened with fit!\n");
      return;
    }
//    if (bEnableThetaModel) {
//      printf("Copied functor has lambda = %f and mPrime = %f\n",localUnmodFunctor.GetLambda(),localUnmodFunctor.GetMPrime());
//    }
    


    double fLocalBkgScale = 1.;
    double fLocalBkgScaleUn = 0.;

    if (fitMethod == 0) { // temporary designation for new method

      if (bkgChoice != 0) {
        fLocalBkgScale = localPi0Fit->GetParameter("Bkg_Scale");
        fLocalBkgScaleUn = localPi0Fit->GetParError(localPi0Fit->GetParNumber("Bkg_Scale"));
        // Important: Applying scale to bkg histogram 
        if ( !fitBkgSub ) {
          printf("Scaling BKG by %f \n",fLocalBkgScale);
          localBackgroundHist->Scale(fLocalBkgScale);
          // Note that the th1 inside the functor is 
          // a clone, not affected by this
        }
      }


      double lMinX = localPi0Fit->GetXmin();
      double lMaxX = localPi0Fit->GetXmax();
      localPi0Peak = new TF1(Form("%s_Peak",localPi0Fit->GetName()),fFunctor,lMinX,lMaxX,fFunctor->GetNParams());
      localPi0Bkg = new TF1(Form("%s_Bkg",localPi0Fit->GetName()),fFunctor,lMinX,lMaxX,fFunctor->GetNParams());

      CopyTF1Details(localPi0Peak, localPi0Fit);
      CopyTF1Details(localPi0Bkg,  localPi0Fit);

      if (useEta ) {
//        localEtaPeak = (TF1 *) localPi0Fit->Clone(Form("%s_EtaPeak",localPi0Fit->GetName()));
        localEtaPeak = new TF1(Form("%s_EtaPeak",localPi0Fit->GetName()),fFunctor,lMinX,lMaxX,fFunctor->GetNParams());
        CopyTF1Details(localEtaPeak,localPi0Fit);
        localEtaPeak->SetRange(kEtaDrawMin,kEtaDrawMax);
      }
    } else {
      localPi0Peak = (TF1 *) localPi0Fit->Clone(Form("%s_Peak",localPi0Fit->GetName()));
      localPi0Bkg  = (TF1 *) localPi0Fit->Clone(Form("%s_Bkg",localPi0Fit->GetName()));

      if (useEta ) {
        printf("debug should not be used anymore\n");
        localEtaPeak = (TF1 *) localPi0Fit->Clone(Form("%s_EtaPeak",localPi0Fit->GetName()));
        localEtaPeak->SetRange(kEtaDrawMin,kEtaDrawMax);
      }
    }

    TH1D * localTotalFit = (TH1D * ) hInvarMasspTBinRotBkg[i]->Clone(Form("%s_TotalFit",fitTarget->GetName()));
    TH1D * localTotalBkg = (TH1D * ) hInvarMasspTBinRotBkg[i]->Clone(Form("%s_TotalBkg",fitTarget->GetName()));

    int nPars = localPi0Fit->GetNpar();
    // Finding nPeakPars:  Need to include extra parameters for more complex peaks.
    //   Check for N or A, the first bkg parameters

    // N.B. This can be done in a more simpler manner if the functor is accessed.

    int nPeakPars = 0;
    for (int j = 0; j < nPars; j++) {
      TString parName = localPi0Peak->GetParName(j);
      if (parName.EqualTo("N") || parName.EqualTo("A") ) {
        nPeakPars = j;
        break;
      }
    }

    // Removing Pi0 Peak from Bkg
    localPi0Bkg->SetParameter("Y",0);
//    }
    if (useEta) {
      // Removing everything from Eta Function
      localEtaPeak->SetParameter("Y",0);
      localEtaPeak->SetParameter("N",0);
      // FIXME temp
      localEtaPeak->SetParameter("p6",0);
      localEtaPeak->SetParameter("p7",0);
      localEtaPeak->SetParameter("p8",0);
      localEtaPeak->SetParameter("p9",0);
      // FIXME mixing background histogram scaling parameter
      //for (int j = 0; j < nPars; j++) {
      //  localEtaPeak->SetParameter(j,0);
     // }

      // Removing the Eta peak from the Residual Background Fit.
      localPi0Bkg->SetParameter("Y_eta",0);
      // FIXME need to remove bkg from residual bkg fit if using new code?    
  //    if (!fitBkgSub) {
//        localPi0Bkg->SetParameter("Bkg_Scale",0);
//      }

      // Adding back the eta peak
      printf("Debug (i = %d) setting the eta peak to yield %f\n",i,localPi0Fit->GetParameter("Y_eta"));
      localEtaPeak->SetParameter("Y_eta",localPi0Fit->GetParameter("Y_eta"));
      localEtaPeak->SetParameter("mu_eta",localPi0Fit->GetParameter("mu_eta"));
      localEtaPeak->SetParameter("sigma_eta",localPi0Fit->GetParameter("sigma_eta"));
    }
    // Removing the background and eta from the Pi0 Peak function
    for (int j = nPeakPars; j < nPars; j++) {
      localPi0Peak->SetParameter(j,0);
    }

    //localTotalFit, localTotalBkg
    if (!fitBkgSub) {
      // Maybe if not subtracting them, we should be replacing the histograms data with the function??? 
      localTotalFit->Reset();
      localTotalBkg->Reset();
    }
    printf("Before adding local pi0bkg function, bkg hist has int = %f\n",localTotalBkg->Integral());//FIXME
    AddFunctionToHist(localTotalFit,localPi0Fit);
    AddFunctionToHist(localTotalBkg,localPi0Bkg); // Need to do this BEFORE setting bkg scale to 0
    printf(" After adding local pi0bkg function, bkg hist has int = %f\n",localTotalBkg->Integral());//FIXME
    // Removing Combinatorial Bkg Term;
    if (bkgChoice != 0 && !fitBkgSub) {
      localPi0Bkg->SetParameter("Bkg_Scale",0.); // localPi0Scale is now just residual background
    }
    fPi0Fit.push_back(localPi0Fit);
    if(localUnmodPi0Fit != 0) fUnmodPi0Fit.push_back(localUnmodPi0Fit);
    fPi0Peak.push_back(localPi0Peak);
    fPi0Bkg.push_back(localPi0Bkg);    
    fEtaPeak.push_back(localEtaPeak); // note: some might be 0

    printf("function has name %s\n",localPi0Bkg->GetName());
//    PrintParameters(localPi0Bkg,true,"Test");

    hTotalFit.push_back(localTotalFit);
    hTotalBkg.push_back(localTotalBkg);

    printf("\tExtracting parameters\n");
    // Extract Yields, means, and sigmas
    Pi0YieldArr.push_back(localPi0Fit->GetParameter(0));
    Pi0YieldArrUn.push_back(localPi0Fit->GetParError(0));
    Pi0MassArr.push_back(localPi0Fit->GetParameter(1));
    Pi0MassArrUn.push_back(localPi0Fit->GetParError(1));
    Pi0SigmaArr.push_back(localPi0Fit->GetParameter(2));
    Pi0SigmaArrUn.push_back(localPi0Fit->GetParError(2));

    // FIXME add additional parameters of interest (lambda, etc)

    //normalizede yield
    printf("\tNormalizing Yields\n");
    Pi0NormYieldArr.push_back(localPi0Fit->GetParameter(0)/(high_pt - low_pt));
    Pi0NormYieldArrUn.push_back(localPi0Fit->GetParError(0)/(high_pt - low_pt));

    if (localPi0Fit->GetNDF() != 0) {
      Pi0ChiSquareArr.push_back(localPi0Fit->GetChisquare()/localPi0Fit->GetNDF());
    } else {
      Pi0ChiSquareArr.push_back(1.);
    }
    double fMean = 0.05;
    double fSigma = 0.05;
    double fPi0MassCutLow = 0.1;
    double fPi0MassCutHigh = 0.3;
    double fSBMassCut0 = -0.05;
    double fSBMassCut1 = -0.05;
    double fSBMassCut2 = -0.05;
    double fSBMassCut3 = -0.05;
    double fSBMassCut4 = -0.05;
    switch (iFixedMassWindows) {
      default:
      case 0: // Fresh Mass Windows
        fMean = localPi0Fit->GetParameter(1);
        fSigma = localPi0Fit->GetParameter(2);
        fPi0MassCutLow  = localPi0Fit->GetParameter(1) - nSigma * localPi0Fit->GetParameter(2);
        if (nSigmaR < 0) fPi0MassCutHigh = localPi0Fit->GetParameter(1) + nSigma * localPi0Fit->GetParameter(2);
        else fPi0MassCutHigh = localPi0Fit->GetParameter(1) + nSigmaR * localPi0Fit->GetParameter(2);

      break;
      case 1: // Windows used in GA Correlation 3 (Trains 55,56)
        printf("DEBUG: Using fixed, preselected mass windows\n");
  //    Double_t fPi0MassFixedValue_3[4][9]
        //fMean  = fPi0MassFixedValue_3[iThetaModelCent][i]; // The windows used in MB
       // fSigma = fPi0SigmaFixedValue_3[iThetaModelCent][i];
        // Windows used in GA
        fMean  = fPi0MassFixedValue_5[iThetaModelCent][i];
        fSigma = fPi0SigmaFixedValue_5[iThetaModelCent][i];
        fPi0MassCutLow = fMean - nSigma * fSigma;
        if (nSigmaR < 0) fPi0MassCutHigh = fMean + nSigma * fSigma;
        else fPi0MassCutHigh = fMean + nSigmaR * fSigma;
      break;
      case 2: // Windows used in MB Correlation 3 (Trains 44ish)
        printf("DEBUG: Using fixed, preselected mass windows\n");
  //    Double_t fPi0MassFixedValue_3[4][9]
        //fMean  = fPi0MassFixedValue_3[iThetaModelCent][i]; // The windows used in MB
       // fSigma = fPi0SigmaFixedValue_3[iThetaModelCent][i];
        // Windows used in GA
        fMean  = fPi0MassFixedValue_4[iThetaModelCent][i];
        fSigma = fPi0SigmaFixedValue_4[iThetaModelCent][i];
        fPi0MassCutLow = fMean - nSigma * fSigma;
        if (nSigmaR < 0) fPi0MassCutHigh = fMean + nSigma * fSigma;
        else fPi0MassCutHigh = fMean + nSigmaR * fSigma;



    }

    fSBMassCut0 = fMean + 3 * fSigma;
    fSBMassCut4 = 0.5;
    double fSBSplit = 0.25*(fSBMassCut4 - fSBMassCut0);
    fSBMassCut1 = fSBMassCut0 + fSBSplit;
    fSBMassCut2 = fSBMassCut0 + 2.*fSBSplit;
    fSBMassCut3 = fSBMassCut0 + 3.*fSBSplit;

    Pi0MassCutLow.push_back(fPi0MassCutLow);
    Pi0MassCutHigh.push_back(fPi0MassCutHigh);

    SBMassCut0.push_back(fSBMassCut0);
    SBMassCut1.push_back(fSBMassCut1);
    SBMassCut2.push_back(fSBMassCut2);
    SBMassCut3.push_back(fSBMassCut3);
    SBMassCut4.push_back(fSBMassCut4);



    // Fit Background subtraction

    printf("\tEstimating Background\n");
    printf("FIXME: (pT = %f)  Mass Window = [  %f  -   %f   ]\n",mean_pt,fPi0MassCutLow,fPi0MassCutHigh);
    int iPi0MassCutLow = localTotalBkg->FindBin(fPi0MassCutLow);
    int iPi0MassCutHigh = localTotalBkg->FindBin(fPi0MassCutHigh);
    printf("Mass Window bins [ %d - %d ]\n",iPi0MassCutLow,iPi0MassCutHigh);

    double fLocalBkg = 0;
    double fLocalBkgUn = 0;
    // FIXME should this be localTotalBkg being integrated?

    fLocalBkg = localTotalBkg->IntegralAndError(iPi0MassCutLow,iPi0MassCutHigh,fLocalBkgUn);
    //fLocalBkg = localPi0Bkg->Integral(fPi0MassCutLow,fPi0MassCutHigh);
    //fLocalBkgUn = localPi0Bkg->IntegralError(fPi0MassCutLow,fPi0MassCutHigh);

    // Dividing background by bin width to get real count (or could use "width" on yield integral);
    double binWidth = hInvarMasspTBin[i]->GetXaxis()->GetBinWidth(hInvarMasspTBin[i]->GetXaxis()->FindBin(0.135));
    // FIXME check if this is needed for integral. unless it's pt bin
//    fLocalBkg = fLocalBkg / binWidth;
//    fLocalBkgUn = fLocalBkgUn / binWidth;
//    printf("Found bkg = %f \\pm %f.  Found yield = %f \\pm %f\n",fLocalBkg,fLocalBkgUn,fitBkgSub,fLocalIntegral,fLocalIntegralUn);


    Pi0BkgArr.push_back(fLocalBkg);
    Pi0BkgArrUn.push_back(fLocalBkgUn);

    double fLocalMCBkg = 0;
    double fLocalMCBkgUn = 0;
    if (haveMCStatus) {
//hInvarMassPtBinMCNoPeak[i]
      fLocalMCBkg = hInvarMassPtBinMCNoPeak[i]->IntegralAndError(hInvarMassPtBinMCNoPeak[i]->GetXaxis()->FindBin(fPi0MassCutLow),hInvarMassPtBinMCNoPeak[i]->GetXaxis()->FindBin(fPi0MassCutHigh),fLocalMCBkgUn);

      printf("LocalMCBkg = %f\n",fLocalMCBkg);
      MCBkgArr.push_back(fLocalMCBkg);
      MCBkgArrUn.push_back(fLocalMCBkgUn);
    }

    printf("\tSubtracting fit to produce residual.\n");
    TH1D * local_hResid = (TH1D *) hInvarMasspTBin[i]->Clone(Form("%s_Resid",hInvarMasspTBin[i]->GetName()));
    local_hResid->Add(localPi0Fit,-1);

    hInvarMasspTBinResid.push_back(local_hResid);

    printf("\tSubtracting background, integrating\n");
    TH1D * local_hBkgSub = (TH1D *) hInvarMasspTBin[i]->Clone(Form("%s_BkgSub",hInvarMasspTBin[i]->GetName()));

    // Should this be the background histo?
  //  local_hBkgSub->Add(localPi0Bkg,-1);
    local_hBkgSub->Add(localTotalBkg,-1);


    local_hBkgSub->SetTitle(Form("#splitline{Background Fit Subtracted}{%s}",local_hBkgSub->GetTitle()));
    hInvarMassBkgSubpTBin.push_back(local_hBkgSub);

    double fLocalIntegralUn = 0;
    double fLocalIntegral = local_hBkgSub->IntegralAndError(local_hBkgSub->GetXaxis()->FindBin(fPi0MassCutLow),local_hBkgSub->GetXaxis()->FindBin(fPi0MassCutHigh),fLocalIntegralUn);
   // double fLocalIntegral = local_hBkgSub->IntegralAndError(local_hBkgSub->GetXaxis()->FindBin(Pi0MassCutLow[i]),local_hBkgSub->GetXaxis()->FindBin(Pi0MassCutHigh[i]),fLocalIntegralUn);

    double fTotalIntegralUn = 0;
    double fTotalIntegral = hInvarMasspTBin[i]->IntegralAndError(hInvarMasspTBin[i]->GetXaxis()->FindBin(fPi0MassCutLow),hInvarMasspTBin[i]->GetXaxis()->FindBin(fPi0MassCutHigh),fTotalIntegralUn);

    printf("Found bkg = %.2e.  Found yield = %.2e.  Total = %.2e.  Test Total = %.2e\n",fLocalBkg,fLocalIntegral,fTotalIntegral,fLocalBkg+fLocalIntegral);

    Pi0IntegralArr.push_back(fLocalIntegral);
    Pi0IntegralArrUn.push_back(fLocalIntegralUn);
    Pi0IntegralTotalArr.push_back(fTotalIntegral);
    Pi0IntegralTotalArrUn.push_back(fTotalIntegralUn);
    Pi0NormIntegralArr.push_back(fLocalIntegral / (high_pt - low_pt));
    Pi0NormIntegralArrUn.push_back(fLocalIntegralUn / (high_pt - low_pt));

    double fLocalMCIntegral = 0, fLocalMCIntegralUn = 0;
    if (haveMCStatus) {
      // Could use Pi0 Mass Window from Fit, or fixed
      fLocalMCIntegral = hInvarMassPtBinMCId[i][2]->IntegralAndError(hInvarMassPtBinMCId[i][2]->GetXaxis()->FindBin(fPi0MassCutLow),hInvarMassPtBinMCId[i][2]->GetXaxis()->FindBin(fPi0MassCutHigh),fLocalMCIntegralUn);
      //fLocalMCIntegral = hInvarMassPtBinMCId[i][2]->IntegralAndError(hInvarMassPtBinMCId[i][2]->GetXaxis()->FindBin(fPi0MassCutLow),hInvarMassPtBinMCId[i][2]->GetXaxis()->FindBin(fPi0MassCutHigh),fLocalMCIntegralUn);
      //fLocalMCIntegral = hInvarMassPtBinMCId[i][2]->IntegralAndError(fLocalMCIntegralUn);

      Pi0MCIntegralArr.push_back(fLocalMCIntegral);
      Pi0MCIntegralArrUn.push_back(fLocalMCIntegralUn);
      Pi0MCNormIntegralArr.push_back(fLocalMCIntegral / (high_pt - low_pt));
      Pi0MCNormIntegralArrUn.push_back(fLocalMCIntegralUn / (high_pt - low_pt));
    }
//    printf("Estimating Yield-to-Bkg Ratio\n");
    double fLocalPeakSigRatio = 0;
    double fLocalPeakSigRatioUn = 0;

    double fLocalYieldBkgRatio = 0;
    double fLocalYieldBkgRatioUn = 0;

    double fLocalYieldTotalRatio = 0;
    double fLocalYieldTotalRatioUn = 0;

    if (fLocalBkg >= 1e-7) {
      fLocalPeakSigRatio = fLocalIntegral / TMath::Sqrt(fLocalBkg);
      fLocalYieldBkgRatio = fLocalIntegral / fLocalBkg;
      if (fLocalIntegral != 0) {
        fLocalPeakSigRatioUn = fLocalPeakSigRatio * TMath::Sqrt(0.25*TMath::Power(fLocalBkgUn/fLocalBkg,2) + TMath::Power(fLocalIntegralUn/fLocalIntegral,2));
        fLocalYieldBkgRatioUn = fLocalYieldBkgRatio * TMath::Sqrt(TMath::Power(fLocalBkgUn/fLocalBkg,2) + TMath::Power(fLocalIntegralUn/fLocalIntegral,2));
      }
    }
//    printf("FIXME: Y/B = %f \\pm %f\n",fLocalYieldBkgRatio, fLocalYieldBkgRatioUn);
    if (fTotalIntegral > 0 && fLocalIntegral > 0) {
      fLocalYieldTotalRatio = fLocalIntegral / fTotalIntegral;
      fLocalYieldTotalRatioUn = fLocalYieldTotalRatio * TMath::Sqrt(TMath::Power(fLocalIntegralUn/ fLocalIntegral,2) + TMath::Power(fTotalIntegralUn/fTotalIntegral,2));

    }

    //std::vector<double> Pi0YieldBkgRatioArr = {};

    Pi0PeakSigRatioArr.push_back(fLocalPeakSigRatio);
    Pi0PeakSigRatioArrUn.push_back(fLocalPeakSigRatioUn);

    Pi0YieldBkgRatioArr.push_back(fLocalYieldBkgRatio);
    Pi0YieldBkgRatioArrUn.push_back(fLocalYieldBkgRatioUn);

    Pi0YieldTotalRatioArr.push_back(fLocalYieldTotalRatio);
    Pi0YieldTotalRatioArrUn.push_back(fLocalYieldTotalRatioUn);

    if (haveMCStatus) {
      double fLocalMCPi0PeakSigRatio = 0, fLocalMCPi0PeakSigRatioUn = 0;  
      double fLocalMCPi0YieldBkgRatio = 0, fLocalMCPi0YieldBkgRatioUn = 0;
      double fLocalMCPi0YieldTotalRatio = 0, fLocalMCPi0YieldTotalRatioUn = 0;
      double fLocalRecMCPi0YieldRatio = 0, fLocalRecMCPi0YieldRatioUn = 0;

      if (fLocalMCBkg >= 1e-7) {
        fLocalMCPi0PeakSigRatio = fLocalMCIntegral / TMath::Sqrt(fLocalMCBkg);
        fLocalMCPi0YieldBkgRatio = fLocalMCIntegral / fLocalMCBkg;
        if (fLocalMCIntegral != 0) {
          fLocalMCPi0PeakSigRatioUn = fLocalMCPi0PeakSigRatio * TMath::Sqrt(0.25*TMath::Power(fLocalMCBkgUn/fLocalMCBkg,2) + TMath::Power(fLocalMCIntegralUn/fLocalMCIntegral,2));
          fLocalMCPi0YieldBkgRatioUn = fLocalMCPi0YieldBkgRatio * TMath::Sqrt(TMath::Power(fLocalMCBkgUn/fLocalMCBkg,2) + TMath::Power(fLocalMCIntegralUn/fLocalMCIntegral,2));
        }
      }

      if (fTotalIntegral > 0 && fLocalMCIntegral > 0) {
        fLocalMCPi0YieldTotalRatio = fLocalMCIntegral / fTotalIntegral;
        fLocalMCPi0YieldTotalRatioUn = fLocalMCPi0YieldTotalRatio * TMath::Sqrt(TMath::Power(fLocalMCIntegralUn/ fLocalMCIntegral,2) + TMath::Power(fTotalIntegralUn/fTotalIntegral,2));


        fLocalRecMCPi0YieldRatio = fLocalIntegral / fLocalMCIntegral;
        fLocalRecMCPi0YieldRatioUn = fLocalRecMCPi0YieldRatioUn * TMath::Sqrt(TMath::Power(fLocalIntegralUn/fLocalIntegral,2) + TMath::Power(fLocalMCIntegralUn/fLocalMCIntegral,2));

      }

      MCPi0PeakSigRatioArr.push_back(fLocalMCPi0PeakSigRatio);
      MCPi0PeakSigRatioArrUn.push_back(fLocalMCPi0PeakSigRatioUn);

      MCPi0YieldBkgRatioArr.push_back(fLocalMCPi0YieldBkgRatio);
      MCPi0YieldBkgRatioArrUn.push_back(fLocalMCPi0YieldBkgRatioUn);

      MCPi0YieldTotalRatioArr.push_back(fLocalMCPi0YieldTotalRatio);
      MCPi0YieldTotalRatioArrUn.push_back(fLocalMCPi0YieldTotalRatioUn);

      RecMCPi0YieldRatioArr.push_back(fLocalRecMCPi0YieldRatio);
      RecMCPi0YieldRatioArrUn.push_back(fLocalRecMCPi0YieldRatioUn);

    }

  }
  printf("Setting line features.  nPt = %d , nSkipPoints = %d\n",nPtBins,nSkipPoints);

  for (int i = 0; i < nPtBins-nSkipPoints; i++) {
    fPi0Fit[i]->SetLineColor(kViolet);
    fPi0Fit[i]->SetLineStyle(2);

    fPi0Peak[i]->SetLineColor(kPi0FitColor);
    fPi0Peak[i]->SetNpx(kFitNpx);
    fPi0Peak[i]->SetLineStyle(1);

    fPi0Bkg[i]->SetLineWidth(4);
    fPi0Bkg[i]->SetLineStyle(9);

    if (fEtaPeak[i]) {
      fEtaPeak[i]->SetLineColor(kEtaPeakLineColor);
      fEtaPeak[i]->SetLineStyle(kEtaPeakLineStyle);
      fEtaPeak[i]->SetNpx(kFitNpx);
    }

    if (bDrawUnmod && fUnmodPi0Fit.size() > i) {
      fUnmodPi0Fit[i]->SetLineColor(kUnmodPi0FitColor);
      fUnmodPi0Fit[i]->SetLineStyle(kUnmodPi0FitLineStyle);
      fUnmodPi0Fit[i]->SetLineWidth(kUnmodPi0FitLineWidth);
    }

    hTotalFit[i]->SetLineColor(kTotalFitColor);
    hTotalFit[i]->SetMarkerColor(kTotalFitColor);
    hTotalFit[i]->SetLineWidth(kTotalFitLineWidth);

    hTotalBkg[i]->SetLineColor(kTotalBkgColor);
    hTotalBkg[i]->SetMarkerColor(kTotalBkgColor);
    hTotalBkg[i]->SetLineWidth(kTotalBkgLineWidth);
    hTotalBkg[i]->SetLineStyle(kTotalBkgLineStyle);
    if (haveMCStatus) {
      for (Int_t k = 0; k < nMCId; k++) {
        hInvarMassPtBinMCId[i][k]->SetLineColor(kMCIdColor[k]);
        hInvarMassPtBinMCId[i][k]->SetMarkerColor(kMCIdColor[k]);
        hInvarMassPtBinMCId[i][k]->SetMarkerStyle(kMCIdStyle[k]);
      }

      hInvarMassPtBinMCNoPeak[i]->SetLineColor(kMCNoPeakColor);
      hInvarMassPtBinMCNoPeak[i]->SetMarkerColor(kMCNoPeakColor);
      hInvarMassPtBinMCNoPeak[i]->SetMarkerStyle(kMCNoPeakMarkerStyle);

      hInvarMassPtBinMCNoEta[i]->SetLineColor(kMCNoEtaColor);
      hInvarMassPtBinMCNoEta[i]->SetMarkerColor(kMCNoEtaColor);
      hInvarMassPtBinMCNoEta[i]->SetMarkerStyle(kMCNoEtaMarkerStyle);

      hInvarMasspTBinRotBkgMCPi0[i]->SetLineColor(kMCRotBkgPi0Color);
      hInvarMasspTBinRotBkgMCPi0[i]->SetMarkerColor(kMCRotBkgPi0Color);
      hInvarMasspTBinRotBkgMCPi0[i]->SetMarkerStyle(kMCRotBkgPi0MarkerStyle);

      hInvarMasspTBinRotBkgMCEta[i]->SetLineColor(kMCRotBkgEtaColor);
      hInvarMasspTBinRotBkgMCEta[i]->SetMarkerColor(kMCRotBkgEtaColor);
      hInvarMasspTBinRotBkgMCEta[i]->SetMarkerStyle(kMCRotBkgEtaMarkerStyle);
    }
  }
}


/** 
  * Fits the real pi0 peak, using whatever the input peak function fit
  */
void PionID::FitMCTruthPi0() {
  cout<<"=============================================="<<endl;
  cout<<"Fitting the Real Pi0 peak (from MC Truth info)"<<endl;
  cout<<"=============================================="<<endl;
  
  Double_t fMCPi0FitXMax = 0.33;

  // fMCPi0Fit
  for (int i = 0; i < nPtBins-nSkipPoints; i++) {
    TH1D * hLocalMCPi0 = hInvarMassPtBinMCId[i][2];
    TF1  * localPi0Fit = 0;

    double fLowPt = Pi0PtBins[i];
    double fMeanPt = (Pi0PtBins[i+1] + Pi0PtBins[i]) / 2.;
    double fThetaModelLambda = 0;
    double fThetaModelMPrime = 0; 

    if (bEnableThetaModel) {
      GetThetaModelParameters(fLowPt, fOpeningAngleCut, fThetaModelLambda, fThetaModelMPrime);
//      GetThetaModelParameters(fMeanPt, fOpeningAngleCut, fThetaModelLambda, fThetaModelMPrime);
      printf("MCFit ptbin %d Got Model Parameters Lambda = %f (1/MeV) and M' = %f (MeV)\n",i,fThetaModelLambda,fThetaModelMPrime);
    } else {
      printf("Theta Model not Enabled\n");
    }
    ParamBkg_Functor * fFunctor = 0;

    //hLocalMCPi0->Draw("e"); // literally just to get the axes working
    hLocalMCPi0->GetXaxis()->SetRangeUser(0.,fMCPi0FitXMax); // FIXME WHY IS THIS NOT WORKING
    printf("DEBUG giving MC Pi0 fit histogram with range %f %f\n",hLocalMCPi0->GetXaxis()->GetXmin(),hLocalMCPi0->GetXaxis()->GetXmax());
// Adding in the pi0 bin (mid-bin pt?) and theta cut
    localPi0Fit = fitPi0Peak(hLocalMCPi0,Form("MCPi0Mass_Fit_%d",i),&fFunctor,0,fitPeakMethod,/* FitBkgMethod = */0,/* Bkg Histogram =*/0, fThetaModelLambda, fThetaModelMPrime); 
//    localPi0Fit = fitPi0Peak(hLocalMCPi0,Form("MCPi0Mass_Fit_%d",i),&fFunctor,0,fitPeakMethod,/* FitBkgMethod = */0,/* Bkg Histogram =*/0); 
  
    if (!hLocalMCPi0) {
      fprintf(stderr,"Missing an MC Pi0 Histogram\n");
      continue;
    }
    
    if (localPi0Fit != 0) {
      fMCPi0Fit.push_back(localPi0Fit);
      if (localPi0Fit->GetNDF() != 0) {
        MCPi0ChiSquareArr.push_back(localPi0Fit->GetChisquare()/localPi0Fit->GetNDF());
      } else {
        MCPi0ChiSquareArr.push_back(1.);
      }
    }
  }


  // Save these as the MC Pre Analysis files. 
  // Only one centrality bin saved [whichever one being analyzed]
  vector<TF1 *> fMCPreAnalysisPi0Fits_OneBin = {};
  //for (int i = 0; i < 7 ; i++) { //load 4-5,5-7,7-9,9-11,11-14,14-17,17-20 FIXME magic number
  for (int i = 0; i < nPtBins-nSkipPoints; i++) {
    fMCPreAnalysisPi0Fits_OneBin.push_back(fMCPi0Fit[i]);
  }
  fMCPreAnalysisPi0Fits.push_back(fMCPreAnalysisPi0Fits_OneBin);
  cout<<"Done Fitting Real Pi0 Peaks"<<endl;
}

void PionID::DrawMassPlots() {
  cout<<"Creating primary output plots"<<endl;
  int cNY = TMath::FloorNint(TMath::Sqrt(nPtBins-nSkipPoints));
  int cNX = TMath::CeilNint(((float) nPtBins-nSkipPoints)/cNY);

  if (iThetaModelCent == 1) kMagicScale = 1.41; //1.28


  // Making a of raw mass plots
  TCanvas *cAllMass = new TCanvas("cMassAll","cMassAll",2400,1400);
  cAllMass->Divide(cNX,cNY,0.001,0.0012);
    //leg4->AddEntry(hInvarMasspTBin[i+iFirstRealBin],Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]),"pe");
  for (int i = 0; i < nPtBins-nSkipPoints; i++)
  {
    TLegend * legAllMass = new TLegend(0.40,0.70,0.87,0.8);
    legAllMass->SetTextSize(1.5*legAllMass->GetTextSize());
   // legAllMass->Clear();
    cAllMass->cd(i+1);
//    hInvarMasspTBin[i]->SetTitle(Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    hInvarMasspTBin[i]->Draw();
    legAllMass->SetHeader(Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
//    legAllMass->AddEntry(hInvarMasspTBin[i],Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]),"pe");
    legAllMass->Draw("SAME");
  }

  // Making the Specifc 2x2 plots for bins 7-9,  9-11, 11-14, 14-17
  //   Bins 3, 4, 5, 6

  TCanvas * c2x2 = new TCanvas("c2x2","c2x2",2400,2800);
  c2x2->Divide(2,2,0.001,0.0012);

  TCanvas * c2x1 = new TCanvas("c2x1","c2x1",2400,1400);
  c2x1->Divide(2,1,0.001,0.0012);
  TCanvas * c1x2 = new TCanvas("c1x2","c1x2",1200,2800);
  c1x2->Divide(1,2,0.001,0.0012);

  TCanvas * cMassFit = new TCanvas("cMassFit","cMassFit",600,700);

  TCanvas * cMassFitMCInfo = 0;
  if (haveMCStatus) {
    cMassFitMCInfo = new TCanvas("cMassFitMCInfo","cMassFitMCInfo",600,700);
  }

  TCanvas * c3 = new TCanvas("c3","c3",1200,900);
  // c3->Divide(cNX,cNY,0,0);
  //    c3->Divide(cNX,cNY,0.005,0.005);
  c3->Divide(cNX,cNY,0.001,0.0012);
  printf("Debug (cNX,cNY) = (%d,%d)\n",cNX,cNY);

  for (int i = 0; i < nPtBins-nSkipPoints; i++)
  {
    c3->cd(i+1);
    //  hInvarMasspTBin[i]->GetXaxis()->SetRangeUser(0.0,0.5);
    //    hInvarMasspTBin[i]->Draw();
    // For convenience;
//    cout<<"A0"<<endl;
    TH1D * hMassSE = hInvarMasspTBin[i]; //..ELI
    SetTH1Histo(hInvarMasspTBin[i],"M_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",1);
    hInvarMasspTBin[i+iFirstRealBin]->SetLineColor(kInvarMassColor);
    hInvarMasspTBin[i+iFirstRealBin]->SetMarkerColor(kInvarMassColor);
    hInvarMasspTBin[i+iFirstRealBin]->SetMarkerStyle(kInvarMassStyle); // 4
    hInvarMasspTBin[i+iFirstRealBin]->SetMarkerSize(0.7); //0.7

    if (bZoomInvarMass) hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetRangeUser(0,0.8);
    

    hInvarMasspTBin[i+iFirstRealBin]->SetFillColorAlpha(0,0);

    // The next two lines only work in the latest version of root
//    hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetMaxDigits(2);
//    hInvarMasspTBin[i+iFirstRealBin]->GetYaxis()->SetMaxDigits(2);
 //   hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetTitleSize(0.05);
 //   hInvarMasspTBin[i+iFirstRealBin]->GetYaxis()->SetTitleSize(0.05);
//    hInvarMasspTBin[i+iFirstRealBin]->GetYaxis()->SetTitleOffset(1.1); //1.3
 //   hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetTitleOffset(1.0); //1.1
    //hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetLabelSize(0.04);
 //   hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetLabelSize(0.04);
  


    hInvarMasspTBin[i+iFirstRealBin]->Draw("E");

    TH1D * hMassBkg = 0;
    if ( haveRotBkg) {
      hMassBkg = hInvarMasspTBinRotBkg[i+iFirstRealBin];

      hInvarMasspTBinRotBkg[i+iFirstRealBin]->SetLineColor(kCombBkgColor);
      hInvarMasspTBinRotBkg[i+iFirstRealBin]->SetMarkerColor(kCombBkgColor);
      hInvarMasspTBinRotBkg[i+iFirstRealBin]->SetMarkerStyle(kCombBkgMarkerStyle);
      hInvarMasspTBinRotBkg[i+iFirstRealBin]->SetLineWidth(kCombBkgLineWidth);
      if (drawBkg) hInvarMasspTBinRotBkg[i+iFirstRealBin]->Draw("SAME");

    }

    if (drawMCInfo && haveMCStatus) {

      for (Int_t k = 0; k < nMCId; k++) {
        hInvarMassPtBinMCId[i+iFirstRealBin][k]->Draw("SAME");
      }
      hInvarMassPtBinMCNoPeak[i+iFirstRealBin]->Draw("E SAME");
      hInvarMassPtBinMCNoEta[i+iFirstRealBin]->Draw("E SAME");
    }
    if (drawFits) {
      // Setting Draw Ranges
      if (hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->GetXmax()>0.750) {
        fitDrawRangeMax = hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->GetXmax();
      }
      fPi0Bkg[i]->SetRange(fitDrawRangeMin,fitDrawRangeMax);

      fPi0Peak[i]->SetRange(fitDrawRangeMin,fitDrawRangeMax);
      fPi0Fit[i]->SetRange(fitDrawRangeMin,fitDrawRangeMax);

      fPi0Bkg[i]->SetLineColor(kPi0BkgColor);
      fPi0Bkg[i]->SetMarkerColor(kPi0BkgColor);

      fPi0Fit[i]->SetLineColor(kPi0FitColor);
      fPi0Fit[i]->SetMarkerColor(kPi0FitColor);


      fPi0Peak[i]->Draw("SAME");

      // FIXME turn this off after devel done
      if (drawResBkg) {
        fPi0Bkg[i]->Draw("SAME");
      }

      if (fEtaPeak[i]) {
        printf("debug drawint eta for %s\n",hInvarMasspTBin[i+iFirstRealBin]->GetName());
        printf("debug drawing eta peak that has yield = %f, mean = %f, sigma = %f\n",fEtaPeak[i]->GetParameter("Y_eta"),fEtaPeak[i]->GetParameter("mu_eta"),fEtaPeak[i]->GetParameter("sigma_eta"));
        fEtaPeak[i]->Draw("SAME");
      }

      if (bDrawUnmod && fUnmodPi0Fit.size() > i) {
        fUnmodPi0Fit[i]->Draw("SAME");
      }

      if (fitBkgSub) hTotalBkg[i]->Draw("HIST SAME");
      hTotalFit[i]->Draw("HIST SAME");

      hInvarMasspTBin[i+iFirstRealBin]->Draw("E SAME");
    }

    TLegend* leg4;
    leg4 = new TLegend(0.53,0.45,0.7,0.75); //..Bkg subtracted //0.59,0.45,0.7,0.75
    leg4->AddEntry(hInvarMasspTBin[i+iFirstRealBin],Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]),"pe");
  //  leg4->AddEntry(hInvarMasspTBin[i],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]),"");
  //  leg4->AddEntry(fPi0Fit[i],"S+B Fit","l");

    leg4->AddEntry(hTotalFit[i],"S+B Fit","l");

    if (bDrawUnmod && fUnmodPi0Fit.size() > i) {
      leg4->AddEntry(fUnmodPi0Fit[i],"Unmodulated Fit","l");
    }

//    leg4->AddEntry(fPi0Bkg[i],"Background","l");
    if (fitBkgSub) leg4->AddEntry(hTotalBkg[i],"Total Bkg","l");
    if (bkgChoice != 0) leg4->AddEntry(hInvarMasspTBinRotBkg[i],"Combinatorial Bkg","l");

    if (drawResBkg) {
      if (bkgChoice != 0) {
        leg4->AddEntry(fPi0Bkg[i],"Residual Bkg. Fit","l");
      } else {
        leg4->AddEntry(fPi0Bkg[i],"Background Fit","l");
      }

    }

    leg4->AddEntry(fPi0Peak[i],"#pi^{0} Signal Fit","l");
    if (fEtaPeak[i]) leg4->AddEntry(fEtaPeak[i],"#eta Signal Fit","l");


    leg4->SetTextColor(kBlack);
    leg4->SetTextSize(0.035);
    leg4->SetBorderSize(0);
    leg4->SetFillColorAlpha(10, 0);
    leg4->Draw("same");

    // Making the individual plots
    cMassFit->cd();
    hInvarMasspTBin[i+iFirstRealBin]->Draw("E HIST");

    // (New) Drawing fitbkg subtracted part

    int   xMaxBin = hInvarMasspTBin[i+iFirstRealBin]->GetMaximumBin();
    float yMin = 0;
    float yMax = hInvarMasspTBin[i+iFirstRealBin]->GetBinContent(xMaxBin) + hInvarMasspTBin[i+iFirstRealBin]->GetBinErrorUp(xMaxBin);

    // FIXME estimate a good yMin to show background subtracted part
    yMin = -yMax * kMagicNegativeScale;

    hInvarMasspTBin[i+iFirstRealBin]->GetYaxis()->SetRangeUser(yMin,yMax*kMagicScale);

    if (haveRotBkg && bkgChoice != 0) {
      if (drawBkg) hInvarMasspTBinRotBkg[i+iFirstRealBin]->Draw("HIST SAME");
    }


    if (drawMCInfo && haveMCStatus) {
      for (Int_t k = 0; k < nMCId; k++) {
        hInvarMassPtBinMCId[i+iFirstRealBin][k]->Draw("E SAME");
      }
      hInvarMassPtBinMCNoPeak[i+iFirstRealBin]->Draw("E SAME");
      hInvarMassPtBinMCNoEta[i+iFirstRealBin]->Draw("E SAME");
      /*
      hInvarMasspTBinMCNoMatch[i]->Draw("E SAME");
      hInvarMasspTBinMCPi0Match[i]->Draw("E SAME");
      hInvarMasspTBinMCEtaMatch[i]->Draw("E SAME");
      hInvarMasspTBinMCGammaPCMatch[i]->Draw("E SAME");
      hInvarMasspTBinMCSinglePartMatch[i]->Draw("E SAME");
      hInvarMasspTBinMCSharedEtaAncMatch[i]->Draw("E SAME");
      hInvarMasspTBinMCSharedAncMatch[i]->Draw("E SAME");*/
    }

    if (drawFits) {
      fPi0Peak[i]->Draw("SAME");

      if (bDrawUnmod && fUnmodPi0Fit.size() > i) {
        fUnmodPi0Fit[i]->Draw("SAME");
      }
      if (fitBkgSub) hTotalBkg[i]->Draw("HIST SAME");
      hTotalFit[i]->Draw("HIST SAME");

      // FIXME turn this off after devel done
      if (drawResBkg) {
        fPi0Bkg[i]->Draw("SAME");
      }
    }



    // draw mass window lines here
    TLine * lMassLineLow = 0;
    TLine * lMassLineHigh = 0;
    if (bDrawMassWindowLines) {
      double fMassCutLow = Pi0MassCutLow[i];
      double fMassCutHigh = Pi0MassCutHigh[i];
      lMassLineLow = new TLine(fMassCutLow,yMin,fMassCutLow,yMax*kMagicScale);
      lMassLineHigh = new TLine(fMassCutHigh,yMin,fMassCutHigh,yMax*kMagicScale);

      lMassLineLow->SetLineStyle(kMassWindowLineStyle);
      lMassLineLow->SetLineColor(kMassWindowLineColor);
      lMassLineLow->SetLineWidth(kMassWindowLineWidth);
      lMassLineHigh->SetLineStyle(kMassWindowLineStyle);
      lMassLineHigh->SetLineColor(kMassWindowLineColor);
      lMassLineHigh->SetLineWidth(kMassWindowLineWidth);

      lMassLineLow->Draw("SAME");
      lMassLineHigh->Draw("SAME");
      leg4->AddEntry(lMassLineHigh,"#pi^{0} mass window","l");
    }

    leg4->Draw("same");
    if (fEtaPeak[i]) {
      fEtaPeak[i]->Draw("SAME");
    }

    DrawWIP(hInvarMasspTBin[i],0.25,0.75,0.23,0.15);
    //DrawAlicePerf(hInvarMasspTBin[i],0.25,0.75,0.23,0.13);
  //  DrawAlicePerf(hInvarMasspTBin[i],0.28,0.78,0.67,0.13);
//    DrawAlicePerf(hInvarMasspTBin[i],0.12,0.78,0.67,0.13);

    hInvarMasspTBin[i+iFirstRealBin]->Draw("E SAME");


    gPad->SetTickx();
    gPad->SetTicky();

    cMassFit->SetTickx();
    cMassFit->SetTicky();

    cMassFit->Print(Form("%s/MassFit_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cMassFit->Print(Form("%s/MassFit_%.0f_%.0f.png",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cMassFit->Print(Form("%s/MassFit_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));

    // Copying to the 2x2 histogram
    if ( i >= 2 && i < 6) {
      printf("Attempting to make 2x2\n");
      c2x2->cd(i-1); // 1, 2, 3, 4
      cMassFit->DrawClonePad();
    }
    // Copying to the 2x1, 1x2 histograms
    if (i == 3 || i == 4) {
      c2x1->cd(i-2);
      cMassFit->DrawClonePad();
      c1x2->cd(i-2);
      cMassFit->DrawClonePad();
    }

    // Making the individual MCInfo Plots
    if (haveMCStatus) {
      cMassFitMCInfo->cd();

      hInvarMasspTBin[i+iFirstRealBin]->Draw("E");
      for (Int_t k = 0; k < nMCId; k++) {
        hInvarMassPtBinMCId[i+iFirstRealBin][k]->Draw("E SAME");
      }
     // hInvarMassPtBinMCNoPeak[i+iFirstRealBin]->Draw("E SAME");
     // hInvarMassPtBinMCNoEta[i+iFirstRealBin]->Draw("E SAME");
      TLegend * legMCInfo = new TLegend(0.71,0.63,0.99,0.96); //..Bkg subtracted

      legMCInfo->AddEntry(hInvarMasspTBin[i],Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]),"pe");
      for (Int_t k = 0; k < nMCId; k++) {
        legMCInfo->AddEntry(hInvarMassPtBinMCId[i+iFirstRealBin][k],sMCIdTitles[k].Data(),"pe");
      }
      legMCInfo->AddEntry(hInvarMassPtBinMCNoPeak[i+iFirstRealBin],"MC Bkg","pe");
      gPad->SetLogy(1);

      legMCInfo->Draw("SAME");

      hInvarMasspTBin[i+iFirstRealBin]->GetYaxis()->UnZoom();

      cMassFitMCInfo->Print(Form("%s/MassFitMC_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
      cMassFitMCInfo->Print(Form("%s/MassFitMC_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1]+iFirstRealBin));
    }
  }

  cAllMass->Print(Form("%s/AllMass.pdf",sOutputDir.Data()));
  cAllMass->Print(Form("%s/AllMass.C",sOutputDir.Data()));

  c3->Print(Form("%s/MassFits.pdf",sOutputDir.Data()));
  c3->Print(Form("%s/MassFits.C",sOutputDir.Data()));

  c2x2->Print(Form("%s/MassFits_2_by_2.pdf",sOutputDir.Data()));
  c2x2->Print(Form("%s/MassFits_2_by_2.C",sOutputDir.Data()));

  c2x1->Print(Form("%s/MassFits_2_by_1.pdf",sOutputDir.Data()));
  c2x1->Print(Form("%s/MassFits_2_by_1.C",sOutputDir.Data()));

  c1x2->Print(Form("%s/MassFits_1_by_2.pdf",sOutputDir.Data()));
  c1x2->Print(Form("%s/MassFits_1_by_2.C",sOutputDir.Data()));


}

void PionID::DrawResultGraphs() {
  cout<<"Drawing Graphs"<<endl;

  int cNY = TMath::FloorNint(TMath::Sqrt(nPtBins-nSkipPoints));
  int cNX = TMath::CeilNint(((float) nPtBins-nSkipPoints)/cNY);
  if (MassBkgArray.size()) {
    TCanvas * cBkgInt = new TCanvas("cBkgInt");
    cBkgInt->Divide(cNX,cNY,0,0);

    TCanvas * cSpecificBin = new TCanvas("cSpecificBin");
    cSpecificBin->cd();
  for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      cBkgInt->cd(i+1);
      if (MassBkgArray[i]) {
        MassBkgArray[i]->Draw();
      } else { 
        printf("Missing bkgInt %d\n",i);
      }
    }
  }

  // Adding to nSkip Points here will skip the last points
  //nSkipPoints++;

  Pi0Yield = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0YieldArr[0],0,&Pi0YieldArrUn[0]);
  Pi0Yield->SetName("Pi0Yield");
  Pi0Yield->SetTitle("#pi_{0} Yield (Not Normalized)");
  Pi0Yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0Yield->SetMarkerStyle(kOpenSquare);
  Pi0Spectrum = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0NormYieldArr[0],0,&Pi0NormYieldArrUn[0]);
  Pi0Spectrum->SetName("Pi0Spectrum");
  Pi0Spectrum->SetTitle("#pi_{0} Spectrum");
  Pi0Spectrum->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0Spectrum->GetYaxis()->SetTitle("dN_{#pi}/d#it{p}_{T} (GeV/#it{c})^{-1}");
  Pi0Spectrum->SetMarkerStyle(kOpenSquare);

  Pi0IntYield = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0IntegralArr[0],0,&Pi0IntegralArrUn[0]);
  Pi0IntYield->SetName("Pi0IntYield");
  Pi0IntYield->SetTitle("#pi_{0}^{Cand.} Raw Yield (Signal + Background)");
  Pi0IntYield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0IntYield->SetMarkerStyle(kFullSquare);

  Pi0IntTotal = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0IntegralTotalArr[0],0,&Pi0IntegralTotalArrUn[0]);
  Pi0IntTotal->SetName("Pi0IntTotal");
  Pi0IntTotal->SetTitle("#pi_{0} Raw Yield (Not Normalized)");
  Pi0IntTotal->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0IntTotal->SetMarkerStyle(kFullSquare);

  Pi0IntSpectrum = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0NormIntegralArr[0],0,&Pi0NormIntegralArrUn[0]);
  Pi0IntSpectrum->SetName("Pi0IntSpectrum");
  Pi0IntSpectrum->SetTitle("#pi_{0} Spectrum (Integral)");
  Pi0IntSpectrum->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0IntSpectrum->GetYaxis()->SetTitle("dN_{#pi}/d#it{p}_{T} (GeV/#it{c})^{-1}");
  Pi0IntSpectrum->SetMarkerStyle(kFullSquare);

  if (haveMCStatus) {
    Pi0MCIntYield = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0MCIntegralArr[0],0,&Pi0MCIntegralArrUn[0]);
    Pi0MCIntYield->SetName("Pi0MCIntYield");
    Pi0MCIntYield->SetTitle("#pi_{0} (MC) Raw Yield (Not Normalized)");
    Pi0MCIntYield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    Pi0MCIntYield->SetMarkerStyle(kFullSquare);
    Pi0MCIntSpectrum = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0MCNormIntegralArr[0],0,&Pi0MCNormIntegralArrUn[0]);
    Pi0MCIntSpectrum->SetName("Pi0MCIntSpectrum");
    Pi0MCIntSpectrum->SetTitle("#pi_{0} Spectrum (MC, Integral)");
    Pi0MCIntSpectrum->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    Pi0MCIntSpectrum->GetYaxis()->SetTitle("dN_{#pi}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    Pi0MCIntSpectrum->SetMarkerStyle(kFullSquare);
  }


  Pi0Bkg = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0BkgArr[0],0,&Pi0BkgArrUn[0]);
  Pi0Bkg->SetName("Pi0Bkg");
  Pi0Bkg->SetTitle("#pi_{0} Background Estimate");
  Pi0Bkg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0Bkg->SetMarkerStyle(kFullSquare);
  Pi0Bkg->SetLineColor(kGreen);

  Pi0PeakSigRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0PeakSigRatioArr[0],0,&Pi0PeakSigRatioArrUn[0]);
  Pi0PeakSigRatio->SetName("Pi0PeakSigRatio");
  Pi0PeakSigRatio->SetTitle("#pi_{0} Peak Significance");
  Pi0PeakSigRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0PeakSigRatio->GetYaxis()->SetTitle("Yield/#sqrt{Background}");
  Pi0PeakSigRatio->SetMarkerStyle(kFullTriangleUp);

  Pi0YieldBkgRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0YieldBkgRatioArr[0],0,&Pi0YieldBkgRatioArrUn[0]);
  Pi0YieldBkgRatio->SetName("Pi0YieldBkgRatio");
  Pi0YieldBkgRatio->SetTitle("#pi_{0} Signal-to-Background Ratio");
  Pi0YieldBkgRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0YieldBkgRatio->GetYaxis()->SetTitle("Yield/Background");
  Pi0YieldBkgRatio->SetMarkerStyle(kFullTriangleUp);

  Pi0YieldTotalRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0YieldTotalRatioArr[0],0,&Pi0YieldTotalRatioArrUn[0]);
  Pi0YieldTotalRatio->SetName("Pi0YieldTotalRatio");
  Pi0YieldTotalRatio->SetTitle("#pi_{0} Signal-to-Total Ratio");
  Pi0YieldTotalRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0YieldTotalRatio->GetYaxis()->SetTitle("Yield/Total");
  Pi0YieldTotalRatio->SetMarkerStyle(kFullSquare);
  Pi0YieldTotalRatio->SetMarkerSize(kYieldTotalRatioMarkerSize);

//  TGraphErrors * MCBkg = 0;
//  TGraphErrors * MCPi0PeakSigRatio = 0;
 // TGraphErrors * MCPi0YieldBkgRatio = 0;
//  TGraphErrors * MCPi0YieldTotalRatio = 0;
  if (haveMCStatus) {
    MCBkg = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&MCBkgArr[0],0,&MCBkgArrUn[0]);
    MCBkg->SetName("MCBkg");
    MCBkg->SetTitle("#pi_{0} Background Estimate (MC)");
    MCBkg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    MCBkg->GetYaxis()->SetTitle("Background");
    MCBkg->SetMarkerStyle(kFullSquare);
    MCBkg->SetLineColor(kGreen);

    MCPi0PeakSigRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&MCPi0PeakSigRatioArr[0],0,&MCPi0PeakSigRatioArrUn[0]);
    MCPi0PeakSigRatio->SetName("MCPi0PeakSigRatio");
    MCPi0PeakSigRatio->SetTitle("#pi_{0} Peak Significance (MC)");
    MCPi0PeakSigRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    MCPi0PeakSigRatio->GetYaxis()->SetTitle("Yield/#sqrt{Background}");
    MCPi0PeakSigRatio->SetMarkerStyle(kFullTriangleUp);

    MCPi0YieldBkgRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&MCPi0YieldBkgRatioArr[0],0,&MCPi0YieldBkgRatioArrUn[0]);
    MCPi0YieldBkgRatio->SetName("MCPi0YieldBkgRatio");
    MCPi0YieldBkgRatio->SetTitle("#pi_{0} Signal-to-Background Ratio (MC)");
    MCPi0YieldBkgRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    MCPi0YieldBkgRatio->GetYaxis()->SetTitle("Yield/Background");
    MCPi0YieldBkgRatio->SetMarkerStyle(kFullTriangleUp);

    MCPi0YieldTotalRatio = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&MCPi0YieldTotalRatioArr[0],0,&MCPi0YieldTotalRatioArrUn[0]);
    MCPi0YieldTotalRatio->SetName("MCPi0YieldTotalRatio");
    MCPi0YieldTotalRatio->SetTitle("#pi_{0} Signal-to-Total Ratio (MC)");
    MCPi0YieldTotalRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    MCPi0YieldTotalRatio->GetYaxis()->SetTitle("Yield/Total");
    MCPi0YieldTotalRatio->SetLineColor(kMCYieldTotalColor);
    MCPi0YieldTotalRatio->SetMarkerColor(kMCYieldTotalColor);
    MCPi0YieldTotalRatio->SetMarkerStyle(kFullSquare);
    MCPi0YieldTotalRatio->SetMarkerSize(kYieldTotalRatioMarkerSize);

    RecMCPi0YieldRatio = new TGraphErrors(nPtBins - nSkipPoints,&ptPointsForTGraph[0],&RecMCPi0YieldRatioArr[0],0,&RecMCPi0YieldRatioArrUn[0]);
    RecMCPi0YieldRatio->SetName("RecMCPi0YieldRatioArr");
    RecMCPi0YieldRatio->SetTitle("#pi_{0} Reconstructed-to-MCTruth Yield Ratio");
    RecMCPi0YieldRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    RecMCPi0YieldRatio->GetYaxis()->SetTitle("Reconstructed Yield / MC Truth Yield");
    RecMCPi0YieldRatio->SetLineColor(kRecMCYieldRatioColor);
    RecMCPi0YieldRatio->SetMarkerColor(kRecMCYieldRatioColor);
    RecMCPi0YieldRatio->SetMarkerStyle(kRecMCYieldRatioMarkerStyle);

  }


  Pi0Mass = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0MassArr[0],0,&Pi0MassArrUn[0]);
  Pi0Mass->SetName("Pi0Mass");
  Pi0Mass->SetTitle("#pi_{0} Mass");
  Pi0Mass->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0Mass->GetYaxis()->SetTitle("M_{#pi_{0}} (GeV/#it{c}^2)");
  Pi0Mass->SetMarkerStyle(kFullSquare);
  Pi0Sigma = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0SigmaArr[0],0,&Pi0SigmaArrUn[0]);
  Pi0Sigma->SetName("Pi0Sigma");
  Pi0Sigma->SetTitle("#pi_{0} Sigma");
  Pi0Sigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0Sigma->GetYaxis()->SetTitle("#sigma (GeV/#it{c}^2)");
  Pi0Sigma->SetMarkerStyle(kFullCircle);

  Pi0ChiSquare = new TGraphErrors(nPtBins - nSkipPoints,&ptPointsForTGraph[0],&Pi0ChiSquareArr[0],0,0);
  Pi0ChiSquare->SetName("Pi0ChiSquare");
  Pi0ChiSquare->SetTitle("#pi_{0} Fit ChiSquare over NDF");
  Pi0ChiSquare->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  Pi0ChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF");
  Pi0ChiSquare->SetMarkerStyle(kFullSquare);

  if (haveMCStatus) {
    MCPi0ChiSquare = new TGraphErrors(nPtBins - nSkipPoints,&ptPointsForTGraph[0],&MCPi0ChiSquareArr[0],0,0);
    MCPi0ChiSquare->SetName("MCPi0ChiSquare");
    MCPi0ChiSquare->SetTitle("MC #pi_{0} Fit ChiSquare over NDF");
    MCPi0ChiSquare->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    MCPi0ChiSquare->GetYaxis()->SetTitle("#chi^{2}/NDF");
    MCPi0ChiSquare->SetMarkerStyle(kFullSquare);
  }

  // Fit Pi0 Mass, Sigma
 // TF1 * pi0MassFit = fitPi0Mass(Pi0Mass,"Pi0MassFit");
//  pi0MassFit = fitPi0Mass_KinkedLine(Pi0Mass,"Pi0MassFit");
  pi0MassFit = fit_KinkedLine(Pi0Mass,"Pi0MassFit",0);
  pi0MassFit->SetLineColor(kViolet);

//  TF1 * pi0SigmaFit = fitPi0Sigma(Pi0Sigma,"Pi0SigmaFit");
//  pi0SigmaFit = fitPi0Sigma_KinkedLine(Pi0Sigma,"Pi0SigmaFit");
  pi0SigmaFit = fit_KinkedLine(Pi0Sigma,"Pi0SigmaFit",1);
  pi0SigmaFit->SetLineColor(kBlue);



  TCanvas * cPi0Yield = new TCanvas("cPi0Yield");
//  Pi0Yield->Draw("ALP");
  Pi0IntSpectrum->Draw("ALP");
 // Pi0Spectrum->Draw("LP SAME");
//  TLegend * lPi0Spectrum = new TLegend(0.6,0.6,0.8,0.8);
//  lPi0Spectrum->AddEntry(Pi0IntSpectrum,"Integral","lp");
//  lPi0Spectrum->AddEntry(Pi0Spectrum,"Parameter","lp");
  cPi0Yield->SetLogy();
//  lPi0Spectrum->Draw("SAME");
  cPi0Yield->Print(Form("%s/Pi0Spectra.pdf",sOutputDir.Data()));

  TCanvas * cPi0Mass = new TCanvas("cPi0Mass");
  Pi0Mass->Draw("ALP");
  pi0MassFit->Draw("SAME");
  cPi0Mass->Print(Form("%s/Pi0Mass.pdf",sOutputDir.Data()));
  cPi0Mass->Print(Form("%s/Pi0Mass.C",sOutputDir.Data()));
  TCanvas * cPi0Sigma = new TCanvas("cPi0Sigma");
  Pi0Sigma->Draw("ALP");
  pi0SigmaFit->Draw("SAME");
  cPi0Sigma->Print(Form("%s/Pi0Sigma.pdf",sOutputDir.Data()));
  cPi0Sigma->Print(Form("%s/Pi0Sigma.C",sOutputDir.Data()));
  // Plot with both
  TCanvas * cPi0MassSigma = new TCanvas("cPi0MassSigma");
  cPi0MassSigma->Divide(1,2,0.0,0.0);
  cPi0MassSigma->cd(1);
  Pi0Mass->Draw("ALP");
  pi0MassFit->Draw("SAME");
  cPi0MassSigma->cd(2);
  Pi0Sigma->Draw("ALP");
  pi0SigmaFit->Draw("SAME");
  cPi0MassSigma->Print(Form("%s/Pi0MassSigma.pdf",sOutputDir.Data()));
  cPi0MassSigma->Print(Form("%s/Pi0MassSigma.C",sOutputDir.Data()));


  TCanvas * cPi0ChiSquare = new TCanvas("cChiSquare");
  Pi0ChiSquare->Draw("ALP");
  cPi0ChiSquare->Print(Form("%s/Pi0ChiSquare.pdf",sOutputDir.Data()));
  cPi0ChiSquare->Print(Form("%s/Pi0ChiSquare.C",sOutputDir.Data()));

  if (haveMCStatus) {
    MCPi0ChiSquare->Draw("ALP");
    cPi0ChiSquare->Print(Form("%s/MCPi0ChiSquare.pdf",sOutputDir.Data()));
    cPi0ChiSquare->Print(Form("%s/MCPi0ChiSquare.C",sOutputDir.Data()));
  }

  TCanvas * cPeakSigRatio = new TCanvas("cPeakSigRatio");
  Pi0PeakSigRatio->Draw("ALP");
  cPeakSigRatio->Print(Form("%s/Pi0PeakSig.pdf",sOutputDir.Data()));
  cPeakSigRatio->Print(Form("%s/Pi0PeakSig.C",sOutputDir.Data()));

  TCanvas * cYieldBkgRatio = new TCanvas("cYieldBkgRatio");
  Pi0YieldBkgRatio->Draw("ALP");
  cYieldBkgRatio->Print(Form("%s/Pi0YieldBkgRatio.pdf",sOutputDir.Data()));
  cYieldBkgRatio->Print(Form("%s/Pi0YieldBkgRatio.C",sOutputDir.Data()));

  TCanvas * cYieldTotalRatio = new TCanvas("cYieldTotalRatio","cYieldTotalRatio",kCanvasWidthWide,kCanvasHeightWide);
  Pi0YieldTotalRatio->Draw("ALP");
  TLegend * leg = new TLegend(0.20,0.20,0.55,0.45);
  //TLegend * leg = new TLegend(0.40,0.70,0.85,0.85);
  if (haveMCStatus) {
    Pi0YieldTotalRatio->GetYaxis()->SetRangeUser(0.,1.);
    MCPi0YieldTotalRatio->Draw("SAME LP");
    leg->AddEntry(Pi0YieldTotalRatio,"Analysis Result","lp");
    leg->AddEntry(MCPi0YieldTotalRatio,"MC Truth","lp");
    leg->Draw("SAME");
  }
  // FIXME add an indicator of the sigma range
  //TLegend * lSigmaIntegration = new TLegend(0.3,0.15,0.5,0.4);
  //lSigmaIntegration->SetHeader(Form("Yield/Total in [m_{#pi^{0}} - %.1f#sigma,m_{#pi^{0}} + %.1f#sigma]",nSigma,nSigma));
  //lSigmaIntegration->Draw("SAME");
  cYieldTotalRatio->Print(Form("%s/Pi0YieldTotalRatio.pdf",sOutputDir.Data()));
  cYieldTotalRatio->Print(Form("%s/Pi0YieldTotalRatio.C",sOutputDir.Data()));

  if (haveRotBkg) {
    TCanvas * cRotBkgSub = new TCanvas("cRotBkgSub");
    //cRotBkgSub->Divide(cNX,cNY,0.005,0.005);

    for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      cRotBkgSub->Clear();
      //cRotBkgSub->cd(i+1);
      hInvarMassPtBinRotSub[i+iFirstRealBin]->SetLineColor(kBlack);
      hInvarMassPtBinRotSub[i+iFirstRealBin]->SetMarkerStyle(0);
      hInvarMassPtBinRotSub[i+iFirstRealBin]->Draw();

      if (drawFits) {
        fPi0Bkg[i]->Draw("SAME");
        fPi0Fit[i]->Draw("SAME");

        if (fEtaPeak[i])  fEtaPeak[i]->Draw("SAME");
        if (haveMCStatus) hInvarMassPtBinMCId[i][2]->Draw("SAME");
      }
      hInvarMassPtBinRotSub[i+iFirstRealBin]->Draw("SAME");
      
      cRotBkgSub->Print(Form("%s/MassFit_BkgSub_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
      cRotBkgSub->Print(Form("%s/MassFit_BkgSub_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    }
    //cRotBkgSub->Print(Form("%s/MassFits_BkgSub.pdf",sOutputDir.Data()));
    //cRotBkgSub->Print(Form("%s/MassFits_BkgSub.C",sOutputDir.Data()));
  }

  // note to self: for generalization, need to be flexible if one or more does not exist
  if (haveRotBkg) {
    TCanvas * cBkgCmp = new TCanvas("cBkgCmp");
    cBkgCmp->Divide(cNX,cNY,0.005,0.005);
    for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      cBkgCmp->cd(i+1);
      hInvarMasspTBin[i]->Draw();
//      hInvarMasspTBinRotBkg[i]->Draw();
      if (haveMCStatus) {
        if (hInvarMassPtBinMCNoPeak[i]) hInvarMassPtBinMCNoPeak[i]->Draw("SAME");
        if (hInvarMassPtBinMCId[i][0]) hInvarMassPtBinMCId[i][0]->Draw("SAME");
        // Old Method (merging same energy, same position information)
        //if (hInvarMasspTBinRotBkgMCPi0[i]) hInvarMasspTBinRotBkgMCPi0[i]->Draw("SAME");
        //if (hInvarMasspTBinRotBkgMCEta[i]) hInvarMasspTBinRotBkgMCEta[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCEtaEnergyPair[i]) hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCEtaPosPair[i]) hInvarMasspTBinRotBkgMCEtaPosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCPi0EnergyPair[i]) hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCPi0PosPair[i]) hInvarMasspTBinRotBkgMCPi0PosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]) hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMinusMCAll[i]) hInvarMasspTBinRotBkgMinusMCAll[i]->Draw("SAME");
      }
      hInvarMasspTBinRotBkg[i]->Draw("SAME");
    }
    TLegend * lCmpBkg = new TLegend(0.1,0.1,0.9,0.9);
    lCmpBkg->AddEntry(hInvarMasspTBin[0],"Full Spectrum","lp");
    lCmpBkg->AddEntry(hInvarMasspTBinRotBkg[0],"Calc. Bkg.","lp");
    lCmpBkg->AddEntry(hTotalBkg[0],"Background Fit","lp");
    if (haveMCStatus) {
      if (hInvarMassPtBinMCNoPeak[0]) lCmpBkg->AddEntry(hInvarMassPtBinMCNoPeak[0],"MC w/o #pi^{0} peak","lp");
      if (hInvarMassPtBinMCId[0][0]) lCmpBkg->AddEntry(hInvarMassPtBinMCId[0][0],"MC unmatched bkg (true comb.)","lp");
      //if (hInvarMasspTBinRotBkgMCPi0[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCPi0[0],"Pos Swapped #pi^{0}","lp");
      //if (hInvarMasspTBinRotBkgMCEta[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCEta[0],"Pos Swapped #eta","lp");

//      if (hInvarMasspTBinRotBkgMCPi0EnergyPair[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCPi0EnergyPair[0],"PS #pi^{0} (Same Energy)","lp");
//      if (hInvarMasspTBinRotBkgMCPi0PosPair[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCPi0PosPair[0],"PS #pi^{0} (Same Position)","lp");
//      if (hInvarMasspTBinRotBkgMCEtaEnergyPair[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCEtaEnergyPair[0],"PS #eta (Same Energy)","lp");
//      if (hInvarMasspTBinRotBkgMCEtaPosPair[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMCEtaPosPair[0],"PS #eta (Same Position)","lp");
//      if (hInvarMasspTBinRotBkgMinusMCPi0PosPair[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMinusMCPi0PosPair[0],"PosSwap Minus PS #pi^{0} (Same Positions)","lp");
//      if (hInvarMasspTBinRotBkgMinusMCAll[0]) lCmpBkg->AddEntry(hInvarMasspTBinRotBkgMinusMCAll[0],"PosSwap Minus PS #pi^{0} and #eta","lp");
    }
    cBkgCmp->cd(11);
    lCmpBkg->Draw();
    cBkgCmp->Print(Form("%s/BkgCmp.pdf",sOutputDir.Data()));
    cBkgCmp->Print(Form("%s/BkgCmp.C"  ,sOutputDir.Data()));

    // Now Individual Plots 
    lCmpBkg->SetY2(0.85);
    lCmpBkg->SetX2(0.95);
    lCmpBkg->SetY1(0.50);
    lCmpBkg->SetX1(0.65);
    cBkgCmp->Clear();
    cBkgCmp->Divide(1,1);
    for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      hInvarMasspTBin[i]->SetTitle(TString::Format("%.1f #leq #it{p}_{T} < %.1f GeV/#it{c}",Pi0PtBins[i],Pi0PtBins[i+1]));
      hInvarMasspTBin[i]->Draw();
      if (haveMCStatus) {
        if (hInvarMassPtBinMCNoPeak[i]) hInvarMassPtBinMCNoPeak[i]->Draw("SAME");
        if (hInvarMassPtBinMCId[i][0]) hInvarMassPtBinMCId[i][0]->Draw("SAME");
        //if (hInvarMasspTBinRotBkgMCPi0[i]) hInvarMasspTBinRotBkgMCPi0[i]->Draw("SAME");
        //if (hInvarMasspTBinRotBkgMCEta[i]) hInvarMasspTBinRotBkgMCEta[i]->Draw("SAME");
        
//        if (hInvarMasspTBinRotBkgMCEtaEnergyPair[i]) hInvarMasspTBinRotBkgMCEtaEnergyPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCEtaPosPair[i]) hInvarMasspTBinRotBkgMCEtaPosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCPi0EnergyPair[i]) hInvarMasspTBinRotBkgMCPi0EnergyPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMCPi0PosPair[i]) hInvarMasspTBinRotBkgMCPi0PosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]) hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->Draw("SAME");
//        if (hInvarMasspTBinRotBkgMinusMCAll[i]) hInvarMasspTBinRotBkgMinusMCAll[i]->Draw("SAME");

        // FIXME plot background fit is present?
        hTotalBkg[i]->Draw("SAME");  
      }
      hInvarMasspTBinRotBkg[i]->Draw("SAME");
      lCmpBkg->Draw("SAME");
      //lCmpBkg->DrawBox(0.65,0.5,0.95,0.85);
      cBkgCmp->Print(Form("%s/BkgCmp_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
      cBkgCmp->Print(Form("%s/BkgCmp_%.0f_%.0f.png",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
      cBkgCmp->Print(Form("%s/BkgCmp_%.0f_%.0f.C"  ,sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    }
  }


/*
FitPeakMethod: 6
#             0 Gaussian
#             1 ExpDecay/Gaussian
#             2 Brent-Wigner
#             3 Crystal Ball
#             4 Crystal Ball (Right Side)
#             5 ExpDecay/Gaussian (Right Side)
#             6 ExpDecay/Gaussian (Both Sides)
#             7 Voigt Profile
*/

  // Draw MC Pi0 Peak Plots
  if (haveMCStatus) {
    TCanvas * cMCPi0 = new TCanvas("cMCPi0");
    TLegend * lMCPi0 = new TLegend(0.44,0.33,0.85,0.75);
    cMCPi0->Divide(cNX,cNY,0.005,0.005);
    for (int i = 0; i < nPtBins - nSkipPoints; i++) {
      cMCPi0->cd(i+1);
      hInvarMassPtBinMCId[i][2]->Draw();
      hInvarMassPtBinMCId[i][2]->GetXaxis()->SetRangeUser(0.,0.4);
      if (fMCPi0Fit.size() > i) {
        fMCPi0Fit[i]->Draw("SAME");
        fMCPi0Fit[i]->SetNpx(150);
        lMCPi0->AddEntry(fMCPi0Fit[i],sPeakNames[fitPeakMethod],"l");
      }
    }
    lMCPi0->Draw("SAME");
    cMCPi0->Print(Form("%s/MCPi0.pdf",sOutputDir.Data()));
    cMCPi0->Print(Form("%s/MCPi0.C"  ,sOutputDir.Data()));

    // Draw individual plots
    cMCPi0->Clear();
    for (int i = 0; i < nPtBins - nSkipPoints; i++) {
      lMCPi0->Clear();
      hInvarMassPtBinMCId[i][2]->Draw();
      hInvarMassPtBinMCId[i][2]->GetXaxis()->SetRangeUser(0.,0.4);
      hInvarMassPtBinMCId[i][2]->GetYaxis()->SetTitleOffset(1.5*hInvarMassPtBinMCId[i][2]->GetYaxis()->GetTitleOffset());
      if (fMCPi0Fit.size() > i) {
        fMCPi0Fit[i]->Draw("SAME");
        lMCPi0->AddEntry(fMCPi0Fit[i],sPeakNames[fitPeakMethod],"l");
        // Add the parameters
        for (int j = 0; j < fMCPi0Fit[i]->GetNpar(); j++ ) {
          lMCPi0->AddEntry((TObject*)0,Form("%s = %f #pm %f",fMCPi0Fit[i]->GetParName(j),fMCPi0Fit[i]->GetParameter(j),fMCPi0Fit[i]->GetParError(j)),"");
        }
      }
      lMCPi0->Draw("SAME");
      cMCPi0->Print(Form("%s/MCPi0_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
      cMCPi0->Print(Form("%s/MCPi0_%.0f_%.0f.C"  ,sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    }
  }

  // Ploting the Normal Background Subtraction results
  TCanvas * cBkgSub = new TCanvas("cBkgSub","cBkgSub",600,700);
  TLegend * legBkgSub = 0; 
  if (haveMCStatus) legBkgSub = new TLegend(0.44,0.33,0.60,0.5);
  else {
    if (iThetaModelCent == 0) legBkgSub = new TLegend(0.53,0.20,0.64,0.52);
    else if (iThetaModelCent == 1) legBkgSub = new TLegend(0.35,0.48,0.64,0.74);
    else if (iThetaModelCent == 2) legBkgSub = new TLegend(0.35,0.41,0.64,0.67);
    else legBkgSub = new TLegend(0.35,0.41,0.64,0.74);
  }
  legBkgSub->SetTextColor(kBlack);
  legBkgSub->SetTextSize(0.035);
  legBkgSub->SetBorderSize(0);
  legBkgSub->SetFillColorAlpha(10, 0);
  //cBkgSub->Divide(cNX,cNY,0.005,0.005);
  for (int i = 0; i < nPtBins - nSkipPoints; i++) {
    //cBkgSub->cd(i+1);
    cBkgSub->Clear();
    legBkgSub->Clear();
//    hInvarMassBkgSubpTBin[i]->SetLineColor(kBlack);
//    hInvarMassBkgSubpTBin[i]->SetMarkerColor(kBlack);
    
    hInvarMassBkgSubpTBin[i]->SetTitle("");
    hInvarMassBkgSubpTBin[i]->GetXaxis()->SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})");
    hInvarMassBkgSubpTBin[i]->SetLineColor(kBkgSubColor);
    hInvarMassBkgSubpTBin[i]->SetMarkerColor(kBkgSubColor);
    hInvarMassBkgSubpTBin[i]->SetMarkerStyle(kBkgSubStyle);
    hInvarMassBkgSubpTBin[i]->SetMarkerSize(0.7);
    hInvarMassBkgSubpTBin[i]->SetFillColorAlpha(0,0.7);

    hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetTitleOffset(0.83);
    hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetLabelSize(0.03);
    hInvarMasspTBin[i+iFirstRealBin]->GetXaxis()->SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})");

    hInvarMasspTBin[i+iFirstRealBin]->Draw("E HIST");

    hTotalBkg[i]->Draw("E SAME");
    hInvarMassBkgSubpTBin[i]->Draw("E SAME");
    

    legBkgSub->AddEntry((TObject*)0,Form("%0.0f < #it{p}_{T}^{#gamma#gamma} < %0.0f GeV/#it{c}",Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]),"");
    legBkgSub->AddEntry(hInvarMasspTBin[i+iFirstRealBin],"Raw invariant mass","pe");

    legBkgSub->AddEntry(hTotalBkg[i],"BG fit","pe");
    // Adding in background fit without angle correction


    legBkgSub->AddEntry(hInvarMassBkgSubpTBin[i],"BG subtracted","pe");

    if (drawFits) {
      fPi0Peak[i]->Draw("SAME");
      legBkgSub->AddEntry(fPi0Peak[i],"#pi^{0} Signal Fit","l");
    }
    if (fEtaPeak[i]) {
      fEtaPeak[i]->Draw("SAME");
      legBkgSub->AddEntry(fEtaPeak[i],"#eta Signal Fit","l");
    }
    if (haveMCStatus) {
      hInvarMassPtBinMCId[i][2]->Draw("SAME");
      legBkgSub->AddEntry(hInvarMassPtBinMCId[i][2],"True #pi^{0}#rightarrow#gamma#gamma Cluster Pairs","p");
    }
    //hInvarMassBkgSubpTBin[i]->Draw("SAME");
    int   xMaxBin = hInvarMasspTBin[i+iFirstRealBin]->GetMaximumBin();
    float yMin = 0;
    float yMax = hInvarMasspTBin[i+iFirstRealBin]->GetBinContent(xMaxBin) + hInvarMasspTBin[i+iFirstRealBin]->GetBinErrorUp(xMaxBin);
     
    // FIXME estimate a good yMin to show background subtracted part
    yMin = -yMax * kMagicNegativeScale;


    // Draw those lines
    TLine * lMassLineLow = 0;
    TLine * lMassLineHigh = 0;
    TLine * lSBMassLineLow = 0; // Fixed to SB23
    TLine * lSBMassLineHigh = 0;
    if (bDrawMassWindowLines) {
      double fMassCutLow = Pi0MassCutLow[i];
      double fMassCutHigh = Pi0MassCutHigh[i];
      double fSBMassCutLow = SBMassCut1[i];
      double fSBMassCutHigh = SBMassCut3[i];
      lMassLineLow = new TLine(fMassCutLow,yMin,fMassCutLow,yMax*kMagicScale);
      lMassLineHigh = new TLine(fMassCutHigh,yMin,fMassCutHigh,yMax*kMagicScale);

      lMassLineLow->SetLineStyle(kMassWindowLineStyle);
      lMassLineLow->SetLineColor(kMassWindowLineColor);
      lMassLineLow->SetLineWidth(kMassWindowLineWidth);
      lMassLineHigh->SetLineStyle(kMassWindowLineStyle);
      lMassLineHigh->SetLineColor(kMassWindowLineColor);
      lMassLineHigh->SetLineWidth(kMassWindowLineWidth);

      lMassLineLow->Draw("SAME");
      lMassLineHigh->Draw("SAME");
      legBkgSub->AddEntry(lMassLineHigh,"#pi^{0} mass window","l");

      lSBMassLineLow = new TLine(fSBMassCutLow,yMin,fSBMassCutLow,hInvarMasspTBin[i+iFirstRealBin]->GetBinContent(hInvarMasspTBin[i+iFirstRealBin]->FindFixBin(fSBMassCutLow)));
      lSBMassLineHigh = new TLine(fSBMassCutHigh,yMin,fSBMassCutHigh,hInvarMasspTBin[i+iFirstRealBin]->GetBinContent(hInvarMasspTBin[i+iFirstRealBin]->FindFixBin(fSBMassCutHigh)));
      //lSBMassLineLow = new TLine(fSBMassCutLow,yMin,fSBMassCutLow,yMax*kMagicScale);
      //lSBMassLineHigh = new TLine(fSBMassCutHigh,yMin,fSBMassCutHigh,yMax*kMagicScale);

      lSBMassLineLow->SetLineStyle(kSBMassWindowLineStyle);
      lSBMassLineLow->SetLineColor(kSBMassWindowLineColor);
      lSBMassLineLow->SetLineWidth(kSBMassWindowLineWidth);
      lSBMassLineHigh->SetLineStyle(kSBMassWindowLineStyle);
      lSBMassLineHigh->SetLineColor(kSBMassWindowLineColor);
      lSBMassLineHigh->SetLineWidth(kSBMassWindowLineWidth);

      lSBMassLineLow->Draw("SAME");
      lSBMassLineHigh->Draw("SAME");
      legBkgSub->AddEntry(lSBMassLineHigh,"Sideband mass window","l");
    }


    legBkgSub->Draw("SAME");
    if (bEnablePerformance) {
      if (iThetaModelCent == 0) DrawAlicePerf(hInvarMasspTBin[i],0.28,0.73,0.23,0.14);
      else if (iThetaModelCent == 1) DrawAlicePerf(hInvarMasspTBin[i],0.28,0.76,0.23,0.14);
      else if (iThetaModelCent == 2) DrawAlicePerf(hInvarMasspTBin[i],0.28,0.73,0.23,0.14);
      else DrawAlicePerf(hInvarMasspTBin[i],0.28,0.76,0.23,0.14);
    } else {
      DrawWIP(hInvarMasspTBin[i],0.25,0.67,0.23,0.18);
    }


    gPad->SetTickx();
    gPad->SetTicky();

    cBkgSub->Print(Form("%s/MassFit_FitSub_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cBkgSub->Print(Form("%s/MassFit_FitSub_%.0f_%.0f.png",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cBkgSub->Print(Form("%s/MassFit_FitSub_%.0f_%.0f.eps",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cBkgSub->Print(Form("%s/MassFit_FitSub_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
  }
  //cBkgSub->Print(Form("%s/MassFits_FitSub.pdf",sOutputDir.Data()));
  //cBkgSub->Print(Form("%s/MassFits_FitSub.C",sOutputDir.Data()));

  // Printing out mu, sigma tables
//    TGraphErrors * Pi0Mass = new TGraphErrors(nPtBins-nSkipPoints,&ptPointsForTGraph[0],&Pi0MassArr[0],0,&Pi0MassArrUn[0]);

}


void PionID::PrintResultsTables() {

  TCanvas * cResultsCanvas = new TCanvas("ResultsCanvas","ResultsCanvas",900,500);

  TPaveText *pt = new TPaveText(.05,.1,.95,.8);


  int kNoGammaBins = 9;

  if (sLabel.Length() > 0) pt->AddText(sLabel.Data());
  if (sLabel2.Length() > 0) pt->AddText(sLabel2.Data());

  pt->AddText("The entries after the first five are just filler");

  TString sPi0Mass = "Pi0Mass:   {";

  for (int i = 0; i < nPtBins - nSkipPoints; i++ ) {
    if (i == 0) sPi0Mass += Form(" %f",Pi0MassArr[i]);
    else sPi0Mass += Form(", %f", Pi0MassArr[i]);
  }
  for (int i = nPtBins - nSkipPoints; i < kNoGammaBins; i++ ){
    sPi0Mass += ", 0.15";
  }
  sPi0Mass += "}";

  pt->AddText(sPi0Mass.Data());

  TString sPi0Sigma = "Pi0Sigma: {";
  for (int i = 0; i < nPtBins - nSkipPoints; i++ ) {
    if (i == 0) sPi0Sigma += Form(" %f",Pi0SigmaArr[i]);
    else sPi0Sigma += Form(", %f", Pi0SigmaArr[i]);
  }
  for (int i = nPtBins - nSkipPoints; i < kNoGammaBins; i++ ){
    sPi0Sigma += ", 0.01";
  }
  sPi0Sigma += "}";

  pt->AddText(sPi0Sigma.Data());

  pt->Draw();

  cResultsCanvas->Print(Form("%s/Results.pdf",sOutputDir.Data()));

  printf("Pt Bins:         {");
  for (int i = 0; i < nPtBins; i++) {
    if (i) printf(", %.2f",Pi0PtBins[i]);
    else printf(" %.2f",Pi0PtBins[i]);
  }
  printf(" };\n");

  printf("Pi0 Mean Pt:     {");
  for (int i = 0; i < nPtBins - nSkipPoints; i++ ) {
    if (i == 0) printf(" %f",ptPointsForTGraph[i]);
    else printf(", %f", ptPointsForTGraph[i]);
  }
  printf(" };\n");

  printf("Pi0Mass:         {");
  for (int i = 0; i < nPtBins - nSkipPoints; i++ ) {
    if (i == 0) printf(" %f",Pi0MassArr[i]);
    else printf(", %f", Pi0MassArr[i]);
  }
  printf(" };\n");
  printf("Pi0Sigma:        {");
  for (int i = 0; i < nPtBins - nSkipPoints; i++ ) {
    if (i == 0) printf(" %f",Pi0SigmaArr[i]);
    else printf(", %f", Pi0SigmaArr[i]);
  }
  printf(" };\n");

  printf("Pi0 Peak Purity: {");
  for (int i = 0; i < nPtBins - nSkipPoints; i++) {
    if (i == 0) printf(" %f",Pi0YieldTotalRatioArr[i]);
    else printf(", %f", Pi0YieldTotalRatioArr[i]);
  }
  printf("};\n");

  // Printing out parametrizations
  printf("Pi0 Mass Fit Parameters: \n");
  PrintParameters(pi0MassFit,0,"Pi0MassFit");
  printf("Pi0 Sigma Fit Parameters: \n");
  PrintParameters(pi0SigmaFit,0,"Pi0SigmaFit");


}

/* For MC, make plots showing potential PSCorrection with Same-Pos Pi0 pairs
 */
void PionID::DrawPosSwapMCSub() {

  // Building the ratios
  for (int i = 0; i < nPtBins - nSkipPoints; i++)
  {
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio.push_back((TH1D *) hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->Clone(Form("%s_Ratio",hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]->GetName())));
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio[i]->Divide(hInvarMassPtBinMCNoPeak[i+iFirstRealBin]);
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio[i]->GetYaxis()->SetTitle("(Pos Swap - MC PosSwapped #pi^{0} peak (SP) ) / True MC Background");
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio[i]->SetTitle("(Pos Swap - MC PosSwapped #pi^{0} peak (SP) ) / True MC Background");
  }


  TCanvas * cPSC_MC = new TCanvas("cPSC_MC","cPSC_MC",fCanvasWidth,fCanvasHeight);
  cPSC_MC->Divide(1,2);

  // Build the Legend
  TLegend * lPSC_MC = new TLegend(0.55,0.15,0.95,0.46);
  lPSC_MC->AddEntry(hInvarMassPtBinMCNoPeak[0],"True MC Background","lp");
  lPSC_MC->AddEntry(hInvarMasspTBinRotBkgMinusMCPi0PosPair[0],"Pos Swap - MC PosSwapped #pi^{0} peak (SP)","lp");

  for (int i = 0; i < nPtBins-nSkipPoints; i++)
  {
    cPSC_MC->cd(1);
    hInvarMasspTBinRotBkgMinusMCPi0PosPair[i+iFirstRealBin]->Draw("E");
    hInvarMassPtBinMCNoPeak[i+iFirstRealBin]->Draw("E SAME");
    lPSC_MC->Draw("SAME");

    cPSC_MC->cd(2);
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio[i]->Draw();
    hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio[i]->GetYaxis()->SetRangeUser(0.5,1.5);

    cPSC_MC->Print(Form("%s/PosSwap_MCCorr_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    cPSC_MC->Print(Form("%s/PosSwap_MCCorr_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
  }

}

void PionID::DrawDemoPlots() {
  cout<<"Drawing the Demo plots."<<endl;
  TCanvas * c = new TCanvas("cDDP","cDDP",1200,900);
  TLegend * l = new TLegend(0.46,0.47,0.9,0.86);

  Int_t kPSColor = kAzure+10;
  Int_t kDemoColor = kYellow;

  Double_t arbitrary_scale = 0.3; //FIXME

  for (int i = 0; i < nPtBins-nSkipPoints; i++)
  {
    l->Clear();
    hInvarMasspTBin[i+iFirstRealBin]->Draw("E");
    l->AddEntry(hInvarMasspTBin[i+iFirstRealBin],"Same Event Pairs","lp");

    hInvarMassPtBinMCNoPeak[i+iFirstRealBin]->Draw("E SAME");
    l->AddEntry(hInvarMassPtBinMCNoPeak[i+iFirstRealBin],"True MC Background","lp");
    hInvarMasspTBinRotBkg[i+iFirstRealBin]->Draw("HIST SAME");
    l->AddEntry(hInvarMasspTBinRotBkg[i+iFirstRealBin],"PosSwap Background","lp");

    if (fPSFinalPeak!=0) {
      hInvarMasspTBinPSCorr[i]->SetLineColor(kPSColor);
      hInvarMasspTBinPSCorr[i]->SetMarkerColor(kPSColor);
      hInvarMasspTBinPSCorr[i]->Draw("SAME");
      l->AddEntry(hInvarMasspTBinPSCorr[i],"Convolved Signal Model","lp");
      TH1D * hLocalSub = (TH1D *) hInvarMasspTBinRotBkg[i]->Clone(Form("RotBkgPSSub_%d",i));
      hLocalSub->Rebin(2); // FIXME
      hLocalSub->Scale(1./2.);
      hLocalSub->SetLineColor(kDemoColor);
      hLocalSub->SetMarkerColor(kDemoColor);
      hLocalSub->Add(hInvarMasspTBinPSCorr[i],-1*arbitrary_scale);
      hLocalSub->Draw("SAME");
      hInvarMasspTRotBkgPSSub.push_back(hLocalSub);
      l->AddEntry(hInvarMasspTRotBkgPSSub[i],"PosSwap - Convolved Signal","lp");
    }
    // Draw the true MC Pi0 Swapped Piece
    hInvarMasspTBinRotBkgMCPi0PosPair[i]->SetLineColor(kBlue);
    hInvarMasspTBinRotBkgMCPi0PosPair[i]->Draw("SAME");
    l->AddEntry(hInvarMasspTBinRotBkgMCPi0PosPair[i],"True Same-Pos Pi0","lp");

    l->Draw("SAME");
    c->Print(Form("%s/DemoPlot_%.0f_%.0f.pdf",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
    c->Print(Form("%s/DemoPlot_%.0f_%.0f.C",sOutputDir.Data(),Pi0PtBins[i+iFirstRealBin],Pi0PtBins[i+1+iFirstRealBin]));
  }
}

void PionID::SaveResults() {
  cout<<"Saving Results"<<endl;

  printf("Attempting saving to output file now.\n");

  //TString outFilePath = "SecondAnalysis.root";
  TString outFilePath = Form("%s/%s",sOutputDir.Data(),sOutputFileName.Data());
  TFile * outFile = TFile::Open(outFilePath,"RECREATE");
  if (!outFile) return;

  printf("So far, so good\n");

  if (haveMCStatus) {
    printf("Evidently, I have MC information, so MC histograms will be stored\n");
  }
//  if (fClusEnergy) {
//  //  fClusEnergy->SetName("ClusEnergySpectrum");
//  //  outFile->Add(fClusEnergy);
//    TH1F * ClusEnergySpectrum = (TH1F *) fClusEnergy->Clone("ClusEnergySpectrum");
//    outFile->Add(ClusEnergySpectrum);
//  }
//  printf("Added a clone of fClusEnergy\n");
//  printf("Did not add a clone of fClusEnergy\n");

  if (hGammaE) {
    outFile->Add(hGammaE);
  }

  outFile->Add(Pi0Spectrum);
  outFile->Add(Pi0IntSpectrum);
  outFile->Add(Pi0Yield);
  outFile->Add(Pi0IntYield);
  outFile->Add(Pi0IntTotal);
  if (haveMCStatus) {
    outFile->Add(Pi0MCIntYield);
    outFile->Add(Pi0MCIntSpectrum);
  }
  outFile->Add(Pi0Mass);
  outFile->Add(pi0MassFit);
  outFile->Add(Pi0Sigma);
  outFile->Add(pi0SigmaFit);
  outFile->Add(Pi0ChiSquare);
  if (haveMCStatus) outFile->Add(MCPi0ChiSquare);
  outFile->Add(Pi0Bkg);
  outFile->Add(Pi0PeakSigRatio);
  outFile->Add(Pi0YieldBkgRatio);
  outFile->Add(Pi0YieldTotalRatio);
  if (haveMCStatus) {
    outFile->Add(MCBkg);
    outFile->Add(MCPi0PeakSigRatio);
    outFile->Add(MCPi0YieldBkgRatio);
    outFile->Add(MCPi0YieldTotalRatio);
    outFile->Add(RecMCPi0YieldRatio);

    // Save individual MC th2s? 
    for (Int_t i = 0; i < nMCId; i++) {
      if(fInvarMassPtMCId[i]) outFile->Add(fInvarMassPtMCId[i]);
    }

    if(fInvarMassPtMCNoPeak) outFile->Add(fInvarMassPtMCNoPeak);
    if(fInvarMassPtMCNoEta) outFile->Add(fInvarMassPtMCNoEta);
  }

  outFile->Add(hPairOpeningAngle);
  outFile->Add(hPairOpeningAngleRotBkg);
  outFile->Add(hPairLambdaBkgOpeningAngle);
  outFile->Add(hPairPtBkgOpeningAngle);
 

  for (int i = 0; i < nPtBins-nSkipPoints; i++) {
    outFile->Add(hInvarMasspTBin[i+iFirstRealBin]);
    if (hInvarMasspTBinResid[i+iFirstRealBin]) outFile->Add(hInvarMasspTBinResid[i+iFirstRealBin]);
    if (hInvarMasspTBinRotBkg[i+iFirstRealBin]) outFile->Add(hInvarMasspTBinRotBkg[i+iFirstRealBin]);
    if (hInvarMassPtBinRotSub[i+iFirstRealBin]) outFile->Add(hInvarMassPtBinRotSub[i+iFirstRealBin]);
    if (hInvarMassBkgSubpTBin[i]) outFile->Add(hInvarMassBkgSubpTBin[i]);

    if (haveMCStatus) {
      for (Int_t k = 0; k < nMCId; k++) {
        outFile->Add(hInvarMassPtBinMCId[i+iFirstRealBin][k]);
      }
      if (hInvarMassPtBinMCNoPeak[i+iFirstRealBin]) outFile->Add(hInvarMassPtBinMCNoPeak[i+iFirstRealBin]);
      if (hInvarMassPtBinMCNoEta[i+iFirstRealBin]) outFile->Add(hInvarMassPtBinMCNoEta[i+iFirstRealBin]);
    }

    if (fPi0Fit[i]) outFile->Add(fPi0Fit[i]);
    if (fPi0Peak[i]) outFile->Add(fPi0Peak[i]);
    if (fPi0Bkg[i]) outFile->Add(fPi0Bkg[i]);
    if (fEtaPeak[i]) outFile->Add(fEtaPeak[i]);
    if (haveMCStatus) {
      if (fMCPi0Fit[i]) outFile->Add(fMCPi0Fit[i]);
    }
  }
  if (hOpeningAngleCorrection) outFile->Add(hOpeningAngleCorrection);
  if (h2DOpeningAngleCorr) outFile->Add(h2DOpeningAngleCorr);
  if (hPairPtOpAngle) outFile->Add(hPairPtOpAngle);
  if (hPairPtOpAngleBkg) outFile->Add(hPairPtOpAngleBkg);

  if (fInvarMasspT) outFile->Add(fInvarMasspT);
  if (fInvarMasspTRotBkg) outFile->Add(fInvarMasspTRotBkg);
  if (fInvarMasspTRotBkgAngleScaled) outFile->Add(fInvarMasspTRotBkgAngleScaled);
  if (fInvarMassPtRotBkgGAPatch) outFile->Add(fInvarMassPtRotBkgGAPatch);
  if (fInvarMassPtRotBkgNoPatch) outFile->Add(fInvarMassPtRotBkgNoPatch);
  if (haveMCStatus) {
    if (fInvarMasspTRotBkgMCPi0) outFile->Add(fInvarMasspTRotBkgMCPi0);
    if (fInvarMasspTRotBkgMCEta) outFile->Add(fInvarMasspTRotBkgMCEta);
    if (fInvarMasspTRotBkgMCPi0EnergyPair) outFile->Add(fInvarMasspTRotBkgMCPi0EnergyPair);
    if (fInvarMasspTRotBkgMCEtaEnergyPair) outFile->Add(fInvarMasspTRotBkgMCEtaEnergyPair);
    if (fInvarMasspTRotBkgMCPi0PosPair) outFile->Add(fInvarMasspTRotBkgMCPi0PosPair);
    if (fInvarMasspTRotBkgMCEtaPosPair) outFile->Add(fInvarMasspTRotBkgMCEtaPosPair);
    if (fInvarMasspTRotBkgMinusMCPi0Pos) outFile->Add(fInvarMasspTRotBkgMinusMCPi0Pos);
    if (fInvarMasspTRotBkgMinusMCAll) outFile->Add(fInvarMasspTRotBkgMinusMCAll);

    for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      if (hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]) outFile->Add(hInvarMasspTBinRotBkgMinusMCPi0PosPair[i]);
    }
  }

  if (fPSEnergyPairDist) outFile->Add(fPSEnergyPairDist);
  if (fPSFinalPeak) outFile->Add(fPSFinalPeak);
  if (hInvarMasspTBinPSCorr.size() != 0) {
    for (int i = 0; i < nPtBins-nSkipPoints; i++) {
      if (hInvarMasspTBinPSCorr[i]) outFile->Add(hInvarMasspTBinPSCorr[i]);
    }
  }
 
  /* Save space by not copying these
  if (hPtEPAnglePionAcc)       outFile->Add(hPtEPAnglePionAcc); 
  if (hPtEPAngleMCPion)        outFile->Add(hPtEPAngleMCPion);
  if (hPtEPAngleTrueRecMCPion) outFile->Add(hPtEPAngleTrueRecMCPion);
  if (hHistTrackPsiEPPtCent)   outFile->Add(hHistTrackPsiEPPtCent);

  if (hPtRPAnglePionAcc)       outFile->Add(hPtRPAnglePionAcc);
  if (hPtRPAngleMCPion)        outFile->Add(hPtRPAngleMCPion);
  if (hPtRPAngleTrueRecMCPion) outFile->Add(hPtRPAngleTrueRecMCPion);
  if (hHistTrackPsiRPPtCent)   outFile->Add(hHistTrackPsiRPPtCent);
 */
  printf("Made it this far\n");
  if (hPtEPAnglePionAcc) {
    for (int i = 0; i <  hPtEPAnglePionAcc_Proj.size(); i++) {
      outFile->Add(hPtEPAnglePionAcc_Proj[i]);
      if (haveMCStatus) {
        if (hPtEPAngleMCPion) outFile->Add(hPtEPAngleMCPion_Proj[i]);
        if (hPtEPAngleTrueRecMCPion) outFile->Add(hPtEPAngleTrueRecMCPion_Proj[i]);
      }

      if (hPtRPAnglePionAcc) outFile->Add(hPtRPAnglePionAcc_Proj[i]);
      if (haveMCStatus) {
        if (hPtRPAngleMCPion) outFile->Add(hPtRPAngleMCPion_Proj[i]);
        if (hPtRPAngleTrueRecMCPion) outFile->Add(hPtRPAngleTrueRecMCPion_Proj[i]);
      }
    }

    if (hHistTrackPsiEPPt) outFile->Add(hHistTrackPsiEPPt);
    if (hHistTrackPsiRPPt) outFile->Add(hHistTrackPsiRPPt);

    if (gTrack_Bv) outFile->Add(gTrack_Bv);
    if (gTrack_V2) outFile->Add(gTrack_V2);
    if (gTrack_V3) outFile->Add(gTrack_V3);
    if (gTrack_V4) outFile->Add(gTrack_V4);
    if (gTrack_V6) outFile->Add(gTrack_V6);
    if (gTrigger_Bv) outFile->Add(gTrigger_Bv);
    if (gTrigger_V2) outFile->Add(gTrigger_V2);
    if (gTrigger_V3) outFile->Add(gTrigger_V3);
    if (gTrigger_V4) outFile->Add(gTrigger_V4);
    if (gTrigger_V6) outFile->Add(gTrigger_V6);

  }

  printf("Done.  Writing to output file %s\n",outFile->GetName());

  outFile->Write();
}


void PionID::Run() {
  cout<<"Beginning Pi0 Candidate Analysis"<<endl;
  
  string PLine     ="+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
  string EmptyLine ="+                                                               +";
  string TitleLine ="+                      Activating Phase 1                       +";

  cout<<PLine<<endl<<EmptyLine<<endl<<TitleLine<<endl<<EmptyLine<<endl<<PLine<<endl;
  

  if (fDebugLevel > 1) PrintSettings();

  if (LoadHistograms()) {
    fprintf(stderr,"There was an error loading the histograms. Exitting ...\n");
    return;
  }
  
  if (haveMCStatus) { // For MC, don't load preanalysis from another file
    // use preanalysis from the same file
    bUseMCPreAnalysis = false;
  }

  if (bUseMCPreAnalysis) {
    LoadMCPreAnalysis();
  }

  if (bUseThetaLookUpTable) {
    LoadThetaModelParameters();
  }

  // FIXME put a switch here
  OpeningAngleAnalysis();

  SetPtBins();
  MakeBasicPlots();

  InitializeWSignal();
  if (fUScaleMatrix && bEnablePSScaleMethod) {
    DrawPSScaleCorrectionPlots();
  }

  if (fPSMassPtMap && bEnablePSDirectMethod) {
    DrawPSDirectCorrectionPlots();
  }

  if (hPtEPAnglePionAcc && !haveMCStatus) { // don't do this for MC (a few bugs)
    ProduceDeltaPsiPlots();
    MeasureVn();
  }

  DoProjections();

  // Trying out doing the MCTruthPi0 before the main analysis
  // and (in the case of MC) storing the results of the Pi0 analysis
  // as if they were a preloaded MC fit
  if (haveMCStatus) {
    FitMCTruthPi0();
    bUseMCPreAnalysis = true; // could make another config variable to disable this
  }

  Pi0MassAnalysis();

  AnalyzeMatchedTracks();

//  if (haveMCStatus) {
//    FitMCTruthPi0();
//  }

  DrawMassPlots();
  DrawResultGraphs();
  PrintResultsTables();

  if (bkgType == 4 && !bUsingClusPairRot) DrawPosSwapMCSub();
  if (fPSMassPtMap && bEnablePSDirectMethod) {
    DrawDemoPlots();
  }

  SaveResults();

}

