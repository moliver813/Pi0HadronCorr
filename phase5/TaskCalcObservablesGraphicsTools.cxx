#ifndef TASKCALCOBSERVABLESGRAPHICSTOOLS_CXX
#define TASKCALCOBSERVABLESGRAPHICSTOOLS_CXX

// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
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

void FindCommonMinMax(std::vector<TH1D *> fHists, Double_t * fMin, Double_t * fMax) {
  Int_t nHists = fHists.size();

  Double_t fFinalMin = fHists[0]->GetBinContent(fHists[0]->GetMinimumBin());
  Double_t fFinalMax = fHists[0]->GetBinContent(fHists[0]->GetMaximumBin());

  for (Int_t i = 1; i < nHists; i++) {
    Double_t fLocalMin = fHists[i]->GetBinContent(fHists[i]->GetMinimumBin());
    Double_t fLocalMax = fHists[i]->GetBinContent(fHists[i]->GetMaximumBin());

    if (fLocalMin < fFinalMin) fFinalMin = fLocalMin;
    if (fLocalMax > fFinalMax) fFinalMax = fLocalMax;
  }
  // FIXME add extra space here
  Double_t fInitialRange = fFinalMax - fFinalMin;
  fFinalMax += 0.1 * fInitialRange;
  fFinalMin -= 0.1 * fInitialRange;

  *fMin=fFinalMin;
  *fMax=fFinalMax;
}


void ZoomCommonMinMax(std::vector<TH1D *> fHists) {
 
  double fMin = 0;
  double fMax = 0;

  FindCommonMinMax(fHists,&fMin,&fMax);

  fHists[0]->GetYaxis()->SetRangeUser(fMin,fMax);
}


#endif

