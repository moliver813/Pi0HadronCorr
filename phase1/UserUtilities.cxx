

#include "UserUtilities.h"

/**
  * Normalize a histogram to unity
  */
void NormalizeHist(TH1 * hist) {
  if (hist == 0) return;
  double sum = hist->Integral("width");
  if (sum == 0) return;
  hist->Scale(1./sum);
}

/**
  * Normalize a 1D histogram to its bin sizes (to get a spectrum dN/dx)
  */
void NormalizeHistByBinWidth(TH1 * hist) {
  if (hist == 0) return;
  int nBins = hist->GetNbinsX();
  for (int i = 1; i < nBins; i++) {
    double localValue = hist->GetBinContent(i);
    double localError = hist->GetBinError(i);
    double localBinWidth = hist->GetBinWidth(i);
    if (localBinWidth == 0) continue;
    localValue /= localBinWidth;
    localError /= localBinWidth;
    hist->SetBinContent(i,localValue);
    hist->SetBinError(i,localError);
  }
}


/**
  * Normalize a 2D Histogram by col
  */
TH2F * Normalize2DHistByCol(TH2F * inHist) {
  TH2F * normHist = (TH2F *) inHist->Clone(Form("%sColNorm",inHist->GetName()));
  normHist->SetTitle(Form("%s (Normalized-by-column)",inHist->GetTitle()));
  TH2F * scaleHist = (TH2F *) inHist->Clone(Form("%sScale",inHist->GetName()));
  scaleHist->SetTitle(Form("%s (Scale)",inHist->GetTitle()));
  scaleHist->Sumw2(0);
  int nbinsx = inHist->GetNbinsX();
  int nbinsy = inHist->GetNbinsY();

  for (int i = 0; i < nbinsx+1; i++) {
    double integral = inHist->Integral(i,i,0,nbinsy+1,"");

  //  double integral = inHist->Integral(i,i,0,nbinsy+1,"width");
//    double integral = inHist->Integral(i,i+1,0,nbinsy+1,"width");

    integral = integral ? 1./integral : 0; // lol
    for (int j = 0; j < nbinsy+1; j++) {
      scaleHist->SetBinContent(i,j,integral);
      scaleHist->SetBinError(i,j,0);
    }
  }

  normHist->Multiply(scaleHist);
  delete scaleHist;
  return normHist;
}

/**
  * Normalize a 2D Histogram by row
  */
TH2F * Normalize2DHistByRow(TH2F * inHist) {
  TH2F * normHist = (TH2F *) inHist->Clone(Form("%sRowNorm",inHist->GetName()));
  normHist->SetTitle(Form("%s (Normalized-by-row)",inHist->GetTitle()));
  TH2F * scaleHist = (TH2F *) inHist->Clone(Form("%sScale",inHist->GetName()));
  scaleHist->SetTitle(Form("%s (Scale)",inHist->GetTitle()));
  scaleHist->Sumw2(0);
  int nbinsx = inHist->GetNbinsX();
  int nbinsy = inHist->GetNbinsY();

  for (int i = 0; i < nbinsy+1; i++) {

    double integral = inHist->Integral(0,nbinsx+1,i,i,"");
    integral = integral ? 1./integral : 0;

    for (int j = 0; j < nbinsx+1; j++) {
      scaleHist->SetBinContent(j,i,integral);
      scaleHist->SetBinError(j,i,0);
    }
  }

  normHist->Multiply(scaleHist);
  delete scaleHist;
  return normHist;
}


