

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

