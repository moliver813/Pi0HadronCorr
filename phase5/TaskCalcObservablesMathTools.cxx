#ifndef TASKCALCOBSERVABLESMATHTOOLS_CXX
#define TASKCALCOBSERVABLESMATHTOOLS_CXX

#include <TGraphErrors.h>

#include <vector>
#include <numeric>

using namespace std;

/**
  * Calculate the systematic uncertainty from whatever is varied among the fGraphArray
  * If given a non-null pointer, will save the results there, such that that graph has the systematic errors but original points.
  * Returns a new TGraphErrors with the systematic uncertainties, but the point values are the mean of these variations
  */
//TGraphErrors * ProduceSystematicFromGraphs(vector<TGraphErrors*> fGraphArray) {

TGraphErrors * ProduceSystematicFromGraphs(vector<TGraphErrors*> fGraphArray, TGraphErrors * fInputGraph) {
  TGraphErrors * fMeanWithSysError = 0;

  int nVar = fGraphArray.size();
  if (nVar == 0) {
    fprintf(stderr,"Empty variant array!\n");
    return fMeanWithSysError;
  }

  if (fInputGraph != 0) {
    fMeanWithSysError = (TGraphErrors *) fInputGraph->Clone(Form("%s_SysErrorWithMean",fInputGraph->GetName()));
  } else { // just copy the first graph
    fMeanWithSysError = (TGraphErrors *) fGraphArray[0]->Clone(Form("%s_SysErrorWithMean",fGraphArray[0]->GetName()));
  }

  int nPoints = fMeanWithSysError->GetN();
  for (int i = 0; i < nPoints; i++) {
    double fMean = 0;
    double fVar = 0;

    double fXValue = fMeanWithSysError->GetX()[i];
    double fXError = fMeanWithSysError->GetEX()[i];

    vector<double> fValues = {};
    printf("\ndebug:Sys (bin %d) : ",i);
    for (int iVar = 0; iVar < nVar; iVar++) {
      fValues.push_back(fGraphArray[iVar]->GetY()[i]);
      printf("%.2e ",fValues[iVar]);
    }
    printf("\n");

    // Learned this from stackoverflow
    double sum = std::accumulate(fValues.begin(),fValues.end(),0.0);
    fMean = sum / nVar;

    double sumSq = std::inner_product(fValues.begin(),fValues.end(),fValues.begin(),0.0);
    fVar = std::sqrt(sumSq / nVar - fMean * fMean);

    fMeanWithSysError->SetPoint(i,fXValue,fMean);
    fMeanWithSysError->SetPointError(i,fXError,fVar);
    if (fInputGraph) {
      printf("Setting val +- error to %f %f\n",fMean,fVar);
      fInputGraph->SetPointError(i,fInputGraph->GetEX()[i],fVar);
      // FIXME testing, 
      printf(" Setting x,y = %f,%f\n",fXValue,fMean);
 //     fInputGraph->SetPoint(i,fXValue,fMean);
    }
  }




  return fMeanWithSysError;
}


TH1D * ProduceSystematicFromHists(vector<TH1D*> fHistArray, TH1D * fInputHist) {
  TH1D * fMeanWithSysError = 0;

  int nVar = fHistArray.size();



  return fMeanWithSysError;
}

#endif
