#ifndef TASKCALCOBSERVABLESMATHTOOLS_CXX
#define TASKCALCOBSERVABLESMATHTOOLS_CXX

#include <TGraphErrors.h>


#include "TaskCalcObservables.h"


#include <vector>
#include <numeric>

using namespace std;

/**
  * Calculate the systematic uncertainty from whatever is varied among the fGraphArray
  * If given a non-null pointer, will save the results there, such that that graph has the systematic errors but original points.
  * Returns a new TGraphErrors with the systematic uncertainties, but the point values are the mean of these variations
  */
//TGraphErrors * ProduceSystematicFromGraphs(vector<TGraphErrors*> fGraphArray) {


double CalculateVariance(vector<double> fValues) {
    int nValues = fValues.size();
    // Learned this from stackoverflow
    double sum = std::accumulate(fValues.begin(),fValues.end(),0.0);
    double fMean = sum / nValues;
    double sumSq = std::inner_product(fValues.begin(),fValues.end(),fValues.begin(),0.0);
    //printf("Raw SumSq/nVar = %f, fMean^2 = %f\n",sumSq/nValues,fMean*fMean);
    return std::sqrt(std::abs(sumSq / nValues - fMean * fMean));
}


TGraphErrors * TaskCalcObservables::ProduceSystematicFromGraphs(vector<TGraphErrors*> fGraphArray, TGraphErrors * fInputGraph) {
  TGraphErrors * fMeanWithSysError = 0;

  TCanvas * canv = new TCanvas();


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

  TGraphErrors * fMinValueGraph = (TGraphErrors *) fMeanWithSysError->Clone(Form("%s_MinValue",fMeanWithSysError->GetName()));
  TGraphErrors * fMaxValueGraph = (TGraphErrors *) fMeanWithSysError->Clone(Form("%s_MaxValue",fMeanWithSysError->GetName()));



  int nPoints = fMeanWithSysError->GetN();
  for (int i = 0; i < nPoints; i++) {
    double fMean = 0;
    double fVar = 0;

    double fXValue = fMeanWithSysError->GetX()[i];
    double fXError = fMeanWithSysError->GetEX()[i];

    vector<double> fValues = {};
    vector<double> fValueErrors = {};
    //printf("\ndebug:Sys (bin %d) : ",i);
    for (int iVar = 0; iVar < nVar; iVar++) {

      // Temporary outlier removal
     // if (TMath::Abs(fGraphArray[iVar]->GetY()[i]) < 4) fValues.push_back(fGraphArray[iVar]->GetY()[i]);
      fValues.push_back(fGraphArray[iVar]->GetY()[i]);
      fValueErrors.push_back(fGraphArray[iVar]->GetEY()[i]);
      //printf("%.2e ",fValues[iVar]);
 //     printf("%.2f #pm %.2f \n",fValues[iVar],fGraphArray[iVar]->GetEY()[i]);
    }
    //printf("\n");

    int nValues = fValues.size();

    // Learned this from stackoverflow
    double sum = std::accumulate(fValues.begin(),fValues.end(),0.0);
    fMean = sum / nValues;
    //fMean = sum / nVar;
    double sumSq = std::inner_product(fValues.begin(),fValues.end(),fValues.begin(),0.0);
   // printf("Raw SumSq/nVar = %f, fMean^2 = %f\n",sumSq/nValues,fMean*fMean);

    //bool bEnableWeightedMean = false; // Disabling for remaking preliminary
    bool bEnableWeightedMean = true;

    double sumOfWeights = 0;
    double weightedSumSq = 0;
    double weightedMean = 0;
    
    if (bEnableWeightedMean) {
      double weightedSum = 0;
      for (int iVar = 0; iVar < nVar; iVar++) {
        double value = fValues[iVar];
        double weight = fValueErrors[iVar];
        if (weight != 0) weight = 1./weight;

        sumOfWeights += weight;
        weightedSum += weight * value;
        weightedSumSq += weight * value * value;
      }
      if (sumOfWeights != 0) {
        weightedMean = weightedSum / sumOfWeights;
        weightedSumSq = weightedSumSq / sumOfWeights;
      }
    }



    //printf("SumSq/nVar = %f, fMean^2 = %f\n",sumSq/nVar,fMean*fMean);
    //printf("SumOfWeights = %.2e, WeightedMean = %.2e, WeightedSumSq = %.2e\n",sumOfWeights,weightedMean,weightedSumSq);


    //fVar = std::sqrt(sumSq / nVar - fMean * fMean);
    //fVar = std::sqrt(std::abs(sumSq / nVar - fMean * fMean));
    if (bEnableWeightedMean) {
      fVar = std::sqrt(std::abs(weightedSumSq - weightedMean * weightedMean));
    } else {
      fVar = std::sqrt(std::abs(sumSq / nValues - fMean * fMean));
    }

    fMeanWithSysError->SetPoint(i,fXValue,fMean);
    fMeanWithSysError->SetPointError(i,fXError,fVar);
    if (fInputGraph) {
      printf("Setting val +- error to %f +- %f\n",fMean,fVar);
      //fInputGraph->SetPointError(i,fXError,fVar);
      fInputGraph->SetPointError(i,fInputGraph->GetEX()[i],fVar);
      // FIXME testing, 
      printf(" Setting x,y = %f,%f\n",fXValue,fMean);
 //     fInputGraph->SetPoint(i,fXValue,fMean);
    }

    double min = *std::min_element(fValues.begin(),fValues.end());
    double max = *std::max_element(fValues.begin(),fValues.end());

    printf("Found min,max = %.3f,%.3f\n",min,max);   

    fMinValueGraph->SetPoint(i,fXValue,min);
    fMaxValueGraph->SetPoint(i,fXValue,max);

    fMinValueGraph->SetPointError(i,0,0);
    fMaxValueGraph->SetPointError(i,0,0);

  }

  float legX1 = 0.45;
  float legY1 = 0.2;
  float legX2 = 0.45;
  float legY2 = 0.2;

  fInputGraph->SetMarkerStyle(kFullSquare);
  fMeanWithSysError->SetMarkerStyle(kOpenSquare);
  fMeanWithSysError->SetMarkerColor(kRed+1);
  fMeanWithSysError->SetLineStyle(kRed+1);

  fMinValueGraph->SetMarkerStyle(kOpenCircle);
  fMinValueGraph->SetMarkerColor(kBlue+1);
  fMinValueGraph->SetLineColor(kBlue+1);
  fMaxValueGraph->SetMarkerStyle(kOpenCircle);
  fMaxValueGraph->SetMarkerColor(kBlue+4);
  fMaxValueGraph->SetLineColor(kBlue+4);

  TMultiGraph * mg = new TMultiGraph();
  TLegend * leg = new TLegend(legX1,legY1,legX2,legY2);
  mg->Add(fInputGraph);

  mg->Add(fMinValueGraph);
  mg->Add(fMaxValueGraph);

  mg->Add(fMeanWithSysError);

  mg->Draw("ALP");
  leg->SetHeader(Form("%s",fInputGraph->GetName()),"c");
  leg->AddEntry(fInputGraph,"Central Value");

  leg->AddEntry(fMinValueGraph,"Min. of Variants");
  leg->AddEntry(fMaxValueGraph,"Max. of Variants");

  leg->AddEntry(fMeanWithSysError,"Average of variations, with errors from variance");
  leg->Draw("SAME");

  canv->SetLogy(1);

  canv->SetGridx(1);
  canv->SetGridy(1);

  //canv->BuildLegend();
  canv->Print(Form("%s/QA/ErrofFromGraphs_QA_%s.pdf",fOutputDir.Data(),fInputGraph->GetName()));
  canv->Print(Form("%s/QA/ErrofFromGraphs_QA_%s.png",fOutputDir.Data(),fInputGraph->GetName()));
  canv->Print(Form("%s/CFiles/ErrofFromGraphs_QA_%s.C",fOutputDir.Data(),fInputGraph->GetName()));


  return fMeanWithSysError;
}


TH1D * TaskCalcObservables::ProduceSystematicFromHists(vector<TH1D*> fHistArray, TH1D * fInputHist) {
  TH1D * fMeanWithSysError = 0;

  int nVar = fHistArray.size();



  return fMeanWithSysError;
}


vector<vector<TGraphErrors *>> Transpose2DTGraphErrArray(vector<vector<TGraphErrors *>> input) {

  vector<vector<TGraphErrors *>> output = {};
  int nRow = input.size();
  int nCol = input[0].size();
  cout<<"Debug::Transpose"<<endl;
  printf("  nRows=%d, nCol=%d\n",nRow,nCol);
  for (int i = 0; i < nRow; i++) {
    printf("  Row %d has %d entries\n",i,(int)input[i].size());
  }

  for (int j = 0; j < nCol; j++) {
    vector<TGraphErrors *> localVector = {};
    for (int i = 0; i < nRow; i++) {
      localVector.push_back(input[i][j]);
    }
    output.push_back(localVector);
  }

  return output;
}

// outermost (first) index is the source. I want to make that the last index
vector<vector<vector<TH1D *>>> Transpose3DTH1DArray(vector<vector<vector<TH1D *>>> input) {
  vector<vector<vector<TH1D *>>> output = {};
  int nSources = input.size();
  printf("  nSources=%d\n",nSources);
  int nObs = input[0].size();
  printf("  nSources=%d, nObs=%d\n",nSources,nObs);
  int nEPs = input[0][0].size();
  printf("  nSources=%d, nObs=%d, nEPs=%d\n",nSources,nObs,nEPs);

  cout<<"Debug::Transpose3D"<<endl;
  for (int i = 0; i < nObs; i++) {
    vector<vector<TH1D *>> localArray = {};
    for (int j = 0; j < nEPs; j++) {
      vector<TH1D *> localVector = {};
      for (int k = 0; k < nSources; k++) {
        localVector.push_back(input[k][i][j]);
      }
      localArray.push_back(localVector);
    }
    output.push_back(localArray);
  }

  return output;
}


#endif
