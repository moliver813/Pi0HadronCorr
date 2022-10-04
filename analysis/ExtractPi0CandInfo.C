#include "TFile.h"
#include "TList.h"


// 
//void ExtractEventPool(TString listName = "AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_histos", TString outFilePath = "TriggerPool.root",TString inFilePath = "AnalysisResults.root") {


void ExtractPi0CandInfo(TString listName = "AliAnalysisTaskGAPi0Candv_Pi0H_SE_tracks_caloClusters_HadCorr_Cent2_histos", TString outFilePath = "AnalysisResults_Pi0Cands.root", TString inFilePath = "AnalysisResults.root") {


  TFile * inFile = TFile::Open(inFilePath,"READ");
  if (!inFile) {
    fprintf(stderr,"Missing %s\n",inFilePath.Data());
    return;
  }

  //TString outFilePath = "TriggerPool.root";

	TList * list = (TList *) inFile->Get(listName);
	if (!list) {
		fprintf(stderr,"Could not find list %s\n",listName.Data());
		return;
	}




	TList * listClone = (TList *) list->Clone("clone");
	TFile * outFile = new TFile(outFilePath.Data(),"RECREATE");
	outFile->Add(listClone);
	outFile->Write();







	cout << "Success" << endl;
	return ;



}
