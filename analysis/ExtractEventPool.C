#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TGraphErrors.h"

#include "AliEventPoolManager.h"


void ExtractEventPool(TString listName = "AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_histos", TString outFilePath = "TriggerPool.root",TString inFilePath = "AnalysisResults.root") {

	TFile * inFile = TFile::Open(inFilePath,"READ");
	if (!inFile) {
		fprintf(stderr,"Missing %s\n",inFilePath.Data());
		return;
	}

	//TString outFilePath = "TriggerPool.root";

	const char * type = "Pi0H";
	Bool_t SearchForList = false;

  TList * keys = inFile->GetListOfKeys();
  TObject * obj;
  TIter next(keys);
//  TString listName = "AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_histos";
//  TString listName = "AliAnalysisTask_Pi0H_SB_SE_tracks_caloClusters_SB6_histos";
	if (SearchForList) {
		while ((obj = next())) { // extra parentheses to keep root 6 from complaining 
			TString name = obj->GetName();
			printf("Found object: %s\n",name.Data());
			if (name.Contains(type)) {
				listName = name;
			}
		}
	}
	printf("Using list %s\n",listName.Data());
  TList * list = (TList * ) inFile->Get(listName);
  if (!list) {
    fprintf(stderr,"List %s not found!\n",listName.Data());
    return;
  }

	AliEventPoolManager * fPoolManager = 0;
	fPoolManager = (AliEventPoolManager *) list->FindObject("AliEventPoolManager");
	if (!fPoolManager) {
		fprintf(stderr,"EventPoolManager not found!\n");
		return;
	}

	TFile * outFile = TFile::Open(outFilePath,"RECREATE");
	outFile->Add(fPoolManager);
	outFile->Write();
	outFile->Close();

}
