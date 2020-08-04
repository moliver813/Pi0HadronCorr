#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TGraphErrors.h"

#include "AliEventPoolManager.h"

const Int_t cWidth = 600;
const Int_t cHeight = 900;

void SetStyle() {
//  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.07);//0.05 //..ELI
 // gStyle->SetPadBottomMargin(0.18);//0.15
  gStyle->SetPadBottomMargin(0.12);//0.15
  gStyle->SetPadRightMargin(0.04);
//  gStyle->SetPadLeftMargin(0.21);
  gStyle->SetPadLeftMargin(0.14);
  TGaxis::SetMaxDigits(3); //

}

void GetEventPoolInfo(TString listName = "AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_histos",TString inFilePath = "AnalysisResults.root") {

	//TString inFilePath = "AnalysisResults.root";	


	TFile * inFile = TFile::Open(inFilePath,"READ");
	if (!inFile) {
		fprintf(stderr,"Missing %s\n",inFilePath.Data());
		return;
	}
	SetStyle();	

	const char * type = "caloClusters_histos";

  TList * keys = inFile->GetListOfKeys();
  TObject * obj;
  TIter next(keys);
//  TString listName = "AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_histos";
//  TString listName = "AliAnalysisTask_Pi0H_SB_SE_tracks_caloClusters_SB3_histos";
  while ((obj = next())) { // extra parentheses to keep root 6 from complaining 
    TString name = obj->GetName();
    printf("Found object: %s\n",name.Data());
    if (name.Contains(type)) {
      listName = name;
    }
  }
	printf("Using list %s\n",listName.Data());
  TList * list = (TList * ) inFile->Get(listName);

	AliEventPoolManager * fPoolManager = 0;

  if (!list) {
    printf("List %s not found!\n",listName.Data());
		printf("Trying to find the event pool in the main directory\n");
		fPoolManager = (AliEventPoolManager *) inFile->Get("AliEventPoolManager");
  } else {
		fPoolManager = (AliEventPoolManager *) list->FindObject("AliEventPoolManager");
	}
	if (!fPoolManager) {
		fprintf(stderr,"EventPoolManager not found!\n");
		return;
	}

	Int_t nTotalBins = fPoolManager->GetNumberOfAllBins();
	Int_t nMultBins  = fPoolManager->GetNumberOfMultBins();
	Int_t nZVtxBins  = fPoolManager->GetNumberOfZVtxBins();
	Int_t nEPBins    = fPoolManager->GetNumberOfPsiBins();
	Int_t nPtBins    = fPoolManager->GetNumberOfPtBins();
	printf("EventPoolManager has %d bins in total.\n",nTotalBins);
	printf("   %d Mult bins\n",nMultBins);
	printf("   %d nZVt bins\n",nZVtxBins);
	printf("   %d EP   bins\n",nEPBins);
	printf("   %d Pt   bins\n",nPtBins);

	Int_t nTotalEvents = 0;
	for (Int_t i = 0; i < nMultBins; i++) {
		for (Int_t j = 0; j < nZVtxBins; j++) {
			for (Int_t k = 0; k < nEPBins; k++) {
				for (Int_t l = 0; l < nPtBins; l++) {
					AliEventPool * pool = fPoolManager->GetEventPool(i,j,k,l);
					Int_t nEv = pool->GetCurrentNEvents();
					printf("[%2d %2d %2d %2d] = %d\n",i,j,k,l,nEv);
					nTotalEvents+= nEv; 
				}
			}
		}
	}
	printf("Total Events Stored: %d\n",nTotalEvents);
}
