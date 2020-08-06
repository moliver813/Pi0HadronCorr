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
#include "TEfficiency.h"

const float kTol = 1e-6;

const float lX1 = 0.56;
const float lX2 = 0.9;
const float lY1 = 0.20;
const float lY2 = 0.9;

//const int nPtBins = 10; // 3,4 then standard
//std::vector <double> fPtBins = {3,4,5,7,9,11,14,17,20,22,30};

//const int nPtBins = 7;
//std::vector <double> fPtBins = {5,7,9,11,13,15,17,19}; // Normal for this

// So Lazy
//std::vector <double> fPtBinsMid = {3.5,4.5,6,8,10,12.5,15.5,18.5,21,26,35};
//std::vector <double> fPtBinsMid = {6,8,10,12,14,16,18}; // Normal for this

//const int nPtBins = 7;
//std::vector <double> fPtBins =   {4,5,7,9,11,14,17,20}; // Analysis Bins
//std::vector <double> fPtBinsMid = {4.5,6,8,10,12.5,15.5,18.5};

// Analysis Bins
const int nPtBins = 7;
std::vector <double> fPtBins    = {5,7,9,11,14,17,20}; // Analysis Bins
std::vector <double> fPtBinsMid = {6,8,10,12.5,15.5,18.5};

// Finer bins
//const int nPtBins = 10;
//std::vector <double> fPtBins    = {4  ,5  ,6  ,7  ,8  ,9  ,10  ,11,15  ,16 ,17};
//std::vector <double> fPtBinsMid = {4.5,5.5,6.5,7.5,8.5,9.5,10.5,13,15.5,16.5};

const bool bGlobalParam = true; // Attempt Global fits of the parameters

//const int nUseThetaBins = 23; //26
const int nFirstThetaBin = 2; // 5 for same event, 2 for mixed eventsa
//const int nFirstThetaBin = 2; // 5 for same event, 2 for mixed eventsa
//const int nFirstThetaBin = 5; // 5 for same event, 2 for mixed eventsa
const int nSkipFinalBins = 0; // Skip some of the last

int iSkipAlternateBins = 0; // 0 no skipping, 1 skip odd, 2 skip even

int nUseThetaBins = 15;
//int nUseThetaBins = 4;//15 //26
int nLastThetaBin = nFirstThetaBin + nUseThetaBins;

const Int_t iType = 1;  // 1 for same event, 2 for RotBkg, 3 for mixed event, 4 for pos swap
const Bool_t bSelectMC = false; // use only specific type of pair
const	Int_t mcType = 3; //2 for pi0s , 0 for true background
	// 3 for correlated background (not including single particle -> 2 cluster background
// 0 is truely uncorrelated background
const int nRebinMass = 5;// 2 //5

const float kMinMassRange = 0;
const float kMaxMassRange = 0.55;

const int cWidth = 1200;
const int cHeight = 900;

const int kFuncColor1 = kRed; // not used yet
const int kFuncColor2 = kViolet;

const float kPaveTextX1 = 0.65;
const float kPaveTextY1 = 0.3;
const float kPaveTextX2 = 0.95;
const float kPaveTextY2 = 0.6;

void NormalizeHistoHighMass(TH1 * histo);

// This is derived from fitRatio_3
TF1 * fitEfficiency(TEfficiency * hRatio, std::string name, double pT, double theta);

TF1 * fitRatio(TH1F * hRatio, std::string name);
TF1 * fitRatio_2(TH1F * hRatio, std::string name, double pT, double theta);
TF1 * fitRatio_3(TH1F * hRatio, std::string name, double pT, double theta);
TF1 * fitRatio_4(TH1F * hRatio, std::string name, double pT, double theta);
TF1 * fitMass_Flat(TH1F * hRatio, std::string name, double pT, double theta);
TF1 * fitMass_Int_4(TH1F * hist, std::string name, double pT, double theta, double ECut, double norm);

TF2 * FitGlobalModel();
TF2 * CreatePtVsMassGlobalModel();

void SetStyle() {
  gStyle->SetOptStat(0);
//	gStyle->SetOptFit(111);
//  gStyle->SetTitleOffset(0.6,"X"); // 0.7 
//  gStyle->SetLabelSize(0.04,"X");
//  gStyle->SetTitleSize(0.01,"all"); //0.055
//  gStyle->SetTitleSize(0.045,"X"); // 0.05

//  gStyle->SetPadTopMargin(0.07);//0.05 //..ELI
//  gStyle->SetPadBottomMargin(0.18);//0.15
//  gStyle->SetPadRightMargin(0.04);
//  gStyle->SetPadLeftMargin(0.21);
  TGaxis::SetMaxDigits(3); //

}

int GetCustomColor(int x) {
	// r,g,b in [0,255]
	//int r = (13 * x) % 256 ;
	//would be simpler to use GetColor(Float_t r, Float_t g, Float_t b)
	float radius  = 1.0 - 0.3 * TMath::Power(TMath::Sin(0.1*x),2);
	float theta_r = 3.1415 * 0.12 * x - 3.1415; // 0.06
	float theta_g = 3.1415 * 0.06 * x + 1.57;   // 0.04
	float theta_b = 3.1415 * 0.04 * x - 3.1415; // 0.02
	int r = round(radius*(128 + 127.*TMath::Cos(theta_r)));
	int g = round(radius*(128 + 127.*TMath::Cos(theta_g)));
	int b = round(radius*(128 + 127.*TMath::Cos(theta_b)));
//	int r = 128 + round(127.*TMath::Cos(theta_r));
//	int g = 128 + round(127.*TMath::Cos(theta_g));
//	int b = 128 + round(127.*TMath::Cos(theta_b));
//	printf("r,g,b = (%d,%d,%d)\n",r,g,b);
	return TColor::GetColor(r,g,b);
}

void angleEffect(int LambdaBinLow = 1, int LambdaBinHigh = 5, int EnergyBinLow = 1, int EnergyBinHigh = -2, int AsymBinLow = 1, int AsymBinHigh = -1, string inputListName = "AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_TrkClsVeto_Cent0_histos") {
//void angleEffect(int LambdaBinLow = 1, int LambdaBinHigh = 5, int EnergyBinLow = 3, int EnergyBinHigh = -2, int AsymBinLow = 1, int AsymBinHigh = -1) {
	TString inFilePath = "AnalysisResults.root";	
	TFile * inFile = TFile::Open(inFilePath,"READ");
	if (!inFile) {
		fprintf(stderr,"Missing %s\n",inFilePath.Data());
		return;
	}
	SetStyle();	

	//const char * type = "Pi0H";
	const char * type = "Cent0";
	//const char * type = "MB";

  TList * keys = inFile->GetListOfKeys();
  TObject * obj;
  TIter next(keys);
//  TString listName = "AliAnalysisTask_Pi0H_ME_tracks_caloClusters_histos";
	//TString listName = "AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_TrkClsVeto_Cent0_histos";
	TString listName(inputListName);
  //TString listName = "Cent0";

	
  while ((obj = next())) { // extra parentheses to keep root 6 from complaining 
    TString name = obj->GetName();
    printf("Found object: %s\n",name.Data());
  //  if (name.Contains(type)) {
  //    listName = name;
  //  }
  }
	printf("Using list %s\n",listName.Data());
  TList * list = (TList * ) inFile->Get(listName);
  if (!list) {
    fprintf(stderr,"List %s not found!\n",listName.Data());
    return;
  }

	// Temporary, can be removed
	if (true) {
		TH1D * fClusEnergy = (TH1D *) list->FindObject("ClusEnergy");
		if (!fClusEnergy) return;
		TCanvas * cClusEnergy = new TCanvas("cClusEnergy","cClusEnergy");

		TF1 * fitClusEnergy = new TF1("fitClusEnergy","[0]*TMath::Exp(-[1]*x)",0.3,fClusEnergy->GetXaxis()->GetXmax());
		fitClusEnergy->SetParameter(0,fClusEnergy->GetBinContent(fClusEnergy->GetMaximumBin()));
		fitClusEnergy->SetParLimits(0,0.0,10*fitClusEnergy->GetParameter(0));
		fitClusEnergy->SetParameter(1,1.);
		fitClusEnergy->SetParLimits(1,0.,20.);
		
		TF1 * fitClusEnergy2 = new TF1("fitClusEnergy2","[0]*TMath::Power(x,-[1])",1.5,20.);
		//TF1 * fitClusEnergy2 = new TF1("fitClusEnergy2","[0]*TMath::Power(x,-[1])",1.5,fClusEnergy->GetXaxis()->GetXmax());
		fitClusEnergy2->SetLineColor(kCyan);
		fitClusEnergy2->SetParameter(0,fClusEnergy->GetBinContent(fClusEnergy->GetMaximumBin()));
		fitClusEnergy2->SetParLimits(0,0.0,10*fitClusEnergy2->GetParameter(0));
		fitClusEnergy2->SetParameter(1,5.);
		fitClusEnergy2->SetParLimits(1,0.5,8.);

		TF1 * fitClusEnergy3 = new TF1("fitClusEnergy3","[0]*TMath::Exp(-[1]*x)+[2]*TMath::Power(x,-[3])",1.5,40.);
		fitClusEnergy3->SetParameter(0,fClusEnergy->GetBinContent(fClusEnergy->GetMaximumBin()));
		fitClusEnergy3->SetParLimits(0,0.0,10*fitClusEnergy3->GetParameter(0));
		fitClusEnergy3->SetParameter(1,1.);
		fitClusEnergy3->SetParLimits(1,0.,20.);
		fitClusEnergy3->SetParameter(2,fClusEnergy->GetBinContent(fClusEnergy->GetMaximumBin()));
		fitClusEnergy3->SetParLimits(2,0.0,10*fitClusEnergy3->GetParameter(0));
		fitClusEnergy3->SetParameter(3,5.);
		fitClusEnergy3->SetParLimits(3,0.5,8.);
		fitClusEnergy3->SetLineColor(kCyan);

		//fClusEnergy->Fit(fitClusEnergy);	
		//fClusEnergy->Fit(fitClusEnergy2,"R");	
		fClusEnergy->Fit(fitClusEnergy3,"R");	

		fClusEnergy->Draw();

		cClusEnergy->SetLogy();

		cClusEnergy->Print("ClusEnergy.pdf");
		cClusEnergy->Print("ClusEnergy.C");
	}

  Int_t nLambdaBinHigh = LambdaBinHigh;
  Int_t nEnergyBinHigh = EnergyBinHigh;

	THnSparse * Pi0Cands = (THnSparse *) list->FindObject("Pi0Cands");
	if (!Pi0Cands) {
		fprintf(stderr,"Missing Pi0Cands THnSparse!\n");
		return;
	}

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

	double fClusterEnergyCut = Pi0Cands->GetAxis(4)->GetBinLowEdge(EnergyBinLow);

	// Apply Asymmetry Cut
	int nAsymBinHigh = AsymBinHigh;
	if (nAsymBinHigh == -1) nAsymBinHigh = Pi0Cands->GetAxis(5)->GetNbins();
	if (nAsymBinHigh == -2) nAsymBinHigh = Pi0Cands->GetAxis(5)->GetNbins() + 1;
	printf("Asym Range:   [%.2f - %.2f] \t\tBins: %d %d\n",Pi0Cands->GetAxis(5)->GetBinLowEdge(AsymBinLow),Pi0Cands->GetAxis(5)->GetBinUpEdge(nAsymBinHigh),AsymBinLow,nAsymBinHigh);
	Pi0Cands->GetAxis(5)->SetRange(AsymBinLow,nAsymBinHigh);

	// Selecting data or rotBkg or Mixed
	Pi0Cands->GetAxis(6)->SetRange(iType,iType);

	// Selecting MC Type:
	if (bSelectMC) {
		int iMCAxis = 8;
		if (Pi0Cands->GetNdimensions() < 9) iMCAxis = 7; 
		if (mcType <= 2) {
			printf("Will use only MC type %d\n",mcType);
			//Pi0Cands->GetAxis(7)->SetRange(mcType+1,mcType+1);
			Pi0Cands->GetAxis(iMCAxis)->SetRange(mcType+1,mcType+1);
		} else if (mcType == 3) {
			int iCorrBkgMin = 3; // pi0 dalitz 
			int iCorrBkgMax = 9; // other shared ancestor
			// unfortunately, this misses the single particle background
			printf("Will use on correlated background types %d through %d.\n",iCorrBkgMin,iCorrBkgMax);
			Pi0Cands->GetAxis(iMCAxis)->SetRange(iCorrBkgMin+1,iCorrBkgMax+1);
		}
	}

	// Seeing how many pairs survived cuts:
	TH1F * hPt = (TH1F *) Pi0Cands->Projection(1);
	double nClusterPairs = hPt->Integral();
	printf("Found %.1f cluster pairs\n",nClusterPairs);

	// Looping over Theta Bins
	int nThetaBins = Pi0Cands->GetAxis(2)->GetNbins();
	// Changing the global nUseThetaBins variable, if necessary
//	if (nThetaBins < 15) {
//		nUseThetaBins = 4;
//		nLastThetaBin = nFirstThetaBin + nUseThetaBins;
//	}
//	if (iType 


	printf("Using bins %d - %d (out of %d bins)\n",nFirstThetaBin,nLastThetaBin,nThetaBins);
	std::vector<Double_t> fMeanTheta = {};
	std::vector<Double_t> fLowTheta  = {};

	std::vector<TH2F *> hMassPtThetaBins = {};
	std::vector<std::vector<TH1F *>> hMassPtBinThetaBin = {};

	// Saving the ones with cuts .  NVM
	///std::vector<std::vector<TH1F *>> hMassPtBinThetaBin = {};
	// Ratios
	std::vector<std::vector<TEfficiency *>> hMassEffPtBinThetaBin = {};
	std::vector<std::vector<TH1F *>> hMassRatioPtBinThetaBin = {};

	TH2F * hMassPtNoThetaCut = (TH2F *) Pi0Cands->Projection(0,1,"e");
	hMassPtNoThetaCut->SetName("hMassPt_NoThetaCut");
	std::vector<TH1F *> hMassPtBinNoThetaCut = {};
		
	Int_t nRealNThetaBins = Pi0Cands->GetAxis(2)->GetNbins();
	nUseThetaBins = min(nRealNThetaBins-1,nUseThetaBins);


	for (Int_t i = 0; i < nUseThetaBins; i++) {
		Pi0Cands->GetAxis(2)->SetRange(nFirstThetaBin + i,Pi0Cands->GetAxis(2)->FindBin(1.57)); // adding pi/2 as a max angle
//		Pi0Cands->GetAxis(2)->SetRange(nFirstThetaBin + i,nThetaBins); //correct
//		Pi0Cands->GetAxis(2)->SetRange(nFirstThetaBin + i,nFirstThetaBin + i);  // FIXME make it do the low cut only
		TH2F * localMassPt = (TH2F *) Pi0Cands->Projection(0,1,"e");
		localMassPt->SetName(Form("hMassPt_Theta_%d",i));
		hMassPtThetaBins.push_back(localMassPt);
		
		fMeanTheta.push_back(Pi0Cands->GetAxis(2)->GetBinCenter(nFirstThetaBin+i));
		fLowTheta.push_back(Pi0Cands->GetAxis(2)->GetBinLowEdge(nFirstThetaBin+i));
	}
	// UnZoom
	Pi0Cands->GetAxis(2)->SetRange(1,nThetaBins);

	// Style
//	std::vector<Int_t> k


	for (Int_t i = 0; i < nPtBins; i++) {
		Float_t lowPt = fPtBins[i];
		Float_t highPt = fPtBins[i+1];
		//Float_t meanPt = 0.5 * (lowPt + highPt);

		Int_t lowPtBin = Pi0Cands->GetAxis(0)->FindBin(lowPt);
		Int_t highPtBin = Pi0Cands->GetAxis(0)->FindBin(highPt);

		TH1F * localMassNoThetaCut = (TH1F *) hMassPtNoThetaCut->ProjectionX(Form("Mass_Pt_%d_NoThetaCut",i),lowPtBin,highPtBin);
		localMassNoThetaCut->Rebin(nRebinMass);
		hMassPtBinNoThetaCut.push_back(localMassNoThetaCut);

		std::vector<TH1F *> lMassThetaBin = {};
		std::vector<TH1F *> lMassRatioThetaBin = {};
		std::vector<TEfficiency *> lMassEffThetaBin = {};

		for (Int_t j = 0; j < nUseThetaBins; j++) {
			float lowTheta = fLowTheta[j]; 
			TH1F * localMass = (TH1F *) hMassPtThetaBins[j]->ProjectionX(Form("Mass_Pt_%d_Theta_%d",i,j),lowPtBin,highPtBin);
			localMass->SetTitle(Form("Mass (%.1f #leq p_{T} < %.1f)",lowPt,highPt));
			localMass->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/c^{2})");
			localMass->Rebin(nRebinMass);
//			printf("Debug: local integral = %f, no cut integral = %f\n",localMass->Integral(),localMassNoThetaCut->Integral());
			Color_t lColor = GetCustomColor(j);
			localMass->SetLineColor(lColor);
			localMass->SetMarkerColor(lColor);

		//	localMass->SetLineColor(kRed+j-7);
		//	localMass->SetMarkerColor(kRed+j-7);

			localMass->SetMarkerStyle(kOpenSquare);
			lMassThetaBin.push_back(localMass);

			// Old Method
			TH1F * localMassRatio = (TH1F * ) localMass->Clone(Form("MassRatio_Pt_%d_Theta_%d",i,j));
			//localMassRatio->Divide(localMassNoThetaCut);
			localMassRatio->Divide(localMassRatio,localMassNoThetaCut,1.0,1.0,"B");
			//  /Old Method

			TH1F * localCutDist = (TH1F * ) localMass->Clone(Form("MassEff_Pt_%d_Theta_%d",i,j));

			TEfficiency * localMassEff = new TEfficiency(*localMass,*localMassNoThetaCut);
			//TEfficiency * localMassRatio = new TEfficiency(*localCutDist,*localMassNoThetaCut);


			if (j == 0) localMassRatio->SetTitle(Form("Ratio (%.0f #leq p_{T} < %.0f GeV/c)",lowPt,highPt));
			else localMassRatio->SetTitle(Form("Ratio (%.0f #leq p_{T} < %.0f GeV/c) (%.3f #leq #theta_{c})",lowPt,highPt,lowTheta));
//			printf("integral of ratio = %f\n",localMassRatio->Integral());
			lMassRatioThetaBin.push_back(localMassRatio);
			lMassEffThetaBin.push_back(localMassEff);
		}
		hMassPtBinThetaBin.push_back(lMassThetaBin);
		hMassRatioPtBinThetaBin.push_back(lMassRatioThetaBin);
		hMassEffPtBinThetaBin.push_back(lMassEffThetaBin);

		//hMassCutDistPtBinThetaBin.push_back(lMassCutDistPtBinThetaBin);
	}

	std::vector<std::vector<TF1 *>> fMassFits = {};
	std::vector<std::vector<TF1 *>> fMassRatioFits = {};
	std::vector<std::vector<Double_t>> fFitChiSqNDF = {};  // [pt][theta]
	std::vector<std::vector<std::vector<Double_t>>> fFitParams = {};  // [param][pt][theta]
	std::vector<std::vector<std::vector<Double_t>>> fFitParamsUn = {};  // [param][pt][theta]

	// Fit Opening Angle Mass Spectral Ratios
	for (Int_t i = 0; i < nPtBins; i++) {
		std::vector<TF1 *> localFits = {};
		std::vector<Double_t > localChiSqNDF = {};
		Float_t lowPt = fPtBins[i];
		Float_t highPt = fPtBins[i+1];
		printf("==============================================\n");
		printf(" Starting ratio fits for pt bin [ %.2f , %.2f ]\n",lowPt,highPt);
		printf("==============================================\n");
		for (Int_t j = 0; j < nUseThetaBins; j++) {
			float lowTheta = fLowTheta[j]; 
			printf("Fitting ratio for theta cut = %f\n",lowTheta);

			TF1 * localPreFit = fitRatio_3(hMassRatioPtBinThetaBin[i][j],Form("RatioPreFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);

			// Fit the TEfficiency object with parameters from the ratio fit
			TF1 * localFit = localPreFit;
			//TF1 * localFit = fitEfficiency(hMassEffPtBinThetaBin[i][j],Form("RatioFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);


	//		printf("Starting Ratio Fit for pT \\in (%.1f,%.1f) GeV/c, theta_cut = %.3f\n",lowPt,highPt,lowTheta);
// fitRatio_2 is what I had before
//			TF1 * localFit = fitRatio_4(hMassRatioPtBinThetaBin[i][j],Form("RatioFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);
			//TF1 * localFit = fitRatio_3(hMassRatioPtBinThetaBin[i][j],Form("RatioFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);
//			TF1 * localFit = fitRatio_2(hMassRatioPtBinThetaBin[i][j],Form("RatioFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);
	//		TF1 * localFit = fitRatio(hMassRatioPtBinThetaBin[i][j],Form("RatioFit_Pt_%d_Theta_%d",i,j));
		//	localFit->SetLineColor(hMassRatioPtBinThetaBin[i][j]->GetLineColor());
	//		localFit->SetFillColor(hMassRatioPtBinThetaBin[i][j]->GetLineColor());
//			localFit->SetLineColor(kGreen);
			localFit->SetNpx(200);
			localFits.push_back(localFit);
			if (localFit->GetNDF() > 0) {
				localChiSqNDF.push_back(localFit->GetChisquare() / localFit->GetNDF());
			} else localChiSqNDF.push_back(-1);
		}
		fMassRatioFits.push_back(localFits);
		fFitChiSqNDF.push_back(localChiSqNDF);
	}

	// Fitting Mass Spectra
	for (Int_t i = 0; i < nPtBins; i++) {
		std::vector<TF1 *> localFits = {};
		Float_t lowPt = fPtBins[i];
		Float_t highPt = fPtBins[i+1];
		for (Int_t j = 0; j < nUseThetaBins; j++) {
			float lowTheta = fLowTheta[j];
	//		printf("Starting Mass  Fit for pT \\in (%.1f,%.1f) GeV/c, theta_cut = %.3f\n",lowPt,highPt,lowTheta);
	// MODE
//			TF1 * localFit = fitMass_Flat(hMassPtBinThetaBin[i][j],Form("MassFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta);
			TF1 * localFit = fitMass_Int_4(hMassPtBinThetaBin[i][j],Form("MassFit_Pt_%d_Theta_%d",i,j),lowPt,lowTheta,fClusterEnergyCut, nClusterPairs / hMassPtBinThetaBin[i][j]->GetXaxis()->GetBinWidth(4)); 
			// include pt bin size in normalization??

			localFits.push_back(localFit);
		}
		fMassFits.push_back(localFits);
	}


	//TCanvas * cPtTheta = new TCanvas("cPtTheta","cPtTheta",cWidth,cHeight);
	TCanvas * cMass = new TCanvas("cMass","cMass",cWidth,cHeight);
	TLegend * lMass = new TLegend(lX1,lY1,lX2,lY2);
	for (Int_t i = 0; i < nPtBins; i++) {
		// Draw no cut?
		for (Int_t j = 0; j < nUseThetaBins; j++) {
			hMassPtBinThetaBin[i][j]->GetXaxis()->SetRangeUser(kMinMassRange,kMaxMassRange);// 0.55
//		Don't skip first bin, ever		
			if (j>0 && iSkipAlternateBins > 0) {
				if ( (j + iSkipAlternateBins) % 2 == 1)  continue;
			}	
			if (!j) hMassPtBinThetaBin[i][j]->Draw();
			else hMassPtBinThetaBin[i][j]->Draw("SAME");
//			lMass->AddEntry(hMassPtBinThetaBin[i][j],Form("(%.3f #leq #theta_{c})",fLowTheta[j]),"lp");
			lMass->AddEntry(hMassPtBinThetaBin[i][j],Form("#theta_{c} = %.3f",fLowTheta[j]),"lp");
//			fMassFits[i][j]->Draw("SAME");
			lMass->Draw("SAME");
			//cMass->Print(Form("MassTheta_PtBin_%d_Frame_%2d.png",i,j));
		}
//		lMass->Draw("SAME");
//		cMass->SetLogy();
		cMass->Print(Form("MassTheta_PtBin_%d.pdf",i));
		cMass->Print(Form("MassTheta_PtBin_%d.C",i));
		cMass->Clear();
		lMass->Clear();
	}

	TCanvas * cMassRatio = new TCanvas("cMassRatio","cMassRatio",cWidth,cHeight);
	TLegend * lMassRatio = new TLegend(lX1,lY1,lX2,lY2);
	for (Int_t i = 0; i < nPtBins; i++) {
		cMassRatio->cd();
		for (Int_t j = 0; j < nUseThetaBins; j++) {
//			 hMassRatioPtBinThetaBin[i][j]->Draw();
			if (j>0 && iSkipAlternateBins > 0) {
				if ( (j + iSkipAlternateBins) % 2 == 1)  continue;
			}	
			if (j == 0) hMassRatioPtBinThetaBin[i][j]->Draw();
			else hMassRatioPtBinThetaBin[i][j]->Draw("SAME");
			// I had to disable this setting to use TEfficiency
			hMassRatioPtBinThetaBin[i][j]->GetXaxis()->SetRangeUser(0.0,0.5); 
//			hMassRatioPtBinThetaBin[i][j]->GetTotalHistogram()->GetXaxis()->SetRangeUser(0.0,0.5); 
			lMassRatio->AddEntry(hMassRatioPtBinThetaBin[i][j],Form("#theta_{c} = %.3f",fLowTheta[j]),"lp");
//			fIndivFile->Add(hMassRatioPtBinThetaBin[i][j]);

			//attempting to draw ratio fit
			fMassRatioFits[i][j]->Draw("SAME");
		}
		lMassRatio->Draw("SAME");

		cMassRatio->Print(Form("RatioThetaFit_PtBin_%d.pdf",i));
//		cMassRatio->Print(Form("RatioThetaFit_PtBin_%d.root",i));
		cMassRatio->Print(Form("RatioThetaFit_PtBin_%d.C",i));
		cMassRatio->Clear();
		lMassRatio->Clear();
	}
	
	for (Int_t i = 0; i < nPtBins; i++) {
		cMassRatio->cd();
		for (Int_t j = 0; j < nUseThetaBins; j++) {
//			 hMassRatioPtBinThetaBin[i][j]->Draw();
			TF1 * fLocalFit = 0;
			fLocalFit = (TF1 *) hMassRatioPtBinThetaBin[i][j]->GetFunction(Form("RatioPreFit_Pt_%d_Theta_%d",i,j));
			if (fLocalFit != 0) fLocalFit->SetBit(TF1::kNotDraw);
			else printf("i = %d j = %d: Could not find prefit function\n",i,j);

			if (j>0 && iSkipAlternateBins > 0) {
				if ( (j + iSkipAlternateBins) % 2 == 1)  continue;
			}	

			if (j == 0) hMassRatioPtBinThetaBin[i][j]->Draw("");
			else hMassRatioPtBinThetaBin[i][j]->Draw("SAME");
			//if (j == 0) hMassRatioPtBinThetaBin[i][j]->Draw("HIST E");
			//else hMassRatioPtBinThetaBin[i][j]->Draw("HIST E SAME");
			// I had to disable this setting to use TEfficiency
			hMassRatioPtBinThetaBin[i][j]->GetXaxis()->SetRangeUser(0.0,0.5); 
//			hMassRatioPtBinThetaBin[i][j]->GetTotalHistogram()->GetXaxis()->SetRangeUser(0.0,0.5); 
			lMassRatio->AddEntry(hMassRatioPtBinThetaBin[i][j],Form("#theta_{c} = %.3f",fLowTheta[j]),"lp");
		}
		lMassRatio->Draw("SAME");

		cMassRatio->Print(Form("RatioTheta_PtBin_%d.pdf",i));
//		cMassRatio->Print(Form("RatioTheta_PtBin_%d.root",i));
		cMassRatio->Print(Form("RatioTheta_PtBin_%d.C",i));
		cMassRatio->Clear();
		lMassRatio->Clear();
	}

	//Making Parameter Plots
	//Create arrays of parameters in bins of par,pt,theta
	Int_t nFitParams = fMassRatioFits[0][0]->GetNpar();
	for (Int_t i = 0; i < nFitParams; i++) {
		std::vector<std::vector<Double_t>> lFitParams_Param = {};
		std::vector<std::vector<Double_t>> lFitParamsUn_Param = {};
		for (Int_t j = 0; j < nPtBins; j++) {
			std::vector<Double_t> lFitParams_Pt = {};
			std::vector<Double_t> lFitParamsUn_Pt = {};
			for (Int_t k = 0; k < nUseThetaBins; k++) {
				Double_t lParam = fMassRatioFits[j][k]->GetParameter(i);
				Double_t lParamUn = fMassRatioFits[j][k]->GetParError(i);

				lFitParams_Pt.push_back(lParam);
				lFitParamsUn_Pt.push_back(lParamUn);
			}
			lFitParams_Param.push_back(lFitParams_Pt);
			lFitParamsUn_Param.push_back(lFitParamsUn_Pt);
		}
		fFitParams.push_back(lFitParams_Param);
		fFitParamsUn.push_back(lFitParamsUn_Param);
	}
	// Create TGraphs from Arrays
	// ChiSq Over NDF
	std::vector<TGraphErrors *> fChiSqGraphs = {};
	for (Int_t i = 0; i < nPtBins; i++) {
		TGraphErrors * lChiSqGraph = new TGraphErrors(nUseThetaBins,&fLowTheta[0],&fFitChiSqNDF[i][0],0,0);
		lChiSqGraph->SetName(Form("ChiSqNDF_PtBin_%d",i));
		lChiSqGraph->SetTitle(Form("#chi^{2} / NDF (%.0f #leq p_{T} < %.0f GeV/c)",fPtBins[i],fPtBins[i+1]));
		lChiSqGraph->SetMarkerStyle(kOpenSquare);
		lChiSqGraph->GetXaxis()->SetTitle("#theta Cut (Rad)");
		lChiSqGraph->GetYaxis()->SetTitle("#chi^{2} / NDF");
		fChiSqGraphs.push_back(lChiSqGraph);
	}
	


	TCanvas * cChiSq = new TCanvas("cChiSq","cChiSq",cWidth,cHeight);
	for (Int_t i = 0; i < nPtBins; i++) {
		fChiSqGraphs[i]->Draw("AP");
		
		cChiSq->Print(Form("ChiSqNDF_ptBin_%d.pdf",i));
		cChiSq->Print(Form("ChiSqNDF_ptBin_%d.C",i));
	}

	// FIXME did you make a mean theta array??
	std::vector<std::vector<TGraphErrors *>> fParamGraphs = {};
	for (Int_t i = 0; i < nFitParams; i++) {
		std::vector<TGraphErrors *> lParamGraphs_Pt = {};
		for (Int_t j = 0; j < nPtBins; j++) {
			//TGraphErrors * lParamGraph = new TGraphErrors(nUseThetaBins,&fMeanTheta[0],fFitParams[i][j],0,0);
			//TGraphErrors * lParamGraph = new TGraphErrors(nUseThetaBins,&fMeanTheta[0],&fFitParams[i][j][0],0,&fFitParamsUn[i][j][0]);
			TGraphErrors * lParamGraph = new TGraphErrors(nUseThetaBins,&fLowTheta[0],&fFitParams[i][j][0],0,&fFitParamsUn[i][j][0]);
			lParamGraph->SetName(Form("Param_%d_Pt_%d",i,j));	
			lParamGraph->SetMarkerStyle(kFullSquare);
			//lParamGraph->SetTitle(Form("Parameter %d Graph, (Pt Bin %d)",i,j));
			lParamGraph->SetTitle(Form("Parameter %d Graph, (%.1f #leq p_{T} < %.1f GeV/c)",i,fPtBins[j],fPtBins[j+1]));
			lParamGraph->GetXaxis()->SetTitle("#theta Cut (Rad)");
			if ( i == 0 ) lParamGraph->GetYaxis()->SetTitle("m' (MeV/c^{2})");
			else if ( i == 1 ) lParamGraph->GetYaxis()->SetTitle("#lambda (c^{2}/MeV)");
			else lParamGraph->GetYaxis()->SetTitle(Form("Parameter %d Value",i));

			lParamGraphs_Pt.push_back(lParamGraph);
		}
		fParamGraphs.push_back(lParamGraphs_Pt);
	}

	// Global Parameter Arrays
	// In pt bins
	std::vector<Double_t> fMassPrimePar0Array = {};
	std::vector<Double_t> fMassPrimePar1Array = {};
	std::vector<Double_t> fLambdaPar0Array    = {};
	std::vector<Double_t> fLambdaPar1Array    = {};

	std::vector<Double_t> fMassPrimePar0Array_Err = {};
	std::vector<Double_t> fMassPrimePar1Array_Err = {};
	std::vector<Double_t> fLambdaPar0Array_Err    = {};
	std::vector<Double_t> fLambdaPar1Array_Err    = {};
	// Could make this flexible 


	TGraphErrors * fMassPrimePar0Graph;
	TGraphErrors * fMassPrimePar1Graph;
	TGraphErrors * fLambdaPar0Graph;
	TGraphErrors * fLambdaPar1Graph;

	// Fitting the parameters as functions of theta, pt
	if (bGlobalParam) {
		printf("Starting the global parameter fits\n");
		// The default for this is just the two paramters (m' and lambda)
		// Maybe just do this the same way as before (param, pt)
		
		// Fitting m' vs theta_c 
		std::vector<TGraphErrors *> lMassPrimeGraphArray = fParamGraphs[0];
		for (Int_t i = 0; i < nPtBins; i++) {
			printf("Fitting MassPrime vs Theta for PtBin %.1f-%.1f ...\n",fPtBins[i],fPtBins[i+1]);
			TGraphErrors * lMassPrimeGraph = lMassPrimeGraphArray[i];
			TF1 * lFit = new TF1(Form("Par%d_PtBin%d",0,i),"([0]+[1]*x)*1000",0.0,0.05); // return to GeV

			lFit->FixParameter(0,0);

			lMassPrimeGraph->Fit(lFit);

			fMassPrimePar0Array.push_back(lFit->GetParameter(0));
			fMassPrimePar0Array_Err.push_back(lFit->GetParError(0));
			fMassPrimePar1Array.push_back(lFit->GetParameter(1));
			fMassPrimePar1Array_Err.push_back(lFit->GetParError(1));
		}

		// Fitting lambda vs theta_c
		std::vector<TGraphErrors *> lLambdaGraphArray = fParamGraphs[1];
		for (Int_t i = 0; i < nPtBins; i++) {
			printf("Fitting Lambda vs Theta for PtBin %.1f-%.1f ...\n",fPtBins[i],fPtBins[i+1]);
			TGraphErrors * lLambdaGraph = lLambdaGraphArray[i];
			TF1 * lFit = new TF1(Form("Par%d_PtBin%d",1,i),"([0]+[1]/(x+1e-6))/1000",0.001,0.05); //original
			//TF1 * lFit = new TF1(Form("Par%d_PtBin%d",1,i),"([0]/(x*x+1e-6) + [1]/(x+1e-6))/1000",0.001,0.05);

			lFit->FixParameter(0,0);

			lLambdaGraph->Fit(lFit);

			fLambdaPar0Array.push_back(lFit->GetParameter(0));
			fLambdaPar0Array_Err.push_back(lFit->GetParError(0));
			fLambdaPar1Array.push_back(lFit->GetParameter(1));
			fLambdaPar1Array_Err.push_back(lFit->GetParError(1));
		}

		// FIXME try making LambdaPar1Graph 2 dimensional? (theta,theta bin size)


		// FIXME replacing mid bins with normal bins (using the low pt value of each bin)
		fMassPrimePar0Graph = new TGraphErrors(nPtBins,&fPtBins[0],&fMassPrimePar0Array[0],0,&fMassPrimePar0Array_Err[0]);
		fMassPrimePar1Graph = new TGraphErrors(nPtBins,&fPtBins[0],&fMassPrimePar1Array[0],0,&fMassPrimePar1Array_Err[0]);
		fLambdaPar0Graph = new TGraphErrors(nPtBins,&fPtBins[0],&fLambdaPar0Array[0],0,&fLambdaPar0Array_Err[0]);
		fLambdaPar1Graph = new TGraphErrors(nPtBins,&fPtBins[0],&fLambdaPar1Array[0],0,&fLambdaPar1Array_Err[0]);

		// The old way
		//fMassPrimePar0Graph = new TGraphErrors(nPtBins,&fPtBinsMid[0],&fMassPrimePar0Array[0],0,&fMassPrimePar0Array_Err[0]);
		//fMassPrimePar1Graph = new TGraphErrors(nPtBins,&fPtBinsMid[0],&fMassPrimePar1Array[0],0,&fMassPrimePar1Array_Err[0]);
		//fLambdaPar0Graph = new TGraphErrors(nPtBins,&fPtBinsMid[0],&fLambdaPar0Array[0],0,&fLambdaPar0Array_Err[0]);
		//fLambdaPar1Graph = new TGraphErrors(nPtBins,&fPtBinsMid[0],&fLambdaPar1Array[0],0,&fLambdaPar1Array_Err[0]);

		fMassPrimePar0Graph->SetName("MassPrimePar0Graph");
		fMassPrimePar1Graph->SetName("MassPrimePar1Graph");
		fLambdaPar0Graph->SetName("LambdaPar0Graph");
		fLambdaPar1Graph->SetName("LambdaPar1Graph");

		fMassPrimePar0Graph->SetTitle("m' Par0");
		fMassPrimePar1Graph->SetTitle("m' Par1");
		fLambdaPar0Graph->SetTitle("#lambda Par0");
		fLambdaPar1Graph->SetTitle("#lambda Par1");

		fMassPrimePar0Graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		fMassPrimePar1Graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		fLambdaPar0Graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		fLambdaPar1Graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");

		fMassPrimePar0Graph->SetMarkerStyle(kOpenSquare);
		fMassPrimePar0Graph->SetMarkerColor(kBlack);
		fMassPrimePar0Graph->SetLineColor(kBlack);

		fMassPrimePar1Graph->SetMarkerStyle(kFullSquare);
		fMassPrimePar1Graph->SetMarkerColor(kBlack);
		fMassPrimePar1Graph->SetLineColor(kBlack);

		fLambdaPar0Graph->SetMarkerStyle(kOpenTriangleUp);
		fLambdaPar0Graph->SetMarkerColor(kBlack);
		fLambdaPar0Graph->SetLineColor(kBlack);

		fLambdaPar1Graph->SetMarkerStyle(kFullTriangleUp);
		fLambdaPar1Graph->SetMarkerColor(kBlack);
		fLambdaPar1Graph->SetLineColor(kBlack);

		//fMassPrimePar0Graph->GetYaxis()->SetTitle("C (100 Sr)");
		//fMassPrimePar1Graph->GetYaxis()->SetTitle("");


//		for (Int_t i = 0; i < nFitParams; i++) {
//			for (Int_t j = 0; j < nPtBins; j++) {
//			}
//		}
	}

	

	// Make Param vs Theta plots for each param, Pt Bin

	TCanvas * cParams = new TCanvas("cParams","cParams",cWidth,cHeight);
	for (Int_t i = 0; i < nFitParams; i++) {
		for (Int_t j = 0; j < nPtBins; j++) {
			fParamGraphs[i][j]->Draw("AP");

//			TPaveText *pt = new TPaveText(.75,.4,.95,.6);
			TPaveText *pt = new TPaveText(kPaveTextX1,kPaveTextY1,kPaveTextX2,kPaveTextY2,"NDC");

			TF1 * fLocalFit = fParamGraphs[i][j]->GetFunction(Form("Par%d_PtBin%d",i,j));
			if (fLocalFit) {
				printf("Preparing to add textbox with fit info\n");
				TFormula * formula = fLocalFit->GetFormula();
				Int_t nGlobalFitParams = fLocalFit->GetNpar();
				pt->AddText(Form("f(x) = %s",formula->GetExpFormula().Data()));
				pt->AddText(Form("p_{0} = %f #pm %f",fLocalFit->GetParameter(0),fLocalFit->GetParError(0)));
				pt->AddText(Form("p_{1} = %f #pm %f",fLocalFit->GetParameter(1),fLocalFit->GetParError(1)));
				pt->AddText(Form("#chi^{2}/NDF = %f/%d = %f",fLocalFit->GetChisquare(),fLocalFit->GetNDF(),fLocalFit->GetChisquare()/fLocalFit->GetNDF()));
			} else printf("Could not find associated function to a TGraphErrors\n");			

			pt->Draw("SAME");

			cParams->Print(Form("Param_%d_ptBin_%d.pdf",i,j));
			cParams->Print(Form("Param_%d_ptBin_%d.C",i,j));
			cParams->Clear();
		}
	}

	// Making Graph and Fit of meta-parameters
	//fMassPrimePar0Graph
	//fMassPrimePar1Graph
	//fLambdaPar0Graph
	//fLambdaPar1Graph

	TCanvas * cMetaParams = new TCanvas("MetaParams","MetaParams",cWidth,cHeight);
	//gStyle->SetOptStat(1);
	//gStyle->SetOptFit(111);

	fMassPrimePar0Graph->Draw("AP");

//	TF1 * fFitMassPrimePar0 = new TF1("MassPrime Par0 Fit","[0]",0,30);
//	fMassPrimePar0Graph->Fit(fFitMassPrimePar0);
	cMetaParams->Print("MP_MassPrime_Par0.pdf");
	cMetaParams->Print("MP_MassPrime_Par0.C");

	fMassPrimePar1Graph->Draw("AP");
	//TF1 * fTheoryMassPrimePar1 = new TF1("MassPrime_Par1_Prediction","[0]*x",0,30);
	//fTheoryMassPrimePar1->SetParameter(0,0.5); // m' = 0.5 pt * theta
	//TF1 * fTheoryMassPrimePar1 = new TF1("MassPrime_Par1_Theory","[0]*TMath::ATan(x/[1])",0,30); // Why did I do this?
/*	fTheoryMassPrimePar1->SetParameter(0,5.);
	fTheoryMassPrimePar1->SetParLimits(0,1.,20.);
	fTheoryMassPrimePar1->SetParameter(1,10.);
	fTheoryMassPrimePar1->SetParLimits(1,1.,20.);
	fMassPrimePar1Graph->Fit(fTheoryMassPrimePar1);
	fTheoryMassPrimePar1->SetLineColor(kViolet);
	fTheoryMassPrimePar1->Draw("SAME"); */
	//TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","[1]*TMath::Sqrt(2*[0]*(x-1-[0]))",0,30);
	//TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","TMath::Sqrt([0]*(x-[1]-[0]))",0,30);

	// Standard
	//TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","(1+[2]*x)*TMath::Sqrt([0]*(x-[1]-[0]))",0,30);
	// Trying this out
	TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","TMath::Sqrt(  ([0]*([2]*x*x + x-[1]-[0]) )*( [0]*([2]*x*x + x-[1]-[0])>0) )",0,30);


	fFitMassPrimePar1->SetParLimits(0,0.25,5);
	
	// FIXME for the test
	fFitMassPrimePar1->SetParLimits(2,0,5);

	// My theory
	//fFitMassPrimePar1->FixParameter(0,fClusterEnergyCut);
//	//fFitMassPrimePar1->FixParameter(0,2.); // this was fixed at 2.5 for some reason

	//TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","[0] + [1]*x",0,30);
	/*TF1 * fFitMassPrimePar1 = new TF1("MassPrime_Par1_Fit","[0]*TMath::ATan((x-[2])/([1]+0.1))",0,30);
	fFitMassPrimePar1->SetParLimits(0,1.,20.);
	fFitMassPrimePar1->SetParLimits(1,1,20.);
	fFitMassPrimePar1->SetParLimits(2,0,20.); */
	fMassPrimePar1Graph->Fit(fFitMassPrimePar1);
	fFitMassPrimePar1->SetLineColor(kFuncColor2);
	fFitMassPrimePar1->Draw("SAME");


	TPaveText *ptMPrimePar1 = new TPaveText(kPaveTextX1,kPaveTextY1,kPaveTextX2,kPaveTextY2,"NDC");
	TFormula * formulaMPrimePar1 = fFitMassPrimePar1->GetFormula();
	ptMPrimePar1->AddText(Form("f(x) = %s",formulaMPrimePar1->GetExpFormula().Data()));
	for (Int_t i = 0; i < fFitMassPrimePar1->GetNpar(); i++) {
		ptMPrimePar1->AddText(Form("p_{%d} = %f #pm %f",i,fFitMassPrimePar1->GetParameter(i),fFitMassPrimePar1->GetParError(i)));
	}
	ptMPrimePar1->AddText(Form("#chi^{2}/NDF = %f/%d = %f",fFitMassPrimePar1->GetChisquare(),fFitMassPrimePar1->GetNDF(),fFitMassPrimePar1->GetChisquare()/fFitMassPrimePar1->GetNDF()));
	ptMPrimePar1->Draw("SAME");


	cMetaParams->Print("MP_MassPrime_Par1.pdf");
	cMetaParams->Print("MP_MassPrime_Par1.C");

	fLambdaPar0Graph->Draw("AP");
//	TF1 * fFitLambdaPar0 = new TF1("Lambda_Par0_Fit","[0]",0,30);
//	fLambdaPar0Graph->Fit(fFitLambdaPar0);
	cMetaParams->Print("MP_Lambda_Par0.pdf");
	cMetaParams->Print("MP_Lambda_Par0.C");

	fLambdaPar1Graph->Draw("AP"); 
	// This fit is prone to floating pt errors
	//TF1 * fFitLambdaPar1 = new TF1("Lambda_Par1_Fit","[0] + [1] * TMath::Exp(-[2]*x)",0,30);
	// FIXME add a parameter for the pt bin size
	TF1 * fFitLambdaPar1 = new TF1("Lambda_Par1_Fit","[0] + [1] / TMath::Sqrt(TMath::Abs(x-[2]+1e-6))",4,30);
	fFitLambdaPar1->SetParameter(0,0.0); //0.6
	fFitLambdaPar1->SetParameter(1,4.4);
//	fFitLambdaPar1->SetParameter(2,0.3);
	fLambdaPar1Graph->Fit(fFitLambdaPar1);
	fFitLambdaPar1->SetLineColor(kFuncColor2);
	fFitLambdaPar1->Draw("SAME"); 

	TPaveText *ptLPar1 = new TPaveText(kPaveTextX1,kPaveTextY1,kPaveTextX2,kPaveTextY2,"NDC");
	TFormula * formulaLPar1 = fFitLambdaPar1->GetFormula();
	//Int_t nGlobalFitParams = fFitLambdaPar1->GetNpar();
	ptLPar1->AddText(Form("f(x) = %s",formulaLPar1->GetExpFormula().Data()));
	for (Int_t i = 0; i < fFitLambdaPar1->GetNpar(); i++) {
		ptLPar1->AddText(Form("p_{%d} = %f #pm %f",i,fFitLambdaPar1->GetParameter(i),fFitLambdaPar1->GetParError(i)));
	}
	//ptLPar1->AddText(Form("p_{1} = %f #pm %f",fFitLambdaPar1->GetParameter(1),fFitLambdaPar1->GetParError(1)));
	//ptLPar1->AddText(Form("p_{2} = %f #pm %f",fFitLambdaPar1->GetParameter(2),fFitLambdaPar1->GetParError(2)));
	ptLPar1->AddText(Form("#chi^{2}/NDF = %f/%d = %f",fFitLambdaPar1->GetChisquare(),fFitLambdaPar1->GetNDF(),fFitLambdaPar1->GetChisquare()/fFitLambdaPar1->GetNDF()));
	ptLPar1->Draw("SAME");

	cMetaParams->Print("MP_Lambda_Par1.pdf");
	cMetaParams->Print("MP_Lambda_Par1.C");

	// I should write this
	TF2 * fGlobalFit = FitGlobalModel();

	Double_t fLambda_Param00 = 0;
	Double_t fLambda_Param01 = 0;
	Double_t fLambda_Param10 = 1.2;
	Double_t fLambda_Param11 = -1./60.;

	Double_t fMassPrime_Param00 = 0;
	Double_t fMassPrime_Param01 = 0;
	Double_t fMassPrime_Param10 = 0;
	Double_t fMassPrime_Param11 = 0.5; //My Theory

	Double_t fDefaultThetaCut = 0.017;

	TF2 * fGlobalModel = CreatePtVsMassGlobalModel();

	fGlobalModel->SetParameter(0,fDefaultThetaCut);

	fGlobalModel->SetParameter(1,fLambda_Param00);
	fGlobalModel->SetParameter(2,fLambda_Param01);
	fGlobalModel->SetParameter(3,fLambda_Param10);
	fGlobalModel->SetParameter(4,fLambda_Param11);
	fGlobalModel->SetParameter(5,fMassPrime_Param00);
	fGlobalModel->SetParameter(6,fMassPrime_Param01);
	fGlobalModel->SetParameter(7,fMassPrime_Param10);
	fGlobalModel->SetParameter(8,fMassPrime_Param11);


	// Normalizing the raw histograms, to make it easy to compare in later files
	for (Int_t i = 0; i < nPtBins; i++) {
		NormalizeHistoHighMass(hMassPtBinNoThetaCut[i]);
//		Double_t scale = hMassPtBinNoThetaCut[i]->Integral("width");
//		if (scale != 0) hMassPtBinNoThetaCut[i]->Scale(1./scale);
		for (Int_t j = 0; j < nUseThetaBins; j++) {
//			scale = hMassPtBinThetaBin[i][j]->Integral("width");
//			if (scale != 0) hMassPtBinThetaBin[i][j]->Scale(1./scale);
			NormalizeHistoHighMass(hMassPtBinThetaBin[i][j]);
		}
	}

	// Saving Plots to an external file
	TFile * outFile = TFile::Open("AngleAnalysis.root","RECREATE");
	if (!outFile) return;
	// 1-D Plots
	
	// Ratios
	for (Int_t i = 0; i < nPtBins; i++) {
		outFile->Add(hMassPtBinNoThetaCut[i]);
		for (Int_t j = 0; j < nUseThetaBins; j++) {
			outFile->Add(hMassPtBinThetaBin[i][j]);
			outFile->Add(hMassRatioPtBinThetaBin[i][j]);
			outFile->Add(hMassEffPtBinThetaBin[i][j]);
		}
	}

	// ChiSq
	for (Int_t i = 0; i < nPtBins; i++) {
		outFile->Add(fChiSqGraphs[i]);
	}

	// Parameters
	for (Int_t i = 0; i < nFitParams; i++) {
		for (Int_t j = 0; j < nPtBins; j++) {
			outFile->Add(fParamGraphs[i][j]);
		}
	}

	if (bGlobalParam) {
		outFile->Add(fMassPrimePar0Graph);
		outFile->Add(fMassPrimePar1Graph);
		outFile->Add(fLambdaPar0Graph);
		outFile->Add(fLambdaPar1Graph);
		outFile->Add(fFitMassPrimePar1); 
		//outFile->Add(fTheoryMassPrimePar1); 
		outFile->Add(fFitLambdaPar1); 
		//outFile->Add(fTheoryLambdaPar1); 

		// FIXME also save the direct values

		// Writing out the final parameters in a copy and pasteable format
		printf("FinalParams:\n");
		for (int i = 0; i < 3; i++) {
			printf("\t\t\t\tpLambda_1%d = %f; // \\pm %f\n",i,fFitLambdaPar1->GetParameter(i),fFitLambdaPar1->GetParError(i));
		}
		for (int i = 0; i < 3; i++) {
			printf("\t\t\t\tpMPrime_1%d = %f; // \\pm %f\n",i,fFitMassPrimePar1->GetParameter(i),fFitMassPrimePar1->GetParError(i));
		}

		if (fGlobalModel != 0) outFile->Add(fGlobalModel);
	}

	outFile->Write();

}

// One part of the I function, plug in the limit u and factor z
// Subtract in a higher function
// May have issue with infinity
Double_t FuncI_Half_4(Double_t z, Double_t u) {
	Double_t xi = 1./TMath::Sqrt(4*z + 1);

	Double_t value =(1./6.) * (1 - u) * ( 2*u*u - 7*u + 6*z + 11);
//	printf("u^2 - u - z = %f\t\t\n",u*u-u-z);
//	if (u*u - u - z >= kTol) {
//		value += (z+0.5) * TMath::Log( u*u - u - z);
//		printf("good\n");
//	} 
	if (TMath::Abs(u*u - u - z) >= kTol) {
		value += (z+0.5) * TMath::Log(TMath::Abs( u*u - u - z));
	} else return 0;
	value += -xi*(2*z*z + 4*z + 1)*TMath::ATanH(xi*(1-2*u) * (abs(xi*(1-2*u)) < 1. - kTol) );

	return value;
}

Double_t FuncI_Full_4(Double_t z, Double_t u1, Double_t u2) {
	Double_t xi = 1./TMath::Sqrt(4*z + 1);
	printf("FuncI 4 : z = %f, xi = %f , u1 = %f, u2 = %f\n",z,xi,u1,u2);

	Double_t value =  (1./6.) * (1 - u1) * ( 2*u1*u1 - 7*u1 + 6*z + 11) - (1./6.) * (1 - u2) * ( 2*u2*u2 - 7*u2 + 6*z + 11);
	
	printf("FullVal 1 = %f\n",value);

	double logInput = (z + u1 - u1*u1)/(z + u2 - u2*u2);

	printf("  logInput 1 = %f\n",logInput);

	if (logInput < kTol)  return 0;
	value += (z+0.5) * TMath::Log(abs(logInput));
	printf("FullVal 2 = %f\n",value);

	logInput = ((1 + xi*(1 - 2*u1)) * (1 - xi*(1 - 2*u2)))/((1 - xi*(1 - 2*u1)) * (1 + xi*(1 - 2*u2)));

	printf("  logInput 2 = %f\n",logInput);
	if (logInput < kTol)  return 0;
	value += -xi * (2*z*z + 4*z + 1)* 0.5 * TMath::Log(logInput);

	printf("FuncInt_Full_4 returning %e\n",value);
	if (value < 0) return 0;
		
	//FIXME test
	value = TMath::Power(value,1./(2.*4.-1.));

	return value;	
}


// x = mass
// pars 
// 0:theta cut
// 1:pt
// 2:A (scale. Should set to nClusterPairs)
// 3:E_Cut
Double_t FuncInt_4(Double_t * x, Double_t *p) {

	Int_t B = 4;

	double E_Cut = p[3];

	double theta_cut = p[0];
	double pt = p[1];
	double u1 = TMath::Cos(theta_cut);
	double u2 = 0.5*(x[0]*x[0]*E_Cut*E_Cut+2*pt*pt*E_Cut*E_Cut-2.*TMath::Power(E_Cut,4)-x[0]*x[0]*E_Cut*TMath::Sqrt(x[0]*x[0] + pt*pt))/(x[0]*x[0]*E_Cut*E_Cut + pt*pt*E_Cut*E_Cut - TMath::Power(E_Cut,4)); 

	if ( u2 >= u1) return 0;

/*	Double_t value = 2*p[2]*TMath::Power(E_Cut/x[0],2*B-1) * TMath::Power(B-1,-2) * TMath::Power(2,-B);
	printf("   1 Value = %e\n",value);
	value = value / pt; // * TMath::Power(x[0],2*B-1));
	////value = value / (pt * TMath::Power(x[0],2*B-1));
	printf("   2 Value = %e\n",value);
//	value = value * (FuncI_Half_4(x[0]*x[0]/(pt*pt),u1) - FuncI_Half_4(x[0]*x[0]/(pt*pt),u2));
	value = value * FuncI_Full_4(x[0]*x[0]/(pt*pt),u1,u2);
*/
/*
	Double_t value =(1./6.) * (1 - u) * ( 2*u*u - 7*u + 6*z + 11);
	if (TMath::Abs(u*u - u - z) >= kTol) {
		value += (z+0.5) * TMath::Log(TMath::Abs( u*u - u - z));
	} else return 0;
	value += -xi*(2*z*z + 4*z + 1)*TMath::ATanH(xi*(1-2*u) * (abs(xi*(1-2*u)) < 1. - kTol) );
*/
	// FIXME test
	Double_t value = 2*p[2]* TMath::Power(B-1,-2) * TMath::Power(2,-B);
	printf("   1 Value = %e\n",value);
	value = value / pt; // * TMath::Power(x[0],2*B-1));
	////value = value / (pt * TMath::Power(x[0],2*B-1));
	printf("   2 Value = %e\n",value);
//	value = value * (FuncI_Half_4(x[0]*x[0]/(pt*pt),u1) - FuncI_Half_4(x[0]*x[0]/(pt*pt),u2));
	value = value * TMath::Power((E_Cut/x[0])*FuncI_Full_4(x[0]*x[0]/(pt*pt),u1,u2),2*B-1);



	printf("Returning Value = %e\n",value);

	return value; 
}

TF1 * fitMass_Int_4(TH1F * hist, std::string name, double pT = 5, double theta = 0.03, double ECut = 0.5, double norm = 1) {
	TF1 * fit = new TF1(name.c_str(),FuncInt_4,0.010,hist->GetXaxis()->GetXmax(),4);
	fit->SetNpx(250);

	fit->SetLineColor(hist->GetLineColor());

	fit->SetParName(0,"#theta_{cut}");
	fit->SetParName(1,"p_{T}");
	fit->SetParName(2,"A");
	fit->SetParName(3,"E_{c}");

	fit->SetParameter(0,theta);
	fit->SetParameter(1,pT);
	fit->SetParameter(2,1);
	fit->SetParameter(3,ECut);

	bool fixNorm = true; 

	if (fixNorm) {
		fit->FixParameter(2,norm);
	} else {
		float normRangeMin = 0.3-.15;
		float normRangeCenter = 0.35-.15;
		float normRangeMax = 0.4-.15;

		double integral = hist->Integral(hist->FindBin(normRangeMin),hist->FindBin(normRangeMax),"");
		double length = hist->FindBin(normRangeMax) - hist->FindBin(normRangeMin);
		double average = integral / length;	
		fit->SetParameter(2,1);
		double temp = fit->Eval(normRangeCenter);

		fit->SetParameter(2,1.*average/temp);
		fit->SetParLimits(2,0.0,5*average/temp);
	}

	return fit;
}



TF1 * fitMass_Flat(TH1F * hist, std::string name, double pT = 5, double theta = 0.03) {

	TString formuoli = "[2] *(x/TMath::Sqrt(4*x*x+[1]*[1]))*(TMath::ATanH((2.*TMath::Cos([0])-1.)/TMath::Sqrt(4.*x*x/([1]*[1])+1)*(0.999999>(2.*TMath::Cos([0])-1.)/TMath::Sqrt(4.*x*x/([1]*[1])+1))) - TMath::ATanH(-1./TMath::Sqrt(4.*x*x/([1]*[1])+1)*(0.999999>(-1.)/TMath::Sqrt(4.*x*x/([1]*[1])+1))) ) ";

	TF1 * fit = new TF1(name.c_str(),formuoli,0.010,hist->GetXaxis()->GetXmax());
	//TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

	fit->SetLineColor(hist->GetLineColor());

	fit->SetParName(0,"#theta_{C}");
	fit->SetParName(1,"p");

	fit->SetParameter(0,theta);
	//fit->SetParameter(0,0.);
	fit->SetParameter(1,pT);

	fit->SetParLimits(0,0,1.57);
	fit->SetParLimits(1,pT,30);
		
	fit->FixParameter(0,fit->GetParameter(0));
	fit->FixParameter(1,fit->GetParameter(1));

	float normRangeMin = 0.3;
	float normRangeCenter = 0.35;
	float normRangeMax = 0.4;

	double integral = hist->Integral(hist->FindBin(normRangeMin),hist->FindBin(normRangeMax),"");
	double length = hist->FindBin(normRangeMax) - hist->FindBin(normRangeMin);
	double average = integral / length;	
	fit->SetParameter(2,1);
	double temp = fit->Eval(normRangeCenter);

	fit->SetParName(2,"A");
	fit->SetParameter(2,1*average/temp);
	fit->SetParLimits(2,0.0,5*average/temp);

//	hist->Fit(fit,"");
	return fit;
}


// June 2019
TF1 * fitRatio_4(TH1F * hRatio, std::string name, double pT = 5, double theta = 0.03) {

	// Inside Sigmoid, want function:
	//   1. Monotonic
	//   2. Zero of function is a parameter (this zero is the inflection point of the sigmoid)
	//   3. Not fully antisymmetric (the data does not indicate a symmetric sigmoid

//	TString formuoli = "TMath::Erf([0]+[1]*x+[2]*TMath::Log(x+1e-3))";
	TString LogArgument = "( 100*x - [0] )* ( 100*x - [0] > 0 ) + 1e-5";
	TString SigmoidArgument = Form("[1]*TMath::Log(%s)+[2]*100*x + [3]",LogArgument.Data());;
	TString formuoli = Form("TMath::Erf(%s)",SigmoidArgument.Data());

	TF1 * fit = new TF1(name.c_str(),formuoli,0.010,hRatio->GetXaxis()->GetXmax());
	//TF1 * fit = new TF1(name.c_str(),formuoli,hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());

	fit->SetLineColor(hRatio->GetLineColor());

//	fit->SetParameter(0,-31.5);
	//fit->SetParameter(0,0.05);
	fit->SetParameter(0,5.);
//	fit->SetParameter(1,100.);
	fit->SetParameter(1,1);
	fit->SetParameter(2,3);
	fit->SetParameter(3,0.5);

	fit->SetParLimits(0,1,40);
//	fit->SetParLimits(0,0.01,0.015);
	fit->SetParLimits(1,0.,10.);
	fit->SetParLimits(2,0.,10.);
	fit->SetParLimits(3,0.,10.);

	hRatio->Fit(fit,"");
	return fit;
}
/**
 * Fit the ratio as a TEfficiency object, using the formula derived in fitRatio_3
 *
 */
TF1 * fitEfficiency(TEfficiency * hRatio, std::string name, double pT = 5, double theta = 0.03) {

	if (theta == 0) {
		fprintf(stderr,"Theta = 0 here. This should not happen.\n");
	}

	// Inside Sigmoid, want function:
	//   1. Monotonic
	//   2. Zero of function is a parameter (this zero is the inflection point of the sigmoid)
	//   3. Not fully antisymmetric (the data does not indicate a symmetric sigmoid

	TString SigmoidArgument = "[1]*(1e3*x - [0])";
	//TString formuoli = Form("0.5 + 0.5 * TMath::Erf(%s)",SigmoidArgument.Data()); 
	TString formuoli = Form("TMath::Erf(%s) * (%s > 0)",SigmoidArgument.Data(),SigmoidArgument.Data()); // FIXME should this be 0.5 + 0.5 erf(%s)?

	TF1 * fit = new TF1(name.c_str(),formuoli,0.010,hRatio->GetTotalHistogram()->GetXaxis()->GetXmax());
	//TF1 * fit = new TF1(name.c_str(),formuoli,hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());

	fit->SetLineColor(hRatio->GetLineColor());


	//Double_t Param0_Guess = 3000 * theta;
	Double_t Param0_Guess = 300 * theta;
	Double_t Param1_Guess = 0.009 / theta;

	//fit->SetParLimits(0,0,25);
	//Double_t fMinPar0 = theta * 100;
	//Double_t fMaxPar0 = 30. + theta*100;
	//Double_t fMinPar0 = theta * 1000;
	Double_t fMinPar0 = 0;
	Double_t fMaxPar0 = 300. + theta*1000;

	Double_t fStepSizePar0 = 0.0001; // MeV
	Double_t fStepSizePar1 = 0.000001; // 1/MeV	

	fit->SetParError(0,fStepSizePar0);
	fit->SetParError(1,fStepSizePar1);

	fit->SetParLimits(0,fMinPar0,fMaxPar0);
	fit->SetParLimits(1,0.,2.);

	fit->SetParameter(0,Param0_Guess); 
	fit->SetParameter(1,Param1_Guess);


//	fit->SetParameter(2,3);
//	fit->SetParameter(3,0.5);


	printf("Setting Par0 Limits %f - %f\n",fMinPar0,fMaxPar0);
	printf("Pre fit parameters: Par0 = %f MeV\n",fit->GetParameter(0));
	printf("                    Par1 = %f 1/MeV\n",fit->GetParameter(1));


//	fit->SetParLimits(0,0.01,0.015);
	//fit->SetParLimits(1,0.,5.);


//	fit->SetParLimits(2,0.,10.);
//	fit->SetParLimits(3,0.,10.);

//	fit->FixParameter(2,0.);

	//hRatio->Fit(fit,"R",0,0.4);
	hRatio->Fit(fit,"");

	printf("Final Param Results:Par0 = %f MeV\n",fit->GetParameter(0));
	printf("                    Par1 = %f 1/MeV\n",fit->GetParameter(1));

	return fit;

}

// June 2019
TF1 * fitRatio_3(TH1F * hRatio, std::string name, double pT = 5, double theta = 0.03) {

	// Inside Sigmoid, want function:
	//   1. Monotonic
	//   2. Zero of function is a parameter (this zero is the inflection point of the sigmoid)
	//   3. Not fully antisymmetric (the data does not indicate a symmetric sigmoid

//	TString formuoli = "TMath::Erf([0]+[1]*x+[2]*TMath::Log(x+1e-3))";
//	TString LogArgument = "( 100*x - [0] )* ( 100*x - [0] > 0 ) + 1e-5";
//	TString SigmoidArgument = Form("[1]*TMath::Log(%s)+[2]*100*x + [3]",LogArgument.Data());;
	//TString SigmoidArgument = "[1]*(100*x - [0]) + [2] * TMath::Log(100*x + 1e-3)";
	//TString SigmoidArgument = "[1]*(100*x - [0])";
	//TString SigmoidArgument = "[1]*(100*x - [0])";
	TString SigmoidArgument = "[1]*(1e3*x - [0])";
	//TString formuoli = Form("0.5 + 0.5 * TMath::Erf(%s)",SigmoidArgument.Data()); 
	TString formuoli = Form("TMath::Erf(%s) * (%s > 0)",SigmoidArgument.Data(),SigmoidArgument.Data()); // FIXME should this be 0.5 + 0.5 erf(%s)?

	TF1 * fit = new TF1(name.c_str(),formuoli,0.010,hRatio->GetXaxis()->GetXmax());
	//TF1 * fit = new TF1(name.c_str(),formuoli,hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());



	fit->SetLineColor(hRatio->GetLineColor());

	//Double_t Param0_Guess = 3000 * theta;
	// Second version:
	//Double_t Param0_Guess = 1000 * theta * (3./8.) * pT;
	//Double_t Param1_Guess = 0.009 / theta;
	// Based on the final result
	Double_t Param0_Guess = 1000 * (9.072 * TMath::ATan(pT / 20.) * theta); // MPrime model
	Double_t Param1_Guess = ((-2.574 + 4.314 * TMath::Exp(-0.02344 * pT)) / theta) / 1000; //Lambda model
	//fit->SetParLimits(0,0,25);
	//Double_t fMinPar0 = theta * 100;
	//Double_t fMaxPar0 = 30. + theta*100;
	//Double_t fMinPar0 = theta * 1000;

//	fit->SetParameter(0,-31.5);
	//fit->SetParameter(0,0.05);
	//fit->SetParameter(0,300 * theta); 
//	fit->SetParameter(0,300 * theta); 
//	fit->SetParameter(1,100.);
	//fit->SetParameter(1,1);
//	fit->SetParameter(1,0.1);

	//fit->SetParLimits(0,0,25);
	Double_t fMinPar0 = theta * 100;
	//Double_t fMaxPar0 = 30. + theta*100;
//	Double_t fMinPar0 = theta * 1000;
	//Double_t fMinPar0 = 0;
	Double_t fMaxPar0 = 300. + theta*1000;

	fit->SetParLimits(0,fMinPar0,fMaxPar0);
	fit->SetParLimits(1,0.,2.);

	fit->SetParameter(0,Param0_Guess); 
	fit->SetParameter(1,Param1_Guess);

	printf("Setting Par0 Limits %f - %f\n",fMinPar0,fMaxPar0);
	printf("Pre fit parameters: Par0 = %f MeV\n",fit->GetParameter(0));
	printf("                    Par1 = %f 1/MeV\n",fit->GetParameter(1));

	fit->SetParLimits(0,fMinPar0,fMaxPar0);

//	fit->SetParLimits(0,0.01,0.015);
	//fit->SetParLimits(1,0.,5.);
	fit->SetParLimits(1,0.,0.5);


	fit->SetParLimits(2,0.,10.);
//	fit->SetParLimits(3,0.,10.);

	fit->FixParameter(2,0.);

	printf("Histogram has integral %f and NEntries %f\n",hRatio->Integral(),hRatio->GetEntries());
	if (hRatio->GetEntries() == 0) {
		printf("Empty histogram, returning unfit function, skipping fit.\n");
		return fit;
	}

	hRatio->Fit(fit,"R",0,0.4);
	printf("Final Param Results:Par0 = %f MeV\n",fit->GetParameter(0));
	printf("                    Par1 = %f 1/MeV\n",fit->GetParameter(1));
	return fit;
}


TF1 * fitRatio_2(TH1F * hRatio, std::string name, double pT = 5, double theta = 0.03) {
//	TString formuoli = "([0]*x - TMath::ATanH((-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))*(-0.9<-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))))";
//	formuoli        += "/(-TMath::ATanH((-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))*(-0.9<-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))))";

	TString formuoli = "(TMath::ATanH((2.*TMath::Cos([0])-1.)/TMath::Sqrt(4.*x*x/([1]*[1])+1)*(0.999999>(2.*TMath::Cos([0])-1.)/TMath::Sqrt(4.*x*x/([1]*[1])+1))) - TMath::ATanH(-1./TMath::Sqrt(4.*x*x/([1]*[1])+1)*(-0.999999<-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))))";
	formuoli        += "/(0.000001+TMath::ATanH(1./TMath::Sqrt(4.*x*x/([1]*[1])+1)*(0.999999>1./TMath::Sqrt(4.*x*x/([1]*[1])+1)))- TMath::ATanH(-1./TMath::Sqrt(4.*x*x/([1]*[1])+1)*(-0.999999<-1./TMath::Sqrt(4.*x*x/([1]*[1])+1))))";

	TF1 * fit = new TF1(name.c_str(),formuoli,0.010,hRatio->GetXaxis()->GetXmax());
	//TF1 * fit = new TF1(name.c_str(),formuoli,hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());

//	printf("Hello\n");

	fit->SetLineColor(hRatio->GetLineColor());

	fit->SetParName(0,"#theta_{C}");
	fit->SetParName(1,"p");

	fit->SetParameter(0,theta);
	fit->SetParameter(1,pT);

	fit->SetParLimits(0,0,1.57);
	fit->SetParLimits(1,pT,30);
		
//	printf("Val (0.1) = %f\n",fit->Eval(0.1));

//	hRatio->Fit(fit,"Q");
	return fit;
}


TF1 * fitRatio(TH1F * hRatio, std::string name) {
//	TString formuoli = "1./(1.+ TMath::Exp(-[0]*(x-[1])))";
//	TString formuoli = "[0] * x + [1] + [2]*x*x+[3]*x*x*x+[4]*x*x*x*x";
//	TString formuoli = "TMath::Power(0.5 + TMath::ATan([0]*(x-[1]))/TMath::Pi(),2)";
//	TString formuoli = "(0.5 + TMath::ATan([0]/(1-(x-[1])))/TMath::Pi())";

	TString formuoli = "1./(1.+TMath::Exp(-[0]*(x-[1])))";
	TF1 * fit = new TF1(name.c_str(),formuoli,hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());

	fit->SetLineColor(hRatio->GetLineColor());

	fit->SetParameter(0,1);
	fit->SetParameter(1,0.12);

	fit->SetParLimits(0,0,450);
	fit->SetParLimits(1,0,0.5);
		
//	fit->SetParameter(2,1);
//	fit->SetParLimits(2,0,4);


//	printf("Debug: This histo has %f\n",hRatio->GetEntries());

	hRatio->Fit(fit,"Q");
	return fit;
}

// Scale so the region to the right of the eta peak averages to unity
void NormalizeHistoHighMass(TH1 * histo) {
	double NormMin = 0.6;
	double NormMax = histo->GetXaxis()->GetXmax();
	double scale = histo->Integral(histo->GetXaxis()->FindBin(NormMin),histo->GetXaxis()->GetNbins(),"width");
	double integral = histo->Integral("width");
	if (scale > 0.2 * integral) {
		histo->Scale((NormMax - NormMin)/scale);
	} else {
		printf("using alternate scaling\n");
		// if no stats in norm region, use entire histogram
//		scale = histo->Integral("width");
		if (integral != 0.) { 
			histo->Scale(0.1*NormMax/integral);
		}
	}
	
}

// Fit the global parameters
TF2 * FitGlobalModel() {
	TString sGlobalModel = "";
	TF2 * fGlobalModel = 0;

	return fGlobalModel;
}

TF2 * CreatePtVsMassGlobalModel() {
	TF2 * fGlobalModel = 0;
	
	// Note: parameters may need to be rescaled, as here we only use GeV

	// The original, with lambda and mass'
//	TString SigmoidArgument = "[1]*(1e3*x - [0])";

	// Param 0 is the theta cut

	// Linear (?) Fit for Lambda
	TString sLambda_Par0 = "[1] + [2]*y"; 
	TString sLambda_Par1 = "[3] + [4]*y";

	TString sMassPrime_Par0 = "[5] + [6]*y";
	TString sMassPrime_Par1 = "[7] + [8]*y";
	
	// Linear Fit for MassPrime

	TString sLambda = Form("(%s) + (%s)/(TMath::Abs([0])+1e-6)",sLambda_Par0.Data(),sLambda_Par1.Data()); 
	TString sMassPrime = Form("(%s) + (%s)*[0]",sMassPrime_Par0.Data(),sMassPrime_Par1.Data());


	//TString SigmoidArgument = Form("(%s)*(1e3*x - (%s))",sLambda.Data(),sMassPrime.Data());
	TString SigmoidArgument = Form("(%s)*(x - (%s))",sLambda.Data(),sMassPrime.Data());

	printf("Building the global model with Sigmoid Argument: %s\n",SigmoidArgument.Data());

	TString sGlobalModel = Form("TMath::Erf(%s) * (%s > 0)",SigmoidArgument.Data(),SigmoidArgument.Data());
	fGlobalModel = new TF2("GlobalModel",sGlobalModel,0,1,5,30);

	fGlobalModel->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
	fGlobalModel->GetYaxis()->SetTitle("p_{T} (GeV/c)");

	return fGlobalModel;
}


