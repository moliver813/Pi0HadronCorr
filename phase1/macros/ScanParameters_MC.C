
TString sTrigger = "MC";

void SetStyle() {

	gStyle->SetOptStat(0);
}

void ScanParametersCent(int iCentBin, int iPtbin);

void ScanParameters_MC(int iPtBin = 1) { 
	
//	if (bTrigger == 1) {
//		sTrigger = "GA";
//		iPtBin = 4;
//	} 
	printf("Beginning Cut Parameter Scan for %s Triggers\n",sTrigger.Data());
	vector<int> fCentBins = {0,1,2,3};
	for (int i : fCentBins) {
		ScanParametersCent(i,iPtBin);
	}
}

void ScanParametersCent(int iCentBin, int iPtBin) {
	printf("Beginning Parameter Scan for Lambda and Theta, Cent Bin %d\n",iCentBin);
	SetStyle();

	vector<int> fLambdaBins = {1, 2, 3, 4, 5, 6, 7, 8};
	vector<int> fThetaBins  = {3, 4, 5, 6, 7};
	int nLambdaBins = (int) fLambdaBins.size();
	int nThetaBins = (int) fThetaBins.size();

	int kChiSqMaskColor = kRed-4;

	TString sOutputFileName = Form("ParamScan_%s/Cent%d/ParamScan.root",sTrigger.Data(),iCentBin);
	TString sOutputDir = Form("ParamScan_%s/Cent%d",sTrigger.Data(),iCentBin);

	// Load the files
	TString sFileFormat = "output/T38/ThetaLambdaScan/Lambda0%d/Theta0%d/Cent%d/EP-1/SecondAnalysis.root";
	vector<vector<TFile *>> fParamScanFiles = {};
	for (int iLambdaBin : fLambdaBins) {
		vector<TFile * > fLocalList = {};
		for (int iThetaBin : fThetaBins) {
			printf("Loading file for Lambda %d and Theta %d\n",iLambdaBin,iThetaBin);
			TString sParamScanFileName = Form(sFileFormat.Data(),iLambdaBin,iThetaBin,iCentBin);
			printf("  filename: %s\n",sParamScanFileName.Data());
			TFile * fParamScanFile = new TFile(sParamScanFileName,"OPEN");
			fLocalList.push_back(fParamScanFile);
		}
		fParamScanFiles.push_back(fLocalList);
	}


	int iPickPtBin = iPtBin;
	//int iPickPtBin = 1; // could iterate this
	// For some graphs, we want to pick a specific pt bin to look act
	//vector<TString> sPickGraphs = {"Pi0ChiSquare", "Pi0YieldTotalRatio", "Pi0PeakSigRatio", "Pi0IntYield"};
	vector<TString> sPickGraphs = {"Pi0ChiSquare", "Pi0YieldTotalRatio", "Pi0PeakSigRatio", "Pi0IntYield", "RecMCPi0YieldRatioArr"};


	// histogram for highlighing the bins that pass a chisquared test
	float kChiSqCut = 1.75; //2.0;

	TH2D * hChiSqMask = new TH2D(Form("ChiSqMask_PtBin%d_Cent%d",iPickPtBin,iCentBin),Form("#Chi^{2} Test : PtBin %d : CentBin %d;#lambda_{0}^{2} Cut Bin;#theta_{c} Cut Bin",iPickPtBin,iCentBin),nLambdaBins,(float) fLambdaBins[0]-0.5,(float) fLambdaBins[nLambdaBins-1]+0.5,nThetaBins,(float) fThetaBins[0] - 0.5,(float) 0.5 + fThetaBins[nThetaBins-1]);
	hChiSqMask->SetLineColor(kChiSqMaskColor);

	// vector<vector<vector<TH1D *>>> fPickGraphs; 
	vector<TH2D *> fPickHistograms = {};

	for (TString sPickGraph : sPickGraphs) {
		printf("  Starting the pick graph process for %s\n",sPickGraph.Data());
		
		TH2D * hLambdaThetaPickValues = new TH2D(Form("PG_%s_PtBin%d_Cent%d",sPickGraph.Data(),iPickPtBin,iCentBin),Form("Parameter Array : %s : PtBin %d : CentBin %d;#lambda_{0}^{2} Cut Bin;#theta_{c} Cut Bin",sPickGraph.Data(),iPickPtBin,iCentBin),nLambdaBins,(float) fLambdaBins[0]-0.5,(float) fLambdaBins[nLambdaBins-1]+0.5,nThetaBins,(float) fThetaBins[0] - 0.5,(float) 0.5 + fThetaBins[nThetaBins-1]);
		for (int iLambda = 0; iLambda < nLambdaBins; iLambda++) {
			for (int iTheta = 0; iTheta < nThetaBins; iTheta++) {
				TGraphErrors * fPickGraph = (TGraphErrors *) fParamScanFiles[iLambda][iTheta]->Get(sPickGraph.Data());
				double fLocalValue = 0;
				double fLocalError = 0;
				if (!fPickGraph) {
					fprintf(stderr,"Missing graph %s for lambda %d and theta %d\n",sPickGraph.Data(),iLambda,iTheta);
					continue;
				}	// could use else and set -1 or something
				fLocalValue = fPickGraph->GetY()[iPickPtBin];
				fLocalError = fPickGraph->GetEY()[iPickPtBin];
//				hLambdaThetaPickValues->SetBinContent(iLambda,iTheta,fLocalValue);
//				hLambdaThetaPickValues->SetBinError(iLambda,iTheta,fLocalError);
				hLambdaThetaPickValues->SetBinContent(iLambda+1,iTheta+1,fLocalValue);
				hLambdaThetaPickValues->SetBinError(iLambda+1,iTheta+1,fLocalError);

				if (sPickGraph.EqualTo("Pi0ChiSquare")) {
					if (fLocalValue < kChiSqCut) 
						hChiSqMask->SetBinContent(iLambda+1,iTheta+1,1.0);
				}
			}
		}
		fPickHistograms.push_back(hLambdaThetaPickValues);
	}

	TCanvas * cCanvas = new TCanvas("Canvas","Canvas",800,600);
	for (auto hist : fPickHistograms) {
		cCanvas->Clear();
		hist->Draw("COLZ TEXT");
//		hist->Draw("COLZ ARR");

		hChiSqMask->Draw("SAME0 BOX");
		cCanvas->Print(Form("%s/%s.pdf",sOutputDir.Data(),hist->GetName()));
		cCanvas->Print(Form("%s/%s.png",sOutputDir.Data(),hist->GetName()));
		cCanvas->Print(Form("%s/%s.C",sOutputDir.Data(),hist->GetName()));
	}

	TFile * fOutputFile = new TFile(sOutputFileName.Data(),"RECREATE");
	for (auto hist : fPickHistograms) fOutputFile->Add(hist);
	fOutputFile->Write();



}
