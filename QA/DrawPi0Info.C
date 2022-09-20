
const Int_t nTypeList = 15;

Int_t colorList[nTypeList] = {
	kBlack,kRed,kOrange-3,
	kGreen-3,kBlue,kViolet,
	kMagenta-5,kRed-10,kGray,
	kCyan,kAzure,kAzure-8,
	kGreen-3,kSpring+9,kYellow+1
};
Int_t markerList[nTypeList] = {
	kFullSquare,kFullCircle,kFullDiamond,
	kOpenSquare,kOpenCircle,kOpenDiamond,
	kFullStar,23,34,
	39,41,43,
	40,42,44
};

const int kNCentBins = 4;
double fCentBins[5] = {0,10,30,50,90};
vector<TString> fCentBinTitles = {"0-10%%","10-30%%","30-50%%","50-90%%"};

//const int kNTrackPtBins = 8;
//std::vector <double> fTrackPtBins ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};
const int kNTrackPtBins = 5;
std::vector <double> fTrackPtBins ={0.2,0.4,0.8,1.5,2.5,4};

const int kNPi0PtBins = 5;
std::vector<double> fTriggerPtBins = {5,7,9,11,14,17};






void SetHistStyle(TH1 * hist, int i) {
	if (i < 0 || i >= nTypeList) return;
	if (hist != 0) {
		hist->SetLineColor(colorList[i]);
		hist->SetMarkerColor(colorList[i]);
		hist->SetMarkerStyle(markerList[i]);
	}
}

// Normalize to input scale
// if normalizing by number of events, input the number of events as the rescale
TH1F * NormalizeTH1(TH1F * hHist, double rescale = 1) {
	TString sNormName = "Norm";
	//if (rescale == 1.0) sNormName = "NormUnity";
	TH1F * hNorm = (TH1F *) hHist->Clone(Form("%s_%s",hHist->GetName(),sNormName.Data()));
	double scale = hNorm->Integral("width") * rescale;
	if (scale != 0) hNorm->Scale(1.0/scale);
	return hNorm;
}

// First index is the cent bin
// within each subarray, [0] is the main, [1] is the error
vector<vector<double>> SumUpProfile(TProfile2D * fProf2D) {
	vector<vector<double>> fVector = {};
	int nBinsX = fProf2D->GetXaxis()->GetNbins();
	int nBinsY = fProf2D->GetYaxis()->GetNbins();

	TProfile * fProf = fProf2D->ProfileY(Form("%s_pfy",fProf2D->GetName()));
	//TProfile * fProf = fProf2D->ProfileY(Form("%s_pfy_Cent%d",fProf2D->GetName(),CentBin));
	for (int i = 0; i < nBinsY; i++) {
		double localVal = fProf->GetBinContent(i+1);
		double localErr = fProf->GetBinError(i+1);
		vector<double> fLocalVector = {localVal,localErr};
		fVector.push_back(fLocalVector);
	}
	return fVector;
}



// sName is title to add to plots
//void  DrawPi0Info(TString sFilename = "AnalysisResults.root", TString sName = "") {
void  DrawPi0Info(TString sName = "",TString sFilename = "AnalysisResults.root") {



	TFile * file = TFile::Open(sFilename.Data(),"READ");

	// Loop over tasks with Pi0H in it
	// make a directory for each? Or just prefix each output filename

	// Cand also look at PWGJE_QnVectorTender

	if (!file) {
		fprintf(stderr,"File %s not found.\n",sFilename.Data());
		return;
	}


	vector <TObject * > fDirectories = {};
	// Pi0H task vector
	//vector <TObject * > fPi0HTasks = {};
	vector <TList * > fPi0HTasks = {};
	// QnVector task


	TList * keys = file->GetListOfKeys();

	TObject * obj,* obj2;
	TIter next(keys);
	while ((obj = next())) {
		if (!obj->InheritsFrom("TNamed")) continue;
    if (!obj->InheritsFrom("TKey")) continue;



    TKey *kobj = dynamic_cast<TKey*>(obj);
    obj2 = kobj->ReadObj();
		TString sObjName = obj2->GetName();
		if (obj2->InheritsFrom("AliEmcalList") && sObjName.Contains("Pi0H")) {
			fPi0HTasks.push_back((TList*)obj2);

		}
	}


	bool hasName = sName != "";
	//TLegend * legName = new TLegend();
	TLegend * legName = new TLegend(0.7,0.4,1.0,0.5);
	//TLegend * legName = new TLegend(0.4,0.4,0.6,0.6);
	legName->SetHeader(sName.Data(),"c");

	//TPaveText * ptName = new TPaveText(0.45,0.45,0.975,0.6);
	//TPaveText * ptName = new TPaveText(0.45,0.45,0.975,0.6,"NDC");
	TPaveText * ptName = new TPaveText(0.85,0.45,0.975,0.5,"NDC");
	ptName->AddText(sName.Data());

  //TPaveText *pt = new TPaveText(.05,.1,.95,.8);
	//    pt->AddText("A TPaveText can contain severals line of text.");

	for (auto fTask : fPi0HTasks) {
		printf("Found Pi0H Directory %s\n",fTask->GetName());

		TString sTaskName = fTask->GetName();
		system(Form("mkdir -p %s",sTaskName.Data()));

		int iCentBin = 2; // fixed for now

		vector<TH1F *> hPtEPAngleTrack_Proj = {};
		vector<TH1F *> hPtEP3AngleTrack_Proj = {};

		vector<TH1F *> hPtEPAngleTrigger_Proj = {};
		vector<TH1F *> hPtEP3AngleTrigger_Proj = {};


		TString sOutputFilePath = Form("%s/QAOutput.root",sTaskName.Data());
		TFile * fQAOutput = new TFile(sOutputFilePath,"RECREATE");
		
		TString k2DStyle = "COLZ";

		TH1D * fHistBinCheckPt_0 = (TH1D *) fTask->FindObject("fHistBinCheckPt_0");
		fHistBinCheckPt_0->Draw();
		ptName->Draw("SAME");
		gPad->Print(Form("%s/HistBinCheckPt.png",sTaskName.Data()));	

		double nEvents=0;

		vector<double> fNEventsCent = {0.,0.,0.,0.};
		vector<double> fNEventScaleCent = {0.,0.,0.,0.};

		TH1F * fHistCentrality = (TH1F *) fTask->FindObject("fHistCentrality");
		TH1F * fHistCentralityNorm = 0;
		if (fHistCentrality) {
			fHistCentralityNorm = NormalizeTH1(fHistCentrality);
			fHistCentrality->Draw();
			if (hasName)	{
				cout<<"This has a name"<<endl;
			//	legName->AddEntry(fHistCentrality,"Centrality","lp");
		//		legName->Draw("SAME");
				ptName->Draw("SAME");
			}
			gPad->Print(Form("%s/Centrality.pdf",sTaskName.Data()));	
			gPad->Print(Form("%s/Centrality.png",sTaskName.Data()));	




			// Count number of events
			
			//double nEvents = fHistCentrality->GetEntries();
			for (int i = 0; i < kNCentBins; i++) {

				double nCentEntries = fHistCentrality->GetEntries();
				//int fHistCent_min = fHistCentrality->GetXaxis()->FindBin(fCentBins[iCentBin]);
				//int fHistCent_max = fHistCentrality->GetXaxis()->FindBin(fCentBins[iCentBin+1]-0.0001); // binning
				int fHistCent_min = fHistCentrality->GetXaxis()->FindBin(fCentBins[i]);
				int fHistCent_max = fHistCentrality->GetXaxis()->FindBin(fCentBins[i+1]-0.0001); // binning
				double fHistCent_MinEdge = fHistCentrality->GetXaxis()->GetBinLowEdge(fHistCent_min);
				double fHistCent_MaxEdge = fHistCentrality->GetXaxis()->GetBinUpEdge(fHistCent_max);
				printf("== CentBin = %d ==\n",i);
				printf("Getting event count by integrating centrality histogram from %d to %d\n",fHistCent_min,fHistCent_max);
				printf("     This corresponds to %f - %f \n",fHistCent_MinEdge,fHistCent_MaxEdge);
				fNEventsCent[i] = fHistCentrality->Integral(fHistCent_min,fHistCent_max);


			}


			nEvents = fNEventsCent[2]; // temporary

		} else {
			printf("fHistCentrality not found.\n");
		}



		fQAOutput->Add(fHistCentrality);
		fQAOutput->Add(fHistCentralityNorm);

		double fNEventScale = 0;
		if (nEvents != 0) fNEventScale = 1./nEvents;

		TH2F * fEtaPhiPionAcc = (TH2F *) fTask->FindObject("fEtaPhiPionAcc");

		fEtaPhiPionAcc->Draw(k2DStyle);
		ptName->Draw("SAME");
		gPad->Print(Form("%s/EtaPhiPionAccepted.pdf",sTaskName.Data()));	
		gPad->Print(Form("%s/EtaPhiPionAccepted.png",sTaskName.Data()));	

		TH1F * fPionAccEta = (TH1F *) fEtaPhiPionAcc->ProjectionX("PionAccEta");
		fPionAccEta->SetTitle("Accepted #pi^{0} #eta");
		fPionAccEta->Scale(fNEventScale,"width");

		TH1F * fPionAccPhi = (TH1F *) fEtaPhiPionAcc->ProjectionY("PionAccPhi");
		fPionAccPhi->SetTitle("Accepted #pi^{0} #phi");
		fPionAccPhi->Scale(fNEventScale,"width");




		TList * fEventCutOutput = (TList *) fTask->FindObject("EventCutOutput");

		TH1D * fCutStats = (TH1D *) fEventCutOutput->FindObject("fCutStats");
		//TH1D * fCutStats = (TH1D *) fEventCutOutput->FindObject("EventCutOutput/fCutStats");
//		TH1D * fCutStats = (TH1D *) fTask->Get("EventCutOutput/fCutStats");
		if (fCutStats) {
			fCutStats->Draw();
			ptName->Draw("SAME");
			gPad->Print(Form("%s/CutStats.pdf",sTaskName.Data()));	
			gPad->Print(Form("%s/CutStats.png",sTaskName.Data()));	
		} else {
			cout<<"fCutStats not found"<<endl;
		}

		TH1D * fRawCentrality = (TH1D *) fEventCutOutput->FindObject("Centrality_raw");

		if (fRawCentrality) {
			fRawCentrality->Draw();
			ptName->Draw("SAME");
			gPad->Print(Form("%s/RawCentrality.pdf",sTaskName.Data()));	
			gPad->Print(Form("%s/RawCentrality.png",sTaskName.Data()));	
		}

		// Clusters
		TH1F * hClusEnergy = (TH1F *) fTask->FindObject("ClusEnergy");
		double fMinClusEnergy = hClusEnergy->GetXaxis()->GetXmin();
		double fMaxClusEnergy = hClusEnergy->GetXaxis()->GetXmax();
		printf("Clus Energy has range %f - %f with %d bins to begin.\n",fMinClusEnergy,fMaxClusEnergy,hClusEnergy->GetNbinsX());

		//hClusEnergy->Rebin(20);

		// rebinning
		int nNewBins = 23;
		vector<double> fNewClusEnergyBins = {0.0};
		double fStepSize = 0.5;
		for (int i = 1; i <= nNewBins; i++) {
			fNewClusEnergyBins.push_back(fNewClusEnergyBins[i-1] + fStepSize);
			if (i == 10) fStepSize = 1.0;
			if (i == 15) fStepSize = 2.0;
			if (i == 20) fStepSize = 10.;
		}
		printf("Bin edges: ");
		for (double fEdge : fNewClusEnergyBins) printf("%.2f ",fEdge);
		printf("\n");

		TH1F * hClusEnergyRebin = (TH1F *) hClusEnergy->Rebin(nNewBins,"ClusEnergyRebin",&fNewClusEnergyBins[0]);
		
		hClusEnergyRebin->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dE}");


		//hClusEnergy->Draw();
		hClusEnergyRebin->Scale(fNEventScale,"width");
		hClusEnergyRebin->Draw();
		gPad->SetLogy(1);
		gPad->SetGridx(1);
		ptName->Draw("SAME");
		gPad->Print(Form("%s/ClusEnergy.pdf",sTaskName.Data()));	
		gPad->Print(Form("%s/ClusEnergy.png",sTaskName.Data()));	
		gPad->SetLogy(0);
		fQAOutput->Add(hClusEnergy);

		TH1F * hClusEnergyRebinNorm = NormalizeTH1(hClusEnergyRebin);
		TH1F * hClusEnergyNorm = NormalizeTH1(hClusEnergy);




		// Accepted Pi0s

		TH2F * fMassPtPionAcc = (TH2F *) fTask->FindObject("fMassPtPionAcc");
		TH3F * fMassPtCentPionAcc = (TH3F *) fTask->FindObject("fMassPtCentPionAcc");

		vector<TH1F *> fMassPtBins = {};
		if (fMassPtPionAcc) {
			for (int i = 0; i < kNTrackPtBins; i++) {
				double fTriggerMinPt = fTriggerPtBins[i];
				double fTriggerMaxPt = fTriggerPtBins[i+1];
				TString sTriggerName = Form("%.0f_%.0f",fTriggerMinPt,fTriggerMaxPt);
				TString sTriggerTitle = Form("%.0f #leq p_{T} < %.0f",fTriggerMinPt,fTriggerMaxPt);
				int iMinBinY = fMassPtPionAcc->GetYaxis()->FindBin(fTriggerMinPt);
				int iMaxBinY = fMassPtPionAcc->GetYaxis()->FindBin(fTriggerMaxPt  - 0.0001);
				TH1F * fProjX = (TH1F *) fMassPtPionAcc->ProjectionX(Form("Mass_Pt_%s",sTriggerName.Data()),iMinBinY,iMaxBinY);
				fProjX->SetTitle(sTriggerTitle.Data());
				fProjX->GetXaxis()->SetRangeUser(0.1,0.5);
				fProjX->Draw();
				gPad->Print(Form("%s/MassPt_%s.png",sTaskName.Data(),sTriggerName.Data()));
				fMassPtBins.push_back(fProjX);
			}
		}






		// Tracks
		TH2F * hHistNChargedCent = (TH2F *) fTask->FindObject("fHistNChargedCent");
		TH3F * hHistTrackPsiEPPtCent = (TH3F *) fTask->FindObject("fHistTrackPsiEPPtCent");
		TH3F * hHistTrackPsiEP3PtCent = (TH3F *) fTask->FindObject("fHistTrackPsiEP3PtCent");

		TH1F * hNCharged = 0;
		if (hHistNChargedCent) {
			hNCharged =	(TH1F *) hHistNChargedCent->ProjectionX("NCharged",iCentBin+1,iCentBin+1,"e");
			hNCharged->Draw();
			ptName->Draw("SAME");
			gPad->Print(Form("%s/NCharged.png",sTaskName.Data()));
			fQAOutput->Add(hNCharged);
		}


		//HistTrackPsiEPPtCent->GetZaxis()->SetRange(iCentBin+1,iCentBin+1);
		TH1F * hTrackPt = (TH1F *) hHistTrackPsiEPPtCent->ProjectionY("TrackPt",0,-1,iCentBin+1,iCentBin+1,"e");

		TH1F * hTrackPtRebin = (TH1F *) hTrackPt->Rebin(nNewBins,"TrackPtRebin",&fNewClusEnergyBins[0]);
		hTrackPtRebin->SetTitle("Track p_{T}");
		hTrackPtRebin->GetYaxis()->SetTitle(Form("#frac{1}{N_{ev}} #frac{dN}{d%s}",hTrackPt->GetXaxis()->GetTitle()));


		//hTrackPt->Rebin(4);
		hTrackPt->Scale(fNEventScale,"width");
		hTrackPtRebin->Scale(fNEventScale,"width");

		hTrackPt->SetTitle("Track p_{T}");
		TH1F * hTrackEP = (TH1F *) hHistTrackPsiEPPtCent->ProjectionX("TrackEP",0,-1,iCentBin+1,iCentBin+1,"e");
		hTrackEP->Rebin(4);

	//	hTrackPt->Scale(fNEventScale,"width");

		hTrackPtRebin->Draw();
		gPad->SetLogy(1);
		ptName->Draw("SAME");
		gPad->Print(Form("%s/TrackPt.png",sTaskName.Data()));	

		gPad->SetLogy(0);
		hTrackEP->Draw();
		ptName->Draw("SAME");
		gPad->Print(Form("%s/TrackEP.png",sTaskName.Data()));	

		TH1F * hTrackPtNorm = NormalizeTH1(hTrackPt);
		TH1F * hTrackEPNorm = NormalizeTH1(hTrackEP);

		hHistTrackPsiEPPtCent->GetZaxis()->SetRange(iCentBin+1,iCentBin+1);
		TH2F * hHistTrackEPPt = (TH2F *) hHistTrackPsiEPPtCent->Project3D("yx");


    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

 //     printf("Projecting Tracks Event Plane Histograms in bin from %.1f to %.1f\n",hPtEPAnglePionAcc->GetYaxis()->GetBinLowEdge(iMinBin),hPtEPAnglePionAcc->GetYaxis()->GetBinUpEdge(iMaxBin));

      TString sFormat = "%s_Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

//(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
      hHistTrackPsiEPPtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEPAngleTrack_Proj = (TH1F *) hHistTrackPsiEPPtCent->Project3D("xe");
      hLocalPtEPAngleTrack_Proj->SetName(Form(sFormat.Data(),hHistTrackPsiEPPtCent->GetName(),i));
			// FIXME normalizing by nEvents, which might include all centrality bins currently
			hLocalPtEPAngleTrack_Proj->Scale(fNEventScale);

      hLocalPtEPAngleTrack_Proj->Rebin(4);
      hLocalPtEPAngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP} (%s)",sPtRange.Data()));
      hLocalPtEPAngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEPAngleTrack_Proj.push_back(hLocalPtEPAngleTrack_Proj);

			fQAOutput->Add(hLocalPtEPAngleTrack_Proj);
    }

		// Repeat for EP3
    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      TString sFormat = "%s_Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

      hHistTrackPsiEP3PtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEP3AngleTrack_Proj = (TH1F *) hHistTrackPsiEP3PtCent->Project3D("xe");
      hLocalPtEP3AngleTrack_Proj->SetName(Form(sFormat.Data(),hHistTrackPsiEP3PtCent->GetName(),i));
			// FIXME normalizing by nEvents, which might include all centrality bins currently
			hLocalPtEP3AngleTrack_Proj->Scale(fNEventScale);

      hLocalPtEP3AngleTrack_Proj->Rebin(4);
      hLocalPtEP3AngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP3} (%s)",sPtRange.Data()));
      hLocalPtEP3AngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEP3AngleTrack_Proj.push_back(hLocalPtEP3AngleTrack_Proj);

			fQAOutput->Add(hLocalPtEP3AngleTrack_Proj);
    }





		// Triggers (Pi0s) vs Event Plane
		TH3F * hPtEPAnglePionAccCent = (TH3F *) fTask->FindObject("PtEPAnglePionAccCent");



    for (int i = 0; i < kNPi0PtBins; i++) {
      double fMinPt = fTriggerPtBins[i];
      double fMaxPt = fTriggerPtBins[i+1];

      int iMinBin = hPtEPAnglePionAccCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hPtEPAnglePionAccCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      TString sFormat = "%s_Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

			hPtEPAnglePionAccCent->GetYaxis()->SetRange(iMinBin,iMaxBin);

      TH1F * hLocalPtEPAngleTrigger_Proj = (TH1F *) hPtEPAnglePionAccCent->Project3D("xe");
      hLocalPtEPAngleTrigger_Proj->SetName(Form(sFormat.Data(),hPtEPAnglePionAccCent->GetName(),i));
			// FIXME normalizing by nEvents, which might include all centrality bins currently
			hLocalPtEPAngleTrigger_Proj->Scale(fNEventScale);
      hLocalPtEPAngleTrigger_Proj->Rebin(4);
      hLocalPtEPAngleTrigger_Proj->SetTitle(Form("Trigger #Delta#Psi_{EP} (%s)",sPtRange.Data()));
      hLocalPtEPAngleTrigger_Proj->GetYaxis()->SetTitle("N_{Triggers}");
      hPtEPAngleTrigger_Proj.push_back(hLocalPtEPAngleTrigger_Proj);

			fQAOutput->Add(hLocalPtEPAngleTrigger_Proj);

		}


		//TH2F * fMassPtPionAcc = (TH2F *) fTask->FindObject("fMassPtPionAcc");
		fMassPtPionAcc->Draw(k2DStyle);
		ptName->Draw("SAME");
		gPad->Print(Form("%s/PionAcceptedMassPt.pdf",sTaskName.Data()));	
		TH2F * fMassPtPionRej = (TH2F *) fTask->FindObject("fMassPtPionRej");
		fMassPtPionRej->Draw(k2DStyle);
		ptName->Draw("SAME");
		gPad->Print(Form("%s/PionRejectedMassPt.pdf",sTaskName.Data()));	

		TH2F * fHistEventHashVsMixingAngle = (TH2F *) fTask->FindObject("HistEventHashVsMixingAngle");


		TH1F * fHistEventPlane = (TH1F *) fTask->FindObject("fHistEventPlane");
		TH1F * fEPAngleV0M = (TH1F *) fTask->FindObject("EPAngleV0M");
		TH1F * fEPAngleTPCA = (TH1F *) fTask->FindObject("EPAngleTPCA");
		TH1F * fEPAngleTPCC = (TH1F *) fTask->FindObject("EPAngleTPCC");
		TH1F * fEP3AngleV0M = (TH1F *) fTask->FindObject("EP3AngleV0M");
		TH1F * fEP3AngleTPCA = (TH1F *) fTask->FindObject("EP3AngleTPCA");
		TH1F * fEP3AngleTPCC = (TH1F *) fTask->FindObject("EP3AngleTPCC");

		SetHistStyle(fHistEventPlane,0);
		SetHistStyle(fEPAngleV0M,1);
		SetHistStyle(fEPAngleTPCA,2);
		SetHistStyle(fEPAngleTPCC,3);
		SetHistStyle(fEP3AngleV0M,4);
		SetHistStyle(fEP3AngleTPCA,5);
		SetHistStyle(fEP3AngleTPCC,6);

		if (fHistEventHashVsMixingAngle) {
			fHistEventHashVsMixingAngle->SetTitle(fTask->GetName());
			fHistEventHashVsMixingAngle->Draw("COLZ");
			gPad->Print(Form("%s/EventHashVsMixingAngle.pdf",sTaskName.Data()));	
			gPad->Print(Form("%s/EventHashVsMixingAngle.png",sTaskName.Data()));
		}


		fHistEventPlane->Draw();
		fEPAngleV0M->Draw("SAME");
		ptName->Draw("SAME");
	gPad->BuildLegend();
	gPad->Print(Form("%s/EventPlane.png",sTaskName.Data()));

	double fYMax = fHistEventPlane->GetBinContent(fHistEventPlane->GetMaximumBin());
	double fYMin = 0;
	fHistEventPlane->GetYaxis()->SetRangeUser(fYMin,fYMax+0.3*(fYMax-fYMin));

	fHistEventPlane->Draw();
	fEPAngleV0M->Draw("SAME");
	fEPAngleTPCA->Draw("SAME");
	fEPAngleTPCC->Draw("SAME");
	fEP3AngleV0M->Draw("SAME");
	fEP3AngleTPCA->Draw("SAME");
	fEP3AngleTPCC->Draw("SAME");
	ptName->Draw("SAME");
		gPad->BuildLegend();
		gPad->Print(Form("%s/EventPlaneCmp.png",sTaskName.Data()));
		gPad->Print(Form("%s/EventPlaneCmp.pdf",sTaskName.Data()));
		fQAOutput->Add(fHistEventPlane);
		fQAOutput->Add(fEPAngleV0M);
		fQAOutput->Add(fEP3AngleV0M);

		TH1F * fHistEventPlaneNorm = NormalizeTH1(fHistEventPlane);

		TH1F * fEPAngleV0M_Norm = NormalizeTH1(fEPAngleV0M);
		TH1F * fEPAngleTPCA_Norm = NormalizeTH1(fEPAngleTPCA);
		TH1F * fEPAngleTPCC_Norm = NormalizeTH1(fEPAngleTPCC);
		TH1F * fEP3AngleV0M_Norm = NormalizeTH1(fEP3AngleV0M);
		TH1F * fEP3AngleTPCA_Norm = NormalizeTH1(fEP3AngleTPCA);
		TH1F * fEP3AngleTPCC_Norm = NormalizeTH1(fEP3AngleTPCC);

		fQAOutput->Add(fHistEventPlaneNorm);
		fQAOutput->Add(fEPAngleV0M_Norm);
		fQAOutput->Add(fEP3AngleV0M_Norm);
		fQAOutput->Add(fEPAngleTPCA_Norm);
		fQAOutput->Add(fEP3AngleTPCA_Norm);
		fQAOutput->Add(fEPAngleTPCC_Norm);
		fQAOutput->Add(fEP3AngleTPCC_Norm);

		int nRebinScaleVsAngleX = 4;
		int nRebinScaleVsAngleY = 4;

		TH2F * hQ2V0MScaleVsAngle = (TH2F *) fTask->FindObject("Q2V0MScaleVsAngle");
		TProfile * hQ2V0MScale = 0;

		if (hQ2V0MScaleVsAngle) {
			hQ2V0MScaleVsAngle->Rebin2D(nRebinScaleVsAngleX,nRebinScaleVsAngleY);
			hQ2V0MScaleVsAngle->Draw("COLZ");
			hQ2V0MScale = hQ2V0MScaleVsAngle->ProfileX("Q2V0Scale");
			hQ2V0MScale->Draw("SAME");
				ptName->Draw("SAME");
			gPad->SetLogz(1);
			gPad->Print(Form("%s/Q2ScaleVsAngle_V0M.png",sTaskName.Data()));
			//fQAOutput->Add(hQ2V0MScale);
		}

		TH1F * hQ2Scale = 0;
		if (hQ2V0MScaleVsAngle) {
			hQ2Scale = (TH1F *) hQ2V0MScaleVsAngle->ProjectionY("Q2Scale");
			hQ2Scale->SetTitle("|q_{2}|");
			hQ2Scale->Rebin(5);
			fQAOutput->Add(hQ2Scale);
			hQ2Scale->Draw();
			ptName->Draw("SAME");
			gPad->Print(Form("%s/Q2Scale.png",sTaskName.Data()));
			TH1F * hQ2ScaleNorm = NormalizeTH1(hQ2Scale);
			fQAOutput->Add(hQ2ScaleNorm);
		}




		// Event Plane Resolution
    //# Event Plane 2
    //# EP3,4 stored as EP3R_CosD%d_N%d, EP4R_...

		vector<TGraphErrors *> fEP2ResGraphs = {};
		vector<TGraphErrors *> fEP3ResGraphs = {};

    printf("Finding EPRs for Event Plane 2 and 3\n");
    //for i in range(6): # individual N values
      //LocalVals = []x
      //for j in range(3): # individual Dn values
        //Pf2Name="EPR_CosD%d_N%d" % (j+1,i+1)
		for (int i = 0; i < 6; i++) {

			TGraphErrors * fEP2ResGraph = new TGraphErrors(4);
			fEP2ResGraph->SetName(Form("EP2_%d_ResGraph",i+1));
			fEP2ResGraph->SetTitle(Form("EPR_{2,%d};Centrality;EPR_{2,%d}",i+1,i+1));


			TGraphErrors * fEP3ResGraph = new TGraphErrors(4);
			fEP3ResGraph->SetName(Form("EP3_%d_ResGraph",i+1));
			fEP3ResGraph->SetTitle(Form("EPR_{3,%d};Centrality;EPR_{3,%d}",i+1,i+1));

			int nFoundPoints = 0;


			vector<vector<vector<double>>> fLocalVals = {};
			for (int j = 0; j < 3; j++) {
				TString sEPProfileName = Form("EPR_CosD%d_N%d",j+1,i+1);

				TProfile2D * fEPProfile = (TProfile2D *) fTask->FindObject(sEPProfileName.Data());

				if (fEPProfile) {
				} else {
					fprintf(stderr,"Missing fEPProfile %s\n",sEPProfileName.Data());
					continue;
				}

				vector<vector<double>> fLocalArray = SumUpProfile(fEPProfile);
				fLocalVals.push_back(fLocalArray);
			}

			vector<vector<vector<double>>> fEP3LocalVals = {};
			for (int j = 0; j < 3; j++) {
				TString sEPProfileName = Form("EP3R_CosD%d_N%d",j+1,i+1);

				TProfile2D * fEPProfile = (TProfile2D *) fTask->FindObject(sEPProfileName.Data());

				if (fEPProfile) {
				} else {
					fprintf(stderr,"Missing fEPProfile %s\n",sEPProfileName.Data());
					continue;
				}

				vector<vector<double>> fLocalArray = SumUpProfile(fEPProfile);
				fEP3LocalVals.push_back(fLocalArray);
			}


			for (int jCent = 0; jCent < 4; jCent++) {
				printf("    Calculating EPR for CentBin %d\n",jCent);
				double fjCentBinCenter = 0.5 * (fCentBins[jCent+1] + fCentBins[jCent]);
				double fjCentBinWidth  = 0.5 * (fCentBins[jCent+1] - fCentBins[jCent]);

				double fLocalRn = 0;
				double fLocalRn_Un = 0;
				double MeanCosD1 = fLocalVals[0][jCent][0];
				double MeanCosD2 = fLocalVals[1][jCent][0];
				double MeanCosD3 = fLocalVals[2][jCent][0];
				double MeanCosD1_Un = fLocalVals[0][jCent][1];
				double MeanCosD2_Un = fLocalVals[1][jCent][1];
				double MeanCosD3_Un = fLocalVals[2][jCent][1];
				if (MeanCosD3 < 0 || MeanCosD2 < 0 || MeanCosD3 < 0) {
					printf("Unusually circumstance: some <cos> value is negative\n");
				}
				if (MeanCosD3 > 0 && MeanCosD2 > 0 && MeanCosD1 > 0) {
					fLocalRn = TMath::Sqrt((MeanCosD1 * MeanCosD2) / MeanCosD3);
					fLocalRn_Un = fLocalRn * TMath::Sqrt((0.5)*TMath::Power(MeanCosD1_Un/MeanCosD1,2) + (0.5)*TMath::Power(MeanCosD2_Un/MeanCosD2,2) + (0.5)*TMath::Power(MeanCosD3_Un/MeanCosD3,2));
					printf("Found R_{%d,2} = %f   \\pm %f\n",i+1,fLocalRn,fLocalRn_Un);

				}
				fEP2ResGraph->SetPoint(jCent,fjCentBinCenter,fLocalRn);
				fEP2ResGraph->SetPointError(jCent,fjCentBinWidth,fLocalRn_Un);
			}

			for (int jCent = 0; jCent < 4; jCent++) {
				printf("    Calculating EP3R for CentBin %d\n",jCent);
				double fjCentBinCenter = 0.5 * (fCentBins[jCent+1] + fCentBins[jCent]);
				double fjCentBinWidth  = 0.5 * (fCentBins[jCent+1] - fCentBins[jCent]);

				double fLocalRn = 0;
				double fLocalRn_Un = 0;
				double MeanCosD1 = fEP3LocalVals[0][jCent][0];
				double MeanCosD2 = fEP3LocalVals[1][jCent][0];
				double MeanCosD3 = fEP3LocalVals[2][jCent][0];
				double MeanCosD1_Un = fEP3LocalVals[0][jCent][1];
				double MeanCosD2_Un = fEP3LocalVals[1][jCent][1];
				double MeanCosD3_Un = fEP3LocalVals[2][jCent][1];
				if (MeanCosD3 < 0 || MeanCosD2 < 0 || MeanCosD3 < 0) {
					printf("Unusually circumstance: some <cos> value is negative\n");
				}
				if (MeanCosD3 > 0 && MeanCosD2 > 0 && MeanCosD1 > 0) {
					fLocalRn = TMath::Sqrt((MeanCosD1 * MeanCosD2) / MeanCosD3);
					fLocalRn_Un = fLocalRn * TMath::Sqrt((0.5)*TMath::Power(MeanCosD1_Un/MeanCosD1,2) + (0.5)*TMath::Power(MeanCosD2_Un/MeanCosD2,2) + (0.5)*TMath::Power(MeanCosD3_Un/MeanCosD3,2));
					printf("Found R_{%d,3} = %f   \\pm %f\n",i+1,fLocalRn,fLocalRn_Un);

				}
				fEP3ResGraph->SetPoint(jCent,fjCentBinCenter,fLocalRn);
				fEP3ResGraph->SetPointError(jCent,fjCentBinWidth,fLocalRn_Un);
			}



			fEP2ResGraphs.push_back(fEP2ResGraph);
			fEP3ResGraphs.push_back(fEP3ResGraph);
		}
		printf("Done with first step of  analyzing EPR stuff\n");

		TMultiGraph * mgEPR = new TMultiGraph();
		for (int i = 0; i < 6; i++) {
			fEP2ResGraphs[i]->SetLineColor(colorList[i]);
			fEP2ResGraphs[i]->SetMarkerColor(colorList[i]);
			fEP2ResGraphs[i]->SetMarkerStyle(markerList[i]);
			mgEPR->Add(fEP2ResGraphs[i]);
			//if (i == 0) fEP2ResGraphs[i]->Draw("ALP");
			//else fEP2ResGraphs[i]->Draw("LP SAME");
		}
		mgEPR->GetXaxis()->SetTitle(fEP2ResGraphs[0]->GetXaxis()->GetTitle());
		mgEPR->GetYaxis()->SetTitle(fEP2ResGraphs[0]->GetYaxis()->GetTitle());
		mgEPR->Draw("ALP");
		gPad->BuildLegend();
		gPad->Print(Form("%s/EPR.png",sTaskName.Data()));	


		for (auto graph : fEP2ResGraphs) fQAOutput->Add(graph);
		for (auto graph : fEP3ResGraphs) fQAOutput->Add(graph);

		double fEPRes_R2 = 1;
		double fEP3Res_R3 = 1;
		double fEPRes_R4 = 1;
		double fEPRes_R6 = 1;
		double fEPRes_R2_Un = 1;
		double fEP3Res_R3_Un = 1;
		double fEPRes_R4_Un = 1;
		double fEPRes_R6_Un = 1;

		if (fEP2ResGraphs[2-1]->GetY()[iCentBin] > 0) {
			fEPRes_R2 = fEP2ResGraphs[2-1]->GetY()[iCentBin];
			fEPRes_R2_Un = fEP2ResGraphs[2-1]->GetEY()[iCentBin];
		}

		if (fEP3ResGraphs[3-1]->GetY()[iCentBin] > 0) {
			fEP3Res_R3 = fEP3ResGraphs[3-1]->GetY()[iCentBin];
			fEP3Res_R3_Un = fEP3ResGraphs[3-1]->GetEY()[iCentBin];
		}

		if (fEP2ResGraphs[4-1]->GetY()[iCentBin] > 0) {
			fEPRes_R4 = fEP2ResGraphs[4-1]->GetY()[iCentBin];
			fEPRes_R4_Un = fEP2ResGraphs[4-1]->GetEY()[iCentBin];
		}
		if (fEP2ResGraphs[6-1]->GetY()[iCentBin] > 0) {
			fEPRes_R6 = fEP2ResGraphs[6-1]->GetY()[iCentBin];
			fEPRes_R6_Un = fEP2ResGraphs[6-1]->GetEY()[iCentBin];
		}



		TCanvas * cVn = new TCanvas("cVn","cVn");
		cVn->SetTopMargin(0.1);
		cVn->SetBottomMargin(0.1);
		cVn->SetLeftMargin(0.1);
		cVn->SetRightMargin(0.1);


		TGraphErrors * gTrack_ChiSq = new TGraphErrors(kNTrackPtBins);
		gTrack_ChiSq->SetName("Track_ChiSq");
		gTrack_ChiSq->SetTitle("ChiSquare/NDF from V_{n} fit");
		gTrack_ChiSq->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gTrack_ChiSq->GetYaxis()->SetTitle("#chi^{2}/NDF");

		TGraphErrors * gTrack_Bv = new TGraphErrors(kNTrackPtBins);
		gTrack_Bv->SetName("Track_Bv");
		gTrack_Bv->SetTitle("B Value from V_{n} fit");
		gTrack_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gTrack_Bv->GetYaxis()->SetTitle("B");

		TGraphErrors * gTrack_V2 = new TGraphErrors(kNTrackPtBins);
		gTrack_V2->SetName("Track_V2");
		gTrack_V2->SetTitle("Calculated v_{2}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
		//gTrack_V2->SetTitle("Calculated #tilde{v}_{2}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
		gTrack_V2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gTrack_V2->GetYaxis()->SetTitle("v_{2}");
		//gTrack_V2->GetYaxis()->SetTitle("#tilde{v}_{2}");
		gTrack_V2->SetMarkerStyle(kOpenSquare);
		TGraphErrors * gTrack_V2_Raw = (TGraphErrors *) gTrack_V2->Clone("Track_V2_Raw");
		gTrack_V2_Raw->SetTitle("Calculated v_{2}^{Track} (Event Plane method, No EPR Correction))");
		gTrack_V2_Raw->GetYaxis()->SetTitle("v_{2}(raw)");

		TGraphErrors * gTrack_V4 = new TGraphErrors(kNTrackPtBins);
		gTrack_V4->SetName("Track_V4");
		gTrack_V4->SetTitle("Calculated v_{4}^{Track} (Event Plane method)");
		//gTrack_V4->SetTitle("Calculated #tilde{v}_{4}^{Track} (Event Plane method)");
		gTrack_V4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gTrack_V4->GetYaxis()->SetTitle("v_{4}");
		gTrack_V4->SetMarkerStyle(kFullSquare);
		//gTrack_V4->GetYaxis()->SetTitle("#tilde{v}_{4}");
		TGraphErrors * gTrack_V4_Raw = (TGraphErrors *) gTrack_V4->Clone("Track_V4_Raw");
		gTrack_V4_Raw->SetTitle("Calculated v_{4}^{Track} (Event Plane method, No EPR Correction))");
		gTrack_V4_Raw->GetYaxis()->SetTitle("v_{4}(raw)");

		TGraphErrors * gTrack_V6 = new TGraphErrors(kNTrackPtBins);
		gTrack_V6->SetName("Track_V6");
		gTrack_V6->SetTitle("Calculated v_{6}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
		//gTrack_V6->SetTitle("Calculated #tilde{v}_{6}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
		gTrack_V6->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gTrack_V6->GetYaxis()->SetTitle("v_{6}");
		//gTrack_V6->GetYaxis()->SetTitle("#tilde{v}_{6}");
		TGraphErrors * gTrack_V6_Raw = (TGraphErrors *) gTrack_V6->Clone("Track_V6_Raw");
		gTrack_V6_Raw->SetTitle("Calculated v_{6}^{Track} (Event Plane method, No EPR Correction))");
		gTrack_V6_Raw->GetYaxis()->SetTitle("v_{6}(raw)");

		printf("Starting the fitting of tracks vs the event plane\n");
		for (int i = 0; i < kNTrackPtBins; i++) {
			printf("Fitting track pt bin %d\n",i);
			double fMinPt = fTrackPtBins[i];
			double fMaxPt = fTrackPtBins[i+1];

			// FIXME
			TH1F * hTrackEPLocal = 0;
			hTrackEPLocal = hPtEPAngleTrack_Proj[i];
			if (!hTrackEPLocal) {
				printf("Missing the track vs EP histogram\n");
				return;
			}

			TF1 * fitTrackEP = new TF1(Form("Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",hTrackEPLocal->GetXaxis()->GetXmin(),hTrackEPLocal->GetXaxis()->GetXmax());
			
			printf("Made a fit  function. Setting initial parameters ...\n");
			cout << "\tFitting using histogram " << hTrackEPLocal->GetName() <<endl;
			cout << "\tand fit function " << fitTrackEP->GetName() <<endl;
			fitTrackEP->SetParameter(0,hTrackEPLocal->Integral("width") / (TMath::Pi() / 2));
			fitTrackEP->SetParameter(1,0.01);
			fitTrackEP->SetParameter(2,0.001);
			fitTrackEP->SetParameter(3,0.);

			cout << "Setting par limits" <<endl;
			fitTrackEP->SetParLimits(1,0.0,0.5);
			fitTrackEP->SetParLimits(2,0.0,0.25);
			fitTrackEP->SetParLimits(3,0.0,0.25);

			cout << "Setting par names" <<endl;
			fitTrackEP->SetParName(0,"B");
			fitTrackEP->SetParName(1,"v_2");
			fitTrackEP->SetParName(2,"v_4");
			fitTrackEP->SetParName(3,"v_6");
			cout << "About to fit" <<endl;

			if (hTrackEPLocal->GetEntries() == 0) {
				printf("Skipping pt bin because of empty histogram\n");
				continue;
			}

			hTrackEPLocal->Fit(fitTrackEP);
			fitTrackEP->SetLineColor(kCyan);

			hTrackEPLocal->SetMarkerStyle(kFullSquare);
			hTrackEPLocal->Draw();
			fitTrackEP->Draw("SAME");

			double fEPFitChiSquare = fitTrackEP->GetChisquare();
			double fEPFitNDF = fitTrackEP->GetNDF();
			double fEPFitChi2OverNDF = 0;
			if (fEPFitNDF > 0) {
				fEPFitChi2OverNDF = fEPFitChiSquare / fEPFitNDF;
			}

			printf("Found raw v2(track) = %f \\pm %f\n",fitTrackEP->GetParameter(1),fitTrackEP->GetParError(1));
			printf("Found raw v4(track) = %f \\pm %f\n",fitTrackEP->GetParameter(2),fitTrackEP->GetParError(2));
			printf("Found raw v6(track) = %f \\pm %f\n",fitTrackEP->GetParameter(3),fitTrackEP->GetParError(3));

			printf("Using EPR2 = %f, EPR4 = %f, EPR6 = %f\n",fEPRes_R2,fEPRes_R4,fEPRes_R6);

			double raw_v2 = fitTrackEP->GetParameter(1);
			double raw_v4 = fitTrackEP->GetParameter(2);
			double raw_v6 = fitTrackEP->GetParameter(3);
			double raw_v2_un = fitTrackEP->GetParError(1);
			double raw_v4_un = fitTrackEP->GetParError(2);
			double raw_v6_un = fitTrackEP->GetParError(3);


			printf("Found raw v2(track) = %f \\pm %f\n",raw_v2,raw_v2_un);
			printf("Found raw v4(track) = %f \\pm %f\n",raw_v4,raw_v4_un);
			printf("Found raw v6(track) = %f \\pm %f\n",raw_v6,raw_v6_un);


			double v2 = fitTrackEP->GetParameter(1)/fEPRes_R2;
			double v4 = fitTrackEP->GetParameter(2)/fEPRes_R4;
			double v6 = fitTrackEP->GetParameter(3)/fEPRes_R6;



			double epsilon = 1e-5;

			double v2_un = 0;
			double v4_un = 0;
			double v6_un = 0;
			if (TMath::Abs(v2) > epsilon) {
				v2_un = v2 * TMath::Sqrt(TMath::Power(raw_v2_un/raw_v2,2)+TMath::Power(fEPRes_R2_Un/fEPRes_R2,2));
			}
			if (TMath::Abs(v4) > epsilon) {
				v4_un = v4 * TMath::Sqrt(TMath::Power(raw_v4_un/raw_v4,2)+TMath::Power(fEPRes_R4_Un/fEPRes_R4,2));
			}
			if (TMath::Abs(v6) > epsilon) {
				v6_un = v6 * TMath::Sqrt(TMath::Power(raw_v6_un/raw_v6,2)+TMath::Power(fEPRes_R6_Un/fEPRes_R6,2));
			}

			printf("Found v2(track) = %f \\pm %f\n",v2,v2_un);
			printf("Found v4(track) = %f \\pm %f\n",v4,v4_un);
			printf("Found v6(track) = %f \\pm %f\n",v6,v6_un);

			gTrack_ChiSq->SetPoint(i,(fMinPt+fMaxPt)/2.,fEPFitChi2OverNDF);
			gTrack_Bv->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(0));
			gTrack_V2->SetPoint(i,(fMinPt+fMaxPt)/2.,v2);
			gTrack_V4->SetPoint(i,(fMinPt+fMaxPt)/2.,v4);
			gTrack_V6->SetPoint(i,(fMinPt+fMaxPt)/2.,v6);

			gTrack_ChiSq->SetPointError(i,(fMaxPt-fMinPt)/2.,0.0);
			gTrack_Bv->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(0));
			gTrack_V2->SetPointError(i,(fMaxPt-fMinPt)/2.,v2_un);
			gTrack_V4->SetPointError(i,(fMaxPt-fMinPt)/2.,v4_un);
			gTrack_V6->SetPointError(i,(fMaxPt-fMinPt)/2.,v6_un);

			gTrack_V2_Raw->SetPoint(i,(fMinPt+fMaxPt)/2.,raw_v2);
			gTrack_V4_Raw->SetPoint(i,(fMinPt+fMaxPt)/2.,raw_v4);
			gTrack_V6_Raw->SetPoint(i,(fMinPt+fMaxPt)/2.,raw_v6);

			gTrack_V2_Raw->SetPointError(i,(fMaxPt-fMinPt)/2.,raw_v2_un);
			gTrack_V4_Raw->SetPointError(i,(fMaxPt-fMinPt)/2.,raw_v4_un);
			gTrack_V6_Raw->SetPointError(i,(fMaxPt-fMinPt)/2.,raw_v6_un);

			TLegend * fitBox= new TLegend(0.55,0.55,0.95,0.9);
			fitBox->AddEntry(hTrackEP,"Data","pe");
			fitBox->AddEntry(fitTrackEP,"b*(1+2*p_{2}*Cos(2x)+2p_{4}Cos(4x)+2p_{6}*Cos(6x))","l");
			fitBox->AddEntry((TObject *)0,Form("b = %.2e #pm %.2e",fitTrackEP->GetParameter(0),fitTrackEP->GetParError(0)),"");
			fitBox->AddEntry((TObject *)0,Form("p_{2}= %.2e #pm %.2e",fitTrackEP->GetParameter(1),fitTrackEP->GetParError(1)),"");
			fitBox->AddEntry((TObject *)0,Form("p_{4}= %.2e #pm %.2e",fitTrackEP->GetParameter(2),fitTrackEP->GetParError(2)),"");
			fitBox->AddEntry((TObject *)0,Form("p_{6}= %.2e #pm %.2e",fitTrackEP->GetParameter(3),fitTrackEP->GetParError(3)),"");
			fitBox->AddEntry((TObject *)0,Form("#chi^{2}/NDF = %.2f",fEPFitChi2OverNDF),"");
			fitBox->Draw("SAME");

			cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f.pdf",sTaskName.Data(),fMinPt,fMaxPt));
		}
		gTrack_Bv->Draw("ALP");
		ptName->Draw("SAME");
		cVn->Print(Form("%s/EPStudy_Track_B.png",sTaskName.Data()));
		

		TMultiGraph * mgVn = new TMultiGraph();
		mgVn->Add(gTrack_V2);
		mgVn->Add(gTrack_V4);
		//gTrack_V2->Draw("ALP");
		//gTrack_V4->Draw("LP SAME");
		mgVn->GetXaxis()->SetTitle(gTrack_V2->GetXaxis()->GetTitle());
		mgVn->Draw("ALP");
		cVn->BuildLegend();
		ptName->Draw("SAME");
		cVn->Print(Form("%s/EPStudy_Track_v2_v4.png",sTaskName.Data()));

		gTrack_ChiSq->Draw("ALP");
		ptName->Draw("SAME");
		cVn->Print(Form("%s/EPStudy_Track_ChiSq.png",sTaskName.Data()));




		// Study V3
		// need fEP3Res_R3, fEP3Res_R3_un


		TGraphErrors * gEP3Track_ChiSq = new TGraphErrors(kNTrackPtBins);
		gEP3Track_ChiSq->SetName("EP3Track_ChiSq");
		gEP3Track_ChiSq->SetTitle("ChiSquare/NDF from V_{n} fit");
		gEP3Track_ChiSq->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gEP3Track_ChiSq->GetYaxis()->SetTitle("#chi^{2}/NDF");

		TGraphErrors * gEP3Track_Bv = new TGraphErrors(kNTrackPtBins);
		gEP3Track_Bv->SetName("EP3Track_Bv");
		gEP3Track_Bv->SetTitle("B Value from V_{n} fit");
		gEP3Track_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		gEP3Track_Bv->GetYaxis()->SetTitle("B");


		TGraphErrors * gEP3Track_V3 = new TGraphErrors(kNTrackPtBins);
		gEP3Track_V3->SetName("EP3Track_V3");
		gEP3Track_V3->SetTitle("Calculated v_{3}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
		gEP3Track_V3->SetMarkerStyle(kOpenSquare);
		TGraphErrors * gEP3Track_V3_Raw = (TGraphErrors *) gEP3Track_V3->Clone("EP3Track_V3_Raw");
		gEP3Track_V3_Raw->SetTitle("Calculated v_{3}^{Track} (Event Plane method, No EPR Correction))");
		gEP3Track_V3_Raw->GetYaxis()->SetTitle("v_{3}(raw)");



		printf("Starting the fitting of tracks vs the 3rd order event plane\n");
		for (int i = 0; i < kNTrackPtBins; i++) {
			printf("Fitting track pt bin %d\n",i);
			double fMinPt = fTrackPtBins[i];
			double fMaxPt = fTrackPtBins[i+1];

			// FIXME
			TH1F * hTrackEP3Local = 0;
			hTrackEP3Local = hPtEP3AngleTrack_Proj[i];
			if (!hTrackEP3Local) {
				printf("Missing the track vs EP3 histogram\n");
				return;
			}


			TF1 * fitTrackEP3 = new TF1(Form("EP3Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(3*x))",hTrackEP3Local->GetXaxis()->GetXmin(),hTrackEP3Local->GetXaxis()->GetXmax());
			printf("Made a fit  function. Setting initial parameters ...\n");
			cout << "\tFitting using histogram " << hTrackEP3Local->GetName() <<endl;
			cout << "\tand fit function " << fitTrackEP3->GetName() <<endl;
			fitTrackEP3->SetParameter(0,hTrackEP3Local->Integral("width") / (TMath::Pi() / 2));
			fitTrackEP3->SetParameter(1,0.03);

			cout << "Setting par limits" <<endl;
			fitTrackEP3->SetParLimits(1,0.0,0.5);

			cout << "Setting par names" <<endl;
			fitTrackEP3->SetParName(0,"B");
			fitTrackEP3->SetParName(1,"v_3");
			cout << "About to fit" <<endl;

			if (hTrackEP3Local->GetEntries() == 0) {
				printf("Skipping pt bin because of empty histogram\n");
				continue;
			}

			hTrackEP3Local->Fit(fitTrackEP3);
			fitTrackEP3->SetLineColor(kCyan);

			hTrackEP3Local->SetMarkerStyle(kFullSquare);
			hTrackEP3Local->Draw();
			fitTrackEP3->Draw("SAME");

			double fEPFitChiSquare = fitTrackEP3->GetChisquare();
			double fEPFitNDF = fitTrackEP3->GetNDF();
			double fEPFitChi2OverNDF = 0;
			if (fEPFitNDF > 0) {
				fEPFitChi2OverNDF = fEPFitChiSquare / fEPFitNDF;
			}
			printf("Found raw v2(track) = %f \\pm %f\n",fitTrackEP3->GetParameter(1),fitTrackEP3->GetParError(1));
			printf("Using EPR3 = %f\n",fEP3Res_R3);
			double raw_v3 = fitTrackEP3->GetParameter(1);
			double raw_v3_un = fitTrackEP3->GetParError(1);
			double v3 = fitTrackEP3->GetParameter(1)/fEP3Res_R3;

			double epsilon = 1e-5;

			double v3_un = 0;
			if (TMath::Abs(v3) > epsilon) {
				v3_un = v3 * TMath::Sqrt(TMath::Power(raw_v3_un/raw_v3,2)+TMath::Power(fEP3Res_R3_Un/fEP3Res_R3,2));
			}
			printf("Found v3(track) = %f \\pm %f\n",v3,v3_un);
			gEP3Track_ChiSq->SetPoint(i,(fMinPt+fMaxPt)/2.,fEPFitChi2OverNDF);
			gEP3Track_Bv->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP3->GetParameter(0));
			gEP3Track_V3->SetPoint(i,(fMinPt+fMaxPt)/2.,v3);
			gEP3Track_ChiSq->SetPointError(i,(fMaxPt-fMinPt)/2.,0.0);
			gEP3Track_Bv->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP3->GetParError(0));
			gEP3Track_V3->SetPointError(i,(fMaxPt-fMinPt)/2.,v3_un);
			gEP3Track_V3_Raw->SetPoint(i,(fMinPt+fMaxPt)/2.,raw_v3);
			gEP3Track_V3_Raw->SetPointError(i,(fMaxPt-fMinPt)/2.,raw_v3_un);

			TLegend * fitBox= new TLegend(0.55,0.55,0.95,0.9);
			fitBox->AddEntry(hTrackEP3Local,"Data","pe");
			fitBox->AddEntry(fitTrackEP3,"b*(1+2*p_{3}*Cos(3x))","l");
			fitBox->AddEntry((TObject *)0,Form("b = %.2e #pm %.2e",fitTrackEP3->GetParameter(0),fitTrackEP3->GetParError(0)),"");
			fitBox->AddEntry((TObject *)0,Form("p_{3}= %.2e #pm %.2e",fitTrackEP3->GetParameter(1),fitTrackEP3->GetParError(1)),"");
			fitBox->AddEntry((TObject *)0,Form("#chi^{2}/NDF = %.2f",fEPFitChi2OverNDF),"");
			fitBox->Draw("SAME");

			cVn->Print(Form("%s/EP3Study_Track_Pt_%.2f_%.2f.pdf",sTaskName.Data(),fMinPt,fMaxPt));
		}



		gEP3Track_V3->Draw("ALP");
		cVn->BuildLegend();
		ptName->Draw("SAME");
		cVn->Print(Form("%s/EP3Study_Track_v3.png",sTaskName.Data()));



		// If CorrVsManyThings Exists, show some projections
		THnSparseF * CorrVsManyThings = 0;

		TH1D * hCorrDeltaPhi = 0;
		TH1D * hCorrDeltaEta = 0;
		TH1D * hCorrEP = 0;
		TH1D * hCorrCent = 0;

		CorrVsManyThings = (THnSparseF *) fTask->FindObject("CorrVsManyThings");
		if (CorrVsManyThings != 0) {
			cout << "Found CorrVsManyThings THnSparse. Analyzing for QA." << endl;

			hCorrDeltaPhi = (TH1D *) CorrVsManyThings->Projection(0,"e");
			hCorrDeltaPhi->Draw();
			gPad->Print(Form("%s/Corr_DeltaPhi.png",sTaskName.Data()));

			hCorrDeltaEta = (TH1D *) CorrVsManyThings->Projection(1,"e");
			hCorrDeltaEta->Draw();
			gPad->Print(Form("%s/Corr_DeltaEta.png",sTaskName.Data()));

			hCorrEP = (TH1D *) CorrVsManyThings->Projection(6,"e");
			hCorrEP->Draw();
			gPad->Print(Form("%s/Corr_EP.png",sTaskName.Data()));

			hCorrCent = (TH1D *) CorrVsManyThings->Projection(7,"e");
			hCorrCent->Draw();
			gPad->Print(Form("%s/Corr_Cent.png",sTaskName.Data()));
		}

		// If TriggerHist Exists, show some projections
		THnSparseF * TriggerHist = 0;
		TH1D * Trigger_pT = 0;
		TH1D * Trigger_EP = 0;
		TH1D * Trigger_Cent = 0;



		TriggerHist = (THnSparseF *) fTask->FindObject("TriggerHist");
		if (TriggerHist != 0) {
			cout << "Found TriggerHist THnSparse. Analyzing for QA." << endl;
			Trigger_pT = (TH1D *) TriggerHist->Projection(0,"e");
			Trigger_pT->Draw();
			gPad->Print(Form("%s/Trigger_pT.png",sTaskName.Data()));

			Trigger_EP = (TH1D *) TriggerHist->Projection(2,"e");
			Trigger_EP->Draw();
			gPad->Print(Form("%s/Trigger_EP.png",sTaskName.Data()));

			Trigger_Cent = (TH1D *) TriggerHist->Projection(3,"e");
			Trigger_Cent->Draw();
			gPad->Print(Form("%s/Trigger_Cent.png",sTaskName.Data()));
		}



		// If Pi0Cands exists, show some projections
		THnSparse * Pi0Cands = 0;

		TH1D * hPi0CandPt = 0;
		TH1D * hPi0CandMass = 0;
		TH1D * hPi0CandType = 0;
		TH1D * hPi0CandEPAngle = 0;


		Pi0Cands = (THnSparse *) fTask->FindObject("Pi0Cands");
		if (Pi0Cands != 0) {
			cout << "Found Pi0Cands THnSparse. Analyzing for QA." << endl;

			hPi0CandPt = (TH1D *) Pi0Cands->Projection(0,"e");
			hPi0CandPt->Draw();
			gPad->Print(Form("%s/Pi0Cand_pT.png",sTaskName.Data()));

			hPi0CandMass = (TH1D *) Pi0Cands->Projection(1,"e");
			hPi0CandMass->Draw();
			gPad->Print(Form("%s/Pi0Cand_mass.png",sTaskName.Data()));

			hPi0CandType = (TH1D *) Pi0Cands->Projection(6,"e");
			hPi0CandType->Draw();
			gPad->Print(Form("%s/Pi0Cand_Type.png",sTaskName.Data()));

			hPi0CandEPAngle = (TH1D *) Pi0Cands->Projection(7,"e");
			hPi0CandEPAngle->Draw();
			gPad->Print(Form("%s/Pi0Cand_EPAngle.png",sTaskName.Data()));
		}


		fQAOutput->Add(gTrack_ChiSq);
		fQAOutput->Add(gTrack_Bv);
		fQAOutput->Add(gTrack_V2);
		fQAOutput->Add(gTrack_V4);
		fQAOutput->Add(gTrack_V6);
		fQAOutput->Add(gTrack_V2_Raw);
		fQAOutput->Add(gTrack_V4_Raw);
		fQAOutput->Add(gTrack_V6_Raw);

		fQAOutput->Add(gEP3Track_ChiSq);
		fQAOutput->Add(gEP3Track_Bv);
		fQAOutput->Add(gEP3Track_V3);
		fQAOutput->Add(gEP3Track_V3_Raw);

		fQAOutput->Write();
		fQAOutput->Close();

	}

}
