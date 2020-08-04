
Bool_t fActuallyDoFinalFit = true;

TH1D * hGammaE;
TH1D * hPairTheta;
std::vector<TH1D *> MassBkgArray;

void PrintParameters(TF1 * func, bool includeErrors = true, std::string label = "Parameters") {

  int nPar = func->GetNpar();
  double parMin;
  double parMax;
  printf("======================================================================\n");
  printf("||  %s\n",label.c_str());
  printf("======================================================================\n");
  const int nBuff = 15;

  for (int i = 0; i < nPar; i++) {
    func->GetParLimits(i,parMin,parMax);
    const char *parName = func->GetParName(i);
    char parNameBuffer[nBuff] = "              ";
    strncpy(parNameBuffer,parName,min((int) strlen(parName),nBuff));
    int parNameLen = strlen(parName);
    if (includeErrors) {
      printf("%s:\t[ %.3f | \t%e \\pm %e \t | %.3f ]  \n",parNameBuffer,parMin,func->GetParameter(i),func->GetParError(i),parMax);
    //  printf("%s:\t\t\t[%.2f, %.2f \\pm %.2f  ,%.2f]  \n",func->GetParName(i),parMin,func->GetParameter(i),func->GetParError(i),parMax);
    } else {
      printf("%s:\t[ %.3f | \t%e \t |%.3f ]  \n",parNameBuffer,parMin,func->GetParameter(i),parMax);
//      printf("%s:\t[ %.3f | \t%.3f \t |%.3f ]  \n",parNameBuffer,parMin,func->GetParameter(i),parMax);
    }
  }
  if (includeErrors) {
      // If errors are included, also show ChiSq, NDF
      printf("Chi^{2} =  %.3f\t\tNDF = %d\n",func->GetChisquare(),func->GetNDF());
  }
  printf("======================================================================\n");
  return;
}


Double_t testFunction(Double_t * x, Double_t * par) {
  Double_t y = 0;
  Double_t m = x[0];
  
  Double_t Y = par[0];
  Double_t mu = par[1];
  Double_t sigma = par[2];


  y =  Y * TMath::Gaus(m,mu,sigma,1);


  return y;
}

// Zero Out a selected region of a histogram
void zeroRegion(TH1D * hist, double min, double max) {
  int minBin = hist->FindBin(min);
  int maxBin = hist->FindBin(max);
  for (int i = minBin; i < maxBin; i++) {
    hist->SetBinContent(i,0); 
    // do something with error?
    hist->SetBinError(i,0);
  }
}

// Functor for Pi0 Fit with free parameter for background
// 1 Gaus for Peak. Poly(2) for Residual, free scale for CombBkg
class ParamBkg_Functor {
	public:
		ParamBkg_Functor() : nPeakParams(3),nBkgParams(0),fEtaIndex(-1),fPeakChoice(-1),fBkgChoice(-1),fUseEta(0),fBkgHist(0),fThetaModelLambda(-1),fThetaModelMPrime(-1) {

		}
		ParamBkg_Functor(TH1D * BkgHist, bool useEta, int iPeakChoice = 0, int iBkgChoice = 4, double ThetaModelLambda = 0., double ThetaModelMPrime = 0.) {
			//fBkgHist = BkgHist;
			fBkgHist = 0;
			if (BkgHist != 0) {
				fBkgHist = (TH1D * ) BkgHist->Clone(Form("%s_Clone",BkgHist->GetName()));
			}
			fThetaModelLambda = ThetaModelLambda;
			fThetaModelMPrime = ThetaModelMPrime;

			fUseEta = useEta;

			fPeakChoice = iPeakChoice;
			fBkgChoice  = iBkgChoice;

			nPeakParams = 0;
			nBkgParams  = 0; 

			switch (fPeakChoice) {
				case 7: // Voigt Profile
					nPeakParams = 4;
					break;
				case 6: // ExpGausExp (2 Sides)
					nPeakParams = 5;
					break;
				case 4: // Crystal Ball
				case 3: // Crystal Ball
					nPeakParams = 5;
					break;
				case 2: // B-W
					nPeakParams = 3;
					break;
				case 5: // ExpGaus
				case 1: // ExpGaus
					nPeakParams = 4; 
					break;
				case 0:
				default:
					nPeakParams = 3;
			}
			// 
			switch (fBkgChoice) {
				case 8: 
					nBkgParams = 5;
					break;
				case 7:
					nBkgParams = 4;			
					break;
				case 6: 
					nBkgParams = 6;
					break;		
				case 5: 
					nBkgParams = 5;
					break;		
				case 4: 
					nBkgParams = 5;
					break;		
				case 3: 
					nBkgParams = 4;
					break;		
				case 2: 
					nBkgParams = 3;
					break;		
				case 1: 
					nBkgParams = 2;
					break;		
				case 0:
				default:
					nBkgParams = 1;
			}
			if (fBkgHist != 0) {
				printf("Since we have a bkg histogram, increase nBkgParams by one.\n");
				nBkgParams += 1;
			} else {
				printf("Since we don't have a bkg histogram, don't increase nBkgParams.\n");
			}
			fEtaIndex = nPeakParams + nBkgParams; // Start of indices for Eta peak parameters


			printf("TEST: nParams     = %d\n",GetNParams());
			printf("TEST: nPeakParams = %d\n",nPeakParams);
			printf("TEST: nBkgParams  = %d\n",nBkgParams);
			printf("TEST: iBkgScale   = %d\n",GetBkgScaleIndex());
			printf("TEST: etaIndex    = %d\n",fEtaIndex);
			printf("Test: GetEtaIndex = %d\n",GetEtaIndex());
		}
		double operator() (double *x, double *p) {
			double value = 0; 

			// For readability, use these variables
			double mu,sigma,alpha,n;
			double A,B,C,D,N;

			// Start with peak choice
			switch (fPeakChoice) {
				case 8: // Sum of Two Gaussians
					value += p[0] * TMath::Gaus(x[0],p[1],p[2],1) + p[3] * TMath::Gaus(x[0],p[4],p[5],1);
					break;
				case 7: // Voigt Profile 
					// Last arg to Voigt is resolution of Voigt approximationg
					value += p[0] * TMath::Voigt(x[0]-p[1],p[2],p[3],4); 
					break;
				case 6: // Gaus with 2 ExpDecays (ExpGaussExp)
					value += p[0]*(TMath::Gaus(x[0],p[1],p[2],1)+TMath::Exp((x[0]-p[1])/p[3])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] > x[0]) + TMath::Exp(-(x[0]-p[1])/p[4])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] < x[0])); 
					// Left
					//value += p[0]*(TMath::Gaus(x[0],p[1],p[2],1)+TMath::Exp((x[0]-p[1])/p[3])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] > x[0])); 
					// Right
					//value += p[0]*(TMath::Gaus(x[0],p[1],p[2],1)+TMath::Exp(-(x[0]-p[1])/p[3])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] < x[0])); 
					break;

				case 5: // GausWithExpDecay (GaussExp Right)
					value += p[0]*(TMath::Gaus(x[0],p[1],p[2],1)+TMath::Exp(-(x[0]-p[1])/p[3])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] < x[0])); 
					break;

				case 4:
				case 3: // Crystal Ball
					mu = p[1];
					sigma = p[2];
					alpha = TMath::Abs(p[3]);
					n = p[4];

					A = TMath::Power( n / alpha, n) * TMath::Exp(-alpha*alpha / 2);
					B = n / alpha - alpha;
					C = ((n / alpha) / (n - 1)) * TMath::Exp(-alpha*alpha / 2); // Need alpha > 0, n > 1
					D = TMath::Sqrt(TMath::PiOver2()) * (1 + TMath::Erf(alpha/TMath::Sqrt(2))) ;
					N = 1 / (sigma * (C + D)); // Need sigma > 0, C+D > 0

					if (fPeakChoice == 3) {
						value += p[0] * N * (TMath::Gaus(x[0],p[1],p[2],0)*(x[0] - mu >= -sigma * alpha)+A*TMath::Power(B - (x[0]-mu)/sigma,-n) * (x[0] - mu < -sigma * alpha));
					} else if (fPeakChoice == 4) {
						value += p[0] * N * (TMath::Gaus(x[0],p[1],p[2],0)*(x[0] - mu <= sigma * alpha)+A*TMath::Power(B + (x[0]-mu)/sigma,-n) * (x[0] - mu > sigma * alpha));
					}
			
			//		value += p[0]*(TMath::Gaus(x[0],p[1],p[2],0)*(x[0] < p[3])+()*(x[0] >= p[3])); 
					break;
				case 2: // BreitWigner
					value += p[0]*TMath::BreitWigner(x[0],p[1],p[2]*TMath::Sqrt(2.*TMath::Log(2)));
					break;
				case 1: // GausWithExpDecay (GaussExp Left)
					value += p[0]*(TMath::Gaus(x[0],p[1],p[2],1)+TMath::Exp((x[0]-p[1])/p[3])*(TMath::Gaus(0,0,p[2],1) - TMath::Gaus(x[0],p[1],p[2],1))*(p[1] > x[0])); 
					break;
				case 0: // Gaussian
				default:
					value += p[0] * TMath::Gaus(x[0],p[1],p[2],1);
			} 

			double bkgValue = 0;
	
			switch (fBkgChoice) {
				// ExpDecay * Poly(2)
				case 8: 
					bkgValue += (x[0] > p[nPeakParams+2])*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(p[nPeakParams]+p[nPeakParams+3]*x[0] + 0.5 * p[nPeakParams+4]*x[0]*x[0]);
// Why did this miss squared term					bkgValue += (x[0] > p[nPeakParams+2])*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(p[nPeakParams]+p[nPeakParams+3]*x[0] + 0.5 * p[nPeakParams+4]*x[0]*x[0]*x[0]);
//					bkgValue += (1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(p[nPeakParams]+p[nPeakParams+3]*x[0] + 0.5 * p[nPeakParams+4]*x[0]*x[0]*x[0]);
////bad					bkgValue += (1 - TMath::Exp(-(x[0] - p[nPeakParams+1])*TMath::Power(p[nPeakParams],2)))*(p[nPeakParams+2]+p[nPeakParams+3]*x[0] + p[nPeakParams+4]*x[0]*x[0]);
//					bkgValue += p[nPeakParams]*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(1+p[nPeakParams+3]*x[0] + p[nPeakParams+4]*x[0]*x[0]);
					break;
				// ExpDecay * Poly(1)
				case 7: 
					bkgValue += (x[0] > p[nPeakParams+2])*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(p[nPeakParams]+p[nPeakParams+3]*x[0]);
//					bkgValue += (1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(p[nPeakParams]+p[nPeakParams+3]*x[0]);
//					bkgValue += p[nPeakParams]*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])*TMath::Power(p[nPeakParams+1],2)))*(1+p[nPeakParams+3]*x[0]);
//					bkgValue += p[nPeakParams]*(1 - TMath::Exp(-(x[0] - p[nPeakParams+2])/p[nPeakParams+1]))*(1+p[nPeakParams+3]*x[0]);
					break;
				case 6:
					bkgValue += p[nPeakParams]*TMath::Exp(-p[nPeakParams+1]*x[0])*(p[nPeakParams+2]+p[nPeakParams+3]*x[0]+(1./2.)*p[nPeakParams+4]*x[0]*x[0]+(1./6.)*p[nPeakParams+5]*x[0]*x[0]*x[0]);
					break;
				case 5:
					bkgValue += p[nPeakParams]*TMath::Exp(-p[nPeakParams+1]*x[0])*(p[nPeakParams+2]+p[nPeakParams+3]*x[0]+(1./2.)*p[nPeakParams+4]*x[0]*x[0]);
					break;
				case 4:
					bkgValue += (1./24.) * p[nPeakParams + 4] * TMath::Power(x[0],4);
				case 3: // For case 0-3, don't use "break;"; take advantage of case structure
					bkgValue += (1./6.) * p[nPeakParams + 3] * x[0] * x[0] * x[0];
				case 2: 
					bkgValue += 0.5 * p[nPeakParams + 2] * x[0] * x[0];
				case 1:
					bkgValue += p[nPeakParams + 1] * x[0];	
				default:
				case 0: // Poly(0)
					bkgValue += p[nPeakParams];
			}

			if (bkgValue > 0.) value += bkgValue;
			if (bkgValue < 0.) { 
				// Mode one: "punish" the fit
////				 value += TMath::Power(TMath::Abs(bkgValue),5); 
//				 value += (bkgValue*1e6); 
				// Mode two: don't add ;

			}

			// Adding Background term:
			if (fBkgHist) {
				value += p[nPeakParams + nBkgParams - 1] * fBkgHist->Interpolate(x[0]);
			}	// nPeakParams + nBkgParams is index of bkg scale

			if (fUseEta) value += p[fEtaIndex] * TMath::Gaus(x[0],p[fEtaIndex + 1],p[fEtaIndex + 2],0); 
//			if (fUseEta) value += p[fEtaIndex] * TMath::Gaus(x[0],p[fEtaIndex + 1],p[fEtaIndex + 2],1); 
			// Appling Opening Angle Cut model
			if (fThetaModelLambda != 0) {
				value = value * TMath::Erf(fThetaModelLambda*(1e3*x[0] - fThetaModelMPrime)) * (1e3*x[0] - fThetaModelMPrime > 0);
			}

			return value;
		}
		int GetNParams()       { return nPeakParams + nBkgParams + fUseEta * 3; }
		int GetNPeakParams()   { return nPeakParams; }
		int	GetNBkgParams()    { return nBkgParams; }
		int GetEtaIndex()     { if (fUseEta) return fEtaIndex; else return -1; }
		int GetBkgScaleIndex() { if(HasBkg()) return nPeakParams + nBkgParams - 1; else return -1; }

		double GetLambda() { return fThetaModelLambda; }
		double GetMPrime() { return fThetaModelMPrime; }
		void SetLambda(double input) {fThetaModelLambda = input; }
		void SetMPrime(double input) {fThetaModelMPrime = input; }

		TF1 * GetFunc()          { return fFunc; }
		void SetFunc(TF1 * Func) { fFunc = Func; }
		bool HasBkg()            { return fBkgHist!=0; }

	private:
		int nPeakParams;
		int nBkgParams;
		int fEtaIndex;

		TF1 * fFunc;

		int fPeakChoice;
		// 0  --  Gaus
		// 1  --  ExpDecay-Gaus 
		// 2  --  Breit-Wigner
		// 3  --  Crystal Ball
		// 3  --  Crystal Ball (Right Side)
		int fBkgChoice; // Choice of Residual Fit Background
		// 0  --  Poly(0)
		// 1  --  Poly(1)
		// 2  --  Poly(2)
		// 3  --  Poly(3)
		// 4  --  Exp * Poly(2)
		// 5  --  Exp * Poly(3)
		// 6  --  ExpDecay * Poly(2)
		// 7  --  ExpDecay * Poly(3)

		bool fUseEta;
		TH1D * fBkgHist;
//		TH2F * fThetaModelParams;
		double fThetaModelLambda;
		double fThetaModelMPrime;
};

void PreFitPi0(TF1 * fit, ParamBkg_Functor * fitFunct, TH1D * hist, bool useEta, int iPeakChoice, int iBkgChoice, bool bHaveMCPreFit = false) ;

//TF1 * fitPi0Peak(TH1D * hist, std::string name, ParamBkg_Functor ** Func, bool useEta = false, int iPeakChoice = 0, int iBkgChoice = 4, TH1D * hBkgHist = 0) {
//TF1 * fitPi0Peak(TH1D * hist, std::string name, bool useEta = false, int iPeakChoice = 0, int iBkgChoice = 4, TH1D * hBkgHist = 0) {
//TF1 * fitPi0Peak(TH1D * hist, std::string name, ParamBkg_Functor ** Func, bool useEta = false, int iPeakChoice = 0, int iBkgChoice = 4, TH1D * hBkgHist = 0, double fThetaModelLambda = 0., double fThetaModelMPrime = 0., bool enableFit = true) {

/** The big pi0 fitting function using the functor and the modular functions
  *
  */
TF1 * fitPi0Peak(TH1D * hist, std::string name, ParamBkg_Functor ** Func, bool useEta = false, int iPeakChoice = 0, int iBkgChoice = 4, TH1D * hBkgHist = 0, double fThetaModelLambda = 0., double fThetaModelMPrime = 0., bool enableFit = true, TF1 * fMCPreFit = 0) {
	printf("Starting Fit for %s, with",name.c_str());
	if (!useEta) printf("out");
	printf(" Eta Peak, using peak function %d and bkg function %d.\n",iPeakChoice,iBkgChoice);

/*	if (!hist->GetEntries()) {
		printf("Skipping empty hist\n");
		TF1 * fEmptyFit = new TF1(name.c_str(),);
		return 0;
	}*/

	if (hBkgHist != 0) {
		printf("Using background histogram (name,title) = (%s,%s)\n",hBkgHist->GetName(),hBkgHist->GetTitle());
	}

	bool bHaveMCPreFit = (fMCPreFit != 0);

	if (bHaveMCPreFit) {
		printf("Received Pre-Analyzed MC Pi0 Fit %s\n",fMCPreFit->GetName());

	}

	int n = 1;
	while (hist->GetBinContent(n) == 0) {
		n++; // Finding first nonempty bin
		if (n >= hist->GetNbinsX() ) {
			printf("Empty Hist! Returning empty function\n");
			TF1 * emptyFunction = new TF1(name.c_str(),"[0]",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
			emptyFunction->SetParameter(0,0.0);
			return emptyFunction;
		}
	}
	

	// FIXME check that function is defined only on fitting range??
  //double minX = 0.018+hist->GetXaxis()->GetBinLowEdge(n);
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
//  double minX = 0.01+hist->GetXaxis()->GetBinLowEdge(n);
//	double maxX = 0.45;
//	double maxX = 0.7;
//  double maxX = 0.75;
	// FIXME make a way to vary this
	double maxX = hist->GetXaxis()->GetXmax();
	printf("DEBUG MinX = %f (n=%d), MaxX = %f\n",minX,n,maxX);

	//ParamBkg_Functor * fitFunct = new ParamBkg_Functor(hBkgHist,useEta,iPeakChoice,iBkgChoice);
	// FIXME Add in ThetaModel parameters here
	ParamBkg_Functor * fitFunct = new ParamBkg_Functor(hBkgHist,useEta,iPeakChoice,iBkgChoice,fThetaModelLambda,fThetaModelMPrime);
//	ParamBkg_Functor * fitFunct = new ParamBkg_Functor(0,useEta,iPeakChoice,iBkgChoice);

	*Func = fitFunct;

	int nParam = fitFunct->GetNParams();

	// Not necessary to define these variables, but useful for readability and coding.
	int iPi0_Yield = 0;  // Always true
	int iPi0_Mean  = 1;
	int iPi0_Sigma = 2;
	
	int iBkg_Par_0 = fitFunct->GetNPeakParams(); // First Bkg parameter
	int iBkg_Par_BScale = fitFunct->GetBkgScaleIndex(); //Index of background scale parameter

	int iEta_Yield = fitFunct->GetEtaIndex();
	int iEta_Mean  = iEta_Yield + 1;
	int iEta_Sigma = iEta_Yield + 2;

	TF1 * fit = new TF1(name.c_str(),fitFunct,minX,maxX,nParam);

	// Label Parameters
	fit->SetParName(iPi0_Yield,"Y");
	fit->SetParName(iPi0_Mean,"mu");
	fit->SetParName(iPi0_Sigma,"sigma");

	// Most parameter setting is in the PreFit
	// Will move parameter setting and limits to here for tail parameters 
	switch (iPeakChoice) {
		case 7: // Voigt Profile
			fit->SetParName(3,"gamma"); //should be 49 eV in theory? or 2.528e-8 eV?
			fit->SetParLimits(3,0.001,0.05);
			fit->SetParameter("gamma",0.005);	
			break;
		case 4: // Crystal Ball (Right Side)
		case 3: // Crystal Ball
			printf("DEBUG Setting up Crystal Ball parameters\n");
			fit->SetParName(3,"alpha");
			fit->SetParName(4,"n");

			fit->SetParLimits(3,0.01,4);
			fit->SetParLimits(4,1.05,6);

			fit->SetParameter(3,1);
			fit->SetParameter(4,2.);
			printf("DEBUG  ... Finished\n");

			break;
		case 2: // Breit-Wigner
			break;
		case 6: // 2-sided ExpDecay-Gaus (GausExp 2 Sides)
			fit->SetParName(3,"lambda_1");
			fit->SetParameter(3,0.020); //0.017
			fit->SetParLimits(3,0.015,0.035); //0.005,0.035
			fit->SetParName(4,"lambda_2");
			fit->SetParameter(4,0.017);
			fit->SetParLimits(4,0.005,0.06);
			break;
		case 5: // ExpDecay-Gaus (GausExp Right)
		case 1: // ExpDecay-Gaus (GaussExp)
		//fit->SetParName(3,"alpha");
			fit->SetParName(3,"lambda");
			fit->SetParLimits(3,0.005,0.035); //0.1
			fit->SetParameter("lambda",0.015);
			break;
		case 0:
		default:
			break;
	}

	// Fix Parameters to values from MC, if available
	if (bHaveMCPreFit) {
		// fMCPreFit
		
		// Check that MC PreFit has same number of peak parameters
		if (fitFunct->GetNPeakParams() != fMCPreFit->GetNpar() - 1) { // -1 for pedestal fit parameter
			fprintf(stderr,"Error: mismatch in number of parameters in MC Pi0 fit\n");
		}
		int nShapeParameters = fitFunct->GetNPeakParams() - 3; // minus Y,mu,sigma
		printf("Loading %d parameters from MC Fits\n",nShapeParameters);
		for (int i = 0; i < nShapeParameters; i++) {
			printf("Setting parameter %d to MC value %f\n",i+3,fMCPreFit->GetParameter(i+3));
			fit->FixParameter(i+3,fMCPreFit->GetParameter(i+3));
		}
	}


	fit->SetParName(iBkg_Par_0,"N"); // The first background parameter
	fit->SetParName(iBkg_Par_BScale,"Bkg_Scale");


	if (iBkgChoice == 5 || iBkgChoice == 6) {
		fit->SetParName(iBkg_Par_0 + 1,"gamma");
	}

	if (iBkgChoice == 7 || iBkgChoice == 8) { //Optional names
		fit->SetParName(iBkg_Par_0 + 1, "gamma");
		fit->SetParName(iBkg_Par_0 + 2, "d");
		fit->SetParName(iBkg_Par_0 + 3, "B");
		if (iBkgChoice == 8) fit->SetParName(iBkg_Par_0 + 4,"C");
	}


	if (useEta) {
		fit->SetParName(iEta_Yield,"Y_eta");
		fit->SetParName(iEta_Mean,"mu_eta");
		fit->SetParName(iEta_Sigma,"sigma_eta");
	}

//	printf("Initializing parameters for %d peak parameters and %d background parameters,",fitFunct->GetNPeakParam(),fitFunct->GetNBkgParam());
//	if (useEta) {
//		printf(" and using a peak for the eta.\n");
//	} else {
//		printf(" and not using a peak for the eta.\n");
//	}
//	printf("Number of Parameters: %d\n",fitFunct->GetNParams());

	if (enableFit) {
		PreFitPi0(fit, fitFunct, hist, useEta, iPeakChoice, iBkgChoice, bHaveMCPreFit); // Can make alternative PreFit methods
	}

	PrintParameters(fit,false,"Initial Parameters");
	printf("Initiating Final Fit...\n");
	
	if (fActuallyDoFinalFit && enableFit) {
		hist->Fit(fit,"0M","",minX,maxX);
////	hist->Fit(fit,"RQ0");
		printf("Done Fitting %s\n",hist->GetName());
	} else {
		printf("Final Fit was not done.\n");
	}
	PrintParameters(fit,true,"Final Parameters");

	return fit;
}

//void PreFitPi0(TF1 * fit, ParamBkg_Functor * fitFunct, TH1D * hist, bool useEta, int iPeakChoice, int iBkgChoice) {
void PreFitPi0(TF1 * fit, ParamBkg_Functor * fitFunct, TH1D * hist, bool useEta, int iPeakChoice, int iBkgChoice, bool bHaveMCPreFit = false) {
	printf("Initializing Pre-Fit ... \n");
	if (bHaveMCPreFit) printf("  Using MC Prefit for tail parameters.\n");
// This is duplicated for convenience ?? ... 
	// Not necessary to define these variables, but useful for readability and coding.
	int iPi0_Yield = 0;  // Always true
	int iPi0_Mean  = 1;
	int iPi0_Sigma = 2;
	
	int nBkg_Par = fitFunct->GetNBkgParams();
	int iBkg_Par_0 = fitFunct->GetNPeakParams(); // First Bkg parameter
	int iBkg_Par_BScale = fitFunct->GetBkgScaleIndex(); // Index of bkg scale

	int iEta_Yield = fitFunct->GetEtaIndex();
	int iEta_Mean  = iEta_Yield + 1;
	int iEta_Sigma = iEta_Yield + 2;
//

	// Finding first bin with data content
	int iFirstUsedBin = 1;
	int nBins = hist->GetXaxis()->GetNbins();
	while ((hist->GetBinContent(iFirstUsedBin) == 0) && (iFirstUsedBin < nBins+1)) iFirstUsedBin++;


	double minX = fit->GetXmin();
	double maxX = fit->GetXmax();
 
	double fPi0_Yield_Guess = 0.1; 
	double fPi0_Mean_Guess  = 0.135;
	double fPi0_Sigma_Guess = 0.010; //0.007

  double etaRangeMin = 0.5;
  double etaRangeMax = 0.65; //0.7;

	double fEta_Yield_Guess = 0.1;
	double fEta_Mean_Guess = 0.55;
	double fEta_Sigma_Guess = 0.03; //0.04;

	double maxEtaInt = 0;

	double meanHeight = hist->Integral(hist->FindBin(minX), hist->FindBin(maxX),"width") / (maxX - minX);
	double maxHeight  = hist->GetBinContent(hist->GetMaximumBin());
	double minHeight  = hist->GetBinContent(hist->GetMinimumBin());
	float  fMassRange = hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin(); 

	// very rough maximum for slopes of residual background
	double bkg_Slope_Min = (minHeight - maxHeight) / fMassRange;
	double bkg_Slope_Max = -bkg_Slope_Min;


	
	fit->SetParameter(iBkg_Par_0,meanHeight); //May need to change this later
	fit->SetParLimits(iBkg_Par_0,0.,4*maxHeight); 
//	fit->SetParameter(iBkg_Par_BScale,meanHeight); //May need to change this later
//	fit->SetParLimits(iBkg_Par_BScale,0.,4*maxHeight); 
//	fit->SetParameter(3,meanHeight); //May need to change this later
//	fit->SetParLimits(3,0.,4*maxHeight); 

	if (useEta) {
		maxEtaInt = hist->Integral(hist->FindBin(etaRangeMin),hist->FindBin(etaRangeMax),"width");
	
		fEta_Yield_Guess = 0.3 * maxEtaInt;
	}

	double fPi0_Yield_Min = 0;
	double fPi0_Yield_Max = meanHeight * (maxX - minX);
	double fPi0_Mean_Min = 0.1;
	double fPi0_Mean_Max = 0.2;
	double fPi0_Sigma_Min = 0.005; //0.003
	double fPi0_Sigma_Max = 0.019;

	double Bkg_BScale_Min = 0.7;
	double Bkg_BScale_Max = 1.0;

	double fEta_Yield_Min = 0;
	double fEta_Yield_Max = maxEtaInt;
	double fEta_Mean_Min = 0.54;
	double fEta_Mean_Max = 0.6;
	double fEta_Sigma_Min = 0.01;
	double fEta_Sigma_Max = 0.045;

	bool fUsePi0 = true;  // disable pi0 if stats too low in the peak region or if no spectral peak found

	//===========================================================================//
	// The Meat of the Code: Prefitting the peak, background, and eta peak
	//===========================================================================//

	// Moving tail parameters to FitPi0Peak function
	// Voigt Profile
	//if (iPeakChoice == 7) {
	//	fit->SetParameter("gamma",0.005);	
	//}

	//if (iPeakChoice == 1 || iPeakChoice == 5) {
	//	fit->SetParameter("lambda",0.015);
	//}

	//if (iPeakChoice == 3 || iPeakChoice == 4) { // Crystal Ball
	//	fit->SetParameter("alpha",1);
	//	fit->SetParameter("n",2.);
	//	//fit->FixParameter(3,1.5);
	//}

	//switch (iBkgChoice) {
//		default:
//	}
	// cycle over bkg params
	for (int i = 0; i < nBkg_Par -1 ; i++) { 
		fit->SetParameter(iBkg_Par_0 + i,0.);
	}

	// FIXME testing Aug 16
	if (fitFunct->HasBkg()) {
		fit->SetParameter(iBkg_Par_BScale,1.);
		fit->SetParLimits(iBkg_Par_BScale,Bkg_BScale_Min,Bkg_BScale_Max);
	}

	if (iBkgChoice >= 0 && iBkgChoice <= 4) {
		fit->SetParameter(iBkg_Par_0,meanHeight); //May need to change this later
		fit->SetParLimits(iBkg_Par_0,-1.*maxHeight,4*maxHeight); 
	}

	// Exp * Poly(n) 
	if (iBkgChoice == 5 || iBkgChoice == 6) {
		fit->SetParameter(iBkg_Par_0, meanHeight);  // N
		fit->SetParameter(iBkg_Par_0 + 1, 1.); // Lambda 
			fit->SetParLimits(iBkg_Par_0 + 1, -10., 10.);
		fit->SetParameter(iBkg_Par_0 + 2, -1.); 
		fit->SetParameter(iBkg_Par_0 + 3, 3.);
		fit->SetParameter(iBkg_Par_0 + 4, -10.);

		if (iBkgChoice == 6) {
			fit->SetParameter(iBkg_Par_0 + 5, 1);
		}
	}
	
	// ExpDecay * Poly(n) 
	if (iBkgChoice == 7 || iBkgChoice == 8) {
		fit->SetParameter(iBkg_Par_0,meanHeight); // Overall scale
			fit->SetParLimits(iBkg_Par_0,0.,2*meanHeight);
		fit->SetParameter(iBkg_Par_0 + 1,5.7);  // gamma of the exponential //4.8
			fit->SetParLimits(iBkg_Par_0 + 1, 4., 8.); //15.
			// FIXME 
//		fit->FixParameter(iBkg_Par_0+1,10);

//			fit->SetParLimits(iBkg_Par_0 + 1, 0., 15.);
		// the 0 point of the function
//		fit->SetParameter(iBkg_Par_0 + 2, 0.05); 
		fit->SetParameter(iBkg_Par_0 +2, hist->GetXaxis()->GetBinLowEdge(iFirstUsedBin));
			fit->SetParLimits(iBkg_Par_0 + 2, 0, hist->GetXaxis()->GetBinLowEdge(iFirstUsedBin) + 0.03); //0.14
		// FIXME set zero point to real zero point
	//	fit->FixParameter(iBkg_Par_0 +2, hist->GetXaxis()->GetBinLowEdge(iFirstUsedBin));

		// 
		fit->SetParameter(iBkg_Par_0 + 3, -2.); // the slope of the poly(x,1)
			fit->SetParLimits(iBkg_Par_0 + 3, bkg_Slope_Min,bkg_Slope_Max); 
	//		fit->SetParLimits(iBkg_Par_0 + 3, -125, 125);
		
		if (iBkgChoice == 8) {
			fit->SetParameter(iBkg_Par_0 + 4, 0); // the quadratic part of the poly(x,2)
//				fit->SetParLimits(iBkg_Par_0 + 4, -1., 1.);
		}
	}

	printf("BKG inits set\n");

	TH1D * hClone        = (TH1D *) hist->Clone("hClone");

	switch (iPeakChoice) {
		default:
		//	break;
		case 0:
			fit->FixParameter(iPi0_Yield,0.);
			fit->FixParameter(iPi0_Mean,0.135);
			fit->FixParameter(iPi0_Sigma,0.05);

			if (useEta) {
				fit->FixParameter(iEta_Yield,0.);
				fit->FixParameter(iEta_Mean,0.55);
				fit->FixParameter(iEta_Sigma,0.05);
			}

			TH1D * hClonePeakRmv = (TH1D *) hist->Clone("hClonePeakRmv");
	
	///		TH1D * hClonePeakRmv = hist;

			double litMass   = 0.135;
			double bandRange = 0.01; //0.045; // 0.06
		
			hClone->GetXaxis()->SetRangeUser(0.08,0.24); //0.08, 0.2
			double localMax = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());
		
			double peakCenterGuess = localMax;

			TSpectrum *sp = new TSpectrum();
			float spSearchSigma = 5; //1
			float spThreshold   = 0.4;
			sp->Search(hClone,spSearchSigma,"nobackground new",spThreshold);
			printf("   TSpectrum found %d peaks.\n",sp->GetNPeaks());
		
			// TSpectrum often finds a peak at the right edge of the window

			if (sp->GetNPeaks() == 0) { 
				printf("    Not using a Pi0 Peak since none was found.\n");
				fUsePi0 = false;
			} else {

				/* if first one is most significant
					if (sp->GetNPeaks() > 0) {
						peakCenterGuess = sp->GetPositionX()[0];
					}
				*/

				for (int i = 0; i < sp->GetNPeaks(); i++) {
					double localPeak = sp->GetPositionX()[i];
					printf("    TSpectrum found peak at %f (y = %f).\n",localPeak,sp->GetPositionY()[i]);
					if (TMath::Abs(localPeak - litMass) < TMath::Abs(peakCenterGuess - litMass)) {
						peakCenterGuess = localPeak;
					}
				}
				printf("  peakCenterGuess = %f\n",peakCenterGuess);

			}
			fPi0_Mean_Guess = peakCenterGuess;
			fPi0_Yield_Guess = hist->Integral(hist->GetXaxis()->FindBin(fPi0_Mean_Guess - 2.*fPi0_Sigma_Guess),hist->GetXaxis()->FindBin(fPi0_Mean_Guess + 2.*fPi0_Sigma_Guess),"width");
			
			fPi0_Yield_Guess = 0.3 * fPi0_Yield_Guess;
	
			printf("New Pi0 Peak Guess:\n\tY=%f\n\tmu=%f\n\tsigma=%f\n",fPi0_Yield_Guess,fPi0_Mean_Guess,fPi0_Sigma_Guess);
	
			// bandRange is roughly 2 sigma
			float lowerPi0CutRange = peakCenterGuess - 1.0 * bandRange; 
			float upperPi0CutRange = peakCenterGuess + 1.2* bandRange;
			
			if (iBkgChoice == 7 || iBkgChoice == 8) {
				lowerPi0CutRange = 0; // remove the lowest end, fit the exp decay part at the end
			}

	//		fPi0_CutOut = 0.7;
		//	zeroRegion(hClonePeakRmv, peakCenterGuess - 0.8 * bandRange, peakCenterGuess + 1.2 * bandRange);
			// FIXME Aug 17: trying removing entire left part
			zeroRegion(hClonePeakRmv, lowerPi0CutRange, upperPi0CutRange);
			// FIXME also zero out Eta region??
			if (useEta) {
				float f = 0.8;
				zeroRegion(hClonePeakRmv, fEta_Mean_Guess - f * fEta_Sigma_Guess, fEta_Mean_Guess + f * fEta_Sigma_Guess);
			}	
			printf("Done zeroing out regions.  Attempting BKG-Only Fit.\n");

			PrintParameters(fit,false,"PreBKG-Only Fit Parameters");

//			hClonePeakRmv->Fit(fit,"Q0","");
			hClonePeakRmv->Fit(fit,"N","");

			printf("Finished BKG-Only Fit\n");

			hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());


//		hClone->Add(fit,-1);

//			hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 1.*bandRange, peakCenterGuess + 1.*bandRange);
//			double integral = TMath::Abs(hClone->Integral("width"));

			hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 0.5*bandRange,peakCenterGuess + 0.5*bandRange);
			
			printf("Building Gauss Prefit\n");
			TF1 * fitGaus = new TF1("fitGaus","[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
//			TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
			fitGaus->SetParameter(0,fPi0_Yield_Guess);
			fitGaus->SetParLimits(0,fPi0_Yield_Min,fPi0_Yield_Max);
	
			fitGaus->SetParameter(1,fPi0_Mean_Guess);
			fitGaus->SetParLimits(1,fPi0_Mean_Min,fPi0_Mean_Max);

			fitGaus->SetParameter(2,fPi0_Sigma_Guess);
			fitGaus->SetParLimits(2,fPi0_Sigma_Min,fPi0_Sigma_Max);

			//fitGaus->FixParameter(3,0);
			//fitGaus->FixParameter(4,0);
			printf("Attempting Gaussian Prefit\n");

			double statsCheck = hClone->Integral();
			//double minStatsForPi0 = 10; //4000
			double minStatsForPi0 = 5; //4000

			printf("DEBUG: StatsCheck = %f.   Min = %f\n",statsCheck,minStatsForPi0);

			// FIXME disabling this check
			//if (statsCheck <= minStatsForPi0) {fUsePi0 = false;			
			//	printf("DEBUG: Disabling pi0 due to low statistics\n");
			//}

			if (fUsePi0) {
				hClone->Fit(fitGaus,"0");

				float chi2OverNDF = fitGaus->GetChisquare() / fitGaus->GetNDF();
				float maxChi2OverNDF = 300; //2
				printf("Gaussian prefit had chisquared over NDF: %f\n",chi2OverNDF);
				if (chi2OverNDF < maxChi2OverNDF) {
					printf("Pi0 Peak PreFit successful!  Using new params.\n");
					fPi0_Yield_Guess = fitGaus->GetParameter(0);	
					fPi0_Mean_Guess = fitGaus->GetParameter(1);	
					fPi0_Sigma_Guess = fitGaus->GetParameter(2);	
					printf("\tY  = %f\n\tmu = %f\n\tsig= %f\n",fPi0_Yield_Guess,fPi0_Mean_Guess,fPi0_Sigma_Guess);
				}
			}

//			double guessRMS = hClone->GetRMS();
		 // double guessMean = hClone->GetMean();
//			double guessMean = peakCenterGuess;
//			double guessYield = integral;

/*
			printf("Building Gauss Prefit\n");
			TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
			fitGaus->SetParameter(0,sqrt(integral));
			fitGaus->SetParameter(1,guessMean);
				fitGaus->SetParLimits(1,0.1,0.22);

				float MaxRMSGuess = 0.11;
			printf("\tSetting RMS guess %f.  If greater than %f, using half of %f.\n",guessRMS,MaxRMSGuess,MaxRMSGuess);
				if (guessRMS > MaxRMSGuess) guessRMS = 0.5 * MaxRMSGuess;
			fitGaus->SetParameter(2,guessRMS);
				fitGaus->SetParLimits(2,0.005,MaxRMSGuess);
			fitGaus->SetParameter(3,0);
			fitGaus->SetParameter(4,0);

			printf("Attempting Gauss Prefit with::\n");
			PrintParameters(fitGaus,false,"Gauss Prefit Prefit Params");

			hClone->Fit(fitGaus,"QM","",peakCenterGuess - 1*bandRange,peakCenterGuess+1.5*bandRange);

			float chi2OverNDF = fitGaus->GetChisquare() / fitGaus->GetNDF();
			float maxChi2OverNDF = 2;
//			printf("Gaussian prefit had chisquared over NDF: %f\n",chi2OverNDF);

			if (chi2OverNDF < maxChi2OverNDF && fitGaus->GetParameter(0) < 2.*integral) {
				printf("Successful Gaussian Prefit!\n");
				fPi0_Yield_Guess = fitGaus->GetParameter(0);
				fPi0_Yield_Guess *= fPi0_Yield_Guess;
				fPi0_Mean_Guess = fitGaus->GetParameter(1);
				fPi0_Sigma_Guess = TMath::Abs(fitGaus->GetParameter(2));
			} else {
				printf("Not using prefit Gaus");
				fPi0_Yield_Guess = 0.2 * integral;
			}
	*/		
			// Release
			printf("Releasing parameters\n");
			if (fUsePi0) {
				fit->ReleaseParameter(iPi0_Yield);
				fit->ReleaseParameter(iPi0_Mean);
				fit->ReleaseParameter(iPi0_Sigma);
			}

			if (useEta) {
				fit->ReleaseParameter(iEta_Yield);
				fit->ReleaseParameter(iEta_Mean);
				fit->ReleaseParameter(iEta_Sigma);
			}

	}

	// Estimating Eta Parameters
	float Eta_Range_Min = 0.35;  // Defining the eta dominated region. 
	float Eta_Range_Max = 0.75;
	//fPi0_Yield_Guess, fPi0_Mean_Guess, fPi0_Sigma_Guess
	double EtaLowMean = 0;
	double EtaHighMean = 0;

	double EtaMidInt = 0;

	float Eta_Inner_Range = 1.2; // N-Sigma
	float Eta_Outer_Range = 3.4;


	double Bkg_Constant_Guess = 0.;
	double Bkg_Slope_Guess    = 0.;

	if (useEta) {

		int Bin_A = hClone->GetXaxis()->FindBin(fEta_Mean_Guess - Eta_Outer_Range * fEta_Sigma_Guess);
		int Bin_B = hClone->GetXaxis()->FindBin(fEta_Mean_Guess - Eta_Inner_Range * fEta_Sigma_Guess);
		int Bin_C = hClone->GetXaxis()->FindBin(fEta_Mean_Guess + Eta_Inner_Range * fEta_Sigma_Guess);
		int Bin_D = hClone->GetXaxis()->FindBin(fEta_Mean_Guess + Eta_Outer_Range * fEta_Sigma_Guess);

		EtaLowMean = hClone->Integral(Bin_A,Bin_B,"width")
			/ (hClone->GetXaxis()->GetBinUpEdge(Bin_B) - hClone->GetXaxis()->GetBinLowEdge(Bin_A));
		EtaHighMean = hClone->Integral(Bin_C,Bin_D,"width")
			/ (hClone->GetXaxis()->GetBinUpEdge(Bin_D) - hClone->GetXaxis()->GetBinLowEdge(Bin_C));
		
	//	EtaMidInt = hClone->Integral(Bin_B,Bin_C,"width");
	//	fEta_Yield_Guess = EtaMidInt - 0.5 * (EtaLowMean + EtaHighMean) 
	//		* (hClone->GetXaxis()->GetBinUpEdge(Bin_C) - hClone->GetXaxis()->GetBinLowEdge(Bin_B));
	//	fEta_Yield_Max = EtaMidInt;

		printf("Eta Setting:  bins %d %d %d %d\n",Bin_A,Bin_B,Bin_C,Bin_D);


		double EtaAllInt = hClone->Integral(Bin_A,Bin_D,"width");
		
		hClone->GetXaxis()->SetRange(Bin_B,Bin_C);
		double EtaPeakMax = hClone->GetBinContent(hClone->GetMaximumBin()) - 0.5 * (EtaLowMean + EtaHighMean);

		fEta_Mean_Guess = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());
		fEta_Mean_Min = fEta_Mean_Guess - Eta_Inner_Range * fEta_Sigma_Guess ;
		fEta_Mean_Max = fEta_Mean_Guess + Eta_Inner_Range * fEta_Sigma_Guess ;


		printf("  Before PreFit:\n");
		printf("            : MeanLow = %f, MeanHigh = %f\n",EtaLowMean,EtaHighMean);
		printf("            : etaAllInt = %f\n",EtaAllInt);
		printf("            : etaPeakMax= %f\n",EtaPeakMax);

		fEta_Yield_Guess = EtaAllInt - 0.5* (EtaLowMean + EtaHighMean) 
			* (hClone->GetXaxis()->GetBinUpEdge(Bin_D) - hClone->GetXaxis()->GetBinLowEdge(Bin_A));

		fEta_Yield_Guess = EtaPeakMax;
		fEta_Yield_Max = 2.*EtaPeakMax;


		//	if (iBkgChoice == 6 || iBkgChoice == 7) {
			// Attempting the Gaussian prefit of the eta region, including a linear bkg.
		
		double bkg_A_Guess = 0; // The first par
	//	double bkg_A_Guess = fit->GetParameter(iBkg_Par_0); // The first par
		double bkg_B_Guess = 0; // the slope par

		if (iBkgChoice == 7 || iBkgChoice == 8) {
	//		bkg_B_Guess = fit->GetParameter(iBkg_Par_0 + 3); // index of slope par in 6,7
			bkg_B_Guess = 0;
		}

		// Could make this guess fix the eta peak to 0 for non-eta pt bins FIXME

		//TF1 * fEtaGaus = new TF1("fitEtaGaus","[0]*x",Eta_Range_Min,maxX);
		TF1 * fEtaGaus = new TF1("fitEtaGaus","[0]*TMath::Gaus(x,[1],[2],0) + [3] + [4]*x",Eta_Range_Min,maxX);
		//TF1 * fEtaGaus = new TF1("fitEtaGaus","[0]*TMath::Gaus(x,[1],[2],0) + [3] + [4]*x",Eta_Range_Min,maxX,5);

		fEtaGaus->SetParameter(0,fEta_Yield_Guess);    fEtaGaus->SetParLimits(0,fEta_Yield_Min,fEta_Yield_Max);
		fEtaGaus->SetParameter(1,fEta_Mean_Guess);     fEtaGaus->SetParLimits(1,fEta_Mean_Min,fEta_Mean_Max);
		fEtaGaus->SetParameter(2,fEta_Sigma_Guess);    fEtaGaus->SetParLimits(2,fEta_Sigma_Min,fEta_Sigma_Max);

		fEtaGaus->SetParameter(3,bkg_A_Guess);
		fEtaGaus->SetParameter(4,bkg_B_Guess);

		hClone->GetXaxis()->SetRangeUser(Eta_Range_Min,Eta_Range_Max);
		hClone->Fit(fEtaGaus,"0","");

		PrintParameters(fEtaGaus,true,"Eta PreFit:");


		Bkg_Constant_Guess = fEtaGaus->GetParameter(3);
		Bkg_Slope_Guess    = fEtaGaus->GetParameter(4);

//	printf("            : MeanLow = %f, MeanHigh = %f\n",EtaLowMean,EtaHighMean);
//	printf("            : etaAllInt = %f\n",EtaAllInt);
//	printf("            : etaPeakMax= %f\n",EtaPeakMax);

	}
	
//	fEta_Yield_Max = fEta_Yield_Guess * 2.;
	// Returning hClone back to normal (for whatever reason)
	hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

	
	//===========================================================================//
	//===========================================================================//
	printf("Setting Final Limits\n");

	// Setting Final Pi0 Parameter Guesses / Range
	fit->SetParLimits(iPi0_Yield, fPi0_Yield_Min, fPi0_Yield_Max);
	fit->SetParLimits(iPi0_Mean, fPi0_Mean_Min, fPi0_Mean_Max);
	fit->SetParLimits(iPi0_Sigma, fPi0_Sigma_Min, fPi0_Sigma_Max);

	fit->SetParameter(iPi0_Yield, fPi0_Yield_Guess);
	fit->SetParameter(iPi0_Mean,  fPi0_Mean_Guess );
	fit->SetParameter(iPi0_Sigma, fPi0_Sigma_Guess);



	// Setting Final Bkg Parameter Guesses
	//   may need a case-switch structure to do this with precision.
//	fit->SetParLimits(iBkg_Par_BScale,Bkg_BScale_Min,Bkg_BScale_Max);
//	fit->SetParameter(iBkg_Par_BScale,1.0); // Good Initial guess (1 is no additional scaling after norm)


	if (iBkgChoice == 7 || iBkgChoice == 8) {
//		fit->SetParameter(iBkg_Par_0   ,Bkg_Constant_Guess);
//		fit->SetParameter(iBkg_Par_0 + 3 ,Bkg_Slope_Guess); 
//		fit->FixParameter(iBkg_Par_0   ,Bkg_Constant_Guess);
//		fit->FixParameter(iBkg_Par_0 + 3 ,Bkg_Slope_Guess); 
		fit->FixParameter(iBkg_Par_0     ,fit->GetParameter(iBkg_Par_0));
		fit->FixParameter(iBkg_Par_0 + 3 ,fit->GetParameter(iBkg_Par_0 + 3));

	//FIXME temp test
	//	fit->FixParameter(iBkg_Par_0+1,10.);
//		fit->FixParameter(iBkg_Par_0 + 1, 5);
	
	}
	// FIXME	
//	printf("DEBUG, fixing BScale (i = %d) at 1.0\n",iBkg_Par_BScale);
//	fit->FixParameter(iBkg_Par_BScale,1.0);

	// Setting Final Eta Parameter Guesses / Range

	if (!fUsePi0) fit->FixParameter(0,0.); 

	if (useEta) {

		fit->SetParLimits(iEta_Yield, fEta_Yield_Min, fEta_Yield_Max);
		fit->SetParLimits(iEta_Mean, fEta_Mean_Min, fEta_Mean_Max);
		fit->SetParLimits(iEta_Sigma, fEta_Sigma_Min, fEta_Sigma_Max);

		fit->SetParameter(iEta_Yield, fEta_Yield_Guess);
		fit->SetParameter(iEta_Mean,  fEta_Mean_Guess );
		fit->SetParameter(iEta_Sigma, fEta_Sigma_Guess);

//		fit->FixParameter(iEta_Yield, fEta_Yield_Guess); // test
//		fit->FixParameter(iEta_Mean,  fEta_Mean_Guess );  // test
		bool FreezeEta = false; // Freeze ALL eta parameters

		if (FreezeEta) {
			fit->FixParameter(iEta_Yield, fEta_Yield_Guess);
			fit->FixParameter(iEta_Mean,  fEta_Mean_Guess );
			fit->FixParameter(iEta_Sigma, fEta_Sigma_Guess);
		}
	}

	return;
}


// Functor used for derived background function
class MassIntThetaPtFit {
  public:
    MassIntThetaPtFit(TH1D * PairPt, TH1D * Theta, double MinPt, double MaxPt, double MinMass = 0, double MaxMass = 0.5) : hPairPt(PairPt),hTheta(Theta),fMinPt(MinPt),fMaxPt(MaxPt),fMinMass(MinMass),fMaxMass(MaxMass) {
 //     hPairPt->GetXaxis()->SetRangeUser(fMinPt,fMaxPt);
      IntegrateMassBkg(); 
      MassBkgArray.push_back(hMassBkg);
    }
    double operator() (double *x, double *p) {
      return p[0]*TMath::Gaus(x[0],p[1],p[2],1) + p[3] * hMassBkg->Interpolate(x[0]);
//      return p[0]*TMath::Gaus(x[0],p[1],p[2],1);
    }  

  private:
    double fMinPt,fMaxPt,fMinMass,fMaxMass;

    TH1D * hPairPt;
    TH1D * hTheta;
    
    TH1D * hMassBkg; // the integrated mass background spectrum, normalized.

    void IntegrateMassBkg() {

      //double clusterCut = 2; // do something with this
   //   double clusterCut = 0.5; // do something with this
      double clusterCut = 1.; // do something with this

  
      // Maybe change nbins to match input
      hMassBkg = new TH1D(TString::Format("MassBkg_%s",hPairPt->GetName()),"Mass Background",100,fMinMass,fMaxMass);
//      float minPt = hPairPt->GetXaxis()->GetXmin();
//      float maxPt = hPairPt->GetXaxis()->GetXmax();

      hMassBkg->SetTitle(TString::Format("Mass Bkg ( %.2f < p_{T} < %.2f)",fMinPt,fMaxPt));  

//      TH1D * localPairPt = (TH1D * ) hPairPt->Clone(TString::Format("PairPt_%s",hPai));
//      localPairPt->GetXaxis()->SetRangeUser(fMinPt,fMaxPt);

      printf("pT = %.2e\n",hPairPt->GetMean());

      int nBins = hMassBkg->GetNbinsX();
      // Step over Mass bins   
      for (int i = 1; i <= nBins; i++) {
        double mass = hMassBkg->GetXaxis()->GetBinCenter(i);
        // Integrate over Theta, Pt
        double sum = 0; 

        int nThetaSteps = 2000;
        float thetaMin = 0.001;
       // float thetaMin = 0.017;
        
//        thetaMin = 0;
 //       float thetaMax = TMath::Pi();  
          // Can calculate this as function of mass
        float thetaMax = TMath::ACos(1 - 0.5 * mass* mass / (clusterCut * clusterCut));
 
        float thetaRange = thetaMax - thetaMin;
        if (thetaRange <= 0) {
          continue;
        }

        int nPtSteps = 1000;
        float ptRange = fMaxPt - fMinPt;
    
        double pT,theta;
        theta = thetaMin;

        for (int j = 0; j < nThetaSteps; j++) {
        // Integrate over pT
          pT = fMinPt; 
          for (int k = 0; k < nPtSteps; k++) {

            double E_1 = GetE_1(mass,pT,theta);
            double E_2 = GetE_2(mass,pT,theta);
          
            if (E_1 < clusterCut || E_2 < clusterCut) {
              pT += ptRange / nPtSteps;
              continue;
            }
            
            double f1 = hGammaE->GetBinContent(hGammaE->FindBin(E_1));
            double f2 = hGammaE->GetBinContent(hGammaE->FindBin(E_2));


           // double f3 = hPairTheta->Interpolate(theta);
            double f3 = hPairTheta->GetBinContent(hPairTheta->FindBin(theta));


            if (f1*f2*f3 == 0) { 
              pT += ptRange / nPtSteps;
              continue;  
            }        
 //           double value = Jacobian(mass,pT,theta) * f1 * f2;
            double value = Jacobian(mass,pT,theta) * f1 * f2 * f3;
            sum += value * (thetaRange / nThetaSteps) * (ptRange / nPtSteps) ;
 //           sum += value * thetaRange / nThetaSteps;
            

            pT += ptRange / nPtSteps;
          }
   //       printf("   theta = %.3e, i =%d, sum = %.3e\n",theta,i,sum);
          
          theta += thetaRange / nThetaSteps;
        }
  //      printf("mass = %.3e, i =%d, sum = %.3e\n",mass,i,sum);
        hMassBkg->SetBinContent(i,sum);
      }

      double integral = hMassBkg->Integral("width");
      if (integral) hMassBkg->Scale(1./integral);
    }

    double Jacobian(double m, double pT, double theta) {
      Double_t thetaMinCut = 0.001;
      if (theta < thetaMinCut) return 0;  // Need to avoid collinear singularities
      double massCutOff = 0.98 * pT * TMath::Tan(theta/2.);
      if ( m >= massCutOff) return 0; // Avoid 0 or complex denominators
      return (0.5 / TMath::Cos(theta/2)) * m / TMath::Sqrt(TMath::Power(pT * TMath::Sin(theta/2),2) - TMath::Power(m*TMath::Cos(theta/2),2));
    }

    double GetE_1(double m, double pT, double theta) {
      double massCutOff = 0.99 * pT * TMath::Tan(theta/2.);
 //     double x = TMath::Power(pT*TMath::Sin(theta/2),2) - TMath::Power(m*TMath::Cos(theta/2),2);  
      if ( m >= massCutOff) return 0; // Avoid 0 or complex denominators
      return (0.5 / TMath::Cos(theta/2)) * (pT - (1./TMath::Sin(theta/2)*TMath::Sqrt(TMath::Power(pT*TMath::Sin(theta/2),2) - TMath::Power(m*TMath::Cos(theta/2),2))));
    }
    double GetE_2(double m, double pT, double theta) {
      double massCutOff = 0.99 * pT * TMath::Tan(theta/2.);
      if ( m >= massCutOff) return 0; // Avoid 0 or complex denominators
      return (0.5 / TMath::Cos(theta/2)) * (pT + (1./TMath::Sin(theta/2)*TMath::Sqrt(TMath::Power(pT*TMath::Sin(theta/2),2) - TMath::Power(m*TMath::Cos(theta/2),2))));
    }

};



TF1 * fitPi0Peak_MassIntBkg(TH1D * hist, std::string name, TH1D * ptHist, TH1D * thetaHist, double minPt, double maxPt) {
	TH1D * h = (TH1D * ) hist->Clone("h");

  // Finding first non-empty bin.
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.4;
  minX = 0.08;
//  minX = 0.0;
  maxX = 0.30;
   // MassIntThetaPtFit(TH1D * PairPt, TH1D * Theta, double MinPt, double MaxPt, double MinMass = 0, double MaxMass = 0.5) : hPairPt(PairPt),hTheta(Theta),fMinPt(MinPt),fMaxPt(MaxPt),fMinMass(MinMass),fMaxMass(MaxMass) {
  MassIntThetaPtFit * intFit = new MassIntThetaPtFit(ptHist,thetaHist,minPt,maxPt,minX,maxX);

  
  TF1 * fit = new TF1(name.c_str(),intFit,minX,maxX,4);



  double litMass = 0.135;


  double totalIntegral = hist->Integral("width");
	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();


  fit->SetParameter(0,nearPeakIntegral);
  fit->SetParameter(1,meanMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(3,0.5);

  fit->SetParLimits(0,0,2*totalIntegral);
  fit->SetParLimits(1,0.1,0.2);
  fit->SetParLimits(2,0,0.1);


  // Estimating scale for A (background parameter)
  // bkg is normed to one.



  fit->FixParameter(0,0);
  fit->FixParameter(1,litMass);
  fit->FixParameter(2,0.05);

  fit->SetParameter(3,0.5*totalIntegral);
  fit->SetParLimits(3,0,2*totalIntegral);

  hist->Fit(fit,"RQ0");

  return fit;
}


Double_t GammaEFit(Double_t * x, Double_t * par) {
  Double_t y = 0;
  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A

  //      par[4] = p_T
  //      par[5] = theta

  Double_t m = x[0];
//  Double_t gammaThr = 10; // GeV

  Double_t Y = par[0];
  Double_t mu = par[1];
  Double_t sigma = par[2];

  Double_t A = par[3];

  
  Double_t pT = par[4];
  Double_t theta = par[5];

  Double_t E_1 = 0;
  Double_t E_2 = 0;

 // Double_t thetaMinCut = 0.3;
 // Double_t thetaMinCut = 0.01;
  Double_t thetaMinCut = 0.02;
  if (theta < thetaMinCut) return 0;  // Need to avoid collinear singularities

  
//  if ( m < pT) {
 
  double massCutOff = 0.95 * pT * TMath::Tan(theta/2.);
//  if ( m < massCutOff) {
  if ( m >= massCutOff) {
  
 //   E_1 = 0.5 * (1./TMath::Cos(theta/2.)) * (pT - TMath::Sqrt(TMath::Abs(pT*pT - m*m)));
 //   E_2 = 0.5 * (1./TMath::Cos(theta/2.)) * (pT + TMath::Sqrt(TMath::Abs(pT*pT - m*m)));
    E_1 = 0.5 * (1./TMath::Sin(theta/2.)) * (pT - TMath::Sqrt(TMath::Abs(pT*pT - TMath::Power(m/TMath::Tan(theta/2.),2))));
    E_2 = 0.5 * (1./TMath::Sin(theta/2.)) * (pT + TMath::Sqrt(TMath::Abs(pT*pT - TMath::Power(m/TMath::Tan(theta/2.),2))));
 //   E_1 = 0.5 * (1./TMath::Sin(theta/2.)) * (pT - TMath::Sqrt(TMath::Abs(pT*pT - TMath::Power(m/TMath::Tan(theta/2.),2))));
 //   E_2 = 0.5 * (1./TMath::Sin(theta/2.)) * (pT + TMath::Sqrt(TMath::Abs(pT*pT - TMath::Power(m/TMath::Tan(theta/2.),2))));
  }

  // Distribution of gamma Energy
  Double_t f1 = 0;
  Double_t f2 = 0;

  if (!hGammaE) {
    return 0 ; //Missing the GammaE histogram
  }
//  Double_t ClusterCut = 2; //
  Double_t ClusterCut = 1; // 
 // Double_t ClusterCut = 0.5; // 
  if (E_1 > ClusterCut) f1 = hGammaE->Interpolate(E_1);
  if (E_2 > ClusterCut) f2 = hGammaE->Interpolate(E_2);
//  if (E_1 > ClusterCut) f1 = hGammaE->GetBinContent(hGammaE->FindBin(E_1));
//  if (E_2 > ClusterCut) f2 = hGammaE->GetBinContent(hGammaE->FindBin(E_2));

//  if (E_1 + E_2 > 0.001)  {
//    printf("Debug m,pt = (%.2e,%.2e) E_1,E_2 = (%.2e  %.2e) f1,f2=(%.2e,%.2e)\n",m,pT,E_1,E_2,f1,f2);
//  //  printf("E_1 + E_2 = %.2e\n",E_1+E_2);
//  }

 // y = Y*TMath::Gaus(m,mu,sigma,1) + A*TMath::Power(TMath::Cos(theta/2.),-2) * (m / TMath::Sqrt(TMath::Abs(pT*pT - m*m))) * f1 * f2;
  y = Y*TMath::Gaus(m,mu,sigma,1) + 0.5*A*TMath::Power(TMath::Sin(theta/2.),-2) * (m / TMath::Sqrt(TMath::Abs(pT*pT - TMath::Power(m/TMath::Tan(theta/2.),2)))) * f1 * f2;
  return y;

}

Double_t GammaEFitThetaIntegral(Double_t * x, Double_t * par) {
  Double_t y = 0;

  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A

  //      par[4] = p_T
  //      par[5] = theta range
  Double_t thetaR = par[5];  

  Double_t x1[1] = {x[0]};
  Double_t theta = 0; 

  int nSteps = 20;
 // int nSteps = 6000;
  for (int i = 0; i < nSteps; i++) {
    Double_t par1[6] = {0,0,1,par[3],par[4],theta};  
    Double_t deltaY = GammaEFit(x1,par1);
    y += deltaY * thetaR/nSteps;
    theta += thetaR / nSteps;
//    printf("Integral debug: (theta, delta Y) = %.1f,%.1f\n",theta,deltaY);

  }

 // return par[0]*TMath::Gaus(x[0],par[1],par[2],1) + y/ TMath::Pi();
  return par[0]*TMath::Gaus(x[0],par[1],par[2],1) + y/ thetaR;
}


Double_t analyticFit(Double_t * x, Double_t * par) {
  Double_t y = 0;
  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A
  //      par[4] = B
  //      par[5] = T1
  //      par[6] = T2

  //      par[7] = p_T
  //      par[8] = theta

  Double_t m = x[0];
  Double_t gammaThr = 10; // GeV

  Double_t Y = par[0];
  Double_t mu = par[1];
  Double_t sigma = par[2];

  Double_t A = par[3];
  Double_t B = par[4];
  Double_t T1 = par[5];
  Double_t T2 = par[6];

//  Double_t pT = par[7]; 
//  Double_t theta = par[8];
  
  Double_t pT = 5.;
  Double_t theta = 0;

//  T1 = 1;
 // T2 = 1;


//  return m;
//  return theta;

  Double_t E_1 = 0;
  Double_t E_2 = 0;

  if ( m < pT) {
    E_1 = 0.5 * (1./TMath::Cos(theta/2.)) * (pT - TMath::Sqrt(TMath::Abs(pT*pT - m*m)));
    E_2 = 0.5 * (1./TMath::Cos(theta/2.)) * (pT + TMath::Sqrt(TMath::Abs(pT*pT - m*m)));
  }
  // Distribution of gamma Energy
  Double_t f1 = (E_1 > 5.) * ( A * TMath::Exp(-E_1/T1) + B * TMath::Exp(-TMath::Abs(E_1 - gammaThr)/T2));
  Double_t f2 = (E_2 > 5.) * ( A * TMath::Exp(-E_2/T1) + B * TMath::Exp(-TMath::Abs(E_2 - gammaThr)/T2));

 // printf("Debug E_1,E_2,f1,f2 = %.1f\t%.1f\t%.1f\t%.1f\n",E_1,E_2,f1,f2);

  y = Y*TMath::Gaus(m,mu,sigma,1) + TMath::Power(TMath::Cos(theta/2.),-2) * (m / TMath::Sqrt(TMath::Abs(pT*pT - m*m))) * f1 * f2;
  return y;

}



TF1 * fitPi0Peak_1(TH1D * hist, std::string name, bool useEta = false) {
  TH1D * h = (TH1D * ) hist->Clone("h");
  
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.4;
  
  // ========================
  // Unique to technique
  // ========================
  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x + [4]*x + [5]";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"A");
  fit->SetParName(4,"B");
  fit->SetParName(5,"C");
  // ========================



  double litMass = 0.135;
  double integral = hist->Integral();

	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.2);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);
  
  delete h;



  hist->Fit(fit,"0Q","",minX,maxX);
  return fit;
}


TF1 * fitPi0Peak_2(TH1D * hist, std::string name, bool useEta = false) {

	TH1D * h = (TH1D * ) hist->Clone("h");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.45;
 
	printf("D 1\n");
 
  // ========================
  // Unique to technique
  // ========================
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*(x+[5])+[4])*(x+[5])/(x+[5]+[6])";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*x+[4])*(x+[5])/(x+[6])";
//  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x/TMath::Sqrt(TMath::Abs([4]*[4] - x*x))";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1) + ([3]*x + [4]) / (1 + TMath::Exp(-(x+[6])/[5]) )";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"A");
  fit->SetParName(4,"B");

  float maxVal = hist->GetBinContent(hist->GetMaximumBin());


	fit->FixParameter(0,0);
	fit->FixParameter(1,0);
	fit->FixParameter(2,0);
	
  fit->SetParLimits(3,0,1.1*maxVal);
	fit->SetParLimits(4,0.51,20.);

  fit->SetParameter(3,maxVal / 2.);
  fit->SetParameter(4,5.);


	hist->Fit(fit,"R0","",0.22,0.5);


  // ========================


  double litMass = 0.135;

	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.2);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);


//	if (!hist) {
//		fprintf(stderr,"Something is wrong!\n");
//		return linFit;	
//	}

//  hist->Fit(fit,"R");
  hist->Fit(fit,"0Q","",minX,maxX);

	//delete h;
	//delete linFit;
  return fit;
}


TF1 * fitPi0Peak_3(TH1D * hist, std::string name, bool useEta = false) {
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.4;
  
  // ========================
  // Unique to technique
  // ========================
  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"A");
  fit->SetParName(4,"B");
  fit->SetParName(5,"C");
  fit->SetParName(6,"D");


  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");
  hClone->GetXaxis()->SetRangeUser(0.2,0.5);
  hClone->Fit(fit,"R0Q");  // Get ABCD

  // ========================



  double litMass = 0.135;
  double integral = hist->Integral(hist->FindBin(minX),hist->FindBin(0.4),"width");
  printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);
  
  fit->SetParLimits(0,0,integral);
  fit->SetParLimits(1,0.12,0.22);
  fit->SetParLimits(2,0,0.1);

  fit->SetParameter(0,0.2*integral);
  fit->SetParameter(1,litMass);
  fit->SetParameter(2,0.025);


  hist->Fit(fit,"0Q","",minX,maxX);
  return fit;
}


// Gaussian + exp*polynomial
TF1 * fitPi0Peak_4(TH1D * hist, std::string name, bool useEta = false) {
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.5;
  
  // ========================
  // Unique to technique
  // ========================
    TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+[7]*x*x)";
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x)";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");


  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"N");
  fit->SetParName(4,"A");
//  fit->SetParName(5,"B");
//  fit->SetParName(6,"C");
//  fit->SetParName(7,"D");

  fit->SetParameter(3,1);
  fit->SetParameter(4,1);
  fit->SetParameter(5,1);
  fit->SetParameter(6,1);
  fit->SetParameter(7,1);


  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");
  hClone->GetXaxis()->SetRangeUser(0.175,0.5);
  hClone->Fit(fit,"0Q");  // Get ABCD Guesses
  
  double litMass = 0.135;

  // Use sidebands to estimate yield
  double bandRange = 0.1; 
  double sbLow  = hClone->Integral(hClone->FindBin(litMass - 1.5*bandRange),hClone->FindBin(litMass - 0.5*bandRange),"width");
  double sbHigh = hClone->Integral(hClone->FindBin(litMass + 0.5*bandRange),hClone->FindBin(litMass + 1.5*bandRange),"width");


  double peakIntegral = TMath::Abs(hClone->Integral(hClone->FindBin(litMass-0.5*bandRange),hClone->FindBin(litMass+0.5*bandRange),"width"));

  // ========================
  double integral = peakIntegral - 0.5 * (sbLow + sbHigh);
  printf("Debug: SBLow = %.3e, SBHigh = %.3e, peakRegion = %.3e, guessIntegral = %.3e\n",sbLow,sbHigh,peakIntegral,integral);
  // Getting RMS
  TF1 * line = new TF1("temp","[0]+[1]*x",0,0.5);
  double slope = (sbHigh - sbLow) / (2 * bandRange * bandRange); // result of rise/run
  double yInt = sbLow/bandRange - slope * (litMass - bandRange);
  line->SetParameter(0,yInt);
  line->SetParameter(1,slope);
  hClone->Add(line,-1);
  hClone->GetXaxis()->SetRangeUser(litMass - 0.5*bandRange, litMass + 0.5*bandRange);
  double guessRMS = hClone->GetRMS();
  printf("Debug: guessRMS = %.3e, secondGuessIntegral = %.3e\n",guessRMS,hClone->Integral("width"));
  double guessMean = hClone->GetMean();

  //printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);
  
  fit->SetParLimits(0,0,integral);
  fit->SetParLimits(1,0.12,0.22);
  fit->SetParLimits(2,0.005,0.1);

  fit->SetParameter(0,0.2*integral);
  fit->SetParameter(1,guessMean);
  fit->SetParameter(2,guessRMS);

  hist->Fit(fit,"0Q","",minX,maxX);
 // hist->Fit(fit,"R0Q","");
////  hist->Fit(fit,"0Q","",minX,0.5);
  return fit;
}


// Gaussian + exp*polynomial, and using a histogram clone with peak zeroed out for bkg fit;
// Third Degree Polynomial
TF1 * fitPi0Peak_5(TH1D * hist, std::string name, bool useEta = false) {
  printf("\tStarting Fit Method 5\n");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
    
  if (n >= hist->GetNbinsX()){
    printf("Empty Hist!\n");
  }

  //double minX = 0.018+hist->GetXaxis()->GetBinLowEdge(n);
  double minX = 0.01+hist->GetXaxis()->GetBinLowEdge(n);
//  double maxX = 0.47;
  double maxX = 0.75;
  
	double etaRangeMin = 0.5;
	double etaRangeMax = 0.7;
	double maxEtaInt = 0;// changed later

	TString formuoli = "";
  // ========================
  // Unique to technique
  // ========================
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+[7]*x*x+[8]*x*x*x)";
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x+(1./6.)*[8]*x*x*x)";
	if (!useEta) {
		formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x+(1./6.)*[8]*x*x*x)";
	} else {
		formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x+(1./6.)*[8]*x*x*x) + [9]*TMath::Gaus(x,[10],[11],1)";
	}
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");


  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"N");
  fit->SetParName(4,"A");
  fit->SetParName(5,"B");
  fit->SetParName(6,"C");
  fit->SetParName(7,"D");
  fit->SetParName(8,"E");

	if (useEta) {
		fit->SetParName(9,"Y_eta");
		fit->SetParName(10,"mu_eta");
		fit->SetParName(11,"sigma_eta");
	}

	double meanHeight = hist->Integral(hist->FindBin(minX),hist->FindBin(maxX),"width")/(maxX - minX);
	double maxHeight  = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParameter(3,meanHeight);
		fit->SetParLimits(3,0,4*maxHeight);
	// First Guesses
	
	// Guesses based on obversvation
  fit->SetParameter(4,1);
	fit->SetParLimits(4,0,4);

  fit->SetParameter(5,-1);
  fit->SetParameter(6,3);
  fit->SetParameter(7,-10);
  fit->SetParameter(8,300);

	if (useEta) {
		maxEtaInt = hist->Integral(hist->FindBin(etaRangeMin),hist->FindBin(etaRangeMax),"width");

		fit->SetParameter(9,0.3*maxEtaInt);
		fit->SetParameter(10,.55);
		fit->SetParameter(11,.05);

		fit->SetParLimits(9,0.,maxEtaInt);
		fit->SetParLimits(10,.5,0.7);
		fit->SetParLimits(11,.01,0.07); //0.01,0.1
	}



  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");

  TH1D * hClonePeakRmv = (TH1D *) hist->Clone("hClonePeakRmv");

  printf("e1\n");

  double litMass = 0.135;
  //double bandRange = 0.08; 
  double bandRange = 0.06; 
  //Pick out the best guess for the mass peak, the local maxima
  hClone->GetXaxis()->SetRangeUser(0.08,0.2);  
  double localMax = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());

  double peakCenterGuess = localMax; // or litMass

  // Use TSpectrum tool to find local max.  Pick closest peak to literature mass
  TSpectrum *sp  = new TSpectrum();
  sp->Search(hClone,1,"nobackground new");
  //sp->Search(hClone);
  for (int i = 0; i < sp->GetNPeaks(); i++) {
    double localPeak = sp->GetPositionX()[i];
    printf("    TSpectrum found peak at %f.\n",localPeak);
    //looking for peak closest to litmass;
    if (TMath::Abs(localPeak - litMass) < TMath::Abs(peakCenterGuess - litMass)) {
      peakCenterGuess = localPeak;
    }
  }
  printf("  peakCenterGuess = %f\n",peakCenterGuess);

  zeroRegion(hClonePeakRmv,peakCenterGuess - 0.8*bandRange, peakCenterGuess + 0.8*bandRange);
  printf("e2\n");
	hClonePeakRmv->Draw();
	fit->Draw("SAME");
 // hClonePeakRmv->Fit(fit,"Q0","",minX,maxX); // 0 \pm 0 Bins are automatically ignored
  hClonePeakRmv->Fit(fit,"Q0",""); // 0 \pm 0 Bins are automatically ignored

  printf("e3\n");
  
  hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  hClone->Add(fit,-1);
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 1.*bandRange,peakCenterGuess + 1.*bandRange);
  double integral = TMath::Abs(hClone->Integral("width"));
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 0.5*bandRange,peakCenterGuess + 0.5*bandRange);
  double guessRMS = hClone->GetRMS();
 // double guessMean = hClone->GetMean();
  double guessMean = peakCenterGuess;
	double guessYield = integral;

  printf("Got guess Yield = %f.\n",guessYield);

  TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
  fitGaus->SetParameter(0,sqrt(integral));
  fitGaus->SetParameter(1,guessMean);
		fitGaus->SetParLimits(1,0.1,0.22);

		float MaxRMSGuess = 0.11;
  printf("\tSetting RMS guess %f.  If greater than %f, using half of %f.\n",guessRMS,MaxRMSGuess,MaxRMSGuess);
		if (guessRMS > MaxRMSGuess) guessRMS = 0.5 * MaxRMSGuess;
  fitGaus->SetParameter(2,guessRMS);
		fit->SetParLimits(2,0.005,MaxRMSGuess);
  fitGaus->SetParameter(3,0);
  fitGaus->SetParameter(4,0);

  hClone->Fit(fitGaus,"QM","",peakCenterGuess - 1*bandRange,peakCenterGuess+1.5*bandRange);
  
	printf("DEBUG: NDF = %d\n",fitGaus->GetNDF());
  float maxChi2OverNDF = 2;
  float chi2OverNDF = 100 * maxChi2OverNDF; //temp
	if (fitGaus->GetNDF() > 0) chi2OverNDF = fitGaus->GetChisquare() / fitGaus->GetNDF();
  printf("Gaussian prefit had chisquared over NDF: %f\n",chi2OverNDF);
  if (chi2OverNDF < maxChi2OverNDF && fitGaus->GetParameter(0) < 2.*integral) {
    guessYield = fitGaus->GetParameter(0);
    guessYield *= guessYield;
    guessMean = fitGaus->GetParameter(1);
    guessRMS = TMath::Abs(fitGaus->GetParameter(2));
  } else {
    printf("Not using prefit Gaus");
    guessYield = 0.2*integral;
  }
  printf("e4\n");


	// Eta estimates
	if (useEta) {
		fit->SetParameter(9,0.3*maxEtaInt);
		fit->SetParameter(10,.55);
		fit->SetParameter(11,.05);

		fit->SetParLimits(9,0.,maxEtaInt);
		fit->SetParLimits(10,.5,0.7);
		fit->SetParLimits(11,.01,0.1);
	}


  //printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);

  printf("\tSetting yield guess %f.  If negative, taking abs value\n",guessYield);
  guessYield = TMath::Abs(guessYield);
  fit->SetParameter(0,guessYield);
  printf("\tSetting mean guess  %f\n",guessMean);
  fit->SetParameter(1,guessMean);
  printf("\tSetting rms guess   %f\n",guessRMS);
  fit->SetParameter(2,guessRMS);

  fit->ReleaseParameter(0);
  fit->ReleaseParameter(1);
  fit->ReleaseParameter(2);



  printf("\tSetting yield limits\n");
  fit->SetParLimits(0,0,TMath::Abs(5*guessYield));
  printf("\tSetting mean limits\n");
  fit->SetParLimits(1,0.1,0.22);
  printf("\tSetting sigma limits\n");
  fit->SetParLimits(2,0.01,0.1);


  printf("\tFitting\n");
  hist->Fit(fit,"0","",minX,maxX);
  printf("\tDone Fitting\n");
  return fit;
}


// Gaussian + exp*polynomial, and using a histogram clone with peak zeroed out for bkg fit;
// 2nd Degree Polynomial
TF1 * fitPi0Peak_6(TH1D * hist, std::string name, bool useEta = false) {
  printf("\tStarting Fit Method 6\n");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
    
  if (n >= hist->GetNbinsX()){
    printf("Empty Hist!\n");
  }

  double minX = 0.018+hist->GetXaxis()->GetBinLowEdge(n);
//  double maxX = 0.7; // Normal
  double maxX = 0.75;

  double etaRangeMin = 0.5;
  double etaRangeMax = 0.7;
  double maxEtaInt = 0;// changed later

  // ========================
  // Unique to technique
  // ========================
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+[7]*x*x+[8]*x*x*x)";
    TString formuoli = "";
  // ========================
/*	if (!useEta) {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x)";
	} else {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x) + [8]*TMath::Gaus(x,[9],[10],1)";
	}*/

	if (!useEta) {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Exp(-[4]*x)*TMath::Power([5]+[6]*x+(1./2.)*[7]*x*x,2)";
	} else {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Exp(-[4]*x)*TMath::Power([5]+[6]*x+(1./2.)*[7]*x*x,2) + [8]*TMath::Gaus(x,[9],[10],1)";
	}

/*	if (!useEta) {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1) + TMath::Exp(-[3]*x)*TMath::Power([4]+[5]*x+(1./2.)*[6]*x*x + (1./6.)*[7]*x*x*x,2)";
	} else {
		formuoli = "[0]*TMath::Gaus(x,[1],[2],1) + TMath::Exp(-[3]*x)*TMath::Power([4]+[5]*x+(1./2.)*[6]*x*x + (1./6.)*[7] * x * x * x,2) + [8]*TMath::Gaus(x,[9],[10],1)";
	}
*/

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"N");
  fit->SetParName(4,"A");
  fit->SetParName(5,"B");
  fit->SetParName(6,"C");
  fit->SetParName(7,"D");

  if (useEta) {
    fit->SetParName(8,"Y_eta");
    fit->SetParName(9,"mu_eta");
    fit->SetParName(10,"sigma_eta");
  }

	double meanHeight = hist->Integral(hist->FindBin(minX),hist->FindBin(maxX),"width")/(maxX - minX);
	double maxHeight  = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParameter(3,meanHeight);
		fit->SetParLimits(3,0,4*maxHeight);
	// First Guesses
	
	// Guesses based on obversvation
  fit->SetParameter(4,1);
	fit->SetParLimits(4,0,4);

  fit->SetParameter(5,-1);
  fit->SetParameter(6,3);
  fit->SetParameter(7,-10);

	if (useEta) {
    maxEtaInt = hist->Integral(hist->FindBin(etaRangeMin),hist->FindBin(etaRangeMax),"width");

    fit->SetParameter(8,0.3*maxEtaInt);
    fit->SetParameter(9,.55);
    fit->SetParameter(10,.05);

    fit->SetParLimits(8,0.,maxEtaInt);
    fit->SetParLimits(9,.5,0.7);
    fit->SetParLimits(10,.01,0.1);
  }

  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");

  TH1D * hClonePeakRmv = (TH1D *) hist->Clone("hClonePeakRmv");

  printf("e1\n");

  double litMass = 0.135;
  //double bandRange = 0.08; 
  double bandRange = 0.06; 
  //Pick out the best guess for the mass peak, the local maxima
  hClone->GetXaxis()->SetRangeUser(0.08,0.2);  
  double localMax = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());

  double peakCenterGuess = localMax; // or litMass

  // Use TSpectrum tool to find local max.  Pick closest peak to literature mass
  TSpectrum *sp  = new TSpectrum();
  sp->Search(hClone);
  for (int i = 0; i < sp->GetNPeaks(); i++) {
    double localPeak = sp->GetPositionX()[i];
    printf("    TSpectrum found peak at %f.\n",localPeak);
    //looking for peak closest to litmass;
    if (TMath::Abs(localPeak - litMass) < TMath::Abs(peakCenterGuess - litMass)) {
      peakCenterGuess = localPeak;
    }
  }
  printf("  peakCenterGuess = %f\n",peakCenterGuess);

  zeroRegion(hClonePeakRmv,peakCenterGuess - 0.8*bandRange, peakCenterGuess + 0.8*bandRange);
  printf("e2\n");
	hClonePeakRmv->Draw();
	fit->Draw("SAME");
 // hClonePeakRmv->Fit(fit,"Q0","",minX,maxX); // 0 \pm 0 Bins are automatically ignored
  hClonePeakRmv->Fit(fit,"Q0",""); // 0 \pm 0 Bins are automatically ignored

  printf("e3\n");
  
  hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  hClone->Add(fit,-1);
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 1.*bandRange,peakCenterGuess + 1.*bandRange);
  double integral = TMath::Abs(hClone->Integral("width"));
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 0.5*bandRange,peakCenterGuess + 0.5*bandRange);
  double guessRMS = hClone->GetRMS();
  double guessMean = hClone->GetMean();
	double guessYield = integral;


  TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
  fitGaus->SetParameter(0,sqrt(integral));
  fitGaus->SetParameter(1,guessMean);
		fitGaus->SetParLimits(1,0.1,0.22);

		float MaxRMSGuess = 0.11;
		if (guessRMS > MaxRMSGuess) guessRMS = 0.5 * MaxRMSGuess;
  fitGaus->SetParameter(2,guessRMS);
		fit->SetParLimits(2,0.005,MaxRMSGuess);
  fitGaus->SetParameter(3,0);
  fitGaus->SetParameter(4,0);

  hClone->Fit(fitGaus,"QM","",peakCenterGuess - 1*bandRange,peakCenterGuess+1.5*bandRange);
  guessYield = fitGaus->GetParameter(0);
	guessYield *= guessYield;
  guessMean = fitGaus->GetParameter(1);
  guessRMS = TMath::Abs(fitGaus->GetParameter(2));

  printf("e4\n");

  // Eta estimates
  if (useEta) {
    fit->SetParameter(8,0.3*maxEtaInt);
    fit->SetParameter(9,.55);
    fit->SetParameter(10,.05);

    fit->SetParLimits(8,0.,maxEtaInt);
    fit->SetParLimits(9,.5,0.7);
    fit->SetParLimits(10,.01,0.1);
  }



  //printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);

  printf("\tSetting yield guess %f.  If negative, taking abs value\n",guessYield);
  guessYield = TMath::Abs(guessYield);
  fit->SetParameter(0,guessYield);
  printf("\tSetting mean guess  %f\n",guessMean);
  fit->SetParameter(1,guessMean);
  printf("\tSetting rms guess   %f\n",guessRMS);
  fit->SetParameter(2,guessRMS);

  fit->ReleaseParameter(0);
  fit->ReleaseParameter(1);
  fit->ReleaseParameter(2);



  printf("\tSetting yield limits\n");
  fit->SetParLimits(0,0,TMath::Abs(5*guessYield));
  printf("\tSetting mean limits\n");
  fit->SetParLimits(1,0.1,0.22);
  printf("\tSetting sigma limits\n");
  fit->SetParLimits(2,0.01,0.1);


  printf("\tFitting\n");
  hist->Fit(fit,"0","",minX,maxX);
  printf("\tDone Fitting\n");
  return fit;
}

// ExponentialDecay-Gaussian + exp*polynomial, and using a histogram clone with peak zeroed out for bkg fit;
// 2nd Degree Polynomial
TF1 * fitPi0Peak_7(TH1D * hist, std::string name, bool useEta = false) {
  printf("\tStarting Fit Method 7\n");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
    
  if (n >= hist->GetNbinsX()){
    printf("Empty Hist!\n");
  }

 
  double minX = 0.01+hist->GetXaxis()->GetBinLowEdge(n); 
  double maxX = 0.45;
  // ========================
  // Unique to technique
  // ========================
    TString formuoli ="[0]*(TMath::Gaus(x,[1],[2],1)+TMath::Exp((x-[1])/[3])*(TMath::Gaus(0,0,[2],1) - TMath::Gaus(x,[1],[2],1))*([1] > x))+[4]*TMath::Exp(-[5]*x)*([6]+[7]*x+(1./2.)*[8]*x*x + [9]*(1./6.)*x*x*x)";
 //   TString formuoli ="[0]*(TMath::Gaus(x,[1],[2],1)+TMath::Exp((x-[1])/[3])*(TMath::Gaus(0,0,[2],1) - TMath::Gaus(x,[1],[2],1))*([1] > x))+[4]*TMath::Exp(-[5]*x)*([6]+[7]*x+(1./2.)*[8]*x*x)";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");
  fit->SetParName(3,"lambda");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(4,"N");
  fit->SetParName(5,"A");
  fit->SetParName(7,"B");
  fit->SetParName(8,"C");
  fit->SetParName(9,"D");

	double meanHeight = hist->Integral(hist->FindBin(minX),hist->FindBin(maxX),"width")/(maxX - minX);
	double maxHeight  = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParameter(4,meanHeight);
		fit->SetParLimits(4,0,4*maxHeight);
	// First Guesses
	
	// Guesses based on obversvation
  fit->SetParameter(3,0.03);

  fit->SetParameter(5,1);
	fit->SetParLimits(5,0,4);

  fit->SetParameter(6,-1);
  fit->SetParameter(7,3);
  fit->SetParameter(8,-10);


  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);
  fit->FixParameter(3,0.03);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");

  TH1D * hClonePeakRmv = (TH1D *) hist->Clone("hClonePeakRmv");

  printf("e1\n");

  double litMass = 0.135;
  //double bandRange = 0.08; 
  double bandRange = 0.06; 
  //Pick out the best guess for the mass peak, the local maxima
  hClone->GetXaxis()->SetRangeUser(0.08,0.2);  
  double localMax = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());

  double peakCenterGuess = localMax; // or litMass

  // Use TSpectrum tool to find local max.  Pick closest peak to literature mass
  TSpectrum *sp  = new TSpectrum();
  sp->Search(hClone);
  for (int i = 0; i < sp->GetNPeaks(); i++) {
    double localPeak = sp->GetPositionX()[i];
    printf("    TSpectrum found peak at %f.\n",localPeak);
    //looking for peak closest to litmass;
    if (TMath::Abs(localPeak - litMass) < TMath::Abs(peakCenterGuess - litMass)) {
      peakCenterGuess = localPeak;
    }
  }
  printf("  peakCenterGuess = %f\n",peakCenterGuess);

  zeroRegion(hClonePeakRmv,peakCenterGuess - 0.8*bandRange, peakCenterGuess + 0.8*bandRange);
  printf("e2\n");
	hClonePeakRmv->Draw();
	fit->Draw("SAME");
 // hClonePeakRmv->Fit(fit,"Q0","",minX,maxX); // 0 \pm 0 Bins are automatically ignored
  hClonePeakRmv->Fit(fit,"Q0",""); // 0 \pm 0 Bins are automatically ignored

  printf("e3\n");
  
  hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  hClone->Add(fit,-1);
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 1.*bandRange,peakCenterGuess + 1.*bandRange);
  double integral = TMath::Abs(hClone->Integral("width"));
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 0.5*bandRange,peakCenterGuess + 0.5*bandRange);
  double guessRMS = hClone->GetRMS();
  double guessMean = hClone->GetMean();
	double guessYield = integral;


  TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
  fitGaus->SetParameter(0,sqrt(integral));
  fitGaus->SetParameter(1,guessMean);
		fitGaus->SetParLimits(1,0.1,0.22);

		float MaxRMSGuess = 0.11;
		if (guessRMS > MaxRMSGuess) guessRMS = 0.5 * MaxRMSGuess;
  fitGaus->SetParameter(2,guessRMS);
		fit->SetParLimits(2,0.005,MaxRMSGuess);
//  fitGaus->SetParameter(3,0);
  fitGaus->SetParameter(4,0);

  hClone->Fit(fitGaus,"QM","",peakCenterGuess - 1*bandRange,peakCenterGuess+1.5*bandRange);
  guessYield = fitGaus->GetParameter(0);
	guessYield *= guessYield;
  guessMean = fitGaus->GetParameter(1);
  guessRMS = TMath::Abs(fitGaus->GetParameter(2));

  printf("e4\n");

  //printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);

  printf("\tSetting yield guess %f.  If negative, taking abs value\n",guessYield);
  guessYield = TMath::Abs(guessYield);
  fit->SetParameter(0,guessYield);
  printf("\tSetting mean guess  %f\n",guessMean);
  fit->SetParameter(1,guessMean);
  printf("\tSetting rms guess   %f\n",guessRMS);
  fit->SetParameter(2,guessRMS);

  fit->ReleaseParameter(0);
  fit->ReleaseParameter(1);
  fit->ReleaseParameter(2);
  fit->ReleaseParameter(3);



  printf("\tSetting yield limits\n");
  fit->SetParLimits(0,0,TMath::Abs(5*guessYield));
  printf("\tSetting mean limits\n");
  fit->SetParLimits(1,0.1,0.22);
  printf("\tSetting sigma limits\n");
  fit->SetParLimits(2,0.01,0.1);
  printf("\tSetting lambda limits\n");
  fit->SetParLimits(3,0.01,0.2);


  printf("\tFitting\n");
  hist->Fit(fit,"0","",minX,maxX);
  printf("\tDone Fitting\n");
  return fit;
}


// Breit-Wigner + exp*(3rd degree polynomial)
TF1 * fitPi0Peak_8(TH1D * hist, std::string name, bool useEta = false) {
  printf("\tStarting Fit Method 8\n");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
    
  if (n >= hist->GetNbinsX()){
    printf("Empty Hist!\n");
  }

  //double minX = 0.018+hist->GetXaxis()->GetBinLowEdge(n);
  double minX = 0.009+hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.65;
  
  // ========================
  // Unique to technique
  // ========================
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+[7]*x*x+[8]*x*x*x)";
   // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x+(1./6.)*[8]*x*x*x)";
    TString formuoli ="[0]*TMath::BreitWigner(x,[1],[2]*TMath::Sqrt(2.*TMath::Log(2)))+[3]*TMath::Exp(-[4]*x)*([5]+[6]*x+(1./2.)*[7]*x*x+(1./6.)*[8]*x*x*x)";
  // ========================

  TF1 * fit = new TF1(name.c_str(),formuoli,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma"); 


  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"N");
  fit->SetParName(4,"A");
  fit->SetParName(5,"B");
  fit->SetParName(6,"C");
  fit->SetParName(7,"D");
  fit->SetParName(8,"E");

	double meanHeight = hist->Integral(hist->FindBin(minX),hist->FindBin(maxX),"width")/(maxX - minX);
	double maxHeight  = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParameter(3,meanHeight);
		fit->SetParLimits(3,0,4*maxHeight);
	// First Guesses
	
	// Guesses based on obversvation
  fit->SetParameter(4,1);
	fit->SetParLimits(4,0,4);

  fit->SetParameter(5,-1);
  fit->SetParameter(6,3);
  fit->SetParameter(7,-10);
  fit->SetParameter(8,300);


  fit->FixParameter(0,0);
  fit->FixParameter(1,0);
  fit->FixParameter(2,0.05);

  TH1D * hClone = (TH1D *) hist->Clone("hClone");

  TH1D * hClonePeakRmv = (TH1D *) hist->Clone("hClonePeakRmv");

  printf("e1\n");

  double litMass = 0.135;
  //double bandRange = 0.08; 
  double bandRange = 0.04; 
  //Pick out the best guess for the mass peak, the local maxima
  hClone->GetXaxis()->SetRangeUser(0.08,0.2);  
  double localMax = hClone->GetXaxis()->GetBinCenter(hClone->GetMaximumBin());

  double peakCenterGuess = localMax; // or litMass

  // Use TSpectrum tool to find local max.  Pick closest peak to literature mass
  double fPeakSearchSigma = 2; // default is 2
	double fPeakSearchThreshold = 0.05;
  TSpectrum *sp  = new TSpectrum();
  sp->Search(hClone,fPeakSearchSigma,"",fPeakSearchThreshold);

  for (int i = 0; i < sp->GetNPeaks(); i++) {
    double localPeak = sp->GetPositionX()[i];
    printf("    TSpectrum found peak at %f.\n",localPeak);
    //looking for peak closest to litmass;
    //FIXME look for most significant peak???
    if (TMath::Abs(localPeak - litMass) < TMath::Abs(peakCenterGuess - litMass)) {
      peakCenterGuess = localPeak;
    }
  }
  printf("  peakCenterGuess = %f\n",peakCenterGuess);

  zeroRegion(hClonePeakRmv,peakCenterGuess - 0.8*bandRange, peakCenterGuess + 0.8*bandRange);
  printf("e2\n");
	hClonePeakRmv->Draw();
	fit->Draw("SAME");
 // hClonePeakRmv->Fit(fit,"Q0","",minX,maxX); // 0 \pm 0 Bins are automatically ignored
  hClonePeakRmv->Fit(fit,"Q0",""); // 0 \pm 0 Bins are automatically ignored

  printf("e3\n");
  
  hClone->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  hClone->Add(fit,-1);
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 1.*bandRange,peakCenterGuess + 1.*bandRange);
  double integral = TMath::Abs(hClone->Integral("width"));  // replace with sum of all positive values?
  hClone->GetXaxis()->SetRangeUser(peakCenterGuess - 0.5*bandRange,peakCenterGuess + 0.5*bandRange);
  double guessRMS = hClone->GetRMS();
  double guessMean = hClone->GetMean();
	double guessYield = integral;


  TF1 * fitGaus = new TF1("fitGaus","[0]*[0]*TMath::Gaus(x,[1],[2],1) + [3] + [4]*x",minX,maxX);
  fitGaus->SetParameter(0,sqrt(integral));
  fitGaus->SetParameter(1,guessMean);
		fitGaus->SetParLimits(1,0.1,0.22);
		float MaxRMSGuess = 0.11;
		if (guessRMS > MaxRMSGuess) guessRMS = 0.5 * MaxRMSGuess;
  fitGaus->SetParameter(2,guessRMS);
		fit->SetParLimits(2,0.005,MaxRMSGuess);
  fitGaus->SetParameter(3,0);
  fitGaus->SetParameter(4,0);

  fitGaus->FixParameter(4,0);

  hClone->Fit(fitGaus,"QM","",peakCenterGuess - 1*bandRange,peakCenterGuess+1.5*bandRange);
  guessYield = fitGaus->GetParameter(0);
	guessYield *= guessYield;
  guessMean = fitGaus->GetParameter(1);
  guessRMS = TMath::Abs(fitGaus->GetParameter(2));

  printf("e4\n");

  //printf("Debug integrating %.2f to %.2f = %.2f\n",minX,0.4,integral);

  printf("\tSetting yield guess %f.  If negative, taking abs value\n",guessYield);
  guessYield = TMath::Abs(guessYield);
  fit->SetParameter(0,guessYield);
  printf("\tSetting mean guess  %f\n",guessMean);
  fit->SetParameter(1,guessMean);
  printf("\tSetting rms guess   %f\n",guessRMS);
  fit->SetParameter(2,guessRMS);

  fit->ReleaseParameter(0);
  fit->ReleaseParameter(1);
  fit->ReleaseParameter(2);



  printf("\tSetting yield limits\n");
  fit->SetParLimits(0,0,TMath::Abs(5*guessYield));
  printf("\tSetting mean limits\n");
  fit->SetParLimits(1,0.1,0.22);
  printf("\tSetting sigma limits\n");
  fit->SetParLimits(2,0.01,0.1);


  printf("\tFitting\n");
  hist->Fit(fit,"0","",minX,maxX);
  printf("\tDone Fitting\n");
  return fit;
}


/*
// Gaussian + Poly(2) and using the background histogram with a free parameter for the background scale;
// Third Degree Polynomial
TF1 * fitPi0Peak_9(TH1D * hist, std::string name, bool useEta = false, TH1D * bkgHist) {
  printf("\tStarting Fit Method 9\n");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
    
  if (n >= hist->GetNbinsX()){
    printf("Empty Hist!\n");
  }

  //double minX = 0.018+hist->GetXaxis()->GetBinLowEdge(n);
  double minX = 0.01+hist->GetXaxis()->GetBinLowEdge(n);
//  double maxX = 0.47;
  double maxX = 0.75;
  
	double etaRangeMin = 0.5;
	double etaRangeMax = 0.7;
//	double maxEtaInt = 0;// changed later
	
	ParamBkg_Functor * fitFunc = new ParamBkg_Functor(bkgHist,useEta); 
	TF1 * fit = new TF1(name.c_str(),fitFunc,minX,maxX,9);


	return fit;
}*/



TF1 * fitPi0Peak_analytic(TH1D * hist, std::string name, double pT = 3) {

	TH1D * h = (TH1D * ) hist->Clone("h");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.45;
 
	printf("D 1\n");
 
  // ========================
  // Unique to technique
  // ========================
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*(x+[5])+[4])*(x+[5])/(x+[5]+[6])";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*x+[4])*(x+[5])/(x+[6])";
//  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x/TMath::Sqrt(TMath::Abs([4]*[4] - x*x))";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1) + ([3]*x + [4]) / (1 + TMath::Exp(-(x+[6])/[5]) )";
  // ========================

//  TF1 * fit = new TF1(name.c_str(),formuoli,minX,maxX);

  printf("Trying to create function %s on %.1f,%.1f\n",name.c_str(),minX,maxX);

//  TF1 * fit = new TF1(name.c_str(),analyticFit,minX,maxX,7);
  TF1 * fit = new TF1(name.c_str(),analyticFit,minX,maxX,9);
 // TF1 * fit = new TF1(name.c_str(),analyticFit,minX,maxX,9,1);

  printf("D 2\n");

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
//  fit->SetParName(3,"A");
//  fit->SetParName(4,"B");
//  fit->SetParName(5,"C");
//  fit->SetParName(6,"D");

  float maxVal = hist->GetBinContent(hist->GetMaximumBin());

  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A
  //      par[4] = B
  //      par[5] = T1
  //      par[6] = T2

  //      par[7] = p_T
  //      par[8] = theta



  fit->SetParLimits(3,0.0,maxVal);
  fit->SetParLimits(4,0.0,maxVal);


  fit->SetParameter(3,0.5*maxVal);
  fit->SetParameter(4,0.5*maxVal);



  fit->SetParLimits(5,0.01,5);
  fit->SetParLimits(6,0.01,5);

  fit->SetParameter(5,1);
  fit->SetParameter(6,1);

//	fit->SetParameter(3,0);

//	fit->SetParLimits(4,-0.5,5);
//	fit->SetParameter(4,1);

//	fit->SetParLimits(5,0.1,10);
//	fit->SetParameter(5,1);

	
  // ========================
	//Prefit
  // ========================

//	TF1 * linFit = new TF1("linFit","[0]*x+[1]",0.22,0.5);
//	hist->Fit(linFit,"RQ0");


	fit->FixParameter(0,0);
	fit->FixParameter(1,0);
	fit->FixParameter(2,0);
//	h->GetXaxis()->SetRangeUser(0.22,0.5);
	
//	hist->Fit(fit,"R0","",0.22,0.5);

//	h->Add(fit,-1);
//	printf("D 3\n");

//	fit->FixParameter(3,fit->GetParameter(3));
//	fit->FixParameter(4,fit->GetParameter(4));
	
//	fit->FixParameter(5,0.01);
//	fit->FixParameter(6,0);

  fit->FixParameter(7,pT);
  fit->FixParameter(8,0.0);

  // ========================

  printf("D 3\n");

  double litMass = 0.135;

	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.2);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);



//	if (!hist) {
//		fprintf(stderr,"Something is wrong!\n");
//		return linFit;	
//	}

  hist->Fit(fit,"R");

  printf("D 4\n");

  return fit;
}


TF1 * fitPi0Peak_test(TH1D * hist, std::string name, double pT = 3) {

	TH1D * h = (TH1D * ) hist->Clone("h");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.45;
 
	printf("D 1\n");
 
  // ========================
  // Unique to technique
  // ========================
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*(x+[5])+[4])*(x+[5])/(x+[5]+[6])";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*x+[4])*(x+[5])/(x+[6])";
//  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x/TMath::Sqrt(TMath::Abs([4]*[4] - x*x))";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1) + ([3]*x + [4]) / (1 + TMath::Exp(-(x+[6])/[5]) )";
  // ========================

//  TF1 * fit = new TF1(name.c_str(),formuoli,minX,maxX);

  TF1 * fit = new TF1(name.c_str(),testFunction,minX,maxX,3);
 // TF1 * fit = new TF1(name.c_str(),analyticFit,minX,maxX,9,1);

  printf("D 2\n");

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
//  fit->SetParName(3,"A");
//  fit->SetParName(4,"B");
//  fit->SetParName(5,"C");
//  fit->SetParName(6,"D");

  float maxVal = hist->GetBinContent(hist->GetMaximumBin());


//	fit->SetParameter(3,0);

//	fit->SetParLimits(4,-0.5,5);
//	fit->SetParameter(4,1);

//	fit->SetParLimits(5,0.1,10);
//	fit->SetParameter(5,1);

	
  // ========================
	//Prefit
  // ========================

//	TF1 * linFit = new TF1("linFit","[0]*x+[1]",0.22,0.5);
//	hist->Fit(linFit,"RQ0");


	fit->FixParameter(0,0);
	fit->FixParameter(1,0);
	fit->FixParameter(2,0);
//	h->GetXaxis()->SetRangeUser(0.22,0.5);


//	hist->Fit(fit,"R0","",0.22,0.5);

//	h->Add(fit,-1);
//	printf("D 3\n");

//	fit->FixParameter(3,fit->GetParameter(3));
//	fit->FixParameter(4,fit->GetParameter(4));
	
//	fit->FixParameter(5,0.01);
//	fit->FixParameter(6,0);

//  fit->FixParameter(7,pT);
//  fit->FixParameter(8,0.0);

  // ========================

  printf("D 3\n");

  double litMass = 0.135;

	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.2);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);

//	if (!hist) {
//		fprintf(stderr,"Something is wrong!\n");
//		return linFit;	
//	}

  hist->Fit(fit,"R");

  printf("D 4\n");

  return fit;
}



TF1 * fitPi0Peak_GammaE(TH1D * hist, std::string name, double pT = 3) {

	TH1D * h = (TH1D * ) hist->Clone("h");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.45;
 
 
  // ========================
  // Unique to technique
  // ========================
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*(x+[5])+[4])*(x+[5])/(x+[5]+[6])";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*x+[4])*(x+[5])/(x+[6])";
//  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x/TMath::Sqrt(TMath::Abs([4]*[4] - x*x))";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1) + ([3]*x + [4]) / (1 + TMath::Exp(-(x+[6])/[5]) )";
  // ========================

//  TF1 * fit = new TF1(name.c_str(),formuoli,minX,maxX);

  printf("Trying to create function %s on %.1f,%.1f\n",name.c_str(),minX,maxX);

  TF1 * fit = new TF1(name.c_str(),GammaEFit,minX,maxX,6);
  fit->SetNpx(500);

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"A");
  fit->SetParName(4,"pT");
  fit->SetParName(5,"theta");
//  fit->SetParName(6,"D");

  float maxVal = hist->GetBinContent(hist->GetMaximumBin());

  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A

  //      par[4] = p_T
  //      par[5] = theta



  fit->FixParameter(4,pT);
  fit->FixParameter(5,0.1); // Reasonable Theta

  fit->SetParameter(0,0);
  fit->SetParameter(3,1);

//  fit->SetParLimits(3,0.0,1.1*maxVal);
 // fit->SetParameter(3,0.5*maxVal);
  printf("Val at 0.3 = %.1e\n",fit->Eval(0.3));

  
//  fit->Draw("AP"); 
  double fitIntegral = fit->Integral(minX,maxX);
  printf("DEBUG integral = %.1e\n",fitIntegral);

  double litMass = 0.135;

	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();

  double histIntegralOverFunction = nearPeakIntegral / fitIntegral;

  if (fitIntegral > 0) {
    fit->SetParameter(3,histIntegralOverFunction);
    fit->SetParLimits(3,0,3*histIntegralOverFunction);
  }  


//	TF1 * linFit = new TF1("linFit","[0]*x+[1]",0.22,0.5);
//	hist->Fit(linFit,"RQ0");


	

//	h->Add(fit,-1);
//	printf("D 3\n");

//	fit->FixParameter(3,fit->GetParameter(3));
//	fit->FixParameter(4,fit->GetParameter(4));
	
//	fit->FixParameter(5,0.01);
//	fit->FixParameter(6,0);


 // double sampleValue = fit->:

  // ========================
	//Prefit
  // ========================


	fit->FixParameter(0,0);
	fit->FixParameter(1,0);
	fit->FixParameter(2,0);
//	h->GetXaxis()->SetRangeUser(0.22,0.5);

//	hist->Fit(fit,"R0","",0.22,0.5);
  // ========================



  


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.2);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);


  fit->FixParameter(0,0);
	fit->FixParameter(1,0);
	fit->FixParameter(2,0);
//	if (!hist) {
//		fprintf(stderr,"Something is wrong!\n");
//		return linFit;	
//	}

//  fit->SetParameter(3,nearPeakIntegral/(0.25 - x));

//  hist->Fit(fit,"R");

  delete h;

  printf("D 4\n");

  return fit;
}

TF1 * fitPi0Peak_GammaEIntegral(TH1D * hist, std::string name, double pT = 3) {

	TH1D * h = (TH1D * ) hist->Clone("h");
  int n = 1;
  while (hist->GetBinContent(n) == 0) n++;
  
  double minX = hist->GetXaxis()->GetBinLowEdge(n);
  double maxX = 0.45;
 // double maxX = 20;
 
	printf("D 1\n");
 
  // ========================
  // Unique to technique
  // ========================
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*(x+[5])+[4])*(x+[5])/(x+[5]+[6])";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+([3]*x+[4])*(x+[5])/(x+[6])";
//  TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x*x*x + [4]*x*x + [5]*x+[6]";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1)+[3]*x/TMath::Sqrt(TMath::Abs([4]*[4] - x*x))";
 // TString formuoli ="[0]*TMath::Gaus(x,[1],[2],1) + ([3]*x + [4]) / (1 + TMath::Exp(-(x+[6])/[5]) )";
  // ========================

//  TF1 * fit = new TF1(name.c_str(),formuoli,minX,maxX);

  double fitMinX = 0.0;

  printf("Trying to create function %s on %.1f,%.1f\n",name.c_str(),fitMinX,maxX);

  TF1 * fit = new TF1(name.c_str(),GammaEFitThetaIntegral,fitMinX,maxX,6);

  printf("D 2\n");

  fit->SetParName(0,"Y");
  fit->SetParName(1,"mu");
  fit->SetParName(2,"sigma");

  // ========================
  // Unique to technique
  // ========================
  fit->SetParName(3,"A");
  fit->SetParName(4,"pT");
  fit->SetParName(5,"theta");
//  fit->SetParName(6,"D");

  float maxVal = hist->GetBinContent(hist->GetMaximumBin());

  // Pars par[0] = yield
  //      par[1] = mu
  //      par[2] = sigma

  //      par[3] = A

  //      par[4] = p_T
  //      par[5] = theta

  fit->FixParameter(4,pT);
  fit->FixParameter(5,TMath::Pi());

  fit->FixParameter(0,0);
  fit->FixParameter(1,0.135);
  fit->FixParameter(2,0.03); 
  fit->SetParameter(3,1);

// hist->Fit(fit,"Q");

  // Initial Normalization

  double fitSample = fit->Eval(0.3);
  printf("Val at 0.3 = %.1e\n",fitSample);



 // fit->SetParLimits(3,0.0,1.1*maxVal);
 // fit->SetParameter(3,0.5*maxVal);
//  double fitIntegral = fit->Integral(minX,maxX);
 // printf("DEBUG integral = %.1e\n",fitIntegral);

  double farMin = 0.25;
  double farMax = 0.5;
  h->GetXaxis()->SetRangeUser(farMin,farMax);
  double farIntegral = h->Integral("width");
  double farMean = farIntegral / (farMax - farMin);


  if (fitSample > 0) {
    double ratioHistOverFit = farMean / fitSample;
    fit->SetParameter(3,ratioHistOverFit);
    fit->SetParLimits(3,0,8*ratioHistOverFit);

  }


	h->GetXaxis()->SetRangeUser(minX,0.25);
	double nearPeakIntegral = h->Integral("width");
  double meanMass = h->GetMean();
  double rms = h->GetRMS();





 // if (fitIntegral > 0) {
 // double histIntegralOverFunction = nearPeakIntegral / fitIntegral;
 //   fit->SetParameter(3,histIntegralOverFunction);
 //   fit->SetParLimits(3,0,3*histIntegralOverFunction);
 // }  



	
//	h->GetXaxis()->SetRangeUser(0.22,0.5);
	
//	hist->Fit(fit,"R0","",0.22,0.5);

//	h->Add(fit,-1);
//	printf("D 3\n");

//	fit->FixParameter(3,fit->GetParameter(3));
//	fit->FixParameter(4,fit->GetParameter(4));
	
//	fit->FixParameter(5,0.01);
//	fit->FixParameter(6,0);



  // ========================
	//Prefit
  // ========================

//	TF1 * linFit = new TF1("linFit","[0]*x+[1]",0.22,0.5);
//	hist->Fit(linFit,"RQ0");


//	fit->FixParameter(0,0);
//	fit->FixParameter(1,0);
//	fit->FixParameter(2,0.1);




  // ========================


  double litMass = 0.135;

//	h->GetXaxis()->SetRangeUser(minX,0.25);
//	double nearPeakIntegral = h->Integral("width");
//  double meanMass = h->GetMean();
//  double rms = h->GetRMS();


  fit->SetParLimits(0,0,2*nearPeakIntegral);
  fit->SetParLimits(1,minX,0.3);
  fit->SetParLimits(2,0,0.15);

  fit->SetParameter(0,nearPeakIntegral);
//  fit->SetParameter(0,0);
  fit->SetParameter(1,meanMass);
//  fit->SetParameter(1,litMass);
  fit->SetParameter(2,rms);
//  fit->SetParameter(2,0.05);


  fit->FixParameter(0,0);
  fit->FixParameter(1,meanMass);
  fit->FixParameter(2,rms);
//	if (!hist) {
//		fprintf(stderr,"Something is wrong!\n");
//		return linFit;	
//	}

  hist->Fit(fit,"R");

  printf("Val at 0.3 = %.1e\n",fit->Eval(0.3));
//  fit->SetParameter(3,0);


  return fit;
}




TF1 * fitPi0Mass(TGraphErrors * massGraph, std::string name) {
  
  //TString formuoli ="[0]*x + [1]";
  TString formuoli ="([0]*x*x + [1]*x + [2])/(x+[3])";
 

  TF1 * fit = new TF1(name.c_str(),formuoli,massGraph->GetXaxis()->GetXmin(),massGraph->GetXaxis()->GetXmax());
  double xMin = 5;
  double xMax = 20;
  fit->SetParameter(0,1);
  fit->SetParameter(1,1);
  fit->SetParameter(2,1);
  fit->SetParameter(3,1);
  
  fit->SetParLimits(3,0,40);


  massGraph->Fit(fit,"0","",xMin,xMax);
  
  return fit;
}



TF1 * fitPi0Mass_KinkedLine(TGraphErrors * massGraph, std::string name) {
  TString formuoli ="(x < [0]) * ([2]*x+[1]-[2]*[0]) + (x >= [0]) * ([3]*x+[1]-[3]*[0])";
  TF1 * fit = new TF1(name.c_str(),formuoli,massGraph->GetXaxis()->GetXmin(),massGraph->GetXaxis()->GetXmax());
  double xMin = 5;
  double xMax = 20;

  fit->SetParName(0,"d");  // kink location
  fit->SetParName(1,"e");  // kink height
  fit->SetParName(2,"m1"); // First Slope
  fit->SetParName(3,"m2"); // Second Slope

  double D_Guess = 10;
  double D_Min = 7;
  double D_Max = 12;
  double E_Guess = massGraph->Eval(D_Guess);

  double M1_Guess = 0;
  double M2_Guess = 0.003;

  fit->SetParameter(0,D_Guess);
  fit->SetParameter(1,E_Guess);
  fit->SetParameter(2,0);
  fit->SetParameter(3,1);

  fit->SetParLimits(0,D_Min,D_Max);

  massGraph->Fit(fit,"0","",xMin,xMax);
  return fit;
}


TF1 * fitPi0Sigma(TGraphErrors * sigmaGraph, std::string name) {
  
  //TString formuoli ="[0]*x + [1]";
  //TString formuoli ="[0]+[1]*x+[2]*x*x+[3]*x*x*x";
  TString formuoli ="([0]*x*x + [1]*x + [2])/(x+[3])";

  TF1 * fit = new TF1(name.c_str(),formuoli,sigmaGraph->GetXaxis()->GetXmin(),sigmaGraph->GetXaxis()->GetXmax());
  double xMin = 5;
  double xMax = 28;
  fit->SetParameter(0,1);
  fit->SetParameter(1,1);
  fit->SetParameter(2,1);
  fit->SetParameter(3,1);
  
  fit->SetParLimits(3,0,40);


  sigmaGraph->Fit(fit,"0","",xMin,xMax);
  
  return fit;
}


TF1 * fitPi0Sigma_KinkedLine(TGraphErrors * massGraph, std::string name) {
  TString formuoli ="(x < [0]) * ([2]*x+[1]-[2]*[0]) + (x >= [0]) * ([3]*x+[1]-[3]*[0])";
  TF1 * fit = new TF1(name.c_str(),formuoli,massGraph->GetXaxis()->GetXmin(),massGraph->GetXaxis()->GetXmax());
  double xMin = 3.3;
  double xMax = 20;

  fit->SetParName(0,"d");  // kink location
  fit->SetParName(1,"e");  // kink height
  fit->SetParName(2,"m1"); // First Slope
  fit->SetParName(3,"m2"); // Second Slope

  double D_Guess = 10;
  double D_Min = 7;
  double D_Max = 12;
  double E_Guess = massGraph->Eval(D_Guess);

  double M1_Guess = 0;
  double M2_Guess = 0.00083;

  fit->SetParameter(0,D_Guess);
  fit->SetParameter(1,E_Guess);
  fit->SetParameter(2,0);
  fit->SetParameter(3,1);

  fit->SetParLimits(0,D_Min,D_Max);

  massGraph->Fit(fit,"0","",xMin,xMax);
  return fit;
}


// Type 0 -> Mass
// Type 1 -> Sigma
TF1 * fit_KinkedLine(TGraphErrors * massGraph, std::string name, Int_t type) {
  TString formuoli ="(x < [0]) * ([2]*x+[1]-[2]*[0]) + (x >= [0]) * ([3]*x+[1]-[3]*[0])";
  TF1 * fit = new TF1(name.c_str(),formuoli,massGraph->GetXaxis()->GetXmin(),massGraph->GetXaxis()->GetXmax());
	double xMin = 0;
	double xMax = 20;

	if (type == 0) { //Mass
		xMin = 5;
	} else if (type == 1) { //Sigma
		xMin = 3.3;
	}

  fit->SetParName(0,"d");  // kink location
  fit->SetParName(1,"e");  // kink height
  fit->SetParName(2,"m1"); // First Slope
  fit->SetParName(3,"m2"); // Second Slope

  double D_Guess = 10;
  double D_Min = 7;
  double D_Max = 12;
  double E_Guess = massGraph->Eval(D_Guess);

  double M1_Guess = 0;
 	double M2_Guess = 0.01;
	if (type == 0) { // Mass
		M2_Guess = 0.003;
	} else { // Sigma
 		M2_Guess = 0.00083;
	} 

  fit->SetParameter(0,D_Guess);
  fit->SetParameter(1,E_Guess);
  fit->SetParameter(2,0);
  fit->SetParameter(3,1);

  fit->SetParLimits(0,D_Min,D_Max);

  massGraph->Fit(fit,"0","",xMin,xMax);
  return fit;
}
