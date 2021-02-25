#ifndef TASKEVENTPLANEMATHTOOLS_CXX
#define TASKEVENTPLANEMATHTOOLS_CXX

// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>                                  
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>                                                  
#include <TSpectrum.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TLine.h>            
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>              
#include <TSystem.h>
#include <TList.h>
#include <TLatex.h>                                           
#include <TColor.h>
#include <TGaxis.h>               
#include <TStyle.h>
#include <TEnv.h>
#include <TMath.h>
#include <TMinuit.h>
#include <THnSparse.h>         
#include <TDatime.h>

using namespace std;
#include "TaskEventPlane.h"

class TaskEventPlane;

//class RPF_Functor;
//class RPF_Functor_Single;

// I think this should be true, in theory ??
// It depends on if you are scaling the background function by (2c/pi) B * beta_R
// The factor of 2 simply cancels out
// Since I don't currently apply that scale, I should leave this as false
// I believe Raymond's code does apply the (2c/pi) scale, so maybe I will too
const bool kDoubleMidPlane = true; 

const bool kV3BugFix = true;

//const Int_t kFitLineColor = kViolet-5;

// Variables that change in Evt Plane Bins
const Double_t kPhi_S[3] = {0., 0.25, 0.5}; // To be multiplied by pi
  // what about 3/4?
const Double_t kC_S[3]   = {1./6.,1./12.,1./6.};

// FIXME Calculate these 
// FIXME set these based on Cent through config
//                         {R_1,R_2,R_3,R_4,R_5,R_6}
//const Double_t kEPRes[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
// More reasonable guesses
//const Double_t kEPRes[6] = {0.0,0.8,0.0,0.4,0.0,0.1};
// Cent 10-30%
const Double_t kEPRes[6] = {0.0,0.885,0.605,0.245,0.0,0.1};


// Some very absolute cuts;

// May want to copy these into the python code.
// Setting the funct->V2A_Min(),Max()
/*
double kV2T_AbsCut = 0.2;
double kV2A_AbsCut = 0.2;
double kV3_AbsCut = 0.1;
double kV4T_AbsCut = 0.2;
double kV4A_AbsCut = 0.2;
*/


class RPF_Functor {
  public:
    RPF_Functor() {
      // Nothing to be done
    }

    void SetEPRes(int i, double input) { fEPRes[i] = input; }

    void SetInitV2T(double input) {fInitV2T = input;}
    void SetInitV2A(double input) {fInitV2A = input;}
    void SetInitV3(double input) {fInitV3 = input;}
    void SetInitV4T(double input) {fInitV4T = input;}
    void SetInitV4A(double input) {fInitV4A = input;}
    void SetInitV5(double input) {fInitV5 = input;}
    void SetInitV6T(double input) {fInitV6T = input;}
    void SetInitV6A(double input) {fInitV6A = input;}

    void SetFixedV2T(double input) {fFixedV2T = input;}
    void SetFixedV2A(double input) {fFixedV2A = input;}
    void SetFixedV3(double input) {fFixedV3 = input;}
    void SetFixedV4T(double input) {fFixedV4T = input;}
    void SetFixedV4A(double input) {fFixedV4A = input;}
    void SetFixedV5(double input) {fFixedV5 = input;}
    void SetFixedV6T(double input) {fFixedV6T = input;}
    void SetFixedV6A(double input) {fFixedV6A = input;}

    void SetV2TRange(double min, double max) { fV2T_Min = min; fV2T_Max = max; }
    void SetV2ARange(double min, double max) { fV2A_Min = min; fV2A_Max = max; }
    void SetV3Range(double min, double max) { fV3_Min = min; fV3_Max = max; }
    void SetV4TRange(double min, double max) { fV4T_Min = min; fV4T_Max = max; }
    void SetV4ARange(double min, double max) { fV4A_Min = min; fV4A_Max = max; }
    void SetV5Range(double min, double max) { fV5_Min = min; fV5_Max = max; }
    void SetV6TRange(double min, double max) { fV6T_Min = min; fV6T_Max = max; }
    void SetV6ARange(double min, double max) { fV6A_Min = min; fV6A_Max = max; }

    double GetInitV2T() { return fInitV2T; }
    double GetInitV2A() { return fInitV2A; }
    double GetInitV3()  { return fInitV3;  }
    double GetInitV4T() { return fInitV4T; }
    double GetInitV4A() { return fInitV4A; }
    double GetInitV5()  { return fInitV5;  }
    double GetInitV6T() { return fInitV6T; }
    double GetInitV6A() { return fInitV6A; }

    double GetFixedV2T() { return fFixedV2T; }
    double GetFixedV2A() { return fFixedV2A; }
    double GetFixedV3()  { return fFixedV3;  }
    double GetFixedV4T() { return fFixedV4T; }
    double GetFixedV4A() { return fFixedV4A; }
    double GetFixedV5()  { return fFixedV5;  }
    double GetFixedV6T() { return fFixedV6T; }
    double GetFixedV6A() { return fFixedV6A; }

    double GetV2T_Min()  { return fV2T_Min; }
    double GetV2T_Max()  { return fV2T_Max; }
    double GetV2A_Min()  { return fV2A_Min; }
    double GetV2A_Max()  { return fV2A_Max; }

    double GetV3_Min()   { return fV3_Min; }
    double GetV3_Max()   { return fV3_Max; }

    double GetV4T_Min()  { return fV4T_Min; }
    double GetV4T_Max()  { return fV4T_Max; }
    double GetV4A_Min()  { return fV4A_Min; }
    double GetV4A_Max()  { return fV4A_Max; }

    double GetV5_Min()   { return fV5_Min; }
    double GetV5_Max()   { return fV5_Max; }

    double GetV6T_Min()  { return fV6T_Min; }
    double GetV6T_Max()  { return fV6T_Max; }
    double GetV6A_Min()  { return fV6A_Min; }
    double GetV6A_Max()  { return fV6A_Max; }

    double operator() (double *x, double *p) {
      Double_t fLocalDPhi = x[0]; 

      Int_t iEP = 0;
      if (x[0] > 0.5 * TMath::Pi()) iEP = 1;
      if (x[0] > 1.5 * TMath::Pi()) iEP = 2;
      fLocalDPhi -= iEP*TMath::Pi(); // Post Bin Fix


      // The Parameters
/*      Double_t B    = p[0];
      Double_t VT_2 = p[1];
      Double_t VA_2 = p[2];
      Double_t VT_3VA_3 = p[3];
      Double_t VT_4 = p[4];
      Double_t VA_4 = p[5];
*/
      // Idea: use parameter 0 for the window size
      Double_t B    = p[1];
      Double_t VT_2 = p[2];
      Double_t VA_2 = p[3];
      Double_t VT_3VA_3 = p[4];
      Double_t VT_4 = p[5];
      Double_t VA_4 = p[6];
      Double_t VT_5VA_5 = p[7];
      Double_t VT_6 = p[8];
      Double_t VA_6 = p[9];
      

      Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
      Double_t fC_S   = TMath::Pi() * kC_S[iEP];

      Double_t value = 1;
      value += 2. * funct_RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4,VT_6) * VA_2 * TMath::Cos(2.*fLocalDPhi);
      value += 2. * funct_RPF_VR_3(fPhi_S,fC_S,VT_2,VT_3VA_3,VT_4,VT_6) * TMath::Cos(3.*fLocalDPhi);
      value += 2. * funct_RPF_VR_4(fPhi_S,fC_S,VT_2,VT_4,VT_6) * VA_4 * TMath::Cos(4.*fLocalDPhi);

      value += 2. * VT_5VA_5 * TMath::Cos(5.*fLocalDPhi);
      value += 2. * funct_RPF_VR_6(fPhi_S,fC_S,VT_2,VT_4,VT_6) * VA_6 * TMath::Cos(4.*fLocalDPhi);

      value = 2 * fC_S / TMath::Pi() * B * funct_RPF_BR(fPhi_S,fC_S,VT_2,VT_4) * value;

      // FIXME check that this is right for mid plane
      if (iEP == 1 && kDoubleMidPlane) {
        value = 2.0 * value;
      }
      return value;
    }

    double funct_RPF_Beta_3(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0, Double_t vt_6 = 0) {
      double j = 2; // w.r.t. 2nd order event plane
      double n = 3; // this is denominator for v3
      double R_jk_j = fEPRes[1]; // R2
      double C_j1_0_j = 1; // C_j_0_j = 1

      double value = 1 + 2 * vt_2 / (j * 1 * fC_S) * TMath::Sin(n*fC_S) * R_jk_j * TMath::Cos(j*1*fPhi_S);
      return value;
    }

    // Implementation of these functions with access to the proper EPR
    double funct_RPF_BR(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0, Double_t vt_6 = 0) {

      Double_t value = 1 + 2. * vt_2 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * fEPRes[1];
      value +=             2. * vt_4 * TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * fEPRes[3];
      value +=             2. * vt_6 * TMath::Cos(6.*fPhi_S) * TMath::Sin(4.*fC_S) / (6.*fC_S) * fEPRes[5];
      return value;
    }

    double funct_RPF_VR_2(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0, Double_t vt_6 = 0) {

      Double_t value = vt_2;
      value +=        TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * fEPRes[1];
      value += vt_4 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * fEPRes[1];
      value += (vt_2 + vt_6) * TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * fEPRes[3];
      value += vt_4 * TMath::Cos(6.*fPhi_S) * TMath::Sin(6.*fC_S) / (6.*fC_S) * fEPRes[5];
      //value += vt_6 * TMath::Cos(8.*fPhi_S) * TMath::Sin(8.*fC_S) / (8.*fC_S) * fEPRes[7]; // if R8 found

      value = value / funct_RPF_BR(fPhi_S,fC_S,vt_2,vt_4,vt_6) ;

      return value;
    }

    double funct_RPF_VR_3(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t v3, Double_t vt_4 = 0, Double_t vt_6 = 0) {

      double value = v3  * fEPRes[2]/ funct_RPF_Beta_3(fPhi_S,fC_S,vt_2,vt_4,vt_6);
      return value;
    }

    double funct_RPF_VR_4(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0, Double_t vt_6 = 0) {

      Double_t value = vt_4;
      value +=        TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * fEPRes[3];
      value += vt_2 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * fEPRes[1];
      value += vt_2 * TMath::Cos(6.*fPhi_S) * TMath::Sin(6.*fC_S) / (6.*fC_S) * fEPRes[5];
      //value += vt_4 * TMath::Cos(8.*fPhi_S) * TMath::Sin(8.*fC_S) / (8.*fC_S) * fEPRes[7];
      //value += vt6 * TMath::Cos(10.*fPhi_S) * TMath::Sin(10.*fC_S) / (10.*fC_S) * fEPRes[9];

      value = value / funct_RPF_BR(fPhi_S,fC_S,vt_2,vt_4,vt_6) ;

      return value;
    }

    double funct_RPF_VR_6(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0, Double_t vt_6 = 0.) {
      Double_t value = vt_6;
      value +=        TMath::Cos(6.*fPhi_S) * TMath::Sin(6.*fC_S) / (6.*fC_S) * fEPRes[5];
      value += vt_4 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * fEPRes[1];
      value += vt_2 * TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * fEPRes[3];
      //value += vt_2 * TMath::Cos(8.*fPhi_S) * TMath::Sin(8.*fC_S) / (8.*fC_S) * fEPRes[7];
      //value += vt_4 * TMath::Cos(10.*fPhi_S) * TMath::Sin(10.*fC_S) / (10.*fC_S) * fEPRes[9];
      // FIXME
      return value;
    }


    void DebugPrint() {
      printf("Debugging RPF Functor\n");
//      printf(" Parameters: B v2t v2a v3  v4t v4a\n");
//      printf("             %f %f %f %f %f\n");
      TString labels[3] = {"In-Plane","Mid-Plane","Out-of-Plane"};
      double v2example = 0.05;
      double v4example = 0.025;
      for (int i = 0; i < 3; i++) {
        printf("==========================\n");
        printf("  %s\n",labels[i].Data());
        printf("==========================\n");
        Double_t fPhi_S = TMath::Pi() * kPhi_S[i];
        Double_t fC_S   = TMath::Pi() * kC_S[i];
        printf("Fix Parameters: Phi_s = pi*%f C_s = pi*%f\n",kPhi_S[i],kC_S[i]);
        printf("RPF_BR(phi,C,v_2=%f,v_4=%f)=%f\n",v2example,v4example,funct_RPF_BR(fPhi_S,fC_S,v2example,v4example));
        printf("RPF_VR_2(phi,C,v_2=%f,v_4=%f)=%f\n",v2example,v4example,funct_RPF_VR_2(fPhi_S,fC_S,v2example,v4example));
        printf("RPF_VR_4(phi,C,v_2=%f,v_4=%f)=%f\n",v2example,v4example,funct_RPF_VR_4(fPhi_S,fC_S,v2example,v4example));

      }
      printf("==========================\n");
    }

    static const int kTotalNumberOfRn = 6;
  protected:
//    TF1 * fFunc;
    double fEPRes[kTotalNumberOfRn] = {0.0,0.8,0.0,0.4,0.0,0.1};

    double fInitV2T = -1;
    double fInitV2A = -1;
    double fInitV3  = -1;
    double fInitV4T = -1;
    double fInitV4A = -1;
    double fInitV5  = -1;
    double fInitV6T = -1;
    double fInitV6A = -1;

    double fFixedV2T = -1;
    double fFixedV2A = -1;
    double fFixedV3  = -1;
    double fFixedV4T = -1;
    double fFixedV4A = -1;
    double fFixedV5  = -1;
    double fFixedV6T = -1;
    double fFixedV6A = -1;
    
    double fV2T_Min = -1;
    double fV2T_Max = -1;
    double fV2A_Min = -1;
    double fV2A_Max = -1;
    double fV3_Min = -1;
    double fV3_Max = -1;
    double fV4T_Min = -1;
    double fV4T_Max = -1;
    double fV4A_Min = -1;
    double fV4A_Max = -1;

    double fV5_Min = -1;
    double fV5_Max = -1;
    double fV6T_Min = -1;
    double fV6T_Max = -1;
    double fV6A_Min = -1;
    double fV6A_Max = -1;
};

class RPF_Functor_Single : public RPF_Functor {
  public:
    
    double operator() (double *x, double *p) {
      Double_t fLocalDPhi = x[0];

      Int_t iEP = (Int_t) p[0]; // 0,1, or 2, or -1 for All 
      Double_t B    = p[1];
      Double_t VT_2 = p[2];
      Double_t VA_2 = p[3];
      Double_t VT_3VA_3 = p[4];
      Double_t VT_4 = p[5];
      Double_t VA_4 = p[6];
      Double_t VT_5VA_5 = p[7];
      Double_t VT_6 = p[8];
      Double_t VA_6 = p[9];

      if (iEP == -1) {
        Double_t value = 1;
        // With EP resolution. Note: R3 already implicitly included in v3  (Not any more)
        // Should the event plane resolutions event be included here?
       // value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi) * fEPRes[1];
       // value += 2. * VT_3VA_3    * TMath::Cos(3.*fLocalDPhi) * fEPRes[2]; 
       // value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi) * fEPRes[3];
        value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi) * fEPRes[1];
        value += 2. * VT_3VA_3    * TMath::Cos(3.*fLocalDPhi) * fEPRes[2]; 
        value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi) * fEPRes[3];
        return B / 3. * value;
      }

      Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
      Double_t fC_S   = TMath::Pi() * kC_S[iEP];

      Double_t value = 1;

      value += 2. * funct_RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4)  *  VA_2      * TMath::Cos(2.*fLocalDPhi);
      value += 2. * funct_RPF_VR_3(fPhi_S,fC_S,VT_2,VT_3VA_3,VT_4,VT_6) * TMath::Cos(3.*fLocalDPhi);
      value += 2. * funct_RPF_VR_4(fPhi_S,fC_S,VT_2,VT_4)  *  VA_4      * TMath::Cos(4.*fLocalDPhi);

      value = 2. * fC_S / TMath::Pi() * B * funct_RPF_BR(fPhi_S,fC_S,VT_2,VT_4) * value;
  
      // FIXME check that this is right for mid plane
      if (iEP == 1 && kDoubleMidPlane) {
        value = 2.0 * value;
      }
      return value;
    }

  private:

};


Double_t RPF_BR(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0) {
	
// from MA:
// den =  1 + 2*v2_t*cos(2*phi)*sin(2*c)/(2*c) + 2*v4_t*cos(4*phi)*sin(4*c)/(4*c)

	Double_t value = 1 + 2. * vt_2 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * kEPRes[1];
	value +=             2. * vt_4 * TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * kEPRes[3]; return value;
}

Double_t RPF_VR_2(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0) {

	Double_t value = vt_2;
	value +=        TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * kEPRes[1];
	value += vt_4 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * kEPRes[1];
	value += vt_2 * TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * kEPRes[3];
	value += vt_4 * TMath::Cos(6.*fPhi_S) * TMath::Sin(6.*fC_S) / (6.*fC_S) * kEPRes[5];

	value = value / RPF_BR(fPhi_S,fC_S,vt_2,vt_4) ;

	return value;
}

Double_t RPF_VR_4(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0) {
	Double_t value = vt_4;
	value +=        TMath::Cos(4.*fPhi_S) * TMath::Sin(4.*fC_S) / (4.*fC_S) * kEPRes[3];
	value += vt_2 * TMath::Cos(2.*fPhi_S) * TMath::Sin(2.*fC_S) / (2.*fC_S) * kEPRes[1];
	value += vt_2 * TMath::Cos(6.*fPhi_S) * TMath::Sin(6.*fC_S) / (6.*fC_S) * kEPRes[5];

	value = value / RPF_BR(fPhi_S,fC_S,vt_2,vt_4) ;

	return value;

}

Double_t TestFunction(Double_t * x, Double_t * par ) {
	Double_t value = 0;

	value = par[0];

	return value;
}

/**
 *  RPF Function for single bin.  Intended to be used for plotting/using final results
 */
//Double_t TaskEventPlane::RPFFunction_Single(Double_t * x, Double_t * par) {
Double_t RPFFunction_Single(Double_t * x, Double_t * par) {
	Double_t fLocalDPhi = x[0];

	Double_t B    = par[0];     //1;
	Double_t VT_2 = par[1];     //0.1;
	Double_t VA_2 = par[2];     //0.2;

	Double_t VT_3VA_3 = par[3]; //0.05;

	Double_t VT_4 = par[4];     //0.05;
	Double_t VA_4 = par[5];     //0.1;

	Int_t iEP = (Int_t) par[6]; // 0,1, or 2, or -1 for All 
	if (iEP == -1) {
	
		Double_t value = 1;
//		value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi);
//		value +=    2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi); 
//		value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi);
  // With EP resolution. Note: R3 already implicitly included in v3
		value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi) * kEPRes[1];
		value +=    2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi); 
		value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi) * kEPRes[3];
		
		return B * value;
	}

	Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
	Double_t fC_S   = TMath::Pi() * kC_S[iEP];

	Double_t value = 1;
	value += 2. * RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4) * VA_2 * TMath::Cos(2.*fLocalDPhi);
	value += 2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi);
	value += 2. * RPF_VR_4(fPhi_S,fC_S,VT_2,VT_4) * VA_4 * TMath::Cos(4.*fLocalDPhi);

  // FIXME check that this is right for mid plane
	if (iEP == 1 && kDoubleMidPlane) {
		value = 2.0 * value;
	}

	value = B * RPF_BR(fPhi_S,fC_S,VT_2,VT_4) * value;

	return value;
}

// FIXME add to TaskEventPlane object, store nbins, single hist range as variables
// what is this hist range stuff I mentioned?


// This code isn't used anymore. see the functors
//Double_t TaskEventPlane::RPFFunction(Double_t * x, Double_t * par ) {
Double_t RPFFunction(Double_t * x, Double_t * par ) {
	// Determine which subset we are in
	//   Determine phi_c, Calculate localDPhi
	Double_t fLocalDPhi = x[0]; 

	// Primary Mode: Nearside in DPhi Only
	Int_t iEP = 0;
	if (x[0] > 0.5 * TMath::Pi()) iEP = 1;
	if (x[0] > 1.5 * TMath::Pi()) iEP = 2;
	// I wonder if I could cheat by not doing the modulus step
	// Temporary thing for 45 bins
//	Int_t iEP = 0;
//	if (x[0] + TMath::Pi()/2. > 2.*TMath::Pi()*(22./45.)) iEP = 1;
//	if (x[0] + TMath::Pi()/2. > 2.*TMath::Pi()*(22.+22.)/45.) iEP = 2;
	
	fLocalDPhi -= iEP*TMath::Pi(); // Post Bin Fix
//	fLocalDPhi -= iEP*2.*TMath::Pi() * 22./45.;
//	printf("x[0] = %f, iEP = %d, fLocalDPhi = %f\n",x[0],iEP,fLocalDPhi);
//	if (iEP == 1) fLocalDPhi -= TMath::Pi(); // can do in one line
//	if (iEP == 2) fLocalDPhi -= 2. * TMath::Pi();

//	if (iEP == 1) fLocalDPhi -= TMath::Pi(); // can do in one line
//	if (iEP == 2) fLocalDPhi -= 2. * TMath::Pi();

	// Secondary Option (that I originally wrote and did not want to delete): All
/*	Int_t iEP = 0;
	if (x[0] > 1.5 * TMath::Pi()) iEP = 1;
	if (x[0] > 3.5 * TMath::Pi()) iEP = 2;
	
	Double_t fLocalDPhi = x[0]; //FIXME
	// I wonder if I could cheat by not doing the modulus step
	if (iEP == 1) fLocalDPhi -= 2. * TMath::Pi(); // can do in one line
	if (iEP == 2) fLocalDPhi -= 4. * TMath::Pi();
*/

	// The Parameters
	Double_t B    = par[0];     //1;
	Double_t VT_2 = par[1];     //0.1;
	Double_t VA_2 = par[2];     //0.2;

	Double_t VT_3VA_3 = par[3]; //0.05;

	Double_t VT_4 = par[4];     //0.05;
	Double_t VA_4 = par[5];     //0.1;

	Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
	Double_t fC_S   = TMath::Pi() * kC_S[iEP];

	Double_t value = 1;
	value += 2. * RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4) * VA_2 * TMath::Cos(2.*fLocalDPhi);
	value += 2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi);
	value += 2. * RPF_VR_4(fPhi_S,fC_S,VT_2,VT_4) * VA_4 * TMath::Cos(4.*fLocalDPhi);

	// FIXME is this right for Mid plane angle range????
	if (iEP == 1 && kDoubleMidPlane) {
		value = 2.0 * value;
	}

// FIXME
//	value = B * value;
	value = B * RPF_BR(fPhi_S,fC_S,VT_2,VT_4) * value;

	return value;
}

//void TaskEventPlane::RPF_Prefit(TF1 * fit,TH1D * fHist, RPF_Functor * funct ) {
void RPF_Prefit(TF1 * fit,TH1D * fHist, RPF_Functor * funct ) {

  double fAverageValue = fHist->Integral() / fHist->GetNbinsX();
  fit->SetParameter(1,fAverageValue);
 // fit->SetParLimits(1,0.8*fAverageValue,1.2*fAverageValue);
 // fit->SetParLimits(1,0.8*fAverageValue,1.2*fAverageValue);
/*
	fit->SetParLimits(2,-0.5,0.5); // v2t
	fit->SetParLimits(3,-0.5,0.5); // v2a
	fit->SetParLimits(4,-0.3,0.3); // v3tv3a
	fit->SetParLimits(5,-0.5,0.5); // v4t
	fit->SetParLimits(6,-0.5,0.5); // v4a
*/
/*
	fit->SetParLimits(2,-kV2T_AbsCut,kV2T_AbsCut); // v2t
	fit->SetParLimits(3,-kV2A_AbsCut,kV2A_AbsCut); // v2a
	fit->SetParLimits(4,-kV3_AbsCut,kV3_AbsCut); // v3tv3a
	fit->SetParLimits(5,-kV4T_AbsCut,kV4T_AbsCut); // v4t
	fit->SetParLimits(6,-kV4A_AbsCut,kV4A_AbsCut); // v4a
*/


//	fit->SetParLimits(1,0,1);
//	fit->SetParLimits(1,-0.03,0.05);
//	fit->SetParLimits(2,-0.03,0.05);
//	fit->SetParLimits(3,-0.05,0.05);
//	fit->SetParLimits(4,-0.02,0.04);
//	fit->SetParLimits(5,-0.02,0.04);

  // Guesses
  /*fit->SetParameter(1,0.05);
  fit->SetParameter(2,0.15);
  fit->SetParameter(3,0.04);
  fit->SetParameter(4,0.05);
  fit->SetParameter(5,0.15);*/
  // Guesses
  fit->SetParameter(2,0.05);
  fit->SetParameter(3,0.10);
  fit->SetParameter(4,0.005);
  fit->SetParameter(5,0.05);
  fit->SetParameter(6,0.10);

/*
  if (funct->GetFixedV2T() > -1) fit->FixParameter(1,funct->GetFixedV2T());
  if (funct->GetFixedV2A() > -1) fit->FixParameter(2,funct->GetFixedV2A());
  if (funct->GetFixedV4T() > -1) fit->FixParameter(4,funct->GetFixedV4T());
  if (funct->GetFixedV4A() > -1) fit->FixParameter(5,funct->GetFixedV4A());
  if (funct->GetFixedV3() > -1) fit->FixParameter(3,funct->GetFixedV3());
  if (funct->GetV2T_Min() > -1) fit->SetParLimits(1,funct->GetV2T_Min(),funct->GetV2T_Max());
  if (funct->GetV2A_Min() > -1) fit->SetParLimits(2,funct->GetV2A_Min(),funct->GetV2A_Max());
  if (funct->GetV4T_Min() > -1) fit->SetParLimits(4,funct->GetV4T_Min(),funct->GetV4T_Max());
  if (funct->GetV4A_Min() > -1) fit->SetParLimits(5,funct->GetV4A_Min(),funct->GetV4A_Max());
  */


  if (funct->GetInitV2T() > -1) fit->SetParameter(2,funct->GetInitV2T());
  if (funct->GetInitV2A() > -1) { 
    fit->SetParameter(3,funct->GetInitV2A());
    printf("Debug in Functor, setting V2A to %f\n",funct->GetInitV2A());
  }
  if (funct->GetInitV3() > -1) fit->SetParameter(4,funct->GetInitV3());
  if (funct->GetInitV4T() > -1) fit->SetParameter(5,funct->GetInitV4T());
  if (funct->GetInitV4A() > -1) fit->SetParameter(6,funct->GetInitV4A());

  if (funct->GetV2T_Min() > -1) fit->SetParLimits(2,funct->GetV2T_Min(),funct->GetV2T_Max());
  if (funct->GetV2A_Min() > -1) fit->SetParLimits(3,funct->GetV2A_Min(),funct->GetV2A_Max());
  if (funct->GetV3_Min() > -1) fit->SetParLimits(4,funct->GetV3_Min(),funct->GetV3_Max());
  if (funct->GetV4T_Min() > -1) fit->SetParLimits(5,funct->GetV4T_Min(),funct->GetV4T_Max());
  if (funct->GetV4A_Min() > -1) fit->SetParLimits(6,funct->GetV4A_Min(),funct->GetV4A_Max());

  if (funct->GetFixedV2T() > -1) fit->FixParameter(2,funct->GetFixedV2T());
  if (funct->GetFixedV2A() > -1) fit->FixParameter(3,funct->GetFixedV2A());
  if (funct->GetFixedV3() > -1) fit->FixParameter(4,funct->GetFixedV3());
  if (funct->GetFixedV4T() > -1) fit->FixParameter(5,funct->GetFixedV4T());
  if (funct->GetFixedV4A() > -1) fit->FixParameter(6,funct->GetFixedV4A());

  
  // Higher order parameters (5, 6)
  // FIXME set the default vn???
  // FIXME preprocess the flowVNModes in main
  switch (fit->GetNpar()) {
    case 10: // V6A
      fit->SetParameter(9,0.0);
      fit->SetParLimits(9,-0.13,0.13);
 //     fit->FixParameter(9,0.0);
      //if (iFlowV6AMode == 0) fit->FixParameter(7,0.0);
      /*if (iFlowV6AMode == 2) {
        if (funct->GetV4T_Min() > -1) fit->SetParLimits(5,funct->GetV4T_Min(),funct->GetV4T_Max());
        if (funct->GetV4A_Min() > -1) fit->SetParLimits(6,funct->GetV4A_Min(),funct->GetV4A_Max());
      }*/
    case 9: // V6T
      fit->SetParameter(8,0.001);
      fit->SetParLimits(8,-0.09,0.09);
//      fit->FixParameter(8,0.0);
      //if (iFlowV6TMode == 0) fit->FixParameter(8,0.0);
    case 8: // V_5
      fit->SetParameter(7,0.0);
      fit->SetParLimits(7,-0.-9,0.09);
      //fit->FixParameter(7,0.0);
      //if (iFlowV5Mode == 0) fit->FixParameter(7,0.0);
      break;
    case 7:
    default:
      break;
  }


  if (funct->GetFixedV5() > -1) fit->FixParameter(7,funct->GetFixedV5());
  if (funct->GetFixedV6T() > -1) fit->FixParameter(8,funct->GetFixedV6T());
  if (funct->GetFixedV6A() > -1) fit->FixParameter(9,funct->GetFixedV6A());

  // FIXME test
 // fit->FixParameter(3,0.0);
  //fit->FixParameter(4,0.0);
  //fit->FixParameter(5,0.0); 

// FIXME loosening cuts. The final analysis with a realistic trigger distribution should stay within normal limits
//	fit->SetParLimits(1,-10.,10.);
//	fit->SetParLimits(2,-10.,10.);
//	fit->SetParLimits(3,-10.,10.);
//	fit->SetParLimits(4,-10.,10.);
//	fit->SetParLimits(5,-10.,10.);

	return;
}

//TF1 * TaskEventPlane::FitRPF(TH1D * fHist, RPF_Functor * fFit, TString fName, Double_t fV2T_Fixed) {
TF1 * FitRPF(TH1D * fHist, RPF_Functor * fFit, TString fName, Double_t fV2T_Fixed) {
	//Int_t nPar = 7;
	Int_t nPar = 10;
	Double_t Min = fHist->GetXaxis()->GetXmin();
	Double_t Max = fHist->GetXaxis()->GetXmax();

	//TF1 * fit = new TF1(Form("%s_Fit",fName.Data()),TaskEventPlane::RPFFunction,Min,Max,nPar); 
	TF1 * fit = new TF1(Form("%s_Fit",fName.Data()),fFit,Min,Max,nPar); 
	fit->SetNpx(300);	

// ====
// Test
//	nPar = 1;
//	TF1 * fit = new TF1(Form("%s_Fit",fName.Data()),TestFunction,Min,Max,nPar); 
// ====

	/*fit->SetParName(0,"B");
	fit->SetParName(1,"v^{t}_{2}");
	fit->SetParName(2,"v^{a}_{2}");
	fit->SetParName(3,"v^{t}_{3}v^{a}_{3}");
	fit->SetParName(4,"v^{t}_{4}");
	fit->SetParName(5,"v^{a}_{4}");*/

  fit->SetParName(0,"EventPlanePar");
	fit->SetParName(1,"B");
	fit->SetParName(2,"v^{t}_{2}");
	fit->SetParName(3,"v^{a}_{2}");
	fit->SetParName(4,"v^{t}_{3}v^{a}_{3}");
	fit->SetParName(5,"v^{t}_{4}");
	fit->SetParName(6,"v^{a}_{4}");
  switch (fit->GetNpar()) {
    case 10:
      fit->SetParName(9,"v^{a}_{6}");
    case 9:
      fit->SetParName(8,"v^{t}_{6}");
    case 8:
      fit->SetParName(7,"v^{t}_{5}v^{a}_{5}");
      break;
    case 7:
    default:
      break;
  }


  fit->FixParameter(0,0.0); // Parameter not used yet

	RPF_Prefit(fit,fHist,fFit);

 // fFit->DebugPrint();

	// Add a switch for this
//	fit->FixParameter(1,0.01);
//	fit->FixParameter(3,0.);
//	fit->FixParameter(4,0.);
//	fit->FixParameter(5,0.);

  if (fV2T_Fixed > -1) {
    printf("Fixing v^{t}_{2} = %f\n",fV2T_Fixed);
//    fit->FixParameter(1,fV2T_Fixed);
    fit->FixParameter(2,fV2T_Fixed);
  }

	fHist->Fit(fit,"0M");

	//fit->SetLineColor(kFitLineColor);

	return fit;
}

/**
  * Code to merge the near side bins of the event plane bins
  */
//TH1D * TaskEventPlane::MergeEvtPlanesForRPF(vector<TH1D *> fHists, TString fName) {
TH1D * MergeEvtPlanesForRPF(vector<TH1D *> fHists, TString fName) {
	Int_t nHists = fHists.size();
	if (nHists==0) return 0;
	// Assume all histograms have same binning
	Int_t nSingleBins = fHists[0]->GetNbinsX();

	Bool_t bNBinsOdd = (Bool_t) nSingleBins % 2;
	if (bNBinsOdd) printf("DPhi Histos have an odd number of bins\n");
	else           printf("DPhi Histos have an even number of bins\n");

	Int_t nUsedBins = nSingleBins / 2; // factor of two for removing awayside

//	Int_t nZeroBin = fHists[0]->FindBin(0.);
//	printf(" for zero bin %d, edges %f , %f\n",nZeroBin,fHists[0]->GetXaxis()->GetBinLowEdge(nZeroBin),fHists[0]->GetXaxis()->GetBinUpEdge(nZeroBin));

	Int_t nFinalBins = nHists * nUsedBins;
//	Double_t fSingleMin = fHists[0]->GetXaxis()->GetXmin();
//	Double_t fSingleMax = fHists[0]->GetXaxis()->GetXmax();
//	Double_t fSingleRange = fSingleMax - fSingleMin; // removeable

	Double_t fSingleMin   = fHists[0]->GetXaxis()->GetBinLowEdge(1);
	Double_t fSingleMax   = fHists[0]->GetXaxis()->GetBinUpEdge(nUsedBins);
	Double_t fSingleRange = fSingleMax - fSingleMin; // removeable

	printf("nSingleBins = %d, nUsedBins = %d, nFinalBins = %d\n",nSingleBins,nUsedBins,nFinalBins);
	printf("Final Range pi*(%f) - pi*(%f), for a range of pi*(%f)\n",fSingleMin/TMath::Pi(),fSingleMax/TMath::Pi(),fSingleRange/TMath::Pi());

	Double_t fFinalMin = fSingleMin; // removeable
	Double_t fFinalMax = fSingleMin + nHists * fSingleRange; //removeable

	printf("output histo range: pi*(%f) - pi*(%f) for a range of pi*(%f)\n",fFinalMin/TMath::Pi(),fFinalMax/TMath::Pi(),(fFinalMax-fFinalMin)/TMath::Pi());

	TH1D * fOutputHist = new TH1D(fName.Data(),fName.Data(),nFinalBins,fFinalMin,fFinalMax);
	fOutputHist->GetXaxis()->SetTitle(fHists[0]->GetXaxis()->GetTitle());

	for (Int_t i = 0; i < nHists; i++) {
		for (Int_t j = 1; j <= nUsedBins; j++) {
			Int_t iFinalBin = j + i * nUsedBins;
			fOutputHist->SetBinContent(iFinalBin,fHists[i]->GetBinContent(j));
			fOutputHist->SetBinError(iFinalBin,fHists[i]->GetBinError(j));
		}
	}

	return fOutputHist;
}




#endif
