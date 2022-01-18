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

#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

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





// Multivariate Distribution Sampler
// Copied from https://root.cern/doc/master/multidimSampling_8C.html

// Adaptions: just take the Covariance matrix through the parameters
// I already have covariance matrices, no need to recalculate from
// the correlation matrices

bool debug = false;

// Define the GausND strcture
struct GausND {

   TVectorD X;
   TVectorD Mu;
   TMatrixDSym CovMat;

   GausND( int dim ) :
      X(TVectorD(dim)),
      Mu(TVectorD(dim)),
      CovMat(TMatrixDSym(dim) )
   {}
   double operator() (double *x, double *p) {
      // 4 parameters
      int dim = X.GetNrows();
      int k = 0;
      for (int i = 0; i<dim; ++i) { X[i] = x[i] - p[k]; k++; }
// Old code assuming 2nd N parameters are the sigmas
//  and the final N(N-1)/2 parameters are the correlation matrix
/*
      for (int i = 0; i<dim; ++i) {
         CovMat(i,i) = p[k]*p[k];
         k++;
      }
      for (int i = 0; i<dim; ++i) {
         for (int j = i+1; j<dim; ++j) {
            // p now are the correlations N(N-1)/2
               CovMat(i,j) = p[k]*sqrt(CovMat(i,i)*CovMat(j,j));
               CovMat(j,i) = CovMat(i,j);
               k++;
         }
      }
*/
      for (int i = 0; i<dim; ++i) {
         // p are the sigmas of the parameters
         CovMat(i,i) = p[k]*p[k];
         k++;
      }
      for (int i = 0; i<dim; ++i) {
         for (int j = i+1; j<dim; ++j) {
            //// p now are the correlations N(N-1)/2
            // Now p are the off-diagonal covariance elements
               CovMat(i,j) = p[k];
               CovMat(j,i) = CovMat(i,j);
               k++;
         }
      }
      if (debug) {
         X.Print();
         CovMat.Print();
      }

      double det = CovMat.Determinant();
      if (det <= 0) {
         //Fatal("GausND","Determinant is <= 0 det = %e",det);
         //CovMat.Print();
         //return 0;

         //FIXME
         //printf("Experiment: det %e < 0, taking abs\n",det);
         det = -det;
      }
      double norm = std::pow( 2. * TMath::Pi(), dim/2) * sqrt(det);
      // compute the gaussians
      CovMat.Invert();
      double fval  = std::exp( - 0.5 * CovMat.Similarity(X) )/ norm;

      if (debug) {
         std::cout << "det  " << det << std::endl;
         std::cout << "norm " << norm << std::endl;
         std::cout << "fval " << fval << std::endl;
      }

      return fval;
   }
};






















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
    double GetEPRes(int i) { return fEPRes[i]; }

    double GetNumTriggers() { return fNumTriggers; }
    void SetNumTriggers(double input) { fNumTriggers = input; }

    double GetAverageValue() { return fAverageValue; }
    void SetAverageValue(double input) {fAverageValue = input; }

    void SetInitB(double input) {fInitB = input;}
    void SetInitV1(double input) {fInitV1 = input;}
    void SetInitV2T(double input) {fInitV2T = input;}
    void SetInitV2A(double input) {fInitV2A = input;}
    void SetInitV3(double input) {fInitV3 = input;}
    void SetInitV4T(double input) {fInitV4T = input;}
    void SetInitV4A(double input) {fInitV4A = input;}
    void SetInitV5(double input) {fInitV5 = input;}
    void SetInitV6T(double input) {fInitV6T = input;}
    void SetInitV6A(double input) {fInitV6A = input;}

    void SetFixedB(double input) {fFixedB = input;}
    void SetFixedV1(double input) {fFixedV1 = input;}
    void SetFixedV2T(double input) {fFixedV2T = input;}
    void SetFixedV2A(double input) {fFixedV2A = input;}
    void SetFixedV3(double input) {fFixedV3 = input;}
    void SetFixedV4T(double input) {fFixedV4T = input;}
    void SetFixedV4A(double input) {fFixedV4A = input;}
    void SetFixedV5(double input) {fFixedV5 = input;}
    void SetFixedV6T(double input) {fFixedV6T = input;}
    void SetFixedV6A(double input) {fFixedV6A = input;}

    void SetBRange(double min, double max) { fB_Min = min; fB_Max = max; }
    void SetV1Range(double min, double max) { fV1_Min = min; fV1_Max = max; }
    void SetV2TRange(double min, double max) { fV2T_Min = min; fV2T_Max = max; }
    void SetV2ARange(double min, double max) { fV2A_Min = min; fV2A_Max = max; }
    void SetV3Range(double min, double max) { fV3_Min = min; fV3_Max = max; }
    void SetV4TRange(double min, double max) { fV4T_Min = min; fV4T_Max = max; }
    void SetV4ARange(double min, double max) { fV4A_Min = min; fV4A_Max = max; }
    void SetV5Range(double min, double max) { fV5_Min = min; fV5_Max = max; }
    void SetV6TRange(double min, double max) { fV6T_Min = min; fV6T_Max = max; }
    void SetV6ARange(double min, double max) { fV6A_Min = min; fV6A_Max = max; }

    double GetInitB()  { return fInitB;  }
    double GetInitV1()  { return fInitV1;  }
    double GetInitV2T() { return fInitV2T; }
    double GetInitV2A() { return fInitV2A; }
    double GetInitV3()  { return fInitV3;  }
    double GetInitV4T() { return fInitV4T; }
    double GetInitV4A() { return fInitV4A; }
    double GetInitV5()  { return fInitV5;  }
    double GetInitV6T() { return fInitV6T; }
    double GetInitV6A() { return fInitV6A; }

    double GetFixedB()  { return fFixedB;  }
    double GetFixedV1()  { return fFixedV1;  }
    double GetFixedV2T() { return fFixedV2T; }
    double GetFixedV2A() { return fFixedV2A; }
    double GetFixedV3()  { return fFixedV3;  }
    double GetFixedV4T() { return fFixedV4T; }
    double GetFixedV4A() { return fFixedV4A; }
    double GetFixedV5()  { return fFixedV5;  }
    double GetFixedV6T() { return fFixedV6T; }
    double GetFixedV6A() { return fFixedV6A; }

    double GetB_Min()   { return fB_Min; }
    double GetB_Max()   { return fB_Max; }

    double GetV1_Min()   { return fV1_Min; }
    double GetV1_Max()   { return fV1_Max; }

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
      Double_t VT_1VA_1 = p[2];
      Double_t VT_2 = p[3];
      Double_t VA_2 = p[4];
      Double_t VT_3VA_3 = p[5];
      Double_t VT_4 = p[6];
      Double_t VA_4 = p[7];
      Double_t VT_5VA_5 = p[8];
      Double_t VT_6 = p[9];
      Double_t VA_6 = p[10];
      

      Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
      Double_t fC_S   = TMath::Pi() * kC_S[iEP];

      Double_t value = 1;
      value += 2. * VT_1VA_1 * TMath::Cos(fLocalDPhi);
      value += 2. * funct_RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4,VT_6) * VA_2 * TMath::Cos(2.*fLocalDPhi);
      value += 2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi);
      //value += 2. * funct_RPF_VR_3(fPhi_S,fC_S,VT_2,VT_3VA_3,VT_4,VT_6) * TMath::Cos(3.*fLocalDPhi);
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

    double fNumTriggers = -1; // Number of triggers (useful for changing normalization)
    double fAverageValue = -1; // The average value of the input histogram. May be used to renormalize the B parameter

    double fInitB  = -1;
    double fInitV1  = -1;
    double fInitV2T = -1;
    double fInitV2A = -1;
    double fInitV3  = -1;
    double fInitV4T = -1;
    double fInitV4A = -1;
    double fInitV5  = -1;
    double fInitV6T = -1;
    double fInitV6A = -1;

    double fFixedB  = -1; 
    double fFixedV1  = -1;
    double fFixedV2T = -1;
    double fFixedV2A = -1;
    double fFixedV3  = -1;
    double fFixedV4T = -1;
    double fFixedV4A = -1;
    double fFixedV5  = -1;
    double fFixedV6T = -1;
    double fFixedV6A = -1;
    
    double fB_Min = -1;
    double fB_Max = -1;
    double fV1_Min = -1;
    double fV1_Max = -1;
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

/**
  * This functor is to be used when splitting into individual event plane
  * bins. 
  */

class RPF_Functor_Single : public RPF_Functor {
  public:
    
    double operator() (double *x, double *p) {
      Double_t fLocalDPhi = x[0];

      Int_t iEP = (Int_t) p[0]; // 0,1, or 2, or -1 for All 
      Double_t B    = p[1];
      Double_t VT_1VA_1 = p[2];
      Double_t VT_2 = p[3];
      Double_t VA_2 = p[4];
      Double_t VT_3VA_3 = p[5];
      Double_t VT_4 = p[6];
      Double_t VA_4 = p[7];
      Double_t VT_5VA_5 = p[8];
      Double_t VT_6 = p[9];
      Double_t VA_6 = p[10];

      if (iEP == -1) {
        Double_t value = 1;
        // With EP resolution. Note: R3 already implicitly included in v3  (Not any more)
        // Should the event plane resolutions event be included here?
       // value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi) * fEPRes[1];
       // value += 2. * VT_3VA_3    * TMath::Cos(3.*fLocalDPhi) * fEPRes[2]; 
       // value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi) * fEPRes[3];
//        value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi) * fEPRes[1];
//        value += 2. * VT_3VA_3    * TMath::Cos(3.*fLocalDPhi) * fEPRes[2]; 
//        value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi) * fEPRes[3];
        // For all event plane angles, the EPR (with respect to EP2) don't enter in. Each vNa just multiplies the vNt
        value += 2. * VT_1VA_1    * TMath::Cos(fLocalDPhi); 
        value += 2. * VT_2 * VA_2 * TMath::Cos(2.*fLocalDPhi);
        value += 2. * VT_3VA_3    * TMath::Cos(3.*fLocalDPhi); 
        value += 2. * VT_4 * VA_4 * TMath::Cos(4.*fLocalDPhi);
        return B / 3. * value;
      }

      Double_t fPhi_S = TMath::Pi() * kPhi_S[iEP];
      Double_t fC_S   = TMath::Pi() * kC_S[iEP];

      Double_t value = 1;

      value += 2. * VT_1VA_1    * TMath::Cos(fLocalDPhi); 
      value += 2. * funct_RPF_VR_2(fPhi_S,fC_S,VT_2,VT_4)  *  VA_2      * TMath::Cos(2.*fLocalDPhi);
      //value += 2. * funct_RPF_VR_3(fPhi_S,fC_S,VT_2,VT_3VA_3,VT_4,VT_6) * TMath::Cos(3.*fLocalDPhi);
      value += 2. * VT_3VA_3 * TMath::Cos(3.*fLocalDPhi);
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

// Example parameter Values
//          V2FP_0 = 5.39182e-01;//   8.36928e-03
//          V2FP_1 = 3.70690e+00;//   1.94541e-01
//          V2FP_2 = 1.54252e+00;//   9.17354e-02
//          V4FP_0 = 3.89071e-01;//   2.13836e-02
//          V4FP_1 = 3.43959e+00;//   1.90255e-01
//          V4FP_2 = 1.23961e+00;//   8.31609e-02

/**
  * Prepare a [0]*Landau(x,[1],[2],0) fit for Vn params
  */
void PrepLandauFit(TF1 * fit) {
  fit->SetParName(0,"v_{0}");
  fit->SetParName(1,"v_{1}");
  fit->SetParName(2,"v_{2}");

  fit->SetParLimits(0,0.,5);
  fit->SetParLimits(1,0.,30);
  fit->SetParLimits(2,0.,30);

  fit->SetParameter(0,0.15);
  fit->SetParameter(1,3.5);
  fit->SetParameter(2,2.0);
}

/**
  * Print out the Landau fit parameters in a format that
  * can be easily copied to the event generator code
  */
void PrintLandauFit(TF1 * fit, int n) {
  for (int i = 0; i < 3; i++)
  printf("V%dFP_%d = %e ;// e=%e\n",n,i,fit->GetParameter(i),fit->GetParError(i));


}












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
  //fit->SetParameter(1,fAverageValue);
  // FIXME the factor should be pi/2 for this to make sense
  //double fAverageValueNorm = TMath::PiOver2() * fAverageValue;
  double fAverageValueNorm = TMath::Pi() * fAverageValue;

 // fit->SetParLimits(1,0.8*fAverageValue,1.2*fAverageValue);
 // fit->SetParLimits(1,0.8*fAverageValue,1.2*fAverageValue);

  // Dealing with overall normalization
  
  if (funct->GetB_Min() > -1) fit->SetParLimits(1, fAverageValueNorm * funct->GetB_Min(),fAverageValueNorm *funct->GetB_Max());
  if (funct->GetInitB() > -1) fit->SetParameter(1, fAverageValueNorm * funct->GetInitB());
  if (funct->GetFixedB() > -1) fit->FixParameter(1,fAverageValueNorm * funct->GetFixedB());

  // Guesses
/*  fit->SetParameter(2,0.05);
  fit->SetParameter(3,0.10);
  fit->SetParameter(4,0.005);
  fit->SetParameter(5,0.05);
  fit->SetParameter(6,0.10);
*/

  fit->SetParameter(2,0.0); // V1

  fit->SetParameter(3,0.05); // V2T
  fit->SetParameter(4,0.10); // V2A
  fit->SetParameter(5,0.005); // V3
  fit->SetParameter(6,0.05); // V4T
  fit->SetParameter(7,0.10); // V4A

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


  if (funct->GetInitV1() > -1) fit->SetParameter(5,funct->GetInitV1());
  if (funct->GetInitV2T() > -1) fit->SetParameter(3,funct->GetInitV2T());
  if (funct->GetInitV2A() > -1) { 
    fit->SetParameter(4,funct->GetInitV2A());
    printf("Debug in Functor, setting V2A to %f\n",funct->GetInitV2A());
  }
  if (funct->GetInitV3() > -1) fit->SetParameter(5,funct->GetInitV3());
  if (funct->GetInitV4T() > -1) fit->SetParameter(6,funct->GetInitV4T());
  if (funct->GetInitV4A() > -1) fit->SetParameter(7,funct->GetInitV4A());

  if (funct->GetV1_Min() > -1) fit->SetParLimits(2,funct->GetV1_Min(),funct->GetV1_Max());
  if (funct->GetV2T_Min() > -1) fit->SetParLimits(3,funct->GetV2T_Min(),funct->GetV2T_Max());
  if (funct->GetV2A_Min() > -1) fit->SetParLimits(4,funct->GetV2A_Min(),funct->GetV2A_Max());
  if (funct->GetV3_Min() > -1) fit->SetParLimits(5,funct->GetV3_Min(),funct->GetV3_Max());
  if (funct->GetV4T_Min() > -1) fit->SetParLimits(6,funct->GetV4T_Min(),funct->GetV4T_Max());
  if (funct->GetV4A_Min() > -1) fit->SetParLimits(7,funct->GetV4A_Min(),funct->GetV4A_Max());

  if (funct->GetFixedV1() > -1) fit->FixParameter(2,funct->GetFixedV1());
  if (funct->GetFixedV2T() > -1) fit->FixParameter(3,funct->GetFixedV2T());
  if (funct->GetFixedV2A() > -1) fit->FixParameter(4,funct->GetFixedV2A());
  if (funct->GetFixedV3() > -1) fit->FixParameter(5,funct->GetFixedV3());
  if (funct->GetFixedV4T() > -1) fit->FixParameter(6,funct->GetFixedV4T());
  if (funct->GetFixedV4A() > -1) fit->FixParameter(7,funct->GetFixedV4A());

  
  // Higher order parameters (5, 6)
  // FIXME set the default vn???
  // FIXME preprocess the flowVNModes in main
  switch (fit->GetNpar()) {
    case 10: // V6A
      fit->SetParameter(10,0.0);
      fit->SetParLimits(10,-0.13,0.13);
 //     fit->FixParameter(9,0.0);
      //if (iFlowV6AMode == 0) fit->FixParameter(7,0.0);
      /*if (iFlowV6AMode == 2) {
        if (funct->GetV4T_Min() > -1) fit->SetParLimits(5,funct->GetV4T_Min(),funct->GetV4T_Max());
        if (funct->GetV4A_Min() > -1) fit->SetParLimits(6,funct->GetV4A_Min(),funct->GetV4A_Max());
      }*/
    case 9: // V6T
      fit->SetParameter(9,0.001);
      fit->SetParLimits(9,-0.09,0.09);
//      fit->FixParameter(8,0.0);
      //if (iFlowV6TMode == 0) fit->FixParameter(8,0.0);
    case 8: // V_5
      fit->SetParameter(8,0.0);
      fit->SetParLimits(8,-0.-9,0.09);
      //fit->FixParameter(7,0.0);
      //if (iFlowV5Mode == 0) fit->FixParameter(7,0.0);
      break;
    case 7:
    default:
      break;
  }


  if (funct->GetFixedV5() > -1) fit->FixParameter(8,funct->GetFixedV5());
  if (funct->GetFixedV6T() > -1) fit->FixParameter(9,funct->GetFixedV6T());
  if (funct->GetFixedV6A() > -1) fit->FixParameter(10,funct->GetFixedV6A());

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

// May be easier to return the fit results, and get the TF1 via arugment.
//TFitResultPtr FitRPF(TH1D * fHist, RPF_Functor * fFit, TString fName, Double_t fV2T_Fixed, int OverallMode, TF1 ** fOutputFit) {
TFitResultPtr FitRPF(TH1D * fHist, RPF_Functor * fFit, TString fName, double fV2T_Fixed, double fV4T_Fixed, int OverallMode, TF1 ** fOutputFit) {
	//Int_t nPar = 7;
	Int_t nPar = 11;
	Double_t Min = fHist->GetXaxis()->GetXmin();
	Double_t Max = fHist->GetXaxis()->GetXmax();

  if (OverallMode > 0) nPar = 2;

	//TF1 * fit = new TF1(Form("%s_Fit",fName.Data()),TaskEventPlane::RPFFunction,Min,Max,nPar); 
	TF1 * fit = new TF1(Form("%s_Fit",fName.Data()),fFit,Min,Max,nPar); 
	fit->SetNpx(300);	

  fit->SetParName(0,"EventPlanePar");
	fit->SetParName(1,"B");

  // FIXME adjusting overallMode to do ZYAM.

  if (OverallMode == 0) {
    fit->SetParName(2,"vt1va1");
    fit->SetParName(3,"vt2");
    fit->SetParName(4,"va2");
    fit->SetParName(5,"vt3va3");
    fit->SetParName(6,"vt4");
    fit->SetParName(7,"va4");
    // Old names
    /*fit->SetParName(2,"v^{t}_{1}v^{a}_{1}");
    fit->SetParName(3,"v^{t}_{2}");
    fit->SetParName(4,"v^{a}_{2}");
    fit->SetParName(5,"v^{t}_{3}v^{a}_{3}");
    fit->SetParName(6,"v^{t}_{4}");
    fit->SetParName(7,"v^{a}_{4}");*/
    switch (fit->GetNpar()) {
      case 10:
        fit->SetParName(10,"va6");
      case 9:
        fit->SetParName(9,"vt6");
      case 8:
        fit->SetParName(8,"vt5va5");
        break;
      case 7:
      default:
        break;
    }
  }

  fit->FixParameter(0,0.0); // Parameter 0 only used in RPF Single EP

  if (OverallMode == 2) { // Far Eta

    double fAverage = fHist->Integral("width") / (Max - Min);
    fAverage *= 3.;
    fit->FixParameter(1,fAverage);


    fit->FixParameter(2,0.);
    fit->FixParameter(3,0.);
    fit->FixParameter(4,0.);
    fit->FixParameter(5,0.);
    fit->FixParameter(6,0.);
    fit->FixParameter(7,0.);
    switch (fit->GetNpar()) {
      case 10:
        fit->FixParameter(10,0.);
      case 9:
        fit->FixParameter(9,0.);
      case 8:
        fit->FixParameter(8,0.);
        break;
      case 7:
      default:
        break;
    }
    TFitResultPtr otherFitResult = fHist->Fit(fit,"0MS");
    *fOutputFit = fit;
    return otherFitResult;
    //return fit;
  }
  RPF_Prefit(fit,fHist,fFit);

   // fFit->DebugPrint();

  if (fV2T_Fixed > -1) {
    printf("Fixing v^{t}_{2} = %f\n",fV2T_Fixed);
    fit->FixParameter(3,fV2T_Fixed);
  }
  if (fV4T_Fixed > -1) {
    printf("Fixing v^{t}_{4} = %f\n",fV4T_Fixed);
    fit->FixParameter(6,fV4T_Fixed);
  }

	TFitResultPtr fitResults = fHist->Fit(fit,"0MS");
  *fOutputFit = fit;
  cout<<"Just tried to return the fit function"<<endl;
	//return fit;
  return fitResults;
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

/** Raise each point in a TGraph 
  * When dealing with negative values: apply to abs, then conserve sign
  */
void ApplyPowerToTGraph(TGraphErrors * fInput, double exponent) {

  int nPoints = fInput->GetN();
  for (int i = 0; i < nPoints; i++) {
    double initX     = fInput->GetX()[i];
    double initX_Err = fInput->GetEX()[i];
    double initY     = fInput->GetY()[i];
    double initY_Err = fInput->GetEY()[i];

    double finalY = 0;
    double finalY_Err = 0;
    if (initY > 0) {
      finalY = TMath::Power(TMath::Abs(initY),exponent);
      finalY_Err = finalY * exponent * initY_Err/ initY;
    } else {
      finalY = initY;
      finalY_Err = initY_Err;
    }

    fInput->SetPoint(i,initX,finalY);
    fInput->SetPointError(i,initX_Err,finalY_Err);
  }
}

/** Divide TGraphErrors fInput by fDenom
  * Note that fInput is updated with the resulting values
  */
void DivideTGraphs(TGraphErrors * fInput, TGraphErrors * fDenom) {

  // What if numbers don't match?
  
  int nPoints = fInput->GetN();
  for (int i = 0; i < nPoints; i++) {
    double initX     = fInput->GetX()[i];
    double initX_Err = fInput->GetEX()[i];
    double initY     = fInput->GetY()[i];
    double initY_Err = fInput->GetEY()[i];

    double initYD     = fDenom->GetY()[i];
    double initYD_Err = fDenom->GetEY()[i];

    double finalY = 0;
    double finalY_Err = 0;

    if (initYD != 0) {
      finalY = initY / initYD;
      if (initY != 0) {
        finalY_Err = TMath::Abs(finalY) * TMath::Sqrt(TMath::Power(initY_Err/initY,2.) + TMath::Power(initYD_Err/initYD,2.));
      }
    }
    fInput->SetPoint(i,initX,finalY);
    fInput->SetPointError(i,initX_Err,finalY_Err);
  }
}
/** Shift the input TGraph Errors points in Y 
  * Shift by error * fErrorShift
  */
void ShiftTGraphByErr(TGraphErrors * fInput, double fErrorShift) {
  for (int i = 0; i < fInput->GetN(); i++) {
    double initX     = fInput->GetX()[i];
    double initX_Err = fInput->GetEX()[i];
    double initY     = fInput->GetY()[i];
    double initY_Err = fInput->GetEY()[i];

    double finalY = initY + fErrorShift * initY_Err;
    double finalY_Err = initY_Err;

    fInput->SetPoint(i,initX,finalY);
    fInput->SetPointError(i,initX_Err,finalY_Err);
  }
}

/**
  * Scale a TGraph by the given scalar and error
  */
void MultiplyTGraphByScalar(TGraphErrors * fInput, double fScalar, double fError = 0.0) {
  for (int i = 0; i < fInput->GetN(); i++) {
    double initX     = fInput->GetX()[i];
    double initX_Err = fInput->GetEX()[i];
    double initY     = fInput->GetY()[i];
    double initY_Err = fInput->GetEY()[i];

    double finalY = initY * fScalar;
    double finalY_Err = 0;
    if (initY != 0) {
      finalY_Err = TMath::Abs(finalY) * TMath::Sqrt(TMath::Power(initY_Err / initY,2.)+TMath::Power(fError/fScalar,2.));
    }

    printf("finalY = %f, finalY_Err = %f; initY = %f, initY_Err = %f; fScalar = %f, fError = %f\n",finalY,finalY_Err,initY,initY_Err,fScalar,fError);

    fInput->SetPoint(i,initX,finalY);
    fInput->SetPointError(i,initX_Err,finalY_Err);
  }
}


#endif
