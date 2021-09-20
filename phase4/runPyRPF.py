#! /usr/bin/env python3

import sys

from typing import Any, Dict, Tuple, Type, TYPE_CHECKING

from ROOT import gROOT, TCanvas, TGraphErrors, TMultiGraph, TFile, TLegend

import logging

import math

#from reaction_plane_fit import base
#import pachyderm.fit as fit_base
from pachyderm.fit.base import calculate_function_errors, FitFailed, BaseFitResult, FitResult
from pachyderm import histogram

from reaction_plane_fit import base
from reaction_plane_fit import three_orientations
from reaction_plane_fit import plot


# event plane information
# From QnVector Framework, on MB
fEPRes_Set_0 = [[0.765960,  0.619163,  0.509267,  0.348666,  0.318429,  0.187868],
                [0.858157,  0.822691,  0.692985,  0.580624,  0.502229,  0.375755],
                [0.832549,  0.771133,  0.639423,  0.507014,  0.439729,  0.305388],
                [0.704550,  0.445893,  0.380824,  0.196809,  0.211605,  0.084895]]

# From QnVector Framework, on EGA
fEPRes_Set_1 = [[0.767842, 0.622406, 0.525535, 0.390981, 0.352418, 0.236852],
                [0.781089, 0.651664, 0.549970, 0.423973, 0.376066, 0.259773],
                [0.779045, 0.646402, 0.544922, 0.414330, 0.369320, 0.252558],
                [0.762544, 0.610570, 0.515899, 0.376056, 0.343955, 0.226055]]



# Perfect EPR, useful for MC
fEPRes_Set_3 = [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]


#ParamNames=["B","V2T","V2A","V3","V4T","V4A"]
#PyParamNames=["B","v2_t","v2_a","v3","v4_t","v4_a"]
ParamNames=["B","V1","V2T","V2A","V3","V4T","V4A"]
PyParamNames=["B","v1","v2_t","v2_a","v3","v4_t","v4_a"]

PyColor=8
PyMarkerStyle=26

PyRPDepColor=46
PyRPDepMarkerStyle=25

CColor=9
CMarkerStyle=22

# Enabling different RPF modes
enableInclusiveFit=False
enableRPDepFit=False
enableReduxFit=False
useMinos=True

#3
nSkipLast=7
#nSkipLast=8
#MCRescaleValue=1e6
# get MCRescaleValue from CTask


def FixErrorsOnHist(hist):
  NumZeroBins = 0
  for i in range(hist.GetNbinsX()):
    value=hist.GetBinContent(i)
    if (value == 0):
      print("Found a zero bin")
      hist.SetBinError(i,1.)
      NumZeroBins+=1
      #FIXME
  print("Found %d Zero bins (out of %d)" % (NumZeroBins,hist.GetNbinsX()))

# could get a lot of information from the C++ task

#def RunRPFCode(fInputFileAll,fInputFile_EP_0,fInputFile_EP_1,fInputFile_EP_2):
def RunRPFCode(fCTask,fOutputDir,fOutputFile):

#  logging.basicConfig(level="DEBUG")


  fNumTriggers=fCTask.GetNumTriggers()
  print("Found the number of triggers = %f" % (fNumTriggers))
  fObservable=fCTask.GetObservable()
  nObsBins=fCTask.GetNObsBins()

  fFlowTermMode=fCTask.GetFlowTermModeAssoc()
  nFixV40Last=fCTask.GetFixV4Threshold()

  iV1Mode=fCTask.GetFlowV1Mode()
  iV3Mode=fCTask.GetFlowV3Mode()
  iV5Mode=fCTask.GetFlowV5Mode()
  iV6TMode=fCTask.GetFlowV6TMode()
  iV6AMode=fCTask.GetFlowV6AMode()


  fObsBins=[0, 1./6, 2./6, 3./6, 4./6, 5./6, 1.]
  # FIXME do the new zt bins
  if fObservable==0:
    print("Error: trigger pt is the observable, haven't coded bins in yet")
  if fObservable==1:
    print("Using z_T as the observable")
  if fObservable==2:
    fObsBins=[0.2,0.4,0.8,1.5,2.5,4,7,11,17]
    #fObsBins=[0.15,0.4,0.8,1.45,2.5,4.2,6.95,11.4,18.6]
    print("Using associated pt as the observable")
    # assoc pt

  print("This analysis is in observable %d, with %d bins" % (fObservable,nObsBins))
  iCentBin=fCTask.GetCentBin()
  print(" Centrality Bin %d" % (iCentBin))

  iEPRSet = fCTask.GetEPRSet()

  fUseEPRSet = fEPRes_Set_0
  if iEPRSet == 1:
    fUseEPRSet = fEPRes_Set_1
  if iEPRSet == 3:
    fUseEPRSet = fEPRes_Set_3

  res_par = {"R22": fUseEPRSet[iCentBin][1], "R42" : fUseEPRSet[iCentBin][3], "R62": fUseEPRSet[iCentBin][5], "R82": 0.0}
 # res_par = {"R22": fEPRes_Set_0[iCentBin][1], "R42" : fEPRes_Set_0[iCentBin][3], "R62": fEPRes_Set_0[iCentBin][5], "R82": 0.0}


  print("Resolution parameters:")
  print(res_par)

#  enableInclusiveFit=False
#  enableRPDepFit=False
#  enableReduxFit=False
#  useMinos=False

  nRebin=1
  MCRescale=-1
  if (fCTask.GetMCGenMode()):
    MCRescale=fCTask.GetMCRescaleFactor()
#    MCRescale=MCRescaleValue
    nRebin=2*nRebin


  CovMatrices = []

  # Initializing some TGraphs
#  Py_B_TGraph: TGraphErrors
  # Graphs for basic RPF
  Py_ChiSq_TGraph = TGraphErrors(nObsBins)
  Py_ChiSq_TGraph.SetName("Py_ChiSq_TGraph")

  Py_B_TGraph   = TGraphErrors(nObsBins)
  Py_B_TGraph.SetName("Py_B_TGraph")
  Py_V1_TGraph  = TGraphErrors(nObsBins)
  Py_V1_TGraph.SetName("Py_V1_TGraph")
  Py_V2T_TGraph = TGraphErrors(nObsBins)
  Py_V2T_TGraph.SetName("Py_V2T_TGraph") 
  Py_V2A_TGraph = TGraphErrors(nObsBins)
  Py_V2A_TGraph.SetName("Py_V2A_TGraph")
  Py_V3_TGraph  = TGraphErrors(nObsBins)
  Py_V3_TGraph.SetName("Py_V3_TGraph")
  Py_V4T_TGraph = TGraphErrors(nObsBins)
  Py_V4T_TGraph.SetName("Py_V4T_TGraph")
  Py_V4A_TGraph = TGraphErrors(nObsBins)
  Py_V4A_TGraph.SetName("Py_V4A_TGraph")

  #Py_TGraphs= [Py_B_TGraph, Py_V2T_TGraph, Py_V2A_TGraph, Py_V3_TGraph, Py_V4T_TGraph, Py_V4A_TGraph]
  Py_TGraphs= [Py_B_TGraph, Py_V1_TGraph, Py_V2T_TGraph, Py_V2A_TGraph, Py_V3_TGraph, Py_V4T_TGraph, Py_V4A_TGraph]

  # Graphs for RPDep Fit
  Py_RPDep_ChiSq_TGraph = TGraphErrors(nObsBins)
  Py_RPDep_ChiSq_TGraph.SetName("Py_RPDep_ChiSq_TGraph")
  Py_RPDep_B_TGraph   = TGraphErrors(nObsBins)
  Py_RPDep_B_TGraph.SetName("Py_RPDep_B_TGraph")
  Py_RPDep_V1_TGraph  = TGraphErrors(nObsBins)
  Py_RPDep_V1_TGraph.SetName("Py_RPDep_V1_TGraph")
  Py_RPDep_V2T_TGraph = TGraphErrors(nObsBins)
  Py_RPDep_V2T_TGraph.SetName("Py_RPDep_V2T_TGraph")
  Py_RPDep_V2A_TGraph = TGraphErrors(nObsBins)
  Py_RPDep_V2A_TGraph.SetName("Py_RPDep_V2A_TGraph")
  Py_RPDep_V3_TGraph  = TGraphErrors(nObsBins)
  Py_RPDep_V3_TGraph.SetName("Py_RPDep_V3_TGraph")
  Py_RPDep_V4T_TGraph = TGraphErrors(nObsBins)
  Py_RPDep_V4T_TGraph.SetName("Py_RPDep_V4T_TGraph")
  Py_RPDep_V4A_TGraph = TGraphErrors(nObsBins)
  Py_RPDep_V4A_TGraph.SetName("Py_RPDep_V4A_TGraph")
  # Reaction Plane Dependent Signal fit unique parameters
  # Yields
  Py_RPDep_IP_YieldNS = TGraphErrors(nObsBins)
  Py_RPDep_IP_YieldNS.SetName("Py_RPDep_IP_YieldNS")
  Py_RPDep_MP_YieldNS = TGraphErrors(nObsBins)
  Py_RPDep_MP_YieldNS.SetName("Py_RPDep_MP_YieldNS")
  Py_RPDep_OP_YieldNS = TGraphErrors(nObsBins)
  Py_RPDep_OP_YieldNS.SetName("Py_RPDep_OP_YieldNS")

  Py_RPDep_IP_YieldAS = TGraphErrors(nObsBins)
  Py_RPDep_IP_YieldAS.SetName("Py_RPDep_IP_YieldAS")
  Py_RPDep_MP_YieldAS = TGraphErrors(nObsBins)
  Py_RPDep_MP_YieldAS.SetName("Py_RPDep_MP_YieldAS")
  Py_RPDep_OP_YieldAS = TGraphErrors(nObsBins)
  Py_RPDep_OP_YieldAS.SetName("Py_RPDep_OP_YieldAS")

  # Sigmas
  Py_RPDep_IP_SigmaNS = TGraphErrors(nObsBins)
  Py_RPDep_IP_SigmaNS.SetName("Py_RPDep_IP_SigmaNS")
  Py_RPDep_MP_SigmaNS = TGraphErrors(nObsBins)
  Py_RPDep_MP_SigmaNS.SetName("Py_RPDep_MP_SigmaNS")
  Py_RPDep_OP_SigmaNS = TGraphErrors(nObsBins)
  Py_RPDep_OP_SigmaNS.SetName("Py_RPDep_OP_SigmaNS")

  Py_RPDep_IP_SigmaAS = TGraphErrors(nObsBins)
  Py_RPDep_IP_SigmaAS.SetName("Py_RPDep_IP_SigmaAS")
  Py_RPDep_MP_SigmaAS = TGraphErrors(nObsBins)
  Py_RPDep_MP_SigmaAS.SetName("Py_RPDep_MP_SigmaAS")
  Py_RPDep_OP_SigmaAS = TGraphErrors(nObsBins)
  Py_RPDep_OP_SigmaAS.SetName("Py_RPDep_OP_SigmaAS")

  Py_RPDep_TGraphs= [Py_RPDep_B_TGraph, Py_RPDep_V1_TGraph, Py_RPDep_V2T_TGraph, Py_RPDep_V2A_TGraph, Py_RPDep_V3_TGraph, Py_RPDep_V4T_TGraph, Py_RPDep_V4A_TGraph]

  # Load the Vn TGraphs for triggers and tracks
  FlowV2TGraph=fCTask.GetTriggerV2()
  FlowV2AGraph=fCTask.GetTrackV2()

  FlowV4TGraph=fCTask.GetTriggerV4()
  FlowV4AGraph=fCTask.GetTrackV4()


  for iObsBin in range(nObsBins):
    print("Doing the thing for %s bin %d" % (fCTask.GetObservableName(),iObsBin))
    ObsBinCenter=0.5*(fObsBins[iObsBin]+fObsBins[iObsBin+1])
    ObsBinWidth=0.5*(fObsBins[iObsBin+1]-fObsBins[iObsBin])


    LocalLabel="%.1f $\\leq p_{T}^{a} <$ %.1f GeV/$c$" % (fObsBins[iObsBin],fObsBins[iObsBin+1])

    print("Loading functor")
    RPFFunctor = fCTask.GetRPFFunctor(iObsBin)
    print("Success")


    UseLogLikelihood=False
    if (iObsBin > nSkipLast):
      continue
    #  UseLogLikelihood=True

    #  Clone these all first so they are not alterred in the original program

    SigInPlaneHistOrig=fCTask.GetNearEtaDPhiProjEP(iObsBin,0)
    SigMidPlaneHistOrig=fCTask.GetNearEtaDPhiProjEP(iObsBin,1)
    SigOutPlaneHistOrig=fCTask.GetNearEtaDPhiProjEP(iObsBin,2)
    SigInclusiveHistOrig=fCTask.GetNearEtaDPhiProjAll(iObsBin)

    BkgInPlaneHistOrig=fCTask.GetFarEtaDPhiProjEP(iObsBin,0)
    BkgMidPlaneHistOrig=fCTask.GetFarEtaDPhiProjEP(iObsBin,1)
    BkgOutPlaneHistOrig=fCTask.GetFarEtaDPhiProjEP(iObsBin,2)
    BkgInclusiveHistOrig=fCTask.GetFarEtaDPhiProjAll(iObsBin)

    SigInPlaneHist=SigInPlaneHistOrig.Clone("%s_Clone" % (SigInPlaneHistOrig.GetName()))
    SigMidPlaneHist=SigMidPlaneHistOrig.Clone("%s_Clone" % (SigMidPlaneHistOrig.GetName()))
    SigOutPlaneHist=SigOutPlaneHistOrig.Clone("%s_Clone" % (SigOutPlaneHistOrig.GetName()))
    SigInclusiveHist=SigInclusiveHistOrig.Clone("%s_Clone" % (SigInclusiveHistOrig.GetName()))
    BkgInPlaneHist=BkgInPlaneHistOrig.Clone("%s_Clone" % (BkgInPlaneHistOrig.GetName()))
    BkgMidPlaneHist=BkgMidPlaneHistOrig.Clone("%s_Clone" % (BkgMidPlaneHistOrig.GetName()))
    BkgOutPlaneHist=BkgOutPlaneHistOrig.Clone("%s_Clone" % (BkgOutPlaneHistOrig.GetName()))
    BkgInclusiveHist=BkgInclusiveHistOrig.Clone("%s_Clone" % (BkgInclusiveHistOrig.GetName()))

    ListOfHists=[SigInPlaneHist,SigMidPlaneHist,SigOutPlaneHist,SigInclusiveHist,BkgInPlaneHist,BkgMidPlaneHist,BkgOutPlaneHist,BkgInclusiveHist]

    # Rescaling by Number of Triggers for numerical betterness
    SigInPlaneHist.Scale(fNumTriggers)
    SigMidPlaneHist.Scale(fNumTriggers)
    SigOutPlaneHist.Scale(fNumTriggers)
    SigInclusiveHist.Scale(fNumTriggers * 3)
    BkgInPlaneHist.Scale(fNumTriggers)
    BkgMidPlaneHist.Scale(fNumTriggers)
    BkgOutPlaneHist.Scale(fNumTriggers)
    BkgInclusiveHist.Scale(fNumTriggers * 3) # to correct for previous graphical downscale

    # FIXME rescale in case of MC
    # Weighting causes all the histograms to have fractional entries even when not divided by num triggers (which also has fractional weighting)
    # if MC is done without reweighting, this would be unnecessary (and wrong)
    # could define a rescale 
    if (MCRescale > 0):
      for hist in ListOfHists:
        hist.Scale(MCRescale)

    if (nRebin > 1):
      for hist in ListOfHists:
        hist.Rebin(nRebin)

    # Fitting just the background
    rp_fit = three_orientations.BackgroundFit(resolution_parameters=res_par,use_log_likelihood=UseLogLikelihood, signal_region=(0,0.8),background_region = (0.8,1.35),use_minos=useMinos)

    # Fitting the background and signal regions. Same yield parameters across RPs?
    rp_fit_IncSig = three_orientations.InclusiveSignalFit(resolution_parameters=res_par,use_log_likelihood=UseLogLikelihood, signal_region=(0,0.8),background_region = (0.8,1.35),use_minos=useMinos)

    rp_fit_RPSig = three_orientations.SignalFit(resolution_parameters=res_par,use_log_likelihood=UseLogLikelihood, signal_region=(0,0.8),background_region = (0.8,1.35),use_minos=useMinos)

    rp_fit_Redux = three_orientations.BackgroundFit(resolution_parameters=res_par,use_log_likelihood=UseLogLikelihood, signal_region=(0,0.8),background_region = (0.8,1.35),use_minos=useMinos)

    dataBkg = {"background": {"in_plane" : BkgInPlaneHist, "mid_plane" : BkgMidPlaneHist, "out_of_plane" : BkgOutPlaneHist, "inclusive" : BkgInclusiveHist }}
    dataFull = { "background": {"in_plane" : BkgInPlaneHist, "mid_plane" : BkgMidPlaneHist, "out_of_plane" : BkgOutPlaneHist, "inclusive" : BkgInclusiveHist }, "signal" : {"in_plane" : SigInPlaneHist, "mid_plane" : SigMidPlaneHist, "out_of_plane" : SigOutPlaneHist, "inclusive" : SigInclusiveHist }}
 
    print("Done loading the histograms?",flush=True)

    print(BkgInPlaneHist)

    print("Fitting the background dominated region only")
    # draw_fit expects data to be in pachyderm histogram1D format

    # Estimate variables?
    # B is approximately 1/pi times the average value of the histogram in [-pi/2,pi/2]
    # FIXME implement this

    maxValue = 0.0
    for hist in ListOfHists:
      localMaxValue=hist.GetBinContent(hist.GetMaximumBin())
      if (localMaxValue > maxValue):
        maxValue = localMaxValue

    # magic number for safety
    maxValue = 1.05*maxValue
    # the near side has the addition of the 

    
    minv1=fCTask.GetGlobalV1Min()
    minv2t=fCTask.GetGlobalV2TMin()
    minv2a=fCTask.GetGlobalV2AMin()
    minv3=fCTask.GetGlobalV3Min()
    minv4t=fCTask.GetGlobalV4TMin()
    minv4a=fCTask.GetGlobalV4AMin()

    maxv1=fCTask.GetGlobalV1Max()
    maxv2t=fCTask.GetGlobalV2TMax()
    maxv2a=fCTask.GetGlobalV2AMax()
    maxv3=fCTask.GetGlobalV3Max()
    maxv4t=fCTask.GetGlobalV4TMax()
    maxv4a=fCTask.GetGlobalV4AMax()

    functorAverage = RPFFunctor.GetAverageValue()
    functorNumTriggers = RPFFunctor.GetNumTriggers()
    BRescale = 1
    if (functorAverage > -1 and functorNumTriggers > -1): BRescale = functorAverage * functorNumTriggers
    # BRescale might need some factor of pi/2 or pi
    BRescale *= math.pi

    print("Functor has average/trigger = %e and numTrigger = %e" % (functorAverage,functorNumTriggers))

    initB = maxValue / 1.05
    minB = 0
    maxB = maxValue
    fixedB = -1.0
    if (RPFFunctor.GetB_Min() > -1): minB = (BRescale) * RPFFunctor.GetB_Min()
    if (RPFFunctor.GetB_Max() > -1): maxB = (BRescale) * RPFFunctor.GetB_Max()
    if (RPFFunctor.GetInitB() > -1): initB= (BRescale) * RPFFunctor.GetInitB()
    if (RPFFunctor.GetFixedB() > -1): fixedB = (BRescale) * RPFFunctor.GetFixedB()

    # Now if the Functor has updated limits, replace the global limits
    # if they are unset, they are -1
    if (RPFFunctor.GetV1_Min() > -1): minv1=RPFFunctor.GetV1_Min()
    if (RPFFunctor.GetV1_Max() > -1): maxv1=RPFFunctor.GetV1_Max()

    if (RPFFunctor.GetV2T_Min() > -1): minv2t=RPFFunctor.GetV2T_Min()
    if (RPFFunctor.GetV2T_Max() > -1): maxv2t=RPFFunctor.GetV2T_Max()
    if (RPFFunctor.GetV2A_Min() > -1): minv2a=RPFFunctor.GetV2A_Min()
    if (RPFFunctor.GetV2A_Max() > -1): maxv2a=RPFFunctor.GetV2A_Max()

    if (RPFFunctor.GetV3_Min() > -1): minv3=RPFFunctor.GetV3_Min()
    if (RPFFunctor.GetV3_Max() > -1): maxv3=RPFFunctor.GetV3_Max()

    if (RPFFunctor.GetV4T_Min() > -1): minv4t=RPFFunctor.GetV4T_Min()
    if (RPFFunctor.GetV4T_Max() > -1): maxv4t=RPFFunctor.GetV4T_Max()
    if (RPFFunctor.GetV4A_Min() > -1): minv4a=RPFFunctor.GetV4A_Min()
    if (RPFFunctor.GetV4A_Max() > -1): maxv4a=RPFFunctor.GetV4A_Max()

    print("RPF Parameter Limits: ")
    print("     B %f %f" % (minB,maxB))
    print("    v1 %f %f" % (minv1,maxv1))
    print("   v2t %f %f" % (minv2t,maxv2t))
    print("   v2a %f %f" % (minv2a,maxv2a))
    print("    v3 %f %f" % (minv3,maxv3))
    print("   v4t %f %f" % (minv4t,maxv4t))
    print("   v4a %f %f" % (minv4a,maxv4a))


    MyDefaultArgs={"limit_B": [minB,maxB],"limit_v1":[minv1,maxv1],"limit_v2_t": [minv2t, maxv2t],"limit_v2_a": [minv2a, maxv2a],"limit_v3": [minv3, maxv3],"limit_v4_t": [minv4t, maxv4t],"limit_v4_a": [minv4a, maxv4a]}
    #MyDefaultArgs={"limit_B": [0.,maxValue],"limit_v1":[minv1,maxv1],"limit_v2_t": [minv2t, maxv2t],"limit_v2_a": [minv2a, maxv2a],"limit_v3": [minv3, maxv3],"limit_v4_t": [minv4t, maxv4t],"limit_v4_a": [minv4a, maxv4a]}
#    MyDefaultArgs={"limit_B": [0.,maxValue],"limit_v2_t": [0.0, maxv2t],"limit_v2_a": [0.0, maxv2a],"limit_v3": [0.0, maxv3],"limit_v4_t": [0.0, maxv4t],"limit_v4_a": [0.0, maxv4a]}
    #MyDefaultArgs={"limit_B": [0.,maxValue],"limit_v2_a": [0.0, 0.5] }
#    MyDefaultArgs={"limit_B": [0.,1e7],"limit_v2_a": [0.0, 0.5],"fix_v3":True,"v3":0.0}



    # Getting the initial values.
    #FlowV2TValue = FlowV2TGraph.GetY()[] # Trigger pt bin
    #FlowV4TValue = FlowV4TGraph.GetY()[] # Trigger pt bin

    #FlowV3Value = # good luck with this one

    #FlowV2AValue = FlowV2AGraph.GetY()[iObsBin]
    #FlowV4AValue = FlowV4AGraph.GetY()[iObsBin]
    # FIXME could use spline interpolation
    # Need
      #fPtAMin = fTrackPtProjectionSE->GetXaxis()->GetBinLowEdge(iObsBin+1);
      #fPtAMax = fTrackPtProjectionSE->GetXaxis()->GetBinUpEdge(iObsBin+1);
    # or write a function for the c++ task that returns the value
    # may also want one for the error

#    FlowV2AValue = FlowV2AGraph.Eval(pTA_Value)
#    FlowV4AValue = FlowV4AGraph.Eval(pTA_Value)

    FlowV2AError = 0
    FlowV2AValue=fCTask.GetFlowVNAFromObsBin(2,iObsBin)
    FlowV2AError=fCTask.GetFlowVNAeFromObsBin(2,iObsBin)

    print("Found Flow V2 = %f +- %f" % (FlowV2AValue,FlowV2AError));

    FlowV4AError = 0
    FlowV4AValue=fCTask.GetFlowVNAFromObsBin(4,iObsBin)
    FlowV4AError=fCTask.GetFlowVNAeFromObsBin(4,iObsBin)
    print("Found Flow V4 = %f +- %f" % (FlowV4AValue,FlowV4AError));

    # Now apply fFlowTermMode
    # 0 = no action
    # 1 = Fix V2A, V4A to interpolation of flow grpahs
    # 2 = Limit V2A, V4A to flow graph interpolation +- sigma
    
    print("   FTM=%d" % fFlowTermMode)
    # This part is probably mostly overwritten by info from functors
    if (fFlowTermMode == 0):
      # Setting initial values
      MyDefaultArgs['v2_a']=FlowV2AValue
      MyDefaultArgs['v4_a']=FlowV4AValue
    if (fFlowTermMode == 1):
      MyDefaultArgs["fix_v2_a"]=True 
      MyDefaultArgs['v2_a']=FlowV2AValue
      #MyDefaultArgs["fix_v4_a"]=True 
      MyDefaultArgs['v4_a']=FlowV4AValue
    if (fFlowTermMode == 2):
      MyDefaultArgs['v2_a']=FlowV2AValue
      MyDefaultArgs['v4_a']=FlowV4AValue
      MyDefaultArgs['limit_v2_a']=[FlowV2AValue-FlowV2AError,FlowV2AValue+FlowV2AError]
      #MyDefaultArgs['limit_v4_a']=[FlowV4AValue-FlowV4AError,FlowV4AValue+FlowV4AError]
       
    # Update this guess
    MyDefaultArgs['v3']=0.01




    if (RPFFunctor.GetInitV1() > -1):
      MyDefaultArgs['v1'] = RPFFunctor.GetInitV1()
    if (RPFFunctor.GetInitV2T() > -1):
      MyDefaultArgs['v2_t'] = RPFFunctor.GetInitV2T()
    if (RPFFunctor.GetInitV2A() > -1):
      MyDefaultArgs['v2_a'] = RPFFunctor.GetInitV2A()
    if (RPFFunctor.GetInitV3() > -1):
      MyDefaultArgs['v3'] = RPFFunctor.GetInitV3()
    if (RPFFunctor.GetInitV4T() > -1):
      MyDefaultArgs['v4_t'] = RPFFunctor.GetInitV4T()
    if (RPFFunctor.GetInitV4A() > -1):
      MyDefaultArgs['v4_a'] = RPFFunctor.GetInitV4A()

    # Now if the Functor has fixed values, use those to replace the initial values

    if (fixedB > -1):
      MyDefaultArgs['B'] = fixedB
      MyDefaultArgs["fix_B"]=True 


    if (RPFFunctor.GetFixedV1() > -1):
      MyDefaultArgs['v1'] = RPFFunctor.GetFixedV1()
      MyDefaultArgs["fix_v1"]=True 
    else:
      MyDefaultArgs["fix_v1"]=False

    if (RPFFunctor.GetFixedV2T() > -1):
      MyDefaultArgs['v2_t'] = RPFFunctor.GetFixedV2T()
      MyDefaultArgs["fix_v2_t"]=True 
    else:
      MyDefaultArgs["fix_v2_t"]=False
    if (RPFFunctor.GetFixedV2A() > -1):
      MyDefaultArgs['v2_a'] = RPFFunctor.GetFixedV2A()
      MyDefaultArgs["fix_v2_a"]=True 
    else:
      MyDefaultArgs["fix_v2_a"]=False

    if (RPFFunctor.GetFixedV3() > -1):
      MyDefaultArgs['v3'] = RPFFunctor.GetFixedV3()
      MyDefaultArgs["fix_v3"]=True 
    else:
      MyDefaultArgs["fix_v3"]=False

    if (RPFFunctor.GetFixedV4T() > -1):
      MyDefaultArgs['v4_t'] = RPFFunctor.GetFixedV4T()
      MyDefaultArgs["fix_v4_t"]=True 
    else:
      MyDefaultArgs["fix_v4_t"]=False
    if (RPFFunctor.GetFixedV4A() > -1):
      MyDefaultArgs['v4_a'] = RPFFunctor.GetFixedV4A()
      MyDefaultArgs["fix_v4_a"]=True 
    else:
      MyDefaultArgs["fix_v4_a"]=False



    # Switch on V1,V5,V6T,V6A terms if requested
    # iV1Mode iV5Mode iV6TMode iV6AMode
    # v5,v6 not implemented here

    #MyDefaultArgs["v1"]=0.0
    #if (iV1Mode==1):
    #  MyDefaultArgs["fix_v1"]=False
    ##  MyDefaultArgs["v1"]=0.0

    #if (iV3Mode==1):
    #  MyDefaultArgs["fix_v3"]=True
    ##if (iV3Mode==2):

    if (iObsBin >= 4):
      print("doing the fix thing")
      for data in dataFull:
        print("data = %s" % (data))
        # background or signal ...
        for entry in dataFull[data]:
          print(entry)
          FixErrorsOnHist(dataFull[data][entry])

    if (iObsBin >= nFixV40Last):
      MyDefaultArgs["fix_v4_t"]=True 
      MyDefaultArgs['v4_t'] = 0.0
      #if (fFlowTermMode != 1): # Checking if already fixed to a value
      #  MyDefaultArgs["fix_v4_a"]=True 
      #  MyDefaultArgs['v4_a'] = 0.0
      # could also fix off the v6,v5

    MyUserArgs={}
    #MyUserArgs={"v4_t": 0,"fix_v4_t": True,"v4_a": 0,"fix_v4_a": True}
    InclusiveUserArgs={}
    RPDepUserArgs={}
    # Get initial parameters from RP Dep Fit
    ReduxUserArgs={}

    MyUserArgs.update(MyDefaultArgs)
    InclusiveUserArgs.update(MyDefaultArgs)
    RPDepUserArgs.update(MyDefaultArgs)
    ReduxUserArgs.update(MyDefaultArgs)
  

    print("Fitting ObsBin %d with user args:" % (iObsBin))
    print(MyUserArgs)

    # Create a covariance matrix in the CTask before running the fit
    fCTask.CreateCovarianceMatrix(1,iObsBin,fCTask.GetNTotalParMethod(1))
    #fCTask.CreateCovarianceMatrix(1,iObsBin,len(Py_TGraphs))

    # The Fitting is done here
#    (success,data_BkgFit,_) = rp_fit.fit(data=dataBkg, user_arguments=MyUserArgs)
 
    try:
      (success,data_BkgFit,_) = rp_fit.fit(data=dataFull, user_arguments=MyUserArgs)
    except FitFailed:
#    except pachyderm.fit.base.FitFailed:
      print("Caught a Fit failed exception. Continuing")
      fCTask.UpdateNumFreePar(1,iObsBin,1)
      continue
    except RuntimeError:
      print("Caught a run-time error. Continuing")
      fCTask.UpdateNumFreePar(1,iObsBin,1)
      continue
    print("Finished doing the fit, maybe")

    print("Fit result: {fit_result}".format(fit_result = rp_fit.fit_result))

    BkgFitResults = rp_fit.fit_result

    BkgChiSquare = BkgFitResults.minimum_val
    BkgNDOF = BkgFitResults.nDOF

    # Get the covariance matrix

    #covariance_matrix: Dict[Tuple[str, str], float]
    CovMatrix = BkgFitResults.covariance_matrix

    Py_ChiSq_TGraph.SetPoint(iObsBin,ObsBinCenter,BkgChiSquare/BkgNDOF)
    Py_ChiSq_TGraph.SetPointError(iObsBin,ObsBinWidth,0)

    NumFreeParams=0

    # Setting the Event Plane bin parameter to be fixed
    # This parameter only exists in the C++ implementation
    fCTask.UpdateFreeMaskArray(1,iObsBin,0,False);
    for j in range(len(Py_TGraphs)):
      Py_Val = BkgFitResults.values_at_minimum[PyParamNames[j]]
      Py_Err = BkgFitResults.errors_on_parameters[PyParamNames[j]]
      Py_TGraphs[j].SetPoint(iObsBin,ObsBinCenter,Py_Val)
      Py_TGraphs[j].SetPointError(iObsBin,ObsBinWidth,Py_Err)

      # Sending CTask information about the free parameters

      # Sending CTask entries in the covariance matrix for all parameters
      #for k in range(len(Py_TGraphs)):
      # take the parameter name, check if it has an entry in the covariance matrix
      # [parname,parname]?
      # if so, iterate over other parameters
      # don't try using k >= j to avoid double entries, as the TMatrixDSym only enforces
      # asymmetry when reading from disk, where the upper right triangle is the only part 
      # saved to disk, while the full matrix exists in memory
      if (PyParamNames[j],PyParamNames[j]) in CovMatrix: 
        # This parameter was a free one
        print("Updating fCTask with free parameter %s = %e #pm %e" % (PyParamNames[j],Py_Val,Py_Err))
        fCTask.UpdateFreeMaskArray(1,iObsBin,j,True);
        fCTask.UpdateParNamesArray(1,iObsBin,NumFreeParams,PyParamNames[j])
        fCTask.UpdateParMuArray(1,iObsBin,NumFreeParams,Py_Val)
        fCTask.UpdateParSigmaArray(1,iObsBin,NumFreeParams,Py_Err)
        print("Finished updating the diagonal")
        NumFreeParams+=1


        for k in range(len(Py_TGraphs)):
          if (PyParamNames[k],PyParamNames[k]) in CovMatrix:
            fCovValue = CovMatrix[(PyParamNames[j],PyParamNames[k])]
            fCTask.UpdateCovarianceMatrix(1,iObsBin,j+1,k+1,fCovValue)
        print("Finished updating the off-diagonal")
      else:
        # This was not a free parameter
        fCTask.UpdateFreeMaskArray(1,iObsBin,j,False);


      # Storing result for next fit
      InclusiveUserArgs[PyParamNames[j]]=Py_Val
      RPDepUserArgs[PyParamNames[j]]=Py_Val
      # Also Set some reasonable starting values for parameters
      RPDepUserArgs["in_plane_ns_sigma"]=0.36
      RPDepUserArgs["in_plane_as_sigma"]=0.42
      RPDepUserArgs["mid_plane_ns_sigma"]=0.36
      RPDepUserArgs["mid_plane_as_sigma"]=0.42
      RPDepUserArgs["out_of_plane_ns_sigma"]=0.36
      RPDepUserArgs["out_of_plane_as_sigma"]=0.42


    fCTask.UpdateNumFreePar(1,iObsBin,NumFreeParams)

#------------------------------------------------------------------------------------------------------------------
#| 0 | in_plane_ns_amplitude        |  0.52E4   |  0.03E4   |            |            |    0    |  1e+07  |       |
#| 1 | in_plane_as_amplitude        |    530    |    310    |            |            |    0    |  1e+07  |       |
#| 2 | in_plane_ns_sigma            |   0.363   |   0.023   |            |            |  0.02   |   0.7   |       |
#| 3 | in_plane_as_sigma            |   0.42    |   0.17    |            |            |  0.02   |   0.7   |       |
#| 4 | in_plane_signal_pedestal     |    0.0    |    1.0    |            |            |         |         |  yes  |
#| 5 | B                            |  1.891E5  |  0.002E5  |            |            |    0    |  1e+07  |       |
#| 6 | v2_t                         | 0.547E-1  | 0.007E-1  |            |            |  0.001  |   0.2   |       |
#| 7 | v2_a                         | 3.000E-2  | 0.008E-2  |            |            |  0.03   |   0.5   |       |
#| 8 | v4_t                         | 0.005E-4  | 0.714E-4  |            |            |    0    |   0.5   |       |
#| 9 | v4_a                         | 0.005E-4  | 2.563E-4  |            |            |    0    |   0.5   |       |
#| 10| v1                           |  0.000E1  |  0.000E1  |            |            |   -1    |    1    |  yes  |
#| 11| v3                           |  2.5E-3   |  0.6E-3   |            |            |   -1    |    1    |       |
#| 12| mid_plane_ns_amplitude       |  0.88E4   |  0.03E4   |            |            |    0    |  1e+07  |       |
#| 13| mid_plane_as_amplitude       |  0.49E4   |  0.03E4   |            |            |    0    |  1e+07  |       |
#| 14| mid_plane_ns_sigma           |   0.467   |   0.022   |            |            |  0.02   |   0.7   |       |
#| 15| mid_plane_as_sigma           |   0.700   |   0.027   |            |            |  0.02   |   0.7   |       |
#| 16| mid_plane_signal_pedestal    |    0.0    |    1.0    |            |            |         |         |  yes  |
#| 17| out_of_plane_ns_amplitude    |  0.76E4   |  0.03E4   |            |            |    0    |  1e+07  |       |
#| 18| out_of_plane_as_amplitude    |  0.361E4  |  0.031E4  |            |            |    0    |  1e+07  |       |
#| 19| out_of_plane_ns_sigma        |   0.436   |   0.019   |            |            |  0.02   |   0.7   |       |
#| 20| out_of_plane_as_sigma        |   0.45    |   0.04    |            |            |  0.02   |   0.7   |       |
#| 21| out_of_plane_signal_pedestal |    0.0    |    1.0    |            |            |         |         |  yes  |
#------------------------------------------------------------------------------------------------------------------




#    Py_B = BkgFitResults.values_at_minimum["B"]
#    Py_B_Err = BkgFitResults.errors_on_parameters["B"]
#    Py_B_TGraph.SetPoint(iObsBin,ObsBinCenter,Py_B)
#    Py_B_TGraph.SetPointError(iObsBin,ObsBinWidth,Py_B_Err)

#    Py_V2T = BkgFitResults.values_at_minimum["v2_t"]
#    Py_V2T_Err = BkgFitResults.errors_on_parameters["v2_t"]
#    Py_V2T_TGraph.SetPoint(iObsBin,ObsBinCenter,Py_V2T)
#    Py_V2T_TGraph.SetPointError(iObsBin,ObsBinWidth,Py_V2T_Err)


    fit_label="%s \n %s \n %s" % (fCTask.GetLabel(),fCTask.GetLabel2(),LocalLabel)
    fit_label=LocalLabel
    filename="%s/PyRPF_BkgFit_ObsBin%d.pdf" % (fOutputDir,iObsBin)
    plot.draw_fit(rp_fit=rp_fit,data=data_BkgFit,fit_label=fit_label,filename=filename)
    filename="%s/PyRPF_BkgFit_ObsBin%d.png" % (fOutputDir,iObsBin)
    plot.draw_fit(rp_fit=rp_fit,data=data_BkgFit,fit_label=fit_label,filename=filename)
#    filename="%s/PyRPF_BkgFit_PlotAll_ObsBin%d.pdf" % (fOutputDir,iObsBin)
#    plot.draw_fit(rp_fit=rp_fit,data=dataFull,fit_label=fit_label,filename=filename)
#    plot.fit_draw_func(data=data_BkgFit,fit_label=fit_label,filename=filename)


    if (enableInclusiveFit):
      # Settings for inclusive fit
      InclusiveUserArgs["fix_v2_t"]=True
      InclusiveUserArgs["fix_v2_a"]=True
      # could estimate yields based on integrals and the B parameter found earlier
      print(str(InclusiveUserArgs))
      print("Fitting background dominated and signal regions")

      #(success_IncSig,data_IncSig,_) = rp_fit_IncSig.fit(data=dataFull,user_arguments=InclusiveUserArgs)
      try:
        (success_IncSig,data_IncSig,_) = rp_fit_IncSig.fit(data=dataFull,user_arguments=InclusiveUserArgs)
      except RuntimeError:
        print("Caught a run-time error. Continuing")
        continue
      except FitFailed:
        print("Caught a Fit failed exception. Continuing")
        continue


#      try:
#        (success_IncSig,data_IncSig,_) = rp_fit_IncSig.fit(data=dataFull,user_arguments=InclusiveUserArgs)
#      except pachyderm.fit.base.FitFailed:
#        print("Da Fit Failed")

      print("Fit result: {fit_result}".format(fit_result = rp_fit_IncSig.fit_result))
      filename="%s/PyRPF_IncFit_ObsBin%d.pdf" % (fOutputDir,iObsBin)
      plot.draw_fit(rp_fit=rp_fit_IncSig,data=data_IncSig,fit_label=fit_label,filename=filename)
      filename="%s/PyRPF_IncFit_ObsBin%d.png" % (fOutputDir,iObsBin)
      plot.draw_fit(rp_fit=rp_fit_IncSig,data=data_IncSig,fit_label=fit_label,filename=filename)

    if (enableRPDepFit):
      print("Fitting background dominated and signal regions with RP dependent signal")
#      (success_RPSig,data_RPSig,_) = rp_fit_RPSig.fit(data=dataFull,user_arguments=RPDepUserArgs)
      try:
        (success_RPSig,data_RPSig,_) = rp_fit_RPSig.fit(data=dataFull,user_arguments=RPDepUserArgs)
      except RuntimeError:
        print("Caught a run-time error. Continuing")
        continue
      except FitFailed:
        print("Caught a Fit failed exception. Continuing")
        continue
      filename="%s/PyRPF_RPDepF_ObsBin%d.pdf" % (fOutputDir,iObsBin)
      plot.draw_fit(rp_fit=rp_fit_RPSig,data=data_RPSig,fit_label=fit_label,filename=filename)
      filename="%s/PyRPF_RPDepF_ObsBin%d.png" % (fOutputDir,iObsBin)
      plot.draw_fit(rp_fit=rp_fit_RPSig,data=data_RPSig,fit_label=fit_label,filename=filename)

      RPDepBkgFitResults=rp_fit_RPSig.fit_result

      RPDepChiSquare = RPDepBkgFitResults.minimum_val
      RPDepNDOF = RPDepBkgFitResults.nDOF

      Py_RPDep_ChiSq_TGraph.SetPoint(iObsBin,ObsBinCenter,RPDepChiSquare/RPDepNDOF)
      Py_RPDep_ChiSq_TGraph.SetPointError(iObsBin,ObsBinWidth,0)


      # Save fit parameters
      if (enableRPDepFit):
        for j in range(len(Py_RPDep_TGraphs)):
          Py_Val = RPDepBkgFitResults.values_at_minimum[PyParamNames[j]]
          Py_Err = RPDepBkgFitResults.errors_on_parameters[PyParamNames[j]]
          Py_RPDep_TGraphs[j].SetPoint(iObsBin,ObsBinCenter,Py_Val)
          Py_RPDep_TGraphs[j].SetPointError(iObsBin,ObsBinWidth,Py_Err)




      if (enableReduxFit):
        for key in rp_fit_RPSig.fit_result.values_at_minimum:
          if key in BkgFitResults.parameters:
            print("Loading parameter %s" % (str(key)))
            ReduxUserArgs[key]=rp_fit_RPSig.fit_result.values_at_minimum[key]
        print("Fitting Background dominated only with initial parameters from RP Dep Fit (Redux)")

        #(success_Redux,data_Redux,_) = rp_fit_Redux.fit(data=dataFull,user_arguments=ReduxUserArgs)
        try:
          (success_Redux,data_Redux,_) = rp_fit_Redux.fit(data=dataFull,user_arguments=ReduxUserArgs)
        except RuntimeError:
          print("Caught a run-time error, continuing")
          continue


        filename="%s/PyRPF_Redux_ObsBin%d.pdf" % (fOutputDir,iObsBin)
        plot.draw_fit(rp_fit=rp_fit_Redux,data=data_Redux,fit_label=fit_label,filename=filename)


  # End of obs bin loop

  print("Finished the Observable Bin Loop")

  # Get Vn from C++ Code
  C_B_TGraph = fCTask.GetParamGraph(0)
  C_V1_TGraph = fCTask.GetParamGraph(1)
  C_V2T_TGraph = fCTask.GetParamGraph(2)
  C_V2A_TGraph = fCTask.GetParamGraph(3)
  C_V3_TGraph = fCTask.GetParamGraph(4)
  C_V4T_TGraph = fCTask.GetParamGraph(5)
  C_V4A_TGraph = fCTask.GetParamGraph(6)

  #C_TGraphs=[C_B_TGraph, C_V2T_TGraph, C_V2A_TGraph,  C_V3_TGraph, C_V4T_TGraph, C_V4A_TGraph]
  C_TGraphs=[C_B_TGraph, C_V1_TGraph, C_V2T_TGraph, C_V2A_TGraph,  C_V3_TGraph, C_V4T_TGraph, C_V4A_TGraph]

  for graph in C_TGraphs:
    graph.SetLineColor(CColor)
    graph.SetMarkerColor(CColor)
    graph.SetMarkerStyle(CMarkerStyle)
  for graph in Py_TGraphs:
    for i in range(1+nObsBins-nSkipLast):
      graph.RemovePoint(nSkipLast+1)
    graph.SetLineColor(PyColor)
    graph.SetMarkerColor(PyColor)
    graph.SetMarkerStyle(PyMarkerStyle)
#    graph.GetXaxis().SetTitle("z_{T}")
  if (enableRPDepFit):
    for graph in Py_RPDep_TGraphs:
      for i in range(1+nObsBins-nSkipLast):
        graph.RemovePoint(nSkipLast+1)
      graph.SetLineColor(PyRPDepColor)
      graph.SetMarkerColor(PyRPDepColor)
      graph.SetMarkerStyle(PyRPDepMarkerStyle)
  #    graph.GetXaxis().SetTitle("z_{T}")
  

  c1 = TCanvas("c1","c1",900,600)  
  c1.cd()
  
  # FIXME note that B value must be scaled by num_triggers to be compared

  # Comparing CTask and Python Bkg Parameters
  MultiGraphs=[]
  MergedList=tuple(zip(C_TGraphs,Py_TGraphs))
  for i in range(len(MergedList)):
    print("i = %d" % (i))
    (CGraph,PyGraph)=MergedList[i]
    c1.Clear()
    tmg=TMultiGraph()

  #  tmg.Add(CGraph,"lp")
  #  tmg.Add(PyGraph,"lp")
    tmg.Add(CGraph.Clone(),"lp")
    tmg.Add(PyGraph.Clone(),"lp")

 #   gROOT.SetOwnership(CGraph,False)
 #   gROOT.SetOwnership(PyGraph,False)

    tmg.Draw("a")
    tmg.SetName("tmg_%d" % (i))
    tmg.SetTitle(CGraph.GetTitle())
    tmg.GetXaxis().SetTitle(CGraph.GetXaxis().GetTitle())
    PyGraph.SetTitle(CGraph.GetTitle())
    PyGraph.GetXaxis().SetTitle(CGraph.GetXaxis().GetTitle())
    MultiGraphs.append(tmg)

    filename="RPF_Comp_Param_%s" % (ParamNames[i])
    c1.Print("%s/%s.pdf" % (fOutputDir,filename))
    c1.Print("%s/CFiles/%s.C" % (fOutputDir,filename))

  # Comparing Python BkgOnly and InclusiveSignal
  PyMultiGraphs=[]
  MergedList=tuple(zip(Py_TGraphs,Py_RPDep_TGraphs))
  for i in range(len(MergedList)):
    print("i = %d" % (i))
    (PyBkgGraph,PySigGraph)=MergedList[i]
    c1.Clear()
    tmg=TMultiGraph()

    #tmg.Add(PyBkgGraph,"lp")
    #tmg.Add(PySigGraph,"lp")
    tmg.Add(PyBkgGraph.Clone(),"lp")
    tmg.Add(PySigGraph.Clone(),"lp")

    tmg.Draw("a")
    tmg.SetName("tmg_%d" % (i))
    tmg.SetTitle(PyBkgGraph.GetTitle())
    tmg.GetXaxis().SetTitle(PyBkgGraph.GetXaxis().GetTitle())

    PySigGraph.SetTitle(PyBkgGraph.GetTitle())
    PySigGraph.GetXaxis().SetTitle(PyBkgGraph.GetXaxis().GetTitle())

    PyMultiGraphs.append(tmg)

    filename="RPF_CompBkgSig_Param_%s" % (ParamNames[i])
    c1.Print("%s/%s.pdf" % (fOutputDir,filename))
    c1.Print("%s/CFiles/%s.C" % (fOutputDir,filename))


  # Drawing the ChiSquare Graphs

#  print("Saving to file %s" % (fOutputFile))
#  OutFile = TFile(fOutputFile,"UPDATE")
#  print("Opened file %s" % (OutFile.GetName()))

  print("Trying to get the output file from the c++ task")
  OutFile = fCTask.GetOutputFile()



#  for graph in Py_TGraphs:
#    OutFile.Add(graph)
 #   graph.Write()
#  if (enableRPDepFit):
#    for graph in Py_RPDep_TGraphs:
#      OutFile.Add(graph)
  #    graph.Write()

  # This code would directly add matrices consisting of only the free parameters
  #print("Adding Covariance Matrices")
  #for CovMatrix in CovMatrices:
    # 1 is for method2
   # fCTask.InputPyCovMatrix(1,CovMatrix)
    # Can also add vector 2, 3, etc for Incl, RPDep



  # Add the chisq/ndf and  parameter graphs to the CTask
  fCTask.InputPyBkgChiSqGraph(Py_ChiSq_TGraph)
  print("About to add %d PyBkg Graphs to CTask" % (len(Py_TGraphs)))
  for i in range(len(Py_TGraphs)):
    print("Adding graph [ %d ] = %s" % (i,Py_TGraphs[i].GetName()))
    fCTask.InputPyBkgParamGraph(i,Py_TGraphs[i])   

  if (enableRPDepFit):
    print("About to add %d PyRODepBkg Graphs to CTask" % (len(Py_RPDep_TGraphs)))
    fCTask.InputPyRPSChiSqGraph(Py_RPDep_ChiSq_TGraph)
    for i in range(len(Py_RPDep_TGraphs)):
      fCTask.InputPyRPSParamGraph(i,Py_RPDep_TGraphs[i])   
 

  print("about to try deleting recursive objects")

  del MultiGraphs
  del PyMultiGraphs

  # Input the Covariance matrices





  # Should maybe find a nice way to save the covariance matrices
  # TH2D ?
  print("Writing File...")
 # OutFile.Write()
  print("Successfully wrote file! (maybe)")
 # OutFile.Close()
  print("Closed file!")

  print("=======================================================")
  print("Done with the python part")
  print("=======================================================")

  sys.stdout.flush()
 
  fCTask.Run_Part2()

#  for graph in Py_TGraphs:
#    del graph
#  for graph in Py_RPDep_TGraphs:
#    del graph

#  C_V2T_TGraph.Draw()
#  Py_V2T_TGraph.Draw("SAME PL")
#  c1.Print("%s/TestComparison.pdf" % (fOutputDir))
#  c1.Print("%s/CFiles/TestComparison.C" % (fOutputDir))
  




#  fInputFileEvt=[fInputFile_EP_0,fInputFile_EP_1,fInputFile_EP_2]

#  VariableInfo=fInputFileAll.Get("VariableInfo")
#  fObservable=VariableInfo.GetBinContent(1)
#  print("This analysis is in observable %d" % (fObservable))
  


#import sys
## sys.argv.append( '-b' )
#from reaction_plane_fit import three_orientations
#import os
#import ROOT
#import yaml

#from ROOT import gROOT
##gROOT.LoadMacro("blah blah blah")
#
#  print("sys.argv = %d" % len(sys.argv))
#  if (len(sys.argv) > 2):
#    YamlFile=sys.argv[1]
#  else:
#    usage="Usage:  %s [yaml config]" % sys.argv[0]
#    print(usage)
#    exit(0)

#  print("o Load config from yamlFile = %s" % YamlFile)

# Example User Arguments
"""
limit_B: [0, 1e7]
v2_t: 0.08
limit_v2_t: [0.001, 0.15]
v2_a: 0.15
limit_v2_a: [0.03, 0.5]
v3: 0.00015
error_v3: 1e-5
v4_t: 0.03
limit_v4_t: [0, 0.1]
v4_a: 0.05
limit_v4_a: [0, 0.2]


limit_in_plane_ns_sigma: [0.15, 0.7]
limit_out_of_plane_ns_sigma: [0.15, 0.7]
limit_mid_plane_ns_sigma: [0.15, 0.7]

in_plane_ns_sigma: 0.5
limit_in_plane_ns_sigma: [0.15, 0.7]
in_plane_as_sigma: 0.5
limit_in_plane_as_sigma: [0.01, 1.0]
# Mid plane
mid_plane_ns_sigma: 0.5
limit_mid_plane_ns_sigma: [0.15, 0.7]
mid_plane_as_sigma: 0.5
limit_mid_plane_as_sigma: [0.01, 1.0]
# Out of plane
out_of_plane_ns_sigma: 0.5
limit_out_of_plane_ns_sigma: [0.15, 0.7]
out_of_plane_as_sigma: 0.5
limit_out_of_plane_as_sigma: [0.01, 1.0]

v4_t: 0
fix_v4_t: true
v4_a: 0
fix_v4_a: true
"""



