#! /usr/bin/env python3

import sys
sys.argv.append( '-b' ) # batch mode
import os
import ROOT
import yaml
import math

from ROOT import gROOT
# gROOT.LoadMacro("asdf.cxx")

from ROOT import TProfile2D,TProfile
import array

# Sum up the bins of the y-axis, return array of (val,err)
def SumUpProfile(Pf2,CentBin):
  valList = []
  nBinsX = Pf2.GetXaxis().GetNbins()
  nBinsY = Pf2.GetYaxis().GetNbins()
  YProj = Pf2.ProfileY(("%s_pfy_Cent%d" % (Pf2.GetName(),CentBin)))
  for j in range(nBinsY):
#    print(j)
    localVal = YProj.GetBinContent(j+1)
    localErr = YProj.GetBinError(j+1)
#    print("bin %d has total %f \\pm %f" % (j,localVal,localErr))
    valList.append((localVal,localErr))

  
  return valList
#    for i in range (nBinsX):
#      localSum += Pf2.GetBinContent(i,j)
#      localErr 
#  delete YProj


def CalcEventPlaneResolution():
  print("---------------------------------")
  print("| Starting Event Plane Res Calc.|")
  print("---------------------------------")

  InputFileName="AnalysisResults.root"

  RootDirCent0Name="AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_Cent0_histos"
  RootDirCent1Name="AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_Cent1_histos"
  RootDirCent2Name="AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_Cent2_histos"
  RootDirCent3Name="AliAnalysisTaskMBPi0Candv_Pi0H_SE_tracks_caloClusters_Cent3_histos"

  InputFile = ROOT.TFile(InputFileName)
  if ( InputFile == 0 ):
    print("Could not open file")  


  dirlist = InputFile.GetListOfKeys()
  iter = dirlist.MakeIterator()
  key = iter.Next()
  dirs = {}
  td = None
  while key:
    if key.GetClassName() == 'AliEmcalList':  #'TDirectory'
      td = key.ReadObj()
      dirName = td.GetName()
      print("found directory", dirName)
      dirs[dirName] = td
    key = iter.Next()

  for dir in dirs:
    print(dir)
    CentIndex=4+dir.find('Cent')
    CentBin = int(dir[CentIndex:CentIndex+1])
    print("Cent Bin %d" % CentBin)
    localDir = InputFile.Get(dir)
#    Vals = []
    for i in range(6): # individual N values
      LocalVals = []
      for j in range(3): # individual Dn values
        Pf2Name="EPR_CosD%d_N%d" % (j+1,i+1)
        Pf2=localDir.FindObject(Pf2Name)
        ValArray = SumUpProfile(Pf2,CentBin) # returns a list of tuples, one per cent bin
        LocalVals.append(ValArray)
#      print("Cent = %d" % CentBin)
   #   print(LocalVals)
      LocalRn=0
      LocalRn_Un=0
      MeanCosD1=LocalVals[0][CentBin][0]
      MeanCosD2=LocalVals[1][CentBin][0]
      MeanCosD3=LocalVals[2][CentBin][0]
      MeanCosD1_Un=LocalVals[0][CentBin][1]
      MeanCosD2_Un=LocalVals[1][CentBin][1]
      MeanCosD3_Un=LocalVals[2][CentBin][1]
      if (LocalVals[2][CentBin][0] > 0.):
        LocalRn=math.sqrt((MeanCosD1 * MeanCosD2) / MeanCosD3)
        #LocalRn=math.sqrt((LocalVals[0][CentBin][0] * LocalVals[1][CentBin][0]) / LocalVals[2][CentBin][0])
        LocalRn_Un = LocalRn * math.sqrt((0.5)*math.pow(MeanCosD1_Un/MeanCosD1,2) + (0.5)*math.pow(MeanCosD2_Un/MeanCosD2,2) + (0.5)*math.pow(MeanCosD3_Un/MeanCosD3,2))
        #(LocalVals[0][CentBin][0]/LocalVals[0][CentBin][1]) + LocalVals[1][CentBin][0]) / LocalVals[2][CentBin][0])
      print("Found R_{%d} = %f   \\pm %f" % (i+1,LocalRn,LocalRn_Un))




#  print("Seraching for object "+ RootDirCent0Name + ".EPR_CosD1_N1")
#  Pf2_CosD1_N1 = InputFile.Get(RootDirCent0Name + ".EPR_CosD1_N1")
#  if (Pf2_CosD1_N1 == 0):
#    print("Could not get object "+ RootDirCent0Name + ".EPR_CosD1_N1")
#    exit()
#
#  print("Found " + Pf2_CosD1_N1.GetName() + " succesfully")



if __name__ == '__main__':
  CalcEventPlaneResolution()
