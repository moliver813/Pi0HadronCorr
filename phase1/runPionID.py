#! /usr/bin/env python3


##  !/usr/bin/python3

#   ##! /usr/bin/env python
# -----------------
# Code to analyze cluster pair mass vs pt distributions, and produce information on pi0 mass peak and width, yield, and background. 
# Should be able to look at individual EP, Centrality bins
# Goal is to be able to act on output of Pi0Cand mode or Pi0HCorr Mode

# ---------------------
#  Running the macro
# ---------------------
# python GHphase3.py
# python GHphase3.py [yaml-file]

# ---------------------
# import
# ---------------------

import sys
sys.argv.append( '-b' ) # batch mode
import os
import ROOT
import yaml

from ROOT import gROOT
gROOT.LoadMacro("PionID.cxx")

def runPionID():

  print("...................................")

  LocalPath=""
  YamlFile=""

  # Check if a yaml config file has been passed in args:
  print("sys.argv = %d" % len(sys.argv))
  if (len(sys.argv) > 2):
    YamlFile=sys.argv[1]
  else:
    usage="Usage:  %s [yaml config]" % sys.argv[0]
    print(usage)
    exit(0)

  print("o Load config from yamlFile = %s" % YamlFile)

  LocalPath=""
  LocalPath+=YamlFile
  print("o Load yaml file: %s" % LocalPath)

  configurations = open(LocalPath)
  configurations = yaml.load(configurations,Loader=yaml.FullLoader)

  if not configurations: print(" did not get the yaml file!")

  FilePath           = configurations['file_path']
  InputFileName_Pi0  = configurations['file_name_pi0']
  if 'output_dir' in  configurations:
    OutputDir      = configurations['output_dir']
  else:
    OutputDir      = "output"
  if 'output_file' in configurations:
    OutputFileName = configurations['output_file']
  else:
    OutputFileName = "SecondAnalysis.root"

  InputFile_Pi0  = ROOT.TFile(FilePath+"/"+InputFileName_Pi0)

  task = ROOT.PionID()
  if "debug_level" in configurations:
    task.SetDebugLevel(configurations['debug_level'])
  else:
    task.SetDebugLevel(2)
  task.SetPi0CandInputFile(InputFile_Pi0)

  if "enable_performance" in configurations:
    task.EnablePerformance(configurations['enable_performance'])

  # Setting the input EmcalList
  if "list_name" in configurations:
    task.SetListName(configurations["list_name"])

  if "output_dir" in configurations:
    OutputDir = configurations["output_dir"]
    if not os.path.exists(OutputDir):
      print("Making directory for output")
      os.makedirs(OutputDir)
#      os.mkdir(OutputDir)
    if not os.path.exists(OutputDir+"/CFiles"):
      os.makedirs(OutputDir+"/CFiles")
    task.SetOutputDirectory(OutputDir)
  else:
    print("Will use default output directory.")
  task.SetOutputFileName(OutputFileName)

  # Setting Ranges

  if 'EventPlaneBinLow' in configurations:
    EventPlaneBinLow = configurations['EventPlaneBinLow']
  else:
    EventPlaneBinLow = 1
  if 'EventPlaneBinHigh' in configurations:
    EventPlaneBinHigh = configurations['EventPlaneBinHigh']
  else:
    EventPlaneBinHigh = -1
  task.SetEventPlaneBinRange(EventPlaneBinLow,EventPlaneBinHigh)

  if 'LambdaBinLow' in configurations:
    LambdaBinLow = configurations['LambdaBinLow']
  else:
    LambdaBinLow = 1
  if 'LambdaBinHigh' in configurations:
    LambdaBinHigh = configurations['LambdaBinHigh']
  else:
    LambdaBinHigh = -1
  task.SetLambdaBinRange(LambdaBinLow,LambdaBinHigh)

  if 'EnergyBinLow' in configurations:
    EnergyBinLow = configurations['EnergyBinLow']
  else:
    EnergyBinLow = 1
  if 'EnergyBinHigh' in configurations:
    EnergyBinHigh = configurations['EnergyBinHigh']
  else:
    EnergyBinHigh = -1
  task.SetEnergyBinRange(EnergyBinLow,EnergyBinHigh)

  if 'AsymBinLow' in configurations:
    AsymBinLow = configurations['AsymBinLow']
  else:
    AsymBinLow = 1
  if 'AsymBinHigh' in configurations:
    AsymBinHigh = configurations['AsymBinHigh']
  else:
    AsymBinHigh = -1
  task.SetAsymBinRange(AsymBinLow,AsymBinHigh)

  if 'OpeningAngleBinLow' in configurations:
    OpeningAngleBinLow = configurations['OpeningAngleBinLow']
  else:
    OpeningAngleBinLow = 1
  if 'OpeningAngleBinHigh' in configurations:
    OpeningAngleBinHigh = configurations['OpeningAngleBinHigh']
  else:
    OpeningAngleBinHigh = -1
  task.SetOpeningAngleBinRange(OpeningAngleBinLow,OpeningAngleBinHigh)

  # Other Options
  if 'nPtBins' in configurations:
    task.SetPtBins(configurations['nPtBins'])
  if 'fitMinX' in configurations:
    task.SetFitMinX(configurations['fitMinX'])
  if 'fitMaxX' in configurations:
    task.SetFitMaxX(configurations['fitMaxX'])
  if 'nSigma' in configurations:
    task.SetNSigma(configurations['nSigma'])
  if 'nSigmaR' in configurations:
    task.SetNSigmaR(configurations['nSigmaR'])
  if 'FixedMassWindows' in configurations:
    task.SetFixedMassWindows(configurations['FixedMassWindows'])
  if 'MassRebin' in configurations:
    task.SetNRebinMass(configurations['MassRebin'])
  if 'PtBinChoice' in configurations:
    task.SetPtBinChoice(configurations['PtBinChoice'])
  if 'BkgChoice' in configurations:
    task.SetBkgChoice(configurations['BkgChoice'])
    print("Setting Background choice to %s" % (configurations['BkgChoice']))
  if 'FitMethod' in configurations:
    task.SetFitMethod(configurations['FitMethod'])
  if 'FitPeakMethod' in configurations:
    task.SetFitPeakMethod(configurations['FitPeakMethod'])
  if 'FitBkgMethod' in configurations:
    task.SetFitBkgMethod(configurations['FitBkgMethod'])

  if 'DisableFlow' in configurations:
    task.SetDisableFlow(configurations['DisableFlow'])

  if 'MCPreAnalysisFile' in configurations:
    task.SetMCPreAnalysisFile(configurations['MCPreAnalysisFile'])

  if 'RmvMCPi0PS' in configurations:
    task.SetRemoveMCPi0PS(configurations['RmvMCPi0PS'])
  if 'RmvMCEta' in configurations:
    task.SetRemoveMCEta(configurations['RmvMCEta'])

  if 'EnableThetaModel' in configurations:
    task.SetEnableThetaModel(configurations['EnableThetaModel'])
  if 'EnableThetaLookUp' in configurations:
    task.SetUseThetaLookUpTable(configurations['EnableThetaLookUp'])
  if 'ThetaModelRootFile' in configurations:
    task.SetThetaModelRootFile(configurations['ThetaModelRootFile'])
  if 'ThetaModelParamChoice' in configurations:
    task.SetThetaModelChoice(configurations['ThetaModelParamChoice'])
  if 'ThetaModelTrigger' in configurations:
    task.SetThetaModelTrigger(configurations['ThetaModelTrigger'])
  if 'ThetaModelCent' in configurations:
    task.SetThetaModelCent(configurations['ThetaModelCent'])
  if 'EnablePSScaleMethod' in configurations:
    task.SetEnablePSScaleMethod(configurations['EnablePSScaleMethod'])
  if 'EnablePSDirectMethod' in configurations:
    task.SetEnablePSDirectMethod(configurations['EnablePSDirectMethod'])

  if 'WPi0Power' in configurations:
    task.SetWPower(configurations['WPi0Power'])
  if 'WPi0Yield' in configurations:
    task.SetWYield(configurations['WPi0Yield'])
  if 'WPi0Pt0' in configurations:
    task.SetWPt0(configurations['WPi0Pt0'])
  if 'WPi0Mass' in configurations:
    task.SetWMass(configurations['WPi0Mass'])
  if 'WPi0Sigma' in configurations:
    task.SetWSigma(configurations['WPi0Sigma'])

  if 'label' in configurations:
    task.SetLabel(configurations['label'])
  if 'label2' in configurations:
    task.SetLabel2(configurations['label2'])

  # Get this thing going!

  task.Run()

  print("This was yaml file %s" % (YamlFile),flush=True)

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  runPionID()


