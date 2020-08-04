#! /usr/bin/env python3

# ---------------------
# Code to use Sideband intermediate resultsto produce background to subtract from pi0_cand-hadron correlations
# Uses input purity (stored in a root file with tgraphs) to do subtractions
#
#
# ---------------------

# ---------------------
#  Running the macro
# ---------------------
# python GHphase3.py
# python GHphase3.py [yaml-file]

# ---------------------
# import
# ---------------------
import sys
sys.argv.append( '-b' )
import os
import ROOT
import yaml

User=1

from ROOT import gROOT
gROOT.LoadMacro("TaskSideband.cxx+")

# ---------------------
#To Do list:
# ---------------------

#This macro takes the histograms from the
#gamma and pi0 hadron correlation and analyzes them
#---------------------------------------------------------------------------------------------------
def GHsidebandSub():

  print("................................... ")

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
#  if User==0: LocalPath = "/Users/Eliane/Software/ALICE_Yale_Dev/analyses/gammaH/"
#  if User==1: LocalPath = ""
  LocalPath+=YamlFile
  print("o Load yaml file: %s" % LocalPath)
  
  configurations = open(LocalPath)
  configurations = yaml.load(configurations,Loader=yaml.FullLoader)
  #configurations = ROOT.LoadFile(LocalPath)
  if not configurations: print(" did not get the yaml file!")

#  TrainHist     = configurations['Train_histo']
  FilePath          = configurations['file_path']

  InputFileName_Pi0  = configurations['file_name_pi0']
  SidebandFileName_3 = configurations['file_name_SB3']
  SidebandFileName_4 = configurations['file_name_SB4']
  SidebandFileName_5 = configurations['file_name_SB5']
  SidebandFileName_6 = configurations['file_name_SB6']

# File with S/(S+B) not assumed to be in same directory
  FilePathPurity    = configurations['file_path_Purity']

  BkgSelection      = configurations['background_selection']
  FitFunction       = 0
  FitSidebandsSel   = 0
  if "background_selection" in configurations:
    FitFunction = configurations['fit_function']
  if "fit_sidebands_selection" in configurations:
    FitSidebandsSel = configurations['fit_sidebands_selection']

#  InputFile_A   = configurations['file_name_A']
#  InputFile_0   = configurations['file_name_0']
#  InputFile_1   = configurations['file_name_1']
#  InputFile_2   = configurations['file_name_2']
#  WagonNameSE   = configurations['wagon_NameSE']
#  WagonNameME   = configurations['wagon_NameME']
 # InputRootFileName = configurations['root_histos']
  if 'output_dir' in configurations:
    OutputDir      = configurations['output_dir']
    if not os.path.exists(OutputDir):
      os.makedirs(OutputDir)
#      os.mkdir(OutputDir)
    CFilesDir = OutputDir + "/CFiles"
    if not os.path.exists(CFilesDir):
      os.makedirs(CFilesDir)
#      os.mkdir(CFilesDir)
  else:
    OutputDir      = "output"
  if 'output_file' in configurations:
    OutputFileName = configurations['output_file']
  else:
    OutputFileName = "SidebandOutput.root"


  print("  o Load keys from yaml file: ")
  print("  o file_path: %s" % FilePath)
  print("  o Purity File: %s" % FilePathPurity)
  print("  o Pi0 File: %s" % InputFileName_Pi0)
  print("  o Sideband 3 File: %s" % SidebandFileName_3)
  print("  o Sideband 4 File: %s" % SidebandFileName_4)
  print("  o Sideband 5 File: %s" % SidebandFileName_5)
  print("  o Sideband 6 File: %s" % SidebandFileName_6)
  print("  o Output File: %s" % OutputFileName)
  print("  o Bkg Selection: %d" % BkgSelection)
  print("  o Fit Function: %d" % FitFunction)
  print("  o FitSideband Selection: %d" % FitSidebandsSel)
  

  InputFile_Pi0  = ROOT.TFile(FilePath+InputFileName_Pi0)
  SidebandFile_3 = ROOT.TFile(FilePath+SidebandFileName_3)
  SidebandFile_4 = ROOT.TFile(FilePath+SidebandFileName_4)
  SidebandFile_5 = ROOT.TFile(FilePath+SidebandFileName_5)
  SidebandFile_6 = ROOT.TFile(FilePath+SidebandFileName_6)

  
  Pi0PurityFile  = ROOT.TFile(FilePathPurity)

  #..Create the task and run it


#  if (savePlots):
  if (not os.path.isdir("output")):
    os.mkdir("output")
#    os.mkdir("output","-p")
      

  if (not os.path.isdir(OutputDir)):
    os.makedirs(OutputDir)
#    os.mkdir(OutputDir)
  if (not os.path.isdir(OutputDir+"/CFiles")):
    os.makedirs(OutputDir+"/CFiles")
#    os.mkdir(OutputDir+"/CFiles")


  task = ROOT.TaskSideband()
  task.SetDebugLevel(2)
#  task.SetSavePlots(savePlots)
  task.SetPlotOptions("COLZ")

  if 'ptbin' in configurations:
    task.SetPtBin(configurations['ptbin'])
  if 'cent' in configurations:
    task.SetCentralityBin(configurations['cent'])
  if 'label' in configurations:
    task.SetLabel(configurations['label'])
  if 'label2' in configurations:
    task.SetLabel2(configurations['label2'])

  if 'purityChoice' in configurations:
    task.SetPurityChoice(configurations['purityChoice'])
  if 'fixedpurity' in configurations:
    task.SetFixedPurity(configurations['fixedpurity'])

  task.SetBackgroundSelection(BkgSelection)
  task.SetScalingFitFunction(FitFunction)
  task.SetSidebandFitMask(FitSidebandsSel)

  task.SetOutputDir(OutputDir)

#  task.SetIntermediateInputFile(-1,InputFileEP_A)
#  task.SetIntermediateInputFile(0,InputFileEP_0)
#  task.SetIntermediateInputFile(1,InputFileEP_1)
#  task.SetIntermediateInputFile(2,InputFileEP_2)

  task.SetPi0CorrInputFile(InputFile_Pi0)
  
  task.SetSidebandInputFile(0,SidebandFile_3)
  task.SetSidebandInputFile(1,SidebandFile_4)
  task.SetSidebandInputFile(2,SidebandFile_5)
  task.SetSidebandInputFile(3,SidebandFile_6)

  task.SetPi0PurityFile(Pi0PurityFile)

  task.SetOutputFileName(OutputFileName)	


  #task.SetWait(1)                       # 1=let the canvases pop up on the screan
  task.Run()

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  GHsidebandSub()

