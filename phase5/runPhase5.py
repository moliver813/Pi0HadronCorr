#! /usr/bin/env python3

# ---------------------
# Code to use three Event Plane windows and run RPF fit
#
#
# ---------------------

# ---------------------
#  Running the macro
# ---------------------
# python runPhase5.py
# python runPhase5.py [yaml-file]

# ---------------------
# import
# ---------------------
import sys
sys.argv.append( '-b' )
import os
import ROOT
import yaml

User=1

from ROOT import TString
from ROOT import gROOT
gROOT.LoadMacro("TaskCalcObservables.cxx+")

# ---------------------
#To Do list:
# ---------------------



# applies formatting options for display
def formatTitle(inputString):
  st = TString(inputString)
  st.ReplaceAll("_"," ")
  return st.Data()


#This macro takes the histograms from the
#gamma and pi0 hadron correlation and analyzes them
#---------------------------------------------------------------------------------------------------
def GHphase5():

  sys.argv.append( '-b' )
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

  LocalPath = ""
  LocalPath+=YamlFile
  print("o Load yaml file: %s" % LocalPath)
  
  configurations = open(LocalPath)
  configurations = yaml.load(configurations,Loader=yaml.FullLoader)
  #configurations = yaml.load(configurations)
  #configurations = ROOT.LoadFile(LocalPath)
  if not configurations: print(" did not get the yaml file!")

#  TrainHist     = configurations['Train_histo']
  FilePath          = configurations['file_path']
  InputFile_Central = configurations['file_name_central']
  Label             = configurations['label']
  if 'output_dir' in configurations:
    OutputDir       = configurations['output_dir']
  else:
    OutputDir       = "output"
  RPFMethod         = configurations['rpf_method']


  print("  o Load keys from yaml file: ")
  print("  o file_path: %s" % FilePath)
  print("  o file_name_Central: %s" % InputFile_Central)

  

  InFileName_Central = FilePath+InputFile_Central
 
  InputFile_Central = ROOT.TFile(InFileName_Central)

#  if InputRootFileName: InputRootFile = ROOT.TFile(InputRootFileName)

  #..Create the task and run it
  savePlots=1   #..do you want to save the plots to a directory "output/"?
#  observable=0  #..do you want to plot 2DCorr for Ga,Zt, or Xi bins?
#  cent=-1       #..do you want to plot all (-1), central (1),.....

#  if (savePlots):
  if (not os.path.isdir("output")):
    os.mkdir("output")
#    os.mkdir("output","-p")
      

  if (not os.path.isdir(OutputDir)):
    os.makedirs(OutputDir)
#    os.mkdir(OutputDir)
  if (not os.path.isdir(OutputDir+"/CFiles")):
    os.mkdir(OutputDir+"/CFiles")
  if (not os.path.isdir(OutputDir+"/QA")):
    os.mkdir(OutputDir+"/QA")

  OutFileName = "output/FinalObs_%s.root" % Label

  if 'output_file' in configurations:
    OutFileName = configurations['output_file']
#  else:
#    print("Missing output_file name. I don't like this any more")
#    return

  OutputFile = ROOT.TFile(OutFileName,"RECREATE")

  # ================================================================
  # Run the C++ Code 
  # ================================================================
  print("Running the Phase 5 code",flush=True)
#  task = ROOT.TaskEventPlane(observable)
  task = ROOT.TaskCalcObservables()
  task.SetDebugLevel(2)
#  task.SetSavePlots(savePlots)
  task.SetPlotOptions("COLZ")

  task.SetOutputDir(OutputDir)
  task.SetOutputFile(OutputFile)

  if 'CentralityBin' in configurations:
    print("Using Centrality Bin %d" % configurations['CentralityBin'])
    task.SetCentralityBin(configurations['CentralityBin'])

  task.SetCentralInputFile(InputFile_Central)

  task.SetRPFMethod(RPFMethod);


  # Inputting systematics.
  if "systematic_uncertainties" in configurations:
    for key in configurations["systematic_uncertainties"]:
      SysUncertFilePath = configurations["systematic_uncertainties"][key]
      if (not os.path.exists(SysUncertFilePath)):
        print("Systematic Uncertainty %s File %s does not exist!" % (key,SysUncertFilePath))
      else:
        SysUncertFile = ROOT.TFile(SysUncertFilePath)
        task.AddSystematicComparison(SysUncertFile,key)

  # Inputting models.
  if "model_comparisons" in configurations:
    for key in configurations["model_comparisons"]:
      ModelFilePath = configurations["model_comparisons"][key]
      if (not os.path.exists(ModelFilePath)):
        print("Model %s File %s does not exist!" % (key,ModelFilePath))
      else:
        print("Adding model %s" % (key))
        ModelFile = ROOT.TFile(ModelFilePath)
        ModelTitle=formatTitle(key)
        #task.AddModelComparison(ModelFile,key,key)
        task.AddModelComparison(ModelFile,key,ModelTitle)


  if 'PlotStatus' in configurations:
    task.SetPlotStatus(configurations['PlotStatus'])

  if 'PtBin' in configurations:
    task.SetPtBin(configurations['PtBin'])

  if 'EventPlaneResSet' in configurations:
    print("Using Event Plane Resolution Set %d" % configurations['EventPlaneResSet'])
    task.SetEPRSet(configurations['EventPlaneResSet'])

  if 'Preliminary' in configurations:
    task.SetPreliminary(configurations['Preliminary'])

  if 'CentralityBin' in configurations:
    task.SetCentralityBin(configurations['CentralityBin'])

  if 'SigmaRangeNSBinChange' in configurations:
    task.SetSigmaRangeNSBinChange(configurations['SigmaRangeNSBinChange'])
  if 'SigmaRangeASBinChange' in configurations:
    task.SetSigmaRangeASBinChange(configurations['SigmaRangeASBinChange'])

  if 'AwaysideFit' in configurations:
    task.SetAwaysideFitOption(configurations['AwaysideFit'])
  if 'NearsideFit' in configurations:
    task.SetNearsideFitOption(configurations['NearsideFit'])

  if 'SkipPoints' in configurations:
    task.SetNSkipPoints(configurations['SkipPoints'])

  if 'Blinded' in configurations:
    task.SetBlinded(configurations['Blinded'])

  if 'MC' in configurations:
    task.SetMCGenMode(configurations['MC'])

  if 'label' in configurations:
    task.SetLabel(configurations['label'])
  if 'label2' in configurations:
    task.SetLabel2(configurations['label2'])

  if 'title' in configurations:
    task.SetMyTitle(configurations['title'])
  if 'title2' in configurations:
    task.SetMyTitle2(configurations['title2'])

# Example of setting something based on yaml config
#  if 'FixV2TToFirstBin' in configurations:
#    task.SetFixV2TMode(configurations['FixV2TToFirstBin'])
#    print("Enabling Fixed V2T mode: V2 of trigger determined with first zt bin and fixed across the rest")



  #task.SetWait(1)                       # 1=let the canvases pop up on the screan
  task.Run()




#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  GHphase5()

