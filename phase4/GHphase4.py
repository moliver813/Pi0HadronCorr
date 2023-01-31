#! /usr/bin/env python3

# ---------------------
# Code to use three Event Plane windows and run RPF fit
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

from reaction_plane_fit import three_orientations

from runPyRPF import *

User=1

from ROOT import gROOT
gROOT.LoadMacro("TaskEventPlane.cxx+")
#if User==0:
#  gROOT.LoadMacro("/Users/Eliane/Software/ALICE_Yale_Dev/analyses/gammaH/PlotGHcorrelation2.cxx+")
#else:
#  gROOT.LoadMacro("/home/moliver/cern/gammaHadron/wrk/phase2/PlotGHcorrelation2.cxx+")
#from ROOT import PlotGHcorrelation2

# ---------------------
#To Do list:
# ---------------------
#-check systematic on fit range for flow function
#-add number of v_ns as flexible parameter. Right now its hard coded

#This macro takes the histograms from the
#gamma and pi0 hadron correlation and analyzes them
#---------------------------------------------------------------------------------------------------
def GHphase3():

  sys.argv.append( '-b' )
  print("................................... ")
  if User==0: print("o Plot_GH_Correlation loaded for: Eliane")
  if User==1: print("o Plot_GH_Correlation loaded for: Michael")

  LocalPath=""
  YamlFile=""

  IsMCGenMode=0

  # Check if a yaml config file has been passed in args:
  print("sys.argv = %d" % len(sys.argv))
  if (len(sys.argv) > 2):
    YamlFile=sys.argv[1]
  else:
    usage="Usage:  %s [yaml config]" % sys.argv[0]
    print(usage)
    exit(0)

  print("o Load config from yamlFile = %s" % YamlFile)

  if User==0: LocalPath = "/Users/Eliane/Software/ALICE_Yale_Dev/analyses/gammaH/"
  if User==1: LocalPath = ""
  LocalPath+=YamlFile
  print("o Load yaml file: %s" % LocalPath)


  # Info on ALICE Published Flow values
  PubFlow="/home/moliver/cern/gammaHadron/wrk/FlowMeasurements/ALICE_PbPb/"



  
  configurations = open(LocalPath)
  configurations = yaml.load(configurations,Loader=yaml.FullLoader)
  #configurations = yaml.load(configurations)
  #configurations = ROOT.LoadFile(LocalPath)
  if not configurations: print(" did not get the yaml file!")

#  TrainHist     = configurations['Train_histo']
  FilePath          = configurations['file_path']

  if 'IsMCGenMode' in configurations:
    IsMCGenMode=configurations['IsMCGenMode']

  # In MCGen Mode the input files will likely all be the same
  InputFile_MC  = ""
  InputFile_A   = ""
  InputFile_0   = ""
  InputFile_1   = ""
  InputFile_2   = ""

  #  InputFile_MC = configurations['file_name_MC']
  InputFile_A   = configurations['file_name_A']
  if (not IsMCGenMode):
    InputFile_0   = configurations['file_name_0']
    InputFile_1   = configurations['file_name_1']
    InputFile_2   = configurations['file_name_2']

  Label         = configurations['label']
#  WagonNameSE   = configurations['wagon_NameSE']
#  WagonNameME   = configurations['wagon_NameME']
 # InputRootFileName = configurations['root_histos']
  if 'output_dir' in configurations:
    OutputDir      = configurations['output_dir']
  else:
    OutputDir      = "output"

  print("  o Load keys from yaml file: ")
  print("  o file_path: %s" % FilePath)
  if (IsMCGenMode):
    print("  o file_name_A: %s" % InputFile_A);
  else:
    print("  o file_name_A: %s" % InputFile_A)
    print("  o file_name_0: %s" % InputFile_0)
    print("  o file_name_1: %s" % InputFile_1)
    print("  o file_name_2: %s" % InputFile_2)


#  InputFileSE   = ROOT.TFile(InFileNameSE)
#  InputFileME   = ROOT.TFile(InFileNameME)
#  if InputRootFileName: InputRootFile = ROOT.TFile(InputRootFileName)

  #..Create the task and run it
  savePlots=1   #..do you want to save the plots to a directory "output/"?
#  observable=0  #..do you want to plot 2DCorr for Ga,Zt, or Xi bins?
#  cent=-1       #..do you want to plot all (-1), central (1),.....
#  evtPl=-1      #..do you want to plot all (-1), in plane (0), mid plane (1), out of plane (2) .....
#  mergeME=0     #..merge ME plots for all bins to have more statistic
#  normMode=0    #..normalize SE after dividing by ME and adding up z-vertex bins (0), or before both (1)
#  plotAdvProj=1 #..plot advanced histograms such as projections
#  plotComp=0    #..plot a comparision of 2 diffent ME correction methods

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

#  print("  o Start analysis while loading files from: ")
 # if TrainHist==0: print("  o Train output")
#  if TrainHist==1: print("  o Histo Root file")

  OutFileName = "output/RPFSub_%s.root" % Label

  if 'output_file' in configurations:
    OutFileName = configurations['output_file']
  else:
    print("Missing output_file name. I don't like this any more")
    return

  OutputFile = ROOT.TFile(OutFileName,"RECREATE")

  # ================================================================
  # Run the C++ Code 
  # ================================================================
  print("Running the C++ RPF implementation",flush=True)
#  task = ROOT.TaskEventPlane(observable)
  task = ROOT.TaskEventPlane()

  if 'debug_level' in configurations:
    task.SetDebugLevel(configurations['debug_level'])
  else:
    task.SetDebugLevel(2)

  task.SetSavePlots(savePlots)
  task.SetPlotOptions("COLZ")

  task.SetOutputDir(OutputDir)
  task.SetOutputFile(OutputFile)

#  if 'UseZYAM' in configurations:
  if 'OverallMode' in configurations:
    task.SetOverallMode(configurations['OverallMode'])
  if (IsMCGenMode):
    task.SetMCGenMode(1)

#  if (IsMCGenMode):
#    InputFileEP_MC = ROOT.TFile(InFileNameEP_MC)
#    InFileNameEP_MC = FilePath+InputFile_MC
#    task.SetIntermediateInputFile(-1,InputFileEP_MC)
#  else:
  InFileNameEP_A = FilePath+InputFile_A
  InputFileEP_A = ROOT.TFile(InFileNameEP_A)
  task.SetIntermediateInputFile(-1,InputFileEP_A)

  if (not IsMCGenMode):
    InFileNameEP_0 = FilePath+InputFile_0
    InFileNameEP_1 = FilePath+InputFile_1
    InFileNameEP_2 = FilePath+InputFile_2

    InputFileEP_0 = ROOT.TFile(InFileNameEP_0)
    InputFileEP_1 = ROOT.TFile(InFileNameEP_1)
    InputFileEP_2 = ROOT.TFile(InFileNameEP_2)

    task.SetIntermediateInputFile(0,InputFileEP_0)
    task.SetIntermediateInputFile(1,InputFileEP_1)
    task.SetIntermediateInputFile(2,InputFileEP_2)



  if 'EnableDeltaEtaScaling' in configurations:
    print("Setting Delta Eta Scaling configuration")
    task.SetEnableDeltaEtaScaling(configurations['EnableDeltaEtaScaling'])

  if 'NumVariants' in configurations:
    print("Setting number of variants to %d" % (configurations['NumVariants']))
    task.SetNumVariants(configurations["NumVariants"])


  if 'EventPlaneResSet' in configurations:
    print("Using Event Plane Resolution Set %d" % configurations['EventPlaneResSet'])
    task.SetEPRSet(configurations['EventPlaneResSet'])
  if 'CentralityBin' in configurations:
    print("Using Centrality Bin %d" % configurations['CentralityBin'])
    task.SetCentralityBin(configurations['CentralityBin'])

  if 'FlowFinderMode' in configurations:
    print("Using Flow Finder Mode %d" % (configurations['FlowFinderMode']))
    task.SetFlowFinderMode(configurations['FlowFinderMode'])


  if 'NRebin' in configurations:
    print("Setting rebin to %d" % (configurations['NRebin']))
    task.SetNRebin(configurations['NRebin'])

  if 'DisableFlow' in configurations:
    print("Setting disable flow to %d. Note that values > 0 disable all vN" % (configurations['DisableFlow']))
    task.SetDisableFlow(configurations['DisableFlow'])

  if 'FixV2TToFirstBin' in configurations:
    task.SetFixV2TMode(configurations['FixV2TToFirstBin'])
    print("Enabling Fixed V2T mode: V2 of trigger determined with first zt bin and fixed across the rest")
  if 'FixV3To0' in configurations:
    task.SetFixV3To0(configurations['FixV3To0'])
    print("Fixing V3 term to 0")


  if 'FlowTermModeTrigger' in configurations:
    task.SetFlowTermModeTrigger(configurations['FlowTermModeTrigger'])
    print("Setting the Flow Term mode to %s" % (configurations['FlowTermModeTrigger'])) 
  if 'FlowTermModeAssoc' in configurations:
    task.SetFlowTermModeAssoc(configurations['FlowTermModeAssoc'])
    print("Setting the Flow Term mode to %s" % (configurations['FlowTermModeAssoc'])) 

  if 'FlowSource' in configurations:
    task.SetFlowSource(configurations['FlowSource']) 


  if 'FlowV1Mode' in configurations:
    task.SetFlowV1Mode(configurations['FlowV1Mode'])

  if 'FixV1Value' in configurations:
    task.SetFlowV1FixValue(configurations['FixV1Value'])
 
  if 'FlowV3Mode' in configurations:
    task.SetFlowV3Mode(configurations['FlowV3Mode'])
  if 'FlowV5Mode' in configurations:
    task.SetFlowV5Mode(configurations['FlowV5Mode'])
  if 'FlowV6TMode' in configurations:
    task.SetFlowV6TMode(configurations['FlowV6TMode'])
  if 'FlowV6AMode' in configurations:
    task.SetFlowV6AMode(configurations['FlowV6AMode'])

  if 'V2ACalcChoice' in configurations:
    task.SetV2ACalcChoice(configurations['V2ACalcChoice'])
  if 'V3CalcChoice' in configurations:
    task.SetV3CalcChoice(configurations['V3CalcChoice'])

  if 'EP1CorrMode' in configurations:
    task.SetEP1CorrMode(configurations['EP1CorrMode'])

  if 'PtBin' in configurations:
    task.SetPtBin(configurations['PtBin'])

  if 'Observable' in configurations:
    task.SetObservable(configurations['Observable'])

  if 'label' in configurations:
    task.SetLabel(configurations['label'])
  if 'label2' in configurations:
    task.SetLabel2(configurations['label2'])

  # Set Trigger Label
  #####  if User==0: task.SetGammaPi0(0)
#  if User==0: task.SetGammaPi0(1)
#  if User==1: task.SetGammaPi0(1)

#  task.SetUseFindLastGoodBin(useFLGB)


  #task.SetWait(1)                       # 1=let the canvases pop up on the screan
  task.Run_Part1()


  # ================================================================
  # Run Raymond's python RPF
  # ================================================================
  print("Running the Python RPF implementation",flush=True)
  #PyRPFOutFileName = "output/PyRPF_RPFSub_%s.root"
  #PyRPFOutputFile = ROOT.TFile(PyRPFOutFileName,"RECREATE")

  # Have input files

  RunRPFCode(task,OutputDir,OutFileName)
#  RunRPFCode(task,OutputDir,PyRPFOutputFile)



#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  GHphase3()

