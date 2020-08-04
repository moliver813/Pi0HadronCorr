#! /usr/bin/env python3

# ---------------------
#  Running the macro
# ---------------------
# python GHcorr_PlotMacro.py
# python GHcorr_PlotMacro.py [yaml-file]

# ---------------------
# import
# ---------------------
import sys
sys.argv.append( '-b' )
import os
import ROOT
import yaml

# User = 0 for Eliane
# User = 1 for Michael
User=1

from ROOT import gROOT
if User==0:
  gROOT.LoadMacro("/Users/Eliane/Software/ALICE_Yale_Dev/analyses/gammaH/PlotGHcorrelation2.cxx+")
else:
  gROOT.LoadMacro("~/cern/gammaHadron/wrk/phase2/PlotGHcorrelation2.cxx+")
 # gROOT.LoadMacro("/home/moliver/cern/gammaHadron/wrk/phase2/PlotGHcorrelation2.cxx+")
from ROOT import PlotGHcorrelation2

# ---------------------
#To Do list:
# ---------------------
#-check systematic on fit range for flow function
#-add number of v_ns as flexible parameter. Right now its hard coded

#//-add scale as a parameter to function to pass on
#//-add number of v_ns as flexible parameter. Right now its hard coded
#//-check with porjection what the ranges are included/excluded


#This macro takes the histograms from the
#gamma and pi0 hadron correlation and analyzes them
#---------------------------------------------------------------------------------------------------
def GHcorr_PlotMacro():

  print("................................... ")
  if User==0: print("o Plot_GH_Correlation loaded for: Eliane")
  if User==1: print("o Plot_GH_Correlation loaded for: Michael")

  LocalPath=""
  YamlFile=""

  # Check if a yaml config file has been passed in args:
  print("sys.argv = %d" % len(sys.argv))
  if (len(sys.argv) > 2):
    YamlFile=sys.argv[1]
  if (YamlFile==""):
    if User==0:
      YamlFile="configGH_EE.yaml"
    else:
      YamlFile="configGH_MO.yaml"

  print("o Load config from yamlFile = %s" % YamlFile)

  if User==0: LocalPath = "/Users/Eliane/Software/ALICE_Yale_Dev/analyses/gammaH/"
  if User==1: LocalPath = ""
  LocalPath+=YamlFile
  print("o Load yaml file: %s" % LocalPath)
  
  configurations = open(LocalPath)
  configurations = yaml.load(configurations,Loader=yaml.FullLoader)
  #configurations = ROOT.LoadFile(LocalPath)
  if not configurations: print(" did not get the yaml file!")

  GammaPiSB     = configurations['gamma_Pi_SB']
  TrainHist     = configurations['Train_histo']
  FilePath      = configurations['file_path']

	# FIXME in progress: for alidock, replace moliver with alidock?


  InputFileSE   = configurations['file_nameSE']
  InputFileME   = configurations['file_nameME']
  WagonNameSE   = configurations['wagon_NameSE']
  WagonNameME   = configurations['wagon_NameME']
  InputRootFileName = "Projections/%s" % (configurations['root_histos'])
  if (not os.path.isdir("Projections")):
    os.mkdir("Projections")
  InputRootFileMEName = ""
  if 'root_histos_ME' in configurations:
    InputRootFileMEName = configurations['root_histos_ME']
  OutputProjectionLabel = configurations['label']
  if 'output_dir' in configurations:
    OutputDir     = configurations['output_dir']
  else:
    OutputDir = "output"
  # more variables
  savePlots     =configurations['savePlots']  #..do you want to save the plots to a directory "output/"?
  observable    =configurations['observable'] #..do you want to plot 2DCorr for Ga,Zt, or Xi bins?
  cent          =configurations['cent']       #..do you want to plot all (-1), central (1),.....
  evtPl         =configurations['evtPl']      #..do you want to plot all (-1), in plane (0), mid plane (1), out of plane (2) .....
  mergeME       =configurations['mergeME']    #..merge ME plots for all bins to have more statistic
  normMode      =configurations['normMode']    #..normalize SE after dividing by ME and adding up z-vertex bins (0), or before both (1)
  plotAdvProj   =configurations['plotAdvProj'] #..plot advanced histograms such as projections
  plotComp      =configurations['plotComp']    #..plot a comparision of 2 diffent ME correction methods
  useFLGB       =configurations['useFLGB']     # use FindLastGoodBin method to set delta eta limits, or use default delta eta limit
  ptMinBin      =configurations['ptMinBin']
  ptMaxBin      =configurations['ptMaxBin']

	
 

  print("o Load keys from yaml file: ")
  print("  o Train_or_histo: %s" % TrainHist)
  if TrainHist==0: print("  o file_path: %s" % FilePath)
  if TrainHist==0: print("  o file_nameSE: %s" % InputFileSE)
  if TrainHist==0: print("  o file_nameME: %s" % InputFileME)
  if TrainHist==0: print("  o wagon_NameSE: %s" % WagonNameSE)
  if TrainHist==0: print("  o wagon_NameME: %s" % WagonNameME)
  if TrainHist==1: print("  o Input root_histos: %s" % InputRootFileName)
  print("  o output_dir: %s" % OutputDir)
  print("  o output_label: %s" % OutputProjectionLabel)

  print("o Load analysis variables from yaml file:")
  print("  o savePlots: %s" % savePlots) 
  print("  o observable: %s" % observable) 
  print("  o cent: %s" % cent) 
  print("  o evtPl: %s" % evtPl) 
  print("  o mergeME: %s" % mergeME) 
  print("  o normMode: %s" % normMode)  ##this needs to be deleted out of the yaml file since it is no longer used
  print("  o plotAdvProj: %s" % plotAdvProj) 
  print("  o plotComp: %s" % plotComp) 
  print("  o ptBin Range: %s - %s" % (ptMinBin, ptMaxBin)) 
  

  InFileNameSE = FilePath+InputFileSE
  InFileNameME = FilePath+InputFileME
 
  InputFileSE   = ROOT.TFile(InFileNameSE)
  InputFileME   = ROOT.TFile(InFileNameME)

	# if runing in train mode, recreate the file
  FileMode = "READ"
  if TrainHist==0:
    FileMode = "RECREATE"
  # if using pre-projected files, merely open the file
  print("filemode = %s" % FileMode)
  InputRootFile = 0
  InputRootFileME = 0
  if InputRootFileName: InputRootFile = ROOT.TFile(InputRootFileName,FileMode)
#  if InputRootFileName: InputRootFile = ROOT.TFile(InputRootFileName)

  if InputRootFileMEName != "" : InputRootFileME = ROOT.TFile(InputRootFileMEName,FileMode)

  if InputRootFile==0: # FIXME this doesn't work
    print("Input File %s Missing!" % InputRootFileName)
    exit()

  #..Create the task and run it
  #  savePlots=1   #..do you want to save the plots to a directory "output/"?
  #  observable=0  #..do you want to plot 2DCorr for Ga,Zt, or Xi bins?
  #  cent=-1       #..do you want to plot all (-1), 0-10% (0), 10-30% (1) .....
  #  evtPl=-1      #..do you want to plot all (-1), in plane (0), mid plane (1), out of plane (2) .....
  #  mergeME=0     #..merge ME plots for all bins to have more statistic
  #  plotAdvProj=1 #..plot advanced histograms such as projections
  #  plotComp=0    #..plot a comparision of 2 diffent ME correction methods

  if (savePlots):
    if (not os.path.isdir("output")):
      os.mkdir("output")
#    os.mkdir("output","-p")
      

  print("o Start analysis while loading files from: ")
  if TrainHist==0: print("  o Train output")
  if TrainHist==1: print("  o Histo Root file")

  task = ROOT.PlotGHcorrelation2(observable,cent,evtPl) #..set options
  if TrainHist==0: task.SetInputTrainFile(InputFileSE,InputFileME,WagonNameSE,WagonNameME)
  if TrainHist==1:
    task.SetInputRootFile(InputRootFile)
    if (InputRootFileME != 0): task.SetInputRootFileME(InputRootFileME)
  task.SetOutputLabel(OutputProjectionLabel)
	# Make output dir if it does not already exist
  if not os.path.exists(OutputDir):
    os.mkdir(OutputDir)
  if not os.path.exists(OutputDir+"/CFiles"):
    os.mkdir(OutputDir+"/CFiles")
  task.SetOutputDir(OutputDir)
  # Set Trigger Label
  task.SetGammaPi0(GammaPiSB)      
  task.SetSavePlots(savePlots)
 # task.SetPlotOptions("COLZ")
  task.SetPlotOptions("LEGO2 FB BB")
  task.SetMergeMEplots(mergeME)         # do not merge mixed event plots

  if 'mcMode' in configurations:
    print("Setting up MC Mode %d" % (configurations['mcMode']))
    task.SetMCMode(configurations['mcMode'])
  if 'nSigma' in configurations:
    print("Setting NSigma to %f" % (configurations['nSigma'])) 
    task.SetNSigma(configurations['nSigma'])
  if 'SetPlot2DHistos' in configurations:
    task.SetPlot2DHistos(configurations['SetPlot2DHistos'])
  if 'MEDEtaRange' in configurations:
    task.SetMEDetaRangeForNorm(configurations['MEDEtaRange'])
  if 'RebinMEForNorm' in configurations:
    task.SetRebinMEForNorm(configurations['RebinMEForNorm'])
  if 'MENormSideBins' in configurations:
    task.SetMENormSideBins(configurations['MENormSideBins'])
  if 'MinDEtaSignalRange' in configurations:
    task.SetMinDEtaSignalRange(configurations['MinDEtaSignalRange'])
  if 'MaxDEtaSignalRange' in configurations:
    task.SetMaxDEtaSignalRange(configurations['MaxDEtaSignalRange'])


  task.SetPtBinRange(ptMinBin, ptMaxBin);
  task.SetPlotAdvanced(plotAdvProj)     # plot advanced histos such as corrected projections to Delta Eta and Delta phi
  task.SetPlotMEcomparision(plotComp)   # plot a comparision between adding all z-vertices or not
  task.SetUseFindLastGoodBin(useFLGB)
  #task.SetWait(1)                       # 1=let the canvases pop up on the screan
  task.Run()

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  GHcorr_PlotMacro()

