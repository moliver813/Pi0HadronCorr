#! /usr/bin/env python3

# ---------------------
# import
# ---------------------
import sys
import os
import ROOT
import yaml
import argparse
import statistics
import math

from ROOT import gROOT
from ROOT import gPad
from ROOT import gStyle
from ROOT import TMath
from ROOT import TColor
from ROOT import TCanvas
from ROOT import TString
from ROOT import TFile
from ROOT import TLegend
from ROOT import TMultiGraph
from ROOT import TGraphErrors
from ROOT import TH1D

gROOT.SetBatch(True)

bBinomialDivision = True
enableCFiles = True
enableRootFiles = False
proj2D = False

fDefaultTopMargin=0.05
fDefaultLeftMargin=0.15
fDefaultBottomMargin=0.1
fDefaultRightMargin=0.05

LineStyle = 2
LineWidth = 4
#LineColor = 920
LineColor = ROOT.kGray

# Default Color List (if Custom Color not used)
colorList = [
    ROOT.kBlack,ROOT.kRed,ROOT.kOrange-3,
    ROOT.kGreen-3,ROOT.kBlue,ROOT.kViolet,
    ROOT.kMagenta-5,ROOT.kRed-10,ROOT.kGray,
    ROOT.kCyan,ROOT.kAzure,ROOT.kAzure-8,
    ROOT.kGreen-3,ROOT.kSpring+9,ROOT.kYellow+1]
markerList = [
    ROOT.kFullSquare,ROOT.kFullCircle,ROOT.kFullDiamond,
    ROOT.kOpenSquare,ROOT.kOpenCircle,ROOT.kOpenDiamond,
    ROOT.kFullStar,23,34,
    39,41,43,
    40,42,44]



kMarkerSize = 0.6

kDefaultLabelSizeX=0.03
kDefaultLabelSizeY=0.03
kDefaultTitleSizeX=0.035
kDefaultTitleSizeY=0.035
kDefaultTitleOffsetX=1.1
kDefaultTitleOffsetY=1.1

# whether to use the rotating color palette, or just fixed list of colors.
useCustomColor = True

# setting for space between pads
small = 1e-5

nTypeList=15

c_width=500
c_height=500

def setStyle():
  NRGBs = 5
  NCont = 255

  #gErrorIgnoreLevel

  #Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  #Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  #Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  #Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  #TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);    
  #gStyle.SetNumberContours(NCont);

  gStyle.SetPadBorderMode(0)
  gStyle.SetFrameBorderMode(0)
  gStyle.SetOptStat(0)
  gStyle.SetOptFit(0)

#  gStyle.SetLabelSize(0.05,"X")
#  gStyle.SetLabelSize(0.05,"Y")
  gStyle.SetLabelSize(0.035,"XY")
#  gStyle.SetTitleSize(1.0,"X")
#  gStyle.SetTitleSize(0.3,"Y")

  gStyle.SetOptStat(0);
  gStyle.SetTitleSize(0.035,"xy");


def cleanName(inputString):
  st = TString(inputString)
  st.ReplaceAll(" ","_")
  st.ReplaceAll(" ","_");
  st.ReplaceAll(".","_");
  st.ReplaceAll("/","_");
  st.ReplaceAll("%","_");
  st.ReplaceAll("+","_");
  st.ReplaceAll("#","_");
  st.ReplaceAll(",","_");
  st.ReplaceAll("-","_");
  st.ReplaceAll("*","");
  st.ReplaceAll("(","");
  st.ReplaceAll(")","");
  st.ReplaceAll("{","");
  st.ReplaceAll("}","");
  st.ReplaceAll("<=","lte");
  st.ReplaceAll(">=","gte");
  st.ReplaceAll("<","lt");
  st.ReplaceAll(">","gt");
  st.ReplaceAll("=","_");
  return st


def GetCustomColor(x):
  # r,g,b in [0,255]
  omega_r = 2*0.12; phi_r = -3.1415;
  omega_g = 2*0.08;  phi_g = 1.57;
  omega_b = 1.5*0.16;  phi_b = -3.1415/2.;

  theta_r = 3.1415 * omega_r * x + phi_r;
  theta_g = 3.1415 * omega_g * x + phi_g;
  theta_b = 3.1415 * omega_b * x + phi_b;
  r = 128 + round(127.*TMath.Cos(theta_r));
  g = 128 + round(127.*TMath.Cos(theta_g));
  b = 128 + round(127.*TMath.Cos(theta_b));
  #  printf("r,g,b = (%d,%d,%d)\n",r,g,b);
  return TColor.GetColor(r,g,b);

def GetObjType(obj):
  # removing TH2s for now
  if obj.InheritsFrom("TH2"):
    return 0
  if obj.InheritsFrom("TGraph"):
    return 1
  if obj.InheritsFrom("TH1"):
    return 2
  return 0


def GetMinValue (hist):
  MinBin = hist.GetMinimumBin()
  return hist.GetBinContent(MinBin)

def GetMaxValue (hist):
  MaxBin = hist.GetMaximumBin()
  return hist.GetBinContent(MaxBin)

def ExpandRange(Min,Max):
  MagicNumber = 0.1
  NewMin = Min - MagicNumber*(Max - Min)
  NewMax = Max + MagicNumber*(Max - Min)
  return (NewMin,NewMax)

def RelocateTLegend(graphs,legend):
  return legend

def ProduceSystematicFromHists(hists):
  nHists=len(hists)
  if (nHists < 2):
    return 0
  primaryHist = hists[0]
  newHist = primaryHist.Clone("%s_SysErr" % (primaryHist.GetName()))
  newHist.SetTitle("")
  nPoints = primaryHist.GetNbinsX()
  for i in range(nPoints):
    primaryY = primaryHist.GetBinContent(i+1)
    listOfYValues=[]
    for j in range(len(hists)):
      localHist = hists[j]
      listOfYValues.append(localHist.GetBinContent(i+1))
    print(listOfYValues)
    stdDev=statistics.stdev(listOfYValues)
    print(" stdev = %f" % (stdDev))
    newHist.SetBinError(i+1,stdDev)
  return newHist

def ProduceSystematicFromGraphs(graphs):
  nGraphs=len(graphs)
  if (nGraphs < 2):
    return 0
  primaryGraph = graphs[0]
  newGraph=primaryGraph.Clone("%s_SysErr" % (primaryGraph.GetName()))
  newGraph.SetTitle("")
  nPoints=primaryGraph.GetN()

  for i in range(nPoints):
    primaryY=primaryGraph.GetY()[i]
    listOfYValues=[]
    myVar=0
    myMean=0
    for j in range(len(graphs)):
      localGraph=graphs[j]
      localYValue = localGraph.GetY()[i]
      # Check to avoid crash from giant numbers
      if (math.fabs(localYValue) < 1e50):
        listOfYValues.append(localYValue)
    print(listOfYValues)
    stdDev=statistics.stdev(listOfYValues)
#    print(" stdev = %f" % (stdDev))
    newGraph.SetPointError(i,0.5,stdDev)#FIXME
    # set the x error to be some width characteristic of the input data, or maybe just some fraction of the bin width
    # but I don't have access to the bin width here? I should make sure the phase1 outputs have the bin widths as the errors
    #newGraph.SetPointError(i,newGraph.GetEX()[i],stdDev)
      # temp
#      localValue=localGraph.GetY()[i]
#      myMean+=localValue
#      myVar+=localValue*localValue
#    myMean = myMean / nGraphs
    # manual calculation to check
#    for val in listOfYValues:
#      myVar += (val - myMean)*(val - myMean)     
    
#    myMean/=nGraphs
#    myVar-= myMean*myMean
#    myStdDev = math.sqrt(myVar / (nGraphs-1))
    # should actually use the standard error, the one with sqrt(n-1)
    # ? should this be stdev or pstdev? python documentation implies stdev
#    pstdDev=statistics.pstdev(listOfYValues)
#    print("pstdev = %f" % (pstdDev))
#    print(" mydev = %f" % (myStdDev))

  # Maybe also add a tgraph whose y values are the sys uncertainty
  # that could be used to compare different systematic uncertainties
  # either using this program again or using the paramscan program
  return newGraph

#---------------------------------------------------------------------------------------------------
def sysCompare():

  parser = argparse.ArgumentParser(description='Produce Systematic Uncertainties From Varied inputs.')
  parser.add_argument('-d','--directory',required=False,default="",type=str,help="Directory for comparison plots. If none given, plots are not saved")
  parser.add_argument('-l','--list',required=True,type=str,help="List of histograms/graphs to use")
  parser.add_argument('-o','--output',required=True,type=str,help="Output file to produce")
  parser.add_argument('-r','--ratioMode',required=False,type=bool,help="Whether to produce plots of the ratios")
  parser.add_argument('-t','--titles',required=True,type=str,nargs='+',help="Titles for each file")
  parser.add_argument('files',metavar='Files',type=str,nargs='+',help="Files to use")

  args = parser.parse_args()

  setStyle()

  stringListOfHists=args.list
  #list of hists/graphs
  listOfHists=stringListOfHists.split()
  #array of files
  fileNames=args.files

  fileTitles=args.titles

  directory=args.directory

  print("List of hists/graphs to use:"),
  print(listOfHists)

  print("List of files to use:")
  print(fileNames)
  print("That is %d files" % (len(fileNames)))
  print("List of titles:")
  print(fileTitles)
  print("That is %d titles" % (len(fileTitles)))

  if (len(fileNames) != len(fileTitles)):
    print("Error: mismatch in number of files vs titles")
    exit(1)


  print("Primary file:")
  print(fileNames[0])

  # list of files
  files=[]

  titlesToFileNames=zip(fileTitles,fileNames)
  titlesToFileNames=set(titlesToFileNames)
  print(titlesToFileNames)

  fileNamesToTitles=zip(fileNames,fileTitles)
  fileNamesToTitles=dict(fileNamesToTitles)  

  filesToTitles={}
  cleanTitlesToFiles={}
  filesToCleanTitles={}

  for filename in fileNames:
    print("Opening File %s" % (filename))
    tfile = TFile.Open(filename,"READ")
    print("opened file %s" % (tfile.GetName()))
    files.append(tfile)

    filesToTitles[tfile] = fileNamesToTitles[filename]
    # clean name that can be added to root object names
    cleanedTitle = cleanName(fileNamesToTitles[filename])
    filesToCleanTitles[tfile] = cleanedTitle

  primaryFile = files[0]

  if (directory == ""):
    print("No directory given. Output plots will not be saved")
  else:
    print("Directory for output %s" % (directory))

  # create output file here, before changing directory

  canvas=TCanvas("canvas","canvas",c_width,c_height)
  if (directory != ""):
    # check if the directory exists
    if (not os.path.isdir(directory)):
      os.makedirs(directory)
    os.chdir(directory)

  # This is where it gets easier than compare: just need to find the objects in the files
  
  for objName in listOfHists:
    print("Starting the thing for object %s" % (objName))
    #print(" Looking in file %s" % (primaryFile.GetName()))
    
    primaryObj = primaryFile.Get(objName)
    if (not primaryObj):
      print("Could not find object %s in file %s" % (objName,fileNames[0]))
    else:
      print("obj title = %s" % (primaryObj.GetTitle()))
    objectTitle=primaryObj.GetTitle()
    listOfObjs=[]
    for tfile in files:
      localObj=tfile.Get(objName) 
      if (not localObj):
        print("Could not find object %s in file %s" % (objName,tfile))
        exit(1)
      else:
        print("obj title = %s" % (localObj.GetTitle()))
      localObj.SetName("%s_%s" % (objName,filesToCleanTitles[tfile]))
      print("  name set to %s" % (localObj.GetName()))
      localObj.SetTitle(filesToTitles[tfile])
      print("  title set to %s" % (localObj.GetTitle()))
      listOfObjs.append(localObj)


    legX=0.6
    legY=0.22
    legWidth=0.25
    legHeight=0.25

    # legend for comparison
#    leg = TLegend(legX,legY,legX+legWidth,legY+legHeight)
    # legend for SysUncert
#    leg2 = TLegend(legX,legY,legX+legWidth,legY+legHeight)
    # legend for comparison
    leg = TLegend(legWidth,legHeight,legWidth,legHeight)
    # legend for SysUncert
    leg2 = TLegend(legWidth,legHeight,legWidth,legHeight)

    # Set Properties and include in legend
    for i in range(len(listOfObjs)):
      localObj=listOfObjs[i]
      localTitle=fileTitles[i]
      color=ROOT.kBlack
      if (i != 0):
        if (useCustomColor):
          color = GetCustomColor(i)
        else:
          color = colorList[i]
      markerStyle = markerList[i]
      localObj.SetMarkerColor(color)
      localObj.SetLineColor(color)
      localObj.SetMarkerStyle(markerStyle)

      leg.AddEntry(localObj,localTitle,"LP")

    leg2.AddEntry(listOfObjs[0],fileTitles[0],"LP")
    # Two key cases: TGraph or TH1
    # For TGraphs we want to do a multigraph
    # for th1, we have to manually track the y-limits
    iObjType=GetObjType(listOfObjs[0])
    # tgraph : 1
    # th1    : 2
    # other  : 0
    
    canvas.Clear()
    canvas.cd()

    if (iObjType == 1): # TGraph
      # Get object name
      #objName = "TestGraph"

      # Reset the margins
      gPad.SetTopMargin(fDefaultTopMargin)
      gPad.SetLeftMargin(fDefaultLeftMargin)
      gPad.SetBottomMargin(fDefaultBottomMargin)
      gPad.SetRightMargin(fDefaultRightMargin)

      mg = TMultiGraph()
      mg.SetTitle(primaryObj.GetTitle())
      mg.GetXaxis().SetTitle(primaryObj.GetXaxis().GetTitle())
      mg.GetYaxis().SetTitle(primaryObj.GetYaxis().GetTitle())
      for lobj in listOfObjs:
        mg.Add(lobj)
      mg.Draw("ALP")
      #leg.Draw("SAME")
      legtest = gPad.BuildLegend()
      legtest.Draw("SAME")
      if (directory != ""):
        canvas.Print("%s_Cmp.pdf" % (objName))
        canvas.Print("%s_Cmp.png" % (objName))

      # Now produce systematic uncertainties (for TGraphErrors
      sysUncertObj=ProduceSystematicFromGraphs(listOfObjs)
      sysUncertObj.SetTitle("Systematic Uncertainty")
      sysUncertObj.SetFillColor(ROOT.kBlue)
      sysUncertObj.SetFillStyle(3002)
      sysUncertObj.Draw("ALP[]5")

      # reset the primary objects title
      primaryObj.Draw("LP")
      leg2.AddEntry(sysUncertObj,"Systematic Uncertainty","F")
      #mg.Draw("LP")
      #leg.Draw("SAME")
      leg2.Draw("SAME")
      # not using the autobuild legend here. could get a good location
      # from the comparison plot
      #legtest2=gPad.BuildLegend()
      #legtest2.Draw("SAME")
      if (directory != ""):
        canvas.Print("%s_SysUncert.pdf" % (objName))
        canvas.Print("%s_SysUncert.png" % (objName))

      # Now plot them both in a split canvas
      canvas.Clear()
      canvas.Divide(1,2,0.01,0.0)
      canvas.cd(1)
      mg.Draw("ALP X+")
      gPad.SetTopMargin(0.0)
      leg.Draw("SAME")
      gPad.SetBottomMargin(0.0)
      canvas.cd(2)
      sysUncertObj.Draw("ALP[]5")
      gPad.SetTopMargin(0.0)
      gPad.SetBottomMargin(0.1)
      primaryObj.Draw("LP")
      leg2.Draw("SAME")
      if (directory != ""):
        canvas.Print("%s_SysUncert_Cmp.pdf" % (objName))
        canvas.Print("%s_SysUncert_Cmp.png" % (objName))



    # for histograms, also draw a plot with each of them
    # separately? Useful if the fit functions are visible
    if (iObjType == 2): # TH1
      primaryObj.Draw()
      fYMin = GetMinValue(primaryObj)
      fYMax = GetMaxValue(primaryObj)
      for lobj in listOfObjs:
        lobj.Draw("SAME")
        fYMin = min(fYMin,GetMinValue(lobj))
        fYMax = max(fYMax,GetMaxValue(lobj))
      (fYMin,fYMax) = ExpandRange(fYMin,fYMax)
      leg.Draw("SAME")
      primaryObj.GetYaxis().SetRangeUser(fYMin,fYMax)
      if (directory != ""):
        canvas.Print("%s_Cmp.pdf" % (objName))
        canvas.Print("%s_Cmp.png" % (objName))
      sysUncertObj=ProduceSystematicFromHists(listOfObjs)
      sysUncertObj.SetFillColor(ROOT.kBlue)
      sysUncertObj.SetFillStyle(3002)
      #sysUncertObj.Draw("E2")
      sysUncertObj.Draw("E4")
      primaryObj.Draw("SAME")
      leg2.AddEntry(sysUncertObj,"Systematic Uncertainty","F")
      #leg.Draw("SAME")
      leg2.Draw("SAME")
      if (directory != ""):
        canvas.Print("%s_SysUncert.pdf" % (objName))
        canvas.Print("%s_SysUncert.png" % (objName))
      

#    primaryObj.Draw()
#    for lobj in listOfObjs:
#      lobj.Draw("SAME")
#    canvas.Print("Test.pdf")




#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  sysCompare()



exit()
