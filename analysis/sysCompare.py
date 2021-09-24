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
from ROOT import TPaveText
from ROOT import TMultiGraph
from ROOT import TGraphErrors
from ROOT import TH1D

gROOT.SetBatch(True)

# Default behavior is to use the first file's data for the central value in results
CenterValueMode=0

bBinomialDivision = True
enableCFiles = True
enableRootFiles = False
proj2D = False

fDefaultTopMargin=0.05
fDefaultLeftMargin=0.15
fDefaultBottomMargin=0.1
fDefaultRightMargin=0.05

canvSmall=1.0e-5

TitlePaveMinX=0.5
TitlePaveMaxX=0.95
TitlePaveMinY=0.91
TitlePaveMaxY=0.99
TitlePaveLineStyle=2
TitlePaveLineWidth=1

LineStyle = 2
LineWidth = 4
#LineColor = 920
LineColor = ROOT.kGray

AxisLabelSizeX=0.035
AxisLabelSizeY=0.035

Legend1Width=0.3
Legend1Height=0.5


# Default Color List (if Custom Color not used)
colorList = [
    ROOT.kBlack,ROOT.kRed,ROOT.kOrange-3,
    ROOT.kGreen-3,ROOT.kBlue,ROOT.kViolet,
    ROOT.kMagenta-5,ROOT.kRed-10,ROOT.kGray,
    ROOT.kCyan,ROOT.kAzure,ROOT.kAzure-8,
    ROOT.kGreen-3,ROOT.kSpring+9,ROOT.kYellow+1,
    ROOT.kBlack,ROOT.kRed,ROOT.kOrange-3,
    ROOT.kGreen-3,ROOT.kBlue,ROOT.kViolet,
    ROOT.kMagenta-5,ROOT.kRed-10,ROOT.kGray,
    ROOT.kCyan,ROOT.kAzure,ROOT.kAzure-8,
    ROOT.kGreen-3,ROOT.kSpring+9,ROOT.kYellow+1,
    ROOT.kBlack,ROOT.kRed,ROOT.kOrange-3,
    ROOT.kGreen-3,ROOT.kBlue,ROOT.kViolet,
    ROOT.kMagenta-5,ROOT.kRed-10,ROOT.kGray,
    ROOT.kCyan,ROOT.kAzure,ROOT.kAzure-8,
    ROOT.kGreen-3,ROOT.kSpring+9,ROOT.kYellow+1
    ]
markerList = [
    ROOT.kFullSquare,ROOT.kFullCircle,ROOT.kFullDiamond,
    ROOT.kOpenSquare,ROOT.kOpenCircle,ROOT.kOpenDiamond,
    ROOT.kFullStar,23,34,
    39,41,43,
    40,42,44,
    ROOT.kFullSquare,ROOT.kFullCircle,ROOT.kFullDiamond,
    ROOT.kOpenSquare,ROOT.kOpenCircle,ROOT.kOpenDiamond,
    ROOT.kFullStar,23,34,
    39,41,43,
    40,42,44,
    ROOT.kFullSquare,ROOT.kFullCircle,ROOT.kFullDiamond,
    ROOT.kOpenSquare,ROOT.kOpenCircle,ROOT.kOpenDiamond,
    ROOT.kFullStar,23,34,
    39,41,43,
    40,42,44
    ]



kMarkerSize = 1.0

kDefaultLabelSizeX=0.03
kDefaultLabelSizeY=0.03
kDefaultTitleSizeX=0.055
kDefaultTitleSizeY=0.055
#kDefaultTitleSizeX=0.035
#kDefaultTitleSizeY=0.035
kDefaultTitleOffsetX=0.75
kDefaultTitleOffsetY=0.9

# whether to use the rotating color palette, or just fixed list of colors.
useCustomColor = False
ColorMode = 0
# 0 -> Predefined color set
# 1 -> Rotating palette

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

#  gStyle.SetPadBorderMode(0)
#  gStyle.SetFrameBorderMode(0)
  gStyle.SetOptStat(0)
  gStyle.SetOptFit(0)

#  gStyle.SetLabelSize(0.05,"X")
#  gStyle.SetLabelSize(0.05,"Y")
# 0.035
#  gStyle.SetLabelSize(0.035,"XY")
  gStyle.SetTitleSize(AxisLabelSizeX,"X")
  gStyle.SetTitleSize(AxisLabelSizeY,"Y")

  gStyle.SetTitleOffset(kDefaultTitleOffsetX,"X")
  gStyle.SetTitleOffset(kDefaultTitleOffsetY,"Y")

  gStyle.SetOptStat(0);
  gStyle.SetOptTitle(0);
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
  # removing TH3s for now
  if obj.InheritsFrom("TH3"):
    return 0
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


# produce TGraph Errors as ratio
# for now, not doing binomial errors
def DivideTGraphErrors(num,denom,name):
  nPoints = num.GetN()
  nSkipPoints = 0 # count of when we can't divide
  newTGraph = TGraphErrors(nPoints)
  newTGraph.SetName(name)
  for i in range(nPoints):
    x = num.GetX()[i]
    x_err = num.GetEX()[i]
    y1 = num.GetY()[i]
    y2 = denom.GetY()[i]
    y = 0.0
    y_err = 0.0
    if ((y1 != 0) and (y2 != 0)):
      # technically could make this work for y1 = 0
      y = y1 / y2
      y1_err = num.GetEY()[i]
      y2_err = denom.GetEY()[i]
      y_err = y * math.sqrt(math.pow(y1_err/y1,2) + math.pow(y2_err/y2,2))
      newTGraph.SetPoint(i-nSkipPoints,x,y)
      newTGraph.SetPointError(i-nSkipPoints,x_err,y_err)
    else:
      newTGraph.RemovePoint(i-nSkipPoints)
      nSkipPoints+=1
  newTGraph.SetLineColor(num.GetLineColor())
  newTGraph.SetMarkerColor(num.GetMarkerColor())
  newTGraph.SetMarkerStyle(num.GetMarkerStyle())
  
  return newTGraph


# Produce new hist where each point's error is the variance of the point from all the objects
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
    print(" CenterValueMode = %d" % (CenterValueMode))
    if (CenterValueMode == 1):
      primaryY = statistics.mean(listOfYValues)
      print(" mean = %f" % (primaryY))
      newHist.SetBinContent(i+1,primaryY)
      # if CVM=0, then the bin content is already set by clone
    newHist.SetBinError(i+1,stdDev)
  return newHist

# Produces a new hist where each point's error is the geometric mean of the statistical and systematic error
def ProduceTotalErrorHists(primaryHist,sysUncertObj):
  newHist = primaryHist.Clone("%s_TotalErr" % (primaryHist.GetName()))
 # newHist.SetTitle("")
  nPoints = primaryHist.GetNbinsX()
  for i in range(nPoints):
    statErr = primaryHist.GetBinError(i+1)
    sysErr  = sysUncertObj.GetBinError(i+1)
    newErr = math.sqrt(statErr*statErr + sysErr*sysErr)
    newHist.SetBinError(i+1,newErr)
  return newHist


def ProduceSystematicFromGraphs(graphs):
  nGraphs=len(graphs)
  if (nGraphs < 2):
    return 0
  primaryGraph = graphs[0]
  newGraph=primaryGraph.Clone("%s_SysErr" % (primaryGraph.GetName()))
  newGraph.SetTitle("")
  nPoints=primaryGraph.GetN()
  ListOfRangeHists=[]
  print("  Producing systematic uncertainty for object %s with %d points each" % (primaryGraph.GetName(),nPoints))
  for i in range(nPoints):
    #xCenter=primaryGraph.GetX()[i]
    xError=primaryGraph.GetEX()[i]
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
    if (CenterValueMode == 1):
      primaryY = statistics.mean(listOfYValues)
      print(" mean = %f" % (primaryY))
      newGraph.SetPoint(i,primaryGraph.GetX()[i],primaryY)

    stdDev=statistics.stdev(listOfYValues)
    nBinsRange=150
    minValue=min(listOfYValues)
    maxValue=max(listOfYValues)
    valueRange=maxValue-minValue
    localHistName="%sBin%dRange" % (primaryGraph.GetName(),i)
    localHistTitle="%s Bin %d Range" % (primaryGraph.GetName(),i)
    localHist=ROOT.TH1F(localHistName,localHistTitle,nBinsRange,minValue-0.15*valueRange,maxValue+0.15*valueRange);
    for value in listOfYValues:
      localHist.Fill(value)
    ListOfRangeHists.append(localHist)
    newGraph.SetPointError(i,xError,stdDev)#FIXME
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
  return (newGraph,ListOfRangeHists)
  #return newGraph

def ProduceTotalErrorGraphs(primaryGraph,sysUncertObj):
  newGraph = primaryGraph.Clone("%s_TotalErr" % (primaryGraph.GetName()))
  nPoints=primaryGraph.GetN()
  print("Producing the total error plots for %s with %d points." % (newGraph.GetName(),nPoints))
  for i in range(nPoints):
    primaryY = primaryGraph.GetY()[i]
    statYE   = primaryGraph.GetEY()[i]
    systYE   = sysUncertObj.GetEY()[i]
    totalE   = math.sqrt(statYE*statYE + systYE*systYE)
    newGraph.SetPointError(i,primaryGraph.GetEX()[i],totalE)
  return newGraph


#---------------------------------------------------------------------------------------------------
def sysCompare():

  parser = argparse.ArgumentParser(description='Produce Systematic Uncertainties From Varied inputs.')
  parser.add_argument('-d','--directory',required=False,default="",type=str,help="Directory for comparison plots. If none given, plots are not saved")
  parser.add_argument('-l','--list',required=False,type=str,default="",help="List of histograms/graphs to use. If not given, all objects are used.")
  parser.add_argument('-o','--output',required=True,type=str,help="Output file to produce")
#  parser.add_argument('-a','--all',requred=False,type=bool,help="whether to use all valid objects in the first file")
  parser.add_argument('-r','--ratioMode',required=False,type=bool,default=False,help="Whether to produce plots of the ratios")
  parser.add_argument('-c','--centerValueMode',required=False,default=0,type=int,help="Whether to use the first file as the central value (0) or the average (1)")
  parser.add_argument('-C','--ColorMode',type=int,default=1,help="Whether to use predefined colors or a rotating spectrum")
#  parser.add_argument('-c','--centerValueMode',required=False,default=0,type=int,help="Whether to use the first file as the central value (0) or the average (1)")
  parser.add_argument('-t','--titles',required=True,type=str,nargs='+',help="Titles for each file")
  parser.add_argument('-f','--files',metavar='Files',required=True,type=str,nargs='+',help="Files to use")

  parser.add_argument('-D','--DeletePoints',required=False,default=0,type=int,help="Delete the beginning N points for TGraphs")
#  parser.add_argument('files',metavar='Files',type=str,nargs='+',help="Files to use")

  parser.add_argument('-T','--OverallTitle',required=False,type=str,default="",help="Title to appear in legend headers")
  parser.add_argument('-y','--LogY',required=False,type=int,default=False,help="Whether to force the plots to have log Y.")


  args = parser.parse_args()

  setStyle()

  stringListOfHists=args.list
  #list of hists/graphs
  listOfHists=stringListOfHists.split()

#  bUseAllObjs=False
#  if (len(listOfHists) == 0):
#    bUseAllObjs=True
  bUseAllObjs = (len(listOfHists) == 0)

  bRatioMode=args.ratioMode
  if (bRatioMode):
    print("Ratio Mode is enabled.")

  #array of files
  fileNames=args.files
  fileTitles=args.titles

  directory=args.directory
  outputFileName=args.output

  global CenterValueMode
  CenterValueMode=args.centerValueMode
  LogYMode=False
  LogYMode=args.LogY

  global ColorMode
  ColorMode=args.ColorMode

  numDelete=args.DeletePoints

  print("Center Value Mode = %d" % (CenterValueMode))

  print("List of hists/graphs to use:"),
  print(listOfHists)

  if (bUseAllObjs):
    print("Will use all valid objects found in the first file")
  else:
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
  if (numDelete > 0):
    print("Will delete the first %d points from TGraphs" % (numDelete))

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
  outputFile=TFile.Open(outputFileName,"RECREATE")


  OverallTitle=args.OverallTitle



  canvas=TCanvas("canvas","canvas",c_width,c_height)
  if (directory != ""):
    # check if the directory exists
    if (not os.path.isdir(directory)):
      os.makedirs(directory)
    os.chdir(directory)

  # This is where it gets easier than compare: just need to find the objects in the files
  
  if (bUseAllObjs):
    print("Building list of objects")
    #primaryObj
    ListOfKeys=primaryFile.GetListOfKeys()
    for key in ListOfKeys:
      print("Adding item %s" % (key.GetName()))
      listOfHists.append(key.GetName())

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
    legWidth=0.2#0.45
    legHeight=0.4#0.225

    # legend for comparison
#    leg = TLegend(legX,legY,legX+legWidth,legY+legHeight)
    # legend for SysUncert
#    leg2 = TLegend(legX,legY,legX+legWidth,legY+legHeight)
    # legend for comparison
    leg = TLegend(legWidth,legHeight,legWidth,legHeight)

    if (OverallTitle != ""):
      leg.SetHeader(OverallTitle,"c")

    # legend for SysUncert
    leg2 = TLegend(legWidth,legHeight,legWidth,legHeight)

    # Set Properties and include in legend
    for i in range(len(listOfObjs)):
      localObj=listOfObjs[i]
      localTitle=fileTitles[i]
      color=ROOT.kBlack
      if (i != 0):
        if (ColorMode==1):
          color = GetCustomColor(i)
        else:
          color = colorList[i]
      if (i >= len(markerList)):
        # This check code is redundant for now, in case I 
        # want to switch so something fancier than just looping
        # marker styles
        markerStyle = markerList[i % len(markerList)]
      else:
        markerStyle = markerList[i]
      localObj.SetMarkerColor(color)
      localObj.SetLineColor(color)
      localObj.SetMarkerStyle(markerStyle)
      localObj.SetTitle(localTitle)

      localObj.SetFillColor(0)

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

    if (iObjType == 0): # Not coded yet
      print("Object is of a type I don't have comparison code for (yet)")
      continue
    if (iObjType == 1): # TGraph
      # Get object name
      #objName = "TestGraph"

      # Reset the margins
      gPad.SetTopMargin(fDefaultTopMargin)
      gPad.SetLeftMargin(fDefaultLeftMargin)
      gPad.SetBottomMargin(fDefaultBottomMargin)
      gPad.SetRightMargin(fDefaultRightMargin)

      mg = TMultiGraph()

      for j in range(numDelete):
        print("Deleting a point from object %s" % (primaryObj.GetName()))
        primaryObj.RemovePoint(0) # remove the first point

      mg.SetTitle(primaryObj.GetTitle())
      mg.GetXaxis().SetTitle(primaryObj.GetXaxis().GetTitle())
      mg.GetYaxis().SetTitle(primaryObj.GetYaxis().GetTitle())
      # FIXME is this going to double deleting points from the primary object
      for lobj in listOfObjs:
        print("Object starts with %d points" % (lobj.GetN()))
        for j in range(numDelete):
          print("Deleting a point from object %s" % (lobj.GetName()))
          lobj.RemovePoint(0) # remove the first point
        print("Object ends with %d points" % (lobj.GetN()))
        mg.Add(lobj)
      mg.Draw("ALP")
      mg.GetXaxis().SetLabelSize(AxisLabelSizeX)
      mg.GetYaxis().SetLabelSize(AxisLabelSizeY)


      # Draw a title for the TGraphs
      tp = TPaveText(TitlePaveMinX,TitlePaveMinY,TitlePaveMaxX,TitlePaveMaxY,"NDC")
      tp.SetLineStyle(TitlePaveLineStyle)
      tp.SetLineWidth(TitlePaveLineWidth)
      if (objectTitle==""):
        tp.AddText(objName)
        tp.Draw("SAME")
      elif (objectTitle=="Graph"):
        tp.AddText("test")
      else:
        tp.AddText(objectTitle)
        tp.Draw("SAME")


      #leg.Draw("SAME")
      legtest = gPad.BuildLegend(Legend1Width,Legend1Height,Legend1Width,Legend1Height)
      #legtest = gPad.BuildLegend()
      legtest.Draw("SAME")
      if (directory != ""):
        canvas.Print("%s_Cmp.pdf" % (objName))
        canvas.Print("%s_Cmp.png" % (objName))
        canvas.Print("%s_Cmp.C" % (objName))

      # Now produce systematic uncertainties (for TGraphErrors
      (sysUncertObj,ListOfRangeHists)=ProduceSystematicFromGraphs(listOfObjs)

      sysUncertObj.SetName("%s_SysErr" % (objName))
      #sysUncertObj.SetName(objName)

      
      totalUncertObj=ProduceTotalErrorGraphs(primaryObj,sysUncertObj)
      totalUncertObj.SetName(objName)
      #totalUncertObj.SetName("%s_TotalErr" % (objName))

      

      sysUncertObj.SetTitle("Systematic Uncertainty")
      sysUncertObj.SetFillColor(ROOT.kBlue)
      sysUncertObj.SetFillStyle(3002)
      sysUncertObj.SetMarkerColor(ROOT.kBlue-9)
      sysUncertObj.Draw("ALP[]5")

      # reset the primary objects title
      primaryObj.Draw("LP")
      if (CenterValueMode == 0):
        leg2.AddEntry(sysUncertObj,"Systematic Uncertainty","F")
      else:
        leg2.AddEntry(sysUncertObj,"Central Value + Systematic Uncertainty","PF")

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
      outputFile.Add(sysUncertObj)
      outputFile.Add(totalUncertObj)
      for hist in ListOfRangeHists:
        outputFile.Add(hist)

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
        canvas.Print("%s_SysUncert_Cmp.C" % (objName))

      canvas.Clear()
      canvas.Divide(1,2,0.01,0.0)
      ratioMg = TMultiGraph()
      legRatio=TLegend(2*legWidth,legHeight,2*legWidth,legHeight)
      #RatioArray=[]
      # Make, draw ratios
      if (bRatioMode):
        max_ratio_y=1.
        min_ratio_y=1.
        nPointsFirst=primaryObj.GetN()
        bDivisionPossible=True
        for lobj in listOfObjs:
          if (nPointsFirst != lobj.GetN()):
            bDivisionPossible=False
        if(bDivisionPossible):
          numObjects=0
          for lobj in listOfObjs:
            numObjects = numObjects + 1
            if (lobj == primaryObj):
              print("Avoided TGraph over itself using object")
              continue
            if (lobj.GetName() == primaryObj.GetName()):
              print("Avoided TGraph over itself using name")
              continue
            if (numObjects == 1):
              print("Avoided TGraph over itself using index")
              continue
            ratioName = "%s_Ratio" % (lobj.GetName())
            lRatio = DivideTGraphErrors(lobj,primaryObj,ratioName)
            lRatio.SetTitle("%s / %s" % (lobj.GetTitle(),fileTitles[0]))
            ratioMg.Add(lRatio)
            legRatio.AddEntry(lRatio,lRatio.GetTitle(),"LP")
           # RatioArray.append(lRatio)
  
        ratioMg.GetYaxis().SetTitle("Ratio over (%s)" % (fileTitles[0]))
  
        canvas.cd(1)
        mg.Draw("ALP X+")
        gPad.SetTopMargin(0.0)
        leg.Draw("SAME")
        gPad.SetBottomMargin(0.0)
        canvas.cd(2)
        ratioMg.Draw("ALP")
        legRatio.Draw("SAME")
        if (directory != ""):
          canvas.Print("%s_Ratio.pdf" % (objName))
          canvas.Print("%s_Ratio.png" % (objName))
          canvas.Print("%s_Ratio.C" % (objName))

    # for histograms, also draw a plot with each of them
    # separately? Useful if the fit functions are visible
    if (iObjType == 2): # TH1


      # Reset the margins
      gPad.SetTopMargin(fDefaultTopMargin)
      gPad.SetLeftMargin(fDefaultLeftMargin)
      gPad.SetBottomMargin(fDefaultBottomMargin)
      gPad.SetRightMargin(fDefaultRightMargin)




      primaryObj.Draw()
 #     if (LogYMode):
 #       canvas.SetLogy(1)
      primaryObj.GetYaxis().UnZoom()
      primaryObj.GetXaxis().SetLabelSize(AxisLabelSizeX)
      primaryObj.GetYaxis().SetLabelSize(AxisLabelSizeY)

      primaryObj.GetXaxis().SetTitleSize(kDefaultTitleSizeX)
      primaryObj.GetYaxis().SetTitleSize(kDefaultTitleSizeY)

      primaryObj.GetXaxis().SetTitleOffset(kDefaultTitleOffsetX)
      primaryObj.GetYaxis().SetTitleOffset(kDefaultTitleOffsetY)

      primaryObj.SetMarkerSize(kMarkerSize)

      fYMin = GetMinValue(primaryObj)
      fYMax = GetMaxValue(primaryObj)
      for lobj in listOfObjs:
        if (lobj != primaryObj):
          lobj.SetMarkerSize(kMarkerSize)
          # h.GetFunction(“myFunction”)->SetBit(TF1::kNotDraw);

          lobj.Draw("SAME")
          fYMin = min(fYMin,GetMinValue(lobj))
          fYMax = max(fYMax,GetMaxValue(lobj))
      
      if (LogYMode != 0):
        (fYMin,fYMax) = ExpandRange(fYMin,fYMax)

      #leg.Draw("SAME")

      legHist1=gPad.BuildLegend()
      if (OverallTitle != ""):
        legHist1.SetHeader(OverallTitle,"c")
      #legHist1.RemoveEntry(0)
      #legHist.ListOfPrimitives()
      legHist1.Draw("SAME")

      primaryObj.GetYaxis().SetRangeUser(fYMin,fYMax)
      if (LogYMode == 1):
        canvas.SetLogy(1)
      primaryObj.GetYaxis().UnZoom()
      # FIXME temp
     # primaryObj.GetXaxis().SetRangeUser(0,25)


      # Draw a title for the TH1's
      tp = TPaveText(TitlePaveMinX,TitlePaveMinY,TitlePaveMaxX,TitlePaveMaxY,"NDC")
      tp.SetLineStyle(TitlePaveLineStyle)
      tp.SetLineWidth(TitlePaveLineWidth)
      if (objectTitle==""):
        tp.AddText(objName)
      else:
        tp.AddText(objectTitle)
      tp.Draw("SAME")

      if (directory != ""):
        canvas.Print("%s_Cmp.pdf" % (objName))
        canvas.Print("%s_Cmp.png" % (objName))
        canvas.Print("%s_Cmp.C" % (objName))
      sysUncertObj=ProduceSystematicFromHists(listOfObjs)
      #sysUncertObj.SetName(objName) # This sets the name to be that of the original object
      sysUncertObj.SetName("%s_SysErr" % (objName)) # This adds the label SysErr

      totalUncertObj=ProduceTotalErrorHists(primaryObj,sysUncertObj)
      totalUncertObj.SetName(objName)
      #totalUncertObj.SetName("%s_TotalErr" % (objName))

      sysUncertObj.SetFillColor(ROOT.kBlue)
      sysUncertObj.SetFillStyle(3002)
      sysUncertObj.SetMarkerColor(ROOT.kBlue-9)
      #sysUncertObj.Draw("E2")
      sysUncertObj.Draw("E4")
      primaryObj.Draw("SAME")

      if (CenterValueMode == 0):
        leg2.AddEntry(sysUncertObj,"Systematic Uncertainty","F")
      else:
        leg2.AddEntry(sysUncertObj,"Central Value + Systematic Uncertainty","PF")
      leg2.Draw("SAME")

      if (directory != ""):
        canvas.Print("%s_SysUncert.pdf" % (objName))
        canvas.Print("%s_SysUncert.png" % (objName))
      outputFile.Add(sysUncertObj)
      outputFile.Add(totalUncertObj)

      if (bRatioMode):
        max_ratio_y=1.
        min_ratio_y=1.
        nBinsFirst = primaryObj.GetNbinsX()
        bDivisionPossible=False
        for hobj in listOfObjs:
          if (nBinsFirst == hobj.GetNbinsX()):
            bDivisionPossible=True
        if (bDivisionPossible):
          RatioArray=[]

        #canvas.Divide(1,2,canvSmall,canvSmall)
        #canvas.cd(1)
        #gPad.SetBottomMargin(small)
        #primaryObj.Draw()
        # Build the ratios.
        # Draw the thing
        for lobj in listOfObjs:
          lRatio = lobj.Clone("%s_Ratio" % (lobj.GetName())) 
          if (bBinomialDivision):
            lRatio.Divide(lRatio,primaryObj,1.0,1.0,"B")
          else:
            lRatio.Divide(primaryObj)
          #lRatio.GetFunction(“myFunction”).SetBit(TF1::kNotDraw)
          localRatioMin=lRatio.GetBinContent(lRatio.GetMinimumBin())
          localRatioMax=lRatio.GetBinContent(lRatio.GetMaximumBin())
          if (localRatioMin < min_ratio_y):
            min_ratio_y = localRatioMin
          if (localRatioMax > max_ratio_y):
            max_ratio_y = localRatioMax
          RatioArray.append(lRatio)
          
        # magic adjustments to ratio min/max
        if (min_ratio_y != 0):
          if (max_ratio_y > 1./min_ratio_y):
            min_ratio_y = 1./max_ratio_y
          else:
            max_ratio_y = 1./min_ratio_y
        if (max_ratio_y > 20.):
          max_ratio_y=1.3
        min_ratio_y=1.0/1.3
        #min_ratio_y = pow(min_ratio_y,1.6)
        #max_ratio_y = pow(max_ratio_y,1.6)


        legRatio=TLegend(legWidth,legHeight,legWidth,legHeight)
        RatioArray[0].Draw("HIST E")
        RatioArray[0].GetYaxis().SetRangeUser(min_ratio_y,max_ratio_y)

        for lRatio in RatioArray:
          lRatio.Draw("SAME HIST E")
          legRatio.AddEntry(lRatio,lRatio.GetTitle(),"LP")
          
        legRatio.Draw("SAME")
        if (directory != ""):
          canvas.Print("%s_Ratio.pdf" % (objName))
          canvas.Print("%s_Ratio.png" % (objName))
          canvas.Print("%s_Ratio.C" % (objName))


#    primaryObj.Draw()
#    for lobj in listOfObjs:
#      lobj.Draw("SAME")
#    canvas.Print("Test.pdf")

  outputFile.Write()


#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  sysCompare()



exit()

