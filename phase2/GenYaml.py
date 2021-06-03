#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Produce yaml file for Phase 2.')
parser.add_argument('label',metavar='Label',type=str,help='Label for Analysis')
parser.add_argument('-o','--observable',required=False,type=int,help='0 for Trigger Pt, 1 for z_t, 2 for Xi')
parser.add_argument('-c','--centrality',required=False,type=int,help='-1 for all, 0 [0-10], 1 [10-30], 2 [30-50], 3 [50-90]')
parser.add_argument('-e','--eventplane',required=False,type=int,default=-1,help='-1 for all, 0 [in], 1 [mid], 2 [out]')
parser.add_argument('-p','--ptbin',required=False,type=int,default=4,help='Trigger Pt Bin. {5-7,7-9,9-11,11-14,14-17')
parser.add_argument('-s','--sefile',required=False,type=str,default='',help='file with Same Event Data')
parser.add_argument('-m','--mefile',required=False,type=str,default='',help='file with Same Event Data')
parser.add_argument('-v','--seWagon',required=False,type=str,default='AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_T09_C07_histos',help='Wagon Name for the Same Event Task')
parser.add_argument('-w','--meWagon',required=False,type=str,default='AliAnalysisTaskGAMT_Pi0H_ME_tracks_caloClusters_T09_C07_histos',help='Wagon Name for the Mixed Event Task')
parser.add_argument('-x','--OldMix',nargs='?',const=True,required=False,type=bool,help='Whether to use old ME task name')

args = parser.parse_args()
#print args.accumulate(args.

Label = args.label
SE_file = "t_23_GA_ST_Corr/AnalysisResults.root"
ME_file = "t_24_GAMB_MT_Corr/merge/test_merge_4_files/AnalysisResults.root"
if args.sefile != '':
  print("# using SE file %s" % args.sefile)
  SE_file = args.sefile
else:
  print ("# using default SE file")
if args.mefile != '':
  print ("# using ME file %s" % args.mefile)
  ME_file = args.mefile
else:
  print("# using default ME file")

PtMin=args.ptbin
PtMax=args.ptbin

Observable = args.observable
Cent       = args.centrality # -1 for all , 0 [0-10], 1 [10-30], 2 [30-50], 3 [50-90]
EventPlane = args.eventplane # -for all, 0 [in], 1 [mid], 2 [out]

#Observable = 1  # 0 for Trigger Pt, 1 for z_t, 2 for Xi
#Cent       = 2 # -1 for all , 0 [0-10], 1 [10-30], 2 [30-50], 3 [50-90]
#EventPlane = 0 # -for all, 0 [in], 1 [mid], 2 [out]

#Label += "_PtBin_%d_%d" % (PtMin,PtMax)
Label += "_PtBin%d_Cent%d_EP%d" % (PtMin,Cent,EventPlane)



# The Rest

output_txt = "## Script Generated yaml file\n"
output_txt += "Train_histo: 0\n\n"

# Analysis Settings
output_txt += "#Analysis Settings\n"

# Whether it is using gammas, pi0s or Sidebands
output_txt += "gamma_Pi_SB: 1\n"


output_txt += "savePlots: 1\n"   #..do you want to save the plots to a directory "output/"?
output_txt += "observable: %d\n" % Observable  #..do you want to plot 2DCorr for Ga,Zt, or Xi bins?
output_txt += "cent: %d\n" % Cent      #..do you want to plot all (-1), central (1),.....
output_txt += "evtPl: %d\n" % EventPlane  #..do you want to plot all (-1), in plane (0), mid plane (1), out of plane (2) .....
output_txt += "mergeME: 0\n"     #..merge ME plots for all bins to have more statistic
output_txt += "normMode: 0\n"    #..normalize SE after dividing by ME and adding up z-vertex bins (0), or before both (1)
output_txt += "plotAdvProj: 1\n" #..plot advanced histograms such as projections
output_txt += "plotComp: 1\n"    #..plot a comparision of 2 diffent ME correction methods
output_txt += "useFLGB: 0\n"     #..use FindLastGoodBin method to set delta eta limits, or use default delta eta limit
output_txt += "ptMinBin: %d\n" % PtMin  #..MinPt Bin (for Xi,Zt)
output_txt += "ptMaxBin: %d\n" % PtMax  #..MaxPt Bin (for Xi,Zt) # [5-22]
 #ptMaxBin: -1    #..MaxPt Bin (for Xi,Zt)

##Name to the path were the train output files are located - will not change
#output_txt += "file_path: \"/home/moliver/cern/gammaHadron/wrk/myTrains/\"\n"
# Not needed anymore
output_txt += "file_path: \"\"\n"

##Name of the Train output files for the analysis
output_txt += "file_nameSE: \"%s\"\n" % SE_file
output_txt += "file_nameME: \"%s\"\n" % ME_file

## Wagon Names
#output_txt += "wagon_NameSE: \"AliAnalysisTaskGASEv_Pi0H_SE_tracks_caloClusters_T09_C07_histos\"\n"
output_txt += "wagon_NameSE: \"%s\"\n" % args.seWagon
if (args.OldMix):
  output_txt += "wagon_NameME: \"AliAnalysisTaskGAMEv_Pi0H_ME_tracks_caloClusters_T09_C07_histos\"\n" # Mix Event
else:
  output_txt += "wagon_NameME: \"%s\"\n" % args.meWagon
#  output_txt += "wagon_NameME: \"AliAnalysisTaskGAMT_Pi0H_ME_tracks_caloClusters_T09_C07_histos\"\n" # Mix Trigger

## Label for projection root files, if created
output_txt += "label: \"%s\"\n" % Label
#output_txt += "label: \"%s_Observable%d\"\n" % (Label,Observable)
#output_txt += "label: \"Pi0_t21_PtBin_1\"\n"

## Root file with already projected histograms
output_txt += "root_histos: \"./Projections/%s_Observable%d_Cent%d_EvtPlane%d_Histograms.root\"\n" % (Label,Observable,Cent,EventPlane)
#output_txt += "root_histos: \"./%s_Observable%d_Cent%d_EvtPlane%d_Histograms.root\"\n" % (Label,Observable,Cent,EventPlane)
#output_txt += "root_histos: \"./Pi0_t21_PtBin_1_Observable1_Cent-1_EvtPlane-1_Histograms.root\"\n"

## Directory for storing output plots, output root files
output_txt += "output_dir: \"output/%s_Observable%d\"\n" % (Label,Observable)
#output_txt += "output_dir: \"output/Pi0_t21_PtBin_1\"\n"

print(output_txt)

