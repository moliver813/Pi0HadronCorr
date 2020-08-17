#!/bin/bash

echo "The purpose of this script is to take a yaml file, and produce several yaml files that very the fit function used"

if [ "$#" -eq 0 ]; then
	exit 1
fi

INPUTYAML=$1
LABEL=$2
echo "Will start with yaml file $INPUTYAML under label $LABEL"


# put these output yaml files in a new dir, with intermediate
# dirs matching the input yaml path
#OUTDIR

# Fit Ranges
FITRANGEMAXS=(0 0.35 0.4 0.45 0.5 0.6 0.65 0.7 0.75 0.8)
# Peak functions
PEAKNAMES=("Gaus" "GausExpLeft" "Brent-Wigner" "CrystalBallLeft" "CrystalBallRight" "GausExpRight" "ExpGausExp" "Voigt")
BKGNAMES=("Poly(0)"  "Poly(1)" "Poly(2)" "Poly(3)" "Poly(4)" "Exp*Poly(2)"  "Exp*Poly(3)" "ExpDecay*Poly(2)" "ExpDecay*Poly(3)")

#PEAKFUNCTIONS=(0 1 2 3 4 5 6 7)
#BKGFUNCTIONS=(0 1 2 3 4 5 6 7)
PEAKFUNCTIONS=(0 1 2 3 4 5 6 7)
BKGFUNCTIONS=(0 1 2 3 4 5 6 7 8)
# 0 is the default, each i refers to a different value
FITRANGEOPTS=(0 1 2 3 4 5 6 7 8 9)

#             0 Gaussian
#             1 ExpDecay/Gaussian
#             2 Brent-Wigner
#             3 Crystal Ball
#             4 Crystal Ball (Right Side)
#             5 ExpDecay/Gaussian (Right Side)
#             6 ExpDecay/Gaussian (Both Sides)
#             7 Voigt Profile

#// Choice of Residual Fit Background

#             0-4 Poly(n)
#             5 Exp * Poly(2)
#             6 Exp * Poly(3)
#             7 ExpDecay * Poly(2)
#             8 ExpDecay * Poly(3)


for PEAKF in ${PEAKFUNCTIONS[@]}
do
	for BKGF in ${BKGFUNCTIONS[@]}
	do
		for FITOPT in ${FITRANGEOPTS[@]}
		do
			FITMAX=${FITRANGEMAXS[$FITOPT]}
			PEAKNAME=${PEAKNAMES[$PEAKF]}
			echo ""
			echo "Peak function $PEAKF (name = $PEAKNAME), background function $BKGF, peak opt $FITOPT (value = $FITMAX)"
			echo "  input file $INPUTYAML"
	#		OUTPUTFILE=Func$PEAKF$BKGF/`echo $INPUTYAML | sed -e "s/\.yaml/_Peak$PEAKF""_Bkg$BKGF\.yaml/"`
	#		OUTDIR=FuncFunc$PEAKF$BKGF/`echo $INPUTYAML | sed -e "s/\/.*\.yaml//g"`
			OUTDIR=FuncScan/$LABEL/Func$PEAKF$BKGF$FITOPT/`ls $INPUTYAML | sed -e "s/config.*\.yaml//g"`
			OUTPUTFILE=FuncScan/$LABEL/Func$PEAKF$BKGF$FITOPT/`ls $INPUTYAML | sed -e "s/\.yaml/_Peak$PEAKF""_Bkg$BKGF""_FitRange$FITOPT\.yaml/"`
			echo "  output directory  $OUTDIR"
			echo "  output file $OUTPUTFILE"
			mkdir -p $OUTDIR/output
			OUTDIRWITHQUOTES=\"$OUTDIR/output/\"
			cat $INPUTYAML | sed -e "s@output_dir.*@output_dir: \"$OUTDIR/output\"@g" -e "s/FitPeakMethod.*/FitPeakMethod: $PEAKF/" -e "s/FitBkgMethod.*/FitBkgMethod: $BKGF/" > $OUTPUTFILE
			echo "fitMinX: 0" >> $OUTPUTFILE
			echo "fitMaxX: $FITMAX" >> $OUTPUTFILE
			echo "Label2: \"Peak $PEAKF ($PEAKNAME), Bkg $BKGF, PeakFitOpt $FITOPT\"" >> $OUTPUTFILE
		done
	done
done



# How to deal with the output. Rewrite the outputdir


#BkgChoice: 0
# bkgChoice = 0 No background
#             1 Mixed Event background (Free Scale)
#             2 Rotational background (Free Scale)
#             3 Mixed Event background (Fixed Scale) 
#             4 Rotational background (Fixed Scale)
#             5 Pos. Swap background (Free Scale)
#             6 Pos. Swap background (Fixed Scale)


#FitPeakMethod: 6
#             0 Gaussian
#             1 ExpDecay/Gaussian
#             2 Brent-Wigner
#             3 Crystal Ball
#             4 Crystal Ball (Right Side)
#             5 ExpDecay/Gaussian (Right Side)
#             6 ExpDecay/Gaussian (Both Sides)
#             7 Voigt Profile

#FitBkgMethod: 4
#             0-4 Poly(n)
#             4 Exp * Poly(2)
#             5 Exp * Poly(3)
#             6 ExpDecay * Poly(2)
#             7 ExpDecay * Poly(3)

