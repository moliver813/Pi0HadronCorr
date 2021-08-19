#!/bin/bash

templateFile="config_SideBand_Corr4EGA_Obs2_Template.yaml"

#SBBkgSels=(0 1 2 3 4 5 6)
SBBkgSels=(0 1 2 3)

CENTBINS=(0 1 2 3)

EPBINS=(-1 0 1 2)

# 0 = Pur=0 (so, no correction applied)
# 1 = Primary (the one from the phase 1 file)
# 2 = 1
# 3 = Purity - Error
# 4 = Purity + Error

PURCHOICES=(1 3 4)

# 0 = No scaling
# 1 = LinearScaling
# 2 = Quadratic Scaling
# 3 = Fractional Scaling
# 4 = Average of LinearScaling and first sideband
FITFUNCS=(0 1 2 3 4)

# Acceptable combinations of SBSel and fit functions
# 0 : 0 1 2 3 4 5 6
# 1 : 0 1 2 3
# 2 : 0
# 3 : 0
# 4 : 0 1 2 3

#PTBIN=4

#PTLOW=5.0
#PTHIGH=7.0
#PTLOW=7.0
#PTHIGH=9.0
#PTLOW=9.0
#PTHIGH=11.0
#PTLOW=11.0
#PTHIGH=14.0
#PTLOW=14.0
#PTHIGH=17.0

PTLOWS=(0 5 7 9 11 14)
PTHIGHS=(0 7 9 11 14 17)


PTBINS=(3 4 5)

for PTBIN in ${PTBINS[@]}
do

	PTLOW=${PTLOWS[$PTBIN]}
	PTHIGH=${PTHIGHS[$PTBIN]}

	mkdir -p Obs2/PtBin$PTBIN

	for CENTBIN in ${CENTBINS[@]}
	do
		CENTLABEL=$CENTBIN
		if [ $CENTBIN = 0 ]
		then
			CENTLABEL="0-10%"
		elif [ $CENTBIN = 1 ]
		then
			CENTLABEL="10-30%"
		elif [ $CENTBIN = 2 ]
		then
			CENTLABEL="30-50%"
		elif [ $CENTBIN = 3 ]
		then
			CENTLABEL="50-80%"
		fi

		for PURCHOICE in ${PURCHOICES[@]}
		do
			mkdir -p Obs2/PtBin$PTBIN/Cent$CENTBIN/PurChoice$PURCHOICE
			for FITFUNC in ${FITFUNCS[@]}
			do
				# Filtering out unacceptable combinations
				if [ $FITFUNC = 0 ]
				then
					SBBkgSels=(0 1 2 3 4 5 6)
				elif [ $FITFUNC = 1 ]
				then
					SBBkgSels=(0 1 2 3)
				elif [ $FITFUNC = 2 ]
				then
					SBBkgSels=(0)
				elif [ $FITFUNC = 3 ]
				then
					SBBkgSels=(0)
				elif [ $FITFUNC = 4 ]
				then
					SBBkgSels=(0 1 2 3)
				fi

				for BkgSel in ${SBBkgSels[@]}
				do

					for EPBIN in ${EPBINS[@]}
					do
				#		NewFileName=Obs2/PtBin$PTBIN/config_P2_Corr4EGA_Cent$CENTBIN""_EP$EPBIN""_PurChoice$PURCHOICE.yaml
						NewFileName=Obs2/PtBin$PTBIN/Cent$CENTBIN/PurChoice$PURCHOICE/config_P2_Corr4EGA_SBSel$BkgSel""_FitFunc$FITFUNC""_PtBin$PTBIN""_Cent$CENTBIN""_EP$EPBIN""_PurChoice$PURCHOICE.yaml
		#				echo "Creating file $NewFileName"
		#				cat $templateFile | sed -e "s/\$PTBIN/$PTBIN/g" -e "s/\$CENTBIN/$CENTBIN/g" -e "s/\$EPBIN/$EPBIN/g" -e "s/\$PTHIGH/$PTHIGH/g" -e "s/\$PTLOW/$PTLOW/g"  -e "s/\$CENTLABEL/$CENTLABEL/g" -e "s/\$PURITYCHOICE/$PURCHOICE/g" > $NewFileName
						echo "Creating file $NewFileName"
						cat $templateFile | sed -e "s/\$BKGSEL/$BkgSel/g" -e "s/\$FITFUNCTION/$FITFUNC/g" -e "s/\$PTBIN/$PTBIN/g" -e "s/\$CENTBIN/$CENTBIN/g" -e "s/\$EPBIN/$EPBIN/g" -e "s/\$PTHIGH/$PTHIGH/g" -e "s/\$PTLOW/$PTLOW/g"  -e "s/\$CENTLABEL/$CENTLABEL/g" -e "s/\$PURITYCHOICE/$PURCHOICE/g" > $NewFileName
					done
				done
			done
		done
	done
done
