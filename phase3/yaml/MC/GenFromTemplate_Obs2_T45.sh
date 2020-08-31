#!/bin/bash

templateFile="config_SideBand_Obs2_Template.yaml"


# 0 = Pur=0 (so, no correction applied)
# 1 = Primary (the one from the phase 1 file)
# 2 = 1
# 3 = Purity - Error
# 4 = Purity + Error
PURITYCHOICES=(0 1 2 3 4)

#CENTBINS=(0 1 2 3)
CENTBINS=(0)

EPBIN=-1
#EPBINS=(-1 0 1 2)

MCMODES=(0 1 2)

PTBIN=1

PTLOW=5.0
PTHIGH=7.0

#PTLOW=7.0
#PTHIGH=9.0

#PTLOW=11.0
#PTHIGH=14.0

#PTLOW=14.0
#PTHIGH=17.0

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
	for PUR_CHOICE in ${PURITYCHOICES[@]}
	do
		for MCMODE in ${MCMODES[@]}
		do
			NewFileName=Obs2/PtBin$PTBIN/config_P2_T53_McMode$MCMODE""_PurChoice$PUR_CHOICE""_Cent$CENTBIN""_EP$EPBIN"".yaml
			echo "Creating file $NewFileName"
			cat $templateFile | sed -e "s/\$MCMODE/$MCMODE/g" -e "s/\$PURITYCHOICE/$PUR_CHOICE/g" -e "s/\$PTBIN/$PTBIN/g" -e "s/\$CENTBIN/$CENTBIN/g" -e "s/\$EPBIN/$EPBIN/g" -e "s/\$PTHIGH/$PTHIGH/g" -e "s/\$PTLOW/$PTLOW/g"  -e "s/\$CENTLABEL/$CENTLABEL/g" > $NewFileName
		done
	done
done

