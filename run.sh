#!/bin/bash 

# https://yutarine.blogspot.com/2018/02/shell-command-bc02.html
# [コマンド] bc(数値計算ソフト)の結果出力で、小数点前の0(ゼロ)を省略させないようにする方法

#PAR="mu"
#PAR="el"
#PAR="gm"
PAR="pr"

OUTROOT="p0.5GeV_${PAR}_water_BisY"

# ./wavy run.mac *.root
# 100cm x 100cm x 30layears is ok w. sx -n 5
#bsub -q sx -n 5 -J ${PAR} "./wavy  g4mac_to_z.mac  gps_${OUTROOT}_toZ.root  > log_${OUTROOT}_toZ.log 2>&1"
#bsub -q sx -n 5 -J ${PAR} "./wavy  g4mac_to_y.mac  gps_${OUTROOT}_toY.root  > log_${OUTROOT}_toY.log 2>&1"
bsub -q sx -n 5 -J g4_${PAR} "./wavy  g4mac_to_x.mac  gps_${OUTROOT}_toX.root  > log_${OUTROOT}_toX.log 2>&1"
bsub -q sx -n 5 -J g4_${PAR} "./wavy  g4mac_to_y.mac  gps_${OUTROOT}_toY.root  > log_${OUTROOT}_toY.log 2>&1"
bsub -q sx -n 5 -J g4_${PAR} "./wavy  g4mac_to_z.mac  gps_${OUTROOT}_toZ.root  > log_${OUTROOT}_toZ.log 2>&1"


<<KILL
DIR="output"

mkdir ./$DIR
#./wavy vis.mac name
rm -rf ./$DIR/*

NPINTS=1000;
for ((itX=0 ; itX<$NPINTS  ; itX++)); do

	MACRO="./${output}/g4mac_RUN${itX}.mac"
   cp ./run.mac $MACRO
   sed -i -e "s;###SEED###;${itX};g" $MACRO
	
   ROOTNAME="./${DIR}/data_RUN${itX}.root"
   LOGDILE="./${DIR}/log_${itX}.log" 

	bsub -q s -J r${itX} "./wavy $MACRO $ROOTNAME ${itX} > $LOGDILE 2>&1"
done
KILL

#hadd ../output_prot/sum.root ../output_prot/*root
#hadd ../output_elec/sum.root ../output_elec/*root
#hadd ../output_muon/sum.root ../output_muon/*root

## 10 mm cube 
#NPINTS=30;
#OFFSET=0.02 # cm
#STEP=`echo "scale=5; (1 - 2*$OFFSET) / $NPINTS" | bc`
#initX=`echo "scale=5; - 1/2 + $OFFSET" | bc | sed 's/^\./0./'`
#initY=`echo "scale=5; - 1/2 + $OFFSET" | bc | sed 's/^\./0./'`
#
#echo "initX=$initX, initY=$initY, STEP=$STEP"
#
#for ((itX=0 ; itX<$NPINTS  ; itX++)); do
#   #   if [ $mod -eq 1 ] || [ $mod -eq 2 ] || [ $mod -eq 4 ] || [ $mod -eq 6 ]; then
#   #      continue
#   #   fi
# 	for ((itY=0 ; itY<$NPINTS ; itY++)); do
#   	#if [ ! -d "0${runs[i]}" ] ; then
#      #fi;
#		xposi=`echo "scale=2; $initX + $STEP*$itX" | bc | sed 's/^\./0./'`		
#		yposi=`echo "scale=2; $initY + $STEP*$itY" | bc | sed 's/^\./0./'`		
#		echo "x=${xposi}, y=${yposi}"
#   	#MACRO="../output/g4_macro_X${xposi}_Y${yposi}.mac"
#		MACRO="../output/g4mac_X${itX}_Y${itY}.mac"
#
#   	#cp ./run.mac ../output/g4_macro_X${xposi}_Y${yposi}.mac 
#   	cp ./run.mac $MACRO 
#   	sed -i -e "s;###XPOSI###;${xposi};g" $MACRO
#   	sed -i -e "s;###YPOSI###;${yposi};g" $MACRO
#
#		ROOTNAME="../output/root_X${itX}_Y${itY}.root"
#		echo "$MACRO"
#      echo "$ROOTNAME"
#      bsub -q s -J X${itX}_Y${itY} "./wavy $MACRO $ROOTNAME > ../output/log_X${itX}_Y${itY}.log 2>&1"
#	done
#done
#
 
