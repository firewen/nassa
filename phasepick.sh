#!/bin/bash

PhaseNetdir=/home/wenj/tmpsoft/PhaseNet-master

eventdir=$1
cd $eventdir

saclst knetwk kstnm f data/*Z.sac | awk '{print $2"."$3}' | while read one
do
	fnm=$(ls data/$one.*.sac)
	#echo $one $fnm
	echo $one $fnm | awk 'BEGIN{print "fname E N Z"}{print $1,$2,$3,$4}' > sac.csv
	
	python $PhaseNetdir/phasenet/predict.py --model=$PhaseNetdir/model/190703-214543 --data_list=sac.csv --data_dir=./ --format=sac --plot_figure --save_prob
        outnm=$(echo $one | awk '{print $1}')
        mv results/result.h5 probs/$outnm.h5
        mv results/picks.csv probs/$outnm.csv
done

cd ..

