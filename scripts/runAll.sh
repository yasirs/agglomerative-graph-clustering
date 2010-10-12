#!/bin/bash
for m in 0.05 0.1 0.2 0.3 0.4 0.5
do
	i=0
	while test $i < 100
	do
		tstamp = './holdout $1 $m'
		thisdir=~/work/agglo/outs/$1\_$m\_$tstamp
		mkdir $thisdir
		./holdout $1 $m
		./get_scores $1
		./writeLabels $1.scores0 $1.gold $1.labels
		cp $1.scores* $1.gold $1.labels $1.jaccard $1.dprod $1.cneighb $1.hyperg $1.edges $thisdir
		R --no-save F=$1 M=$m < addPerfs.r
		i=`expr $i + 1`
	done
done





