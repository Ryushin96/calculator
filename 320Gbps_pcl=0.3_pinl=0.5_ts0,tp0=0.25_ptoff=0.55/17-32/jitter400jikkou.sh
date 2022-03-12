#!/bin/sh
	start_time=`date "+%Y-%m-%d %H:%M:%S"`
	echo $start_time
	
		source /home/fukuchi/intel/bin/iccvars.sh intel64
		sleep 1
		
		icc -O3 -o 400jitter 400jitter.c
		sleep 15
		
		ulimit -s unlimited
		sleep 1

		nohup ./400jitter
		sleep 1
	
	end_time=`date "+%Y-%m-%d %H:%M:%S"`
	echo $end_time

done
