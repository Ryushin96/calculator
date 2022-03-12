#!/bin/sh
		source /opt/intel/bin/iccvars.sh intel64
		sleep 1
		
		icc -O3 -o 400jitter 400jitter.c
		sleep 15
		
		ulimit -s unlimited
		sleep 1

		nohup ./400jitter &
		sleep 1

done
