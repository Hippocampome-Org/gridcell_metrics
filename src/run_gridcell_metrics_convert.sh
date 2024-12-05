#!/bin/bash
less ./run_gridcell_metrics_with_config.sh > ./run_gridcell_metrics_temp1.m &&
sed 's/"//g' ./run_gridcell_metrics_temp1.m > ./run_gridcell_metrics_temp2.m &&
sed "s/'//g" ./run_gridcell_metrics_temp2.m > ./run_gridcell_metrics_temp3.m &&
sed "s/^..//g" ./run_gridcell_metrics_temp3.m > ./run_gridcell_metrics_temp4.m &&
less ./run_gridcell_metrics_temp4.m > ./run_gridcell_metrics.m &&
rm ./run_gridcell_metrics_temp1.m &&
rm ./run_gridcell_metrics_temp2.m &&
rm ./run_gridcell_metrics_temp3.m &&
rm ./run_gridcell_metrics_temp4.m
