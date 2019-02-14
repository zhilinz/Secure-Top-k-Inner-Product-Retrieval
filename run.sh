#!/usr/bin/env bash



#rm -rf build
#mkdir build
#cd build
#cmake ../
#make

declare -a  dataset_list=("./data/MovieLens/" "./data/Yelp/")
#declare -a dataset_list=("./data/MovieLens/" "./data/Yelp/" "./data/Netflix/")
numOfDataSets=${#dataset_list[@]}

k_list="1"

maxip="50000 100000 200000 400000 800000"



bits="1024"

################### SKIP ####################
for (( index=0; index < numOfDataSets; index++ )); do
	for ip in $maxip
 	do
    	./bin/Release/MatrixGen SKIP $k_list ${dataset_list[$index]} $ip $bits
	done
done
################### SKIP ####################



################### IP-PACKING ####################
for (( index=0; index < numOfDataSets; index++ )); do
	for ip in $maxip
 	do
    	./bin/Release/MatrixGen IP-PACKING $k_list ${dataset_list[$index]} $ip $bits
	done
done
################### IP-PACKING ####################



