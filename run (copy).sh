#!/usr/bin/env bash



#rm -rf build
#mkdir build
#cd build
#cmake ../
#make

declare -a  dataset_list=("./data/MovieLens/" "./data/Netflix/")
#declare -a dataset_list=("./data/MovieLens/" "./data/Yelp/" "./data/Netflix/")
numOfDataSets=${#dataset_list[@]}

k_list="1,2,5,10,50"

maxip="100000"

bits="1024"

################### SKIP ####################
for (( index=0; index < numOfDataSets; index++ )); do
    ./bin/Release/MatrixGen SKIP $k_list ${dataset_list[$index]} $maxip $bits
done
################### SKIP ####################


################### AIPE-SS ####################
for (( index=0; index < numOfDataSets; index++ )); do
    ./bin/Release/MatrixGen AIPE-SS $k_list ${dataset_list[$index]} $maxip $bits
done
################### AIPE-SS ####################


################### IP-PACKING ####################
for (( index=0; index < numOfDataSets; index++ )); do
    ./bin/Release/MatrixGen IP-PACKING "10" ${dataset_list[$index]} $maxip $bits
done
################### IP-PACKING ####################


################### AIPE ####################
for (( index=0; index < numOfDataSets; index++ )); do
    ./bin/Release/MatrixGen AIPE "10" ${dataset_list[$index]} $maxip $bits
done
################### AIPE ####################
