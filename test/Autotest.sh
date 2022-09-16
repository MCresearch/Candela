#!/bin/bash
########################
MPICOMPILE=OFF
########################

narg=$#
if ((narg==0)); then
    mpi=$MPICOMPILE
else
    mpi=$1
fi
allcases=`ls | grep ^[0-9]`
echo -e "\033[33m---------------AUTOTEST---------------\033[0m"
for case in $allcases
do
    echo "  ==================================  "
    cd $case

    echo "  Running $case..."
    time_start=`date "+%H:%M:%S.%N"`
    python3 ../tool/check_one.py $mpi
    time_end=`date "+%H:%M:%S.%N"`
    start_seconds=`date --date="$time_start" +%s.%N`
    end_seconds=`date --date="$time_end" +%s.%N`
    during=`echo "$end_seconds-$start_seconds" |bc`
    echo "  Cost time: $during s"

    cd ../
    echo "  ==================================  "
    
done
echo -e "\033[33m--------------------------------------\033[0m"