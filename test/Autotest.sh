#!/bin/bash
allcases=`ls | grep ^[0-9]`
echo -e "\033[33m---------------AUTOTEST---------------\033[0m"
for case in $allcases
do
    echo "  ==================================  "
    cd $case

    echo "  Running $case..."
    time_start=`date "+%Y-%m-%d %H:%H:%S.%N"`
    python3 ../tool/check_one.py
    time_end=`date "+%Y-%m-%d %H:%H:%S.%N"`
    start_seconds=`date --date="$time_start" +%s.%N`
    end_seconds=`date --date="$time_end" +%s.%N`
    during=`echo "$end_seconds-$start_seconds" |bc`
    echo "  Cost time: $during s"

    cd ../
    echo "  ==================================  "
    
done
echo -e "\033[33m--------------------------------------\033[0m"