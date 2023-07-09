#!/bin/bash
if (($# == 0)); then
    allcases=`ls | grep ^[0-9]`
elif (($# > 0)); then
    allcases=$@
fi
for case in $allcases
do
    cd $case
    python3 ../tool/reset.py
    cd ../
done