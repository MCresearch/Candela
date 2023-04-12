#!/bin/bash

allcases=`ls | grep ^[0-9]`
for case in $allcases
do
    cd $case
    python3 ../tool/reset.py
    cd ../
    
done