#!/bin/bash
allcases=`ls | grep ^[0-9]`
echo -e "\033[33m===============AUTOTEST===============\033[0m"
for case in $allcases
do
    echo "========================"
    echo "Running $case..."
    cd $case
    python3 ../tool/check_one.py
    cd ../
    echo "========================"
done