#!/bin/bash

# directory of trajectories
target_dir=../../md_files

# how many processors I need
split_to_nfile=10

# geometry interval 
# 0.02 fs * 20 = 0.4 fs
interval=20

# how many geometry files I need
max=50000
size_each=`echo "$max/$split_to_nfile" | bc`

echo "size_each = $size_each"

count=0
while [ $count -lt $split_to_nfile ]
do

start_geo=`echo "$count*$size_each" | bc` 
stop_geo=`echo "$start_geo+$size_each-1" | bc`
echo "$start_geo $stop_geo"


dir="cal_ssf_$count"
test -d $dir || mkdir $dir

cat > $dir/INPUT << EOF
calculation   ssf           # pair Distribution Function.
geo_in_type   PROFESS       # input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA
geo_directory $target_dir 
geo_1         $start_geo        # IN PROFESS, start from file: ion.*.dat, where * is 'geo_1'
geo_2         $stop_geo         # IN PROFESS, end at file : ion.*.dat, where * is 'geo_2'
geo_interval  $interval         # geometry interval between geo_1 and geo_2

ssf_out       Li_ssf.txt    # output static structure factor name.
ntype         1             # number of different types of atoms.
natom         128           # total number of atoms.
struf_dgx     .4384633635   # delta G in G space, 2pi/a
struf_dgy     .4384633635   # delta G in G space
struf_dgz     .4384633635   # delta G in G space
struf_ng      8            # number of G points
EOF

cp slurm.sh ./$dir/

cd $dir
sbatch slurm.sh
cd ..

let count++
done
