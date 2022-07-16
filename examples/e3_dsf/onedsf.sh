#!/bin/bash
#A script to calculate DSF and get maximam
#w: meV  t:ps
wcut=10 # make sure the maximam position is larger than wcut
a[0]=0.12
a[1]=0.45
a[2]=0.29
a[3]=0.22
a[4]=0.179
a[5]=0.162
a[6]=0.145
i=0;
echo "w max"
for kk in "2.00" #"0.24" "0.31" "0.41" "0.51" "0.60" "0.70"
do  
file="${kk}864isf.txt";
outfile="${kk}864dsf.txt"
nwt=2000;
dw=0.2
tcut=${a[i]};
((i++))
interval=0.00025;
sigma=$tcut;
awk -v nwt=$nwt -v sigma=$sigma -v tcut=$tcut -v dw=$dw -v inte=$interval 'BEGIN{
	pi=3.1415926;
	pre=0.33;
	sigma=pre*sigma;
	for (nw=0;nw<nwt;nw++)
	{
		re[nw]=0;
	}
}{
	pp=1;
	if($1==0) pp=0.5;
	for (nw=0;nw<nwt;nw++)
	{
		if ($1<=tcut)
		{
			gaus=1;
		} 
		else 
		{
			phi=-1*($1-tcut)/sigma;
			if(phi<-100) {gaus=0}
			else {gaus=exp(phi)}
		}
		re[nw]+=pp*$2*cos(nw*dw*$1)*inte*gaus;
	}
}END{	
	for (nw=0;nw<nwt;nw++)
	{
		printf "%g %g\n",nw*dw/1.5192665814,re[nw]/2/pi*100*1.519*2*1.2286;
	}


}' $file>$outfile

cat $outfile|awk -v wcut=$wcut 'BEGIN{
max=0;w=0
}
{
	if($1>wcut)
	{
		if ($2>max)
		{
		max=$2;
		w=$1;
		}
	}	
}
END{
printf "%g %g\n",w,max}'
done
exit 0
