#include "input.h"
#include "bonded.h"
#include "ring_group.h"
#include "compute_gr.h"

void Bonded::Init(Atom* atom, int* total_ring_number, int* ringtot_number, int &frame, double &time)
{
	TITLE("Bonded","Init");

	// set the max number of rings
	int max_ring=RingGroup::nr_max;
	int ringsize = INPUT.ring_size_max;

	assert(ringsize>=3);

	RingGroup *RG1 = new RingGroup[ringsize+1]; //store the rings we want (Lixin's definition) 
	RingGroup *RG2 = new RingGroup[ringsize+1]; //store all rings, only used in Lingzhu's definition
	// starting from Ringsize 6, RG2 will be used to screening rings.
	for(int ig=3; ig<ringsize+1; ++ig)
	{
		RG1[ig].r = new Ring[max_ring];
		RG1[ig].na_of_ring = ig;
		RG2[ig].r = new Ring[max_ring];
		RG2[ig].na_of_ring = ig;
	}

	// set the max size of rings
	int na=INPUT.natom1;

	// set index of atoms in a ring
	short int* index = new short int[ringsize+1];

	// 3 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
			if( atom[k].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;
					if(INPUT.ring_definition==2) RG2[3].push(index,atom); // Ring Group 2 
				if( RG1[3].push(index,atom) ){
				
				if(INPUT.print_ring_min<=3 and 3<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 3 " << " " << i+1 << " " << j+1 << " " << k+1 << endl;	
	}}}}}}}}}

	if(ringsize==3) goto ring_analysis_done;

	// 4 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
			if( atom[l].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;
				if( RG1[4].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[4].push(index,atom); // Ring Group 2 
 				if( RG1[4].push(index,atom) ){
				if(INPUT.print_ring_min<=4 and 4<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 4 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << endl;	
	}}}}}}}}}}}}}

	if(ringsize==4) goto ring_analysis_done;

	// 5 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
			if( atom[m].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				if( RG1[5].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[5].push(index,atom); // Ring Group 2 
 				if( RG1[5].push(index,atom) ){
				if(INPUT.print_ring_min<=5 and 5<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 5 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << endl;	
	}}}}}}}}}}}}}}}}

	if(ringsize==5) goto ring_analysis_done;

	// 6 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
			if( atom[n].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;
				if( RG1[6].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[6].push(index,atom); // Ring Group 2 
				if( RG1[6].newbond(index, RG1, RG2, 6)){
 				if( RG1[6].push(index,atom) ){
				if(INPUT.print_ring_min<=6 and 6<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 6 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " "<< n+1 << endl;	
	}}}}}}}}}}}}}}}}}}}}

	if(ringsize==6) goto ring_analysis_done;

	// 7 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
			if( atom[o].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;
				if( RG1[7].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[7].push(index,atom); // Ring Group 2 
				if( RG1[7].newbond(index, RG1, RG2, 7)){
 				if( RG1[7].push(index,atom) ){
				if(INPUT.print_ring_min<=7 and 7<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 7 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 
					<< endl;
	}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==7) goto ring_analysis_done;

	// 8 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
	for(short int p=0; p<na; ++p){
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
			if( atom[o].bonded(p) == true ){
			if( atom[p].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;
				if( RG1[8].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[8].push(index,atom); //  Ring Group 2 
				if( RG1[8].newbond(index, RG1, RG2, 8)){
 				if( RG1[8].push(index,atom) ){
				if(INPUT.print_ring_min<=8 and 8<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 8 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 
					<< endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==8) goto ring_analysis_done;

	// 9 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
	for(short int p=0; p<na; ++p){
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
			if( atom[o].bonded(p) == true ){
	for(short int q=0; q<na; ++q){
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
			if( atom[p].bonded(q) == true ){
			if( atom[q].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;
				if( RG1[9].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[9].push(index,atom); // Ring Group 2 
				if( RG1[9].newbond(index, RG1, RG2, 9)){
 				if( RG1[9].push(index,atom) ){
				if(INPUT.print_ring_min<=9 and 9<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 9 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==9) goto ring_analysis_done;

	// 10 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
	for(short int p=0; p<na; ++p){
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
			if( atom[o].bonded(p) == true ){
	for(short int q=0; q<na; ++q){
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
			if( atom[p].bonded(q) == true ){
	for(short int r=0; r<na; ++r){
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
			if( atom[q].bonded(r) == true ){
			if( atom[r].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				if( RG1[10].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[10].push(index,atom); // Ring Group 2 
				if( RG1[10].newbond(index, RG1, RG2, 10)){
 				if( RG1[10].push(index,atom) ){
				if(INPUT.print_ring_min<=10 and 10<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 10 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1
					<< endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==10) goto ring_analysis_done;

	// 11 rings
	for(short int i=0; i<na; ++i){
	for(short int i2=0; i2<atom[i].nn; ++i2){short int j=atom[i].neighbour[i2]-1;
	for(short int i3=0; i3<atom[j].nn; ++i3){short int k=atom[j].neighbour[i3]-1;
		if(k!=i and k!=j){
	for(short int i4=0; i4<atom[k].nn; ++i4){short int l=atom[k].neighbour[i4]-1;
		if(l!=i and l!=j and l!=k){
	for(short int i5=0; i5<atom[l].nn; ++i5){short int m=atom[l].neighbour[i5]-1;
		if(m!=i and m!=j and m!=k and m!=l){
	for(short int i6=0; i6<atom[m].nn; ++i6){short int n=atom[m].neighbour[i6]-1;
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
	for(short int i7=0; i7<atom[n].nn; ++i7){short int o=atom[n].neighbour[i7]-1;
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
	for(short int i8=0; i8<atom[o].nn; ++i8){short int p=atom[o].neighbour[i8]-1;
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
	for(short int i9=0; i9<atom[p].nn; ++i9){short int q=atom[p].neighbour[i9]-1;
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
	for(short int i10=0; i10<atom[q].nn; ++i10){short int r=atom[q].neighbour[i10]-1;
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
	for(short int i11=0; i11<atom[r].nn; ++i11){short int s=atom[r].neighbour[i11]-1;
		if(s!=i and s!=j and s!=k and s!=l and s!=m and s!=n and s!=o and s!=p and s!=q and s!=r){
// new algorithm above
			if( atom[s].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				index[10]=s;
				if( RG1[11].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[11].push(index,atom); // Ring Group 2 
				if( RG1[11].newbond(index, RG1, RG2, 11)){
 				if( RG1[11].push(index,atom) ){
				if(INPUT.print_ring_min<=11 and 11<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 11 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1 << " " << s+1 << endl;
	}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==11) goto ring_analysis_done;


	// 12 rings
	for(short int i=0; i<na; ++i){
	for(short int i2=0; i2<atom[i].nn; ++i2){short int j=atom[i].neighbour[i2]-1;
	for(short int i3=0; i3<atom[j].nn; ++i3){short int k=atom[j].neighbour[i3]-1;
		if(k!=i and k!=j){
	for(short int i4=0; i4<atom[k].nn; ++i4){short int l=atom[k].neighbour[i4]-1;
		if(l!=i and l!=j and l!=k){
	for(short int i5=0; i5<atom[l].nn; ++i5){short int m=atom[l].neighbour[i5]-1;
		if(m!=i and m!=j and m!=k and m!=l){
	for(short int i6=0; i6<atom[m].nn; ++i6){short int n=atom[m].neighbour[i6]-1;
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
	for(short int i7=0; i7<atom[n].nn; ++i7){short int o=atom[n].neighbour[i7]-1;
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
	for(short int i8=0; i8<atom[o].nn; ++i8){short int p=atom[o].neighbour[i8]-1;
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
	for(short int i9=0; i9<atom[p].nn; ++i9){short int q=atom[p].neighbour[i9]-1;
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
	for(short int i10=0; i10<atom[q].nn; ++i10){short int r=atom[q].neighbour[i10]-1;
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
	for(short int i11=0; i11<atom[r].nn; ++i11){short int s=atom[r].neighbour[i11]-1;
		if(s!=i and s!=j and s!=k and s!=l and s!=m and s!=n and s!=o and s!=p and s!=q and s!=r){
	for(short int i12=0; i12<atom[s].nn; ++i12){short int t=atom[s].neighbour[i12]-1;
		if(t!=i and t!=j and t!=k and t!=l and t!=m and t!=n and t!=o and t!=p and t!=q and t!=r and t!=s){
			if( atom[t].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				index[10]=s;index[11]=t;
				if( RG1[12].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[12].push(index,atom); // Ring Group 2 
				if( RG1[12].newbond(index, RG1, RG2, 12)){
 				if( RG1[12].push(index,atom) ){
				if(INPUT.print_ring_min<=12 and 12<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 12 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1 << " " << s+1 << " " << t+1 << endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==12) goto ring_analysis_done;


	// 13 rings
	for(short int i=0; i<na; ++i){
	for(short int i2=0; i2<atom[i].nn; ++i2){short int j=atom[i].neighbour[i2]-1;
	for(short int i3=0; i3<atom[j].nn; ++i3){short int k=atom[j].neighbour[i3]-1;
		if(k!=i and k!=j){
	for(short int i4=0; i4<atom[k].nn; ++i4){short int l=atom[k].neighbour[i4]-1;
		if(l!=i and l!=j and l!=k){
	for(short int i5=0; i5<atom[l].nn; ++i5){short int m=atom[l].neighbour[i5]-1;
		if(m!=i and m!=j and m!=k and m!=l){
	for(short int i6=0; i6<atom[m].nn; ++i6){short int n=atom[m].neighbour[i6]-1;
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
	for(short int i7=0; i7<atom[n].nn; ++i7){short int o=atom[n].neighbour[i7]-1;
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
	for(short int i8=0; i8<atom[o].nn; ++i8){short int p=atom[o].neighbour[i8]-1;
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
	for(short int i9=0; i9<atom[p].nn; ++i9){short int q=atom[p].neighbour[i9]-1;
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
	for(short int i10=0; i10<atom[q].nn; ++i10){short int r=atom[q].neighbour[i10]-1;
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
	for(short int i11=0; i11<atom[r].nn; ++i11){short int s=atom[r].neighbour[i11]-1;
		if(s!=i and s!=j and s!=k and s!=l and s!=m and s!=n and s!=o and s!=p and s!=q and s!=r){
	for(short int i12=0; i12<atom[s].nn; ++i12){short int t=atom[s].neighbour[i12]-1;
		if(t!=i and t!=j and t!=k and t!=l and t!=m and t!=n and t!=o and t!=p and t!=q and t!=r and t!=s){
	for(short int i13=0; i13<atom[t].nn; ++i13){short int u=atom[t].neighbour[i13]-1;
		if(u!=i and u!=j and u!=k and u!=l and u!=m and u!=n and u!=o and u!=p and u!=q and u!=r and u!=s and u!=t){
			if( atom[u].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				index[10]=s;index[11]=t;index[12]=u;
				if( RG1[13].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[13].push(index,atom);
				if( RG1[13].newbond(index, RG1, RG2, 13)){
 				if( RG1[13].push(index,atom) ){
				if(INPUT.print_ring_min<=13 and 13<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 13 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1 << " " << s+1 << " " << t+1 << " " << u+1 << endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==13) goto ring_analysis_done;

	// 14 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
	for(short int p=0; p<na; ++p){
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
			if( atom[o].bonded(p) == true ){
	for(short int q=0; q<na; ++q){
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
			if( atom[p].bonded(q) == true ){
	for(short int r=0; r<na; ++r){
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
			if( atom[q].bonded(r) == true ){
	for(short int s=0; s<na; ++s){
		if(s!=i and s!=j and s!=k and s!=l and s!=m and s!=n and s!=o and s!=p and s!=q and s!=r){
			if( atom[r].bonded(s) == true ){
	for(short int t=0; t<na; ++t){
		if(t!=i and t!=j and t!=k and t!=l and t!=m and t!=n and t!=o and t!=p and t!=q and t!=r and t!=s){
			if( atom[s].bonded(t) == true ){
	for(short int u=0; u<na; ++u){
		if(u!=i and u!=j and u!=k and u!=l and u!=m and u!=n and u!=o and u!=p and u!=q and u!=r and u!=s and u!=t){
			if( atom[t].bonded(u) == true ){
	for(short int v=0; v<na; ++v){
		if(v!=i and v!=j and v!=k and v!=l and v!=m and v!=n and v!=o and v!=p and v!=q and v!=r and v!=s and v!=t and v!=u){
			if( atom[u].bonded(v) == true ){
			if( atom[v].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				index[10]=s;index[11]=t;index[12]=u;index[13]=v;
				if( RG1[14].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[14].push(index,atom);
				if( RG1[14].newbond(index, RG1, RG2, 14)){
 				if( RG1[14].push(index,atom) ){
				if(INPUT.print_ring_min<=14 and 14<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 14 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1 << " " << s+1 << " " << t+1 << " " << u+1 << " " << v+1 << endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
	
	if(ringsize==14) goto ring_analysis_done;

	// 15 rings
	for(short int i=0; i<na; ++i){
	for(short int j=i; j<na; ++j){
		if( j!=i ){
			if( atom[i].bonded(j) == true ){
	for(short int k=0; k<na; ++k){
		if(k!=i and k!=j){
			if( atom[j].bonded(k) == true ){
	for(short int l=0; l<na; ++l){
		if(l!=i and l!=j and l!=k){
			if( atom[k].bonded(l) == true ){
	for(short int m=0; m<na; ++m){
		if(m!=i and m!=j and m!=k and m!=l){
			if( atom[l].bonded(m) == true ){
	for(short int n=0; n<na; ++n){
		if(n!=i and n!=j and n!=k and n!=l and n!=m){
			if( atom[m].bonded(n) == true ){
	for(short int o=0; o<na; ++o){
		if(o!=i and o!=j and o!=k and o!=l and o!=m and o!=n){
			if( atom[n].bonded(o) == true ){
	for(short int p=0; p<na; ++p){
		if(p!=i and p!=j and p!=k and p!=l and p!=m and p!=n and p!=o){
			if( atom[o].bonded(p) == true ){
	for(short int q=0; q<na; ++q){
		if(q!=i and q!=j and q!=k and q!=l and q!=m and q!=n and q!=o and q!=p){
			if( atom[p].bonded(q) == true ){
	for(short int r=0; r<na; ++r){
		if(r!=i and r!=j and r!=k and r!=l and r!=m and r!=n and r!=o and r!=p and r!=q){
			if( atom[q].bonded(r) == true ){
	for(short int s=0; s<na; ++s){
		if(s!=i and s!=j and s!=k and s!=l and s!=m and s!=n and s!=o and s!=p and s!=q and s!=r){
			if( atom[r].bonded(s) == true ){
	for(short int t=0; t<na; ++t){
		if(t!=i and t!=j and t!=k and t!=l and t!=m and t!=n and t!=o and t!=p and t!=q and t!=r and t!=s){
			if( atom[s].bonded(t) == true ){
	for(short int u=0; u<na; ++u){
		if(u!=i and u!=j and u!=k and u!=l and u!=m and u!=n and u!=o and u!=p and u!=q and u!=r and u!=s and u!=t){
			if( atom[t].bonded(u) == true ){
	for(short int v=0; v<na; ++v){
		if(v!=i and v!=j and v!=k and v!=l and v!=m and v!=n and v!=o and v!=p and v!=q and v!=r and v!=s and v!=t and v!=u){
			if( atom[u].bonded(v) == true ){
	for(short int w=0; w<na; ++w){
		if(w!=i and w!=j and w!=k and w!=l and w!=m and w!=n and w!=o and w!=p and w!=q and w!=r and w!=s and w!=t and w!=u and w!=v){
			if( atom[v].bonded(w) == true ){
			if( atom[w].bonded(i) == true ){
				index[0]=i;index[1]=j;index[2]=k;index[3]=l;index[4]=m;
				index[5]=n;index[6]=o;index[7]=p;index[8]=q;index[9]=r;
				index[10]=s;index[11]=t;index[12]=u;index[13]=v;index[14]=w;
				if( RG1[15].sumUp(atom, index) ){
					if(INPUT.ring_definition==2) RG2[15].push(index,atom);
				if( RG1[15].newbond(index, RG1, RG2, 15)){
 				if( RG1[15].push(index,atom) ){
				if(INPUT.print_ring_min<=15 and 15<=INPUT.print_ring_max)
				ofs_running << " FRAME " << frame << " A_Ring 15 " << i+1 << " " << j+1 << " " << k+1 
					<< " " << l+1 << " " << m+1 << " " << n+1 
					<< " " << o+1 << " " << p+1 << " " << q+1 
					<< " " << r+1 << " " << s+1 << " " << t+1 << " " << u+1 << " " << v+1 << " " << w+1 << endl;
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

	if(ringsize==15) goto ring_analysis_done;

ring_analysis_done:

	ofs_running << setw(10) << "RingSize" << setw(15) << "ValidRingG1#" << setw(15) << "ValidRingG2#" << endl;
	for(int i=3; i<ringsize+1; ++i)
	{
		ofs_running << setw(10) << i << setw(15) << RG1[i].nr << setw(15) << RG2[i].nr << endl;
	}



	if(INPUT.ion_analysis==0)
	{
		// print out the final results
		ofs_running << setw(10) << "RingSize" << setw(15) << "Ring#" << setw(15) << "TotRing#" << endl;

		for(int i=3; i<ringsize+1; ++i)
		{
			total_ring_number[i]+=RG1[i].nr;
			ofs_running << setw(10) << i << setw(15) << RG1[i].nr << setw(15) << total_ring_number[i] << endl; 
		}
		
		// tot_rn: total_ring_number
		if(INPUT.tot_rn==1) 
		{
			int tmp=0;
			for(int i=3; i<ringsize+1; ++i)
			{
				tmp+=RG1[i].nr;
			}
			assert(tmp<1000);
			if(tmp>0) ringtot_number[tmp]++;
		}
		else if(INPUT.tot_rn==2) // second method, 2016-06-07
		{
			for(int ia=0; ia<INPUT.natom1; ++ia)
			{
				int tmp=0;
				for(int is=3; is<ringsize+1; ++is)
				{
					for(int ir=0; ir<RG1[is].nr; ++ir)
					{
						for(int ib=0; ib<RG1[is].r[ir].na_ring; ++ib)
						{
							const int index=RG1[is].r[ir].index_of_atoms[ib];
							if(index==ia)
							{
								++tmp;
							}
						}
					}
				}
				if(tmp>0) ringtot_number[tmp]++;
			}
		}
	}
	else //ion analysis or multiple jump
	{
		int* rr = new int[ringsize+1];
		int numallring=0;
		for(int is=3; is<ringsize+1; ++is)
		{
			rr[is]=0;
			for(int ir=0; ir<RG1[is].nr; ++ir)
			{
				for(int ia=0; ia<RG1[is].r[ir].na_ring; ++ia)
				{
					const int index=RG1[is].r[ir].index_of_atoms[ia];
					if(atom[index].IamIon)
					{
						++total_ring_number[is];
						++numallring;
						++rr[is];
						if(is>=INPUT.print_ionring_min and is<=INPUT.print_ionring_max) // mohan added 2016-08-11
						{
							ofs_running << " FRAME " << frame << " ION_RING " << is; 
							for(int ib=0; ib<RG1[is].r[ir].na_ring; ++ib)
							{
								ofs_running << " " << RG1[is].r[ir].index_of_atoms[ib]+1;
							}
							ofs_running << endl;
						}
						break;
					}
				}
			}
		}
		if(numallring>0) ++ringtot_number[numallring];
		ofs_running << " Frame " << frame << " time " << time << " number_of_all_rings " << numallring;
		// print again
		for(int is=3; is<ringsize+1; ++is)
		{
			ofs_running << " " << rr[is];
		}
		ofs_running << endl;
		delete[] rr;
	}

	
	delete[] index;

	if(INPUT.cal_gr==true)
	{
		//compute g(r) for each size of ring
		for(int ig=3; ig<ringsize+1; ++ig)
		{
			ofs_running << " Compute g(r) for ring size " << ig << endl;
			CGR.gr_of_ring(RG1[ig],atom);
		}
	}


	for(int ig=3; ig<ringsize+1; ++ig)
	{
		for(int ir=0; ir<RG1[ig].nr; ++ir)
			delete[] RG1[ig].r[ir].index_of_atoms;
		for(int ir=0; ir<RG2[ig].nr; ++ir)
			delete[] RG2[ig].r[ir].index_of_atoms;
		delete[] RG1[ig].r;
		delete[] RG2[ig].r;
	}
	delete[] RG1; 
	delete[] RG2; 

	return;
}
