# How to develop CANDELA?
## 1. Autotest
After we add some new functions into CANDELA, we must add corresponding tests.
Directory **test** contains:

**001_PDF**,...: test cases

**Autotest.sh**: shell script to run auto tests.

**tool**: some tool functions

**geo**: contains geometry files

### a. Preparation for one case:

**compare**:
```
threshold 1e-4
enable_mpi ON
pdf.ref  pdf.txt
other.ref other.result
...
```
threshold: When the relative errors between reference results and calculating results are larger than threshold, the tests will fail.

enable_mpi: Whether this case supports MPI or not.

The remaining lines are reference files and result files. Each pair is in one line.

**INPUT**:
```
calculation  pdf # Pair Distribution Function.
system Al
geo_in_type  LAMMPS
geo_directory ../geo/Al64.dump
geo_1        0
geo_2        20
geo_interval 1
geo_ignore   4

geo_out      pdf.txt # output pdf name.

ntype        1        # number of different types of atoms.
natom        64   # total number of atoms.
natom1       64
rcut         6
dr           0.01     # delta r in real space
```
**pdf.ref / other.ref / ...**:
The reference results are in them.

### b. Run auto test
method 1:
`make CC=g++ TEST=ON`
Autotest.sh will be executed after serial-version **candela** is compiled.
```
......
g++ -O3  -fsanitize=address -fno-omit-frame-pointer -D__DEBUG -c src/DoubleDonor.cpp -o obj/DoubleDonor.o
g++ -O3  -fsanitize=address -fno-omit-frame-pointer -D__DEBUG -c src/qe_input.cpp -o obj/qe_input.o
g++ -O3  -fsanitize=address -fno-omit-frame-pointer -o bin/candela obj/velcor.o ... *.o
---------------AUTOTEST---------------
  ==================================  
  Running 001_PDF...
  One Processor      PASS
  Cost time: .174286415 s
  ==================================
  ......
--------------------------------------
```
`make CC=mpicxx TEST=ON`
Autotest.sh will be executed after parallel-version **candela** is compiled. Multi-processor cases will be tested at the same time.
```
......
mpicxx -O3  -fsanitize=address -fno-omit-frame-pointer -D__DEBUG -c src/DoubleDonor.cpp -o obj/DoubleDonor.o
mpicxx -O3  -fsanitize=address -fno-omit-frame-pointer -D__DEBUG -c src/qe_input.cpp -o obj/qe_input.o
mpicxx -O3  -fsanitize=address -fno-omit-frame-pointer -o bin/candela obj/velcor.o ... *.o
---------------AUTOTEST---------------
  ==================================  
  Running 001_PDF...
  One Processor      PASS
  Multi Processor    PASS
  Cost time: .470076373 s
  ==================================
  ......
--------------------------------------
```
***All developers should run 
`make CC=g++ TEST=ON`
and `make CC=mpicxx TEST=ON`
to make sure some functions are not be destroyed.***

method 2:
After **candela** is compiled.

`make test` or `cd test;sh Autotest.sh`
It will run one-processor tests.

`make test CC=mpicxx` or `cd test;sh Autotest.sh ON`
It will run one-processor and multi-processor tests.


