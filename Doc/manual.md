# CANDELA User’s Manual 
This is the manual for CANDELA, which was developed in Mohan Chen’s research group in Peking University.

## Outline
In Section 1, we briefly introduce the functions of CANDELA. In Section 2, we introduce how to obtain and install CANDELA. In Section 3, a list of examples is provided. In Section 4, some key parameters are explained in detail. Finally, in section 5, the code structure is introduced.

## 1. Introduction
CANDELA derives from D310 package and D310 is the office name of Mohan Chen when he was a postdoc at Princeton University, where the code was created. The code was designed in a way to read in data from either Quantum Espresso, VASP, PROFESS, ABACUS, LAMMPS etc., most of the data are molecular dynamics trajectories. Next, CANDELA is mainly designed for post processing and provides results like radial distribution functions, bond angle distribution functions, diffusion coefficients, and many others. It is convenience for users who use multiple packages and would like to analyze the data in a same manner.

## 2. Download:
The code is stored in the github website and can be downloaded from the website with the following address:
https://github.com/MCresearch/Candela

## 3. Structure of the CANDELA Code

<font color="blue">Blue</font>: commonly used.

<font color="black">Black</font>: base modules.

<font color="red">Red</font>: analysis for liquid water system.

<font color="purple">Purple</font>: analysis for mechanical properties.

(Colors have not been completely painted yet.)
| files | introduction|
|:------|:------------|
|atoms.cpp, atoms.h  |	Define the class of atoms
|<font color="blue">bdf.cpp, bdf.h</font>	|Calculate the bond angle distribution functions
|<font color="blue">bdf_rcut.cpp, bdf_rcut.h</font>	|Calculate the bond angle distribution functions within a radius cutoff
|cell.cpp, cell.h	|Cell Information
|cellFile.cpp, cellFile.h|	Cell Files from various kinds of software
|cellABACUS.cpp|	Read in cell information from ABACUS
|cellLAMMPS.cpp	|Read in cell information from LAMMPS
|cellPROFESS.cpp	|Read in cell information from PROFESS
|cellPWmat.cpp	|Read in cell information from PWmat
|cellQE2.cpp	|Read in cell information from Quantum Espresso
|cellQE.cpp	|Read in cell information from Quantum Espresso
|cellRAW.cpp	|Read in cell information from RAW format
|cellVASP.cpp	|Read in cell information from VASP
|cellXYZ.cpp	|Read in cell information from XYZ format
|const.h	|Define constants
|<font color="red">density2D.cpp, density2D.h</font>	|For water, to plot two dimensional figures regarding the relation between density and length of covalent bonds.
|dielectric.cpp, dielectric.h	|Compute dielectric constants for liquid water.
|distri2D.cpp, 
|distri2D.h	|Compute angle-distance relation for liquid water (2D).
|distri3D_ions.cpp, 
|distri3D_ions.h	|Compute 3D distribution of atoms/MLWFs around ions(hydroxide,Cl), read by VESTA
|dsf.cpp, dsf.h	|Dynamic Structure Factors
|ele_conductivity.cpp,
|ele_conductivity.h	|Compute the electron conductivity.
|ext.cpp, ext.h	|Extend the cells.
|find_ion.cpp, find_ion.h	|Find the ion positions in the interface of water-oil system.
|gfun.cpp, gfun.h	|Define globally used functions.
|HBs.cpp, HBs.h	|Define hydrogen bonds in liquid water.
|<font color="purple">Honeycutt.cpp Honeycutt.h	</font>|Honeycutt analysis (not used yet).
|planarity.cpp, planarity.h	|Compute the hyper structures of OH- ion via the order parameter named planarity.
|ili_3D.cpp, ili_3D.h |	Generate 3D data of Instantaneous Liquid Interface (ILI)
|ili.cpp, ili.h	|Compute Instantaneous Liquid Interface (ILI)
|input.cpp, input.h	|Input data.
|insert.cpp, insert.h	|Randomly insert atoms into existing structures.
|iprof.cpp, iprof.h	|Ionic density profile.
|isf2.cpp, isf2.h	|Compute intermediate scattering function, new by Qianrui.
|isf.cpp, isf.h	|Compute intermediate scattering function, old.
|main.cpp	|Main function of CANDELA package.
|Makefile	|Makfile of CANDELA package.
|math.cpp, math.h	|Math functions.
|matrix3.cpp, matrix3.h	|3-dimensional matrix format.
|mdp2.cpp, mdp2.h	|Compute mean density profile based on mean liquid interfaces.
|mdp3.cpp, mdp3.h	|Analysis based on Mean Density Profile and Instantaneous Liquid Interface
|mdp.cpp, mdp.h	|Compute mean density profile based on instantaneous liquid interface (ILI).
|mj.cpp, mj.h	|Provide proton transfer data.
|movie_hexane.cpp, movie_hexane.h	|Print movie data for hexane project.
|msd.cpp, msd.h	|Compute the mean square displacements.
|msd_multiple.cpp, msd_multiple.h	|Compute the mean square displacements through multiple trajectories.
|pdf2d.cpp, pdf2d.h	|PDF for 2d materials.
|pdf5.cpp, pdf5.h	|Compute PDF for each shell of neighboring atoms.
|pdf_added.cpp, pdf_added.h	|The extra functions used in computing the PDF of systems.
|pdf.cpp, pdf.h	|Compute the radial distribution functions.
|powers.cpp, powers.h	|Compute the power spectra of system (FFT of VAF).
|presolvation.cpp, presolvation.h	|Compute the presolvation structure (yes/no).
|local_pseudopotential.cpp, local_pseudopotential.h	|Change the form of the readin pseudopotential.
|random.h	|Compute the random numbers.
|reorganize.cpp, reorganize.h	|Reorganize the input trajectories through nbin,
|binftream.cpp, binfstream.h	|Read in binary data.
|ssf.cpp, ssf.h	|Compute static structure factors
|ssf_selected.cpp, ssf_selected.h	|Compute static structure factors for selected q points.
|tetra_order.cpp, tetra_order.h	|Tetrahedral order parameter computed for each water molecule.
|vacuum.cpp, vacuum.h	|Add vacuum for a structure.
|vec3.h	Vector3 format.
|velcor.cpp, velcor.h	|Velocity autocorrelation functions.
|vel.cpp, vel.h	|Compute the distribution of velocities.
|void.cpp, void.h	|Create a void in a material.
|wannier1.cpp, wannier1.h	|Compute distributions of Wannier functions, including 1D, 2D, and 3D.
|wannier.cpp, wannier.h	|Compute Infrared Spectra from Wannier functions, plot distribution of dipoles.
|water.cpp, water.h	|Construct a water class for each read-in water.
|waterwire2.cpp, waterwire2.h	|Compute the free energy map of waterwire involving single and double proton transfer.
|waterwire.cpp, waterwire.h	|Compute the free energy map of waterwire involving single and double proton transfer.
|wavefunc.cpp, wavefunc.h	|Wave functions, added by Qianrui
|wfFile.cpp, wfFile.h	|Wave functions, added by Qianrui.
|wfPWmat.cpp, wfPWmat.h	|Wave functions, added by Qianrui.
|wfqe.cpp, wfqe.h	|Wave functions, added by Qianrui
|wfRead.cpp, wfRead.h	|Wave functions, added by Qianrui.
|write.cpp, write.h|	Wave functions, added by Qianrui.
|ww_compress.cpp, ww_compress.h	|Compute the compression of water wires.
|xsf.cpp, xsf.h	|Compute the properties from 3D data in XSF format.


## 4. List of Examples:
1).	[Radial distribution functions for O-O with varying cells](#rdf-O-O)

2).	[Bond angle distribution function for O triplets with varying cells](#adf-O)

3).	[Bond angle distribution function for O triplets with cutoff in varying cells](#bdf-vc)

4).	[Hydrogen bond analysis](#hba)

5).	[Mean square displacement](#msd)

6).	[Mean square displacement computed by multiple sections](#multi-msd)

7).	[Wannier centers and dipoles](#wannier-c-d)

8).	[Tetrahedrality](#tetrahedrality)

9).	[Reorganize](#reorganize)

10). [Radial distribution function of O*H in hydroxide solution](#rdf-OH)

11). [Hydrogen bonds for hydroxide](#hb-H)

12). [Multiple jumps](#multi-j)

13). [Planarity of hydroxide ion and its neighbors](#planarity-H)

14). [Movie](#movie)

15). [Distance 2D](#dis-2D)

16). [XSF 3D](#xsf-3D)

17). [Distance 2D Wannier](#dis-2D-wannier)

18). [PDF O(donate)-O(accept)](#pdf-d-a)

19). [Presolvation](#presolvation)

20). [Hydroxide Wannier](#wannier-HO)

21). [Waterwire2](#waterwire)

22). [Waterwire compress](#waterwire-compression)

23). [Static Structure Factor](#ssf)

24). [Instantaneous Liquid Interface](#ili)

25). [Infrared Spectra](#infs)

26). [Dynamic Structure Factor](#dsf)

27). [Electrical Conductivity](#ele-cond)


 
### 4.1).Radial distribution functions for O-O with varying cells<a id="rdf-O-O"></a>
```
  calculation  pdf # Pair Distribution Function.
  system water
  geo_in_type  QE
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836  # first 5 ps

  geo_out      pdf.txt # output pdf name.

  ntype        2        # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128
  dr           0.01     # delta r in real space
  rcut         6.22     # real space cutoff

  id1 O
  id2 H
  ele1 O
  ele2 O

  struf_dgx   0.05
  struf_ng    480
```

### 4.2）Bond Angle distribution function for O triplets with varying cells<a id="adf-O"></a>
 
 ```
calculation  bdf_rcut
  system water
  geo_in_type  QE  #PROFESS/VASP/QE/ABINIT/MESIA/XYZ/PIMD
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836 ! ignore first 5 ps

  ntype        2         # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128

  dr           0.01      # delta r in real space

  id1 O
  id2 H

  rcut1        4.0
  bdf_rcut     3.154
  bdf_dtheta   0.5     # d(theta) for degree between (0,180]

  ele1 O
  ele2 O

  func 1

  factor 0.9834679045
  x0 2.25
  y0 40
  nx 50
  ny 70
  dx 0.025
  dy 2.0
 ```
### 4.3）Bond Angle distribution function for O triplets with cutoff in varying cells<a id="bdf-vc"></a>
Notes: if func_b\==1, do nothing; func_b\==2, both bonded; func_b\==3, both are not boned; func_b\==4, one bonded, the other not bonded.
```
  calculation  bdf_rcut       #
  system water
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/MESIA/XYZ/PIMD
  geo_directory  ../../SCAN_volume/water.pos
  cell_file ../../SCAN_volume/water.cel
  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836 ! first 5 ps

  ntype        2         # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128

  dr           0.01      # delta r in real space

  id1 O
  id2 H

  rcut1        4.0
  bdf_rcut     3.154     # cutoff for O triplets
  bdf_out      3.154.dat  # output file name
  bdf_dtheta   0.5     # d(theta) for degree between (0,180]

  ele1 O
  ele2 O
  func_b  2

  factor 1.1542012927
  x0 2.25
  y0 40
  nx 50
  ny 70
  dx 0.02
  dy 2.0
```

 
### 4.4）Hydrogen bond analysis<a id="hba"></a>
```
  calculation  hbs  # hydrogen bond analysis
  system water
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/ABACUS
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_out      PBE0_OH_HBs.dat
  geo_1        1   
  geo_2        122517
  geo_interval 1   # geometry selected with this interval
  geo_ignore   20836

  ntype        2         # number of different types of atoms.
  natom        192       # total number of atoms.

  natom1 64 # O
  natom2 128 # H
  id1 O
  id2 H

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30
```
------------------------------------- FOR LAMMPS -----------------------------------------
```
  calculation  hbs  # hydrogen bond analysis
  system water
  geo_in_type  LAMMPS
  geo_directory  ../dump.lammpstrj
  geo_1        1
  geo_2        100 40001
  geo_interval 1
  geo_ignore   0 ! first 5 ps

  ntype        2         # number of different types of atoms.
  natom        96       # total number of atoms.

  natom1 32 # O
  natom2 64 # H
  id1 O
  id2 H

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  cartesian 1
  ```


### 4.5）Mean square displacement<a id="msd"></a>
```
  calculation  msd      # pair Distribution Function.
  system water
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_1        1
  geo_2       122517
  geo_interval 1
  geo_ignore   0

  geo_out      pdf.txt # output pdf name.

  ntype        2         # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128
  dr           0.01      # delta r in real space
  rcut         6.22     # real space cutoff

  id1 O
  id2 H
  ele1 O
  ele2 O

  struf_dgx   0.05
  struf_ng    480

  func 2

  cartesian 1
  msd_dt 0.005 
 ```
### 4.6）Mean square displacement computed by multiple sections<a id="multi-msd"></a>
Example 1 for QE and Example 2 for LAMMPS
```
  calculation  msd_multiple  # mean square displacements
  system hydroxide
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_1        1     
  geo_2        122517 
  geo_interval 1  

  ntype        2         # number of different types of atoms.
  natom        192       # total number of atoms.

  natom1 64 # O
  natom2 128 # H
  id1 O
  id2 H

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  msd_n  5  # number of msd needed
  msd_t0 5 # starting point of msd
  msd_t 12.20 # length of msd (in ps)
  msd_dt0 3 # difference between different different MSD
  msd_dt .00084661 # delta t between 2 snapshots
  msd_natom 192
  msd_stokes 0
  system water
```

```
calculation msd_multiple # mean square displacements
system water
geo_in_type LAMMPS 
geo_directory  ../lmp.nve.64/dump.lammpstrj
geo_1        1
geo_2        200001
geo_interval 1
geo_ignore   1000 ! first 5 ps

ntype 2 # number of different types of atoms.
natom 192 # total number of atoms.
natom1 64 # O
natom2 128 # H
id1 O
id2 H
rcut_oo 3.5
rcut_oh 1.24
acut_hoo 30
msd_n 10 # number of msd needed
msd_t0 0 # starting point of msd
msd_t 100 # length of msd (in ps)
msd_dt0 100 # difference between different different MSD
msd_dt 0.005 # delta t between 2 snapshots (ps)
msd_natom 192
msd_stokes 0


calculation  msd_multiple
geo_in_type PROFESS
geo_directory  ../out/
geo_1 4000
geo_2 400000
geo_interval 4
geo_ignore 0

ntype 2
natom 128
natom1 118
natom2 10
id1 L7
id2 L6
msd_single 1 #1：open msd for only one type 0: close(by default)
msd_type 0   # type id (starts from 0)
msd_natom 10

msd_n 10
msd_t0 1
msd_t 10
msd_dt0 10
msd_dt 0.01
msd_stokes 0
dt_snapshots 0.00025
```

### 4.7）Wannier centers and dipoles<a id="wannier-c-d"></a>

```
  calculation  wannier       # pair Distribution Function.
  system water
  geo_in_type  QE   # input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  wannier_file ../SCAN_volume/water.wfc
  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836

  ntype        2         # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128

  id1 O
  id2 H
  ele1 O

  nbands 256

  dr 0.001
  rcut 1.00

  dz 0.02  # for dipole moment
  rcut1 10.0 # for dipole moment

  nx 70
  ny 70
  x0 0.95
  y0 0.25
  dx 0.0025
  dy 0.005
```
 
### 4.8）Tetrahedrality <a id="tetrahedrality"></a>
```
  calculation  top  # hydrogen bond analysis
  system hydroxide
  geo_in_type  QE   # input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory  ../SCAN_volume/water.pos
  cell_file ../SCAN_volume/water.cel
  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836 ! first 5 ps

  ntype        2         # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128
  id1 O
  id2 H

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  bdf_rcut 3.15
```
 
### 4.9）Reorganize<a id="reorganize"></a>
```
  calculation  reorganize 
  system water          # system is water
  geo_in_type  QE       # PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory  ../SCAN_volume/water.pos
  wannier_file ../SCAN_volume/water.wfc
  cell_file ../SCAN_volume/water.cel

  geo_1        1
  geo_2        122517
  geo_interval 1
  geo_ignore   20836 ! first 5 ps


  ntype        2        # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128
  id1 O
  id2 H

  nbands 256 # number of bands for the system

  nbin 10
```
### 4.10）Radial distribution function of O*H in hydroxide solution
<a id="rdf-OH"></a>
```
  calculation  pdf       # pair Distribution Function.
  system hydroxide
  geo_in_type  QE   # PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory OH_PBE0_vdW.pos
  geo_1        1
  geo_2        95306
  geo_interval 1

  geo_out      pdf.txt # output pdf name.

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.
  natom1       64
  natom2       127
  dr           0.002      # delta r in real space
  rcut         6.2     # real space cutoff

  id1 O
  id2 H
  ele1 O
  ele2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509
```

### 4.11）Hydrogen bonds for hydroxide <a id="hb-H"></a>
```
  calculation  hbs  # hydrogen bond analysis
  system hydroxide
  geo_in_type  QE   # PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory OH_PBE0_vdW.pos
  geo_out      PBE0_OH_HBs.dat
  geo_1        1
  geo_2        95306
  geo_interval 1   # pick up geometry with this interval
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30
```

### 4.12）Multiple jumps <a id="multi-j"></a>
```
  calculation  mj  # multiple jump
  system hydroxide
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/ABACUS
  geo_2         87

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  func 2
```
 
### 4.13）Planarity of hydroxide ion and its neighbors<a id="planarity-H"></a>
Note: If delta is not set, then the criterion is not used.
Nacc set to 4 suggests that the hyper coordinated structures are selected,
Ndon set to 1 suggests that only those water molecules who donate 1 HB are selected.
```
  calculation  hyper  # instantaneous liquid interfaces
  system hydroxide
  geo_in_type  QE   # PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory OH-_PBE0_vdW.pos
  geo_out      hyper.dat
  geo_1        1     #
  geo_2        95306 #
  geo_interval 1   # pick up geometry with this interval
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut  6.0
  dr 0.1
  nacc 4
  delta -0.1
```
 
### 4. 14）Movie <a id="movie"></a>
```
  calculation  movie  # instantaneous liquid interfaces
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/ABACUS
  geo_directory OH-_PBE0_vdW.pos
  geo_out      movie.xyz
  geo_1        1   
  geo_2        62338 
  geo_interval 62338 
  ntype        2        # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  snatom 5 49 15 16 58 62

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509 # celldm3 in Angstrom

  func 2
```
 
### 4.15) Distance 2D <a id="dis-2D"></a>
```
  calculation  dist  
  system hydroxide
  geo_in_type  QE  
  geo_directory OH-_PBE0_vdW.pos
  geo_1        1       
geo_2         95306 
  geo_interval 1   

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut  12.0
  rcut1  4.0
  dr 0.1

  nx 61
  ny 21

  u1 41
  u2 41
  u3 61

  func 1
  ele1 O
```
### 4.16) XSF 3D<a id="xsf-3D"></a>
Nacc40-plot the *.xsf file (3D pot)
``` 
calculation  dist  
  system hydroxide
  geo_in_type  QE   
  geo_directory OH-_PBE0_vdW.pos
  geo_1        1     
  geo_2        95306 
  geo_interval 1   
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut  12.0
  rcut1  4.0
  dr 0.1

  nx 61
  ny 41

  nacc 40
  func 2

  u1 41
  u2 41
  u3 61

  ele1 O
```
Add the following header in the front part of dist3D.dat file, change the name to 3d.xsf

```
CRYSTAL
PRIMVEC
    8.0 0 0
    0 8.0 0
    0  0 12.0
PRIMCOORD
         2 1
O    4.0 4.0 6.0
H    4.0 4.0 7.0
 BEGIN_BLOCK_DATAGRID_3D
 3D_PWSCF
 DATAGRID_3D_UNKNOWN
       41 41 61
 0.000  0.000 0.000
8 0 0
0 8 0
0 0 12
```

### 4.17) Distance 2D wannier<a id="dis-2D-wannier"></a>
```
  calculation  dist  # instantaneous liquid interfaces
  system hydroxide
  geo_in_type  QE   #PROFESS/VASP/QE/ABINIT/ABACUS
  geo_directory OH_PBE0_vdW.pos
  wannier_file  OH_PBE0_vdW.wfc
  nbands 256
  geo_1        1     
  geo_2        95306 
  geo_interval 1   
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut  3.0
  rcut1  1.0
  dr 0.1

  nx 60
  ny 20

  #nacc 30
  func 4
```

### 4.18) PDF O(donate)-O(accept) <a id="pdf-d-a"></a>
```
  calculation  pdf
  system hydroxide
  geo_in_type  QE   
  geo_directory OH-_PBE0_vdW.pos
  geo_1        1
  geo_2       10000 95306

  geo_interval 1

  geo_out      pdf.txt # output pdf name.

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.
  natom1       64
  natom2       127
  dr           0.02      # delta r in real space
  rcut         6.22     

  id1 O
  id2 H
  ele1 O
  ele2 O

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  nacc 4
  func 2 # g(O_donate, O_accept)
```

 
### 4.19) Presolvation<a id="presolvation"></a>
```
  calculation  pre  # hydrogen bond analysis
  system hydroxide
  geo_in_type  QE   
  geo_directory OH_PBE0_vdW.pos
  geo_1        1     
  geo_2       95306 
  geo_interval 1   

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut 5
  dr 0.01
``` 
### 4.20) Hydroxide Wannier<a id="wannier-HO"></a>
Note: add “nacc 5” can select those hydroxide ions with 5 accepted HBs
```
  calculation  wannier  
  system hydroxide
  geo_in_type  QE   
  geo_directory  OH-_PBE0_vdW.pos
  wannier_file  OH_PBE0_vdW.wfc
  geo_1        1     
  geo_2        95306   
  geo_interval 1   
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.
  nbands       256

  natom1 64 # O
  natom2 127 # H
  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  dr 0.005
  rcut 1.00
  ele1 O
```
 
### 4.21) Waterwire 2<a id="waterwire"></a>
```
  calculation  waterwire2  # hydrogen bond analysis
  system hydroxide
  geo_in_type  QE   
  geo_directory OH_PBE0_vdW.pos
  geo_1        1     
  geo_2        95306 
  geo_interval 1   
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H

  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut 6.5
  dr 0.01

  nx 60
  ny 40
  x0 -1.5
  y0 4.5
  dx 0.05
  dy 0.05

  func 2
#  factor 2.4783510188 #100ps
#  factor 1.2391755094 #50ps only if output=10
#  factor .6195877547 #50ps only if output=5
#  factor .3717526527 #30ps
  factor .4956702036 #40ps
#  func_b 31
``` 
Waterwire 2: 41
```
  calculation  waterwire2  
  system hydroxide
  geo_in_type  QE   
  geo_directory OH_PBE0_vdW.pos
  geo_1        1     
  geo_2        95306 
  geo_interval 1   
  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H

  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  rcut 6.5
  dr 0.01

  nx 60
  ny 40
  x0 -1.5
  y0 4.5
  dx 0.05
  dy 0.05
  func 2
  factor .4956702036 #40ps
  func_b 41
```

### 4.22) Waterwire compression<a id="waterwire-compression"></a>
```
  calculation  ww_compress  # hydrogen bond analysis
  system hydroxide
  geo_in_type  QE   
  geo_directory  OH_PBE0_vdW.pos
  geo_1        1     
  geo_2        95306 
  geo_interval 1   

  ntype        2         # number of different types of atoms.
  natom        191       # total number of atoms.

  natom1 64 # O
  natom2 127 # H

  id1 O
  id2 H

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 12.444655509

  rcut_oo 3.5
  rcut_oh 1.24
  acut_hoo 30

  dr 0.05
  rcut 8.0
```

### 4.23) Static structure factor<a id="ssf"></a>
A. ssf_selected
The name of the code for calculating static structure factors is “‘ssf_selected.cpp” in the source directory. The code can be parallerized.

1. Prepare the input file
```
calculation   ssf_selected  # command of calculating static structure factor with selected k
geo_in_type  PROFESS  #input type of geometry file: PROFESS/VASP/QE/ABINIT/ABACUS
geo_directory ../../md_files    
geo_1         0            
geo_2         999          
geo_interval  2             
ssf_out       Li_ssf.txt      # output static structure factor name.
ntype         1            # number of different types of atoms.
natom         6750        # total number of atoms.
struf_dgx     .1080405600   # delta G in G space, 2pi/a, a in the unit of Angstrom
struf_dgy     .1080405600   # delta G in G space
struf_dgz     .1080405600   # delta G in G space
struf_ng      60            # number of G points with delta G described above
```

2. Prepare the “INPUT” file as described above and run CANDELA.exe. The code will generate an output file named “SSF.input0”. There number of lines in the file “SSF.input0” depends on the largest G vector, each line represents one particular G vector with four numbers, namely, Gx, Gy, Gz, and |G|.
3. Use the following command to sort all G vectors in an ascending order:
sort -n -k4 SSF.input0 > SSF.input1
Next, delete some of the G vectors in this SSF.input1 file and only reserve those G vectors that you want to calculate. The principle here is that you can delete G vectors that aretoo close in values, thus the reduced number of G points will save the computational time for calculating the static structure factor.
4. Add the number of lines (the number of G vectors you want to compute) as the first line in SSF.input1 and run CANDELA again. The code will generate a new file named “SSF.input”. Next, delete both SSF.input0 and SSF.input1 and only keep the “SSF.input” file.
5. Add the number of k points in the SSF.input file as the first line and now you can submit the final job to calculate the static structure factor.

 
B. ssf 
The code for calculating the static structure factors (SSF) is ssf.cpp in the source directory.SSF is calculated partly with parallel_script.sh and averaged through collect_ssf.cpp.

1.Modify the script slurm.sh or files with other suffix names (like .bsub) to submit jobs.
2.Modify the parallel_script.sh
a. modify some parameters
target_dir: directory of geometry
split_to_nfile: how many parts I need
interval: geometry interval
max,min: beginning and ending geometry
addcount: restart from which part

b. modify INPUT in parallel_script.sh
```
calculation   ssf           # pair Distribution Function.
geo_in_type   PROFESS       # input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA
geo_directory $target_dir 
geo_1         $start_geo        # IN PROFESS, start from file: ion.*.dat, where * is 'geo_1'
geo_2         $stop_geo         # IN PROFESS, end at file : ion.*.dat, where * is 'geo_2'
geo_interval  $interval    # geometry interval between geo_1 and geo_2

ssf_out       Al_ssf.txt    # output static structure factor name.
ntype         1             # number of different types of atoms.
natom         864           # total number of atoms.
struf_dgx     .2586375438   # delta G in G space, 2pi/a
struf_dgy     .2586375438   # delta G in G space
struf_dgz     .2586375438   # delta G in G space
struf_ng      18            # number of G points
```
Bold parameters should be modified.
c. Modify slurm.sh in parallel_script.sh to the right name.
3.Type: sh parallel_script  and wait till it’s over.
4. Modify collect_ssf.cpp
split_to_nfile: how many parts in total
line:　how many lines of each ssf_out file
5. Compile collect_ssf.cpp and use it to get Final_ssf.txt
 
### 4.24) Instantaneous Liquid Interface<a id="ili"></a>
```
  calculation  ili  # instantaneous liquid interfaces
  geo_in_type  QE   
  geo_directory ../60.pos
  geo_out      ili.dat
  geo_1        1     
  geo_2        7441   
  geo_interval 1   
  ntype        3         # number of different types of atoms.
  natom        407       # total number of atoms.

  natom1 96 # O
  natom2 275 # H
  natom3 36 # C
  id1 O
  id2 H
  id3 C

  celldm1 12.444655509
  celldm2 12.444655509
  celldm3 26.9617536455 # celldm3 in Angstrom

  zeta 2.4
  d 3
  nx 30
  ny 30

  ntry 20
  z0 24
  dz -0.5

  ele1 O
  ref_den 0.016
```

### 4.25) Infrared Spectra<a id="infs"></a>

Step1: use ‘calculation wannier’ to generate all_dipole.dat and all_vdipole.dat files

Step2: obtain the above two files and change calculation to `infrared`

```
  calculation  infrared wannier  # step 1: wannier step2: infrared
  system water          # system is water
  geo_in_type  QE       
  geo_directory  ./reorganize10/new.pos
  wannier_file ./reorganize10/new.wfc
  cell_file ./reorganize10/new.cel
  geo_1        1       # starting geometry index
  geo_2        39411   # ending geometry index
  geo_interval 1       # interval between starting and ending geometries
  geo_ignore   0       # ignore the first xxx geometries

  ntype        2        # number of different types of atoms.
  natom        192      # total number of atoms.
  natom1       64
  natom2       128
  id1 O
  id2 H

  ele1 O  # please do not change this value

  nbands 256 # number of bands for the system

  factor 0.05 // only for func_b=2 in infrared mode
  bdf_rcut 0.1 // only for func_b=2 in infrared mode

  dr 0.001 # for MLWF distribution
  rcut 1.00 # for MLWF distribution

  dz 0.02  # for dipole moment distribution (in Debye)
  rcut1 10.0 # for dipole moment range (in Deybe)

  tcor 600 # correlation with 'tcor' snapshots for  correlation function
 dt_snapshots # the delta_t of two snapshots, unit is ps
  func_b  2
```

### 4.26) Dynamic Structure Factor<a id="dsf"></a>

Step 1 ：Calculate intermediate scattering function.
```
calculation     isf2 # Pair Distribution Function.
geo_in_type     PROFESS
geo_directory   ../out/
isf_outfile 2.00864isf.txt #outfile for isf

geo_1           10000  #(make sure: interval*(nt+nT)<=geo_2-geo_1)
geo_2           100000
geo_interval    1

natom           864      # total number of atoms.
natom1      864
ntype       1
isf_nt1     10000  #The range of time for ISF.
isf_nt2     80000  #The number of time for averaging ISF.
dt_snapshots    0.00025 #The time between two continuous snapshots.

isf_target_q    2.00 # The target q of ISF. (in A^-1)
isf_ngx     20      # The number of dg in x-direction for searching target_q
isf_ngy     20
isf_ngz     20
isf_dgx     0.25863754383  #(in A^-1) The minimal step of q.
isf_dgy     0.25863754383
isf_dgz     0.25863754383
```
Step 2 : Calculate  DSF and get maximam.
	   Modify tcut,dw,nwt,interval in onedsf.sh
	   tcut : The cut time for ISF. After this cut time, the ISF will decay exponentially. 
	   dw : the minimal energy step for DSF (in eV)
	   nwt : the maximal number of dw for DSF
	   interval : dt_snapshots in INPUT

`sh onedsf.sh` or `python onedsf.py`(recommended)
 
### 3.27) Electric Conductivity<a id="ele-cond"></a>

Step 1: Create directories containing different snapshots. (The order of folders should be 1,2,3…)
```
calculation write
write_cartesian 0
geo_in_type PWmat
geo_directory MOVEMENT
geo_1 1
geo_2 9000
geo_ignore 7100
geo_interval 100
headfile    headqe.txt
tailfile    tailqe.txt
geo_out     al.md.in
ntype   1
id1 Al
natom   64
natom1  64
```
Step 2: Calculate electric conductivity.
```
calculation      ele_conductivity
wf_in_type       PWmat/QE1/QE2/ABACUS
multi_directory path_to_multi_files #the directory of multiple snapshots files. (only valid when nscf > 0)
wfdirectory      path_to_wfs #(if nscf==0, then we read wf in path_to_wfs; else if nscf=Nfolders(Nfolders>0), then we read wf in path_to_multi_files/[1,2,…,Nfolders]/path_to_wfs)

wcut    12    #maximal energy (in eV) for conductivity
dw      0.9   #minimla energy step (in eV) for conductivity
fwhm    0.1   #Full Width at Half Maximam of Gauss or Lorentz function. Only useful when smear is 1 or 2
（or   n_fwhm  2    
       fwhm    0.1 0.2 #use it if we want to print results with different fwhm at the same time.)
smear   0     #0(default) for linear approx.; 1 for Gauss approx.; 2 for Lorentz approx.
nscf			  #0(default, only calculate one snapshot) or the number of snapshots used
temperature 10000 # in K
error_con    0/1 #if print the std of multi snapshots (only valid when nscf > 0)
tpk          1  #a factor convert student distribution to normal distribution
 ```


## Key words
The units of input atomic coordinates depend on the package you used for molecular dynamics simulations. For example, the output files of QE are all in Bohr, and the output files of LAMMPS are typically in Angstroms. We design specific interfaces for different simulation packages, and all data once read into CANDELA, are transformed into data with Angstroms.


