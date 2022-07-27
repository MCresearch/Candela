# An introduction to Candela
- [Overview](#Overview)
- [Basic usage and variables](#Basic-usage-and-variables)
    - [Compilation](#Compilation)
    - [Using Candela](#Using-Candela)
- [List of supported analyses](#List-of-supported-analyses)
    - [Water-related analyses](#Water-related-analyses)
        - [pdf](#pdf)
        - [hbs](#hbs)
        - [bdf_rcut](#bdf_rcut)
        - [dist](#dist)
        - [dist2](#dist2)
        - [hb_stat](#hb_stat)
        - [hb_stat2](#hb_stat2)
        - [hyper](#hyper)
        - [incremental_pdf](#incremental_pdf)
        - [mj](#mj)
        - [movie/movie2/OH_movie](#movie/movie2/OH_movie)
        - [reorganize](#reorganize)
        - [ssf](#ssf)
        - [tetra_order](#tetra_order)
        - [wannier](#wannier)
        - [orientation_tcf](#orientation_tcf)
        - [msd](#msd)
        - [msd_multiple](#msd_multiple)
        - [special_msd](#special_msd)
        - [angular_jump](#angular_jump)
        - [hb_correlation/hb_correlation2](#hb_correlation/hb_correlation2)
        - [nonhb_correlation/nonhb_correlation2/nonhb_correlation3](#nonhb_correlation/nonhb_correlation2/nonhb_correlation3)
        
- [Analysis program template](#Analysis-program-template)
    - [Basic classes in Candela](#Basic-classes-in-Candela)
    - [The structure of a typical analysis C++ file](#The-structure-of-a-typical-analysis-C++-file)
## Overview
Candela is short for Collection of ANalysis DEsigned for Large-scale Atomic simulations. It is developed by [MCresearch](https://github.com/MCresearch) to conduct analyses on MD trajectory in different formats. The analyses herein are mostly targeted at liquid water and warm condensed matter system. 

Section [Basic usage and variables](#Basic-usage-and-variables) introduces the usage of Candela;

Section [List of supported analyses](#List-of-supported-analyses) introduces the implemented analyses in Candela for now;

Section [Analysis program template](#Analysis-program-template) introduces the structure of a typical analysis subprogram, in case you may need to add new analysis program to Candela.

## Basic usage and variables 
### Compilation
Type `make` in the `Candela` directory to compile the executable `Candela.exe` in `\bin` directory. The default compiler is g++, which could be switched to other compiliers. 

### Using Candela
First, you have to prepare an `INPUT` file containing the input parameters of Candela. The following shows an example of the `INPUT` file:

```bash
  calculation  pdf       # pair Distribution Function (radial distribution function).
  system hydronium
  geo_in_type  QE   # input type of geometry file: PROFESS/VASP/QE/ABINIT/MESIA
  geo_directory ../H_PBE.pos	
  geo_1        1
  geo_2        10 
  geo_interval 1 

  geo_out      pdf.dat # output pdf name.

  ntype        2         # number of different types of atoms.
  natom        190       # total number of atoms.
  natom1       63
  natom2       127
  dr           0.01      # delta r in real space 
  rcut         6.22     # real space distance to evaluate the structure factor

  id1 O
  id2 H
  ele1 O
  ele2 O

  celldm1 12.444655509 
  celldm2 12.444655509 
  celldm3 12.444655509 
```

In each line the first item specifies the name of the parameter, and you may fill in the value of the parameter after it. The following chart displyas the general parameters and their purposes.

 Name  | Type  | Example                                                      | Discription                                                      |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| **calculation** | String | pdf | Name of the analysis
| **system** | String | hydronium | Name of the system
|  **geo_in_type** | String  | QE | Format of the MD trajectory file
| **geo_directory** | String | ../H_PBE.pos | directory of the MD trajectory file
| **geo_1** | Int | 1 | Starting index of snapshot in the MD trajectory file. In most cases it has to be 1.
| **geo_2** | Int | 10 | Final index of snapshot in the MD trajectory file
| **geo_interval** | Int | 1 | The program analyse 1 snapshot for every `geo_interval` snapshots
| **ntype** | Int | 2 | Number of atom species. The maximum number of spicies is mostly 3 or 4.
| **natom** | Int | 190 | Number of atoms
| **natom1** | Int | 63 | Number of atoms of the 1st species
| **natom2** | Int | 127 | Number of atoms of the 2nd species
| **id1** | String | O | Element name of the 1st species of atoms
| **id2** | String | H | Element name of the 2nd species of atoms
| **celldm1** | float | 12.444655509 | Length of the cell in the 1st dimension
| **celldm2** | float | 12.444655509 | Length of the cell in the 2nd dimension
| **celldm3** | float | 12.444655509 | Length of the cell in the 3rd dimension
| **dt_snapshots** | float | 4.83776865e-5 | A single MD timestep of `cp.x` in ps. The value is needed for trajectory over 100 ps.

For `geo_in_type`, currently supported formats include pw.x (`QE2`), cp.x (`QE`), VASP (`VASP`), ABACUS (`ABACUS`), PWmat (`PWMAT`) and so on.

The physical quantities used in the program are all in Angstrom/ps/eV.

Other parameters are needed for specific `calculation`, which you may refer to [List of supported analyses](#List-of-supported-analyses) where a detailed summary of suppoted `calculation` of analyses is given. You may refer to the example `INPUT` files in `/examples` as well for different kinds of analyses.

## List of supported analyses
### Water-related analyses
- **pdf** <a id="pdf"></a>
Pair distribution function, also radial distribution function. The following parameters are needed for `pdf` calculation:

    - `rcut` Radius cutoff of pair distribution function calculation.

    - `dr` Grid length of radius.

    - `id1` Name of the first type of atom of the pdf. For example, for $g_{OH}(r)$, `id1` is O.

    - `id2` Name of the second type of atom of the pdf. For example, for $g_{OH}(r)$, `id2` is H.
- **hbs** <a id="hbs"></a> Hydrogen bond (HB) analysis. Multiple properties are calculated, among which the most used ones are HB number statistics (at the end of `running0.log`) and proton transfer (PT) time record (`trans.dat`, only for `system` is `hydroxide` or `hydronium`). Another important function in this module is 

    ```cpp
    void HBs::setup_water(const Cell &cel, Water *water)
    ```

    The function is called in EVERY analysis involving water to set up the `water` instance, recording multiple information of each water molecule, including O, H indexes, O-H distance, accepted/donated water index, accepted/donated angle and so on. The following parameters are additionlly needed for setting up `water` instance.

    - `rcut_oo` The radial cutoff for hydrogen-bonded water molecules. The default value is 3.5 (Angstroms), following Luzar et al., 379, 55, Nature (1996).

    - `acut_hoo` The H-O-O angular cutoff for hydrogen-bonded water molecules. The default value is 30 (degrees), following Luzar et al., 379, 55, Nature (1996).

    - `rcut_oh` The radial cutoff to calculate neighboring hydrogen atoms of oxygen atom. Default value is 1.24 (Angstroms).

- **bdf_rcut**<a id=bdf_rcut></a>
The analysis calculates the OOO angle distribution of the central water and its two neighboring water molecules. (`bdf_rcut1` has the same function.) It also calculates the free energy map of OOO angle and O-O distance. (Figure. 3D-3E in Chen, et al., 114, 10846, Proc. Nat. Acad. Sci. (2017))

- **dist**<a id=dist></a>
The analysis calculates the 2D/3D spatial distribution of water molecules around hydroxide ion. The output file could be put into Vesta to generate the 3D distribution plot. (Figure. 4a-c of Chen et al. 10, 413, Nat. Chem. (2018)) The following parameters are also needed for the analysis:
    - `nx`/`ny` number of grids on the x/y direction of the 2D distribution.
    - `u1`/`u2`/`u3` number of grids on the x/y/z direction of the 3D distribution.
    - `rcut` radius cutoff of the distribution. See the code of `src/dist.cpp` for a more detailed understanding.
    - `rcut1` radius cutoff of the distribution. See the code of `src/dist.cpp` for a more detailed understanding.


- **dist2**<a id=dist2></a>
The analysis calculates the 3D spatial distribution of first-shell (`func=2`) and second-shell (`func=3`) water molecules or the spatial distribution of Wannier function centers (`func=1`). The Wannier function centers distribution only support `qe` format (cp.x) for now. The `wannier_file` and `nbands` (number of Wannier centers) have to be specified when involving calculation of Wannier centers. The output file could be put into Vesta to generate the 3D distribution plot as well. The radial cutoff and number of grids on 3 dimensions are all specified with `rcut` and `u1`.

- **hb_stat**<a id=hb_stat></a>
The analysis calculates distribution of HB number, length and lifetime etc.

- **hb_stat2**<a id=hb_stat2></a>
The analysis calculates the O-O-O angular distribution of water molecules with different number of acceptence/donation number.

- **hyper**<a id=hyper></a>
The analysis calculates the planarity distribution of hydroxide accepted water molecules. See the definition in Chen et al. 10, 413, Nat. Chem. (2018)

- **incremental_pdf**<a id=incremental_pdf></a>
The analysis calculates the pdf of water molecules in order of distance or topological neighbors (order of hydrogen bond). INPUT parameter `nshell` is needed to specify the number of shells considered.

- **mj** <a id=mj></a>
The analysis takes in `trans.dat` file from `hbs` analysis containing PT time and ion indexes and classify the PTs into different classes including `single`, `double`, `triple`, `quadraple` and `rattle`. For detailed classification standard and discussions refer to Chen et al. 10, 413, Nat. Chem. (2018) and Liu et al. J. Chem. Phys. 157, 024503 (2022).

- **movie/movie2/OH_movie** <a id=movie/movie2/OH_movie></a>
These three analyses are all designed to search certain configurations in the snapshots and output them in `xyz` format for VESTA plot. Please refer to the codes for exactly what configurations are extracted.

- **reorganize** <a id=reorganize></a>
The analysis read in atomic position file in one format and output the positions in `qe` and `xyz` format.

- **ssf** <a id=ssf></a>
The analysis calculates static structure factor (SSF) of single element. For water it only calculates SSF between oxygen atoms $S_{OO}(k)$. The following parameters are also needed for the calculation:
    - `struf_dgx/struf_dgy/struf_dgz` The parameters specify $\Delta k$ on three directions. Typically it has to be an integer multiple of $2\pi/\mathrm{celldm}$, where $\mathrm{celldm}$ is the cell length on corresponding direction. 
    - `struf_ng` Number of `dg` on each direction.
    - `ssf_out` SSF output file name.

- **tetra_order** <a id=tetra_order></a>
The analysis calculates the tetrahedral order of water and its closest neighbors. See, for example, Errington, J. et al. 409, 318. Nat. (2001) for a precise definition.

- **wannier** <a id=wannier></a>
The analysis calculates the distribution of distance from Wannier centers to oxygen atom (stored in `distribution_MLWF.dat`) and distribution of dipole moment of water molecules (stored in `distribution_dipole.dat`).
Again, input parameters `wannier_file` and `nband` are needed in the `INPUT` file. 

- **orientation_tcf** <a id=orientation_tcf></a>
The analysis calculates the 2nd order angular time correlation function (TCF) of water molecule or ion. For a detailed definition and example see Liu et al. J. Chem. Phys. 157, 024503 (2022).

- **msd** <a id=msd></a>
The analysis calculates the mean square displacement (MSD) by averaging square displacement over all water molecules. By Einstein's relation, MSD of molecule in liquid is linearly dependent on time by which we could deduce the diffusivity. The following parameters are needed for calculation of MSD:
    - `msd_dt` Time step between two printed snapshots.

- **msd_multiple** <a id=msd_multiple></a>
The analysis calculates the MSD by averaging over both water molecules and different time intervals. The trajectory is split into multiple segments with a length of `msd_t` ps for each segment and the initial snapshots were displaced by an interval of `msd_dt0` ps. MSDs are calculated for each time segment and averaged to render the final MSD. The method is mainly used for water ions, since in AIMD there are usually only one water ion and time average is necessary to generate MSD linearly dependent on time. See more details in Liu et al. J. Chem. Phys. 157, 024503 (2022). Except for `msd_dt` and `msd_dt0`, `msd_dt` is also needed in the `INPUT` file.

- **special_msd**<a id=special_msd></a>
The analysis calculates and outputs square displacement of each single water molecule.

- **hb_stat3/hb_stat4/hbs_near_angularjump/hbs_near_doubledonor**<a id=angular_jump></a>
These four programs focus on the analyses of angular jump of hydrogen bonds (see Laage, D. 311, 832. Sci. (2006)). However, the programs are not well organized yet and might contain bugs.

- **hb_correlation/hb_correlation2**<a id=hb_correlation/hb_correlation2></a>
The two analyses calculate the hydrogen bond time correlation function in two different ways. For a definition of HB TCF, see Luzar et al., 379, 55, Nature (1996). The following parameters are needed for the calculation:
    - `msd_dt` time interval between two printed snapshots;
    - `msd_t` maximum time length of the TCF;
    - `nPT` maximum number of recorded HB combinations.

- **nonhb_correlation/nonhb_correlation2/nonhb_correlation3**<a id=nonhb_correlation/nonhb_correlation2/nonhb_correlation3></a>
`nonhb_correlation` analysis calculate the TCF of two water molecules staying together within a certain distance. `nonhb_correlation2` and `nonhb_correlation3` calculate the TCF of water molecule converting between H-bonded water molecule and non-H-bonded water molecule.

## Analysis program template
Before adding new analysis subprogram to Candela, we recommand you to go through the existing analysis program list in case the required analysis is already fulfilled. This part introduces the basic structure of an analysis subroutine, by learning which you could write new analysis subroutines on your own.
### Basic classes in Candela
The following table summarizes the basic classes and their properties and methods that are most used in Candela. There are mainly four kinds of classes: `Vector3<T>`, `Cell`, `Atom` and `Water`. `Cellfile` is a class inherited from `Cell` and mainly contains methods reading data from MD trajectory in different formats. All of the classes are straightforwardly defined and directly correspond to physical entities.

 Class Name (source file)  | Properties  | Property description      |Methods | Method description                                                      | 
| :---------------- | :--------------------- | :------------------------------- | :------------- |:---------------------------- |
| **`Vector3<T>`** (`vec3.h`) | `T x, y, z` | component of the vector on three directions | `Vector<T> +-*/ T` | $(x_1+-\times/ T, y_1+-\times/T, z_1+-\times/T)$
| | | | `Vector cos/sin (Vector)` | $(\cos(x), \cos(y), \cos(z))$ or $(\sin(x), \sin(y), \sin(z))$
| | | | `T norm (Vector)` | $\sqrt{x^2+y^2+z^2}$
| | | | `T norm2 (Vector)` | $x^2+y^2+z^2$
| | | | `void normalize (Vector)` | $(x, y, z)/\sqrt{x^2+y^2+z^2}$
| | | | `Vector +-*^/ Vector` | `^` stands for outer product
| | | | `T dot(Vector, Vector)` | $x_1\cdot x_2+y_1\cdot y_2+z_1\cdot z_2$
| | | | `Vector cross(Vector, Vector)` | $(y_1\cdot z_2-y_2\cdot z_1, z_1\cdot x_2-x_1\cdot z_2, x_1\cdot y_2-y_1\cdot x_2)$
| `Atom` (atoms.h/atoms.cpp) | `string id` | name of the element | `void read_pos/read_pos2/read_pos3/read_pos4` | read in positions of atom. The input arguements are not listed in detailed here.
| | `int na` | number of atoms of the element | `void read_vel` |  read in velocities of atom, only for cp.x (`QE`) format.
| | `Vector3<double> * pos` | Cartesian positions of atom | |  
| | `Vector3<double> * posd` | direct positions of atom | |  
| | `Vector3<double> * vel` | velocities of atom| |  
| | `double mass` | relative atomic mass | |  
| | `double charge` | charge | | 
| `Cell` (cell.h/cell.cpp) | `int nat` | total number of atoms in the cell | `void direct2cartesian` | convert positions in direct coordination into Cartesian coordination
| | `int ntypes` |  number of elements| `void cartesian2direct` | convert positions in Cartesian coordination into direct coordination
| | `Atom* atom` | atoms in the cell | `void read_wannier_centers` | read in Wannier center positions. Now only cp.x (`QE`) format is supported.
| | `Vector3<double>* wan_centers` | Wannier centers in the cell | `void read_eig ()` | read in eigen value of the solved KS equation
| | `double snapshot_time` | simulation time of the current snapshot | | 
| | `int snapshot_index` | simulation index of the current snapshot | |
| | `double* eig` | eigen value of the solved KS equation | |
| `Cell CellFile` (`cellFile.h/cpp, cellFile_PROFESS/VASP/QE/QE2.cpp, etc.`)| `static ifstream ifs_pos_kept` | ifstream opening the file containing atomic positions | `static bool ReadGeometry()` | read in atomic positions|
| | `static ifstream ifs_wan_kept` | ifstream opening the file containing Wannier center positions | `static bool ReadGeometry_PROFESS/VASP/QE/QE2 ...` | read in atomic positions in different formats
| | `static ifstream ifs_cel_kept` | ifstream opening the file containing cell arrays. Only useful for cp.x (`QE`) format. | `static bool WriteGeometry()` | write atomic positions in certain format. Now cp.x (`QE`) and xyz `XYZ` are supported.
| | `static ifstream ifs_eig_kept` | ifstream opening the file containing eigenvalues | `static bool ReadVelocity()` | read in atomic velocities. Only useful for cp.x (`QE`) format.
| | `static ifstream ifs_vel_kept` | ifstream opening the file containing atomic velocities | |
| `Water` (`water.h/cpp, HBs.h/cpp`) | `static int nions` | number of ions in the cell | `void HBs::setup_water()` | set up the water array in the cell |
||`int indexO`| the index is usually the same as the index of the water||
||`int naccept`|number of accepted water molecules ||
||`int ndonate`|number of donated water molecules ||
||`int* indexH`|H index of the water molecule ||
||`double* disH`|distance of H atoms from O atom in Angstroms||
||`int* acceptO`|index of accepted O atoms||
||`int* acceptH`|index of accepted H atoms||
||`double* accept_angle`|HOO angle of the accepted HB||
||`double* accept_disO`|distance of the accepted O atom from the central O atom||
||`int* donateO`|index of donated O atoms||
||`int* donateH`|index of donated H atoms||
||`double* donate_angle`|HOO angle of the donated HB||
||`double* donate_disO`|distance of the donated O atom from the central O atom||
||`double dipole[3]`|dipole moment of the water molecule||

Generally, the relationship between the classes is: `Cell` contains different kinds of `Atom`, and `Atom` could make up to form `Water`. The properties of these classes in vector form are stored in `Vector3<T>` instances.

### The structure of a typical analysis C++ file
Here we give a brief description of `example.cpp` contained in `/examples/example` to show how a typical analysis could be done under the frame of Candela.

First, the class `Example` and its properties and methods are defined in `example.h` as follows:

```cpp
#ifndef EXAMPLE_H
#define EXAMPLE_H
#include "Cellfile.h"
class Example
{
    public:
    Example();
    ~Example();

    void Routine();
    void calc(CellFile &cel);
    void output();

    double* example_arr;
    int example_len;
};

#endif
```

The function `Routine()` is to be called in the `main.cpp` to implement the analysis. Function `calc()` takes in all information of a cell of a single snapshot and do the calculations. Function `output()` writes out the results in files.

In `example.cpp`, function `Routine()` is implemented as the follows:

```cpp
void Example::Routine()
{
    // Initialize necessary arrays.
    this->example_arr = new double[this->example_len]; 

    int count_geometry_number = 0;
    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

        if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean(); // renxi added 20200614
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;
        this->calc(cel);
        cel.clean();
	}//igeo

    this->output();
    delete[] this->example_arr;
    return;
}
```

In the function, it first initialize the memory of array used in this analysis. Then, it iteratively do the following things to each snapshots in the MD trajectory:
1. judge whether the snapshot has to be analyzed (`if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) `). If so, label the snapshot as non-ignorable (`cel.read_and_used=true;`);
2. read in the atomic positions and velocities, Wannier centers etc. if necessary (`if( !CellFile::ReadGeometry( cel ) ) continue;`);
3. judge whether the snapshot is ignorable, if so, clean the cell and move on to the next iteration (`cel.clean(); continue;`); else, go to step 4;
4. do analysis to the cell of the snapshot (`this->calc(cel);`) and clean the cell.

Function `calc(Cellfile &cel)` is implemented as the follows:

```cpp

void Example::calc(CellFile &cel)
{
    int ito=-1;
	int ith=-1;
	int itc=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
	}
	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
    
    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);
    for (int iwater1=0; iwater1<cel.atom[ito].na; iwater1++)
    {
        // The analysis begins here...
    }
    return;
}
```

In the function, each element is first assigned with an index (`ito`, `ith`, etc.). Then memories are allocated for an array of `water` class instances, corresponding to the water molecules in the current snapshot cell (`Water *water = new Water[cel.atom[ito].na];`). The water array is fulfilled later (`HBs::setup_water(cel, water);`), and you could do your calculation to each of the water molecules.

And of course, you have to add your analysis in `main.cpp` as well to let it run.

That's all we would like to say in this documentation. Please feel free to contact us for any help. Happy using Candela and coding!