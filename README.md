[User's Manual](Doc/manual.md)

[Developer's Manual](Doc/develop.md)

- [Overview](#overview)
- [Basic usage and variables](#basic-usage-and-variables)
    - [Compilation](#compilation)
    - [Using Candela](#using-candela)
- [List of supported analyses](#list-of-supported-analyses)
    - [PDF](#pdf)
    - [SSF](#ssf)
    - [MSD](#msd)

## Overview
Candela is short for Collection of ANalysis DEsigned for Large-scale Atomic simulations. It is developed by [MCresearch](https://github.com/MCresearch) to conduct analyses on MD trajectory in different formats. The analyses herein are mostly targeted at liquid water and warm condensed matter system. 

Section [Basic usage and variables](#Basic-usage-and-variables) introduces the usage of Candela;

Section [List of supported analyses](#List-of-supported-analyses) introduces the implemented analyses in Candela for now;

Section [Analysis program template](#Analysis-program-template) introduces the structure of a typical analysis subprogram, in case you may need to add new analysis program to Candela.

## Basic usage and variables 
### Compilation
Type `make` in the `Candela` directory to compile the executable `Candela.exe` in `\bin` directory. The default compiler is g++, which could be switched to other compiliers by changing `CC` in `Makefile`. 

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

In each line the first item specifies the name of the parameter, and you may fill in the value of the parameter after it. The following chart displays the general parameters and their purposes.

 Name  | Type  | Example                                                      | Discription                                                      |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| **calculation** | String | pdf | Name of the analysis
| **system** | String | hydronium | Name of the system
|  **geo_in_type** | String  | QE | Format of the MD trajectory file
| **geo_directory** | String | ../H_PBE.pos | directory of the MD trajectory file
| **geo_1** | Int | 1 | Starting index of snapshot in the MD trajectory file.  It is recommended to set it to 1.
| **geo_2** | Int | 10 | Final index of snapshot in the MD trajectory file
| **geo_interval** | Int | 1 | The program analyse 1 snapshot for every `geo_interval` snapshots
| **ntype** | Int | 2 | Number of atom species. It depends on the MD format. The maximum number of spicies is mostly 3 or 4.
| **natom** | Int | 190 | Number of atoms
| **natom1** | Int | 63 | Number of atoms of the 1st species
| **natom2** | Int | 127 | Number of atoms of the 2nd species
| **id1** | String | O | Element name of the 1st species of atoms
| **id2** | String | H | Element name of the 2nd species of atoms
| **celldm1** | float | 12.444655509 | Length of the cell in the 1st dimension
| **celldm2** | float | 12.444655509 | Length of the cell in the 2nd dimension
| **celldm3** | float | 12.444655509 | Length of the cell in the 3rd dimension
| **dt_snapshots** | float | 4.83776865e-5 | A single MD timestep of `cp.x` in ps. The value is needed for trajectory over 100 ps.

For `geo_in_type`, currently supported formats include ABACUS (`ABACUS`), pw.x (`QE2`), cp.x (`QE`), VASP (`VASP`), PWmat (`PWMAT`) and so on.

The physical quantities used in the program are mostly in Angstrom/ps/eV unit.

Other parameters are needed for specific `calculation`, which you may refer to [List of supported analyses](#List-of-supported-analyses) where a detailed summary of suppoted `calculation` of analyses is given. You may refer to the example `INPUT` files in `/examples` as well for different kinds of analyses.

## List of supported analyses

- **pdf** <a id="pdf"></a>
Pair distribution function, also radial distribution function. See `test/001_PDF` and `test/001_PDF2` for examples of `pdf` calculation. The following parameters are needed for `pdf` calculation:

    - `rcut` Radius cutoff of pair distribution function calculation.

    - `dr` Grid length of radius.

    - `id1` Name of the first type of atom of the pdf. For example, for $g_{OH}(r)$, `id1` is O.

    - `id2` Name of the second type of atom of the pdf. For example, for $g_{OH}(r)$, `id2` is H.

- **ssf** <a id=ssf></a>
The analysis calculates static structure factor (SSF) of single element. For water it only calculates SSF between oxygen atoms $S_{OO}(k)$. The following parameters are also needed for the calculation:
    - `struf_dgx/struf_dgy/struf_dgz` The parameters specify $\Delta k$ on three directions. Typically it has to be an integer multiple of $2\pi/\mathrm{celldm}$, where $\mathrm{celldm}$ is the cell length on corresponding direction. 
    - `struf_ng` Number of `dg` on each direction.
    - `ssf_out` SSF output file name.


- **msd** <a id=msd></a>
The analysis calculates the mean square displacement (MSD) by averaging square displacement over all water molecules. By Einstein's relation, MSD of molecule in liquid is linearly dependent on time by which we could deduce the diffusivity. For system with more than one type of element, the program will treat every type of element as one element. See `003_MSD` for example of `msd` calculation. The following parameters are needed for calculation of MSD:
    - `msd_dt` Time step between two printed snapshots.