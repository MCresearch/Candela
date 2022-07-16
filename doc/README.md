# An introduction to Candela
- [Overview](#Overview)
- [Basic usage and variables](#Basic-usage-and-variables)
    - [Compilation](#Compilation)
    - [Using Candela](#Using-Candela)
- [List of supported analyses](#List-of-supported-analyses)
    - [Water-related analyses](#Water-related-analyses)
        - [pdf](#pdf)
- [Analysis program template](#Analysis-template)
## Overview
Candela is short for Collection of ANalysis DEsigned for Large-scale Atomic simulations. It is developed by [MCresearch](https://github.com/MCresearch) to conduct analyses on MD trajectory in different formats. The analyses herein are mostly targeted at liquid water and warm condensed matter system. 

Section [Basic usage and variables](#Basic-usage-and-variables) introduces the usage of Candela;

Section [List of supported analyses](#List-of-supported-analyses) introduces the implemented analyses in Candela for now;

Section [Analysis program template](#Analysis-template) introduces the structure of a typical analysis subprogram, in case you may need to add new analysis program to Candela.

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
| **celldm1** | float | 12.444655509 | Length of the cell in the 1st dimention
| **celldm2** | float | 12.444655509 | Length of the cell in the 2nd dimention
| **celldm3** | float | 12.444655509 | Length of the cell in the 3rd dimention

For `geo_in_type`, currently supported formats include pw.x (`QE2`), cp.x (`QE`), VASP (`VASP`), ABACUS (`ABACUS`), PWmat (`PWMAT`) and so on.

The physical quantities used in the program are all in Angstrom/ps/eV.

Other parameters are needed for specific `calculation`, which you may refer to [List of supported analyses](#List-of-supported-analyses) where a detailed summary of suppoted `calculation` of analyses is given.

## List of supported analyses
### Water-related analyses
- **pdf** <a id="pdf"></a>
Pair distribution function, also radial distribution function. The following parameters are needed for `pdf` calculation:

    - `rcut` Radius cutoff of pair distribution function calculation.

    - `dr` Grid length of radius.

    - `id1` Name of the first type of atom of the pdf. For example, for $g_{OH}(r)$, `id1` is O.

    - `id2` Name of the second type of atom of the pdf. For example, for $g_{OH}(r)$, `id2` is H.
