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
| **celldm1** | float | 12.444655509 | Length of the cell in the 1st dimension
| **celldm2** | float | 12.444655509 | Length of the cell in the 2nd dimension
| **celldm3** | float | 12.444655509 | Length of the cell in the 3rd dimension
| **dt_snapshots** | float | 4.83776865e-5 | A single MD timestep of `cp.x` in ps. The value is needed for trajectory over 100 ps.

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
- **hbs** <a id="hbs"></a> Hydrogen bond (HB) analysis. Multiple properties are calculated, among which the most used ones are HB number statistics (at the end of `running0.log`) and proton transfer (PT) time record (only for `system` is `hydroxide` or `hydronium`). Another important function in this module is 

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
