.. Candela documentation master file, created by
   sphinx-quickstart on Mon Mar  6 21:43:38 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
Welcome to Candela's documentation!
===================================

-------------------------
Overview
-------------------------

Candela is short for Collection of ANalysis DEsigned for Large-scale Atomic simulations. 
It is developed by `MCresearch <https://github.com/MCresearch>`_ to conduct analyses on MD trajectory in different formats. 
Right now the program only supports analysis of pair distribution function (PDF), static structure factor (SSF) and mean square displacement (MSD).

Section :ref:`Basic-usage-and-vars` introduces the usage of Candela;

Section :ref:`List-of-supported-analyses` introduces the implemented analyses in Candela for now;

Section :ref:`Develop-guide` introduces the structure of a typical analysis subprogram, 
in case you may need to add new analysis program to Candela.
The test method is also introduced.

.. _Basic-usage-and-vars:

Basic usage and variables
===================================

Compilation
-----------

Type `make` in the `Candela` directory to compile the executable `Candela.exe` in `\bin` directory. 
The default compiler is g++, which could be switched to other compiliers by changing `CC` in `Makefile`. 

Using Candela
-------------

First, you have to prepare an `INPUT` file containing the input parameters of Candela. 
The following shows an example of the `INPUT` file::


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


In each line the first item specifies the name of the parameter, 
and you may fill in the value of the parameter after it. 
The following chart displays the general parameters and their purposes.

+---------------+-----------------+------------------+--------------------------------+
|Name           |       Type      |    Example       |   Discription                  |
+---------------+-----------------+------------------+--------------------------------+
| *calculation* |    String       |   pdf            | Name of analysis               |
+---------------+-----------------+------------------+--------------------------------+
| *system*      |    String       |   water          |Name of the system              |
+---------------+-----------------+------------------+--------------------------------+
| *geo_in_type* |    String       |   QE             |directory of the                |
|               |                 |                  |MD trajectory file              |
+---------------+-----------------+------------------+--------------------------------+
| *geo_1*       |    Int          |   1              |Starting index of snapshot      |
|               |                 |                  |in the MD trajectory file.      |
|               |                 |                  |In most cases it has to be 1.   |
+---------------+-----------------+------------------+--------------------------------+
| *geo_2*       |    Int          |   10             |End index of snapshot           |
|               |                 |                  |in the MD trajectory file.      |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *geo_interval*|    Int          |   1              | The program analyse 1 snapshot |
|               |                 |                  | for every `geo_interval`       |
|               |                 |                  | snapshots                      |
+---------------+-----------------+------------------+--------------------------------+
| *ntype*       |    Int          |   2              | Number of atom species.        |
|               |                 |                  |                                |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *natom*       |    Int          |   192            | Number of atoms.               |
|               |                 |                  |                                |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *natom1*      |    Int          |   64             | Number of atoms of type 1.     |
|               |                 |                  |                                |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *natom2*      |    Int          |   128            | Number of atoms of type 2.     |
|               |                 |                  |                                |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *id1*         |    String       |   O              | Element name of the            |
|               |                 |                  | 1st species of atoms.          |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *id2*         |    String       |   H              | Element name of the            |
|               |                 |                  | 2nd species of atoms.          |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *celldm1*     |    Float        |   12.444655509   | Length of the cell             |
|               |                 |                  | in the 1st dimension.          |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *celldm2*     |    Float        |   12.444655509   | Length of the cell             |
|               |                 |                  | in the 2nd dimension.          |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *celldm3*     |    Float        |   12.444655509   | Length of the cell             |
|               |                 |                  | in the 3rd dimension.          |
|               |                 |                  |                                |
+---------------+-----------------+------------------+--------------------------------+
| *dt_snapshots*|    Float        |  4.83776865e-5   |A single MD timestep of `cp.x`  |
|               |                 |                  |in ps. The value is needed      |
|               |                 |                  |for trajectory over 100 ps.     |
+---------------+-----------------+------------------+--------------------------------+

For `geo_in_type`, currently supported formats include pw.x (`QE2`), 
cp.x (`QE`), VASP (`VASP`), ABACUS (`ABACUS`), PWmat (`PWMAT`) and so on.

The physical quantities used in the program are all in Angstrom/ps/eV unit.

Other parameters are needed for specific `calculation`, 
which you may refer to :ref:`List-of-supported-analyses` where a detailed summary of supported `calculation` of analyses is given. 
You may refer to the example `INPUT` files in `/examples` as well for different kinds of analyses.

.. _List-of-supported-analyses:

List of supported analyses
===================================

PDF (Pair distribution function)
--------------------------------

Pair distribution function, also radial distribution function. 
The following parameters are needed for `pdf` calculation:

    - `rcut` Radius cutoff of pair distribution function calculation.

    - `dr` Grid length of radius.

    - `id1` Name of the first type of atom of the pdf. For example, for :math:`g_{OH}(r)`, `id1` is O.

    - `id2` Name of the second type of atom of the pdf. For example, for :math:`g_{OH}(r)`, `id2` is H.

SSF (Static structure factor)
--------------------------------

The analysis calculates static structure factor (SSF) of single element. For water it only calculates SSF between oxygen atoms :math:`S_{OO}(k)`. 
The following parameters are also needed for the calculation:
   
    - `struf_dgx/struf_dgy/struf_dgz` The parameters specify :math:`\Delta k` on three directions. Typically it has to be an integer multiple of :math:`2\pi/\mathrm{celldm}`, where :math:`\mathrm{celldm}` is the cell length on corresponding direction. 
    
    - `struf_ng` Number of `dg` on each direction.
    
    - `ssf_out` SSF output file name.

MSD (Mean square displacement)
--------------------------------

The analysis calculates the mean square displacement (MSD) by averaging square displacement over all atoms. 
By Einstein's relation, MSD of molecule in liquid is linearly dependent on time by which we could deduce the diffusivity. 
For system with more than one type of element, the program will treat every type of element as one element. 
See `003_MSD` for example of `msd` calculation.
The following parameters are needed for calculation of MSD:
   
    - `msd_dt` Time step between two printed snapshots.

.. _Develop-guide:

Developing guide
===================================

Before adding new analysis subprogram to Candela, 
we recommand you to go through the existing analysis subprogram list in case the needed analysis is already fulfilled. 
This part introduces the basic structure of an analysis subroutine, 
by learning which you could write new analysis subroutines on your own.

Analysis program template
--------------------------

To write new subprograms in Candela, it is necessary to understand and take use of the basic classes implemented. 
They will greatly help you in developing your own codes.
The following table summarizes the basic classes and their properties and methods that are most used in Candela. 
There are mainly four kinds of classes: :code:`Vector3<T>`, `Cell`, `Atom` and `Water`. `Cellfile` is a class inherited from `Cell` and mainly contains methods reading data from MD trajectory in different formats. All of the classes are straightforwardly defined and directly correspond to physical entities.