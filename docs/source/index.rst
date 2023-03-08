.. Candela documentation master file, created by
   sphinx-quickstart on Mon Mar  6 21:43:38 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
Welcome to Candela's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Overview
   :ref:`Overview`

.. toctree::
   :maxdepth: 3
   :caption: Basic usage and variables
   :ref:`Basic-usage-and-variables`
   
   :ref:`Compilation`
   :ref:`Using-Candela`
   :ref:`List-of-supported-analyses`

.. toctree::
   :maxdepth: 2
   :caption: Develop guide
   :ref:`Develop-guide`

   :ref:`Analysis-program-template`
   :ref:`structure-of-file`
   :ref:`Adding-test`

-------------------------
.. _Overview:
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

===================================
.. _Basic-usage-and-vars:

Basic usage and variables
===================================

-----------
.. _Compilation:

Compilation
-----------

Type `make` in the `Candela` directory to compile the executable `Candela.exe` in `\bin` directory. 
The default compiler is g++, which could be switched to other compiliers by changing `CXX` in `Makefile`. 

-------------
.. _Using-Candela:

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


===================================
.. _List-of-supported-analyses:

List of supported analyses
===================================

--------------------------------
.. _PDF:

PDF (Pair distribution function)
--------------------------------

Pair distribution function, also radial distribution function. 
The following parameters are needed for `pdf` calculation:

    - `rcut` Radius cutoff of pair distribution function calculation.

    - `dr` Grid length of radius.

    - `id1` Name of the first type of atom of the pdf. For example, for :math:`g_{OH}(r)`, `id1` is O.

    - `id2` Name of the second type of atom of the pdf. For example, for :math:`g_{OH}(r)`, `id2` is H.


--------------------------------
.. _SSF:

SSF (Static structure factor)
--------------------------------

The analysis calculates static structure factor (SSF) of single element. For water it only calculates SSF between oxygen atoms :math:`S_{OO}(k)`. 
The following parameters are also needed for the calculation:
   
    - `struf_dgx/struf_dgy/struf_dgz` The parameters specify :math:`\Delta k` on three directions. Typically it has to be an integer multiple of :math:`2\pi/\mathrm{celldm}`, where :math:`\mathrm{celldm}` is the cell length on corresponding direction. 
    
    - `struf_ng` Number of `dg` on each direction.
    
    - `ssf_out` SSF output file name.


--------------------------------
.. _MSD:

MSD (Mean square displacement)
--------------------------------

The analysis calculates the mean square displacement (MSD) by averaging square displacement over all atoms. 
By Einstein's relation, MSD of molecule in liquid is linearly dependent on time by which we could deduce the diffusivity. 
For system with more than one type of element, the program will treat every type of element as one element. 
See `003_MSD` for example of `msd` calculation.
The following parameters are needed for calculation of MSD:
   
    - `msd_dt` Time step between two printed snapshots.


===================================
.. _Develop-guide:

Developing guide
===================================

Before adding new analysis subprogram to Candela, 
we recommand you to go through the existing analysis subprogram list in case the needed analysis is already fulfilled. 
This part introduces the basic structure of an analysis subroutine, 
by learning which you could write new analysis subroutines on your own.


--------------------------
.. _Analysis-program-template:

Analysis program template
--------------------------

To write new subprograms in Candela, it is necessary to understand and take use of the basic classes implemented. 
They will greatly help you in developing your own codes.
The following table summarizes the basic classes and their properties and methods that are most used in Candela. 
There are mainly four kinds of classes: :code:`Vector3<T>`, :code:`Cell`, :code:`Atom` and :code:`Water`. 
`Cellfile` is a class inherited from :code:`Cell` and mainly contains methods reading data from MD trajectory in different formats. 
All of the classes are straightforwardly defined and directly correspond to physical entities.

+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|Class Name (source file)   | Properties  |Property description              |Methods                        | Method description                                |
|                           |             |                                  |                               |                                                   |
+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|:code:`Vector3<T>`         |`T x, y, z`  |component of the vector           |:code:`Vector<T> +-*/ T`       |:math:`+-\times \div`                              |
|(`vec3.h`)                 |             |on three directions               +-------------------------------+---------------------------------------------------+
|                           |             |                                  |:code:`Vector cos/sin (Vector)`|:math:`\sin/ \cos`                                 |
|                           |             |                                  +-------------------------------+---------------------------------------------------+
|                           |             |                                  |                               |                                                   |
|                           |             |                                  |:code:`T norm (Vector)`        |:math:`x^2+y^2+z^2`                                |
|                           |             |                                  +-------------------------------+---------------------------------------------------+
|                           |             |                                  |                               |                                                   |
|                           |             |                                  |:code:`Vector +-*^/ Vector`    |^ for outer product                                |
|                           |             |                                  +-------------------------------+---------------------------------------------------+
|                           |             |                                  |                               |                                                   |
|                           |             |                                  |:code:`T dot(Vector, Vector)`  |inner product                                      |
|                           |             |                                  +-------------------------------+---------------------------------------------------+
|                           |             |                                  |                               |                                                   |
|                           |             |                                  |:code:`T cross(Vector, Vector)`|outer product                                      |
+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|:code:`Atom`               |`string id`  |name of the element               |:code:`void read_pos/read_pos2`|read in positions of atom.                         |
|(`atoms.h/atoms.cpp`)      |             |                                  |:code:`/read_pos3/read_pos4`   |Input arguements are not listed here.              |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`int na`     |number of atoms of the element    |:code:`void read_vel`          |read in velocities of atom. Currently only         |
|                           |             |                                  |                               |available in `QE`.                                 |  
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`Vector3<dou`|positions of atoms                |                               |                                                   |
|                           |`ble>*pos`   |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`Vector3<dou`|direct positions of atoms         |                               |                                                   |
|                           |`ble>*posd`  |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`Vector3<dou`|velocities of atoms               |                               |                                                   |
|                           |`ble>*vel`   |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`double mass`|mass of atom                      |                               |                                                   |
|                           |             |                                  |                               |                                                   |
+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|:code:`Cell`               |`int na`     |number of atoms of the element    |:code:`void direct2cartesian`  |convert positions in direct coordination into      |
|(`cell.h/cell.cpp`)        |             |                                  |                               |Cartesian coordination                             |  
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`int ntypes` |number of elements                |:code:`void cartesian2direct`  |convert positions in Cartesian coordination into   |
|                           |             |                                  |                               |direct coordination                                |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`Atom* atom` |atoms in the cell                 |:code:`void read_wannier`      |read in Wannier center positions.                  |
|                           |             |                                  |:code:`_centers`               |Now only cp.x (`QE`) format is supported.          |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`Vector3`    |Wannier centers in the cell       |:code:`void read_eig ()`       |read in eigen value of the solved KS equation.     |
|                           |`<double>*wa`|                                  |                               |Now only pw.x (`QE2`) format is supported.         |
|                           |`n_centers`  |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`double `    |simulation time of the            |                               |                                                   |
|                           |`snapshot`   |current snapshot                  |                               |                                                   |
|                           |`_time`      |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`double `    |simulation snapshot index of the  |                               |                                                   |
|                           |`snapshot`   |current snapshot                  |                               |                                                   |
|                           |`_index`     |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`double* eig`|eigen value of the solved         |                               |                                                   |
|                           |             |KS equation                       |                               |                                                   |
|                           |             |                                  |                               |                                                   |
+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|:code:`Cell CellFile`      |`static ifst`|ifstream opening the file         |:code:`static`                 |read in atomic positions                           |
|(`cellFile.h/cpp`)         |`ream ifs_`  |containing atomic positions       |:code:`bool ReadGeometry()`    |                                                   |
|                           |`pos_kept`   |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`static ifst`|ifstream opening the file         |:code:`static bool`            |read in atomic positions of different type         |
|                           |`ream ifs_`  |containing Wannier centers        |:code:`ReadGeometry_PROFESS`   |                                                   |
|                           |`wan_kept`   |                                  |:code:`/VASP/QE/QE2 ...`       |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`static ifst`|ifstream opening the file         |                               |                                                   |
|                           |`ream ifs_`  |containing lattice vectors        |                               |                                                   |
|                           |`cel_kept`   |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`static ifst`|ifstream opening the file         |                               |                                                   |
|                           |`ream ifs_`  |containing eigenvalues            |                               |                                                   |
|                           |`eig_kept`   |                                  |                               |                                                   |
|                           +-------------+----------------------------------+-------------------------------+---------------------------------------------------+
|                           |`static ifst`|ifstream opening the file         |                               |                                                   |
|                           |`ream ifs_`  |containing velocities of atom     |                               |                                                   |
|                           |`vel_kept`   |                                  |                               |                                                   |
+---------------------------+-------------+----------------------------------+-------------------------------+---------------------------------------------------+

--------------------------------------------
: _structure-of-file:

The structure of a typical analysis C++ file
--------------------------------------------

Here we give a brief description of `example.cpp` contained in `/examples/example` to show how a typical analysis of water system could be done under the frame of Candela.

First, the class `Example` and its properties and methods are defined in `example.h` as follows:

.. code-block :: cpp

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

The function `Routine()` is to be called in the `main.cpp` to implement the analysis. Function `calc()` takes in all information of a cell of a single snapshot and do the calculations. Function `output()` writes out the results in files.

In `example.cpp`, function `Routine()` is realized as the follows:

.. code-block :: cpp

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

The function first initialize the memory of array used in this analysis. Then, it iteratively do the following things to each snapshots in the MD trajectory:
   - judge whether the snapshot has to be analyzed (`if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) `). If so, label the snapshot as non-ignorable (`cel.read_and_used=true;`);
   - read in the atomic positions and velocities, Wannier centers etc. if necessary (`if( !CellFile::ReadGeometry( cel ) ) continue;`);
   - judge whether the snapshot is ignorable, if so, clean the cell and move on to the next iteration (`cel.clean(); continue;`); else, go to step 4;
   - do some calculations to the cell of the snapshot (`this->calc(cel);`) and clean the cell.

And of course, you have to add your analysis to `main.cpp` as well to let it run.

--------------------------------------------
: _Adding-test:

Adding test for the subprogram
--------------------------------------------

After adding new functions into CANDELA, we must add corresponding tests.
Directory **test** contains:

**001_PDF**,...: test cases

**Autotest.sh**: shell script to run auto tests.

**tool**: some tool functions

**geo**: contains geometry files

--------------------------------------------
Preparation for one case
--------------------------------------------

**compare**:

.. code-block:: sh

   threshold 1e-4

   enable_mpi ON
   
   pdf.ref  pdf.txt
   
   other.ref other.result
   
   ...

threshold: When the relative errors between reference results and calculating results are larger than threshold, the tests will fail.

enable_mpi: Whether this case supports MPI or not.

The remaining lines are reference files and result files. Each pair is in one line.

**INPUT**:

.. code-block:: sh

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

**pdf.ref / other.ref / ...**:
The reference results are in them.

--------------------------------------------
Run auto test
--------------------------------------------

**method 1**:

.. code-block::

   make CXX=g++ TEST=ON # serial 
   make CXX=mpicxx TEST=ON # parallel

Autotest.sh will be executed after serial-version **candela** is compiled.


`make CXX=mpicxx TEST=ON`

Autotest.sh will be executed after parallel-version **candela** is compiled. Multi-processor cases will be tested at the same time.

**All developers should run 

`make CXX=g++ TEST=ON`
and `make CXX=mpicxx TEST=ON`

to make sure some functions are not be destroyed.**

**method 2**:

After **candela** is compiled.

`make test` or `cd test;sh Autotest.sh`

It will run one-processor tests.

`make test CXX=mpicxx` or `cd test;sh Autotest.sh ON`

It will run one-processor and multi-processor tests.




That's all we would like to say in this documentation. Please feel free to contact us for any help. 

Happy using Candela and coding!
