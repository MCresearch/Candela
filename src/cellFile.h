#ifndef CELLFILE_H
#define CELLFILE_H

#include "gfun.h"
#include "vec3.h"
#include "cell.h"

// each file containing a set of type of atoms.
class CellFile : public Cell 
{
	public: 
	
	CellFile();
	~CellFile();

	static bool ReadGeometry( Cell &cel );
	static bool CheckGeometry( Cell &cel );
	static void WriteGeometry( Cell &cel, bool cartesian=false );

	static bool ReadVelocity( Cell &cel );

	private:

	// set this function as static because
	// not each cell need to read the geometry
	static bool ReadGeometry_PROFESS( Cell &cel );
	static bool ReadGeometry_VASP( Cell &cel );
	static bool ReadGeometry_QE( Cell &cel, ifstream &ifs, ifstream &ifs_cel, ifstream &ifs_wan, 
								ifstream &ifs_eig, ifstream &ifs_pos_ili, ifstream &ifs_vel);
	static bool ReadGeometry_QE2( Cell &cel, ifstream &ifs );//qianrui add 2020-1-6
	static bool ReadGeometry_PWmat( Cell &cel, ifstream &ifs );//qianrui add 2020-1-6
	static bool ReadGeometry_XYZ( Cell &cel, ifstream &ifs);
	static bool ReadGeometry_LAMMPS( Cell &cel, ifstream &ifs );
	//static bool ReadGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24
	static bool ReadGeometry_ABACUS( Cell &cel, ifstream &ifs);
	static bool ReadGeometry_RAW( Cell &cel, ifstream &ifs ); // mohan add 2019-03-17

	// check geometry
	static bool CheckGeometry_PROFESS( Cell &cel );
	static bool CheckGeometry_VASP( Cell &cel );
	static bool CheckGeometry_QE( Cell &cel );
	static bool CheckGeometry_XYZ( Cell &cel );
	static bool CheckGeometry_LAMMPS( Cell &cel );
	static bool CheckGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24
	static bool CheckGeometry_RAW( Cell &cel ); // mohan add 2019-03-17

	// write geometry
	static void WriteGeometry_PROFESS( Cell &cel, bool cartesian );
	static void WriteGeometry_VASP( Cell &cel, bool cartesian );
	static void WriteGeometry_QE( Cell &cel );
	static void WriteGeometry_XYZ( Cell &cel );
	static void WriteGeometry_LAMMPS( Cell &cel ); // mohan add 2013-09-20
	static void WriteGeometry_ABACUS( Cell &cel ); // mohan add 2015-07-24
	static void WriteGeometry_RAW( Cell &cel ); // mohan add 2019-03-17

	// read velocity
	static bool ReadVelocity_PROFESS( Cell &cel );
	static bool ReadVelocity_VASP( Cell &cel );

	private:

	static bool file_open;
	static ifstream ifs_pos_kept;
	static ifstream ifs_cel_kept;
	static ifstream ifs_wan_kept;
	static ifstream ifs_eig_kept;
	static ifstream ifs_pos_ili_kept;
	static ifstream ifs_vel_kept;

};


#endif
