#include "wfFile.h"
#include "input.h"

bool WfFile::wf_open=false;
FILE* WfFile::fp_kept;
WfFile::WfFile()
{
}

WfFile::~WfFile()
{
}

bool WfFile::ReadWf()
{
	TITLE("WfFile","ReadWf");


	if(INPUT.wf_in_type=="PWmat")
    {
        if(wf_open == false)
        {
            stringstream ss1;
            ss1 << INPUT.wf_directory;
            cout << " ReadWaveFunction: " << ss1.str() << endl;
            wf_kept=open(ss1.str().c_str(),"r");
            stringstream ss2;
            ss2 << INPUT.wf_directory;
            cout << " ReadWaveFunction: " << ss2.str() << endl;
            wf_kept=open(ss2.str().c_str(),"r");
            stringstream ss3;
            ss3 << INPUT.wf_directory;
            cout << " ReadWaveFunction: " << ss3.str() << endl;
            wf_kept=open(ss3.str().c_str(),"r");
			
            if(!fp_kept)
            {
              	cout << "Could not open the wave function file." << endl;
                exit(0);
            }
            wf_open = true;
        }
        return ReadGeometry_PWmat(cel,fp_kept);
    }
    /*else if(INPUT.geo_in_type=="QE")
    {
        if(wf_open == false)
        {
            stringstream ss;
            ss << INPUT.geo_directory;
                        cout << " ReadGeometry : " << ss.str() << endl;
                        fp_kept.open(ss.str().c_str());
                        cout << " File name is " << ss.str() << endl;
                        cout << " WfFile::wf_open = " << wf_open << endl;
                        if(!fp_kept)
                        {
                                cout << "Could not open the file." << endl;
                                exit(0);
                        }
            wf_open = true;
        }
        return ReadGeometry_QE2(cel,fp_kept);
    }*/
	else
	{
		cout << " Warning! We can only accept PWmat/QE input Wavefunction." << endl;
		cout << " geo_in_type = " << INPUT.geo_in_type << endl;
		exit(0);
	}
	return false;
}	



