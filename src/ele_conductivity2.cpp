#include "input.h"
#include "const.h"
#include "gfun.h"
#include "ele_conductivity.h"
#include "matrixmultip.h"

#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR      1.839939223835727e7

/**
 * @brief Calculate conducitities use time correlation funtions
 * 
 * @note In this method smear is set to 1.
 *       n_fwhm, nscf are not used
 *       smearinvw is not used and it equivalent to smearinvw = 0
 *        
 */
void Ele_Conductivity::method2()
{
	double times = 16;

	double wcut = INPUT.wcut;
	double dw_in = INPUT.dw;
	double fwhm_in = INPUT.fwhm[0];
	
	int nw = ceil(wcut / dw_in);
    double dw = dw_in / P_Ry2eV; // converge unit in eV to Ry
    double sigma = fwhm_in / TWOSQRT2LN2 / P_Ry2eV;
    double dt = M_PI / (dw * nw) / times; // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20) / sigma / dt);
    cout << "nw: " << nw << " ; dw: " << dw * P_Ry2eV << " eV" << endl;
    cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << endl;
    assert(nw >= 1);
    assert(nt >= 1);

    double *ct11 = new double[nt];
    double *ct12 = new double[nt];
    double *ct22 = new double[nt];
    ZEROS(ct11, nt);
    ZEROS(ct12, nt);
    ZEROS(ct22, nt);
	Wavefunc WF;
	WfRead wfr(WF);
	wfr.Init();
	int nk=INPUT.nkpoint;
	assert(nk>0);
	int allcount=-1;
	//kpioint loop
	for(int ik=0;ik<nk;ik++)
	{
		allcount++;
#ifdef __MPI
        if((allcount)%NPROC!=RANK)
		{
			wfr.ignore(ik);
			wfr.clean();
			continue;
		}
#endif
		wfr.readWF(ik);
		int nband=WF.nband;
		for(int ib=0;ib<nband;ib++)
		{
			WF.checknorm(ik,ib);//check if wavefunction is normalized to 1.
		}
		std::cout<<"kpoint "<<ik+1<<endl;
		std::cout<<"nband: "<<nband<<endl;
		jjcorr_ks(ik, nt, dt, WF, ct11,ct12,ct22);
		wfr.clean();
	}

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

	if (RANK == 0)
    {
        calcondw(nt, dt, fwhm_in, wcut, dw_in, ct11, ct12, ct22);
    }
	delete[] ct11;
    delete[] ct12;
    delete[] ct22;
	
	wfr.cleanclass();
}

void Ele_Conductivity::jjcorr_ks(const int ik, const int nt, const double dt,  Wavefunc& wf,
            double* ct11, double* ct12, double* ct22)
{
	const int ndim = 3;
    const int npw = wf.ngtot;

    const int nbands = wf.nband;
    const double ef = INPUT.fermiE / P_Ry2eV;
	double *pij2 = new double[(nbands-1) * nbands / 2];
	getpij2(wf, pij2);
    for (int it = 0; it < nt; ++it)
    {
        double tmct11 = 0;
        double tmct12 = 0;
        double tmct22 = 0;
        double *enb = wf.eigE;
		int ijb = 0;
        for (int ib = 0; ib < nbands; ++ib)
        {
            double ei = enb[ib] / P_Ry2eV;
            double fi = wf.occ[ib];
            for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
            {
                double ej = enb[jb] / P_Ry2eV;
                double fj = wf.occ[jb];
                double tmct = sin((ej - ei) * (it)*dt) * (fi - fj) * pij2[ijb];
                tmct11 += tmct;
                tmct12 += -tmct * ((ei + ej) / 2 - ef);
                tmct22 += tmct * pow((ei + ej) / 2 - ef, 2);
            }
        }
        ct11[it] += tmct11;
        ct12[it] += tmct12;
        ct22[it] += tmct22;
    }
	delete[] pij2;
}

void Ele_Conductivity::calcondw(const int nt, const double dt, const double fwhmin, const double wcut, 
            const double dw_in, double *ct11, double *ct12, double *ct22)
{
	double factor = FACTOR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / P_Ry2eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / P_Ry2eV;
    ofstream ofscond("je-je.txt");
    ofscond << setw(8) << "#t(a.u.)" << setw(15) << "c11(t)" << setw(15) << "c12(t)" << setw(15) << "c22(t)" << setw(15)
            << "decay" << endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << setw(8) << (it)*dt << setw(15) << -2 * ct11[it] << setw(15) << -2 * ct12[it] << setw(15)
                << -2 * ct22[it] << setw(15) << exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) << endl;
    }
    ofscond.close();
    double *cw11 = new double[nw];
    double *cw12 = new double[nw];
    double *cw22 = new double[nw];
    double *kappa = new double[nw];
    ZEROS(cw11, nw);
    ZEROS(cw12, nw);
    ZEROS(cw22, nw);
    for (int iw = 0; iw < nw; ++iw)
    {
        for (int it = 0; it < nt; ++it)
        {
            cw11[iw] += -2 * ct11[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw12[iw] += -2 * ct12[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw22[iw] += -2 * ct22[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond << setw(8) << "## w(eV) " << setw(20) << "sigma(Sm^-1)" << setw(20) << "kappa(W(mK)^-1)" << setw(20)
            << "L12/e(Am^-1)" << setw(20) << "L22/e^2(Wm^-1)" << endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= double(2) / ndim / INPUT.vol * factor; // unit in Sm^-1
        cw12[iw]
            *= double(2) / ndim / INPUT.vol * factor * 2.17987092759e-18 / 1.6021766208e-19; // unit in Am^-1
        cw22[iw] *= double(2) / ndim / INPUT.vol * factor
                    * pow(2.17987092759e-18 / 1.6021766208e-19, 2); // unit in Wm^-1
        kappa[iw] = (cw22[iw] - pow(cw12[iw], 2) / cw11[iw]) / INPUT.temperature;
        ofscond << setw(8) << (iw + 0.5) * dw * P_Ry2eV << setw(20) << cw11[iw] << setw(20) << kappa[iw]
                << setw(20) << cw12[iw] << setw(20) << cw22[iw] << endl;
    }
    cout << setprecision(6) << "DC electrical conductivity: " << cw11[0] - (cw11[1] - cw11[0]) * 0.5 << " Sm^-1"
         << endl;
    cout << setprecision(6) << "Thermal conductivity: " << kappa[0] - (kappa[1] - kappa[0]) * 0.5 << " W(mK)^-1" << endl;
    ;
    ofscond.close();

    delete[] cw11;
    delete[] cw12;
    delete[] cw22;
    delete[] kappa;
}