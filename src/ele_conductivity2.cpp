#include "input.h"
#include "const.h"
#include "gfun.h"
#include "ele_conductivity.h"
#include "matrixmultip.h"


/**
 * @brief Calculate conducitities use time correlation funtions
 * 
 * @note In this method n_fwhm, nscf are not used
 *       smearinvw is not used and it equivalent to smearinvw = 0
 *        
 */
void Ele_Conductivity::method2()
{
	double wcut = INPUT.wcut;
	double dw_in = INPUT.dw;
	double fwhm_in = INPUT.fwhm[0];
	
	int nw = ceil(wcut / dw_in);
    double dw = dw_in / P_Ry2eV; // converge unit in eV to Ry
    double dt = INPUT.cond_dt;
    double sigma;
    int nt;
    if(INPUT.smear == 1)
    {
        sigma = fwhm_in / TWOSQRT2LN2 / P_Ry2eV;
        nt = ceil(sqrt(32) / sigma / dt);
    }
    else if(INPUT.smear )
    {
        sigma = fwhm_in / 2 / P_Ry2eV;
        nt = ceil(16 / sigma / dt);
    }
    else
    {
        std::cout<<"Wrong smear"<<endl;
        exit(0);
    }

    cout << "nw: " << nw << " ; dw: " << dw_in << " eV" << endl;
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
		std::cout<<"kpoint "<<ik+1<<endl;
		std::cout<<"nband: "<<nband<<endl;
		jjcorr_ks(ik, nt, dt, wfr, ct11,ct12,ct22);
		wfr.clean();
	}

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

	if (RANK == 0)
    {
        calcondw(nt, dt, sigma, wcut, dw_in, ct11, ct12, ct22);
    }
	delete[] ct11;
    delete[] ct12;
    delete[] ct22;
	
	wfr.cleanclass();
}

void Ele_Conductivity::jjcorr_ks(const int ik, const int nt, const double dt,  WfRead& wfr,
            double* ct11, double* ct12, double* ct22)
{
	const int ndim = 3;
    const int npw = wfr.wfpt->ngtot;

    const int nbands = wfr.wfpt->nband;
    const double ef = INPUT.fermiE / P_Ry2eV;

    const int nbb = (nbands-1) * nbands / 2;
	double *vmatrix = new double [nbb];
	if(INPUT.readvmatrix)
	{
		wfr.readvmatrix(ik, vmatrix);
	}
	else
	{
        wfr.calvmatrix(vmatrix);
	}
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ct11[:nt], ct12[:nt], ct22[:nt])
#endif	
    for (int it = 0; it < nt; ++it)
    {
        double tmct11 = 0;
        double tmct12 = 0;
        double tmct22 = 0;
        double *enb = wfr.wfpt->eigE;
		int ijb = 0;
        for (int ib = 0; ib < nbands; ++ib)
        {
            double ei = enb[ib] / P_Ry2eV;
            double fi = wfr.wfpt->occ[ib];
            for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
            {
                double ej = enb[jb] / P_Ry2eV;
                double fj = wfr.wfpt->occ[jb];
                double tmct = sin((ej - ei) * (it)*dt) * (fi - fj) * vmatrix[ijb];
                tmct11 += tmct;
                tmct12 += -tmct * ((ei + ej) / 2 - ef);
                tmct22 += tmct * pow((ei + ej) / 2 - ef, 2);
            }
        }
        ct11[it] += tmct11;
        ct12[it] += tmct12;
        ct22[it] += tmct22;
    }
    delete[] vmatrix;
}

void Ele_Conductivity::calcondw(const int nt, const double dt, const double sigma, const double wcut, 
            const double dw_in, double *ct11, double *ct12, double *ct22)
{
	double factor = 4*P_QE*P_QE/P_BOHR/1e-10/P_HBAR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / P_Ry2eV; // converge unit in eV to Ry
    ofstream ofscond("je-je.txt");
    ofscond << setw(8) << "#t(a.u.)" << setw(15) << "c11(t)" << setw(15) << "c12(t)" << setw(15) << "c22(t)" << setw(15)
            << "decay" << endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << setw(8) << (it)*dt << setw(15) << -2 * ct11[it] << setw(15) << -2 * ct12[it] << setw(15)
                << -2 * ct22[it];
        if(INPUT.smear == 1)
                ofscond << setw(15) << exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) << endl;
        else
                ofscond << setw(15) << exp(- sigma * it * dt) << endl;
    }
    ofscond.close();
    double *cw11 = new double[nw];
    double *cw12 = new double[nw];
    double *cw22 = new double[nw];
    double *kappa = new double[nw];
    ZEROS(cw11, nw);
    ZEROS(cw12, nw);
    ZEROS(cw22, nw);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:cw11[:nw], cw12[:nw], cw22[:nw])
#endif	
    for (int iw = 0; iw < nw; ++iw)
    {
        for (int it = 1; it < nt; ++it)
        {
            double tmp;
            if(INPUT.smear == 1)
            {
                tmp = -2 * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            }
            else
            {
                tmp = -2 * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-sigma * it * dt) / (iw + 0.5) / dw * dt;
            }
            cw11[iw] += ct11[it] * tmp;
            cw12[iw] += ct12[it] * tmp;
            cw22[iw] += ct22[it] * tmp;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond << setw(8) << "## w(eV) " << setw(20) << "sigma(Sm^-1)" << setw(20) << "kappa(W(mK)^-1)" << setw(20)
            << "L12/e(Am^-1)" << setw(20) << "L22/e^2(Wm^-1)" << endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= double(2) / ndim / INPUT.vol * factor; // unit in Sm^-1
        cw12[iw]
            *= double(2) / ndim / INPUT.vol * factor * P_Ry2eV; // unit in Am^-1
        cw22[iw] *= double(2) / ndim / INPUT.vol * factor
                    * pow(P_Ry2eV, 2); // unit in Wm^-1
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