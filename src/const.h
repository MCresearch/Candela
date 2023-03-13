#ifndef CONST_H
#define CONST_H

constexpr double P_BOHR = 0.52917720859; //in QE 
constexpr double P_HBAR = 1.05457266e-34;
constexpr double P_ME =  9.10956e-31;
constexpr double P_QE = 1.6021766341e-19;
constexpr double P_KB = 1.380649e-23;

constexpr double P_HA = P_HBAR*P_HBAR/(P_BOHR * P_BOHR * 1e-20)/P_ME/P_QE;
constexpr double P_EV2K = P_QE/P_KB;
constexpr double P_Ry2eV = P_HA / 2;
constexpr double P_t_AU = P_HBAR/P_Ry2eV/P_QE;
 

#define MASS_H  1.00794
#define MASS_D  2.01410
#define MASS_T  3.016
#define MASS_Li  6.941
#define MASS_Be  9.0121831
#define MASS_B  10.811
#define MASS_C  12.011
#define MASS_O  15.9994
#define MASS_F  18.9984
#define MASS_Na  22.98976928
#define MASS_Mg  24.3050
#define MASS_Al  26.981539
#define MASS_Si  28.0855

#endif
