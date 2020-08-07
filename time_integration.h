 

#ifndef TIME_INT_H
#define TIME_INT_H

/* solve a volterra integral of first kind to integrate the time-dependence (see derivation_time_solver.pdf)
 * int Nt: number of time-steps per decimation step 
 * int t: current time-step
 * int size: size of the inner matrices (since this program solves both Eq. 165 and 180 in Simon Langs's diploma thesis this is left open)
 * double dt: time step
 * dKInv: Inverse Matrix (see Eq. 10 in derivation_time_solver.pdf)
 * K: input memory kernels (saved for times t>Nt/2)
 * V: output variable (Eq. 180: the memory kernel K, Eq. 165: phi)
 */
//void solve_volterra_integral(int Nt, int t, int q, int size, double dt, double * dKInv, double * K, double * dK, double * V , double * dV);
//void solve_volterra_integral2(int Nt, int t, int q, int size, double dt, double * dKInv, double * K, double * dK, double * V , double * dV);

/* solve an integro-differential equation to integrate the time-dependence (see derivation_time_solver_brownian2.pdf)
 * int Nt: number of time-steps per decimation step 
 * int t: current time-step
 * int size: size of the inner matrices (since this program solves both Eq. 165 and 180 in Simon Langs's diploma thesis this is left open)
 * double dt: time step
 * dKInv: Inverse Matrix (see Eq. 10 in derivation_time_solver.pdf)
 * K: input memory kernels (saved for times t>Nt/2)
 * V: output variable (Eq. 180: the memory kernel K, Eq. 165: phi)
 */
//void solve_volterra_integral_brownian(int Nt, int t, int q, int size, double dt, double * D, double * dKInv, double * K, double * dK, double * V , double * dV);
//void solve_volterra_integral_brownian2(int Nt, int t, int q, int size, double dt, double * D, double * dKInv, double * K, double * dK, double * V , double * dV);



void solve_volterra_integral_M_to_K(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n,double * mathcalDInv,double * dKInv,double * tmp_save, double * K, double * dK, double * V , double * dV) ;
void solve_volterra_integral_M_to_Kr(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n,double * mathcalDInv,double * dKInv,double * tmp_save, double * K, double * dK, double * V , double * dV);
//void solve_volterra_integral_M_to_Kr_predictive(int M,int Nt, int t, int q,double dt, double q0, double dq, double * Q, double * n,double * mathcalDInv,double * dKInv, double * K, double * dK, double * V , double * dV);


void solve_volterra_integral_K_to_Meff(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q,double * n, double * dKInv, double * tmp_save,double * K, double * dK, double * V , double * dV) ;
void solve_volterra_integral_Kr_to_Meff(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q,double * n, double * dKInv, double * tmp_save1, double * tmp_save2,double * K, double * dK, double * V , double * dV) ;

void solve_volterra_integral_Meff_to_Phi(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q,double * n, double * dKInv, double * tmp_save,double * K, double * dK, double * V , double * dV) ;

//void solve_volterra_integral_K_to_Phi(int Nt, int t, int q, int size, double dt, double ** phi, double * dKInv, double * K, double * dK, double * V , double * dV);

#endif
