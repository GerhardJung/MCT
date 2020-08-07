
#ifndef MCT_CONFINED
#define MCT_CONFINED

void read_c_phi_modes(int Nt, int N, int M, double q0, double dq, double ** c, double ** phi, double ** dphi);
void read_n_modes(int M, double L, double * n, double * v,const char* density_file);

// contraction
void contract_K(int q, int  M,double q0,double dq,double L, double n0, double * Q,double * v,double ** c, double * Kc, double * K);

void init_values_moments(int Nt, int N, int M, double n0, double  L, double q0, double dq, double dt,double D0, double * Q,double * v, double * n, double ** c,   double * mathcalD, double ** DInv, double ** YS, double ** YM, double ** YSQ, double ** YMQ, double ** phi, double **m ,double **K ,double **Kr ,double **Kc,double **Krc , double **Meff , double **Meff_long , double **dm ,double **dK , double **dKr , double **dKc  , double **dKrc , double **dMeff , double **dMeff_long);

// decimation step to go to better time resolution
void decimation_moments(int Nt, int N, int size,double ** dinput,double ** input);
void decimation_values(int Nt, int N, int size,double ** input);

void calc_dKInv(int N, int size, double dt, double *mathcalD, double **Dif, double **DInv, double **DphiInv, double ** phi, double ** dK , double ** dKInv, int integro);

void initInv(int N, int M, double n0, double q0, double dq, double D0, double *Q, double * n, double ** phi,  double **D, double **DInv, double **DphiInv, double *mathcalD, double *mathcalDInv);

void calc_matrixInv(double *in, int size, double * out);

void enlarge_array_values(int N, int size, int Nt, double *** array, double ** array_value);
void enlarge_array_moments(int N, int size, int Nt, double *** array, double ** array_moment, double ** array_value);

#endif
