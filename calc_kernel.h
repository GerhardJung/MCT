 
#ifndef CALC_KERNEL
#define CALC_KERNEL

#define PI 3.14159265359
#define ABS(value)  ( (value) >=0 ? (value) : -(value) )
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
     
     
/* function: Initilize some precalculation for mkernel calculation
 */
void calc_mkernel_init(int t, int N, int M, double n0, double L,  double * Q, double ** c, double ** phi, double ** YS, double ** YM, double ** YSQ, double ** YMQ);

/* function: calculate memory kernel for slit geometry (Eq. 199 in Simon Lang's diploma thesis). 
 * This is the optimized function which precalculates some sums to improve the scaling to O(N^3*M^4) (instead of O(N^3*M^6)).
 * It also uses all the symmetry to calculate as few as possible.
 * 
 * input: 
 * number of time steps Nt
 * current time step to be calculated t
 * number of wavevector modes N
 * current wavevector to be calculated q
 * total average density n0
 * (accessible) slit size L=D-1
 * number of index modes M
 * offset wavevector q0
 * stepsize wavevector dq
 * YS/YM/YSQ/YMQ precalculated sums for mkernel calculation (see notes)
 * output memory kernel m
 */
void calc_mkernel(int t, int N, int q,  int M, double n0, double L, double q0, double dq,  double ** c, double *Q, double *n, double ** YS, double ** YM, double ** YSQ, double ** YMQ, double ** phi, double * m);


/* function: Initilize some precalculation for mkernel calculation
 */
void calc_gkernel_init(int t, int N, int M, double n0, double L,  double * Q, double * v, double ** c, double ** phi, double ** YS, double ** YM, double ** YSQ, double ** YMQ);
/* function: calculate memory kernel for slit geometry (Eq. 205 in Simon Lang's diploma thesis). 
 * This is the optimized function which precalculates some sums to improve the scaling to O(N^3*M^4) (instead of O(N^3*M^6)).
 * It also uses all the symmetry to calculate as few as possible.
 * 
 * input: 
 * number of time steps Nt
 * current time step to be calculated t
 * number of wavevector modes N
 * current wavevector to be calculated q
 * total average density n0
 * (accessible) slit size L=D-1
 * number of index modes M
 * offset wavevector q0
 * stepsize wavevector dq
 * YS/YM/YSQ/YMQ precalculated sums for mkernel calculation (see notes)
 * output memory kernel m
 */
void calc_gkernel(int t, int N, int q,  int M, double n0, double L, double q0, double dq, double ** c, double *Q,double ** YS, double ** YM, double ** YSQ, double ** YMQ, double ** phi, double * m);



#endif
