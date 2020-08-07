/* program to evaluate MCT memory kernels for the slit geometry
 * author: Gerhard Jung (UIBK)
 * date: 18.03.2019 
 */

#include "calc_kernel.h"
#include "mct_slit.h"
#include "nrutil.h"
#include "def.h"

#include <math.h>
#include <stdio.h>
#include <time.h> 

void calc_mkernel_init(int t, int N, int M, double n0, double L,double * Q, double ** c, double ** phi, double ** YS, double ** YM, double ** YSQ, double ** YMQ) {
  
  int Mtot = 2*M+1;
  
  // calculate YS
  for (int q=0; q<N; q++) {
    for (int sigma=0; sigma<Mtot; sigma++) {
      for (int sigma1=0; sigma1<Mtot; sigma1++) {
	for (int sigma2=0; sigma2<Mtot; sigma2++) {
	  YS[q][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] = 0.0;
	  YSQ[q][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] = 0.0;
	  for (int kappa=0; kappa<Mtot; kappa++) {
	    if (sigma-sigma2+M>=0 && sigma-sigma2+M<Mtot) {
	      YS[q][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] += c[q][kappa*Mtot+sigma-sigma2+M]*phi[q][t*Mtot*Mtot+kappa*Mtot+sigma1];
	      YSQ[q][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] += c[q][kappa*Mtot+sigma-sigma2+M]*Q[sigma-sigma2+3*M]*phi[q][t*Mtot*Mtot+kappa*Mtot+sigma1];
	    }
	  }
	  //printf("%f %f\n",YS[q][sigma2+sigma1*Mtot+sigma2*Mtot*Mtot],YSQ[q][sigma2+sigma1*Mtot+sigma2*Mtot*Mtot]);
	}
      }
    }
  }
  // calculate YM
  for (int q1=0; q1<N; q1++) {
    for (int q2=0; q2<N; q2++) {
      for (int sigma=0; sigma<Mtot; sigma++) {
	for (int sigma1=0; sigma1<Mtot; sigma1++) {
	  for (int sigma2=0; sigma2<Mtot; sigma2++) {
	    YM[q1*N+q2][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] = 0.0;
	    YMQ[q1*N+q2][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] = 0.0;
	    for (int kappa=0; kappa<Mtot; kappa++) {
	      if (sigma-kappa+M>=0 && sigma-kappa+M<Mtot) {
		YM[q1*N+q2][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] += c[q1][sigma1*Mtot+sigma-kappa+M]*phi[q2][t*Mtot*Mtot+kappa*Mtot+sigma2];
		YMQ[q1*N+q2][sigma2+sigma1*Mtot+sigma*Mtot*Mtot] += c[q1][sigma1*Mtot+sigma-kappa+M]*Q[sigma-kappa+3*M]*phi[q2][t*Mtot*Mtot+kappa*Mtot+sigma2];
	      }
	    }
	    //printf("%f %f\n",YM[q1*N+q2][sigma2+sigma1*Mtot+sigma2*Mtot*Mtot],YMQ[q1*N+q2][sigma2+sigma1*Mtot+sigma2*Mtot*Mtot]);
	  }
	}
      }
    }
  }
  
  
}

/**********************************************************************************************************************************************************/

void calc_mkernel(int t, int N, int q,  int M, double n0, double L, double q0, double dq,   double ** c, double *Q, double *n, double ** YS, double ** YM, double ** YSQ, double ** YMQ, double ** phi, double * m) {
  //time_t start = time(0);
  //printf("Start kernel calculation ... Optimized! \n");
  
  int Mtot = 2*M+1;
  
  double prefactor = 1.0/(2.0*PI*PI)*n0/(L*L*L*L)*dq*dq;
  //double prefactor = 1.0/(2.0*PI*PI)*n0*n0*n0/(L*L*L*L*L*L*L*L)*dq*dq*chi_scale;
  
  double q_val = q0 + q*dq; 
  
  for (int a=0; a<2; a++) {
    for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
      int b=a;
      int nu=mu;
#else
      for (int b=0; b<2; b++) {
	for (int nu=0; nu<Mtot; nu++) {
#endif

	  //printf("mu %d nu %d\n",mu,nu);
	  int index = nu + b*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot;
	  // initialize m vector
	  m[index] = 0.0;
	  
	  if (a==b && nu<mu) m[index]=m[mu + b*Mtot + nu*2*Mtot + 2*a*Mtot*Mtot];
	  else if (a==1 && b==0) m[index]=m[mu + 1*Mtot + nu*2*Mtot + 2*0*Mtot*Mtot];
	  else if (a==b && mu > M) m[index]=m[2*M-nu + b*Mtot + (2*M-mu)*2*Mtot + 2*a*Mtot*Mtot];
	  else if (a!=b && mu > M) m[index]=-m[2*M-nu + b*Mtot + (2*M-mu)*2*Mtot + 2*a*Mtot*Mtot];
	  else {
	    // loop over all wavevectors q1,q2
	    for (int q1=0; q1<N;q1++) {
	      double q1_val = q0+q1*dq;
	      for (int q2=ABS(q-q1); q2<MIN(q+q1+1,N);q2++) {
		double q2_val = q0+q2*dq;
		
		double local_result=0.0;
		
		double prefactor_q = q1_val*q2_val/sqrt(4.0*q_val*q_val*q1_val*q1_val-(q_val*q_val+q1_val*q1_val-q2_val*q2_val)*(q_val*q_val+q1_val*q1_val-q2_val*q2_val));
		double qq1 = (q_val*q_val + q1_val*q1_val - q2_val*q2_val)/(2.0*q_val);
		double qq2 = (q_val*q_val + q2_val*q2_val - q1_val*q1_val)/(2.0*q_val);
		
		// due to the cutoff in the indices mu (etc...) also the difference mu - mu1 (etc...) has to be cut off (DELETED)
		//int mu_min = MAX(mu-M,0);
		//int mu_max = MIN(mu+M+1,2*M+1);
		//int nu_min = MAX(nu-M,0);
		//int nu_max = MIN(nu+M+1,2*M+1);
  #ifdef DIAGONAL_APPROXIMATION1
		// exactly the same as in the previous execer since J_mu = 1
		for (int mu1=0; mu1<Mtot; mu1++) {
		  if ( (mu-mu1+M >= 0) && (mu-mu1+M<Mtot) ) {
		    double val = 0.0;
		    if (a==0) val = qq1*c[q1][mu1*Mtot+mu1]+qq2*c[q2][(mu-mu1+M)*Mtot+(mu-mu1+M)];
		    else val = Q[mu1+2*M]*c[q1][mu1*Mtot+mu1]+Q[mu-mu1+3*M]*c[q2][(mu-mu1+M)*Mtot+(mu-mu1+M)];
		    //if (q==0 && mu==3) printf("Q1 %f Q2 %f\n",Q[mu1+2*M],Q[mu-mu1+3*M]);
		    local_result += val*val*phi[q1][t*Mtot*Mtot+mu1*Mtot+mu1]*phi[q2][t*Mtot*Mtot+(mu-mu1+M)*Mtot+(mu-mu1+M)];
		  }
		}
		//if (q==0 && mu==3) printf("local_result %f\n",local_result);
  #else
		for (int gamma=0; gamma<4; gamma++) {
		  double local_result_gamma = 0.0;
		  double prefactor_gamma = 1.0;
		  if (a==0 && b==0 ) {
		    if (gamma==0) prefactor_gamma = qq1*qq1;
		    else if (gamma==1 || gamma==2) prefactor_gamma = qq1*qq2;
		    else prefactor_gamma = qq2*qq2;
		  } else if (a==0 && b==1) {
		    if (gamma==0 || gamma==2) prefactor_gamma = qq1;
		    else prefactor_gamma = qq2;
		  } else if (a==1 && b==0) {
		    if (gamma==0 || gamma==1) prefactor_gamma = qq1;
		    else prefactor_gamma = qq2;
		  }
		  
		  //printf("%f \n",prefactor_gamma);
		  
		  //set pointers outside to reduce if-clause in the lowest loops
		  double ** pa = NULL;
		  double ** pb = NULL;
		  if (a==0) {
		    pa = YS;
		  } else {
		    pa = YSQ;
		  }
		  if (b==0) {
		    if (gamma==0||gamma==3) pb = YM;
		    else pb = YS;
		  } else {
		    if (gamma==0||gamma==3) pb = YMQ;
		    else pb = YSQ;
		  }
		  
		  // some constants
		  int index_mu = mu*Mtot*Mtot;
		  int index_nu = nu*Mtot*Mtot;
		  int indexq1q2 = q2+q1*N;
		  int indexq2q1 = q1+q2*N;
		  
		  // the most important loops
		  if (gamma==0) {
		    for (int sigma1=0; sigma1<Mtot; sigma1++) {
		      for (int sigma2=0; sigma2<Mtot; sigma2++) {
			local_result_gamma += pa[q1][sigma2+sigma1*Mtot+index_mu]*pb[indexq1q2][sigma2+sigma1*Mtot+index_nu];  
			}
		    }
		  } else if (gamma==1) {
		    for (int sigma1=0; sigma1<Mtot; sigma1++) {
		      for (int sigma2=0; sigma2<Mtot; sigma2++) {
			local_result_gamma += pa[q1][sigma2+sigma1*Mtot+index_mu]*pb[q2][sigma1+sigma2*Mtot+index_nu];
			}
		    }
		  } else if (gamma==2) {
		    for (int sigma1=0; sigma1<Mtot; sigma1++) {
		      for (int sigma2=0; sigma2<Mtot; sigma2++) {
			local_result_gamma += pa[q2][sigma2+sigma1*Mtot+index_mu]*pb[q1][sigma1+sigma2*Mtot+index_nu];   
			}
		    }
		  } else {
		    for (int sigma1=0; sigma1<Mtot; sigma1++) {
		      for (int sigma2=0; sigma2<Mtot; sigma2++) {
			local_result_gamma += pa[q2][sigma2+sigma1*Mtot+index_mu]*pb[indexq2q1][sigma2+sigma1*Mtot+index_nu];
			}
		    }
		  }
		  
		  local_result += local_result_gamma*prefactor_gamma;
		}
  #endif
		
		m[index] += local_result*prefactor_q;
		
	      } //q2-loop
	    } // q1-loop
	    m[index]*= prefactor;
	  }

	  //if (a==b && nu<mu) printf("calc. %f sym1. %f\n",m[index],m[mu + b*Mtot + nu*2*Mtot + 2*a*Mtot*Mtot]);
	  //if (a==1 && b==0) printf("calc. %f sym2. %f\n",m[index],m[mu + 1*Mtot + nu*2*Mtot + 2*0*Mtot*Mtot]);
	  //if (a==b && mu > M) printf("calc. %.15f sym3. %.15f\n",m[index],m[2*M-nu + b*Mtot + (2*M-mu)*2*Mtot + 2*a*Mtot*Mtot]);
	  //if (a!=b && mu > M) printf("calc. %.15f sym4. %.15f\n",m[index],-m[2*M-nu + b*Mtot + (2*M-mu)*2*Mtot + 2*a*Mtot*Mtot]);
	    //printf("mu %d nu %d: %f\n",mu,nu,m[index][q]);
	  
	  //printf("%d mu %d nu %d: %f\n",t,mu,nu,m[index]);
#ifndef DIAGONAL_APPROXIMATION1
	}
      }
#endif
    }
  }
 
  
#ifdef DIAGONAL_APPROXIMATION1

#else

// calc DInv
  double * tmp_D;
  double * DInv;
  tmp_D = dvector(0,4*Mtot*Mtot-1);
  DInv = dvector(0,4*Mtot*Mtot-1);
  for (int a=0; a<2; a++) {
    for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
      int nu=mu;
      int b=a;
#else
      for (int b=0; b<2; b++) {
        for (int nu=0; nu<Mtot; nu++) {
#endif
          if (a==b) tmp_D[nu+b*Mtot+2*Mtot*mu+2*Mtot*Mtot*a] =  n[mu -nu + 3*M]/n[3*M];
          else tmp_D[nu+b*Mtot+2*Mtot*mu+2*Mtot*Mtot*a] =  0.0;
#ifndef DIAGONAL_APPROXIMATION1
        }
      }
#endif
    }
  }
  calc_matrixInv(tmp_D,2*Mtot,DInv);

  // save m
  double * mloc = dvector(0, 4*Mtot*Mtot-1);
  
  for (int a=0; a<2; a++) {
    for (int mu=0; mu<Mtot; mu++) {
      for (int b=0; b<2; b++) {
        for (int nu=0; nu<Mtot; nu++) {
          int index = nu + b*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot;
          mloc[index] = m[index];
        }
      }
    }
  }
  
  // include DInv
  for (int mu=0; mu<2*Mtot; mu++) {
    for (int nu=0; nu<2*Mtot; nu++) {
      m[nu+mu*2*Mtot] = 0.0;
      for (int sigma=0; sigma<2*Mtot; sigma++) {
 /*       for (int kappa=0; kappa<2*Mtot; kappa++) {
          m[nu+mu*2*Mtot] += DInv[sigma + mu*2*Mtot]*mloc[kappa+sigma*2*Mtot]*DInv[nu + kappa*2*Mtot];
        }
  */
	m[nu+mu*2*Mtot] += mloc[sigma+mu*2*Mtot]*DInv[nu + sigma*2*Mtot];
      }
    }
  }
	  
  free_dvector(mloc,0, 4*Mtot*Mtot-1);
  free_dvector(tmp_D,0, 4*Mtot*Mtot-1);
  free_dvector(DInv,0, 4*Mtot*Mtot-1);
#endif
  
  //time_t end = time(0);
  //double time = difftime(end, start);
  //printf("Finished kernel calculation ... Optimized! Time: %f \n",time);
}

