
#include "time_integration.h"
#include "mct_slit.h"
#include "nrutil.h"
#include "def.h"

#include <stdio.h>

void solve_volterra_integral_M_to_K(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n,double * mathcalDInv,double * dKInv,double * tmp_save, double * K, double * dK, double * V , double * dV){
  int Mtot = 2*M+1;
  double * tmp;
  tmp = dvector(0,4*Mtot*Mtot-1);
  int ibar = t/2;
  
  // calculate D
  double * tmp_D;
  tmp_D = dvector(0,4*Mtot*Mtot-1);
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
  
  // calculate integration step
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    int nu=mu;
#else
    for (int nu=0; nu<2*Mtot; nu++) {
#endif
      
      if (step==0) { 
	tmp[nu+mu*2*Mtot] = 0.0;
	// include convolution
	for (int j=2; j<=ibar; j++) {
  #ifdef DIAGONAL_APPROXIMATION1
	  int sigma = nu;
	  int kappa = mu;
	  tmp[nu+mu*2*Mtot] -= dt/2.0*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]+K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
  #else
	  for (int sigma=0; sigma<2*Mtot; sigma++) {
	    tmp[nu+mu*2*Mtot] -= dt/2.0*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]+K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
	  }
  #endif
	}

	for (int j=2; j<=t-ibar; j++) {
  #ifdef DIAGONAL_APPROXIMATION1
	  int sigma = nu;
	  int kappa = mu;
	  tmp[nu+mu*2*Mtot] -= dt/2.0*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]+V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
  #else
	  for (int sigma=0; sigma<2*Mtot; sigma++) {
	    tmp[nu+mu*2*Mtot] -= dt/2.0*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]+V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
	  }
  #endif
	}
	tmp_save[nu+mu*2*Mtot] = tmp[nu+mu*2*Mtot];
      } else {
	tmp[nu+mu*2*Mtot] = tmp_save[nu+mu*2*Mtot];
      }
      
      
#ifdef DIAGONAL_APPROXIMATION1
      int sigma = mu;
      tmp[nu+mu*2*Mtot] -= K[sigma+mu*2*Mtot + 4*Mtot*Mtot*t]*tmp_D[nu+sigma*2*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
	tmp[nu+mu*2*Mtot] -= K[sigma+mu*2*Mtot + 4*Mtot*Mtot*t]*tmp_D[nu+sigma*2*Mtot];
      }
#endif
      
#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*2*Mtot] -= dt/2.0*dK[sigma+mu*2*Mtot+4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
	tmp[nu+mu*2*Mtot] -= dt/2.0*dK[sigma+mu*2*Mtot+4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
      }
#endif
      
#ifdef DIAGONAL_APPROXIMATION1
      int sigma = nu;
      int kappa = mu;
      tmp[nu+mu*2*Mtot] -= dt/2.0*(K[sigma+mu*2*Mtot+(t)*4*Mtot*Mtot]+K[sigma+mu*2*Mtot+(t-1)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+4*Mtot*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
	tmp[nu+mu*2*Mtot] -= dt/2.0*(K[sigma+mu*2*Mtot+(t)*4*Mtot*Mtot]+K[sigma+mu*2*Mtot+(t-1)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+4*Mtot*Mtot];
      }
#endif
      
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    V[mu+mu*2*Mtot + 4*Mtot*Mtot*t] = dKInv[mu+mu*2*Mtot]*tmp[mu+mu*2*Mtot];
#else
    for (int nu=0; nu<2*Mtot; nu++) {
      V[nu+mu*2*Mtot + 4*Mtot*Mtot*t] = 0.0;
      for (int sigma=0; sigma<2*Mtot; sigma++) {
        V[nu+mu*2*Mtot + 4*Mtot*Mtot*t] += dKInv[sigma+mu*2*Mtot]*tmp[nu+sigma*2*Mtot];
      }
    }
#endif
  }
  
  free_dvector(tmp,0,4*Mtot*Mtot-1);
  free_dvector(tmp_D,0,4*Mtot*Mtot-1);
}

/**********************************************************************************************************************************************************************************************************/

void solve_volterra_integral_M_to_Kr(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n,double * mathcalDInv,double * dKInv,double * tmp_save, double * K, double * dK, double * V , double * dV){
  int Mtot = 2*M+1;
  double * tmp;
  tmp = dvector(0,4*Mtot*Mtot-1);
  int ibar = t/2;
  
  // calculate D
  double * tmp_D;
  tmp_D = dvector(0,4*Mtot*Mtot-1);
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
  
  // calculate integration step
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<2*Mtot; nu++) {
#endif
      
      if (step==0) { 
	tmp[nu+mu*2*Mtot] = 0.0;
	// include convolution
	for (int j=2; j<=ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int sigma = mu;
	  tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
#else
	  for (int sigma=0; sigma<2*Mtot; sigma++) {
	    tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
	  }
#endif
	}

	for (int j=2; j<=t-ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int sigma = mu;
	  int kappa = mu;
	  tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
#else
	  for (int sigma=0; sigma<2*Mtot; sigma++) {
	    tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
	  }
#endif
	}
	tmp_save[nu+mu*2*Mtot] = tmp[nu+mu*2*Mtot];
      } else {
	tmp[nu+mu*2*Mtot] = tmp_save[nu+mu*2*Mtot];
      }
      
      // include previous timesteps
      tmp[nu+mu*2*Mtot] += alpha1*V[nu+mu*2*Mtot+(t-1)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] += alpha2*V[nu+mu*2*Mtot+(t-2)*4*Mtot*Mtot];

#ifdef DIAGONAL_APPROXIMATION1
      int sigma = mu;
      int kappa = mu;
      tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t-ibar)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar)*4*Mtot*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
          tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t-ibar)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar)*4*Mtot*Mtot];
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
        int sigma = mu;
        tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-1)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+4*Mtot*Mtot];
#else
        for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-1)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+4*Mtot*Mtot];
        }
#endif
    
      tmp[nu+mu*2*Mtot] += beta0*dt*tmp_D[nu+mu*2*Mtot];
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    V[mu+mu*2*Mtot + 4*Mtot*Mtot*t] = dKInv[mu+mu*2*Mtot]*tmp[mu+mu*2*Mtot];
#else
    for (int nu=0; nu<2*Mtot; nu++) {
      V[nu+mu*2*Mtot + 4*Mtot*Mtot*t] = 0.0;
      for (int sigma=0; sigma<2*Mtot; sigma++) {
        V[nu+mu*2*Mtot + 4*Mtot*Mtot*t] += dKInv[sigma+mu*2*Mtot]*tmp[nu+sigma*2*Mtot];
      }
    }
#endif
  }
  
  free_dvector(tmp,0,4*Mtot*Mtot-1);
  free_dvector(tmp_D,0,4*Mtot*Mtot-1);
}

/**********************************************************************************************************************************************************************************************************/

/*void solve_volterra_integral_M_to_Kr_predictive(int M,int Nt, int t, int q,double dt, double q0, double dq, double * Q,double * n,double * mathcalDInv, double * dKInv, double * K, double * dK, double * V , double * dV){
  int Mtot = 2*M+1;
  double * tmp, * val_t1, * val_t2;
  tmp = dvector(0,4*Mtot*Mtot-1);
  val_t1 = dvector(0,4*Mtot*Mtot-1);
  val_t2 = dvector(0,4*Mtot*Mtot-1);
  int ibar = t/2;
  int ibar1= (t+1)/2;
  double qval = q0 + q*dq;
  
  // calculate D
  double * tmp_D;
  tmp_D = dvector(0,4*Mtot*Mtot-1);
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
  
  // calculate integration step 1
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<2*Mtot; nu++) {
#endif
      // include previous timesteps
      tmp[nu+mu*2*Mtot] = alpha1*V[nu+mu*2*Mtot+(t-1)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] += alpha2*V[nu+mu*2*Mtot+(t-2)*4*Mtot*Mtot];

#ifdef DIAGONAL_APPROXIMATION1
      int sigma = mu;
      int kappa = mu;
      tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t-ibar)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar)*4*Mtot*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t-1)*4*Mtot*Mtot];
          tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t-ibar)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar)*4*Mtot*Mtot];
      }
#endif
      
      // include convolution
      for (int j=1; j<=ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int sigma = mu;
        tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
#else
        for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
        }
#endif
      }

      for (int j=2; j<=t-ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int sigma = mu;
        int kappa = mu;
        tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
#else
        for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t-j)*4*Mtot*Mtot]);
        }
#endif
      }
    
      tmp[nu+mu*2*Mtot] += beta0*dt*tmp_D[nu+mu*2*Mtot];
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    val_t1[mu+mu*2*Mtot ] = dKInv[mu+mu*2*Mtot]*tmp[mu+mu*2*Mtot];
#else
    for (int nu=0; nu<2*Mtot; nu++) {
      val_t1[nu+mu*2*Mtot ] = 0.0;
      for (int sigma=0; sigma<2*Mtot; sigma++) {
        val_t1[nu+mu*2*Mtot] += dKInv[sigma+mu*2*Mtot]*tmp[nu+sigma*2*Mtot];
      }
    }
#endif
  }
  
  // interpolate for predictive move
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<2*Mtot; nu++) {
#endif
      K[nu+mu*2*Mtot+(t+1)*4*Mtot*Mtot] = K[nu+mu*2*Mtot+(t)*4*Mtot*Mtot] + (K[nu+mu*2*Mtot+(t)*4*Mtot*Mtot] - K[nu+mu*2*Mtot+(t-1)*4*Mtot*Mtot])*dt;
      V[nu+mu*2*Mtot + t*4*Mtot*Mtot] = val_t1[nu+mu*2*Mtot];
      #ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  
  // calculate integration step 2
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<2*Mtot; nu++) {
#endif
      // include previous timesteps
      tmp[nu+mu*2*Mtot] = alpha1*V[nu+mu*2*Mtot+(t)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] += alpha2*V[nu+mu*2*Mtot+(t-1)*4*Mtot*Mtot];

#ifdef DIAGONAL_APPROXIMATION1
      int sigma = mu;
      int kappa = mu;
      tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t)*4*Mtot*Mtot];
      tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t+1-ibar1)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar1)*4*Mtot*Mtot];
#else
      for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] += beta0*dt*dK[sigma+mu*2*Mtot + 4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(t)*4*Mtot*Mtot];
          tmp[nu+mu*2*Mtot] -= beta0*dt*K[sigma+mu*2*Mtot + (t+1-ibar1)*4*Mtot*Mtot]*V[nu+sigma*2*Mtot+(ibar1)*4*Mtot*Mtot];
      }
#endif
      
      // include convolution
      for (int j=1; j<=ibar1; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int sigma = mu;
        tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t+1-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t+1-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
#else
        for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] -= beta0*dt*(K[sigma+mu*2*Mtot+(t+1-j+1)*4*Mtot*Mtot]-K[sigma+mu*2*Mtot+(t+1-j)*4*Mtot*Mtot])*dV[nu+sigma*2*Mtot+(j)*4*Mtot*Mtot];
        }
#endif
      }

      for (int j=2; j<=t+1-ibar1; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int sigma = mu;
        int kappa = mu;
        tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t+1-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t+1-j)*4*Mtot*Mtot]);
#else
        for (int sigma=0; sigma<2*Mtot; sigma++) {
          tmp[nu+mu*2*Mtot] -= beta0*dt*dK[sigma+mu*2*Mtot+(j)*4*Mtot*Mtot]*(V[nu+sigma*2*Mtot+(t+1-j+1)*4*Mtot*Mtot]-V[nu+sigma*2*Mtot+(t+1-j)*4*Mtot*Mtot]);
        }
#endif
      }
    
      tmp[nu+mu*2*Mtot] += beta0*dt*tmp_D[nu+mu*2*Mtot];
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    val_t2[mu+mu*2*Mtot ] = dKInv[mu+mu*2*Mtot]*tmp[mu+mu*2*Mtot];
#else
    for (int nu=0; nu<2*Mtot; nu++) {
      val_t2[nu+mu*2*Mtot ] = 0.0;
      for (int sigma=0; sigma<2*Mtot; sigma++) {
        val_t2[nu+mu*2*Mtot ] += dKInv[sigma+mu*2*Mtot]*tmp[nu+sigma*2*Mtot];
      }
    }
#endif
  }
  
    // interpolate
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      V[nu+mu*2*Mtot + t*4*Mtot*Mtot] = (V[nu+mu*2*Mtot + (t-1)*4*Mtot*Mtot] + 3.0*val_t1[nu+mu*2*Mtot] + val_t2[nu+mu*2*Mtot])/5.0;
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  
  free_dvector(tmp,0,4*Mtot*Mtot-1);
  free_dvector(val_t1,0,4*Mtot*Mtot-1);
  free_dvector(val_t2,0,4*Mtot*Mtot-1);
  free_dvector(tmp_D,0,4*Mtot*Mtot-1);
}*/

/**********************************************************************************************************************************************************/

void solve_volterra_integral_K_to_Meff(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n, double * dKInv, double * tmp_save, double * K, double * dK, double * V , double * dV) {
  int Mtot = 2*M+1;
  double * tmp;
  tmp = dvector(0,Mtot*Mtot-1);
  int ibar = t/2;
  double qval = q0 + q*dq;
  
  // calculate inverse of D
  double * tmp_D;
  double * DInv;
  tmp_D = dvector(0,Mtot*Mtot-1);
  DInv = dvector(0,Mtot*Mtot-1);
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      tmp_D[mu*Mtot+nu] =  n[mu -nu + 3*M]/n[3*M] * (qval * qval + Q[mu+2*M]*Q[nu+2*M] );
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  calc_matrixInv(tmp_D,Mtot,DInv);
  
  // calculate integration step
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      
      if (step==0) {
	tmp[nu+mu*Mtot] = 0.0;
	// include convolution
	for (int j=2; j<=ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int sigma = nu;
	  int kappa=mu;
	  tmp[nu+mu*Mtot] -= dt/2.0*dV[sigma+mu*Mtot+(j)*Mtot*Mtot]*(K[kappa+sigma*Mtot+(t-j+1)*Mtot*Mtot]+K[kappa+sigma*Mtot+(t-j)*Mtot*Mtot])*DInv[nu+kappa*Mtot];
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    for (int kappa=0; kappa<Mtot; kappa++) {
	      tmp[nu+mu*Mtot] -= dt/2.0*dV[sigma+mu*Mtot+(j)*Mtot*Mtot]*(K[kappa+sigma*Mtot+(t-j+1)*Mtot*Mtot]+K[kappa+sigma*Mtot+(t-j)*Mtot*Mtot])*DInv[nu+kappa*Mtot];
	    }
	  }
#endif
	}

	for (int j=2; j<=t-ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int sigma = nu;
	  int kappa=mu;
	  tmp[nu+mu*Mtot] -= dt/2.0*(V[sigma+mu*Mtot+(t-j+1)*Mtot*Mtot]+V[sigma+mu*Mtot+(t-j)*Mtot*Mtot])*dK[kappa+sigma*Mtot+(j)*Mtot*Mtot]*DInv[nu+kappa*Mtot];
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    for (int kappa=0; kappa<Mtot; kappa++) {
	      tmp[nu+mu*Mtot] -= dt/2.0*(V[sigma+mu*Mtot+(t-j+1)*Mtot*Mtot]+V[sigma+mu*Mtot+(t-j)*Mtot*Mtot])*dK[kappa+sigma*Mtot+(j)*Mtot*Mtot]*DInv[nu+kappa*Mtot];
	    }
	  }
#endif
	} 
	tmp_save[nu+mu*Mtot] = tmp[nu+mu*Mtot];
      } else {
	tmp[nu+mu*Mtot] = tmp_save[nu+mu*Mtot];
      }
      
#ifdef DIAGONAL_APPROXIMATION1
      int sigma=mu;
      int kappa=mu;
      tmp[nu+mu*Mtot] += -DInv[sigma+mu*Mtot]*K[kappa+sigma*Mtot+t*Mtot*Mtot]*DInv[nu+kappa*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        for (int kappa=0; kappa<Mtot; kappa++) {
          tmp[nu+mu*Mtot] += -DInv[sigma+mu*Mtot]*K[kappa+sigma*Mtot+t*Mtot*Mtot]*DInv[nu+kappa*Mtot];
        }
      }
#endif
         
#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] -= dt/2.0*V[sigma+mu*Mtot+(t-1)*Mtot*Mtot]*dK[kappa+sigma*Mtot+Mtot*Mtot]*DInv[nu+kappa*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        for (int kappa=0; kappa<Mtot; kappa++) {
          tmp[nu+mu*Mtot] -=  dt/2.0*V[sigma+mu*Mtot+(t-1)*Mtot*Mtot]*dK[kappa+sigma*Mtot+Mtot*Mtot]*DInv[nu+kappa*Mtot];
        }
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      int sigma = nu;
      int kappa=mu;
      tmp[nu+mu*Mtot] -= dt/2.0*dV[sigma+mu*Mtot+*Mtot*Mtot]*(K[kappa+sigma*Mtot+(t)*Mtot*Mtot]+K[kappa+sigma*Mtot+(t-1)*Mtot*Mtot])*DInv[nu+kappa*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
	for (int kappa=0; kappa<Mtot; kappa++) {
	  tmp[nu+mu*Mtot] -= dt/2.0*dV[sigma+mu*Mtot+Mtot*Mtot]*(K[kappa+sigma*Mtot+(t)*Mtot*Mtot]+K[kappa+sigma*Mtot+(t-1)*Mtot*Mtot])*DInv[nu+kappa*Mtot];
	}
      }
#endif
      
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    V[mu+mu*Mtot + Mtot*Mtot*t] = dKInv[mu+mu*Mtot]*tmp[mu+mu*Mtot];
#else
    for (int nu=0; nu<Mtot; nu++) {
      V[nu+mu*Mtot + Mtot*Mtot*t] = 0.0;
      for (int sigma=0; sigma<Mtot; sigma++) {
        V[nu+mu*Mtot + Mtot*Mtot*t] += tmp[sigma+mu*Mtot]*dKInv[nu+sigma*Mtot];
      }
    }
#endif
  }
  
  free_dvector(tmp,0,Mtot*Mtot-1);
  free_dvector(tmp_D,0,Mtot*Mtot-1);
  free_dvector(DInv,0,Mtot*Mtot-1);
  
}

/**********************************************************************************************************************************************************/

void solve_volterra_integral_Kr_to_Meff(int M,int Nt, int t, int q,double dt, double q0, double dq, int step, double * Q, double * n,double * dKInv,double * tmp_save1,double * tmp_save2, double * K, double * dK, double * V , double * dV) {
  int Mtot = 2*M+1;
  double * tmp, * val_t1, * val_t2;
  tmp = dvector(0,Mtot*Mtot-1);
  val_t1 = dvector(0,Mtot*Mtot-1);
  val_t2 = dvector(0,Mtot*Mtot-1);
  int ibar = t/2;
  int ibar1= (t+1)/2;
  double qval = q0 + q*dq;
  
  // calculate inverse of D
  double * tmp_D;
  double * DInv;
  tmp_D = dvector(0,Mtot*Mtot-1);
  DInv = dvector(0,Mtot*Mtot-1);
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      tmp_D[mu*Mtot+nu] =  n[mu -nu + 3*M]/n[3*M] * (qval * qval + Q[mu+2*M]*Q[nu+2*M] );
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  calc_matrixInv(tmp_D,Mtot,DInv);
  
  // calculate integration step 1
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      if (step==0) {
	tmp[nu+mu*Mtot] = 0.0;
	for (int j=2; j<=ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  tmp[nu+mu*Mtot] -= beta0*dt*(V[sigma+mu*Mtot + (t-j+1)*Mtot*Mtot]-V[sigma+mu*Mtot + (t-j)*Mtot*Mtot])*dK[nu+sigma*Mtot+(j)*Mtot*Mtot];
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    tmp[nu+mu*Mtot] -= beta0*dt*(V[sigma+mu*Mtot + (t-j+1)*Mtot*Mtot]-V[sigma+mu*Mtot + (t-j)*Mtot*Mtot])*dK[nu+sigma*Mtot+(j)*Mtot*Mtot];
	  }
#endif
	}
	for (int j=2; j<=t-ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + (j)*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t-j+1)]-K[nu+sigma*Mtot + Mtot*(t-j)*Mtot]);
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + (j)*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t-j+1)]-K[nu+sigma*Mtot + Mtot*(t-j)*Mtot]);
	  }
#endif
	}
	tmp_save1[nu+mu*Mtot] = tmp[nu+mu*Mtot];
      } else {
	tmp[nu+mu*Mtot] = tmp_save1[nu+mu*Mtot];
      }
      
      
      
#ifdef DIAGONAL_APPROXIMATION1
      int sigma=mu;
      tmp[nu+mu*Mtot] -=  DInv[sigma+mu*Mtot]*(K[nu+sigma*Mtot+t*Mtot*Mtot]-4.0/3.0*K[nu+sigma*Mtot+(t-1)*Mtot*Mtot]+1.0/3.0*K[nu+sigma*Mtot+(t-2)*Mtot*Mtot]);
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] -=  DInv[sigma+mu*Mtot]*(K[nu+sigma*Mtot+t*Mtot*Mtot]-4.0/3.0*K[nu+sigma*Mtot+(t-1)*Mtot*Mtot]+1.0/3.0*K[nu+sigma*Mtot+(t-2)*Mtot*Mtot]);
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] += beta0*dt*V[sigma+mu*Mtot+(t-1)*Mtot*Mtot]*dK[nu+sigma*Mtot + Mtot*Mtot];
      tmp[nu+mu*Mtot] -= beta0*dt*V[sigma+mu*Mtot + (t-ibar)*Mtot*Mtot]*K[nu+sigma*Mtot+(ibar)*Mtot*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] += beta0*dt*V[sigma+mu*Mtot+(t-1)*Mtot*Mtot]*dK[nu+sigma*Mtot + Mtot*Mtot];
        tmp[nu+mu*Mtot] -= beta0*dt*V[sigma+mu*Mtot + (t-ibar)*Mtot*Mtot]*K[nu+sigma*Mtot+(ibar)*Mtot*Mtot];
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t)]-K[nu+sigma*Mtot + Mtot*(t-1)*Mtot]);
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t)]-K[nu+sigma*Mtot + Mtot*(t-1)*Mtot]);
      }
#endif
      
      if (mu==nu) tmp[nu+mu*Mtot] += beta0*dt;
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  // include inverse
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    val_t1[mu+mu*Mtot] = dKInv[mu+mu*Mtot]*tmp[mu+mu*Mtot];
#else
    for (int nu=0; nu<Mtot; nu++) {
      val_t1[nu+mu*Mtot] = 0.0;
      for (int sigma=0; sigma<Mtot; sigma++) {
        val_t1[nu+mu*Mtot] += tmp[sigma+mu*Mtot]*dKInv[nu+sigma*Mtot];
      }
    }
#endif
  }
  
  
  // interpolate for predictive move
    for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      K[nu+mu*Mtot+(t+1)*Mtot*Mtot] = K[nu+mu*Mtot+(t)*Mtot*Mtot] + (K[nu+mu*Mtot+(t)*Mtot*Mtot] - K[nu+mu*Mtot+(t-1)*Mtot*Mtot])*dt;
      V[nu+mu*Mtot + t*Mtot*Mtot] = val_t1[nu+mu*Mtot];
      #ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  
  
  // calculate integration step 2
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      
      if (step ==0) {
	tmp[nu+mu*Mtot] = 0.0;
	for (int j=2; j<=ibar1; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  tmp[nu+mu*Mtot] -= beta0*dt*(V[sigma+mu*Mtot + (t+1-j+1)*Mtot*Mtot]-V[sigma+mu*Mtot + (t+1-j)*Mtot*Mtot])*dK[nu+sigma*Mtot+(j)*Mtot*Mtot];
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    tmp[nu+mu*Mtot] -= beta0*dt*(V[sigma+mu*Mtot + (t+1-j+1)*Mtot*Mtot]-V[sigma+mu*Mtot + (t+1-j)*Mtot*Mtot])*dK[nu+sigma*Mtot+(j)*Mtot*Mtot];
	  }
#endif
	}
	for (int j=3; j<=t+1-ibar1; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + (j)*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t+1-j+1)]-K[nu+sigma*Mtot + Mtot*(t+1-j)*Mtot]);
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + (j)*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t+1-j+1)]-K[nu+sigma*Mtot + Mtot*(t+1-j)*Mtot]);
	  }
#endif
	}
      	tmp_save2[nu+mu*Mtot] = tmp[nu+mu*Mtot];
      } else {
	tmp[nu+mu*Mtot] = tmp_save2[nu+mu*Mtot];
      }
      
      
      
#ifdef DIAGONAL_APPROXIMATION1
      int sigma=mu;
      tmp[nu+mu*Mtot] -=  DInv[sigma+mu*Mtot]*(K[nu+sigma*Mtot+(t+1)*Mtot*Mtot]-4.0/3.0*K[nu+sigma*Mtot+(t)*Mtot*Mtot]+1.0/3.0*K[nu+sigma*Mtot+(t-1)*Mtot*Mtot]);
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] -=  DInv[sigma+mu*Mtot]*(K[nu+sigma*Mtot+(t+1)*Mtot*Mtot]-4.0/3.0*K[nu+sigma*Mtot+(t)*Mtot*Mtot]+1.0/3.0*K[nu+sigma*Mtot+(t-1)*Mtot*Mtot]);
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] += beta0*dt*V[sigma+mu*Mtot+(t)*Mtot*Mtot]*dK[nu+sigma*Mtot + Mtot*Mtot];
      tmp[nu+mu*Mtot] -= beta0*dt*V[sigma+mu*Mtot + (t+1-ibar1)*Mtot*Mtot]*K[nu+sigma*Mtot+(ibar1)*Mtot*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] += beta0*dt*V[sigma+mu*Mtot+(t)*Mtot*Mtot]*dK[nu+sigma*Mtot + Mtot*Mtot];
        tmp[nu+mu*Mtot] -= beta0*dt*V[sigma+mu*Mtot + (t+1-ibar1)*Mtot*Mtot]*K[nu+sigma*Mtot+(ibar1)*Mtot*Mtot];
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t+1)]-K[nu+sigma*Mtot + Mtot*(t)*Mtot]);
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t+1)]-K[nu+sigma*Mtot + Mtot*(t)*Mtot]);
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + 2*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t)]-K[nu+sigma*Mtot + Mtot*(t-1)*Mtot]);
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        tmp[nu+mu*Mtot] -= beta0*dt*dV[sigma+mu*Mtot + 2*Mtot*Mtot]*(K[nu+sigma*Mtot+ Mtot*Mtot*(t)]-K[nu+sigma*Mtot + Mtot*(t-1)*Mtot]);
      }
#endif
      
      if (mu==nu) tmp[nu+mu*Mtot] += beta0*dt;
      
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  // include inverse
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    val_t2[mu+mu*Mtot] = dKInv[mu+mu*Mtot]*tmp[mu+mu*Mtot];
#else
    for (int nu=0; nu<Mtot; nu++) {
      val_t2[nu+mu*Mtot] = 0.0;
      for (int sigma=0; sigma<Mtot; sigma++) {
        val_t2[nu+mu*Mtot] += tmp[sigma+mu*Mtot]*dKInv[nu+sigma*Mtot];
      }
    }
#endif
  }
  
  
  // interpolate
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      V[nu+mu*Mtot + t*Mtot*Mtot] = (V[nu+mu*Mtot + (t-1)*Mtot*Mtot] + 3.0*val_t1[nu+mu*Mtot] + val_t2[nu+mu*Mtot])/5.0;
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  
 
    free_dvector(tmp,0,Mtot*Mtot-1);
     free_dvector(val_t1,0,Mtot*Mtot-1);
      free_dvector(val_t2,0,Mtot*Mtot-1);
  free_dvector(tmp_D,0,Mtot*Mtot-1);
  free_dvector(DInv,0,Mtot*Mtot-1);
  
}

/**********************************************************************************************************************************************************/

void solve_volterra_integral_Meff_to_Phi(int M,int Nt, int t, int q,double dt, double q0, double dq, int step,  double * Q, double * n, double * dKInv,double * tmp_save, double * K, double * dK, double * V , double * dV) {
  int Mtot = 2*M+1;
  double * tmp;
  tmp = dvector(0,Mtot*Mtot-1);
  int ibar = t/2;
  double qval = q0 + q*dq;
  
  // calculate inverse of D
  double * tmp_D;
  tmp_D = dvector(0,Mtot*Mtot-1);
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      tmp_D[mu*Mtot+nu] =  n[mu -nu + 3*M]/n[3*M] * (qval * qval + Q[mu+2*M]*Q[nu+2*M] );
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }
  
  // calculate integration step
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
  int nu=mu;
#else
    for (int nu=0; nu<Mtot; nu++) {
#endif
      

      if (step==0) {
	tmp[nu+mu*Mtot] = 0.0;
	// include convolution
	for (int j=2; j<=ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int sigma=mu;
	  int kappa=mu;
	  tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*(K[kappa+sigma*Mtot+ Mtot*Mtot*(t-j+1)]-K[kappa+sigma*Mtot+ Mtot*Mtot*(t-j)])*dV[nu+kappa*Mtot+(j)*Mtot*Mtot];
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    for (int kappa=0; kappa<Mtot; kappa++) {
	      tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*(K[kappa+sigma*Mtot+ Mtot*Mtot*(t-j+1)]-K[kappa+sigma*Mtot+ Mtot*Mtot*(t-j)])*dV[nu+kappa*Mtot+(j)*Mtot*Mtot];
	    }
	  }
#endif
	  
	}

	for (int j=2; j<=t-ibar; j++) {
#ifdef DIAGONAL_APPROXIMATION1
	int sigma=mu;
	int kappa=mu;
	tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*dK[kappa+sigma*Mtot+(j)*Mtot*Mtot]*(V[nu+kappa*Mtot + Mtot*Mtot*(t-j+1)]-V[nu+kappa*Mtot + Mtot*Mtot*(t-j)]);
#else
	  for (int sigma=0; sigma<Mtot; sigma++) {
	    for (int kappa=0; kappa<Mtot; kappa++) {
	      tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*dK[kappa+sigma*Mtot+(j)*Mtot*Mtot]*(V[nu+kappa*Mtot + Mtot*Mtot*(t-j+1)]-V[nu+kappa*Mtot + Mtot*Mtot*(t-j)]);
	    }
	  }
#endif
	}
        tmp_save[nu+mu*Mtot] = tmp[nu+mu*Mtot];
      } else {
	tmp[nu+mu*Mtot] = tmp_save[nu+mu*Mtot];
      }

      // include previous timesteps
      tmp[nu+mu*Mtot] += alpha1*V[nu+mu*Mtot+(t-1)*Mtot*Mtot];
      tmp[nu+mu*Mtot] += alpha2*V[nu+mu*Mtot+(t-2)*Mtot*Mtot];
      
      // include additional terms
#ifdef DIAGONAL_APPROXIMATION1
      int sigma=mu;
      int kappa=mu;
      tmp[nu+mu*Mtot] += beta0*dt*tmp_D[sigma+mu*Mtot]*dK[kappa+sigma*Mtot+Mtot*Mtot]*V[nu+kappa*Mtot+(t-1)*Mtot*Mtot];
      tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*K[kappa+sigma*Mtot+Mtot*Mtot*(t-ibar)]*V[nu+kappa*Mtot+(ibar)*Mtot*Mtot];
      //tmp[a] -= beta0*dt/2.0*D*K[a+size*(ibar)]*V[a+(t-ibar)*size];
      tmp[nu+mu*Mtot] += beta0*dt*tmp_D[sigma+mu*Mtot]*K[kappa+sigma*Mtot+Mtot*Mtot*(t)]*V[nu+kappa*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
        for (int kappa=0; kappa<Mtot; kappa++) {
          tmp[nu+mu*Mtot] += beta0*dt*tmp_D[sigma+mu*Mtot]*dK[kappa+sigma*Mtot+Mtot*Mtot]*V[nu+kappa*Mtot+(t-1)*Mtot*Mtot];
          tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*K[kappa+sigma*Mtot+Mtot*Mtot*(t-ibar)]*V[nu+kappa*Mtot+(ibar)*Mtot*Mtot];
          //tmp[a] -= beta0*dt/2.0*D*K[a+size*(ibar)]*V[a+(t-ibar)*size];
          tmp[nu+mu*Mtot] += beta0*dt*tmp_D[sigma+mu*Mtot]*K[kappa+sigma*Mtot+Mtot*Mtot*(t)]*V[nu+kappa*Mtot];
        }
      }
#endif

#ifdef DIAGONAL_APPROXIMATION1
      int sigma=mu;
      int kappa=mu;
      tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*(K[kappa+sigma*Mtot+ Mtot*Mtot*(t)]-K[kappa+sigma*Mtot+ Mtot*Mtot*(t-1)])*dV[nu+kappa*Mtot+Mtot*Mtot];
#else
      for (int sigma=0; sigma<Mtot; sigma++) {
	for (int kappa=0; kappa<Mtot; kappa++) {
	  tmp[nu+mu*Mtot] -= beta0*dt*tmp_D[sigma+mu*Mtot]*(K[kappa+sigma*Mtot+ Mtot*Mtot*(t)]-K[kappa+sigma*Mtot+ Mtot*Mtot*(t-1)])*dV[nu+kappa*Mtot+Mtot*Mtot];
	}
      }
#endif
      
#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  // include inverse
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
    V[mu+mu*Mtot + Mtot*Mtot*t] = dKInv[mu+mu*Mtot]*tmp[mu+mu*Mtot];
    //if(q==6 && mu==M) printf("%.15e %.15e %.15e\n",tmp[mu+mu*Mtot],dKInv[mu+mu*Mtot] ,V[mu+mu*Mtot + Mtot*Mtot*t]);

#else
    for (int nu=0; nu<Mtot; nu++) {
      V[nu+mu*Mtot + Mtot*Mtot*t] = 0.0;
      for (int sigma=0; sigma<Mtot; sigma++) {
        V[nu+mu*Mtot + Mtot*Mtot*t] += dKInv[sigma+mu*Mtot]*tmp[nu+sigma*Mtot];
      }
    }
#endif
  }
  
  free_dvector(tmp,0,Mtot*Mtot-1);
  free_dvector(tmp_D,0,Mtot*Mtot-1);

}

/**********************************************************************************************************************************************************/

/*void solve_volterra_integral_K_to_Phi(int Nt, int t, int q, int size, double dt, double ** phi, double * dKInv, double * K, double * dK, double * V , double * dV) {
  
  double * tmp;
  tmp = dvector(0,size*size-1);
  
  // calculate inverse of phi
  double * phiInv;
  phiInv = dvector(0,size*size-1);
  calc_matrixInv(&phi[q][0],size,phiInv);

  // calculate integration step
  for (int a=0; a<size; a++) {
#ifdef DIAGONAL_APPROXIMATION1
    int b=a;
#else
    for (int b=0; b<size; b++) {
#endif
      // include previous timesteps
      tmp[b+a*size] = alpha1*V[b+a*size+(t-1)*size*size];
      tmp[b+a*size] += alpha2*V[b+a*size+(t-2)*size*size];
      
#ifdef DIAGONAL_APPROXIMATION1
      int c = a;
      int d = a;
      tmp[b+a*size] -= beta0*dt*dt/2.0*dK[c+a*size+size*size]*phiInv[d+c*size]*V[b+d*size+(t-1)*size*size];
#else
      for (int c=0; c<size; c++) {
        for (int d=0; d<size; d++) {
          tmp[b+a*size] -= beta0*dt*dt/2.0*dK[c+a*size+size*size]*phiInv[d+c*size]*V[b+d*size+(t-1)*size*size];
        }
      }
#endif
    
      // include convolution
      for (int j=1; j<t-Nt/2+1; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int c = a;
        int d = a;
        tmp[b+a*size] -= beta0*dt*dt/2.0*(K[c+a*size + size*size*(t-j)]+K[c+a*size + size*size*(t-j+1)])*phiInv[d+c*size]*dV[b+d*size+(j)*size*size];
#else
        for (int c=0; c<size; c++) {
          for (int d=0; d<size; d++) {
            tmp[b+a*size] -= beta0*dt*dt/2.0*(K[c+a*size + size*size*(t-j)]+K[c+a*size + size*size*(t-j+1)])*phiInv[d+c*size]*dV[b+d*size+(j)*size*size];
          }
        }
#endif
      }
      
      for (int j=2; j<=Nt/2; j++) {
#ifdef DIAGONAL_APPROXIMATION1
        int c = a;
        int d = a;
        tmp[b+a*size] -= beta0*dt*dt/2.0*dK[c+a*size+(j)*size*size]*phiInv[d+c*size]*(V[b+d*size + size*size*(t-j)]+V[b+d*size + size*size*(t-j+1)]);
#else
        for (int c=0; c<size; c++) {
          for (int d=0; d<size; d++) {
            tmp[b+a*size] -= beta0*dt*dt/2.0*dK[c+a*size+(j)*size*size]*phiInv[d+c*size]*(V[b+d*size + size*size*(t-j)]+V[b+d*size + size*size*(t-j+1)]);
          }
        }
#endif
      }


#ifndef DIAGONAL_APPROXIMATION1
    }
#endif
  }

  
  for (int a=0; a<size; a++) {
#ifdef DIAGONAL_APPROXIMATION1
    V[a+a*size + size*size*t] = dKInv[a+a*size]*tmp[a+a*size];
#else
    for (int b=0; b<size; b++) {
      V[b+a*size + size*size*t] = 0.0;
      for (int c=0; c<size; c++) {
        V[b+a*size + size*size*t] += dKInv[c+a*size]*tmp[b+c*size];
      }
    }
#endif
  }
  
  free_dvector(tmp,0,size*size-1);
  free_dvector(phiInv,0,size*size-1);
}   */
