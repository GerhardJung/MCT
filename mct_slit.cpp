/* program to evaluate MCT memory kernels for the slit geometry
 * author: Gerhard Jung (UIBK)
 * date: 18.03.2019 
 */

#include "mct_slit.h"
#include "calc_kernel.h"
#include "time_integration.h"
#include "nrutil.h"
#include "svd.h"
#include "def.h"

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <sstream>
#include <omp.h>
#include <sys/time.h> 


int main(int argc, char *argv[]) {
  
  if (argc < 4) {
    printf("COMMAND-LINE ARGUMENTS:\nN [Number of wavevectors]\nM [Number of mode indices]\nL [Accesible channel widths]\nPhi [Volume fraction]\n");
    exit(0);
  }
  
  int Nt = 1024;				// number of time steps per decimation step
  int D = 100;				// number of decimation steps
  
  int N = atoi(argv[1]);				// number of wavevectors
  int M = atoi(argv[2]);				// number of modes
  double L = atof(argv[3]);			// channel widths (accesible)
  double phi_frac = atof(argv[4]);		// packing fraction
  double n0 = 6.0*phi_frac/PI*(L+1.0);	// density
  printf("n0 %.10f\n",n0);
  
  double qmax= 30.0;			// maximum wavevector
  double q0 = 0.0;			// minimum wavevector
  double dq = (qmax-q0)/(N-1);	// discretization
  q0 = 0.303*dq;			// redefinition
  qmax = q0+(N-1)*dq;		// redefinition
  printf("Wavenumbers: qmin = %.4f ,   qmax = %.4f ,   dq = %.4f \n",q0,qmax, dq);
  
  double D0 = 1.0;
  
  double dt = 1.0e-9;
  double eps = 1e-7;
  double maxeps_save = 100;
  int maxstep = 100; 
  const char* density_file = "rho.dat";
  int Dswap=18;
  int Dlinear1=22;
  int Dlinear2=24;
  int Dlinear3=31;
  int Nt_max = 8192;
  
  double * Q;
  double * n; double * v;
  double ** c; double ** phi; double ** phiold; double ** m; double ** K; double ** Kr; double ** Kc; double ** Krc; double ** Meff; double ** Meff_long;
  // moments
  double ** dphi; double ** dm; double ** dK;  double ** dKr; double ** dKc; double ** dKrc; double ** dMeff; double ** dMeff_long;
  
  // precalculate inner sums (in O(N^2*M^4))
  double ** YS; double ** YM;
  double ** YSQ; double ** YMQ;
  
  // inverse for time integration
  double ** M_to_K_Inv; double ** M_to_Kr_Inv; double ** K_to_Meff_Inv; double ** Kr_to_Meff_Inv; double ** Meff_to_Phi_Inv; double ** K_to_Phi_Inv;
  double ** M_to_K_save; double ** M_to_Kr_save; double ** K_to_Meff_save; double ** Kr_to_Meff_save1; double ** Kr_to_Meff_save2; double ** Meff_to_Phi_save;
  
  double * mathcalD;
  double * mathcalDInv;
  double ** DInv; 
  double ** Dif; 
  double ** DphiInv;
  
  FILE * out_coherent = fopen("coherent_scattering_function.dat","w");

  int Mtot = 2*M+1;
  int Mtot3 = 6*M +1; //Mtot3 is neccessary to get as much information about the modes as possible
  Q = dvector(0,Mtot3-1);
  n = dvector(0,Mtot3-1);
  v = dvector(0,Mtot3-1);
  c = dmatrix(0,N-1,0, Mtot*Mtot-1);
  phi = dmatrix(0,N-1,0, (Nt+1)*Mtot*Mtot-1);
  phiold = dmatrix(0,N-1,0, Mtot*Mtot-1);
  m = dmatrix(0,N-1,0, (Nt+1)*4*Mtot*Mtot-1);
  K = dmatrix(0,N-1,0, (Nt+1)*4*Mtot*Mtot-1);
  Kr = dmatrix(0,N-1,0, (Nt+1)*4*Mtot*Mtot-1);
  Kc = dmatrix(0,N-1,0, (Nt+1)*Mtot*Mtot-1);
  Krc = dmatrix(0,N-1,0, (Nt+1)*Mtot*Mtot-1);
  Meff = dmatrix(0,N-1,0, (Nt+1)*Mtot*Mtot-1);
  Meff_long = dmatrix(0,N-1,0, (Nt+1)*Mtot*Mtot-1);
  
  dphi = dmatrix(0,N-1,0, (Nt/2+2)*Mtot*Mtot-1);
  dm = dmatrix(0,N-1,0, (Nt/2+2)*4*Mtot*Mtot-1);
  dK = dmatrix(0,N-1,0, (Nt/2+2)*4*Mtot*Mtot-1);
  dKr = dmatrix(0,N-1,0, (Nt/2+2)*4*Mtot*Mtot-1);
  dKc = dmatrix(0,N-1,0, (Nt/2+2)*Mtot*Mtot-1);
  dKrc = dmatrix(0,N-1,0, (Nt/2+2)*Mtot*Mtot-1);
  dMeff = dmatrix(0,N-1,0, (Nt/2+2)*Mtot*Mtot-1);
  dMeff_long = dmatrix(0,N-1,0, (Nt/2+2)*Mtot*Mtot-1);
  
  YS=dmatrix(0,N-1,0, Mtot*Mtot*Mtot-1);
  YM=dmatrix(0,N*N-1,0, Mtot*Mtot*Mtot-1);
  YSQ=dmatrix(0,N-1,0, Mtot*Mtot*Mtot-1);
  YMQ=dmatrix(0,N*N-1,0, Mtot*Mtot*Mtot-1);
  
  M_to_K_Inv = dmatrix(0,N-1,0, 4*Mtot*Mtot-1);
  M_to_Kr_Inv = dmatrix(0,N-1,0, 4*Mtot*Mtot-1);
  K_to_Meff_Inv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Kr_to_Meff_Inv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Meff_to_Phi_Inv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  
  M_to_K_save = dmatrix(0,N-1,0, 4*Mtot*Mtot-1);
  M_to_Kr_save = dmatrix(0,N-1,0, 4*Mtot*Mtot-1);
  K_to_Meff_save = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Kr_to_Meff_save1 = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Kr_to_Meff_save2 = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Meff_to_Phi_save = dmatrix(0,N-1,0, Mtot*Mtot-1);
  //K_to_Phi_Inv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  
  mathcalD = dvector(0, 4*Mtot*Mtot-1);
  mathcalDInv = dvector(0, 4*Mtot*Mtot-1);
  DInv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  Dif = dmatrix(0,N-1,0, Mtot*Mtot-1);
  DphiInv = dmatrix(0,N-1,0, Mtot*Mtot-1);
  
  // initialize dicretizaed modes
  for (int mu=0; mu<Mtot3; mu++) {
    Q[mu] = 2.0*PI/L*((double) mu - 3*M);
  }
  
  // initialize and read static structure factor ,direct correlation function and density modes
  printf("Initialization...\n");
  read_c_phi_modes(Nt,N,M,q0,dq,c,phi,dphi);
  read_n_modes(M,L,n,v,density_file);
  initInv(N, M, n0,q0,dq,D0, Q, n, phi,Dif,DInv,DphiInv,mathcalD,mathcalDInv);
  printf("Initialization... Done!\n");
#ifdef NONERGODICITY 
  int t=0, d=0;
  
#else // (for full time-dependence)
  // initialize values and moments for t<=Nt/2
  init_values_moments(Nt, N, M, n0, L,q0,dq,dt,D0,Q, v, n, c, mathcalD, DInv, YS, YM,YSQ, YMQ, phi, m ,K ,Kr ,Kc , Krc, Meff, Meff_long, dm ,dK ,dKr ,dKc,dKrc,dMeff, dMeff_long);
  
  // initialize inverse for time integration
  calc_dKInv(N, 2*Mtot, dt,mathcalD,Dif, DInv, DphiInv,phi, dm , M_to_K_Inv, 0);
  calc_dKInv(N, 2*Mtot, dt,mathcalD,Dif, DInv, DphiInv,phi, dm , M_to_Kr_Inv, 1);
  calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dKc , K_to_Meff_Inv, 2);
  calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dKrc , Kr_to_Meff_Inv, 3);
  calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dMeff , Meff_to_Phi_Inv, 4);
  //calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv,phi, dKc , K_to_Phi_Inv, 5);
  
  //printf("%.15e %.15e %.15e %.15e\n",M_to_K_Inv[6][M+M*2*Mtot],M_to_K_Inv[6][M+Mtot+(M+Mtot)*2*Mtot] ,K_to_Meff_Inv[6][M+M*Mtot],Meff_to_Phi_Inv[6][M+M*Mtot]);

  for (int d=0; d<D; d++) { // the main decimation loop: starting simulation
    
    printf("Start Decimation Loop %d\n",d);
    
    for (int t=Nt/2; t<Nt; t++) { // the time dependence loop for one decimation step
      
      if (t%100==0) printf("Start Time Step %d\n",t);
      
      struct timeval start_tot, end_tot; 
	gettimeofday(&start_tot, NULL); 
      
      // initialize new phi
      for (int q=0; q<N; q++) { 
        for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
          int nu = mu;
#else
          for (int nu=0; nu<Mtot; nu++) {
#endif
            phi[q][t*Mtot*Mtot+mu*Mtot+nu]=phi[q][(t-1)*Mtot*Mtot+mu*Mtot+nu];
#ifndef DIAGONAL_APPROXIMATION1
          }
#endif
        }
      }
#endif //NONERGODICITY
      
      int step = 0;
      int continue_it = 1;
      do { // solve the self-consistent equation (see derivation_time_solver.pdf) 
	
	
#ifndef DIAGONAL_APPROXIMATION1
#ifdef NONERGODICITY
        calc_mkernel_init(t, N, M, n0, L, Q, c, phi,YS, YM, YSQ, YMQ);
#else
        calc_mkernel_init(t, N, M, n0, L, Q, c, phi,YS, YM, YSQ, YMQ);
#endif
#endif
	
	// the main computations (performed independently for all wavevectors)
	  
#ifdef NONERGODICITY
        
        //printf("Starting Iteration!\n");
        
        for (int q=0; q<N; q++) {
          calc_mkernel(t, N, q,  M, n0,L, q0, dq, c, Q,n,YS,YM,YSQ,YMQ, phi, &m[q][0]);
        }
        
#ifdef DIAGONAL_APPROXIMATION1
        // first need to extrapolate q=0, mu = M, a=1 value, since it is inherently 0!!
        m[0][M + Mtot + M*2*Mtot + 2*Mtot*Mtot]=2.0*m[1][M + Mtot + M*2*Mtot + 2*Mtot*Mtot] - m[2][M + Mtot + 2*M*Mtot + 2*Mtot*Mtot];
        // calc inverse
        for (int q=0; q<N; q++) {	
          for (int mu=0; mu<Mtot; mu++) {
            for (int a=0; a<2; a++) {
              K[q][mu + a*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot] = 1.0/(m[q][mu + a*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot]);
            }
          }
        }
#else
        // Invert m (for Eq. 217)
        double ** mtmp;
        double ** mInvtmp;
        int size = 2*Mtot;
        mtmp = dmatrix(0, size-1,0,size-1);
        mInvtmp = dmatrix(0, size-1,0,size-1);
        for (int q=0; q<N; q++) {
          for (int mu=0; mu<2*Mtot; mu++) {
            for (int nu=0; nu<2*Mtot; nu++) {
	      mtmp[mu][nu] =  0.0;
	      for (int sigma=0; sigma<2*Mtot; sigma++) {
		mtmp[mu][nu] +=  mathcalDInv[mu*2*Mtot+sigma]*m[q][sigma*2*Mtot+nu];
	      }
              //if (q==0) printf("%d %d %f\n",mu,nu,mtmp[mu][nu]);
            }
          }
          invcmp2(mtmp, size, mInvtmp);
          for (int a=0; a<size; a++) {
            for (int b=0; b<size; b++) {
              K[q][b + a*size] = mInvtmp[a][b];
              //if (q==0) printf("%f\n",K[q][b + a*size]);
            }
          }	
        }
        free_dmatrix(mtmp,0, size-1,0,size-1);
        free_dmatrix(mInvtmp,0, size-1,0,size-1);
        
        //printf("M inverted!\n");

        
#endif //DIAGONAL_APPROXIMATION1
#else //TIME_DEPENDENCE
	
	struct timeval start, end; 
	gettimeofday(&start, NULL); 
	
	#pragma omp parallel for
        for (int q=0; q<N; q++) {
          calc_mkernel(t, N, q,  M, n0,L, q0, dq, c, Q,n,YS,YM,YSQ,YMQ, phi, &m[q][4*Mtot*Mtot*t]);
        }
        
        gettimeofday(&end, NULL); 
	double time_taken; 
	time_taken = (end.tv_sec - start.tv_sec) * 1e6; 
	time_taken = (time_taken + (end.tv_usec -  
                              start.tv_usec)) * 1e-6; 
	printf("Finished kernel calculation ... Optimized! Time: %f \n",time_taken);

		
	struct timeval start_int, end_int; 
	gettimeofday(&start_int, NULL); 
	
	#pragma omp parallel for
        for (int q=0; q<N; q++) { 
	  if (d<=Dswap)  solve_volterra_integral_M_to_K(M,Nt, t, q, dt,q0,  dq, step, Q,n,mathcalDInv,&M_to_K_Inv[q][0],&M_to_K_save[q][0],&m[q][0], &dm[q][0], &K[q][0], &dK[q][0]); 
          solve_volterra_integral_M_to_Kr(M,Nt, t, q, dt,q0,  dq, step, Q,n,mathcalDInv,&M_to_Kr_Inv[q][0],&M_to_Kr_save[q][0],&m[q][0], &dm[q][0], &Kr[q][0], &dKr[q][0]);  
	  //solve_volterra_integral_M_to_Kr_predictive(M,Nt, t, q, dt,q0,  dq, Q,n,mathcalDInv,&M_to_Kr_Inv[q][0],&m[q][0], &dm[q][0], &Kr[q][0], &dKr[q][0]);
        }
        
        gettimeofday(&end_int, NULL); 
	double time_taken_int; 
	time_taken_int = (end_int.tv_sec - start_int.tv_sec) * 1e6; 
	time_taken_int = (time_taken_int + (end_int.tv_usec -  start_int.tv_usec)) * 1e-6; 
	printf("Finished time integration Time: %f \n",time_taken_int);
	
#endif	//NONERGODICITY

#ifndef NONERGODICITY	
	gettimeofday(&start_int, NULL); 
#endif
		
        for (int q=0; q<N; q++) {
          if (d<=Dswap)  contract_K(q, M, q0, dq,L,n0,Q, v,c, &K[q][4*Mtot*Mtot*t], &Kc[q][Mtot*Mtot*t]);
          contract_K(q, M, q0, dq,L,n0,Q, v,c, &Kr[q][4*Mtot*Mtot*t], &Krc[q][Mtot*Mtot*t]);
        } //q-loop
        

        
	#pragma omp parallel for
        for (int q=0; q<N; q++) { 
          if (d <= Dswap) solve_volterra_integral_K_to_Meff(M,Nt, t, q, dt,q0,  dq, step, Q,n,&K_to_Meff_Inv[q][0],&K_to_Meff_save[q][0],&Kc[q][0], &dKc[q][0], &Meff[q][0], &dMeff[q][0]);
          if (d > Dswap - 3) solve_volterra_integral_Kr_to_Meff(M,Nt, t, q, dt,q0,  dq, step,Q,n,&Kr_to_Meff_Inv[q][0],&Kr_to_Meff_save1[q][0],&Kr_to_Meff_save2[q][0],&Krc[q][0], &dKrc[q][0], &Meff_long[q][0], &dMeff_long[q][0]);
        }
        
        if (d>Dswap) {
          for (int q=0; q<N; q++) { 
            for (int mu=0; mu<Mtot; mu++) {
              for (int nu=0; nu<Mtot; nu++) {
                Meff[q][nu+mu*Mtot+t*Mtot*Mtot] = Meff_long[q][nu+mu*Mtot+t*Mtot*Mtot];
              }
            }
          }
        }
        if (d<=Dswap - 3) {
	  for (int q=0; q<N; q++) { 
            for (int mu=0; mu<Mtot; mu++) {
              for (int nu=0; nu<Mtot; nu++) {
                Meff_long[q][nu+mu*Mtot+t*Mtot*Mtot] = Meff[q][nu+mu*Mtot+t*Mtot*Mtot];
              }
            }
          }  
	}
         
        // barrier (not necesarry if no shared memory)
	#pragma omp parallel for
        for (int q=0; q<N; q++) {
          //store phi to calculate convergence
          for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
            int nu = mu;
#else
            for (int nu=0; nu<Mtot; nu++) {
#endif
              phiold[q][mu*Mtot+nu]=phi[q][t*Mtot*Mtot+mu*Mtot+nu];
#ifndef DIAGONAL_APPROXIMATION1
            }
#endif
          }
          
#ifdef NONERGODICITY
          // evaluate Eq. 260 in Simon Lang's Diploma thesis (phi stored in t=1
#ifdef DIAGONAL_APPROXIMATION1
          for (int mu=0; mu < Mtot; mu++) {
            phi[q][mu*Mtot+mu] = phi[q][Mtot*Mtot+mu*Mtot+mu]/(1.0+Kc[q][mu*Mtot+mu]/phi[q][Mtot*Mtot+mu*Mtot+mu]);
          }
#else
          int Mtot = 2*M+1;
          // calculate the inverse of the new phi
          double * KphiInv;
          double ** phiinvtmp;
          double ** phitmp;
          double * phiInv;

          KphiInv = dvector(0, Mtot*Mtot-1);
          phiinvtmp = dmatrix(0, Mtot-1,0, Mtot-1);
          phitmp = dmatrix(0, Mtot-1,0, Mtot-1);
          phiInv = dvector(0, Mtot*Mtot-1);
          
          calc_matrixInv(&phi[q][Mtot*Mtot],Mtot,phiInv);
          
          //printf("Phi inverted!\n");
        
          for (int mu=0; mu<Mtot; mu++) {
            for (int nu=0; nu<Mtot; nu++) {
              KphiInv[mu*Mtot+nu] = 0.0;
              //if (q==0) printf("%f\n",phiInv[mu*Mtot+nu]);
              if (mu==nu) KphiInv[mu*Mtot+nu] += 1.0;
              for (int sigma=0; sigma<Mtot; sigma++) {
                KphiInv[mu*Mtot+nu] += Kc[q][mu*Mtot+sigma]*phiInv[sigma*Mtot+nu];
              }
            }
          }
          
          for (int mu=0; mu<Mtot; mu++) {
            for (int nu=0; nu<Mtot; nu++) {
              phiinvtmp[mu][nu] = KphiInv[nu + mu*Mtot];
              
              
            }
          }
          invcmp2(phiinvtmp, Mtot, phitmp);
          for (int mu=0; mu<Mtot; mu++) {
            for (int nu=0; nu<Mtot; nu++) {
              phi[q][mu*Mtot + nu] = 0.0;
              for (int sigma=0; sigma<Mtot; sigma++) {
                phi[q][ mu*Mtot+ nu] += phitmp[mu][sigma]*phi[q][Mtot*Mtot+sigma*Mtot+nu];
              }
            }
          }
          
          //printf("New phi determined!\n");

#endif //DIAGONAL_APPROXIMATION1
#else //FULL TIME-DEPENDENCE

	  

          solve_volterra_integral_Meff_to_Phi(M,Nt, t, q, dt,q0,  dq, step,Q,n,&Meff_to_Phi_Inv[q][0],&Meff_to_Phi_save[q][0],&Meff[q][0], &dMeff[q][0], &phi[q][0], &dphi[q][0]);     
          

#endif	  
          
          // symmetrize
          for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
            int nu=mu;
#else
            for (int nu=0; nu<Mtot; nu++) {
#endif
              if (nu < mu) phi[q][mu*Mtot + nu+t*Mtot*Mtot]= phi[q][nu*Mtot+mu+t*Mtot*Mtot];
              if (mu > M) phi[q][mu*Mtot + nu+t*Mtot*Mtot]= phi[q][(2*M-mu)*Mtot+ 2*M- nu+t*Mtot*Mtot];
          
#ifndef DIAGONAL_APPROXIMATION1
            }
#endif
          }
          
        } // end-q-loop
#ifndef NONERGODICITY        
	gettimeofday(&end_int, NULL); 
	time_taken_int = (end_int.tv_sec - start_int.tv_sec) * 1e6; 
	time_taken_int = (time_taken_int + (end_int.tv_usec -  start_int.tv_usec)) * 1e-6; 
	printf("Finished time integration Time 2: %f \n",time_taken_int);
#endif
	
// barrier
	
        // check whether self-iteration loop can be exited 
        double maxeps = 0.0;
        double max = 0.0;
        int max_mu= -1;
        int max_nu= -1;
        int max_q = -1;
        for (int q=0; q<N; q++) { 
          for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
            int nu = mu;
#else
            for (int nu=0; nu<Mtot; nu++) {
#endif
              maxeps = MAX(maxeps, ABS((phiold[q][mu*Mtot+nu]-phi[q][t*Mtot*Mtot+mu*Mtot+nu])));
              if (max < phi[q][t*Mtot*Mtot+mu*Mtot+nu]) {
                max_mu = mu;
                max_nu = nu;
                max_q=q;
              }
              max = MAX(max, ABS(phi[q][t*Mtot*Mtot+mu*Mtot+nu]));
      #ifndef DIAGONAL_APPROXIMATION1
            }
      #endif
          }
        }
        if (t%1==0) printf("It. Step %d (t=%f): Maxeps = %.20f, Max = %.20f (q=%d, mu/nu=%d/%d) \n", step,t*dt, maxeps,max,max_q,max_mu,max_nu);
        
        //printf("%f %f %f %f\n",phi[6][t*Mtot*Mtot+3*Mtot+3],phi[6][t*Mtot*Mtot+4*Mtot+4],phi[6][t*Mtot*Mtot+5*Mtot+5],phi[6][t*Mtot*Mtot+6*Mtot+6]);
        
        if (maxeps < eps || maxeps > maxeps_save) {
	  
	  continue_it = 0;
	  
#ifdef NONERGODICITY
	  FILE * out_state = fopen("state.dat","w");
	  if (max > 0.01 && maxeps < maxeps_save) {
	    fprintf(out_state,"Glass state\n");
	  } else {
	    fprintf(out_state,"Liquid state\n");
	  }
	  fclose(out_state);
#endif
	  
	} else {
	  
	  maxeps_save = maxeps;
	}
        
        // update step
        step++;
	
	if (step > maxstep) exit(0);
	
      } while (continue_it); // solve the self-consistent equation, finish if solution is self-consistent
      
      maxeps_save = 100;
#ifndef NONERGODICITY      
      gettimeofday(&end_tot, NULL); 
	double time_taken_tot; 
	time_taken_tot = (end_tot.tv_sec - start_tot.tv_sec) * 1e6; 
	time_taken_tot = (time_taken_tot + (end_tot.tv_usec -  
                              start_tot.tv_usec)) * 1e-6; 
	printf("Finished full step! Time: %f \n",time_taken_tot);
#endif
      
      if (t % 1 == 0) {
        fprintf(out_coherent,"%.15e ",t*dt);
        int q=6;
        int mu=3;
        int nu = mu;

        fprintf(out_coherent,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",phi[q][mu*Mtot+nu+t*Mtot*Mtot],Meff[q][mu*Mtot+nu+t*Mtot*Mtot],Meff_long[q][mu*Mtot+nu+t*Mtot*Mtot],Kc[q][mu*Mtot+nu+t*Mtot*Mtot],Krc[q][mu*Mtot+nu+t*Mtot*Mtot],K[q][mu*2*Mtot+nu+t*4*Mtot*Mtot],K[q][mu*2*Mtot+nu+2*Mtot*Mtot+t*4*Mtot*Mtot],K[q][mu*2*Mtot+nu+Mtot+t*4*Mtot*Mtot],K[q][mu*2*Mtot+nu+Mtot+2*Mtot*Mtot+t*4*Mtot*Mtot],Kr[q][mu*2*Mtot+nu+t*4*Mtot*Mtot],Kr[q][mu*2*Mtot+nu+2*Mtot*Mtot+t*4*Mtot*Mtot],Kr[q][mu*2*Mtot+nu+Mtot+t*4*Mtot*Mtot],Kr[q][mu*2*Mtot+nu+Mtot+2*Mtot*Mtot+t*4*Mtot*Mtot],m[q][mu*2*Mtot+nu+t*4*Mtot*Mtot],m[q][mu*2*Mtot+nu+2*Mtot*Mtot+t*4*Mtot*Mtot],m[q][mu*2*Mtot+nu+Mtot+t*4*Mtot*Mtot],m[q][mu*2*Mtot+nu+Mtot+2*Mtot*Mtot+t*4*Mtot*Mtot]);
      }
    
      fflush(out_coherent);
      
      // necessary for predictive step
      if (t==Nt/2 || t==Nt/2+1) {
        int indext = t;
        for (int q=0; q<N; q++) {
          for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
            int nu=mu;
            dMeff_long[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Meff_long[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Meff_long[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	    dMeff[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Meff[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Meff[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	    dphi[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(phi[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+phi[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	    if (d<=Dswap)  dKc[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Kc[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Kc[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	    dKrc[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Krc[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Krc[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
#else
            for (int nu=0; nu<Mtot; nu++) {
              dMeff_long[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Meff_long[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Meff_long[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	      dMeff[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Meff[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Meff[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	      dphi[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(phi[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+phi[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	      if (d<=Dswap)  dKc[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Kc[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Kc[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
	      dKrc[q][mu*Mtot+nu+indext*Mtot*Mtot] = 1.0/2.0*(Krc[q][(indext-1)*Mtot*Mtot+mu*Mtot+nu]+Krc[q][(indext)*Mtot*Mtot+mu*Mtot+nu]);
            }
#endif
          }
          
          for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
            int nu=mu;
            dm[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(m[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+m[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
	    if (d<=Dswap)  dK[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(K[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+K[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
	    dKr[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(Kr[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+Kr[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
#else
            for (int nu=0; nu<2*Mtot; nu++) {
              dm[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(m[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+m[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
	      if (d<=Dswap)  dK[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(K[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+K[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
	      dKr[q][mu*2*Mtot+nu+indext*4*Mtot*Mtot] = 1.0/2.0*(Kr[q][(indext-1)*4*Mtot*Mtot+mu*2*Mtot+nu]+Kr[q][(indext)*4*Mtot*Mtot+mu*2*Mtot+nu]);
            }
#endif
          }
        }
      }

#ifndef NONERGODICITY //(only for full time-dependence)
      
    } // the time dependence loop for one decimation step 

    
    if (d != Dlinear1 && d != Dlinear2 && d != Dlinear3  ) {
      // perform decimation, calculate new moments
      printf("Perform Decimation...\n");
      decimation_moments(Nt, N, Mtot, dphi,phi);
      if (d<=Dswap) decimation_moments(Nt, N, Mtot, dKc,Kc);
      decimation_moments(Nt, N, Mtot, dKrc,Krc);
      decimation_moments(Nt, N, Mtot, dMeff,Meff);
      decimation_moments(Nt, N, Mtot, dMeff_long,Meff_long);
      if (d<=Dswap) decimation_moments(Nt, N, 2*Mtot, dK,K);
      decimation_moments(Nt, N, 2*Mtot, dKr,Kr);
      decimation_moments(Nt, N, 2*Mtot, dm,m);
      
      decimation_values(Nt, N, Mtot, phi);
      if (d<=Dswap) decimation_values(Nt, N, Mtot, Kc);
      decimation_values(Nt, N, Mtot, Krc);
      decimation_values(Nt, N, Mtot, Meff);
      decimation_values(Nt, N, Mtot, Meff_long);
      if (d<=Dswap) decimation_values(Nt, N, 2*Mtot, K);
      decimation_values(Nt, N, 2*Mtot, Kr);
      decimation_values(Nt, N, 2*Mtot, m);
      
      dt *= 2.0;
      
      printf("Update Inverse...\n");
      calc_dKInv(N, 2*Mtot, dt,mathcalD, Dif, DInv, DphiInv, phi,dm , M_to_K_Inv, 0);
      calc_dKInv(N, 2*Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dm , M_to_Kr_Inv, 1);
      calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dKc , K_to_Meff_Inv, 2);
      calc_dKInv(N, Mtot, dt,mathcalD, Dif,DInv, DphiInv, phi,dKrc , Kr_to_Meff_Inv, 3);
      calc_dKInv(N, Mtot, dt,mathcalD, Dif,DInv, DphiInv, phi,dMeff , Meff_to_Phi_Inv, 4);
      //calc_dKInv(N, Mtot, dt,mathcalD,Dif, DInv, DphiInv, phi,dKc , K_to_Phi_Inv, 5);
    } else {
      printf("Enlarge Arrays...\n");
      // just enlarge arrays and double Nt
      enlarge_array_values(N, Mtot, Nt, &phi,phi);
      enlarge_array_values(N, 2*Mtot, Nt, &m,m);
      if (d<=Dswap)  enlarge_array_values(N, 2*Mtot, Nt, &K,K);
      enlarge_array_values(N, 2*Mtot, Nt, &Kr,Kr);
      if (d<=Dswap)  enlarge_array_values(N, Mtot, Nt, &Kc,Kc);
      enlarge_array_values(N, Mtot, Nt, &Krc,Krc);
      enlarge_array_values(N, Mtot, Nt, &Meff,Meff);
      enlarge_array_values(N, Mtot, Nt, &Meff_long,Meff_long);
      
      enlarge_array_moments(N, Mtot, Nt, &dphi, dphi, phi);
      enlarge_array_moments(N, 2*Mtot, Nt, &dm,dm, m);
      if (d<=Dswap)  enlarge_array_moments(N, 2*Mtot, Nt, &dK,dK, K);
      enlarge_array_moments(N, 2*Mtot, Nt, &dKr,dKr, Kr);
      if (d<=Dswap)  enlarge_array_moments(N, Mtot, Nt, &dKc,dKc, Kc);
      enlarge_array_moments(N, Mtot, Nt, &dKrc,dKrc, Krc);
      enlarge_array_moments(N, Mtot, Nt, &dMeff,dMeff, Meff);
      enlarge_array_moments(N, Mtot, Nt, &dMeff_long,dMeff_long, Meff_long);
      
      Nt *= 2;
    }
    
  //printf("%.15e %.15e %.15e %.15e\n",M_to_K_Inv[6][M+M*2*Mtot],M_to_K_Inv[6][M+Mtot+(M+Mtot)*2*Mtot] ,K_to_Meff_Inv[6][M+M*Mtot],Meff_to_Phi_Inv[6][M+M*Mtot]);
    
  } // the main decimation loop: finished simulation
  
#endif //NONERGODICITY
  
}

/**********************************************************************************************************************************************************/

/*IO functions
 */
void read_c_phi_modes(int Nt, int N, int M, double q0, double dq, double ** c, double ** phi, double ** dphi) {
  	double **sbest; double **cbest; double **s;
	double *qbest;
	double dummy_d;
	int Mtot = 2*M+1;
	int MIn = 5;
	int Mtot_read = 2*MIn+1;
	int Nbest = 1024;
	double *cf2;
	s = dmatrix(0,N-1,0, Mtot*Mtot-1);
	qbest = dvector(0,Nbest-1);
	sbest = dmatrix(0,Nbest-1,0, Mtot_read*Mtot_read-1);
	cbest = dmatrix(0,Nbest-1,0, Mtot_read*Mtot_read-1);
	
	FILE * in;
	in = fopen("c_modes.dat","r");
	if (in == NULL) {
		perror("Failed: read_c_modes ");
		exit (EXIT_FAILURE);
	}
	fscanf(in, "%*[^\n]\n", NULL);
	fscanf(in, "%*[^\n]\n", NULL);
	fscanf(in, "%*[^\n]\n", NULL);
	for (int q=0; q<Nbest; q++) {
		fscanf(in,"%lf ", &dummy_d);
		for (int mu=0; mu<Mtot_read; mu++) {
		  for (int nu=0; nu<Mtot_read; nu++) {
		    fscanf(in,"%lf ",&dummy_d);
		    cbest[q][Mtot_read*mu+nu] = dummy_d;
		  }
		}
		fscanf(in,"\n");
	}
	fclose(in);

	in = fopen("s_modes.dat","r");
	if (in == NULL) {
		perror("Failed: read_s_modes ");
		exit (EXIT_FAILURE);
	}
	fscanf(in, "%*[^\n]\n", NULL);
	fscanf(in, "%*[^\n]\n", NULL);
	fscanf(in, "%*[^\n]\n", NULL);
	for (int q=0; q<Nbest; q++) {
		fscanf(in,"%lf ", &dummy_d);
		qbest[q] = dummy_d;
		for (int mu=0; mu<Mtot_read; mu++) {
		  for (int nu=0; nu<Mtot_read; nu++) {
		    fscanf(in,"%lf ",&dummy_d);
		    sbest[q][Mtot_read*mu+nu] = dummy_d;
		  }
		}
		fscanf(in,"\n");
	}
	fclose(in);
	
  
	// find nhg and the values of s and c
	for (int iq=0; iq<N; iq++){
		double xp = q0+iq*dq;
		for (int j = 1; j<Nbest; j++) {
			if (qbest[j]>xp && qbest[j-1]<=xp) {
				double t = (xp-qbest[j-1])/(qbest[j]-qbest[j-1]);
				for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
				  int nu = mu;
#else
				  for (int nu=0; nu<Mtot; nu++) {
#endif
				    int mup = mu+MIn-M;
				    int nup = nu+MIn-M;
				    s[iq][Mtot*mu+nu] = sbest[j-1][Mtot_read*mup+nup] + t*(sbest[j][Mtot_read*mup+nup]-sbest[j-1][Mtot_read*mup+nup]);
				    c[iq][Mtot*mu+nu] = cbest[j-1][Mtot_read*mup+nup] + t*(cbest[j][Mtot_read*mup+nup]-cbest[j-1][Mtot_read*mup+nup]);
				    
#ifndef DIAGONAL_APPROXIMATION1
				  }
#endif
				}
				break;
			}
		}
	}

	// initialize phi
for (int q=0; q<N; q++) {
	for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
	  int nu = mu;
#else
	  for (int nu=0; nu<Mtot; nu++) {
#endif
		for (int t=0; t<Nt/2; t++) {

				phi[q][t*Mtot*Mtot+mu*Mtot+nu] = s[q][mu*Mtot+nu];

		}
#ifndef DIAGONAL_APPROXIMATION1
	  }
#endif
	}
    }
	
	// initialize moments of phi
	for (int q=0; q<N; q++) {
		for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
		  int nu = mu;
#else
		  for (int nu=0; nu<Mtot; nu++) {
#endif
			for (int t=1; t<Nt/2; t++) {
				dphi[q][t*Mtot*Mtot+mu*Mtot+nu] = 0.5*(phi[q][t*Mtot*Mtot+mu*Mtot+nu]+phi[q][(t-1)*Mtot*Mtot+mu*Mtot+nu]);
			}
#ifndef DIAGONAL_APPROXIMATION1
		  }
#endif
		}
	}
	
}

/**********************************************************************************************************************************************************/

void read_n_modes(int M, double L, double * n, double * v,const char* density_file) {
  FILE * in;
  in = fopen(density_file,"r");

  if (in == NULL) {
    perror("Failed: read_n_modes ");
    exit (EXIT_FAILURE);
  }
  
  // use all the information we have -> n and v can be calculated for more modes than c,phi and m
  double dummy_d1, dummy_d2;
  int Mtot = 2*M+1;
  int Mtot3 = 6*M+1;
  double Qpre = 2*PI/L; 
  
  for (int mu=0; mu<Mtot3; mu++) {
    n[mu] = 0.0;
    v[mu] = 0.0;
  }
  
  int N = (int) (L*1000+1.5);
 double dz = 0.001;
  
  for (int z=0; z<N; z++) {
    fscanf(in,"%lf %lf\n", &dummy_d1, &dummy_d2);
    //printf("%lf %lf\n",dummy_d1,dummy_d2);
    double factor = dz;
    if (z==0 || z==N-1) factor *= 0.5;
    for (int mu=0; mu<Mtot3; mu++) {
      n[mu] += factor*dummy_d2*cos(Qpre*((double) (mu-3*M))*(dummy_d1-(L/2.0)));
      v[mu] += factor/dummy_d2*cos(Qpre*((double) (mu-3*M))*(dummy_d1-(L/2.0)));
    }
  }
  printf("n0 %.10f\n",n[3*M]);

}

/**********************************************************************************************************************************************************/


void initInv(int N, int M, double n0, double q0, double dq, double D0, double *Q, double * n, double ** phi,  double **Dif, double **DInv, double **DphiInv, double *mathcalD, double *mathcalDInv){  
  int Mtot = 2*M+1;
  
  double * tmp_phiInv = dvector(0,Mtot*Mtot-1);
  
  
#ifdef DIAGONAL_APPROXIMATION1
  for (int q=0; q<N; q++) {
    double qval = q0 + q*dq;
    for (int mu =0; mu< Mtot; mu++) {
      DphiInv[q][mu*Mtot+mu] = D0*( qval*qval + Q[mu+2*M]*Q[mu+2*M] ) / phi[q][mu*Mtot + mu];
      Dif[q][mu*Mtot+mu] = D0*( qval*qval + Q[mu+2*M]*Q[mu+2*M] ) ;
      DInv[q][mu*Mtot+mu] = 1.0/ (D0*( qval*qval + Q[mu+2*M]*Q[mu+2*M] ) );
      mathcalD[mu*2*Mtot+mu] = 1.0;
      mathcalD[(mu+Mtot)*2*Mtot+mu+Mtot] = 1.0;
    }

  }
  
#else
  // calculate inverse (sum rules not applicable due to M cutoff
  for (int q=0; q<N; q++) {
    double qval = q0 + q*dq;
    for (int mu =0; mu< Mtot; mu++) {
      for (int nu=0; nu<Mtot; nu++) {
        Dif[q][mu*Mtot+nu] = D0*n[mu-nu+3*M]*( qval*qval + Q[mu+2*M]*Q[nu+2*M] )/n0;
      }
    }
    
    calc_matrixInv(&phi[q][0],Mtot,tmp_phiInv);
    calc_matrixInv(&Dif[q][0],Mtot,&DInv[q][0]);
    for (int mu =0; mu< Mtot; mu++) {
      for (int nu=0; nu<Mtot; nu++) {
        DphiInv[q][mu*Mtot+nu] = 0.0;
        for (int sigma=0; sigma<Mtot; sigma++) {
          DphiInv[q][mu*Mtot+nu] += D0*n[mu-sigma+3*M]/n0*( qval*qval + Q[mu+2*M]*Q[sigma+2*M] ) * tmp_phiInv[sigma*Mtot + nu];
        }
      }
    }
  }
  
   for (int a=0; a<2; a++) {
    for (int mu=0; mu<Mtot; mu++) {
      for (int b=0; b<2; b++) {
        for (int nu=0; nu<Mtot; nu++) {
          if (a==b) mathcalD[nu+b*Mtot+2*Mtot*mu+2*Mtot*Mtot*a] =  n[mu -nu + 3*M]/n[3*M];
          else mathcalD[nu+b*Mtot+2*Mtot*mu+2*Mtot*Mtot*a] =  0.0;
        }
      }
    }
  }
 calc_matrixInv(mathcalD,2*Mtot,mathcalDInv);
  
#endif


  free_dvector(tmp_phiInv,0,Mtot*Mtot-1);

}

/**********************************************************************************************************************************************************/

void init_values_moments(int Nt, int N, int M, double n0, double  L, double q0, double dq, double dt,double D0, double * Q,double * v, double * n, double ** c,  double * mathcalD,  double ** DInv, double ** YS, double ** YM, double ** YSQ, double ** YMQ, double ** phi, double **m ,double **K ,double **Kr ,double **Kc,double **Krc , double **Meff , double **Meff_long , double **dm ,double **dK , double **dKr , double **dKc  , double **dKrc , double **dMeff , double **dMeff_long){

  int Mtot = 2*M+1;
  
  printf("Start Initialization\n");

  // initialize m, Kc and K 
#ifndef DIAGONAL_APPROXIMATION1
  calc_mkernel_init(0, N, M, n0, L, Q,  c, phi,YS, YM, YSQ, YMQ);
#endif
  for (int q=0; q<N; q++) { 
    // calc_mkernel for one time step
    calc_mkernel(0, N, q,  M, n0,L, q0, dq, c, Q,n,YS,YM,YSQ,YMQ,phi, &m[q][0]);
  }
    
  for (int q=0; q<N; q++) { 
    // since phi is constant for all t it can be transfered to t<Nt/2
    for (int t=1; t<Nt/2; t++) {
      for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
        int nu=mu;
#else
        for (int nu=0; nu<2*Mtot; nu++) {
#endif
          int indext = nu + mu*2*Mtot + t*4*Mtot*Mtot;
          int index = nu +  mu*2*Mtot;
          m[q][indext] = m[q][index];
#ifndef DIAGONAL_APPROXIMATION1
        }
#endif
      }
    }
    
    // K is equal to the negative memory kernel, Kr is diffusive
    for (int t=0; t<Nt/2; t++) {
        for (int mu=0; mu<2*Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
          int nu=mu;
#else
          for (int nu=0; nu<2*Mtot; nu++) {
#endif
            int indext = nu + mu*2*Mtot + t*4*Mtot*Mtot;
            K[q][indext] = 0.0;
            Kr[q][indext] = dt*t*mathcalD[nu+mu*2*Mtot];
#ifdef DIAGONAL_APPROXIMATION1
            int sigma = mu;
            int kappa = mu;
            K[q][indext] -= mathcalD[sigma+mu*2*Mtot]*m[q][kappa+sigma*2*Mtot + 4*Mtot*Mtot*t]*mathcalD[nu+kappa*2*Mtot];
#else
            for (int sigma=0; sigma<2*Mtot; sigma++) {
                K[q][indext] -= m[q][sigma+mu*2*Mtot + 4*Mtot*Mtot*t]*mathcalD[nu+sigma*2*Mtot];
            }
#endif
          //if (q==5 &&mu==2) printf("K: %f\n",K[q][mu+mu*2*Mtot]);
	      
#ifndef DIAGONAL_APPROXIMATION1
        }
#endif
      }
    }
    
    // K is the constracted version of Kc, same for Krc
    for (int t=0; t<Nt/2; t++) {
      contract_K(q, M, q0, dq,L,n0,Q, v,c,&K[q][4*Mtot*Mtot*t], &Kc[q][Mtot*Mtot*t]);
      contract_K(q, M, q0, dq,L,n0,Q, v,c,&Kr[q][4*Mtot*Mtot*t], &Krc[q][Mtot*Mtot*t]);
    }
    
    // init Meff is equal to Kc
    for (int t=0; t<Nt/2; t++) {
      for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
        int nu=mu;
#else
        for (int nu=0; nu<Mtot; nu++) {
#endif
        
          int indext = nu+mu*Mtot+ t*Mtot*Mtot;
          double qval = q0+q*dq;

          Meff[q][indext] = 0.0;
          Meff_long[q][indext] = 0.0;
          //if (q==6 &&mu==M) printf("DInv: %f\n",DInv[q][mu+mu*Mtot]);
#ifdef DIAGONAL_APPROXIMATION1
          int sigma = mu;
          int kappa = mu;
          Meff[q][indext] -= DInv[q][sigma+mu*Mtot]*Kc[q][kappa+sigma*Mtot + Mtot*Mtot*t]*DInv[q][nu+kappa*Mtot];
          Meff_long[q][indext] -= DInv[q][sigma+mu*Mtot]*Kc[q][kappa+sigma*Mtot + Mtot*Mtot*t]*DInv[q][nu+kappa*Mtot];
#else
          for (int sigma=0; sigma<Mtot; sigma++) {
            for (int kappa=0; kappa<Mtot; kappa++) {
              Meff[q][indext] -= DInv[q][sigma+mu*Mtot]*Kc[q][kappa+sigma*Mtot + Mtot*Mtot*t]*DInv[q][nu+kappa*Mtot];
              Meff_long[q][indext] -= DInv[q][sigma+mu*Mtot]*Kc[q][kappa+sigma*Mtot + Mtot*Mtot*t]*DInv[q][nu+kappa*Mtot];
            }
          }
          
#endif
        //  if (q==5 &&mu==2) printf("Kc: %f\n",Kc[q][mu+mu*Mtot]);
        //if (q==5 &&mu==2) printf("Meff: %f\n",Meff[q][mu+mu*Mtot]);
        
#ifndef DIAGONAL_APPROXIMATION1
        }
#endif
  
      }
    }
        
    // init moments of m and K
    for (int t=1; t<Nt/2; t++) {
      for (int a=0; a<2; a++) {
        for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
          int b=a;
          int nu=mu;
#else
          for (int b=0; b<2; b++) {
            for (int nu=0; nu<Mtot; nu++) {
#endif
              int indext = nu + b*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot + t*4*Mtot*Mtot;
              int indext1 = nu + b*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot + (t-1)*4*Mtot*Mtot;
              dm[q][indext] = 0.5*(m[q][indext]+ m[q][indext1]);
              dK[q][indext] = 0.5*(K[q][indext]+ K[q][indext1]);
              dKr[q][indext] = 0.5*(Kr[q][indext]+ Kr[q][indext1]);
#ifndef DIAGONAL_APPROXIMATION1
            }
          }
#endif
        }
      }
    }
    
    // init moments of Kc,Meff
    for (int t=1; t<Nt/2; t++) {
      for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
        int nu=mu;
#else
        for (int nu=0; nu<Mtot; nu++) {
#endif
          int indext = nu + mu*Mtot + t*Mtot*Mtot;
          int indext1 = nu + mu*Mtot + (t-1)*Mtot*Mtot;
          dKc[q][indext] = 0.5*(Kc[q][indext]+ Kc[q][indext1]);
          dKrc[q][indext] = 0.5*(Krc[q][indext]+ Krc[q][indext1]);
          dMeff[q][indext] = 0.5*(Meff[q][indext]+ Meff[q][indext1]);
          dMeff_long[q][indext] = 0.5*(Meff_long[q][indext]+ Meff_long[q][indext1]);
#ifndef DIAGONAL_APPROXIMATION1
        }
#endif
      }
    }
  }
  
    printf("Done Initialization\n");

}

/**********************************************************************************************************************************************************/

void contract_K(int q, int  M,double q0,double dq,double L, double n0, double * Q,double * v,double ** c, double * K, double * Kc){
  
  int Mtot = 2*M+1;
  
  // contract
  for (int mu=0; mu<Mtot; mu++) {
#ifdef DIAGONAL_APPROXIMATION1
      int nu=mu;
#else
	for (int nu=0; nu<Mtot; nu++) {
#endif
      Kc[nu+mu*Mtot] = 0.0;
      for (int a=0; a<2; a++) {
        double vala = q0+q*dq;
        if (a==1) vala = Q[mu+2*M];
#ifdef DIAGONAL_APPROXIMATION1
        int b=a;
#else
        for (int b=0; b<2; b++) {
#endif
          double valb = q0+q*dq;
          if (b==1) valb = Q[nu+2*M];
          Kc[nu+mu*Mtot] += vala*valb*K[nu + b*Mtot + mu*2*Mtot + 2*a*Mtot*Mtot];
        }
      
#ifndef DIAGONAL_APPROXIMATION1
      }
    }
#endif
  }
  
}

/**********************************************************************************************************************************************************/

void decimation_moments(int Nt, int N, int size,double ** dinput,double ** input){
  double * tmp;
  tmp = dvector(0,Nt-1);
  
  for (int q=0; q<N; q++) {
    for (int a=0; a<size; a++) {
#ifdef DIAGONAL_APPROXIMATION1
      int b=a;
#else
	for (int b=0; b<size; b++) {
#endif
	// perform decimation 1
	for (int t=1; t<Nt/4; t++) {
	  tmp[t] = 0.5*(dinput[q][(2*t-1)*size*size+a*size+b]+dinput[q][(2*t)*size*size+a*size+b]);
	}
	
	// perform decimation 2
	for (int t=Nt/4; t<Nt/2; t++) {
	  tmp[t] = 1.0/6.0*(input[q][(2*t-2)*size*size+a*size+b]+4.0*input[q][(2*t-1)*size*size+a*size+b]+input[q][(2*t)*size*size+a*size+b]);
	}
	
	for (int t=1; t<Nt/2; t++) {
	  dinput[q][t*size*size+a*size+b] = tmp[t];
	}
#ifndef DIAGONAL_APPROXIMATION1
      }
#endif
    }
  }
  free_dvector(tmp,0,Nt-1);
}

/**********************************************************************************************************************************************************/

void decimation_values(int Nt, int N, int size,double ** input){
  
  for (int q=0; q<N; q++) {
    for (int a=0; a<size; a++) {
#ifdef DIAGONAL_APPROXIMATION1
      int b=a;
#else
	for (int b=0; b<size; b++) {
#endif
	// perform decimation 1
	for (int t=0; t<Nt/2; t++) {
	  input[q][t*size*size+a*size+b] = input[q][2*t*size*size+a*size+b];
	}
#ifndef DIAGONAL_APPROXIMATION1
      }
#endif
    }
  }
}

/**********************************************************************************************************************************************************/

void calc_dKInv(int N, int size, double dt, double *mathcalD, double **Dif,double **DInv, double **DphiInv, double ** phi,  double ** dK , double ** dKInv, int integro){
  
  double ** dKtmp;
  double ** dKInvtmp;
  
  dKtmp = dmatrix(0, size-1,0,size-1);
  dKInvtmp = dmatrix(0, size-1,0,size-1);
  
  // calculate inverse of phi
  double * phiInv;
  phiInv = dvector(0,size*size-1);

  
  for (int q=0; q<N; q++) {
    if (integro == 5) calc_matrixInv(&phi[q][0],size,phiInv);
    for (int a=0; a<size; a++) {
      for (int b=0; b<size; b++) {
#ifdef DIAGONAL_APPROXIMATION1
        if (a==b) {
#endif
        dKtmp[a][b] = 0.0;
        if (integro == 0) { // integro = 0 -> M_to_K
#ifdef DIAGONAL_APPROXIMATION1
          int c = a;
          dKtmp[a][b] += dt/2.0*dK[q][a*size+b+size*size];
#else
          dKtmp[a][b] += dt/2.0*dK[q][a*size+b+size*size];
#endif
        } else if (integro == 1) { // integro = 1 -> M_to_Kr
#ifdef DIAGONAL_APPROXIMATION1
          int c = a;
          dKtmp[a][b] += beta0*dt*dK[q][a*size+b+size*size];
#else
          dKtmp[a][b] += beta0*dt*dK[q][a*size+b+size*size];
#endif
        } else if (integro == 2) {  // integro = 2 -> K_to_Meff
#ifdef DIAGONAL_APPROXIMATION1
          int c = a;
          dKtmp[a][b] += dt/2.0*dK[q][a*size+c+size*size]*DInv[q][c*size+b];
          //if (q==6&&a==5) printf("dKInv %.10f %.10f %.10f \n",dKtmp[a][b],dK[q][a*size+c+size*size],DInv[q][c*size+b]);
#else
          for (int c=0; c<size; c++) {
            dKtmp[a][b] += dt/2.0*dK[q][a*size+c+size*size]*DInv[q][c*size+b];
          }
#endif
        } else if (integro == 3) {  // integro = 3 -> Kr_to_Meff
          dKtmp[a][b] += beta0*dt*dK[q][b+a*size+size*size];
        } else if (integro == 4) { // integro = 4 -> Meff_to_Phi
#ifdef DIAGONAL_APPROXIMATION1
          int c = a;
          dKtmp[a][b] += beta0*dt*DphiInv[q][a*size+b] + beta0*dt*Dif[q][a*size+c]*dK[q][c*size+b+size*size];
#else
          dKtmp[a][b] += beta0*dt*DphiInv[q][a*size+b];
          for (int c=0; c<size; c++) {
            dKtmp[a][b] +=  beta0*dt*Dif[q][a*size+c]*dK[q][c*size+b+size*size];
          }
#endif
        } else if (integro == 5) { // integro = 5 -> K_to_Phi
#ifdef DIAGONAL_APPROXIMATION1
          int c = a;
          dKtmp[a][b] += beta0*dt*DphiInv[q][a*size+b] + beta0*dt*dt/2.0*dK[q][c+a*size+size*size]*phiInv[b+c*size];
#else
          dKtmp[a][b] += beta0*dt*DphiInv[q][a*size+b];
          for (int c=0; c<size; c++) {
            dKtmp[a][b] +=  beta0*dt*dt/2.0*dK[q][c+a*size+size*size]*phiInv[b+c*size];
          }
#endif
        } else {
          printf("Integro switch does not exist! %d \n",integro);
          exit;
        }
        if (a==b && integro != 3) dKtmp[a][b] += 1.0;
#ifdef DIAGONAL_APPROXIMATION1
        } else {
          dKtmp[a][b] = 0.0;
        }
#endif

      }
    }
    invcmp2(dKtmp, size, dKInvtmp);
    for (int a=0; a<size; a++) {
      for (int b=0; b<size; b++) {
        dKInv[q][b + a*size] = dKInvtmp[a][b];
      }
    }
  }
  
  free_dmatrix(dKtmp,0, size-1,0,size-1);
  free_dmatrix(dKInvtmp,0, size-1,0,size-1);	
    free_dvector(phiInv,0,size*size-1);

}

/**********************************************************************************************************************************************************/

void calc_matrixInv(double *in, int size, double * out){
  double ** tmp;
  double ** Invtmp;
  
  tmp = dmatrix(0, size-1,0,size-1);
  Invtmp = dmatrix(0, size-1,0,size-1);
  
  for (int a=0; a<size; a++) {
    for (int b=0; b<size; b++) {
      tmp[a][b] = in[a*size+b];
    }
  }
  
  invcmp2(tmp, size, Invtmp);
  
  for (int a=0; a<size; a++) {
    for (int b=0; b<size; b++) {
      out[b + a*size] = Invtmp[a][b];
    }
  }

  
  free_dmatrix(tmp,0, size-1,0,size-1);
  free_dmatrix(Invtmp,0, size-1,0,size-1);	
}

/**********************************************************************************************************************************************************/

void enlarge_array_values(int N, int size, int Nt, double *** array , double ** array_value){
  double ** array_loc = dmatrix(0,N-1,0, (2*Nt+1)*size*size-1);
  for (int q=0; q<N; q++) {
    for (int t=0; t<Nt; t++) {
      for (int mu=0; mu<size; mu++) {
	for (int nu=0; nu<size; nu++) {
	  array_loc[q][nu+mu*size+size*size*t] = array_value[q][nu+mu*size+size*size*t];
	}
      }
    }
  }
  double ** save = *array;
  *array = array_loc;
  free_dmatrix(save,0,N-1,0, (Nt+1)*size*size-1);
}

/**********************************************************************************************************************************************************/

void enlarge_array_moments(int N, int size, int Nt, double *** array, double ** array_moment, double ** array_value){
  double ** array_loc = dmatrix(0,N-1,0, (Nt+2)*size*size-1);
  for (int q=0; q<N; q++) {
    for (int t=1; t<Nt/2; t++) {
      for (int mu=0; mu<size; mu++) {
	for (int nu=0; nu<size; nu++) {
	  array_loc[q][nu+mu*size+size*size*t] = array_moment[q][nu+mu*size+size*size*t];
	}
      }
    }
    for (int t=Nt/2; t<Nt; t++) {
      for (int mu=0; mu<size; mu++) {
	for (int nu=0; nu<size; nu++) {
	  array_loc[q][nu+mu*size+size*size*t] = 1.0/2.0*(array_value[q][(t-1)*size*size+mu*size+nu]+array_value[q][t*size*size+mu*size+nu]);
	}
      }
    }
    
  }
  double ** save = *array;
  *array = array_loc;
  free_dmatrix(save,0,N-1,0, (Nt/2+2)*size*size-1);
}
