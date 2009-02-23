#include <stdlib.h>
#include <math.h>

/* Derived from Ramos-Onsins & Rozas's mlcoalsim */


/*  Calculates Equation 23 from Ewens 
 *
 *      N           - number of samples
 *      i           - ranges between 1 and total number of haplotypes
 *      theta       - in this case is the average pairwise number of 
 *                    subsitutions (pi)
 *      qew_        - 
 *
 *  Returns a double
 *  Also updates the qew_ array in progress
 */
double FunEq23Ewens(int N, int i, double theta, double *qew_)
{                  
    long int acceso; 
    int jj;
    double ValorN;  /* log del numerador */
    double ValorD;  /* log del denominador */

    acceso = (long int)(N-1) * (long int)N + (long int)i - (long int)1;    
    ValorN = 0.0;
    ValorD = 0.0;        

    if (qew_[acceso] < 0.0) {   
        if (i==1) {
            /* calculo de qj,1   (i = 1)   Antigua equacion 19  */
            if (N > 2) {             
                for (jj=2; jj<N; jj++)
                    ValorN = ValorN + log((double)jj);  
            }
            ValorN = ValorN + log(theta);
            for (jj=0; jj<N; jj++)
                ValorD  = ValorD + log((double)theta + (double)jj);      
            qew_[acceso] = exp((double)(ValorN - ValorD)); 
        }    
        if (i==N) {          
            /* calculo de qj,j   (n = i)   antigua equacion 20 */
            ValorN = log((double)theta) * (double)N;
            for (jj=0; jj<N; jj++)     
                ValorD  = ValorD + log((double)theta + (double)jj);
            qew_[acceso] = exp((double)(ValorN - ValorD));
	      }
	      if (i>1 && i<N) {    
            /*  recursividad  */
            qew_[acceso] = FunEq23Ewens(N-1,i,  theta,qew_) * ((double)(N-1)/(theta + (double)N-1.0))
                         + FunEq23Ewens(N-1,i-1,theta,qew_) *         (theta/(theta + (double)N-1.0));
        }    
    }  

    return qew_[acceso];
}

/*  Calculates Fu's Fs
 *
 *      Nsample     - number of samples
 *      pi          - average number of pairwise substitutions
 *      NumAlelos   - number of haplotypes
 *
 *  Returns a double
 */
double Fs(int Nsample, double pi, int NumAlelos)
{
    /* Rozas program */
	
    double    SumaP,
              RestaP,
              ValorFs,
              est_var;
    double    *qew;
    long int  i;          /* iterator */       

    if (pi == 0.0 || Nsample < 2) 
        return(-10000);	

    est_var = pi;

    qew  = (double *)malloc((long int)Nsample*(long int)Nsample*sizeof(double));
    
    for (i=0; i<(long int)Nsample*(long int)Nsample ;i++)
    	  qew[i] = -1.0;
            
    SumaP=RestaP=0.0;

    for (i=1; i<NumAlelos; i++) {
        /* calculo q(n,aleloI)   ecuacion 21 (recurrente con eq. 19 y 20) */
        SumaP += FunEq23Ewens(Nsample, i, est_var, qew);
    }

    if (SumaP > 1.-1E-37) {
    	  for (i=NumAlelos; i<=Nsample; i++)
            RestaP += FunEq23Ewens(Nsample, i, est_var, qew);	 	
        if (RestaP < 1E-37)
            return -10000;

        ValorFs = log((double)RestaP) - log((double)1.0-RestaP);
    } else {
        if (SumaP < 1E-37)
            return +10000;
        else
            ValorFs = log((double)1.0 - (double)SumaP) - log((double)SumaP);
    }

    if (fabs(ValorFs) < 1.0E-15)
        ValorFs = 0.0;	    
    
    free(qew);

    return ValorFs;
}


