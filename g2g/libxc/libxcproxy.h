#ifndef LIBXCPROXY_H
#define LIBXCPROXY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <xc.h>
#include <vector>
#include "../scalar_vector_types.h"
#include "print_utils.h"
#include "../fix_compile.h"
//#include <cuda_runtime.h>
#include "../timer.h"

extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg);

template <class T, int width>
class LibxcProxy
{
private:

    // The libxc components
    xc_func_type funcForExchange;
    xc_func_type funcForCorrelation;

    // Functional ids
    int funcIdForExchange;
    int funcIdForCorrelation;
    int nspin;

    // Is inited
    bool inited;

    void printFunctionalInformation (xc_func_type* func);

public:
    LibxcProxy ();
    LibxcProxy (int exchangeFunctionId, int correlationFuncionalId, int nspin);
    ~LibxcProxy ();

    void doGGA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);

    void doGGA (T* dens,
                const int number_of_points,
		const G2G::vec_type<T,width>* grad,
                const G2G::vec_type<T,width>* hess1,
                const G2G::vec_type<T,width>* hess2,
                T* ex,
                T* ec,
                T* y2a);

    void doGGA (T dens,
                T sigma,
                T* v2rho2,
                T v2rhosigma,
                T v2sigma2);

    void doGGA (T* dens,
		T* sigma,
                const int number_of_points,
		T* v2rho2,
		T* v2rhosigma,
		T* v2sigma2);

    void doGGA (T* dens,
                const int number_of_points,
		const T* contracted_grad,
		const G2G::vec_type<T,width>* grad,
                const G2G::vec_type<T,width>* hess1,
                const G2G::vec_type<T,width>* hess2,
                T* ex,
                T* ec,
                T* y2a);

    void coefLR(double* rho,double* sigma,
                double red,double cruz,double* lrCoef);

    // Open LR
    void coefLR(double* rho,double* tra,double* lrCoef);

    void coefZv(double* rho,double* sgm, // inputs
                // first derivatives
                double* vrho,double* vsigma, 
                // second derivatives
                double* v2rho2, double* v2rhosigma,
                double* v2sigma2, 
                // third derivatives
                double* v3rho3, double* v3rho2sigma,
                double* v3rhosigma2, double* v3sigma3);

    void doLDA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);

    void init (int exId, int xcId, int nspin);
    void closeProxy ();
    void printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId);
};

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy()
{
    funcIdForExchange = 0;
    funcIdForCorrelation = 0;
    nspin = 0;
    inited = false;
}

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy (int exchangeFunctionalId, int correlationFuncionalId, int nSpin)
{
//    printf("LibxcProxy::LibxcProxy (%u, %u, %u) \n", exchangeFunctionalId, correlationFuncionalId, nSpin);
/*
    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFuncionalId;
    nspin = nSpin;

    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
    }
*/
    init (exchangeFunctionalId, correlationFuncionalId, nSpin);
}

template <class T, int width>
void LibxcProxy<T, width>::init (int exchangeFunctionalId, int correlationFunctionalId, int nSpin)
{
//    printf("LibxcProxy::init(%u, %u, %u)\n", exchangeFunctionalId, correlationFunctionalId, nSpin);
    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFunctionalId;
    nspin = nSpin;

    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
    }

    inited = true;
}

template <class T, int width>
LibxcProxy <T, width>::~LibxcProxy ()
{
    //xc_func_end (&funcForExchange);
    //xc_func_end (&funcForCorrelation);
//    printf("LibxcProxy::~LibxcProxy()\n");
    closeProxy ();
}

template <class T, int width>
void LibxcProxy <T, width>::closeProxy ()
{
//    printf("LibxcProxy::closeProxy()\n");
    if (inited) {
	xc_func_end (&funcForExchange);
        xc_func_end (&funcForCorrelation);
	inited = false;
    }
}

template <class T, int width>
void LibxcProxy<T, width>::printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId) 
{
    if (!inited) 
    {
	init (exchangeFunctionalId, correlationFunctionalId, 1);
    }

    printFunctionalInformation (&funcForExchange);
    printFunctionalInformation (&funcForCorrelation);

    if (inited) 
    {
	closeProxy ();
    }
}

template <class T, int width>
void LibxcProxy<T, width>::printFunctionalInformation (xc_func_type* func)
{
    printf("The functional '%s' is ", func->info->name);
    switch (func->info->kind) {
	case (XC_EXCHANGE):
	    printf("an exchange functional");
	break;
	case (XC_CORRELATION):
	    printf("a correlation functional");
	break;
	case (XC_EXCHANGE_CORRELATION):
	    printf("an exchange-correlation functional");
	break;
	case (XC_KINETIC):
	    printf("a kinetic energy functional");
	break;
	default:
	    printf("of unknown kind");
	break;
    }

    printf(", it belongs to the '", func->info->name);
    switch (func->info->family) {
	case (XC_FAMILY_LDA):
	    printf("LDA");
        break;
	case (XC_FAMILY_GGA):
	    printf("GGA");
        break;
	case (XC_FAMILY_HYB_GGA):
	    printf("Hybrid GGA");
	break;
	case (XC_FAMILY_MGGA):
	    printf("MGGA");
        break;
	case (XC_FAMILY_HYB_MGGA):
	    printf("Hybrid MGGA");
        break;
	default:
	    printf("unknown");
        break;
    }
    printf("' family and is defined in the reference(s):\n");

    for (int ii = 0; func->info->refs[ii] != NULL; ii++) {
	printf ("[%d] %s\n", ii+1, func->info->refs[ii]->ref);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 1
// Calls the XC_GGA function from libxc for multiple points.
// dens: pointer for the density array
// grad:.
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
// y2a:
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T dens,
    const G2G::vec_type<T, width> &grad,
    const G2G::vec_type<T, width> &hess1,
    const G2G::vec_type<T, width> &hess2,
    T &ex, T &ec, T &y2a)
{
    //printf("LibxcProxy::doGGA cpu simple(...) \n");

    const double rho[1] = {dens};
    // Libxc needs the 'contracted gradient'
    double sigma[1] = {(grad.x * grad.x) + (grad.y * grad.y) + (grad.z * grad.z)};
    double exchange[1];
    double correlation[1];

    // The outputs for exchange
    double vrho [1];
    double vsigma [1];
    double v2rho [1];
    double v2rhosigma[1];
    double v2sigma [1];

    // The outputs for correlation
    double vrhoC [1];
    double vsigmaC [1];
    double v2rhoC [1];
    double v2rhosigmaC [1];
    double v2sigmaC [1];

    // The exchange values
    xc_gga (&funcForExchange, 1,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,NULL);

    // Now the correlation value.
    xc_gga (&funcForCorrelation, 1,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);

    ex = exchange[0];
    ec = correlation[0];

    // Merge the results for the derivatives.
    vrho[0] += vrhoC[0];
    vsigma[0] += vsigmaC[0];
    v2rho[0] += v2rhoC[0];
    v2rhosigma[0] += v2rhosigmaC[0];
    v2sigma[0] += v2sigmaC[0];

    // Now, compute y2a value.
    y2a = vrho[0] - (2 * sigma[0] * v2rhosigma[0]
            + 2 * (hess1.x + hess1.y + hess1.z) * vsigma[0]
            + 4 * v2sigma[0] * (grad.x * grad.x * hess1.x + grad.y * grad.y * hess1.y + grad.z * grad.z * hess1.z + 2 * grad.x * grad.y * hess2.x + 2 * grad.x * grad.z * hess2.y + 2 * grad.y * grad.z * hess2.z));

    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 2
// Calls the XC_GGA function from libxc for multiple points.
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// grad: gradient value in each point
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
// y2a:
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T* dens,
    const int number_of_points,
    const G2G::vec_type<T, width>* grad,
    const G2G::vec_type<T, width>* hess1,
    const G2G::vec_type<T, width>* hess2,
    T* ex,
    T* ec,
    T* y2a)
{
    //printf("LibxcProxy::doGGA cpu multiple (...) \n");

    int array_size = sizeof(double)*number_of_points;
    double* rho = (double*)malloc(array_size);
    for (int i=0; i<number_of_points; i++) {
	rho[i] = (double)dens[i];
    }

    // Libxc needs the 'contracted gradient'
    double* sigma = (double*)malloc(array_size);
    for (int i=0; i< number_of_points; i++) {
	sigma[i] = (double)((grad[i].x * grad[i].x) + (grad[i].y * grad[i].y) + (grad[i].z * grad[i].z));
    }
    double* exchange = (double*)malloc(array_size);
    double* correlation = (double*)malloc(array_size);

    // The outputs for exchange
    double* vrho = (double*)malloc(array_size);
    double* vsigma = (double*)malloc(array_size);
    double* v2rho = (double*)malloc(array_size);
    double* v2rhosigma = (double*)malloc(array_size);
    double* v2sigma = (double*)malloc(array_size);

    // The outputs for correlation
    double* vrhoC = (double*)malloc(array_size);
    double* vsigmaC = (double*)malloc(array_size);
    double* v2rhoC = (double*)malloc(array_size);
    double* v2rhosigmaC = (double*)malloc(array_size);
    double* v2sigmaC = (double*)malloc(array_size);

    // Exchange values
    xc_gga (&funcForExchange, number_of_points,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);

    // Now the correlation value.
    xc_gga (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);

    for (int i=0; i<number_of_points; i++) {
	ex[i] = exchange[i];
        ec[i] = correlation[i];
        // Merge the results for the derivatives.
        vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];

	// Now, compute y2a value.
        y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
	        + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
    		+ 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x +
				    grad[i].y * grad[i].y * hess1[i].y +
				    grad[i].z * grad[i].z * hess1[i].z +
				    2 * grad[i].x * grad[i].y * hess2[i].x +
				    2 * grad[i].x * grad[i].z * hess2[i].y +
				    2 * grad[i].y * grad[i].z * hess2[i].z));
    }

    // Free memory.
    free(rho);
    free(sigma);
    free(exchange);
    free(correlation);

    // The outputs for exchange
    free(vrho);
    free(vsigma);
    free(v2rho);
    free(v2rhosigma);
    free(v2sigma);

    // The outputs for correlation
    free(vrhoC);
    free(vsigmaC);
    free(v2rhoC);
    free(v2rhosigmaC);
    free(v2sigmaC);

    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 3
// Calls the XC_GGA_FXC function from libxc for multiple points.
// rho: pointer for the density array.
// sigma: contracted gradient array.
// number_of_points: the size of all the input arrays.
// v2rho2: second partial derivative of the energy per unit volume in terms of the density.
// v2rhosigma: second partial derivative of the energy per unit volume in terms of the density and sigma.
// v2sigma2: second partial derivative of the energy per unit volume in terms of sigma.
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA (T* rho,
    T* sigma,
    const int number_of_points,
    T* v2rho2,
    T* v2rhosigma,
    T* v2sigma2)
{
    int array_size = sizeof(double)*number_of_points;

    // The outputs for exchange
    double* v2rho2X = (double*)malloc(array_size);
    double* v2rhosigmaX = (double*)malloc(array_size);
    double* v2sigma2X = (double*)malloc(array_size);

    // The outputs for correlation
    double* v2rho2C = (double*)malloc(array_size);
    double* v2rhosigmaC = (double*)malloc(array_size);
    double* v2sigma2C = (double*)malloc(array_size);

    // Exchange values
    xc_gga_fxc (&funcForExchange, number_of_points,
                rho,
                sigma,
                v2rho2X,
                v2rhosigmaX,
                v2sigma2X);

    // Now the correlation value.
    xc_gga_fxc (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                v2rho2C,
                v2rhosigmaC,
                v2sigma2C);

    for (int i=0; i<number_of_points; i++) {
        // Merge the results for the derivatives.
        v2rho2[i] = v2rho2X[i] + v2rho2C[i];
        v2rhosigma[i] = v2rhosigmaX[i] + v2rhosigmaC[i];
        v2sigma2[i] += v2sigma2X[i] + v2sigma2C[i];
    }

    // Free memory
    // The outputs for exchange
    free(v2rho2X);
    free(v2rhosigmaX);
    free(v2sigma2X);

    // The outputs for correlation
    free(v2rho2C);
    free(v2rhosigmaC);
    free(v2sigma2C);

    return;
}

template <class T, int width>
void LibxcProxy <T, width>::coefLR (double *rho,
                double* sgm,
                double red,
                double cruz,
                double* lrCoef)
{
   // The otputs for exchange
   double vrhoX[2], vsigmaX[3],v2rho2X[3],v2rhosigmaX[6],v2sigma2X[6];

   // The ouputs for correlation
   double vrhoC[2], vsigmaC[3],v2rho2C[3],v2rhosigmaC[6],v2sigma2C[6];

   // NOT Refence
   double exc = 0.0f;

   // convert to alfa and beta;
   double dens[2], sigma[3];
          dens[0] = dens[1] = *rho * 0.5f;
          sigma[0] = sigma[1] = sigma[2] = *sgm * 0.25f;

   // Exchange values
   xc_gga(&funcForExchange,1,dens,sigma,&exc,vrhoX,vsigmaX,
          v2rho2X,v2rhosigmaX,v2sigma2X,NULL,NULL,NULL,NULL,
          NULL,NULL,NULL,NULL,NULL);
   
   // Correlation values
   xc_gga(&funcForCorrelation,1,dens,sigma,&exc,vrhoC,vsigmaC,
          v2rho2C,v2rhosigmaC,v2sigma2C,NULL,NULL,NULL,NULL,
          NULL,NULL,NULL,NULL,NULL);

   // Results
   double term1, term2, term3;
   term1 = red * v2rho2X[0] + 2.0f * v2rhosigmaX[0] * cruz;
   term2 = red * v2rho2C[0] + 2.0f * v2rhosigmaC[0] * cruz;
         term2 += cruz * v2rhosigmaC[1];
   term3 = red * v2rho2C[1] + 2.0f * v2rhosigmaC[2] * cruz;
         term3 += cruz * v2rhosigmaC[1] ;
   lrCoef[0] = term1 + term2 + term3;

   term1 = red * v2rhosigmaX[0] * 2.0f + cruz * v2sigma2X[0] * 4.0f;
   term2 = red * (2.0f * v2rhosigmaC[0] + v2rhosigmaC[1]);
   term2 += cruz * (4.0f * v2sigma2C[0] + 2.0f * v2sigma2C[1]);
   term2 += cruz * (2.0f * v2sigma2C[1] + v2sigma2C[3]);
   term3 = red * (2.0f * v2rhosigmaC[3] + v2rhosigmaC[1]);
   term3 += cruz * (4.0f*v2sigma2C[2]+2.0f*v2sigma2C[1]+2.0f*v2sigma2C[4]+v2sigma2C[3]);
   lrCoef[1] = term1 + term2 + term3;
   
   term1 = 2.0f * vsigmaX[0];
   term2 = 2.0f * vsigmaC[0];
   term3 = vsigmaC[1];
   lrCoef[2] = term1 + term2 + term3;

   return;
}

// OPEN LR COEFF
template <class T, int width>
void LibxcProxy <T, width>::coefLR (double* rho,double* tra,double* lrCoef)
{
/*
INPUTS:
   rho[0] = alpha density
   rho[1] = alpha density gradient X
   rho[2] = alpha density gradient Y
   rho[3] = alpha density gradient Z
   rho[4] = beta density
   rho[5] = beta density gradient X
   rho[6] = beta density gradient Y
   rho[7] = beta density gradient Z

   tra[0] = alpha transition density
   tra[1] = alpha transition density gradient X
   tra[2] = alpha transition density gradient Y
   tra[3] = alpha transition density gradient Z
   tra[4] = beta transition density
   tra[5] = beta transition density gradient X
   tra[6] = beta transition density gradient Y
   tra[7] = beta transition density gradient Z
   tra[8] = alpha transition density TWO gradients 
   tra[9] = beta transition density TWO gradiensts

*/
   // The otputs for exchange
   double vrhoX[2], vsigmaX[3],v2rho2X[3],v2rhosigmaX[6],v2sigma2X[6];

   // The ouputs for correlation
   double vrhoC[2], vsigmaC[3],v2rho2C[3],v2rhosigmaC[6],v2sigma2C[6];

   // NOT Refence
   double exc = 0.0f;

   // Varibles to LIBXC
   double dens[2], sigma[3];
   dens[0]  = rho[0]; dens[1] = rho[4];
   sigma[0] = rho[1] * rho[1] + rho[2] * rho[2] + rho[3] * rho[3]; // sigma_aa
   sigma[1] = rho[1] * rho[5] + rho[2] * rho[6] + rho[3] * rho[7]; // sigma_ab
   sigma[2] = rho[5] * rho[5] + rho[6] * rho[6] + rho[7] * rho[7]; // sigma_bb

   // Exchange values
   xc_gga(&funcForExchange,1,dens,sigma,&exc,vrhoX,vsigmaX,
          v2rho2X,v2rhosigmaX,v2sigma2X,NULL,NULL,NULL,NULL,
          NULL,NULL,NULL,NULL,NULL);

   // Correlation values
   xc_gga(&funcForCorrelation,1,dens,sigma,&exc,vrhoC,vsigmaC,
          v2rho2C,v2rhosigmaC,v2sigma2C,NULL,NULL,NULL,NULL,
          NULL,NULL,NULL,NULL,NULL);

   double DUMGRB, DUMGRC;

   // ALPHA
   DUMGRB=tra[1]*rho[1]+tra[2]*rho[2]+tra[3]*rho[3];
   DUMGRC=tra[1]*rho[5]+tra[2]*rho[6]+tra[3]*rho[7];

   // COEFFICIENTS
   lrCoef[0] = v2rho2X[0]*tra[0]+2.0f*v2rhosigmaX[0]*DUMGRB;
   lrCoef[0] += v2rho2C[0]*tra[0]+2.0f*v2rhosigmaC[0]*DUMGRB+
                v2rhosigmaC[1]*DUMGRC;
   lrCoef[1] = 2.0f*v2rhosigmaX[0]*tra[0]+4.0f*v2sigma2X[0]*DUMGRB;
   lrCoef[1] += 2.0f*v2rhosigmaC[0]*tra[0]+4.0f*v2sigma2C[0]*DUMGRB+
                2.0f*v2sigma2C[1]*DUMGRC;
   lrCoef[2] = v2rhosigmaC[1]*tra[0]+2.0f*v2sigma2C[1]*DUMGRB+
               v2sigma2C[3]*DUMGRC;
   lrCoef[3] = 2.0f*vsigmaX[0] + 2.0f*vsigmaC[0];
   lrCoef[4] = v2rho2C[1]*tra[0]+2.0f*v2rhosigmaC[3]*DUMGRB+
               v2rhosigmaC[4]*DUMGRC;
   lrCoef[5] = 2.0f*v2rhosigmaC[2]*tra[0]+4.0f*v2sigma2C[2]*DUMGRB+
               2.0f*v2sigma2C[4]*DUMGRC;
   lrCoef[6] = v2rhosigmaC[1]*tra[0]+2.0f*v2sigma2C[1]*DUMGRB+
               v2sigma2C[3]*DUMGRC;
   lrCoef[7] = vsigmaC[1];

   // BETA
   DUMGRB=tra[5]*rho[5]+tra[6]*rho[6]+tra[7]*rho[7];
   DUMGRC=tra[5]*rho[1]+tra[6]*rho[2]+tra[7]*rho[3];

   // COEFFICIENTS
   lrCoef[8] = v2rho2X[2]*tra[4]+2.0f*v2rhosigmaX[5]*DUMGRB;
   lrCoef[8] += v2rho2C[2]*tra[4]+2.0f*v2rhosigmaC[5]*DUMGRB+
                v2rhosigmaC[4]*DUMGRC;
   lrCoef[9] = 2.0f*v2rhosigmaX[5]*tra[4]+4.0f*v2sigma2X[5]*DUMGRB;
   lrCoef[9] += 2.0f*v2rhosigmaC[5]*tra[4]+4.0f*v2sigma2C[5]*DUMGRB+
                2.0f*v2sigma2C[4]*DUMGRC;
   lrCoef[10] = v2rhosigmaC[4]*tra[4]+2.0f*v2sigma2C[4]*DUMGRB+
                v2sigma2C[3]*DUMGRC;
   lrCoef[11] = 2.0f*vsigmaX[2] + 2.0f*vsigmaC[2];
   lrCoef[12] = v2rho2C[1]*tra[4]+2.0f*v2rhosigmaC[2]*DUMGRB+
                v2rhosigmaC[1]*DUMGRC;
   lrCoef[13] = 2.0f*v2rhosigmaC[3]*tra[4]+4.0f*v2sigma2C[2]*DUMGRB+
                2.0f*v2sigma2C[1]*DUMGRC;
   lrCoef[14] = v2rhosigmaC[4]*tra[4]+2.0f*v2sigma2C[4]*DUMGRB+
                v2sigma2C[3]*DUMGRC;
   lrCoef[15] = vsigmaC[1];

   return;
}

template <class T, int width>
void LibxcProxy <T, width>::coefZv(double* rho,double* sgm, // inputs
            // first derivatives
            double* vrho,double* vsigma, 
            // second derivatives
            double* v2rho2, double* v2rhosigma,
            double* v2sigma2,
            // third derivatives
            double* v3rho3, double* v3rho2sigma,
            double* v3rhosigma2, double* v3sigma3)
{

   // The otputs for exchange
   double vrhoX[2], vsigmaX[3];
   double v2rho2X[3],v2rhosigmaX[6],v2sigma2X[6];
   double v3rho3X[4], v3rho2sigmaX[9], v3sigma3X[10], v3rhosigma2X[12];

   // The ouputs for correlation
   double vrhoC[2], vsigmaC[3];
   double v2rho2C[3],v2rhosigmaC[6],v2sigma2C[6];
   double v3rho3C[4], v3rho2sigmaC[9], v3sigma3C[10], v3rhosigma2C[12];

   double exc = 0.0f; // Not Reference

   // Get alpha and beta inputs
   double dens[2], sigma[3];
          dens[0] = dens[1] = *rho * 0.5f;
          sigma[0] = sigma[1] = sigma[2] = *sgm * 0.25f;

   // Exchange Values
   xc_gga(&funcForExchange,1,dens,sigma,&exc,vrhoX,vsigmaX,
          v2rho2X,v2rhosigmaX,v2sigma2X,v3rho3X,v3rho2sigmaX,
          v3rhosigma2X,v3sigma3X,
          NULL,NULL,NULL,NULL,NULL);
 
   // Correlation Values
   xc_gga(&funcForCorrelation,1,dens,sigma,&exc,vrhoC,vsigmaC,
          v2rho2C,v2rhosigmaC,v2sigma2C,v3rho3C,v3rho2sigmaC,
          v3rhosigma2C,v3sigma3C,
          NULL,NULL,NULL,NULL,NULL);

   // Results: The first values corresponding to exchange values
   // first derivative
   vrho[0] = vrhoX[0]; vrho[1] = vrhoX[1];
   vsigma[0] = vsigmaX[0]; vsigma[1] = vsigmaX[2];

   vrho[2] = vrhoC[0]; vrho[3] = vrhoC[1];
   vsigma[2] = vsigmaC[0]; vsigma[3] = vsigmaC[2]; vsigma[4] = vsigmaC[1];
 
   // second derivative
   v2rho2[0] = v2rho2X[0]; v2rho2[1] = v2rho2X[2];
   v2rhosigma[0] = v2rhosigmaX[0]; v2rhosigma[1] = v2rhosigmaX[5];
   v2sigma2[0] = v2sigma2X[0]; v2sigma2[1] = v2sigma2X[5];

   v2rho2[2] = v2rho2C[0]; v2rho2[3] = v2rho2C[1]; v2rho2[4] = v2rho2C[2];
   v2rhosigma[2] = v2rhosigmaC[0]; v2rhosigma[3] = v2rhosigmaC[2];
   v2rhosigma[4] = v2rhosigmaC[1]; v2rhosigma[5] = v2rhosigmaC[3];
   v2rhosigma[6] = v2rhosigmaC[5]; v2rhosigma[7] = v2rhosigmaC[4];
   v2sigma2[2] = v2sigma2C[0]; v2sigma2[3] = v2sigma2C[2];
   v2sigma2[4] = v2sigma2C[1]; v2sigma2[5] = v2sigma2C[5];
   v2sigma2[6] = v2sigma2C[4]; v2sigma2[7] = v2sigma2C[3];

   // third derivative
   v3rho3[0] = v3rho3X[0]; v3rho3[1] = v3rho3X[3];
   v3rho2sigma[0] = v3rho2sigmaX[0]; v3rho2sigma[1] = v3rho2sigmaX[8];
   v3rhosigma2[0] = v3rhosigma2X[0]; v3rhosigma2[1] = v3rhosigma2X[11];
   v3sigma3[0] = v3sigma3X[0]; v3sigma3[1] = v3sigma3X[9];
   v3rho3[2] = v3rho3C[2]; v3rho3[3] = v3rho3C[1] - v3rho3C[2] + v3rho3C[3];
   v3rho3[4] = v3rho3[3]; v3rho3[5] = v3rho3C[2];
   // ================== dudas ===============================
   v3rho2sigma[2] = v3rho2sigma[3] = v3rho2sigmaC[0]; 
   v3rho2sigma[4] = v3rho2sigma[7] = 2.0f*v3rho2sigmaC[0]; 
   v3rho2sigma[5] = v3rho2sigma[6] = v3rho2sigmaC[0];
   v3rho2sigma[8] = v3rho2sigma[9] = v3rho2sigmaC[0]; 
   v3rho2sigma[10] = 2.0f*v3rho2sigmaC[0]; 
   // ========================================================
   v3rhosigma2[2] = v3rhosigma2[3] = v3rhosigma2[5] = v3rhosigma2C[0];
   v3rhosigma2[4] = v3rhosigma2[6] = v3rhosigma2[12] = 2.0f*v3rhosigma2C[0];
   v3rhosigma2[7] = v3rhosigma2[13] = 4.0f*v3rhosigma2C[0];
   v3rhosigma2[8] = v3rhosigma2[9] = v3rhosigma2[11] = v3rhosigma2C[0];
   v3rhosigma2[10] = 2.0f*v3rhosigma2C[0];
   v3sigma3[2] = v3sigma3[3] = v3sigma3[5] = v3sigma3C[0];
   v3sigma3[4] = v3sigma3[6] = v3sigma3[9] = 2.0f*v3sigma3C[0];
   v3sigma3[7] = v3sigma3[10] = 4.0f*v3sigma3C[0]; 
   v3sigma3[8] = v3sigma3C[0];
   v3sigma3[11] = 8.0f*v3sigma3C[0]; 

   return;
}

template <class T, int width>
void LibxcProxy <T, width>::doGGA (T rho,
               T sigma,
               T* v2rho2,
               T v2rhosigma,
               T v2sigma2)
{
    const double dens = rho;
    const double sig = sigma;

    // The outputs for exchange
    double v2rho2X = 0.0;
    double v2rhosigmaX = 0.0;
    double v2sigma2X = 0.0;

    // The outputs for correlation
    double v2rho2C = 0.0;
    double v2rhosigmaC = 0.0;
    double v2sigma2C = 0.0;

    // Exchange values
    xc_gga_fxc (&funcForExchange,1,
                &dens,
                &sig,
                &v2rho2X,
                &v2rhosigmaX,
                &v2sigma2X);

    // Correlation values
    xc_gga_fxc (&funcForCorrelation,1,
                &dens,
                &sig,
                &v2rho2C,
                &v2rhosigmaC,
                &v2sigma2C);


    *v2rho2 = v2rho2C + v2rho2X;
    return;
}

#ifdef __CUDACC__

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::joinResults
// Kernel for join the partial results that we obtain from
// calling Libxc.
//
template <class T, int width>
__global__ void joinResults(
		    double* ex, double* exchange,
		    double* ec, double* correlation,
		    double* vrho, double* vrhoC,
		    double* vsigma, double* vsigmaC,
		    double* v2rho, double* v2rhoC,
		    double* v2rhosigma, double* v2rhosigmaC,
		    double* v2sigma, double* v2sigmaC,
		    double* y2a,
		    const double* sigma,
		    const G2G::vec_type<double, width>* grad,
		    const G2G::vec_type<double, width>* hess1,
		    const G2G::vec_type<double, width>* hess2,
		    int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
	ex[i] = exchange[i];
	ec[i] = correlation[i];

	// Merge the results for the derivatives.
	vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];
        // Now, compute y2a value.
	y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
            + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
            + 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x +
					    grad[i].y * grad[i].y * hess1[i].y +
					    grad[i].z * grad[i].z * hess1[i].z +
					    2 * grad[i].x * grad[i].y * hess2[i].x +
					    2 * grad[i].x * grad[i].z * hess2[i].y +
					    2 * grad[i].y * grad[i].z * hess2[i].z));
    }
}

/////////////////////////////////////
// Conversion KERNELS
//
// Utils for data type conversion from lio to libxc
__global__ void convertFloatToDouble(const float* input, double* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (double)input[i];
    }
}

__global__ void convertDoubleToFloat(const double* input, float* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (float)input[i];
    }
}

__global__ void convertFloatToDouble(const G2G::vec_type<float,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (double)(input[i].x);
	output[i].y = (double)(input[i].y);
	output[i].z = (double)(input[i].z);
	output[i].w = (double)(input[i].w);
    }
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertDoubleToFloat(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertFloatToDouble(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertFloatToDouble(const G2G::vec_type<double,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    return;
}

__global__ void convertDoubleToFloat(const G2G::vec_type<double,4>* input, G2G::vec_type<float,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (float)(input[i].x);
	output[i].y = (float)(input[i].y);
	output[i].z = (float)(input[i].z);
	output[i].w = (float)(input[i].w);
    }
}

// end Conversion KERNELS
////////////////////////////////

#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - GPU version
// Calls the xc_gga function from libxc
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// contracted_grad: the contracted grad for libxc
// grad: 
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
//
// Note: all the pointer data are pointers in CUDA memory.
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T* dens,
    const int number_of_points,
    const T* contracted_grad,
    const G2G::vec_type<T, width>* grad,
    const G2G::vec_type<T, width>* hess1,
    const G2G::vec_type<T, width>* hess2,
    T* ex,
    T* ec,
    T* y2a)
{
#ifdef __CUDACC__
    //printf("doGGA - GPU \n");
    //printf("Number of points: %u\n", number_of_points);

    // Este flag esta asi ya que a veces lio utiliza precision mixta
    // y solo en tiempo de ejecucion podemos saber que tipos
    // de datos esta utilizando.
    bool full_double = (sizeof(T) == 8);

    // Variables for the Kernels
    int threadsPerBlock = 256;
    int blocksPerGrid = (number_of_points + threadsPerBlock - 1) / threadsPerBlock;

    cudaError_t err = cudaSuccess;

    // All the arrays for libxc must be of double*
    int array_size = sizeof(double) * number_of_points;
    int vec_size = sizeof(G2G::vec_type<double,width>) * number_of_points;

    double* rho = NULL;
    double* sigma = NULL;

    double* ex_double = NULL;
    double* ec_double = NULL;
    double* y2a_double = NULL;
    G2G::vec_type<double, width>* grad_double = NULL;
    G2G::vec_type<double, width>* hess1_double = NULL;
    G2G::vec_type<double, width>* hess2_double = NULL;

    err = cudaMalloc((void**)&rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void**)&sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device sigma!\n");
	exit(EXIT_FAILURE);
    }

    // Si el tipo de datos es float, creamos los arrays para copiar
    // los inputs y convertirlos a floats.
    if (!full_double) {
	err = cudaMalloc((void**)&ex_double, array_size);
        if (err != cudaSuccess)
        {
	    fprintf(stderr, "Failed to allocate device ex_double!\n");
    	    exit(EXIT_FAILURE);
        }
	cudaMemset(ex_double,0,array_size);

        err = cudaMalloc((void**)&ec_double, array_size);
	if (err != cudaSuccess)
        {
	    fprintf(stderr, "Failed to allocate device ec_double!\n");
	    exit(EXIT_FAILURE);
        }
	cudaMemset(ec_double,0,array_size);

        err = cudaMalloc((void**)&y2a_double, array_size);
	if (err != cudaSuccess)
        {
    	    fprintf(stderr, "Failed to allocate device y2a_double!\n");
	    exit(EXIT_FAILURE);
        }
	cudaMemset(y2a_double,0,array_size);

        err = cudaMalloc((void**)&grad_double, vec_size);
	if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device grad_double!\n");
    	    exit(EXIT_FAILURE);
        }

	err = cudaMalloc((void**)&hess1_double, vec_size);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device hess1_double!\n");
	    exit(EXIT_FAILURE);
        }

	err = cudaMalloc((void**)&hess2_double, vec_size);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device hess2_double!\n");
	    exit(EXIT_FAILURE);
        }
    }

    // Preparamos los datos.
    if (full_double) {
	err = cudaMemcpy(rho, dens, array_size, cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to copy data from dens->rho\n");
        }

	err = cudaMemcpy(sigma, contracted_grad, array_size, cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to copy data from contracted_grad->sigma\n");
	}

        // Usamos los datos como vienen ya que son todos doubles.
	ex_double = (double*)ex;
	ec_double = (double*)ec;
        y2a_double = (double*)y2a;
	grad_double = (G2G::vec_type<double,4>*)grad;
        hess1_double = (G2G::vec_type<double,4>*)hess1;
        hess2_double = (G2G::vec_type<double,4>*)hess2;

    } else {
	// Como los inputs son float, los convertimos para libxc
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(dens, rho, number_of_points);
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(contracted_grad, sigma, number_of_points);
        convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(grad, grad_double, number_of_points);
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(hess1, hess1_double, number_of_points);
        convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(hess2, hess2_double, number_of_points);
    }

    // Preparamos los arrays de salida.
    double* exchange = NULL;
    err = cudaMalloc((void **)&exchange, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exchange!\n");
        exit(EXIT_FAILURE);
    }

    double* correlation = NULL;
    err = cudaMalloc((void **)&correlation, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device correlation!\n");
        exit(EXIT_FAILURE);
    }

    cudaMemset(exchange,0,array_size);
    cudaMemset(correlation,0,array_size);

    // The outputs for exchange
    double* vrho = NULL;
    double* vsigma = NULL;
    double* v2rho = NULL;
    double* v2rhosigma = NULL;
    double* v2sigma = NULL;

    err = cudaMalloc((void **)&vrho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigma!\n");
        exit(EXIT_FAILURE);
    }

    // Clear arrays
    cudaMemset(vrho, 0, array_size);
    cudaMemset(vsigma, 0, array_size);
    cudaMemset(v2rho, 0, array_size);
    cudaMemset(v2rhosigma, 0, array_size);
    cudaMemset(v2sigma, 0, array_size);

    // The outputs for correlation
    double* vrhoC = NULL;
    double* vsigmaC = NULL;
    double* v2rhoC = NULL;
    double* v2rhosigmaC = NULL;
    double* v2sigmaC = NULL;

    err = cudaMalloc((void **)&vrhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigmaC!\n");
        exit(EXIT_FAILURE);
    }

    ///////////////////////////////////
    // Clear arrays for correlation
    cudaMemset(vrhoC, 0, array_size);
    cudaMemset(vsigmaC, 0, array_size);
    cudaMemset(v2rhoC, 0, array_size);
    cudaMemset(v2rhosigmaC, 0, array_size);
    cudaMemset(v2sigmaC, 0, array_size);

    /////////////////////////////
    // Call LIBXC for exchange
    try {
        xc_gga (&funcForExchange, number_of_points,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
        return;
    }

    ////////////////////////////////
    // Call LIBXC for correlation
    try {
        // Now the correlation value.
        xc_gga (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Correlation '%d' \n", exception);
        return;
    }

    ////////////////////////
    // Gather the results
    joinResults<T, width><<<blocksPerGrid, threadsPerBlock>>>(
	ex_double, exchange,
	ec_double, correlation,
	vrho, vrhoC,
	vsigma, vsigmaC,
	v2rho, v2rhoC,
	v2rhosigma, v2rhosigmaC,
	v2sigma, v2sigmaC,
	y2a_double,
	sigma,
	grad_double,
	hess1_double,
	hess2_double,
	number_of_points);

    //////////////////////////
    // Convert if necessary
    if (!full_double) {
	convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (ex_double, ex, number_of_points);
	convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (ec_double, ec, number_of_points);
        convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (y2a_double, y2a, number_of_points);
    }

    /////////////////////////
    // Free device memory
    if (rho != NULL) {
	cudaFree(rho);
    }
    if (sigma != NULL) {
	cudaFree(sigma);
    }
    if (exchange != NULL) {
	cudaFree(exchange);
    }
    if (correlation != NULL) {
	cudaFree(correlation);
    }
    if (vrho != NULL) {
        cudaFree(vrho);
    }
    if (vsigma != NULL) {
	cudaFree(vsigma);
    }
    if (v2rho != NULL) {
	cudaFree(v2rho);
    }
    if (v2rhosigma != NULL) {
	cudaFree(v2rhosigma);
    }
    if (v2sigma != NULL) {
	cudaFree(v2sigma);
    }
    if (vrhoC != NULL) {
        cudaFree(vrhoC);
    }
    if (vsigmaC != NULL) {
	cudaFree(vsigmaC);
    }
    if (v2rhoC != NULL) {
	cudaFree(v2rhoC);
    }
    if (v2rhosigmaC != NULL) {
	cudaFree(v2rhosigmaC);
    }
    if (v2sigmaC != NULL) {
	cudaFree(v2sigmaC);
    }

    if (!full_double) {
        if (ex_double != NULL) {
            cudaFree(ex_double);
        }
        if (ec_double != NULL) {
            cudaFree(ec_double);
        }
        if (y2a_double != NULL) {
        cudaFree(y2a_double);
        }
        if (grad_double != NULL) {
            cudaFree((void*)grad_double);
    	}
    	if (hess1_double != NULL) {
        cudaFree((void*)hess1_double);
    	}
        if (hess2_double != NULL) {
    	    cudaFree((void*)hess2_double);
        }
    }
#endif
    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doLDA - GPU version
// Calls the xc_lda function from libxc
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// contracted_grad: the contracted grad for libxc
// grad:.
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
//
// Note: all the pointer data are pointers in CUDA memory.
//

template <class T, int width>
void LibxcProxy <T, width>::doLDA(T dens, const G2G::vec_type<T, width> &grad, const G2G::vec_type<T, width> &hess1, const G2G::vec_type<T, width> &hess2, T &ex, T &ec, T &y2a)
{
    //TODO: not implemented yet!
    return;
}

#endif // LIBXCPROXY_H
