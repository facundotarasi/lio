#include <iostream>
#include <stdio.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

using namespace G2G;
using namespace std; 

void coef_calculator(double pd,double sigma,double pdx,double pdy,double pdz,
                     double td,double tdx,double tdy,double tdz,double* COEF)
{

// LIBXC INITIALIZATION
   const int nspin = XC_POLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   LibxcProxy<double,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);

// OUTPUTS FOR LIBXC
   double* vrho        = (double*)malloc(4*sizeof(double));
   double* vsigma      = (double*)malloc(5*sizeof(double));
   double* v2rho2      = (double*)malloc(5*sizeof(double));
   double* v2rhosigma  = (double*)malloc(8*sizeof(double));
   double* v2sigma2    = (double*)malloc(8*sizeof(double));
   double* v3rho3      = (double*)malloc(6*sizeof(double));
   double* v3rho2sigma = (double*)malloc(11*sizeof(double));
   double* v3rhosigma2 = (double*)malloc(14*sizeof(double));
   double* v3sigma3    = (double*)malloc(12*sizeof(double));

// LIBXC CALCULATE DERIVATIVES
   libxcProxy.coefZv(&pd,&sigma,
                     vrho,vsigma,
                     v2rho2,v2rhosigma,v2sigma2,
                     v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3);

   //for(int i=0; i<12; i++)
     //printf("%d   %lf\n",i,v3sigma3[i]);

/*
   v3rho2sigma[2] = 1062217532.3723156;
   v3rho2sigma[3] = 1062217532.3723156;
   v3rho2sigma[4] = 2124435064.7202091;
   v3rho2sigma[5] = 1162554460.1091642;
   v3rho2sigma[6] = 1162554460.1091642;
   v3rho2sigma[7] = 2325108920.1885753;
   v3rho2sigma[8] = 1062217532.3722907;
   v3rho2sigma[9] = 1062217532.3722907;
   v3rho2sigma[10] = 2124435064.7201595;
*/

   //exit(-1);

   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[20];

   //pd = 8.6792344676238335E-006;//DEBUG
   //sigma = 6.1329471537828528E-010;//DEBUG

   DUMNV[0] = DUMNV[1] = td;

//-- DUMGRV
//1:A*A 2:B*B 3:A*B 4:B*A
   DUMGRV[0]=DUMGRV[1]=tdx*pdx+tdy*pdy+tdz*pdz;
   DUMGRV[2]=DUMGRV[3]=DUMGRV[0];
   //cout << "DUMGRV[0] " << DUMGRV[0] << endl; // bien

   DUMXX[0]=DUMXX[1]=tdx*tdx+tdy*tdy+tdz*tdz;
   DUMXX[2]=DUMXX[3]=DUMXX[0];
   //cout << "DUMXX " << DUMXX[0] << endl; // bien

//---- FOR G(L2):3RD DERIVATIVES :F_CORE
//********F_CORE****************
//********G_ALPHA***************DUM"A",DUM"A"G
//--AA
   C[0]=2.0f*DUMXX[0];
   C[1]=DUMNV[0]*DUMNV[0];
   C[2]=2.0f*DUMNV[0]*DUMGRV[0];
   C[3]=2.0f*DUMGRV[0]*DUMNV[0];
   C[4]=DUMGRV[0]*DUMNV[0];
   C[5]=DUMNV[0]*DUMGRV[0];
   C[6]=4.0f*DUMGRV[0]*DUMGRV[0];
   C[7]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[8]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[9]=DUMGRV[0]*DUMGRV[0];
/*
      cout << "C0 " << C[0] << endl;
      cout << "C1 " << C[1] << endl;
      cout << "C2 " << C[2] << endl;
      cout << "C3 " << C[3] << endl;
      cout << "C4 " << C[4] << endl;
      cout << "C5 " << C[5] << endl;
      cout << "C6 " << C[6] << endl;
      cout << "C7 " << C[7] << endl;
      cout << "C8 " << C[8] << endl;
      cout << "C9 " << C[9] << endl;
*/  // bien

// -- EXCHANGE
   double XDUMA=0.0f;
   double XDUMAG=0.0f;
// A-AA--1
   //XDUMA=XDUMA+C1*EX(IIPT,KRAGA)
   XDUMA=XDUMA+C[0]*v2rhosigma[0];
   //XDUMAG=XDUMAG+C1*TWO*EX(IIPT,KGAGA)
   XDUMAG=XDUMAG+C[0]*2.0f*v2sigma2[0];
   //cout << "XDUMA,XDUMAG " << XDUMA << " " << XDUMAG << endl;//bien
// --2
   //XDUMA=XDUMA+C2*EX(IIPT,KRARARA)
   XDUMA=XDUMA+C[1]*v3rho3[0];
   //XDUMAG=XDUMAG+C2*TWO*EX(IIPT,KRARAGA)
   XDUMAG=XDUMAG+C[1]*2.0f*v3rho2sigma[0];
   //cout << "XDUMA,XDUMAG " << XDUMA << " " << XDUMAG << endl;//bien
// --3
   //XDUMA=XDUMA+C3*EX(IIPT,KRARAGA)
   XDUMA=XDUMA+C[2]*v3rho2sigma[0];
   //XDUMAG=XDUMAG+C3*TWO*EX(IIPT,KRAGAGA)
   XDUMAG=XDUMAG+C[2]*2.0f*v3rhosigma2[0];
// --4
   //XDUMA=XDUMA+C4*EX(IIPT,KRARAGA)
   XDUMA=XDUMA+C[3]*v3rho2sigma[0];
   //XDUMAG=XDUMAG+C4*TWO*EX(IIPT,KRAGAGA)
   XDUMAG=XDUMAG+C[3]*2.0f*v3rhosigma2[0];
// --7
   //XDUMA=XDUMA+C7*EX(IIPT,KRAGAGA)
   XDUMA=XDUMA+C[6]*v3rhosigma2[0];
   //XDUMAG=XDUMAG+C7*TWO*EX(IIPT,KGAGAGA)
   XDUMAG=XDUMAG+C[6]*2.0f*v3sigma3[0];
   //cout << "XDUMA,XDUMAG " << XDUMA << " " << XDUMAG << endl;//bien
// -- CORRELATION
   double CDUMA=0.0f;
   double CDUMAG1=0.0f;
   double CDUMAG2=0.0f;
// A-AA--1     GA
   //CDUMA=CDUMA+C1*EC(IIPT,IRAGA)
   CDUMA=CDUMA+C[0]*v2rhosigma[2];
   //CDUMAG1=CDUMAG1+C1*TWO*EC(IIPT,IGAGA)
   CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[2];
   //CDUMAG2=CDUMAG2+C1*EC(IIPT,IGAGC)
   CDUMAG2=CDUMAG2+C[0]*v2sigma2[4];
// --2        RARA
   //CDUMA=CDUMA+C2*EC(IIPT,IRARARA)
   CDUMA=CDUMA+C[1]*v3rho3[2];
   //CDUMAG1=CDUMAG1+C2*TWO*EC(IIPT,IRARAGA)
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[2];
   //CDUMAG2=CDUMAG2+C2*EC(IIPT,IRARAGC)
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[4];
// --3        RAGA
   //CDUMA=CDUMA+C3*EC(IIPT,IRARAGA)
   CDUMA=CDUMA+C[2]*v3rho2sigma[2];
   //CDUMAG1=CDUMAG1+C3*TWO*EC(IIPT,IRAGAGA)
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[2];
   //CDUMAG2=CDUMAG2+C3*EC(IIPT,IRAGAGC)
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[4];
// --4        RAGA
   //CDUMA=CDUMA+C4*EC(IIPT,IRARAGA)
   CDUMA=CDUMA+C[3]*v3rho2sigma[2];
   //CDUMAG1=CDUMAG1+C4*TWO*EC(IIPT,IRAGAGA)
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[2];
   //CDUMAG2=CDUMAG2+C4*EC(IIPT,IRAGAGC)
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[4];
// --5        RAGC
   //CDUMA=CDUMA+C5*EC(IIPT,IRARAGC)
   CDUMA=CDUMA+C[4]*v3rho2sigma[4];
   //CDUMAG1=CDUMAG1+C5*TWO*EC(IIPT,IRAGAGC)
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[4];
   //CDUMAG2=CDUMAG2+C5*EC(IIPT,IRAGCGC)
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[7];
// --6        RAGC
   //CDUMA=CDUMA+C6*EC(IIPT,IRARAGC)
   CDUMA=CDUMA+C[5]*v3rho2sigma[4];
   //CDUMAG1=CDUMAG1+C6*TWO*EC(IIPT,IRAGAGC)
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[4];
   //CDUMAG2=CDUMAG2+C6*EC(IIPT,IRAGCGC)
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[7];
// --7        GAGA
   //CDUMA=CDUMA+C7*EC(IIPT,IRAGAGA)
   CDUMA=CDUMA+C[6]*v3rhosigma2[2];
   //CDUMAG1=CDUMAG1+C7*TWO*EC(IIPT,IGAGAGA)
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[2];
   //CDUMAG2=CDUMAG2+C7*EC(IIPT,IGAGAGC)
   CDUMAG2=CDUMAG2+C[6]*v3sigma3[4];
// --8        GAGC
   //CDUMA=CDUMA+C8*EC(IIPT,IRAGAGC)
   CDUMA=CDUMA+C[7]*v3rhosigma2[4];
   //CDUMAG1=CDUMAG1+C8*TWO*EC(IIPT,IGAGAGC)
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[4];
   //CDUMAG2=CDUMAG2+C8*EC(IIPT,IGAGCGC)
   CDUMAG2=CDUMAG2+C[7]*v3sigma3[7];
// --9        GAGC
   //CDUMA=CDUMA+C9*EC(IIPT,IRAGAGC)
   CDUMA=CDUMA+C[8]*v3rhosigma2[4];
   //CDUMAG1=CDUMAG1+C9*TWO*EC(IIPT,IGAGAGC)
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[4];
   //CDUMAG2=CDUMAG2+C9*EC(IIPT,IGAGCGC)
   CDUMAG2=CDUMAG2+C[8]*v3sigma3[7];
// --10       GCGC
   //CDUMA=CDUMA+C10*EC(IIPT,IRAGCGC)
   CDUMA=CDUMA+C[9]*v3rhosigma2[7];
   //CDUMAG1=CDUMAG1+C10*TWO*EC(IIPT,IGAGCGC)
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[7];
   //CDUMAG2=CDUMAG2+C10*EC(IIPT,IGCGCGC)
   CDUMAG2=CDUMAG2+C[9]*v3sigma3[11];
   //cout << "CDUMA,CDUMAG1,CDUMAG2 " << CDUMA << " " << CDUMAG1 << " " << CDUMAG2 << endl;// bien
// --BB
   C[0]=2.0f*DUMXX[1];
   C[1]=DUMNV[1]*DUMNV[1];
   C[2]=2.0f*DUMNV[1]*DUMGRV[1];
   C[3]=2.0f*DUMGRV[1]*DUMNV[1];
   C[4]=DUMGRV[1]*DUMNV[1];
   C[5]=DUMNV[1]*DUMGRV[1];
   C[6]=4.0f*DUMGRV[1]*DUMGRV[1];
   C[7]=2.0f*DUMGRV[1]*DUMGRV[1];
   C[8]=2.0f*DUMGRV[1]*DUMGRV[1];
   C[9]=DUMGRV[1]*DUMGRV[1];
/*
      cout << "C0 " << C[0] << endl;
      cout << "C1 " << C[1] << endl;
      cout << "C2 " << C[2] << endl;
      cout << "C3 " << C[3] << endl;
      cout << "C4 " << C[4] << endl;
      cout << "C5 " << C[5] << endl;
      cout << "C6 " << C[6] << endl;
      cout << "C7 " << C[7] << endl;
      cout << "C8 " << C[8] << endl;
      cout << "C9 " << C[9] << endl;
*/ // bien
// A-BB--1      GB
   //CDUMA=CDUMA+C1*EC(IIPT,IRAGB)
   CDUMA=CDUMA+C[0]*v2rhosigma[3];
   //CDUMAG1=CDUMAG1+C1*TWO*EC(IIPT,IGAGB)
   CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[3];
   //CDUMAG2=CDUMAG2+C1*EC(IIPT,IGBGC)
   CDUMAG2=CDUMAG2+C[0]*v2sigma2[6];
// --2        RBRB
   //CDUMA=CDUMA+C2*EC(IIPT,IRARBRB)
   CDUMA=CDUMA+C[1]*v3rho3[4];
   //CDUMAG1=CDUMAG1+C2*TWO*EC(IIPT,IRBRBGA)
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[8];
   //CDUMAG2=CDUMAG2+C2*EC(IIPT,IRBRBGC)
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[10];
// --3        RBGB
   //CDUMA=CDUMA+C3*EC(IIPT,IRARBGB)
   CDUMA=CDUMA+C[2]*v3rho2sigma[6];
   //CDUMAG1=CDUMAG1+C3*TWO*EC(IIPT,IRBGAGB)
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[9];
   //CDUMAG2=CDUMAG2+C3*EC(IIPT,IRBGBGC)
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[12];
   //cout << "CDUMA,CDUMAG1,CDUMAG2 " << CDUMA << " " << CDUMAG1 << " " << CDUMAG2 << endl;// bien
// --4        RBGB
   //CDUMA=CDUMA+C4*EC(IIPT,IRARBGB)
   CDUMA=CDUMA+C[3]*v3rho2sigma[6];
   //CDUMAG1=CDUMAG1+C4*TWO*EC(IIPT,IRBGAGB)
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[9];
   //CDUMAG2=CDUMAG2+C4*EC(IIPT,IRBGBGC)
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[12];
// --5        RBGC
   //CDUMA=CDUMA+C5*EC(IIPT,IRARBGC)
   CDUMA=CDUMA+C[4]*v3rho2sigma[7];
   //CDUMAG1=CDUMAG1+C5*TWO*EC(IIPT,IRBGAGC)
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[10];
   //CDUMAG2=CDUMAG2+C5*EC(IIPT,IRBGCGC)
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[13];
// --6        RBGC
   //CDUMA=CDUMA+C6*EC(IIPT,IRARBGC)
   CDUMA=CDUMA+C[5]*v3rho2sigma[7];
   //CDUMAG1=CDUMAG1+C6*TWO*EC(IIPT,IRBGAGC)
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[10];
   //CDUMAG2=CDUMAG2+C6*EC(IIPT,IRBGCGC)
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[13];
// --7        GBGB
   //CDUMA=CDUMA+C7*EC(IIPT,IRAGBGB)
   CDUMA=CDUMA+C[6]*v3rhosigma2[5];
   //CDUMAG1=CDUMAG1+C7*TWO*EC(IIPT,IGAGBGB)
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[5];
   //CDUMAG2=CDUMAG2+C7*EC(IIPT,IGBGBGC)
   CDUMAG2=CDUMAG2+C[6]*v3sigma3[9];
// --8        GBGC
   //CDUMA=CDUMA+C8*EC(IIPT,IRAGBGC)
   CDUMA=CDUMA+C[7]*v3rhosigma2[6];
   //CDUMAG1=CDUMAG1+C8*TWO*EC(IIPT,IGAGBGC)
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[6];
   //CDUMAG2=CDUMAG2+C8*EC(IIPT,IGBGCGC)
   CDUMAG2=CDUMAG2+C[7]*v3sigma3[10];
// --9        GBGC
   //CDUMA=CDUMA+C9*EC(IIPT,IRAGBGC)
   CDUMA=CDUMA+C[8]*v3rhosigma2[6];
   //CDUMAG1=CDUMAG1+C9*TWO*EC(IIPT,IGAGBGC)
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[6];
   //CDUMAG2=CDUMAG2+C9*EC(IIPT,IGBGCGC)
   CDUMAG2=CDUMAG2+C[8]*v3sigma3[10];
// --10       GCGC
   //CDUMA=CDUMA+C10*EC(IIPT,IRAGCGC)
   CDUMA=CDUMA+C[9]*v3rhosigma2[7];
   //CDUMAG1=CDUMAG1+C10*TWO*EC(IIPT,IGAGCGC)
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[7];
   //CDUMAG2=CDUMAG2+C10*EC(IIPT,IGCGCGC)
   CDUMAG2=CDUMAG2+C[9]*v3sigma3[11];
// --ABORBA
   C[10]=DUMXX[2];
   C[11]=DUMNV[0]*DUMNV[1];
   C[12]=2.0f*DUMNV[0]*DUMGRV[1];
   C[13]=2.0f*DUMGRV[0]*DUMNV[1];
   C[14]=DUMNV[0]*DUMGRV[3];
   C[15]=DUMGRV[2]*DUMNV[1];
   C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
   C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
   C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
   C[19]=DUMGRV[2]*DUMGRV[3];
// A-AB--11   GC
   //CDUMA=CDUMA+C11*EC(IIPT,IRAGC)
   CDUMA=CDUMA+C[10]*v2rhosigma[4];
   //CDUMAG1=CDUMAG1+C11*TWO*EC(IIPT,IGAGC)
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[4];
   //CDUMAG2=CDUMAG2+C11*EC(IIPT,IGCGC)
   CDUMAG2=CDUMAG2+C[10]*v2sigma2[7];
// --12     RARB
   //CDUMA=CDUMA+C12*EC(IIPT,IRARARB)
   CDUMA=CDUMA+C[11]*v3rho3[3];
   //CDUMAG1=CDUMAG1+C12*TWO*EC(IIPT,IRARBGA)
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[5];
   //CDUMAG2=CDUMAG2+C12*EC(IIPT,IRARBGC)
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[7];
// --13     RAGB
   //CDUMA=CDUMA+C13*EC(IIPT,IRARAGB)
   CDUMA=CDUMA+C[12]*v3rho2sigma[3];
   //CDUMAG1=CDUMAG1+C13*TWO*EC(IIPT,IRAGAGB)
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[3];
   //CDUMAG2=CDUMAG2+C13*EC(IIPT,IRAGBGC)
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[6];
// --14     RBGA
   //CDUMA=CDUMA+C14*EC(IIPT,IRARBGA)
   CDUMA=CDUMA+C[13]*v3rho2sigma[5];
   //CDUMAG1=CDUMAG1+C14*TWO*EC(IIPT,IRBGAGA)
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[8];
   //CDUMAG2=CDUMAG2+C14*EC(IIPT,IRBGAGC)
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[10];
// --15     RAGC
   //CDUMA=CDUMA+C15*EC(IIPT,IRARAGC)
   CDUMA=CDUMA+C[14]*v3rho2sigma[4];
   //CDUMAG1=CDUMAG1+C15*TWO*EC(IIPT,IRAGAGC)
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[4];
   //CDUMAG2=CDUMAG2+C15*EC(IIPT,IRAGCGC)
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[7];
// --16     RBGC
   //CDUMA=CDUMA+C16*EC(IIPT,IRARBGC)
   CDUMA=CDUMA+C[15]*v3rho2sigma[7];
   //CDUMAG1=CDUMAG1+C16*TWO*EC(IIPT,IRBGAGC)
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[10];
   //CDUMAG2=CDUMAG2+C16*EC(IIPT,IRBGCGC)
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[13];
// --17     GAGB
   //CDUMA=CDUMA+C17*EC(IIPT,IRAGAGB)
   CDUMA=CDUMA+C[16]*v3rhosigma2[3];
   //CDUMAG1=CDUMAG1+C17*TWO*EC(IIPT,IGAGAGB)
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[3];
   //CDUMAG2=CDUMAG2+C17*EC(IIPT,IGAGBGC)
   CDUMAG2=CDUMAG2+C[16]*v3sigma3[6];
// --18     GAGC
   //CDUMA=CDUMA+C18*EC(IIPT,IRAGAGC)
   CDUMA=CDUMA+C[17]*v3rhosigma2[4];
   //CDUMAG1=CDUMAG1+C18*TWO*EC(IIPT,IGAGAGC)
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[4];
   //CDUMAG2=CDUMAG2+C18*EC(IIPT,IGAGCGC)
   CDUMAG2=CDUMAG2+C[17]*v3sigma3[7];
// --19     GBGC
   //CDUMA=CDUMA+C19*EC(IIPT,IRAGBGC)
   CDUMA=CDUMA+C[18]*v3rhosigma2[6];
   //CDUMAG1=CDUMAG1+C19*TWO*EC(IIPT,IGAGBGC)
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[6];
   //CDUMAG2=CDUMAG2+C19*EC(IIPT,IGBGCGC)
   CDUMAG2=CDUMAG2+C[18]*v3sigma3[10];
// --20     GCGC
   //CDUMA=CDUMA+C20*EC(IIPT,IRAGCGC)
   CDUMA=CDUMA+C[19]*v3rhosigma2[7];
   //CDUMAG1=CDUMAG1+C20*TWO*EC(IIPT,IGAGCGC)
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[7];
   //CDUMAG2=CDUMAG2+C20*EC(IIPT,IGCGCGC)
   CDUMAG2=CDUMAG2+C[19]*v3sigma3[11];
// **BA_TERM=AB_TERM**
// --ABORBA
   C[10]=DUMXX[2];
   C[11]=DUMNV[0]*DUMNV[1];
   C[12]=2.0f*DUMNV[0]*DUMGRV[1];
   C[13]=2.0f*DUMGRV[0]*DUMNV[1];
   C[14]=DUMNV[0]*DUMGRV[3];
   C[15]=DUMGRV[2]*DUMNV[1];
   C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
   C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
   C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
   C[19]=DUMGRV[2]*DUMGRV[3];
// A-AB--11   GC
   //CDUMA=CDUMA+C11*EC(IIPT,IRAGC)
   CDUMA=CDUMA+C[10]*v2rhosigma[4];
   //CDUMAG1=CDUMAG1+C11*TWO*EC(IIPT,IGAGC)
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[4];
   //CDUMAG2=CDUMAG2+C11*EC(IIPT,IGCGC)
   CDUMAG2=CDUMAG2+C[10]*v2sigma2[7];
// --12     RARB
   //CDUMA=CDUMA+C12*EC(IIPT,IRARARB)
   CDUMA=CDUMA+C[11]*v3rho3[3];
   //CDUMAG1=CDUMAG1+C12*TWO*EC(IIPT,IRARBGA)
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[5];
   //CDUMAG2=CDUMAG2+C12*EC(IIPT,IRARBGC)
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[7];
// --13     RAGB
   //CDUMA=CDUMA+C13*EC(IIPT,IRARAGB)
   CDUMA=CDUMA+C[12]*v3rho2sigma[3];
   //CDUMAG1=CDUMAG1+C13*TWO*EC(IIPT,IRAGAGB)
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[3];
   //CDUMAG2=CDUMAG2+C13*EC(IIPT,IRAGBGC)
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[6];
// --14     RBGA
   //CDUMA=CDUMA+C14*EC(IIPT,IRARBGA)
   CDUMA=CDUMA+C[13]*v3rho2sigma[5];
   //CDUMAG1=CDUMAG1+C14*TWO*EC(IIPT,IRBGAGA)
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[8];
   //CDUMAG2=CDUMAG2+C14*EC(IIPT,IRBGAGC)
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[10];
// --15     RAGC
   //CDUMA=CDUMA+C15*EC(IIPT,IRARAGC)
   CDUMA=CDUMA+C[14]*v3rho2sigma[4];
   //CDUMAG1=CDUMAG1+C15*TWO*EC(IIPT,IRAGAGC)
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[4];
   //CDUMAG2=CDUMAG2+C15*EC(IIPT,IRAGCGC)
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[7];
// --16     RBGC
   //CDUMA=CDUMA+C16*EC(IIPT,IRARBGC)
   CDUMA=CDUMA+C[15]*v3rho2sigma[7];
   //CDUMAG1=CDUMAG1+C16*TWO*EC(IIPT,IRBGAGC)
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[10];
   //CDUMAG2=CDUMAG2+C16*EC(IIPT,IRBGCGC)
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[13];
// --17     GAGB
   //CDUMA=CDUMA+C17*EC(IIPT,IRAGAGB)
   CDUMA=CDUMA+C[16]*v3rhosigma2[3];
   //CDUMAG1=CDUMAG1+C17*TWO*EC(IIPT,IGAGAGB)
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[3];
   //CDUMAG2=CDUMAG2+C17*EC(IIPT,IGAGBGC)
   CDUMAG2=CDUMAG2+C[16]*v3sigma3[6];
// --18     GAGC
   //CDUMA=CDUMA+C18*EC(IIPT,IRAGAGC)
   CDUMA=CDUMA+C[17]*v3rhosigma2[4];
   //CDUMAG1=CDUMAG1+C18*TWO*EC(IIPT,IGAGAGC)
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[4];
   //CDUMAG2=CDUMAG2+C18*EC(IIPT,IGAGCGC)
   CDUMAG2=CDUMAG2+C[17]*v3sigma3[7];
// --19     GBGC
   //CDUMA=CDUMA+C19*EC(IIPT,IRAGBGC)
   CDUMA=CDUMA+C[18]*v3rhosigma2[6];
   //CDUMAG1=CDUMAG1+C19*TWO*EC(IIPT,IGAGBGC)
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[6];
   //CDUMAG2=CDUMAG2+C19*EC(IIPT,IGBGCGC)
   CDUMAG2=CDUMAG2+C[18]*v3sigma3[10];
// --20     GCGC
   //CDUMA=CDUMA+C20*EC(IIPT,IRAGCGC)
   CDUMA=CDUMA+C[19]*v3rhosigma2[7];
   //CDUMAG1=CDUMAG1+C20*TWO*EC(IIPT,IGAGCGC)
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[7];
   //CDUMAG2=CDUMAG2+C20*EC(IIPT,IGCGCGC)
   CDUMAG2=CDUMAG2+C[19]*v3sigma3[11];

// ********G_ALPHA END***********DUM"A",DUM"A"G
// ********G_BETA****************DUM"B",DUM"B"G
// *RB OR GB OR GC
// AA
   C[0]=2.0f*DUMXX[0];
   C[1]=DUMNV[0]*DUMNV[0];
   C[2]=2.0f*DUMNV[0]*DUMGRV[0];
   C[3]=2.0f*DUMGRV[0]*DUMNV[0];
   C[4]=DUMGRV[0]*DUMNV[0];
   C[5]=DUMNV[0]*DUMGRV[0];
   C[6]=4.0f*DUMGRV[0]*DUMGRV[0];
   C[7]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[8]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[9]=DUMGRV[0]*DUMGRV[0];

// -- CORRELATION
      double CDUMB=0.0f;
      double CDUMBG1=0.0f;
      double CDUMBG2=0.0f;
// B-AA--1     GA
      //CDUMB=CDUMB+C1*EC(IIPT,IRBGA)
      CDUMB=CDUMB+C[0]*v2rhosigma[5];
      //CDUMBG1=CDUMBG1+C1*TWO*EC(IIPT,IGAGB)
      CDUMBG1=CDUMBG1+C[0]*2.0f*v2sigma2[3];
      //CDUMBG2=CDUMBG2+C1*EC(IIPT,IGAGC)
      CDUMBG2=CDUMBG2+C[0]*v2sigma2[4];
// --2        RARA
      //CDUMB=CDUMB+C2*EC(IIPT,IRARARB)
      CDUMB=CDUMB+C[1]*v3rho3[3];
      //CDUMBG1=CDUMBG1+C2*TWO*EC(IIPT,IRARAGB)
      CDUMBG1=CDUMBG1+C[1]*2.0f*v3rho2sigma[2];
      //CDUMBG2=CDUMBG2+C2*EC(IIPT,IRARAGC)
      CDUMBG2=CDUMBG2+C[1]*v3rho2sigma[4];
// --3        RAGA
      //CDUMB=CDUMB+C3*EC(IIPT,IRARBGA)
      CDUMB=CDUMB+C[2]*v3rho2sigma[5];
      //CDUMBG1=CDUMBG1+C3*TWO*EC(IIPT,IRAGAGB)
      CDUMBG1=CDUMBG1+C[2]*2.0f*v3rhosigma2[3];
      //CDUMBG2=CDUMBG2+C3*EC(IIPT,IRAGAGC)
      CDUMBG2=CDUMBG2+C[2]*v3rhosigma2[4];
// --4        RAGA
      //CDUMB=CDUMB+C4*EC(IIPT,IRARBGA)
      CDUMB=CDUMB+C[3]*v3rho2sigma[5];
      //CDUMBG1=CDUMBG1+C4*TWO*EC(IIPT,IRAGAGB)
      CDUMBG1=CDUMBG1+C[3]*2.0f*v3rhosigma2[3];
      //CDUMBG2=CDUMBG2+C4*EC(IIPT,IRAGAGC)
      CDUMBG2=CDUMBG2+C[3]*v3rhosigma2[4];
// --5        RAGC
      //CDUMB=CDUMB+C5*EC(IIPT,IRARBGC)
      CDUMB=CDUMB+C[4]*v3rho2sigma[7];
      //CDUMBG1=CDUMBG1+C5*TWO*EC(IIPT,IRAGBGC)
      CDUMBG1=CDUMBG1+C[4]*2.0f*v3rhosigma2[6];
      //CDUMBG2=CDUMBG2+C5*EC(IIPT,IRAGCGC)
      CDUMBG2=CDUMBG2+C[4]*v3rhosigma2[7];
// --6        RAGC
      //CDUMB=CDUMB+C6*EC(IIPT,IRARBGC)
      CDUMB=CDUMB+C[5]*v3rho2sigma[7];
      //CDUMBG1=CDUMBG1+C6*TWO*EC(IIPT,IRAGBGC)
      CDUMBG1=CDUMBG1+C[5]*2.0f*v3rhosigma2[6];
      //CDUMBG2=CDUMBG2+C6*EC(IIPT,IRAGCGC)
      CDUMBG2=CDUMBG2+C[5]*v3rhosigma2[7];
// --7        GAGA
      //CDUMB=CDUMB+C7*EC(IIPT,IRBGAGA)
      CDUMB=CDUMB+C[6]*v3rhosigma2[8];
      //CDUMBG1=CDUMBG1+C7*TWO*EC(IIPT,IGAGAGB)
      CDUMBG1=CDUMBG1+C[6]*v3sigma3[3];
      //CDUMBG2=CDUMBG2+C7*EC(IIPT,IGAGAGC)
      CDUMBG2=CDUMBG2+C[6]*v3sigma3[4];
// --8        GAGC
      //CDUMB=CDUMB+C8*EC(IIPT,IRBGAGC)
      CDUMB=CDUMB+C[7]*v3rhosigma2[10];
      //CDUMBG1=CDUMBG1+C8*TWO*EC(IIPT,IGAGBGC)
      CDUMBG1=CDUMBG1+C[7]*2.0f*v3sigma3[6];
      //CDUMBG2=CDUMBG2+C8*EC(IIPT,IGAGCGC)
      CDUMBG2=CDUMBG2+C[7]*v3sigma3[7];
// --9        GAGC
      //CDUMB=CDUMB+C9*EC(IIPT,IRBGAGC)
      CDUMB=CDUMB+C[8]*v3rhosigma2[10];
      //CDUMBG1=CDUMBG1+C9*TWO*EC(IIPT,IGAGBGC)
      CDUMBG1=CDUMBG1+C[8]*2.0f*v3sigma3[6];
      //CDUMBG2=CDUMBG2+C9*EC(IIPT,IGAGCGC)
      CDUMBG2=CDUMBG2+C[8]*2.0f*v3sigma3[7];
// --10       GCGC
      //CDUMB=CDUMB+C10*EC(IIPT,IRBGCGC)
      CDUMB=CDUMB+C[9]*v3rhosigma2[13];
      //CDUMBG1=CDUMBG1+C10*TWO*EC(IIPT,IGBGCGC)
      CDUMBG1=CDUMBG1+C[9]*2.0f*v3sigma3[10];
      //CDUMBG2=CDUMBG2+C10*EC(IIPT,IGCGCGC)
      CDUMBG2=CDUMBG2+C[9]*v3sigma3[11];
// --BB
      C[0]=2.0f*DUMXX[1];
      C[1]=DUMNV[1]*DUMNV[1];
      C[2]=2.0f*DUMNV[1]*DUMGRV[1];
      C[3]=2.0f*DUMGRV[1]*DUMNV[1];
      C[4]=DUMGRV[1]*DUMNV[1];
      C[5]=DUMNV[1]*DUMGRV[1];
      C[6]=4.0f*DUMGRV[1]*DUMGRV[1];
      C[7]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[8]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[9]=DUMGRV[1]*DUMGRV[1];
// -- EXCHANGE
      double XDUMB=0.0f;
      double XDUMBG=0.0f;
// B-BB--1      GB
      //XDUMB=XDUMB+C1*EX(IIPT,KRBGB)
      XDUMB=XDUMB+C[0]*v2rhosigma[1];
      //XDUMBG=XDUMBG+C1*TWO*EX(IIPT,KGBGB)
      XDUMBG=XDUMBG+C[0]*v2sigma2[1];
// --2          RBRB
      //XDUMB=XDUMB+C2*EX(IIPT,KRBRBRB)
      XDUMB=XDUMB+C[1]*v3rho3[1];
      //XDUMBG=XDUMBG+C2*TWO*EX(IIPT,KRBRBGB)
      XDUMBG=XDUMBG+C[1]*2.0f*v3rho2sigma[1];
// --3          RBGB
      //XDUMB=XDUMB+C3*EX(IIPT,KRBRBGB)
      XDUMB=XDUMB+C[2]*v3rho2sigma[1];
      //XDUMBG=XDUMBG+C3*TWO*EX(IIPT,KRBGBGB)
      XDUMBG=XDUMBG+C[2]*2.0f*v3rhosigma2[1];
// --4          RBGB
      //XDUMB=XDUMB+C4*EX(IIPT,KRBRBGB)
      XDUMB=XDUMB+C[3]*v3rho2sigma[1];
      //XDUMBG=XDUMBG+C4*TWO*EX(IIPT,KRBGBGB)
      XDUMBG=XDUMBG+C[3]*2.0f*v3rhosigma2[1];
// --7          GBGB
      //XDUMB=XDUMB+C7*EX(IIPT,KRBGBGB)
      XDUMB=XDUMB+C[6]*v3rhosigma2[1];
      //XDUMBG=XDUMBG+C7*TWO*EX(IIPT,KGBGBGB)
      XDUMBG=XDUMBG+C[6]*2.0f*v3sigma3[1];
// B-BB--1      GB
      //CDUMB=CDUMB+C1*EC(IIPT,IRBGB)
      CDUMB=CDUMB+C[0]*v2rhosigma[6];
      //CDUMBG1=CDUMBG1+C1*TWO*EC(IIPT,IGBGB)
      CDUMBG1=CDUMBG1+C[0]*2.0f*v2sigma2[5];
      //CDUMBG2=CDUMBG2+C1*EC(IIPT,IGBGC)
      CDUMBG2=CDUMBG2+C[0]*v2sigma2[6];
// --2        RBRB
      //CDUMB=CDUMB+C2*EC(IIPT,IRBRBRB)
      CDUMB=CDUMB+C[1]*v3rho3[5];
      //CDUMBG1=CDUMBG1+C2*TWO*EC(IIPT,IRBRBGB)
      CDUMBG1=CDUMBG1+C[1]*2.0f*v3rho2sigma[9];
      //CDUMBG2=CDUMBG2+C2*EC(IIPT,IRBRBGC)
      CDUMBG2=CDUMBG2+C[1]*v3rho2sigma[10];
// --3        RBGB
      //CDUMB=CDUMB+C3*EC(IIPT,IRBRBGB)
      CDUMB=CDUMB+C[2]*v3rho2sigma[9];
      //CDUMBG1=CDUMBG1+C3*TWO*EC(IIPT,IRBGBGB)
      CDUMBG1=CDUMBG1+C[2]*2.0f*v3rhosigma2[11];
      //CDUMBG2=CDUMBG2+C3*EC(IIPT,IRBGBGC)
      CDUMBG2=CDUMBG2+C[2]*v3rhosigma2[12];
// --4        RBGB
      //CDUMB=CDUMB+C4*EC(IIPT,IRBRBGB)
      CDUMB=CDUMB+C[3]*v3rho2sigma[9];
      //CDUMBG1=CDUMBG1+C4*TWO*EC(IIPT,IRBGBGB)
      CDUMBG1=CDUMBG1+C[3]*2.0f*v3rhosigma2[11];
      //CDUMBG2=CDUMBG2+C4*EC(IIPT,IRBGBGC)
      CDUMBG2=CDUMBG2+C[3]*v3rhosigma2[12];
// --5        RBGC
      //CDUMB=CDUMB+C5*EC(IIPT,IRBRBGC)
      CDUMB=CDUMB+C[4]*v3rho2sigma[10];
      //CDUMBG1=CDUMBG1+C5*TWO*EC(IIPT,IRBGBGC)
      CDUMBG1=CDUMBG1+C[4]*2.0f*v3rhosigma2[12];
      //CDUMBG2=CDUMBG2+C5*EC(IIPT,IRBGCGC)
      CDUMBG2=CDUMBG2+C[4]*v3rhosigma2[13];
// --6        RBGC
      //CDUMB=CDUMB+C6*EC(IIPT,IRBRBGC)
      CDUMB=CDUMB+C[5]*v3rho2sigma[10];
      //CDUMBG1=CDUMBG1+C6*TWO*EC(IIPT,IRBGBGC)
      CDUMBG1=CDUMBG1+C[5]*2.0f*v3rhosigma2[12];
      //CDUMBG2=CDUMBG2+C6*EC(IIPT,IRBGCGC)
      CDUMBG2=CDUMBG2+C[5]*v3rhosigma2[13];
// --7        GBGB
      //CDUMB=CDUMB+C7*EC(IIPT,IRBGBGB)
      CDUMB=CDUMB+C[6]*v3rhosigma2[11];
      //CDUMBG1=CDUMBG1+C7*TWO*EC(IIPT,IGBGBGB)
      CDUMBG1=CDUMBG1+C[6]*2.0f*v3sigma3[8];
      //CDUMBG2=CDUMBG2+C7*EC(IIPT,IGBGBGC)
      CDUMBG2=CDUMBG2+C[6]*v3sigma3[9];
// --8        GBGC
      //CDUMB=CDUMB+C8*EC(IIPT,IRBGBGC)
      CDUMB=CDUMB+C[7]*v3rhosigma2[12];
      //CDUMBG1=CDUMBG1+C8*TWO*EC(IIPT,IGBGBGC)
      CDUMBG1=CDUMBG1+C[7]*2.0f*v3sigma3[9];
      //CDUMBG2=CDUMBG2+C8*EC(IIPT,IGBGCGC)
      CDUMBG2=CDUMBG2+C[7]*v3sigma3[10];
// --9        GBGC
      //CDUMB=CDUMB+C9*EC(IIPT,IRBGBGC)
      CDUMB=CDUMB+C[8]*v3rhosigma2[12];
      //CDUMBG1=CDUMBG1+C9*TWO*EC(IIPT,IGBGBGC)
      CDUMBG1=CDUMBG1+C[8]*2.0f*v3sigma3[9];
      //CDUMBG2=CDUMBG2+C9*EC(IIPT,IGBGCGC)
      CDUMBG2=CDUMBG2+C[8]*v3sigma3[10];
// --10       GCGC
      //CDUMB=CDUMB+C10*EC(IIPT,IRBGCGC)
      CDUMB=CDUMB+C[9]*v3rhosigma2[13];
      //CDUMBG1=CDUMBG1+C10*TWO*EC(IIPT,IGBGCGC)
      CDUMBG1=CDUMBG1+C[9]*2.0f*v3sigma3[10];
      //CDUMBG2=CDUMBG2+C10*EC(IIPT,IGCGCGC)
      CDUMBG2=CDUMBG2+C[9]*v3sigma3[11];
// --ABORBA
      C[10]=DUMXX[2];
      C[11]=DUMNV[0]*DUMNV[1];
      C[12]=2.0f*DUMNV[0]*DUMGRV[1];
      C[13]=2.0f*DUMGRV[0]*DUMNV[1];
      C[14]=DUMNV[0]*DUMGRV[3];
      C[15]=DUMGRV[2]*DUMNV[1];
      C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
      C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
      C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
      C[19]=DUMGRV[2]*DUMGRV[3];
// B-AB--11   GC
      //CDUMB=CDUMB+C11*EC(IIPT,IRBGC)
      CDUMB=CDUMB+C[10]*v2rhosigma[7];
      //CDUMBG1=CDUMBG1+C11*TWO*EC(IIPT,IGBGC)
      CDUMBG1=CDUMBG1+C[10]*2.0f*v2sigma2[6];
      //CDUMBG2=CDUMBG2+C11*EC(IIPT,IGCGC)
      CDUMBG2=CDUMBG2+C[10]*v2sigma2[7];
// --12     RARB
      //CDUMB=CDUMB+C12*EC(IIPT,IRARBRB)
      CDUMB=CDUMB+C[11]*v3rho3[4];
      //CDUMBG1=CDUMBG1+C12*TWO*EC(IIPT,IRARBGB)
      CDUMBG1=CDUMBG1+C[11]*2.0f*v3rho2sigma[6];
      //CDUMBG2=CDUMBG2+C12*EC(IIPT,IRARBGC)
      CDUMBG2=CDUMBG2+C[11]*v3rho2sigma[7];
// --13     RAGB
      //CDUMB=CDUMB+C13*EC(IIPT,IRARBGB)
      CDUMB=CDUMB+C[12]*v3rho2sigma[6];
      //CDUMBG1=CDUMBG1+C13*TWO*EC(IIPT,IRAGBGB)
      CDUMBG1=CDUMBG1+C[12]*2.0f*v3rhosigma2[5];
      //CDUMBG2=CDUMBG2+C13*EC(IIPT,IRAGBGC)
      CDUMBG2=CDUMBG2+C[12]*v3rhosigma2[6];
// --14     RBGA
      //CDUMB=CDUMB+C14*EC(IIPT,IRBRBGA)
      CDUMB=CDUMB+C[13]*v3rho2sigma[8];
      //CDUMBG1=CDUMBG1+C14*TWO*EC(IIPT,IRBGAGB)
      CDUMBG1=CDUMBG1+C[13]*2.0f*v3rhosigma2[9];
      //CDUMBG2=CDUMBG2+C14*EC(IIPT,IRBGAGC)
      CDUMBG2=CDUMBG2+C[13]*v3rhosigma2[10];
// --15     RAGC
      //CDUMB=CDUMB+C15*EC(IIPT,IRARBGC)
      CDUMB=CDUMB+C[14]*v3rho2sigma[7];
      //CDUMBG1=CDUMBG1+C15*TWO*EC(IIPT,IRAGBGC)
      CDUMBG1=CDUMBG1+C[14]*2.0f*v3rhosigma2[6];
      //CDUMBG2=CDUMBG2+C15*EC(IIPT,IRAGCGC)
      CDUMBG2=CDUMBG2+C[14]*v3rhosigma2[7];
// --16     RBGC
      //CDUMB=CDUMB+C16*EC(IIPT,IRBRBGC)
      CDUMB=CDUMB+C[15]*v3rho2sigma[10];
      //CDUMBG1=CDUMBG1+C16*TWO*EC(IIPT,IRBGBGC)
      CDUMBG1=CDUMBG1+C[15]*2.0f*v3rhosigma2[12];
      //CDUMBG2=CDUMBG2+C16*EC(IIPT,IRBGCGC)
      CDUMBG2=CDUMBG2+C[15]*v3rhosigma2[13];
// --17     GAGB
      //CDUMB=CDUMB+C17*EC(IIPT,IRBGAGB)
      CDUMB=CDUMB+C[16]*v3rhosigma2[9];
      //CDUMBG1=CDUMBG1+C17*TWO*EC(IIPT,IGAGBGB)
      CDUMBG1=CDUMBG1+C[16]*2.0f*v3sigma3[5];
      //CDUMBG2=CDUMBG2+C17*EC(IIPT,IGAGBGC)
      CDUMBG2=CDUMBG2+C[16]*v3sigma3[6];
// --18     GAGC
      //CDUMB=CDUMB+C18*EC(IIPT,IRBGAGC)
      CDUMB=CDUMB+C[17]*v3rhosigma2[10];
      //CDUMBG1=CDUMBG1+C18*TWO*EC(IIPT,IGAGBGC)
      CDUMBG1=CDUMBG1+C[17]*2.0f*v3sigma3[6];
      //CDUMBG2=CDUMBG2+C18*EC(IIPT,IGAGCGC)
      CDUMBG2=CDUMBG2+C[17]*v3sigma3[7];
// --19     GBGC
      //CDUMB=CDUMB+C19*EC(IIPT,IRBGBGC)
      CDUMB=CDUMB+C[18]*v3rhosigma2[12];
      //CDUMBG1=CDUMBG1+C19*TWO*EC(IIPT,IGBGBGC)
      CDUMBG1=CDUMBG1+C[18]*2.0f*v3sigma3[9];
      //CDUMBG2=CDUMBG2+C19*EC(IIPT,IGBGCGC)
      CDUMBG2=CDUMBG2+C[18]*v3sigma3[10];
// --20     GCGC
      //CDUMB=CDUMB+C20*EC(IIPT,IRBGCGC)
      CDUMB=CDUMB+C[19]*v3rhosigma2[13];
      //CDUMBG1=CDUMBG1+C20*TWO*EC(IIPT,IGBGCGC)
      CDUMBG1=CDUMBG1+C[19]*2.0f*v3sigma3[10];
      //CDUMBG2=CDUMBG2+C20*EC(IIPT,IGCGCGC)
      CDUMBG2=CDUMBG2+C[19]*v3sigma3[11];

// **BA_TERM=AB_TERM**
// B-BA--11   GC
      //CDUMB=CDUMB+C11*EC(IIPT,IRBGC)
      CDUMB=CDUMB+C[10]*v2rhosigma[7];
      //CDUMBG1=CDUMBG1+C11*TWO*EC(IIPT,IGBGC)
      CDUMBG1=CDUMBG1+C[10]*2.0f*v2sigma2[6];
      //CDUMBG2=CDUMBG2+C11*EC(IIPT,IGCGC)
      CDUMBG2=CDUMBG2+C[10]*v2sigma2[7];
// --12     RARB
      //CDUMB=CDUMB+C12*EC(IIPT,IRARBRB)
      CDUMB=CDUMB+C[11]*v3rho3[4];
      //CDUMBG1=CDUMBG1+C12*TWO*EC(IIPT,IRARBGB)
      CDUMBG1=CDUMBG1+C[11]*2.0f*v3rho2sigma[6];
      //CDUMBG2=CDUMBG2+C12*EC(IIPT,IRARBGC)
      CDUMBG2=CDUMBG2+C[11]*v3rho2sigma[7];
// --13     RAGB
      //CDUMB=CDUMB+C13*EC(IIPT,IRARBGB)
      CDUMB=CDUMB+C[12]*v3rho2sigma[6];
      //CDUMBG1=CDUMBG1+C13*TWO*EC(IIPT,IRAGBGB)
      CDUMBG1=CDUMBG1+C[12]*2.0f*v3rhosigma2[5];
      //CDUMBG2=CDUMBG2+C13*EC(IIPT,IRAGBGC)
      CDUMBG2=CDUMBG2+C[12]*v3rhosigma2[6];
// --14     RBGA
      //CDUMB=CDUMB+C14*EC(IIPT,IRBRBGA)
      CDUMB=CDUMB+C[13]*v3rho2sigma[8];
      //CDUMBG1=CDUMBG1+C14*TWO*EC(IIPT,IRBGAGB)
      CDUMBG1=CDUMBG1+C[13]*2.0f*v3rhosigma2[9];
      //CDUMBG2=CDUMBG2+C14*EC(IIPT,IRBGAGC)
      CDUMBG2=CDUMBG2+C[13]*v3rhosigma2[10];
// --15     RAGC
      //CDUMB=CDUMB+C15*EC(IIPT,IRARBGC)
      CDUMB=CDUMB+C[14]*v3rho2sigma[7];
      //CDUMBG1=CDUMBG1+C15*TWO*EC(IIPT,IRAGBGC)
      CDUMBG1=CDUMBG1+C[14]*2.0f*v3rhosigma2[6];
      //CDUMBG2=CDUMBG2+C15*EC(IIPT,IRAGCGC)
      CDUMBG2=CDUMBG2+C[14]*v3rhosigma2[7];
// --16     RBGC
      //CDUMB=CDUMB+C16*EC(IIPT,IRBRBGC)
      CDUMB=CDUMB+C[15]*v3rho2sigma[10];
      //CDUMBG1=CDUMBG1+C16*TWO*EC(IIPT,IRBGBGC)
      CDUMBG1=CDUMBG1+C[15]*2.0f*v3rhosigma2[12];
      //CDUMBG2=CDUMBG2+C16*EC(IIPT,IRBGCGC)
      CDUMBG2=CDUMBG2+C[15]*v3rhosigma2[13];
// --17     GAGB
      //CDUMB=CDUMB+C17*EC(IIPT,IRBGAGB)
      CDUMB=CDUMB+C[16]*v3rhosigma2[9];
      //CDUMBG1=CDUMBG1+C17*TWO*EC(IIPT,IGAGBGB)
      CDUMBG1=CDUMBG1+C[16]*2.0f*v3sigma3[5];
      //CDUMBG2=CDUMBG2+C17*EC(IIPT,IGAGBGC)
      CDUMBG2=CDUMBG2+C[16]*v3sigma3[6];
// --18     GAGC
      //CDUMB=CDUMB+C18*EC(IIPT,IRBGAGC)
      CDUMB=CDUMB+C[17]*v3rhosigma2[10];
      //CDUMBG1=CDUMBG1+C18*TWO*EC(IIPT,IGAGBGC)
      CDUMBG1=CDUMBG1+C[17]*2.0f*v3sigma3[6];
      //CDUMBG2=CDUMBG2+C18*EC(IIPT,IGAGCGC)
      CDUMBG2=CDUMBG2+C[17]*v3sigma3[7];
// --19     GBGC
      //CDUMB=CDUMB+C19*EC(IIPT,IRBGBGC)
      CDUMB=CDUMB+C[18]*v3rhosigma2[12];
      //CDUMBG1=CDUMBG1+C19*TWO*EC(IIPT,IGBGBGC)
      CDUMBG1=CDUMBG1+C[18]*2.0f*v3sigma3[9];
      //CDUMBG2=CDUMBG2+C19*EC(IIPT,IGBGCGC)
      CDUMBG2=CDUMBG2+C[18]*v3sigma3[10];
// --20     GCGC
      //CDUMB=CDUMB+C20*EC(IIPT,IRBGCGC)
      CDUMB=CDUMB+C[19]*v3rhosigma2[13];
      //CDUMBG1=CDUMBG1+C20*TWO*EC(IIPT,IGBGCGC)
      CDUMBG1=CDUMBG1+C[19]*2.0f*v3sigma3[10];
      //CDUMBG2=CDUMBG2+C20*EC(IIPT,IGCGCGC)
      CDUMBG2=CDUMBG2+C[19]*v3sigma3[11];

// ********G_BETA END************DUM"B",DUM"B"G
// *****F_CORE_END**********************
// ********EDGE****************
// --EXCHANGE
// --GA*
      double XDUMAGEA=0.0f;
// --E1 (GA)*RA
      //XDUMAGEA=XDUMAGEA+TWO*TWO*DUMNV(1)*EX(IIPT,KRAGA)
      XDUMAGEA=XDUMAGEA+4.0f*DUMNV[0]*v2rhosigma[0];
// --E2 (GA)*GA
     // XDUMAGEA=XDUMAGEA+TWO*TWO*TWO*DUMGRV(2)*EX(IIPT,KGAGA)
      XDUMAGEA=XDUMAGEA+8.0f*DUMGRV[1]*v2sigma2[0];
// --GB
      double XDUMAGEB=0.0f;
// --E1 (GB)*RB
      //XDUMAGEB=XDUMAGEB+TWO*TWO*DUMNV(2)*EX(IIPT,KRBGB)
      XDUMAGEB=XDUMAGEB+4.0f*DUMNV[1]*v2rhosigma[1];
// --E2 (GB)*GB
      //XDUMAGEB=XDUMAGEB+TWO*TWO*TWO*DUMGRV(2)*EX(IIPT,KGBGB)
      XDUMAGEB=XDUMAGEB+8.0f*DUMGRV[1]*v2sigma2[1];
// --CORRELATION
      double CDUMAGEA=0.0f;
// --GA
// --E1A GA*RA
      //CDUMAGEA=CDUMAGEA+TWO*DUMNV(1)*TWO*EC(IIPT,IRAGA)
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[0]*v2rhosigma[2];
// --E2A GA*GA
      //CDUMAGEA=CDUMAGEA+TWO*DUMGRV(1)*TWO*TWO*EC(IIPT,IGAGA)
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[0]*v2sigma2[2];
// --E3A GA*GC !
      //CDUMAGEA=CDUMAGEA+TWO*DUMGRV(3)*TWO*EC(IIPT,IGAGC)
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[2]*v2sigma2[4];
// --E1B GA*RB
      //CDUMAGEA=CDUMAGEA+TWO*DUMNV(2)*TWO*EC(IIPT,IRBGA)
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[1]*v2rhosigma[5];
// --E2B GA*GB
      //CDUMAGEA=CDUMAGEA+TWO*DUMGRV(2)*TWO*TWO*EC(IIPT,IGAGB)
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[1]*v2sigma2[3];
// --E3B GA*GC
      //CDUMAGEA=CDUMAGEA+TWO*DUMGRV(4)*TWO*EC(IIPT,IGAGC)
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[3]*v2sigma2[4];
      double CDUMAGEB=0.0f;
// --GB
// --E1A GB*RA
     //CDUMAGEB=CDUMAGEB+TWO*DUMNV(1)*TWO*EC(IIPT,IRAGB)
      CDUMAGEB=CDUMAGEB+4.0f*DUMNV[0]*v2rhosigma[3];
// --E2A GB*GA
      //CDUMAGEB=CDUMAGEB+TWO*DUMGRV(1)*TWO*TWO*EC(IIPT,IGAGB)
      CDUMAGEB=CDUMAGEB+8.0f*DUMGRV[0]*v2sigma2[3];
// --E3A GB*GC !
      //CDUMAGEB=CDUMAGEB+TWO*DUMGRV(3)*TWO*EC(IIPT,IGBGC)
      CDUMAGEB=CDUMAGEB+4.0f*DUMGRV[2]*v2sigma2[6];
// --E1B GB*RB
      //CDUMAGEB=CDUMAGEB+TWO*DUMNV(2)*TWO*EC(IIPT,IRBGB)
      CDUMAGEB=CDUMAGEB+4.0f*DUMNV[1]*v2rhosigma[6];
// --E2B GB*GB
      //CDUMAGEB=CDUMAGEB+TWO*DUMGRV(2)*TWO*TWO*EC(IIPT,IGBGB)
      CDUMAGEB=CDUMAGEB+8.0f*DUMGRV[1]*v2sigma2[5];
// --E3B GB*GC
      //CDUMAGEB=CDUMAGEB+TWO*DUMGRV(4)*TWO*EC(IIPT,IGBGC)
      CDUMAGEB=CDUMAGEB+4.0f*DUMGRV[3]*v2sigma2[6];

      double CDUMAGEC=0.0f;
// --GC*
// --E1A GC*RA
      //CDUMAGEC=CDUMAGEC+TWO*DUMNV(1)*EC(IIPT,IRAGC)
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[0]*v2rhosigma[4];
// --E2A GCA*GA
      //CDUMAGEC=CDUMAGEC+TWO*DUMGRV(1)*TWO*EC(IIPT,IGAGC)
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[0]*v2sigma2[4];
// --E3A GCA*GC
      //CDUMAGEC=CDUMAGEC+TWO*DUMGRV(3)*EC(IIPT,IGCGC)
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[2]*v2sigma2[7];
// --E1B GCA*RB
      //CDUMAGEC=CDUMAGEC+TWO*DUMNV(2)*EC(IIPT,IRBGC)
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[1]*v2rhosigma[7];
// --E2B GCA*GB
      //CDUMAGEC=CDUMAGEC+TWO*DUMGRV(2)*TWO*EC(IIPT,IGBGC)
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[1]*v2sigma2[6];
// --E3B GCA*GC
      //CDUMAGEC=CDUMAGEC+TWO*DUMGRV(4)*EC(IIPT,IGCGC)
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[3]*v2sigma2[7];
// ********EDGE END************
// --CONTRUCTION
// ALPHA
      double DUM1A=XDUMA+CDUMA;
      double DUM2A=XDUMAG+CDUMAG1+CDUMAG2;
      double DUM3A=XDUMAGEA+CDUMAGEA+CDUMAGEC;
      COEF[0] = DUM1A;
      COEF[1] = DUM2A;
      COEF[2] = DUM3A;

// FREE MEMORY FROM LIBXC OUTPUT'S
      free(vrho);
      free(vsigma);
      free(v2rho2);
      free(v2rhosigma);
      free(v2sigma2);
      free(v3rho3);
      free(v3rho2sigma);
      free(v3rhosigma2);
      free(v3sigma3);
}
