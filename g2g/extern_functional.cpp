#include <iostream>

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "init.h"
#include "partition.h"

using namespace G2G;


#if USE_LIBXC
#include "functional.h"

extern "C" void g2g_extern_functional_(int& main_id, bool* externFunc,
                              int* HF, double* HF_fac, double* screen)
{
   if ( *externFunc == 0 ) return;
   cout << " " << endl;
   cout << " Extern Functional Module " << endl;

   // Allocate and Free Memory
   if ( fortran_vars.func_id != NULL ) {
      free(fortran_vars.func_id); fortran_vars.func_id = NULL;
   }
   if ( fortran_vars.func_coef != NULL ) {
      free(fortran_vars.func_coef); fortran_vars.func_coef = NULL;
   }
   if ( fortran_vars.HF != NULL ) {
      free(fortran_vars.HF); fortran_vars.HF = NULL;
   }
   if ( fortran_vars.HF_fac != NULL ) {
      free(fortran_vars.HF_fac); fortran_vars.HF_fac = NULL;
   }
   
   switch (main_id) {
      // PBE
      case 101:
           set_pbe(HF,HF_fac,screen); break;

      case 130:
           set_pbe(HF,HF_fac,screen); break;

      // PBE0
      case 406:
           set_pbe0(HF,HF_fac,screen); break;

      // B3LYP
      case 402:
           set_b3lyp(HF,HF_fac,screen); break;

      // CAM-B3LYP
      case 433:
           set_cam_b3lyp(HF,HF_fac,screen); break;

      // LC-WPBE
      case 478:
           set_lc_wpbe(HF,HF_fac,screen); break;

      // LC-BLYP
      case 400:
           set_lc_blyp(HF,HF_fac,screen); break;

      // WB97X
      case 464:
           set_wb97x(HF,HF_fac,screen); break;
   
      default:
           cout << "The Functional id " << main_id << " doesn't implemented yet" << endl;
           exit(-1); break;
   }
   cout << " " << endl;
}
#else
extern "C" void g2g_extern_functional_(int& main_id, bool* externFunc,
                              int* HF, double* HF_fac, double* screen)
{
   fortran_vars.fexc = 1.0f;
   if ( *externFunc == 0 ) return;
   cout << " " << endl;
   cout << " Extern Functional Module " << endl;

   // Allocate and Free Memory
   if ( fortran_vars.HF != NULL ) {
      free(fortran_vars.HF); fortran_vars.HF = NULL;
   }
   if ( fortran_vars.HF_fac != NULL ) {
      free(fortran_vars.HF_fac); fortran_vars.HF_fac = NULL;
   }
   fortran_vars.HF     = (int*   ) malloc(sizeof(double)*3);
   fortran_vars.HF_fac = (double*) malloc(sizeof(double)*3);

   if ( main_id == 406 ) {
      fortran_vars.fexc = 0.75f;
      fortran_vars.HF[0] = HF[0] = 1;
      fortran_vars.HF[1] = HF[1] = 0;
      fortran_vars.HF[2] = HF[2] = 0;
      fortran_vars.HF_fac[0] = HF_fac[0] = 0.25f;
      fortran_vars.HF_fac[1] = HF_fac[1] = 0.0f;
      fortran_vars.HF_fac[2] = HF_fac[2] = 0.0f;
      fortran_vars.screen = *screen = -1.0f;
   } else {
      cout << "In order to use external Functional you need to recompile ";
      cout << "LIO with libxc=1 or 2" << endl;
      fflush(stdout);
      exit(-1);
   }
}
#endif
