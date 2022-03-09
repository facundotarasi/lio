#include <iostream>
#include <omp.h>

#include <stdio.h>
#include <string.h> 

#include "../common.h"
#include "../init.h"
#include "../partition.h"

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_saverho_()
{
   cout << " Saving density and derivatives of Ground State" << endl;
   partition.lr_init();
   fflush(stdout); // NOT BUFFERED
}

//######################################################################
//######################################################################

namespace G2G {

void Partition::lr_init()
{

#pragma omp parallel for schedule(static)
    for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind - cubes.size()]->lr_closed_init();
         } else {
           cubes[ind]->lr_closed_init();
         }
      }
   }
   fflush(stdout);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

template<class scalar_type> void PointGroupCPU<scalar_type>::
               lr_closed_init()
{
   // This routine is temporary; we dont need save the GS density 
   // in memory any more.
}
//######################################################################
//######################################################################
#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
