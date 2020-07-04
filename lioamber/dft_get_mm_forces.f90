#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates the QMMM forces in the QM and the MM regions.
subroutine dft_get_mm_forces(dxyzcl, dxyzqm)
   use garcha_mod, only: natom, ntatom, r, d, pc, Iz, nsol, Pmat_vec
   use basis_data, only: M
   use faint_cpu , only: int1G, intsolG
   use lj_switch , only: ljs_gradients_qmmm
   use excited_data,only: excited_forces, pack_dens_exc

   implicit none
   LIODBLE, intent(inout) :: dxyzqm(3, natom)
   LIODBLE, intent(inout) :: dxyzcl(3, nsol)

   LIODBLE, dimension (:,:), allocatable :: ff, ffcl
   integer :: MM, iatom, jatom, jcrd, igpu

   if (nsol .le. 0) return
   MM = M * (M +1) / 2

   call aint_query_gpu_level(igpu)
   call g2g_timer_sum_start('QM/MM gradients')
   if (igpu .lt. 2 .or. excited_forces) then
      ! The old version of intsolG expected the MM force array to be
      ! padded in front with # QM atoms spots for some reason
      allocate(ff(natom,3), ffcl(ntatom,3))
      ffcl = 0.0D0
      ff   = 0.0D0

      call g2g_timer_start('intsolG')
      if ( excited_forces ) then
         call intsolG(ff, ffcl, natom, ntatom, pack_dens_exc, d, r, pc, Iz)
      else
         call intsolG(ff, ffcl, natom, ntatom, Pmat_vec, d, r, pc, Iz)
      endif
      call g2g_timer_stop('intsolG')

      do jatom = 1, nsol
      do iatom = 1, 3
        dxyzcl(iatom, jatom) = ffcl(natom + jatom, iatom)
      enddo
      enddo
   else
      ! The GPU version of the QM/MM gradients only uses space for the MM
      ! forces in the MM force array
      allocate(ff(natom,3), ffcl(nsol,3))
      ffcl = 0.0D0
      ff   = 0.0D0

      if (igpu.gt.3) call int1G(ff, Pmat_vec, d, r, Iz, natom, ntatom, .true., .false.)
      call g2g_timer_start('aint_qmmm_forces')
      call aint_qmmm_forces(ff, ffcl)
      call g2g_timer_stop('aint_qmmm_forces')

      do jatom = 1, nsol
      do iatom = 1, 3
         dxyzcl(iatom,jatom) = ffcl(jatom,iatom)
      enddo
      enddo
   endif
 
   ! Accumulates the MM forces for the QM region.
   do iatom = 1, natom
   do jcrd = 1, 3
      dxyzqm(jcrd, iatom) = ff(iatom,jcrd) + dxyzqm(jcrd,iatom)
   enddo
   enddo

   call ljs_gradients_qmmm(dxyzqm, dxyzcl, r, natom)

   call g2g_timer_sum_stop('QM/MM gradients')
   call g2g_timer_sum_stop('Forces')

   deallocate (ff,ffcl)
end subroutine dft_get_mm_forces
