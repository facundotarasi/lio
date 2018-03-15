!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc2_FockMao( FockMao, Energy )
!
!  Time is in ps?fs?
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use maskrmm,     only: rmmput_fock, rmmget_fock

   use faint_cpu77, only: intsol, int2, int3mem

   use garcha_mod,  only: M, Md, RMM, kkind, kkinds, cool, cools, igrid2, MEMO

   implicit none
   real*8,intent(out)    :: FockMao(M,M)
   real*8,intent(out)    :: Energy

   real*8   :: Energy_1e
   real*8   :: Energy_Coulomb
   real*8   :: Energy_SolvT,Energy_SolvF

   integer  :: kk, idx0
   integer  :: MM, MMd, igpu
!
!
!  Initializations
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-init')
   call rmmput_fock(FockMao)

   if (allocated(kkind))  deallocate(kkind)
   if (allocated(kkinds)) deallocate(kkinds)
   if (allocated(cool))   deallocate(cool)
   if (allocated(cools))  deallocate(cools)

   call g2g_reload_atom_positions(igrid2)
   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()
   call g2g_timer_stop('RMMcalc2-init')
!
!
! Calculate fixed-parts of fock
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-sol2coul')
   if (igpu.le.1) then
      call intsol(Energy_SolvF,Energy_SolvT,.true.)
   else
      call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
   endif

   call int2()
   if (igpu.gt.2) call aint_coulomb_init()
   if (igpu.eq.5) MEMO = .false.
   call g2g_timer_stop('RMMcalc2-sol2coul')

   if (MEMO) then
      call g2g_timer_start('RMMcalc2-int3mem')
      call int3mem()
      call g2g_timer_stop('RMMcalc2-int3mem')
   endif
!
!
!  Prepare Outputs
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-exit')
   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2
   idx0=3*MM+2*MMd
   Energy_1e=0.0d0
   do kk=1,MM
      Energy_1e=Energy_1e+RMM(kk)*RMM(idx0+kk)
   enddo

!  Energy=0.0d0
   Energy=Energy+Energy_1e
   Energy=Energy+Energy_Coulomb
   Energy=Energy+Energy_SolvT

   call rmmget_fock(FockMao)
   call g2g_timer_stop('RMMcalc2-exit')

end subroutine RMMcalc2_FockMao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!