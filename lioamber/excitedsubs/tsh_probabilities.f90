subroutine tsh_probabilities(C,E,Xexc,Eexc,NCO,M,Mlr,Ndim,Nvirt,Etot,Nstat)
use garcha_mod  , only: natom, Pmat_vec, nucvel, atom_mass
use excited_data, only: TSH, root, gamma_old, nucvel_old, Nesup, Nesup_old
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, Nstat
   LIODBLE, intent(in) :: C(M,Mlr), E(Mlr)
   LIODBLE, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   LIODBLE, intent(inout) :: Etot

   integer :: ii, jj
   LIODBLE :: Knr
   LIODBLE, allocatable :: Zvec(:), Xvec(:)
   LIODBLE, allocatable :: Xmat(:,:), Zmat(:,:), Zsym(:,:)
   LIODBLE, allocatable :: gammaWS(:,:), gammaXC(:,:), gammaCou(:,:)
   LIODBLE, allocatable :: gammaH(:,:), gammaT(:,:), gammaTot(:,:)
   LIODBLE, allocatable :: rhoG(:,:)

   if ( .not. TSH ) return
   if ( root == 0 ) return

!TODO: for the moment, this is the best place to put Energy
   Nesup(1) = Etot ! GS energy
   Etot = Etot + Eexc(root)
   Nesup(2) = Etot ! ES energy

   if ( root > Nstat ) then
      print*, "The root variable is bigger than nstates"
      print*, "Please set root <= nstates"
      stop
   endif
   
   print*, "TSH probabilities"
   allocate(Zvec(Ndim),Xvec(Ndim))
   Xvec = Xexc(:,root) / dsqrt(2.0d0)
   Zvec = Xvec / Eexc(root)

   ! Form Matrix in AO.
   allocate(Zmat(M,M),Xmat(M,M),gammaWS(natom,3))
   call VecToMat(Xvec,Xmat,C,Ndim,NCO,M,Mlr)
   call VecToMat(Zvec,Zmat,C,Ndim,NCO,M,Mlr)

   ! Obtain gamma WS: this calculates -gammaWS
   call gammaWS_calc(Xvec,Zvec,Zmat,E,C,gammaWS,NCO,M,Mlr,Ndim,Nvirt,&
                     natom)

   ! Obtain gamma Core
   allocate(Zsym(M,M),gammaH(natom,3))
   Zsym = Zmat + transpose(Zmat)
   call HVgradcalc(Zsym,gammaH,M,natom,.false.)

   ! Obtain gamma Coulomb
   allocate(rhoG(M,M)); call spunpack_rho('L',M,Pmat_vec,rhoG)
   allocate(gammaCou(natom,3)); gammaCou = 0.0d0
   call g2g_calcgammcou(rhoG,Zsym,gammaCou)
   deallocate(Zsym,rhoG)
   
   ! Obtain gamma XC
   allocate(Zsym(M,M),gammaXC(3,natom)); Zsym = 0.0d0; gammaXC = 0.0d0
   call g2g_calcgradxc(Zmat,Zsym,gammaXC,1)
   deallocate(Zsym)

   ! Obtain gamma T
   allocate(gammaT(natom,3)); gammaT = 0.0d0
   call intSG_Exc(gammaT,Xmat,natom,M)

   ! Obtain Total gamma
   allocate(gammaTot(natom,3))
   gammaTot = gammaWS + gammaH + gammaCou + transpose(gammaXC) - gammaT
!  print*, "gamma Tot"
!  do ii=1,natom
!     print*, ii,gammaTot(ii,1),gammaTot(ii,2),gammaTot(ii,3)
!  enddo
   deallocate(gammaWS,gammaH,gammaCou,gammaXC,gammaT)
   deallocate(Zvec,Xvec,Zmat,Xmat)

   ! Norm of NACMEs: When you perform Dynamic only.
   if (allocated(gamma_old)) then
      Knr = 0.0d0
      do ii=1,natom
      do jj=1,3
         Knr = Knr + gammaTot(ii,jj) / atom_mass(ii)
      enddo
      enddo
      print*, "NORM of NACVs", Knr*Knr*4.0d0

      call coef_propagator(gammaTot,nucvel,Nesup,gamma_old,nucvel_old,Nesup_old,natom)

      ! Save variables
      gamma_old = gammaTot
      nucvel_old = nucvel
      Nesup_old = Nesup
   endif
   
   deallocate(gammaTot)
end subroutine tsh_probabilities


subroutine coef_propagator(g,v,e,gOld,vOld,eOld,natom)
use excited_data, only: dE_accum, lambda, tsh_time_dt, B_old, &
                        tsh_Jstate, tsh_Kstate, tsh_coef,     &
                        excited_forces, CI_found
use fstsh_data, only: coef_Stat, dot_Stat, elec_Coup, elec_Ene, elec_Pha, current_state, tsh_minprob
use fstshsubs,  only: dot_calculation,interpol_EneCoupVel,phase_calculation,coef_evolution,correct_decoherence
use garcha_mod, only: atom_mass, npas
   implicit none

   integer, intent(in) :: natom
   LIODBLE, intent(in) :: g(natom,3), v(3,natom), e(2)
   LIODBLE, intent(in) :: gOld(natom,3), vOld(3,natom), eOld(2)

   integer :: ii, iter, tsh_nucStep, state_before, surf_old, surf_new
   LIODBLE :: Q, Qold, dt_elec, tot_time
   LIODBLE :: coup(2,2), coup_old(2,2)
   LIODBLE, allocatable :: elec_vel(:,:)
   TDCOMPLEX :: cj, ck, tmpc
   LIODBLE   :: norm, number_random, tmpr, cj2, prob, probFinal

   ! Hardcoded variables
   integer :: tsh_Enstep = 20

!  integer :: ii
!  LIODBLE :: Q, Gprob, factor, pop, temp1, temp2
!  LIODBLE :: number_random, kin_e, mass, norm
!  complex(kind=8) :: B, B_tot, B_abs, zero, B1, B2, pot, c_j, c_k
!  complex(kind=8), allocatable :: Uprop(:,:)

   state_before = tsh_Jstate
   current_state = tsh_Jstate

   if ( npas == 1 ) then
      ! Obtain Cdot at actual time
      !   this routine only needs variables at actual electronic time step
      call dot_calculation(coef_Stat, elec_Coup, elec_Pha, dot_Stat, 2)

      ! Save electronic variables
      coef_Stat(3,:) = coef_Stat(2,:); coef_Stat(2,:) = coef_Stat(1,:)
      elec_Coup(3,:,:) = elec_Coup(2,:,:); elec_Coup(2,:,:) = elec_Coup(1,:,:)
      elec_Pha(3,:,:) = elec_Pha(2,:,:); elec_Pha(2,:,:) = elec_Pha(1,:,:)
      dot_stat(3,:) = dot_stat(2,:); dot_stat(2,:) = dot_stat(1,:)
     return
   endif

   !  v = Nuclear Velocity [Bohr/au]
   !  g = Non-Adiabatic Coupling Vector [1/Bohr]
   !  Q = v X h [1/au]
   Q = 0.0d0; Qold = 0.d0
   coup = 0.d0; coup_old = 0.d0
   do ii=1,natom
     Q    = Q    + v(1,ii) * g(ii,1)        + v(2,ii) * g(ii,2)        + v(3,ii) * g(ii,3)
     Qold = Qold + vOld(1,ii) * gOld(ii,1)  + vOld(2,ii) * gOld(ii,2)  + vOld(3,ii) * gOld(ii,3)
   enddo
   coup(1,2) = Q; coup(2,1) = -Q
   coup_old(1,2) = Qold; coup_old(2,1) = -Qold
   if (CI_found) then
      coup = 0.0d0
      coup_old = 0.0d0
   endif

   allocate(elec_vel(3,natom))
   elec_vel = 0.d0
   probFinal = 0.0d0

   ! Time and Step Variables
   ! tsh_time_dt: Nuclear Time step in a.u
   ! tsh_nucStep: Nuclear step
   ! dt_elec    : Electronic Time step in a.u
   ! tsh_Enstep : Number of Electronic steps in each Nuclear step
   !
   dt_elec = tsh_time_dt / real(tsh_Enstep)
   tsh_nucStep = npas - 1

   ! Main Electronic Interpolation
   do iter = 0, tsh_Enstep-1
      tot_time = (dt_elec * iter + tsh_nucStep * tsh_time_dt) * 0.02418884254d0 ! change a.u to femto                                             
      write(*,"(3X,A,I2,A,F8.4,A)") "Electronic Sub-step= ",iter," Time= ", tot_time, " fs."

      ! Energy, Coupling and Velocities Interpolation
      call interpol_EneCoupVel(e,eOld,coup,coup_old,elec_Ene,elec_Coup,2,v,vOld,elec_vel, &
                               natom,tsh_time_dt,dt_elec,iter)

      ! Phase Calculation
      call phase_calculation(elec_Ene,elec_Pha,dt_elec,2,tsh_nucStep)
      
      ! Coefficients Evolution
      call coef_evolution(coef_Stat,elec_Coup,elec_Pha,dot_Stat,2,dt_elec,tsh_nucStep)
      
      ! Probabilities of Hopp Calculates
      cj   = coef_Stat(1,tsh_Jstate)
      tmpc = exp(cmplx(0.0d0,elec_Pha(1,tsh_Kstate,tsh_Jstate),COMPLEX_SIZE/2))
      ck   = conjg(coef_Stat(1,tsh_Kstate))
      tmpr = real(cj*ck*tmpc)
      prob = tmpr*elec_Coup(1,tsh_Kstate,tsh_Jstate)*(-2.0d0)
      norm = abs(coef_Stat(1,tsh_Kstate))**2
      if ( prob < 0.0d0 ) prob = 0.0d0
      cj2 = abs(cj)**2
      probFinal = probFinal + prob * dt_elec/cj2
      call random_number(number_random)
      print*, "probability, random", probFinal, number_random
      print*, "poblacion1", abs(coef_Stat(1,1))**2
      print*, "poblacion2", abs(coef_Stat(1,2))**2
      if ( probFinal > number_random .and. probFinal > tsh_minprob .and. (.not. CI_found) ) then
              write(*,"(4X,A,I2,A,I2)") "HOPP= ", tsh_Jstate, " -> ", tsh_Kstate
              current_state = tsh_Kstate
              CI_found = .true.
              tsh_Jstate = 1
              tsh_Kstate = 2
      endif
      if ( elec_Ene(1,tsh_Jstate) < elec_Ene(1,1) .and. tsh_Jstate /= 1 .and. (.not. CI_found) ) then
              write(*,"(4X,A)") "Forcing the system at Ground State"
              current_state = tsh_Kstate
              CI_found = .true.
              tsh_Jstate = 1
              tsh_Kstate = 2
      endif

      ! Decoherence Term
      call correct_decoherence(coef_Stat,elec_Ene,elec_vel,natom,2,dt_elec)

      ! Save Variables
      coef_Stat(3,:) = coef_Stat(2,:); coef_Stat(2,:) = coef_Stat(1,:)
      elec_Ene(3,:) = elec_Ene(2,:); elec_Ene(2,:) = elec_Ene(1,:)
      dot_stat(3,:) = dot_stat(2,:); dot_stat(2,:) = dot_stat(1,:)
      elec_Coup(3,:,:) = elec_Coup(2,:,:); elec_Coup(2,:,:) = elec_Coup(1,:,:)
      elec_Pha(3,:,:) = elec_Pha(2,:,:); elec_Pha(2,:,:) = elec_Pha(1,:,:)
   enddo ! ENDDO interpolation
   deallocate(elec_vel)

   !  Add this in order to change PES to GS
   if ( (state_before /= current_state) ) then
      print*, "HOPP After propagation"
      CI_found = .True.
      tsh_Jstate = 1
      tsh_Kstate = 2
      excited_forces = .False.
      dot_Stat  = cmplx(0.0d0,0.0d0,COMPLEX_SIZE/2)
      elec_Pha  = 0.0d0; elec_Coup = 0.0d0; elec_Ene = 0.0d0
   endif
end subroutine coef_propagator


subroutine gammaWS_calc(Xv,Zv,Zm,E,C,gamm,NCO,M,Mlr,Ndim,Nvirt,natom)
use garcha_mod  , only: ntatom, r, d
use faint_cpu   , only: intSG
use excited_data, only: fittExcited, Cocc, Cocc_trans, Coef_trans
use extern_functional_data, only: need_HF
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, natom
   LIODBLE, intent(in) :: Xv(Ndim), Zv(Ndim), Zm(M,M)
   LIODBLE, intent(in) :: C(M,Mlr), E(Mlr)
   LIODBLE, intent(out) :: gamm(natom,3)

   integer :: ii, jj, NCOc, pos1, MM, ind
   LIODBLE, allocatable :: F2e(:,:), Fxc(:,:), Ftot(:,:)
   LIODBLE, allocatable :: HZIJ(:,:), scratch(:,:)
   LIODBLE, allocatable :: Wmat(:,:), WmatMO(:,:), Wtot(:)

!  CALCULATE TWO ELECTRON PART
   allocate(F2e(M,M))
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(Zm,F2e,1)
      F2e = (F2e+transpose(F2e))
      call g2g_timer_stop("Fock 2e LR")
   elseif ( fittExcited .and. (.not. need_HF) ) then
      call g2g_timer_start("Fock 2e LR")
      call calc2eFITT(Zm,F2e,M)
      call g2g_timer_stop("Fock 2e LR")
   else
      print*, "Error in 2 Electron Repulsion Integrals"
      print*, "Check HF in the functional and fittExcited"
      stop
   endif

!  CALCULATE XC PART
   allocate(Fxc(M,M)); Fxc = 0.0d0
   call g2g_calculateXC(Zm,Fxc,2)

!  TOTAL FOCK
   allocate(Ftot(M,M)); Ftot = F2e + Fxc + Fxc
   deallocate(F2e,Fxc)

!  CHANGE BASIS of FOCK. AO -> MO(OCC x OCC)
   allocate(scratch(M,NCO),HZIJ(NCO,NCO))
   call dgemm('N','N',M,NCO,M,1.0d0,Ftot,M,Cocc,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,NCO,M,1.0d0,Cocc_trans,NCO,scratch,M, &
              0.0d0,HZIJ,NCO)
   deallocate(scratch,Ftot)

   allocate(WmatMO(Mlr,Mlr)); WmatMO = 0.0d0
!  FORM OCC x OCC Block
   NCOc = NCO + 1
   do ii=1,NCO
   do jj=1,ii
      WmatMO(NCOc-ii,NCOc-jj) = HZIJ(NCOc-ii,NCOc-jj)
   enddo
   enddo
   do ii=1,NCO
      WmatMO(ii,ii) = WmatMO(ii,ii) * 0.5d0
   enddo
   deallocate(HZIJ)

!  FORM OCC x VIRT Block
   do ii=1,NCO
   do jj=1,Nvirt
      pos1 = (ii-1) * Nvirt + jj
      WmatMO(NCOc-ii,NCO+jj) = E(NCOc-ii) * Zv(pos1) + Xv(pos1)
   enddo
   enddo

   ! CHANGE BASIS of Wmat. MO -> AO
   allocate(scratch(M,Mlr),Wmat(M,M))
   call dgemm('N','N',M,Mlr,Mlr,1.0d0,C,M,WmatMO,Mlr,0.0d0,scratch,M)
   call dgemm('N','N',M,M,Mlr,1.0d0,scratch,M,Coef_trans,Mlr,0.0d0, &
              Wmat,M)
   deallocate(scratch,WmatMO)

   ! NACVs
   MM = M * (M + 1) / 2
   allocate(Wtot(MM))
   ind = 1
   do ii=1,M
      Wtot(ind) = Wmat(ii,ii)
      ind = ind + 1
      do jj=ii+1,M
         Wtot(ind) = Wmat(ii,jj) + Wmat(jj,ii)
         ind = ind + 1
      enddo
   enddo
   Wtot = (-1.0d0) * Wtot
   gamm = 0.0d0
   call intSG(gamm, Wtot, r, d, natom, ntatom)
   deallocate(Wtot,Wmat)
   gamm = 2.0d0 * gamm
end subroutine gammaWS_calc

