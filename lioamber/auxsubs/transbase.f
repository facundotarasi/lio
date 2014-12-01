!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! BASETRANSFORM PROCEDURES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function transbase_d(M,Utrp,Mati,Umat) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Utrp(M,M),Umat(M,M)
       real*8,intent(in)      :: Mati(M,M)
       real*8,allocatable     :: Matm(:,:)
       real*8,allocatable     :: Mato(:,:)
       integer                :: kii,kki,kkd,kdd
!
!
! Initialization of Matm(nnd,ndd) and Mato(nii,ndd)
!--------------------------------------------------------------------!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DBLE(0)
       Mato=DBLE(0)
!
!
! First Product Mati(nni,nnd)*Umat(nnd,ndd)
!--------------------------------------------------------------------!
       do kdd=1,M
       do kkd=1,M
       do kki=1,M
!         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Umat(kkd,kdd)
         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Utrp(kdd,kkd)
       enddo
       enddo
       enddo
!
! Second Product Utrp(nii,nni)*Matm(nni,ndd)
!--------------------------------------------------------------------!
       do kdd=1,M
       do kki=1,M
       do kii=1,M
         Mato(kii,kdd)=Mato(kii,kdd)+Utrp(kii,kki)*Matm(kki,kdd)
       enddo
       enddo
       enddo
!
       return;end function
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function transbase_z(M,Utrp,Mati,Umat) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Utrp(M,M),Umat(M,M)
       complex*16,intent(in)  :: Mati(M,M)
       complex*16,allocatable :: Matm(:,:)
       complex*16,allocatable :: Mato(:,:)
       integer                :: kii,kki,kkd,kdd
!
!
! Initialization of Matm(nnd,ndd) and Mato(nii,ndd)
!--------------------------------------------------------------------!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DCMPLX(0,0)
       Mato=DCMPLX(0,0)
!
!
! First Product Mati(nni,nnd)*Umat(nnd,ndd)
!--------------------------------------------------------------------!
       do kdd=1,M
       do kkd=1,M
       do kki=1,M
!         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Umat(kkd,kdd)
         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Utrp(kdd,kkd)
       enddo
       enddo
       enddo
!
! Second Product Utrp(nii,nni)*Matm(nni,ndd)
!--------------------------------------------------------------------!
       do kdd=1,M
       do kki=1,M
       do kii=1,M
         Mato(kii,kdd)=Mato(kii,kdd)+Utrp(kii,kki)*Matm(kki,kdd)
       enddo
       enddo
       enddo
!
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
