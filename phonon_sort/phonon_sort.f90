program phonon_sort
  implicit none
  integer :: iq, jq, kq, ig, jg, kg, nq1, nq2, nq3, density, ngws, ibrav, nq
  real(kind=8) :: ql, dql, qx0, qy0, qz0, weight
  real(kind=8) :: celldm(6), a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), &
                  qpoint(3), q_cart(3), q_crys(3), gws(0:3,200)
  namelist /q_grid/ nq1, nq2, nq3, dql, density
  !
  !default values for the necessary variables 
  nq1 = 60         ! q numbers for each of three axes in the fine grid
  nq2 = 60
  nq3 = 1         
  dql = 0.002      ! starting distance away from Gamma point
  density = 500    ! number density of q-points for catching the branch continuity
  !
  ! read variables from input file "phonon_sort.in"
  read (*, nml = q_grid)
  !
  open (10, file='ifc', status='old')
  read(10,*) iq, jq, ibrav, celldm
  close(10)
  ! a1, a2, a3 in alat unit; b1, b2, b3 in tpi/alat unit
  call latgen_lib(ibrav, celldm, a1, a2, a3, b1, b2, b3)
  !
  call wsinit(gws,ngws,b1, b2, b3)
  !
  do iq = 1, nq1
    do jq = 1, nq2
      do kq = 1, nq3
        qpoint = b1*(iq-1)/nq1 + b2*(jq-1)/nq2 + b3*(kq-1)/nq3
        !
        do ig = -3, 3
          do jg = -3, 3
            do kg = -3, 3
              q_cart = qpoint + ig*b1 + jg*b2 + kg*b3
              call wsweight(q_cart,gws,ngws,weight)
              if(weight > 0) goto 100
            end do
          end do
        end do
        ! convert to crystal frame
100     q_crys(1) = dot_product(a1,q_cart)
        q_crys(2) = dot_product(a2,q_cart)
        q_crys(3) = dot_product(a3,q_cart)
        ql = sqrt(q_cart(1)**2 + q_cart(2)**2 + q_cart(3)**2)
        nq = nint(density*ql)
        if (ql <= dql) then
          nq = 1
          qx0 = q_crys(1)
          qy0 = q_crys(2)
          qz0 = q_crys(3)
        else
          qx0 = dql*q_crys(1)/ql
          qy0 = dql*q_crys(2)/ql
          qz0 = dql*q_crys(3)/ql
        end if
        !
        open (10, file='matdyn.in', status='replace')
        write (10, '(a)') '&input'
        write (10, '(a)') "flfrc = 'ifc'"
        write (10, '(a)') "asr = 'all'"
        write (10, '(a)') 'loto_2d = .true.'
        write (10, '(a)') 'read_lr = .true.'
        write (10, '(a)') 'q_in_band_form = .true.'
        write (10, '(a)') 'q_in_cryst_coord = .true.'
        write (10, '(a)') 'eigen_similarity = .true.'
        write (10, '(a)') 'lphsort = .true.'
        write (10, '(a)') '/'
        write (10, '(i1)') 2
        write (10, '(3f20.15,i8)') qx0, qy0, qz0, nq
        write (10, '(3f20.15,i8)') q_crys(1), q_crys(2), q_crys(3), 1
        close(10)
        !
        call system('./matdyn.x <matdyn.in> matdyn.out')
        !
      end do
    end do
  end do
  !
  call system('rm -rf matdyn.*')
  !
  stop  
end program
!    
!    
SUBROUTINE latgen_lib(ibrav,celldm,a1,a2,a3,b1,b2,b3)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  real(kind=8), INTENT(inout) :: celldm(6)
  real(kind=8), INTENT(inout) :: a1(3), a2(3), a3(3)
  real(kind=8), INTENT(out) :: b1(3), b2(3), b3(3)
  !
  character(len=100) :: errormsg
  INTEGER :: ierr,ir
  real(kind=8) :: term, cbya, tpi, sr2, sr3, term1, term2, singam, sen, omega, alat
  !
  tpi = 2.d0*acos(-1.d0)
  sr2 = sqrt(2.d0)
  sr3 = sqrt(3.d0)
  !
  ! pre-set to zero, in case we quit because of error
  !
  !  user-supplied lattice vectors
  !
  IF (ibrav == 0) THEN
     IF (sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  &
         THEN; errormsg='wrong at for ibrav=0'; ierr= 1; RETURN; ENDIF
     IF (sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  &
         THEN; errormsg='wrong at for ibrav=0'; ierr= 2; RETURN; ENDIF
     IF (sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  &
         THEN; errormsg='wrong at for ibrav=0'; ierr= 3; RETURN; ENDIF

     IF ( celldm(1) /= 0.D0 ) THEN
     !
     ! ... input at are in units of alat => convert them to a.u.
     !
         a1(:) = a1(:) * celldm(1)
         a2(:) = a2(:) * celldm(1)
         a3(:) = a3(:) * celldm(1)
     ELSE
     !
     ! ... input at are in atomic units: define celldm(1) from a1
     !
         celldm(1) = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
     ENDIF
     !
  ELSE
     a1(:) = 0.d0
     a2(:) = 0.d0
     a3(:) = 0.d0
  ENDIF
  !
  IF (celldm (1) <= 0.d0) THEN; errormsg='wrong celldm(1)'; ierr= abs(ibrav); RETURN; ENDIF
  !
  !  index of bravais lattice supplied
  !
  IF (ibrav == 1) THEN
     !
     !     simple cubic lattice
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  ELSEIF (ibrav == 2) THEN
     !
     !     fcc lattice
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  ELSEIF (abs(ibrav) == 3) THEN
     !
     !     bcc lattice
     !
     term=celldm(1)/2.d0
     DO ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     ENDDO
     IF ( ibrav < 0 ) THEN
        a1(1)=-a1(1)
        a2(2)=-a2(2)
        a3(3)=-a3(3)
     ELSE
        a2(1)=-a2(1)
        a3(1)=-a3(1)
        a3(2)=-a3(2)
     ENDIF
     !
  ELSEIF (ibrav == 4) THEN
     !
     !     hexagonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr= ibrav; RETURN; ENDIF
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (abs(ibrav) == 5) THEN
     !
     !     trigonal lattice
     !
     IF (celldm (4) <= -0.5d0 .or. celldm (4) >= 1.0d0) &
             THEN; errormsg='wrong celldm(4)'; ierr=5; RETURN; ENDIF
     !
     term1=sqrt(1.0d0 + 2.0d0*celldm(4))
     term2=sqrt(1.0d0 - celldm(4))
     !
     IF ( ibrav == 5) THEN
        !     threefold axis along c (001)
        a2(2)=sr2*celldm(1)*term2/sr3
        a2(3)=celldm(1)*term1/sr3
        a1(1)=celldm(1)*term2/sr2
        a1(2)=-a1(1)/sr3
        a1(3)= a2(3)
        a3(1)=-a1(1)
        a3(2)= a1(2)
        a3(3)= a2(3)
     ELSEIF ( ibrav == -5) THEN
        !     threefold axis along (111)
        ! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
        ! does not yield the x,y,z axis, but an equivalent rotated triplet:
        !    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
        ! If you prefer the x,y,z axis as cubic limit, you should modify the
        ! definitions of a1(1) and a1(2) as follows:'
        !    a1(1) = celldm(1)*(term1+2.0d0*term2)/3.0d0
        !    a1(2) = celldm(1)*(term1-term2)/3.0d0
        ! (info by G. Pizzi and A. Cepellotti)
        !
        a1(1) = celldm(1)*(term1-2.0d0*term2)/3.0d0
        a1(2) = celldm(1)*(term1+term2)/3.0d0
        a1(3) = a1(2)
        a2(1) = a1(3)
        a2(2) = a1(1)
        a2(3) = a1(2)
        a3(1) = a1(2)
        a3(2) = a1(3)
        a3(3) = a1(1)
     ENDIF
  ELSEIF (ibrav == 6) THEN
     !
     !     tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=6; RETURN; ENDIF
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (ibrav == 7) THEN
     !
     !     body centered tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=7; RETURN; ENDIF
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  ELSEIF (ibrav == 8) THEN
     !
     !     Simple orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; errormsg='wrong celldm(2)'; ierr=8; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=8; RETURN; ENDIF
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF ( abs(ibrav) == 9) THEN
     !
     !     One face (base) centered orthorhombic lattice  (C type)
     !
     IF (celldm (2) <= 0.d0) THEN; errormsg='wrong celldm(2)'; ierr=9; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=9; RETURN; ENDIF
     !
     IF ( ibrav == 9 ) THEN
        !   old PWscf description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) = a1(1) * celldm(2)
        a2(1) = - a1(1)
        a2(2) = a1(2)
     ELSE
        !   alternate description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) =-a1(1) * celldm(2)
        a2(1) = a1(1)
        a2(2) =-a1(2)
     ENDIF
     a3(3) = celldm(1) * celldm(3)
     !
  ELSEIF ( ibrav == 91 ) THEN
     !
     !     One face (base) centered orthorhombic lattice  (A type)
     !
     IF (celldm (2) <= 0.d0) THEN; errormsg='wrong celldm(2)'; ierr=91; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=91; RETURN; ENDIF
     !
     a1(1) = celldm(1)
     a2(2) = celldm(1) * celldm(2) * 0.5d0
     a2(3) = - celldm(1) * celldm(3) * 0.5d0
     a3(2) = a2(2)
     a3(3) = - a2(3)
     !
  ELSEIF (ibrav == 10) THEN
     !
     !     All face centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; errormsg='wrong celldm(2)'; ierr=10; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=10; RETURN; ENDIF
     !
     a2(1) = 0.5d0 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 11) THEN
     !
     !     Body centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) THEN; errormsg='wrong celldm(2)'; ierr=11; RETURN; ENDIF
     IF (celldm (3) <= 0.d0) THEN; errormsg='wrong celldm(3)'; ierr=11; RETURN; ENDIF
     !
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 12) THEN
     !
     !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
     !
     IF (celldm (2) <= 0.d0)   THEN; errormsg='wrong celldm(2)'; ierr=12; RETURN; ENDIF
     IF (celldm (3) <= 0.d0)   THEN; errormsg='wrong celldm(3)'; ierr=12; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; errormsg='wrong celldm(4)'; ierr=12; RETURN; ENDIF
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF (ibrav ==-12) THEN
     !
     !     Simple monoclinic lattice, unique axis: b (more common)
     !
     IF (celldm (2) <= 0.d0)   THEN; errormsg='wrong celldm(2)'; ierr=12; RETURN; ENDIF
     IF (celldm (3) <= 0.d0)   THEN; errormsg='wrong celldm(3)'; ierr=12; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; errormsg='wrong celldm(5)'; ierr=12; RETURN; ENDIF
     !
     sen=sqrt(1.d0-celldm(5)**2)
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(3)=celldm(1)*celldm(3)*sen
     !
  ELSEIF (ibrav == 13) THEN
     !
     !     One face centered monoclinic lattice unique axis c
     !
     IF (celldm (2) <= 0.d0)   THEN; errormsg='wrong celldm(2)'; ierr=13; RETURN; ENDIF
     IF (celldm (3) <= 0.d0)   THEN; errormsg='wrong celldm(3)'; ierr=13; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; errormsg='wrong celldm(4)'; ierr=13; RETURN; ENDIF
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(3) =-a1(1) * celldm(3)
     a2(1) = celldm(1) * celldm(2) * celldm(4)
     a2(2) = celldm(1) * celldm(2) * sen
     a3(1) = a1(1)
     a3(3) =-a1(3)
  ELSEIF (ibrav == -13) THEN
     errormsg='BEWARE: axis for ibrav=-13 changed, see documentation!'
     !
     !     One face centered monoclinic lattice unique axis b
     !
     IF (celldm (2) <= 0.d0)   THEN; errormsg='wrong celldm(2)'; ierr=13; RETURN; ENDIF
     IF (celldm (3) <= 0.d0)   THEN; errormsg='wrong celldm(3)'; ierr=13; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; errormsg='wrong celldm(5)'; ierr=13; RETURN; ENDIF
     !
     sen = sqrt( 1.d0 - celldm(5) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a2(1) =-a1(1)
     a2(2) = a1(2)
     a3(1) = celldm(1) * celldm(3) * celldm(5)
     a3(3) = celldm(1) * celldm(3) * sen
     !
  ELSEIF (ibrav == 14) THEN
     !
     !     Triclinic lattice
     !
     IF (celldm (2) <= 0.d0)   THEN; errormsg='wrong celldm(2)'; ierr=14; RETURN; ENDIF
     IF (celldm (3) <= 0.d0)   THEN; errormsg='wrong celldm(3)'; ierr=14; RETURN; ENDIF
     IF (abs(celldm(4))>=1.d0) THEN; errormsg='wrong celldm(4)'; ierr=14; RETURN; ENDIF
     IF (abs(celldm(5))>=1.d0) THEN; errormsg='wrong celldm(5)'; ierr=14; RETURN; ENDIF
     IF (abs(celldm(6))>=1.d0) THEN; errormsg='wrong celldm(6)'; ierr=14; RETURN; ENDIF
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)
     IF (term < 0.d0) THEN; errormsg='celldm do not make sense, check your data'; ierr=14; RETURN; ENDIF
     term= sqrt(term/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  ELSE IF (ibrav /= 0) THEN
     !
     errormsg='nonexistent bravais lattice'
     ierr=ABS(ibrav)
     RETURN
     !
  ENDIF
  !
  !  calculate unit-cell volume omega
  !
  omega = a1(1)*(a2(2)*a3(3)-a2(3)*a3(2)) - &
         a1(2)*(a2(1)*a3(3)-a2(3)*a3(1)) + &
         a1(3)*(a2(1)*a3(2)-a2(2)*a3(1))
  !
  ! the lattice parameter 
  alat = celldm(1)
  !
  !  Generate the reciprocal lattice vectors from the direct unit cell       
 b1(1) = tpi*(a2(2)*a3(3)-a3(2)*a2(3))/omega
 b1(2) = tpi*(a3(1)*a2(3)-a2(1)*a3(3))/omega
 b1(3) = tpi*(a2(1)*a3(2)-a3(1)*a2(2))/omega
 b2(1) = tpi*(a3(2)*a1(3)-a1(2)*a3(3))/omega
 b2(2) = tpi*(a1(1)*a3(3)-a3(1)*a1(3))/omega
 b2(3) = tpi*(a3(1)*a1(2)-a1(1)*a3(2))/omega
 b3(1) = tpi*(a1(2)*a2(3)-a2(2)*a1(3))/omega
 b3(2) = tpi*(a2(1)*a1(3)-a1(1)*a2(3))/omega
 b3(3) = tpi*(a1(1)*a2(2)-a2(1)*a1(2))/omega
 
!  cart. coord. in units of alat and 2pi/alat
 a1 = a1/alat; a2 = a2/alat; a3 = a3/alat
 b1 = b1*alat/tpi; b2 = b2*alat/tpi; b3 = b3*alat/tpi
 !
  RETURN
  !
END SUBROUTINE latgen_lib
!
subroutine wsinit(rws,nrws,atw1,atw2,atw3)
!-----------------------------------------------------------------------
! find the points which surround the origin in a relatively safe range,
! the vertical planes through the midpoint of lines form the BZ boundaries
!

 implicit none

!Arguments -------------------------------
!scalars
 integer(kind=4),intent(out) :: nrws
!arrays
 real(kind=8),intent(in) :: atw1(3), atw2(3), atw3(3)
 real(kind=8),intent(out) :: rws(0:3,200)

!Local variables-------------------------------
!scalars
 integer(kind=4) :: i,ii,ir,jr,kr,nx
 real(kind=8) :: eps

! *********************************************************************
 
 rws = 0.d0
 eps = 1.d-6
 nx = 2
 ii = 1
 do ir=-nx,nx
   do jr=-nx,nx
     do kr=-nx,nx
       do i=1,3
         rws(i,ii) = atw1(i)*ir + atw2(i)*jr + atw3(i)*kr
       end do
       rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+rws(3,ii)*rws(3,ii)
       rws(0,ii)=0.5d0*rws(0,ii)
       if (rws(0,ii).gt.eps) ii = ii + 1
     end do
   end do
 end do
 nrws = ii - 1
 
 return
end subroutine wsinit
!
subroutine wsweight(r,rws,nrws,weight)
!-----------------------------------------------------------------------
! wsweights assigns this weight:
! - if a point is inside the Wigner-Seitz cell:    weight=1
! - if a point is outside the WS cell:             weight=0
! - if a point q is on the border of the WS cell, it finds the number N 
!   of translationally equivalent point q+G  (where G is a lattice vector)
!   that are also on the border of the cell. Then: weight = 1/N

! I.e. if a point is on the surface of the WS cell of a cubic lattice 
! it will have weight 1/2; on the vertex of the WS it would be 1/8; 
! the K point of an hexagonal lattice has weight 1/3 and so on.

 implicit none

!Arguments -------------------------------
!scalars
 integer(kind=4),intent(in) :: nrws
 real(kind=8),intent(out) :: weight
!arrays
 real(kind=8),intent(in) :: r(3),rws(0:3,nrws)

!Local variables-------------------------------
!scalars
 integer(kind=4) :: ir,nreq
 real(kind=8) :: rrt,ck,eps

! *********************************************************************
 eps=1.d-6
 weight = 0.d0
 nreq = 1
 do ir =1,nrws
   rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
   ck = rrt-rws(0,ir)
   if (ck .gt. eps) return
   if (abs(ck) .lt. eps) nreq = nreq + 1
 end do
 weight = 1.d0/dble(nreq)
 
 return
end subroutine wsweight