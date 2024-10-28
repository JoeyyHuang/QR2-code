PROGRAM Raman_postproc
!-----------------------------------------------------------------------
!!
!! postprocessing code for Raman intensity PLOT and Raman mode ANALYSIS
!!
!! by Jianqi Huang      Dec 2022
!-----------------------------------------------------------------------
  USE MPI
  IMPLICIT NONE
  CHARACTER(LEN = 256) :: dir_raman
  !! specify the directory where the Raman files are generated in the previous run
  CHARACTER(LEN = 6) :: Raman_type
  CHARACTER(LEN = 12) :: fRaman
  !! name of file with q-resolved Raman intensity
  CHARACTER(LEN = 12) :: fmodes
  !! name of file with Raman modes analysis
  CHARACTER(LEN = 2) :: str
  !! '#q' string
  INTEGER :: istep, nstep, iq, iiq, nq, nmodes, mu, nu, nRaman_modes, iRaman_modes, ipair, &
             my_id, num_procs, proc_index, residual_value, q_per_proc, ierr, ios, nphonon_modes
  INTEGER, ALLOCATABLE :: location(:, :)
  REAL(KIND = 8) :: time, twopi, Rs_min, Rs_max, Rs_inc, Lorentz_func, Lgamma, frq, frq2, &
                    Raman_p, Raman_m, Raman_e_p, Raman_e_m, Raman_h_p, Raman_h_m, &
                    Raman_pp, Raman_mm, Raman_pm, Raman_mp, &
                    Raman_ee_pp, Raman_ee_mm, Raman_ee_pm, Raman_ee_mp, &
                    Raman_hh_pp, Raman_hh_mm, Raman_hh_pm, Raman_hh_mp, &
                    Raman_eh_pp, Raman_eh_mm, Raman_eh_pm, Raman_eh_mp, &
                    Raman_he_pp, Raman_he_mm, Raman_he_pm, Raman_he_mp, &
                    Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                    Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m, &
                    Raman_intensity
  REAL(KIND = 8) :: Raman_modes(50), q(3)
  !! at most for 50 Raman modes analysis
  REAL(KIND = 8), ALLOCATABLE :: Raman_tot(:), Raman_tot0(:)
  !! total Raman intensity at master node or slave node
  REAL(KIND = 8), ALLOCATABLE :: Raman_e(:)
  !! Raman intensity only from el-ph in one-phonon case
  REAL(KIND = 8), ALLOCATABLE :: Raman_h(:)
  !! Raman intensity only from ho-ph in one-phonon case
  REAL(KIND = 8), ALLOCATABLE :: Raman_ee(:), Raman_ee0(:)
  !! Raman intensity only from el-ph in two-phonon case
  REAL(KIND = 8), ALLOCATABLE :: Raman_hh(:), Raman_hh0(:)
  !! Raman intensity only from ho-ph in one-phonon case
  REAL(KIND = 8), ALLOCATABLE :: Raman_eh(:), Raman_eh0(:)
  !! Raman intensity from el-ph and ho-ph in sequence in two-phonon case
  REAL(KIND = 8), ALLOCATABLE :: Raman_he(:), Raman_he0(:)
  !! Raman intensity from ho-ph and el-ph in sequence in two-phonon case
  REAL(KIND = 8), ALLOCATABLE :: modes_intensity(:, :, :), modes_intensity0(:, :, :)
  !! modes-resolved Raman intensity at master node or slave node in DRR case
  REAL(KIND = 8), ALLOCATABLE :: mode_intensity(:, :), mode_intensity0(:, :)
  !! modes-resolved Raman intensity at master node or slave node in defect-induced DRR case
  REAL(KIND = 8), ALLOCATABLE :: frq0(:)
  !! phonon frequency in branch-sorted order
  LOGICAL :: lRaman_modes
  !! if .ture., perform Raman modes analysis
  LOGICAL :: lhoph
  !! if .true., consider the hole-phonon coupling
  LOGICAL :: lRay_sca
  !! if .true., consider the Rayleigh scattering around 0 cm-1
  REAL(KIND = 8) :: Ray_thr
  !! the radius threshold around 0 cm-1 insider which regarded as Rayleigh scattering
  LOGICAL :: alive
  !! check if qraman file exist
  LOGICAL, ALLOCATABLE :: mask(:, :, :)
  !! for maxloc function
  !
  NAMELIST / PLOT / dir_raman, Raman_type, lhoph, Rs_min, Rs_max, Rs_inc, Lgamma, lRay_sca, Ray_thr
  NAMELIST / ANALYSIS / lRaman_modes, nphonon_modes, nRaman_modes, Raman_modes
  !
  ! default values for the necessary variables
  dir_raman = '../qraman'
  Raman_type = 'double'
  lhoph = .false.
  lRay_sca = .false.
  Ray_thr = 5.d0
  Rs_min = -4.d3
  Rs_max = 4.d3
  Rs_inc = 1.d-1
  Lgamma = 1.d1
  lRaman_modes = .false.
  nphonon_modes = 6
  nRaman_modes = 4
  Raman_modes(:) = 0.d0
  !
  twopi = 2.d0 * DACOS(-1.d0)
  !
  ! MPI initialization
  CALL MPI_INIT (ierr)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  !
  IF (my_id == 0) THEN
    ! read namelist PLOT from standard input file "raman.in"
    READ(*, NML = PLOT)
    ! number of steps in Raman shift
    nstep = NINT((Rs_max - Rs_min) / Rs_inc) + 1
    !
    INQUIRE(FILE = TRIM(dir_raman), EXIST = alive)
    IF(.NOT. alive) THEN
      WRITE(*, '("ERROR: dir_raman not exist")')
      STOP
    ENDIF
    ! how many fine q points
    IF (Raman_type == 'double' .OR. Raman_type == 'defect') THEN
      OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
      DO WHILE(.true.)
        READ(10, '(a2)', iostat = ios) str
        IF (ios /= 0) EXIT
        IF (str == '#q') THEN
          BACKSPACE(10)
          READ(10, '(4x,i5)') nq
        ENDIF
      ENDDO
      CLOSE(10)
    ENDIF
    !
    ! how many normal modes in phonon
    OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
    READ(10, *)
    READ(10, *)
    IF (Raman_type == 'single') THEN
      DO WHILE(.true.)
        READ(10, '(2x, i3)', IOSTAT = ios) nmodes
        IF (ios /= 0) EXIT
      ENDDO
    ELSEIF (Raman_type == 'double' .OR. Raman_type == 'defect') THEN
      DO WHILE(.true.)
        READ(10, '(a2)', iostat = ios) str
        IF (ios /= 0) THEN
          BACKSPACE(10)
          READ(10, '(7x, i3)') nmodes
          EXIT
        ELSEIF (str == '#q') THEN
          BACKSPACE(10)
          BACKSPACE(10)
          READ(10, '(7x, i3)') nmodes
          EXIT
        ENDIF
      ENDDO
    ENDIF
    CLOSE(10)
    !
    ! read namelist ANALYSIS from standard input file "raman.in"
    READ(*, NML = ANALYSIS)
    !
    ! number of Raman modes for analysis
    IF (nRaman_modes > 50) THEN
      WRITE(*, '("ERROR: too many Raman modes, should be less than 50")')
      STOP
    ENDIF
    !
  ENDIF
  !
  CALL MPI_Bcast(dir_raman, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Raman_type, 6, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(lhoph, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(nmodes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(nq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Rs_min, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Rs_max, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Rs_inc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(nstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Lgamma, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(lRay_sca, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Ray_thr, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(lRaman_modes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(nphonon_modes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(nRaman_modes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(Raman_modes, 50, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !
  !
  SELECT CASE (Raman_type)
  !
  !################################################################
  CASE ('single')   ! single resonance Raman 
  !################################################################
  !
  iq = 1
  !
  ! not consider hole-phonon scattering
  ! *********************************************************************
  IF (.NOT. lhoph) THEN
  ! *********************************************************************
    IF (my_id == 0) THEN
      !
      CALL CPU_TIME(time)
      WRITE(*, '("begin time:",f18.4," secs")') time
      WRITE(*, *)
      !
      ALLOCATE(Raman_tot(nstep))
      ! Raman data initialization
      Raman_tot(:) = 0.d0
      !
      OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
      READ(10, *)
      READ(10, *)
      DO mu = 1, nmodes
        READ(10, '(7x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
        ! calculate Raman data
        DO istep = 1, nstep
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc-frq)**2)
          Raman_tot(istep) = Raman_tot(istep) + Raman_p * Lorentz_func
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc+frq)**2)
          Raman_tot(istep) = Raman_tot(istep) + Raman_m * Lorentz_func
        ENDDO
        !
      ENDDO ! mu
      CLOSE(10)
      !
      ! output Raman data
      OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
      WRITE(20, '(a, 10x, a)') '#Raman_shift/cm-1', 'I_tot'
      DO istep = 1, nstep
        WRITE(20, '(f8.1, 5x, es22.12e3)') Rs_min + (istep-1) * Rs_inc, Raman_tot(istep)
      ENDDO
      CLOSE(20)
      !
      DEALLOCATE(Raman_tot)
      !
      CALL CPU_TIME(time)
      WRITE(*, '("total time:", f18.4, " secs")') time
      WRITE(*, *)
      !
    ENDIF ! my_id
    !
  ! consider hole-phonon scattering
  ! *********************************************************************
  ELSEIF (lhoph) THEN
  ! *********************************************************************
    IF (my_id == 0) THEN
      !
      CALL CPU_TIME(time)
      WRITE(*, '("begin time:",f18.4," secs")') time
      WRITE(*, *)
      !
      ALLOCATE(Raman_tot(nstep))
      ALLOCATE(Raman_e(nstep))
      ALLOCATE(Raman_h(nstep))
      ! Raman data initialization
      Raman_tot(:) = 0.d0
      Raman_e(:) = 0.d0
      Raman_h(:) = 0.d0
      !
      OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
      READ(10, *)
      READ(10, *)
      DO mu = 1, nmodes
        READ(10, '(7x, f10.4, 3(1x, 2(2x, es22.12e3)))') frq, Raman_p, Raman_m, &
                                    Raman_e_p, Raman_e_m, Raman_h_p, Raman_h_m
        ! calculate Raman data
        DO istep = 1, nstep
          ! tot
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc-frq)**2)
          Raman_tot(istep) = Raman_tot(istep) + Raman_p * Lorentz_func
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc+frq)**2)
          Raman_tot(istep) = Raman_tot(istep) + Raman_m * Lorentz_func
          ! e
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc-frq)**2)
          Raman_e(istep) = Raman_e(istep) + Raman_e_p * Lorentz_func
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc+frq)**2)
          Raman_e(istep) = Raman_e(istep) + Raman_e_m * Lorentz_func
          ! h
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc-frq)**2)
          Raman_h(istep) = Raman_h(istep) + Raman_h_p * Lorentz_func
          Lorentz_func = Lgamma / twopi / ((5.d-1*Lgamma)**2 + (Rs_min+(istep-1)*Rs_inc+frq)**2)
          Raman_h(istep) = Raman_h(istep) + Raman_h_m * Lorentz_func
        ENDDO
        !
      ENDDO ! mu
      CLOSE(10)
      !
      ! output Raman data
      OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
      WRITE(20, '(a, 10x, a, 17x, a, 19x, a)') '#Raman_shift/cm-1', 'I_tot', 'I_e', 'I_h'
      DO istep = 1, nstep
        WRITE(20, '(f8.1, 5x, 3es22.12e3)') Rs_min + (istep-1) * Rs_inc, Raman_tot(istep), &
                                                           Raman_e(istep), Raman_h(istep)
      ENDDO
      CLOSE(20)
      !
      DEALLOCATE(Raman_tot)
      DEALLOCATE(Raman_e)
      DEALLOCATE(Raman_h)
      !
      CALL CPU_TIME(time)
      WRITE(*, '("total time:", f18.4, " secs")') time
      WRITE(*, *)
      !
    ENDIF ! my_id
    !
  ENDIF ! lhoph
  !
  !
  !################################################################
  CASE ('double')   ! double resonance Raman 
  !################################################################
  !
  residual_value = MOD(nq, num_procs - 1)
  q_per_proc = (nq - residual_value) / (num_procs - 1)
  !
  ! not consider hole-phonon scattering
  ! *********************************************************************
  IF (.NOT. lhoph) THEN
  ! *********************************************************************
    IF (.NOT. lRaman_modes) THEN
      IF (my_id == 0) THEN
        !
        CALL CPU_TIME(time)
        WRITE(*, '("begin time:",f18.4," secs")') time
        WRITE(*, *)
        !
        ALLOCATE(Raman_tot(nstep))
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot(:) = 0.d0
        !
        ! output Raman data
        OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
        WRITE(20, '(a, 10x, a)') '#Raman_shift/cm-1', 'I_tot'
        DO proc_index = 1, num_procs-1
          CALL MPI_RECV(Raman_tot0, nstep, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          Raman_tot(:) = Raman_tot(:) + Raman_tot0(:)
        ENDDO
        Raman_tot(:) = Raman_tot(:) / nq
        !
        DO istep = 1, nstep
          WRITE(20, '(f8.1, 5x, es22.12e3)') Rs_min + (istep - 1) * Rs_inc, Raman_tot(istep)
        ENDDO
        CLOSE(20)
        !
        DEALLOCATE(Raman_tot)
        DEALLOCATE(Raman_tot0)
        !
        CALL CPU_TIME(time)
        WRITE(*, '("total time:", f18.4, " secs")') time
        WRITE(*, *)
        !
      ELSEIF (my_id <= residual_value) THEN
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
          CONTINUE
        ELSE
          DO iq = 1, (my_id - 1) * (q_per_proc + 1)
            READ(10, *)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          DO mu = 1, nmodes
            DO nu = 1, nmodes
              READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
              ! how to treat the strong intensity around 0 cm-1
              IF (.NOT. lRay_sca) THEN
                IF (ABS(frq - frq2) < Ray_thr) THEN
                  Raman_pm = 0.d0
                  Raman_mp = 0.d0
                ENDIF
              ENDIF
              ! calculate Raman data
              DO istep = 1, nstep
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mp * Lorentz_func
              ENDDO
              !
            ENDDO ! nu
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        !
      ELSE 
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
          CONTINUE
        ELSE
          DO iq = 1, residual_value + (my_id - 1) * q_per_proc
            READ(10, *)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          DO mu = 1, nmodes
            DO nu = 1, nmodes
              READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
              ! how to treat the strong intensity around 0 cm-1
              IF (.NOT. lRay_sca) THEN
                IF (ABS(frq - frq2) < Ray_thr) THEN
                  Raman_pm = 0.d0
                  Raman_mp = 0.d0
                ENDIF
              ENDIF
              ! calculate Raman data
              DO istep = 1, nstep
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mp * Lorentz_func
              ENDDO
              !
            ENDDO ! nu
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        !
      ENDIF ! my_id
      !
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ELSEIF (lRaman_modes) THEN
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO iRaman_modes = 1, nRaman_modes
        IF (my_id == 0) THEN
          !
          IF (iRaman_modes == 1) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("begin time:",f18.4," secs")') time
            WRITE(*, *)
          ENDIF
          !
          ALLOCATE(modes_intensity(nmodes, nmodes, 4))
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          ALLOCATE(location(3,nphonon_modes))
          ALLOCATE(mask(nmodes, nmodes, 4))
          ALLOCATE(frq0(nmodes))
          modes_intensity(:, :, :) = 0.d0
          mask(:, :, :) = .true.
          !
          ! read the specific total Raman intensiy 
          OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'old')
          READ(20, *)
          DO istep = 1, NINT((Raman_modes(iRaman_modes) - Rs_min) / Rs_inc)
            READ(20, *)
          ENDDO
          READ(20, '(13x, es22.12e3)') Raman_intensity
          CLOSE(20)
          !
          DO proc_index = 1, num_procs-1
            CALL MPI_RECV(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            modes_intensity(:, :, :) = modes_intensity(:, :, :) + modes_intensity0(:, :, :)
          ENDDO
          modes_intensity(:, :, :) = modes_intensity(:, :, :) / nq
          !
          ! extract the first nphonon_modes major contribution
          DO ipair = 1, nphonon_modes
            location(:, ipair) = MAXLOC(modes_intensity(:, :, :), mask(:, :, :))
            mask(location(1, ipair), location(2, ipair), location(3, ipair)) = .false.
          ENDDO
          !
          ! record the first nphonon_modes major contribution in the standard output file
          WRITE(*, '(a, i2.2, a, f8.1, a)') 'Raman mode ', iRaman_modes, ' : ', Raman_modes(iRaman_modes), ' cm-1'
          WRITE(*, '(a, i2, a)') 'the top ', nphonon_modes, ' significant branch-pairs:'
          DO ipair = 1, nphonon_modes
            IF (location(3, ipair) == 1) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' + ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 1), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 1) / Raman_intensity, '%'
            IF (location(3, ipair) == 2) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' - ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 2), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 2) / Raman_intensity, '%'
            IF (location(3, ipair) == 3) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' - ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 3), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 3) / Raman_intensity, '%'
            IF (location(3, ipair) == 4) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' + ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 4), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 4) / Raman_intensity, '%'
          ENDDO
          WRITE(*, *)
          WRITE(*, *)
          !
          ! write the q-resolved Raman intensity of the nphonon_modes branch-pairs
          WRITE(fmodes, '(a, i2.2)') 'Raman_mode', iRaman_modes
          OPEN(UNIT = 30, FILE = fmodes, STATUS = 'replace')
          WRITE(30, '(a, 3(6x, a), 10x, a, i2, a)') '#  iq', 'qx', 'qy', 'qz', &
                                                    'Raman intensity of top ', nphonon_modes, ' branch-pairs'
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          DO iq = 1, nq
            READ(10, '(9x, 3f10.6)') q(:)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq0(mu), frq0(nu), modes_intensity(mu, nu, :)
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq0(mu) - frq0(nu)) < Ray_thr) THEN
                    modes_intensity(mu, nu, 3) = 0.d0
                    modes_intensity(mu, nu, 4) = 0.d0
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            !
            WRITE(30, '(i5.5, 3f10.6)', ADVANCE='no') iq, q(:)
            DO ipair = 1, nphonon_modes
              IF (location(3, ipair) == 1) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) - frq0(location(1, ipair)) - frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 1) * Lorentz_func
              ELSEIF (location(3, ipair) == 2) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) + frq0(location(1, ipair)) + frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 2) * Lorentz_func
              ELSEIF (location(3, ipair) == 3) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) - frq0(location(1, ipair)) + frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 3) * Lorentz_func
              ELSEIF (location(3, ipair) == 4) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) + frq0(location(1, ipair)) - frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 4) * Lorentz_func
              ENDIF
            ENDDO ! ipair
            !
            WRITE(30, *)
            !
          ENDDO ! iq
          !
          CLOSE(10)
          CLOSE(30)
          !
          DEALLOCATE(modes_intensity)
          DEALLOCATE(modes_intensity0)
          DEALLOCATE(location)
          DEALLOCATE(mask)
          DEALLOCATE(frq0)
          !
          IF (iRaman_modes == nRaman_modes) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("total time:", f18.4, " secs")') time
            WRITE(*, *)
          ENDIF
          !
        ELSEIF (my_id <= residual_value) THEN
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
            CONTINUE
          ELSE
            DO iq = 1, (my_id - 1) * (q_per_proc + 1)
              READ(10, *)
              DO mu = 1, nmodes
                DO nu = 1, nmodes
                  READ(10, *)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq - frq2) < Ray_thr) THEN
                    Raman_pm = 0.d0
                    Raman_mp = 0.d0
                  ENDIF
                ENDIF
                ! contribution from respective branch combination
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq - frq2)**2)
                modes_intensity0(mu, nu, 1) = modes_intensity0(mu, nu, 1) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq + frq2)**2)
                modes_intensity0(mu, nu, 2) = modes_intensity0(mu, nu, 2) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq + frq2)**2)
                modes_intensity0(mu, nu, 3) = modes_intensity0(mu, nu, 3) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq - frq2)**2)
                modes_intensity0(mu, nu, 4) = modes_intensity0(mu, nu, 4) + Raman_mp * Lorentz_func
                !
              ENDDO ! nu
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ELSE
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
            CONTINUE
          ELSE
            DO iq = 1, residual_value + (my_id - 1) * q_per_proc
              READ(10, *)
              DO mu = 1, nmodes
                DO nu = 1, nmodes
                  READ(10, *)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq - frq2) < Ray_thr) THEN
                    Raman_pm = 0.d0
                    Raman_mp = 0.d0
                  ENDIF
                ENDIF
                ! contribution from respective branch combination
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq - frq2)**2)
                modes_intensity0(mu, nu, 1) = modes_intensity0(mu, nu, 1) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq + frq2)**2)
                modes_intensity0(mu, nu, 2) = modes_intensity0(mu, nu, 2) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq + frq2)**2)
                modes_intensity0(mu, nu, 3) = modes_intensity0(mu, nu, 3) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq - frq2)**2)
                modes_intensity0(mu, nu, 4) = modes_intensity0(mu, nu, 4) + Raman_mp * Lorentz_func
                !
              ENDDO ! nu
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ENDIF ! my_id
        !
      ENDDO ! iRaman_modes
      !
    ENDIF ! lRaman_modes
    !
  ! consider hole-phonon scattering
  ! *********************************************************************
  ELSEIF (lhoph) THEN
  ! *********************************************************************
    IF (.NOT. lRaman_modes) THEN
      IF (my_id == 0) THEN
        !
        CALL CPU_TIME(time)
        WRITE(*, '("begin time:",f18.4," secs")') time
        WRITE(*, *)
        !
        ALLOCATE(Raman_tot(nstep))
        ALLOCATE(Raman_ee(nstep))
        ALLOCATE(Raman_hh(nstep))
        ALLOCATE(Raman_eh(nstep))
        ALLOCATE(Raman_he(nstep))
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot(:) = 0.d0
        Raman_ee(:) = 0.d0
        Raman_hh(:) = 0.d0
        Raman_eh(:) = 0.d0
        Raman_he(:) = 0.d0
        !
        ! output Raman data
        OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
        WRITE(20, '(a, 10x, a, 17x, a, 18x, a, 18x, a, 18x, a)') '#Raman_shift/cm-1', 'I_tot', 'I_ee', 'I_hh', 'I_eh', 'I_he'
        DO proc_index = 1, num_procs-1
          CALL MPI_RECV(Raman_tot0, nstep, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_ee0, nstep, MPI_REAL8, proc_index, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_hh0, nstep, MPI_REAL8, proc_index, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_eh0, nstep, MPI_REAL8, proc_index, 400, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_he0, nstep, MPI_REAL8, proc_index, 500, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          Raman_tot(:) = Raman_tot(:) + Raman_tot0(:)
          Raman_ee(:) = Raman_ee(:) + Raman_ee0(:)
          Raman_hh(:) = Raman_hh(:) + Raman_hh0(:)
          Raman_eh(:) = Raman_eh(:) + Raman_eh0(:)
          Raman_he(:) = Raman_he(:) + Raman_he0(:)
        ENDDO
        Raman_tot(:) = Raman_tot(:) / nq
        Raman_ee(:) = Raman_ee(:) / nq
        Raman_hh(:) = Raman_hh(:) / nq
        Raman_eh(:) = Raman_eh(:) / nq
        Raman_he(:) = Raman_he(:) / nq
        !
        DO istep = 1, nstep
          WRITE(20, '(f8.1, 5x, 5es22.12e3)') Rs_min + (istep - 1) * Rs_inc, Raman_tot(istep), Raman_ee(istep), &
                                                             Raman_hh(istep), Raman_eh(istep), Raman_he(istep)
        ENDDO
        CLOSE(20)
        !
        DEALLOCATE(Raman_tot)
        DEALLOCATE(Raman_ee)
        DEALLOCATE(Raman_hh)
        DEALLOCATE(Raman_eh)
        DEALLOCATE(Raman_he)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
        CALL CPU_TIME(time)
        WRITE(*, '("total time:", f18.4, " secs")') time
        WRITE(*, *)
        !
      ELSEIF (my_id <= residual_value) THEN
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        Raman_ee0(:) = 0.d0
        Raman_hh0(:) = 0.d0
        Raman_eh0(:) = 0.d0
        Raman_he0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
          CONTINUE
        ELSE
          DO iq = 1, (my_id - 1) * (q_per_proc + 1)
            READ(10, *)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          DO mu = 1, nmodes
            DO nu = 1, nmodes
              READ(10, '(10x, 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') frq, frq2, &
                             Raman_pp, Raman_mm, Raman_pm, Raman_mp, &
                             Raman_ee_pp, Raman_ee_mm, Raman_ee_pm, Raman_ee_mp, &
                             Raman_hh_pp, Raman_hh_mm, Raman_hh_pm, Raman_hh_mp, &
                             Raman_eh_pp, Raman_eh_mm, Raman_eh_pm, Raman_eh_mp, &
                             Raman_he_pp, Raman_he_mm, Raman_he_pm, Raman_he_mp
              ! how to treat the strong intensity around 0 cm-1
              IF (.NOT. lRay_sca) THEN
                IF (ABS(frq - frq2) < Ray_thr) THEN
                  Raman_pm = 0.d0
                  Raman_mp = 0.d0
                  Raman_ee_pm = 0.d0
                  Raman_ee_mp = 0.d0
                  Raman_hh_pm = 0.d0
                  Raman_hh_mp = 0.d0
                  Raman_eh_pm = 0.d0
                  Raman_eh_mp = 0.d0
                  Raman_he_pm = 0.d0
                  Raman_he_mp = 0.d0
                ENDIF
              ENDIF
              ! calculate Raman data
              DO istep = 1, nstep
                ! tot
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mp * Lorentz_func
                ! ee
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_mp * Lorentz_func
                ! hh
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_mp * Lorentz_func
                ! eh
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_mp * Lorentz_func
                ! he
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_mp * Lorentz_func
              ENDDO
              !
            ENDDO ! nu
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_ee0, nstep, MPI_REAL8, 0, 200, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_hh0, nstep, MPI_REAL8, 0, 300, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_eh0, nstep, MPI_REAL8, 0, 400, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_he0, nstep, MPI_REAL8, 0, 500, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
      ELSE
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        Raman_ee0(:) = 0.d0
        Raman_hh0(:) = 0.d0
        Raman_eh0(:) = 0.d0
        Raman_he0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
          CONTINUE
        ELSE
          DO iq = 1, residual_value + (my_id - 1) * q_per_proc
            READ(10, *)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          DO mu = 1, nmodes
            DO nu = 1, nmodes
              READ(10, '(10x, 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') frq, frq2, &
                             Raman_pp, Raman_mm, Raman_pm, Raman_mp, &
                             Raman_ee_pp, Raman_ee_mm, Raman_ee_pm, Raman_ee_mp, &
                             Raman_hh_pp, Raman_hh_mm, Raman_hh_pm, Raman_hh_mp, &
                             Raman_eh_pp, Raman_eh_mm, Raman_eh_pm, Raman_eh_mp, &
                             Raman_he_pp, Raman_he_mm, Raman_he_pm, Raman_he_mp
              ! how to treat the strong intensity around 0 cm-1
              IF (.NOT. lRay_sca) THEN
                IF (ABS(frq - frq2) < Ray_thr) THEN
                  Raman_pm = 0.d0
                  Raman_mp = 0.d0
                  Raman_ee_pm = 0.d0
                  Raman_ee_mp = 0.d0
                  Raman_hh_pm = 0.d0
                  Raman_hh_mp = 0.d0
                  Raman_eh_pm = 0.d0
                  Raman_eh_mp = 0.d0
                  Raman_he_pm = 0.d0
                  Raman_he_mp = 0.d0
                ENDIF
              ENDIF
              ! calculate Raman data
              DO istep = 1, nstep
                ! tot
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_tot0(istep) = Raman_tot0(istep) + Raman_mp * Lorentz_func
                ! ee
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_mp * Lorentz_func
                ! hh
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_mp * Lorentz_func
                ! eh
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_mp * Lorentz_func
                ! he
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq - frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq + frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq + frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq - frq2)**2)
                Raman_he0(istep) = Raman_he0(istep) + Raman_he_mp * Lorentz_func
              ENDDO
              !
            ENDDO ! nu
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_ee0, nstep, MPI_REAL8, 0, 200, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_hh0, nstep, MPI_REAL8, 0, 300, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_eh0, nstep, MPI_REAL8, 0, 400, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_he0, nstep, MPI_REAL8, 0, 500, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
      ENDIF ! my_id
      !
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ELSEIF (lRaman_modes) THEN
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO iRaman_modes = 1, nRaman_modes
        IF (my_id == 0) THEN
          !
          IF (iRaman_modes == 1) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("begin time:",f18.4," secs")') time
            WRITE(*, *)
          ENDIF
          !
          ALLOCATE(modes_intensity(nmodes, nmodes, 4))
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          ALLOCATE(location(3,nphonon_modes))
          ALLOCATE(mask(nmodes, nmodes, 4))
          ALLOCATE(frq0(nmodes))
          modes_intensity(:, :, :) = 0.d0
          mask(:, :, :) = .true.
          !
          ! read the specific total Raman intensiy 
          OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'old')
          READ(20, *)
          DO istep = 1, NINT((Raman_modes(iRaman_modes) - Rs_min) / Rs_inc)
            READ(20, *)
          ENDDO
          READ(20, '(13x, es22.12e3)') Raman_intensity
          CLOSE(20)
          !
          DO proc_index = 1, num_procs-1
            CALL MPI_RECV(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            modes_intensity(:, :, :) = modes_intensity(:, :, :) + modes_intensity0(:, :, :)
          ENDDO
          modes_intensity(:, :, :) = modes_intensity(:, :, :) / nq
          !
          ! extract the first nphonon_modes major contribution
          DO ipair = 1, nphonon_modes
            location(:, ipair) = MAXLOC(modes_intensity(:, :, :), mask(:, :, :))
            mask(location(1, ipair), location(2, ipair), location(3, ipair)) = .false.
          ENDDO
          !
          ! record the first nphonon_modes major contribution in the standard output file
          WRITE(*, '(a, i2.2, a, f8.1, a)') 'Raman mode ', iRaman_modes, ' : ', Raman_modes(iRaman_modes), ' cm-1'
          WRITE(*, '(a, i2, a)') 'the top ', nphonon_modes, ' significant branch-pairs:'
          DO ipair = 1, nphonon_modes
            IF (location(3, ipair) == 1) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' + ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 1), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 1) / Raman_intensity, '%'
            IF (location(3, ipair) == 2) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' - ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 2), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 2) / Raman_intensity, '%'
            IF (location(3, ipair) == 3) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' - ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 3), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 3) / Raman_intensity, '%'
            IF (location(3, ipair) == 4) WRITE(*, '(2(a, i3.3), a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' + ', location(2, ipair), ' : ', &
                                 modes_intensity(location(1, ipair), location(2, ipair), 4), &
            1.d2*modes_intensity(location(1, ipair), location(2, ipair), 4) / Raman_intensity, '%'
          ENDDO
          WRITE(*, *)
          WRITE(*, *)
          !
          ! write the q-resolved Raman intensity of the nphonon_modes branch-pairs
          WRITE(fmodes, '(a, i2.2)') 'Raman_mode', iRaman_modes
          OPEN(UNIT = 30, FILE = fmodes, STATUS = 'replace')
          WRITE(30, '(a, 3(6x, a), 10x, a, i2, a)') '#  iq', 'qx', 'qy', 'qz', &
                                                    'Raman intensity of top ', nphonon_modes, ' branch-pairs'
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          DO iq = 1, nq
            READ(10, '(9x, 3f10.6)') q(:)
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq0(mu), frq0(nu), modes_intensity(mu, nu, :)
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq0(mu) - frq0(nu)) < Ray_thr) THEN
                    modes_intensity(mu, nu, 3) = 0.d0
                    modes_intensity(mu, nu, 4) = 0.d0
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            !
            WRITE(30, '(i5.5, 3f10.6)', ADVANCE='no') iq, q(:)
            DO ipair = 1, nphonon_modes
              IF (location(3, ipair) == 1) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) - frq0(location(1, ipair)) - frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 1) * Lorentz_func
              ELSEIF (location(3, ipair) == 2) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) + frq0(location(1, ipair)) + frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 2) * Lorentz_func
              ELSEIF (location(3, ipair) == 3) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) - frq0(location(1, ipair)) + frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 3) * Lorentz_func
              ELSEIF (location(3, ipair) == 4) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + &
                               (Raman_modes(iRaman_modes) + frq0(location(1, ipair)) - frq0(location(2, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 4) * Lorentz_func
              ENDIF
            ENDDO ! ipair
            !
            WRITE(30, *)
            !
          ENDDO ! iq
          !
          CLOSE(10)
          CLOSE(30)
          !
          DEALLOCATE(modes_intensity)
          DEALLOCATE(modes_intensity0)
          DEALLOCATE(location)
          DEALLOCATE(mask)
          DEALLOCATE(frq0)
          !
          IF (iRaman_modes == nRaman_modes) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("total time:", f18.4, " secs")') time
            WRITE(*, *)
          ENDIF
          !
        ELSEIF (my_id <= residual_value) THEN
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
            CONTINUE
          ELSE
            DO iq = 1, (my_id - 1) * (q_per_proc + 1)
              READ(10, *)
              DO mu = 1, nmodes
                DO nu = 1, nmodes
                  READ(10, *)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq - frq2) < Ray_thr) THEN
                    Raman_pm = 0.d0
                    Raman_mp = 0.d0
                  ENDIF
                ENDIF
                ! contribution from respective branch combination
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq - frq2)**2)
                modes_intensity0(mu, nu, 1) = modes_intensity0(mu, nu, 1) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq + frq2)**2)
                modes_intensity0(mu, nu, 2) = modes_intensity0(mu, nu, 2) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq + frq2)**2)
                modes_intensity0(mu, nu, 3) = modes_intensity0(mu, nu, 3) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq - frq2)**2)
                modes_intensity0(mu, nu, 4) = modes_intensity0(mu, nu, 4) + Raman_mp * Lorentz_func
                !
              ENDDO ! nu
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ELSE
          ALLOCATE(modes_intensity0(nmodes, nmodes, 4))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
            CONTINUE
          ELSE
            DO iq = 1, residual_value + (my_id - 1) * q_per_proc
              READ(10, *)
              DO mu = 1, nmodes
                DO nu = 1, nmodes
                  READ(10, *)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            DO mu = 1, nmodes
              DO nu = 1, nmodes
                READ(10, '(10x, 2(2x, f10.4), 1x, 4(2x, es22.12e3))') frq, frq2, Raman_pp, Raman_mm, Raman_pm, Raman_mp
                ! how to treat the strong intensity around 0 cm-1
                IF (.NOT. lRay_sca) THEN
                  IF (ABS(frq - frq2) < Ray_thr) THEN
                    Raman_pm = 0.d0
                    Raman_mp = 0.d0
                  ENDIF
                ENDIF
                ! contribution from respective branch combination
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq - frq2)**2)
                modes_intensity0(mu, nu, 1) = modes_intensity0(mu, nu, 1) + Raman_pp * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq + frq2)**2)
                modes_intensity0(mu, nu, 2) = modes_intensity0(mu, nu, 2) + Raman_mm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq + frq2)**2)
                modes_intensity0(mu, nu, 3) = modes_intensity0(mu, nu, 3) + Raman_pm * Lorentz_func
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq - frq2)**2)
                modes_intensity0(mu, nu, 4) = modes_intensity0(mu, nu, 4) + Raman_mp * Lorentz_func
                !
              ENDDO ! nu
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ENDIF ! my_id
        !
      ENDDO ! iRaman_modes
      !
    ENDIF ! lRaman_modes
    !
  ENDIF ! lhoph
  !
  !
  !################################################################
  CASE ('defect')   ! defect-induced double resonance Raman
  !################################################################
  !
  residual_value = MOD(nq, num_procs - 1)
  q_per_proc = (nq - residual_value) / (num_procs - 1)
  !
  ! not consider hole-phonon scattering
  ! *********************************************************************
  IF (.NOT. lhoph) THEN
  ! *********************************************************************
    IF (.NOT. lRaman_modes) THEN
      IF (my_id == 0) THEN
        !
        CALL CPU_TIME(time)
        WRITE(*, '("begin time:",f18.4," secs")') time
        WRITE(*, *)
        !
        ALLOCATE(Raman_tot(nstep))
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot(:) = 0.d0
        !
        ! output Raman data
        OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
        WRITE(20, '(a, 10x, a)') '#Raman_shift/cm-1', 'I_tot'
        DO proc_index = 1, num_procs-1
          CALL MPI_RECV(Raman_tot0, nstep, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          Raman_tot(:) = Raman_tot(:) + Raman_tot0(:)
        ENDDO
        Raman_tot(:) = Raman_tot(:) / nq
        !
        DO istep = 1, nstep
          WRITE(20, '(f8.1, 5x, es22.12e3)') Rs_min + (istep - 1) * Rs_inc, Raman_tot(istep)
        ENDDO
        CLOSE(20)
        !
        DEALLOCATE(Raman_tot)
        DEALLOCATE(Raman_tot0)
        !
        CALL CPU_TIME(time)
        WRITE(*, '("total time:", f18.4, " secs")') time
        WRITE(*, *)
        !
      ELSEIF (my_id <= residual_value) THEN
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
          CONTINUE
        ELSE
          DO iq = 1, (my_id - 1) * (q_per_proc + 1)
            READ(10, *)
            DO mu = 1, nmodes
              READ(10, *)
            ENDDO
            DO nu = 1, nmodes
              READ(10, *)
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
            ! calculate Raman data
            DO istep = 1, nstep
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
            ENDDO
            !
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
            ! calculate Raman data
            DO istep = 1, nstep
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
            ENDDO
            !
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        !
      ELSE 
        ALLOCATE(Raman_tot0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
          CONTINUE
        ELSE
          DO iq = 1, residual_value + (my_id - 1) * q_per_proc
            READ(10, *)
            DO mu = 1, nmodes
              READ(10, *)
            ENDDO
            DO nu = 1, nmodes
              READ(10, *)
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
            ! calculate Raman data
            DO istep = 1, nstep
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
            ENDDO
            !
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
            ! calculate Raman data
            DO istep = 1, nstep
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
            ENDDO
            !
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        !
      ENDIF ! my_id
      !
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ELSEIF (lRaman_modes) THEN
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO iRaman_modes = 1, nRaman_modes
        IF (my_id == 0) THEN
          !
          IF (iRaman_modes == 1) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("begin time:",f18.4," secs")') time
            WRITE(*, *)
          ENDIF
          !
          ALLOCATE(modes_intensity(nmodes, 2, 2))
          ! nmodes, phonon(1) or defect(2) scattering first, p(1) or m(2) 
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          ALLOCATE(location(3,nphonon_modes))
          ALLOCATE(mask(nmodes, 2, 2))
          ALLOCATE(frq0(nmodes))
          modes_intensity(:, :, :) = 0.d0
          mask(:, :, :) = .true.
          !
          ! read the specific total Raman intensity 
          OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'old')
          READ(20, *)
          DO istep = 1, NINT((Raman_modes(iRaman_modes) - Rs_min) / Rs_inc)
            READ(20, *)
          ENDDO
          READ(20, '(13x, es22.12e3)') Raman_intensity
          CLOSE(20)
          !
          DO proc_index = 1, num_procs-1
            CALL MPI_RECV(modes_intensity0, 4*nmodes, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            modes_intensity(:, :, :) = modes_intensity(:, :, :) + modes_intensity0(:, :, :)
          ENDDO
          modes_intensity(:, :, :) = modes_intensity(:, :, :) / nq
          !
          ! extract the first nphonon_modes major contribution
          DO ipair = 1, nphonon_modes
            location(:, ipair) = MAXLOC(modes_intensity(:, :, :), mask(:, :, :))
            mask(location(1, ipair), location(2, ipair), location(3, ipair)) = .false.
          ENDDO
          !
          ! record the first nphonon_modes major contribution in the standard output file
          WRITE(*, '(a, i2.2, a, f8.1, a)') 'Raman mode ', iRaman_modes, ' : ', Raman_modes(iRaman_modes), ' cm-1'
          WRITE(*, '(a, i2, a)') 'the top ', nphonon_modes, ' significant branch'
          DO ipair = 1, nphonon_modes
            IF (location(2, ipair) == 1 .and. location(3, ipair) == 1) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' + 0 : ', modes_intensity(location(1, ipair), 1, 1), &
            1.d2*modes_intensity(location(1, ipair), 1, 1) / Raman_intensity, '%'
            IF (location(2, ipair) == 1 .and. location(3, ipair) == 2) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' + 0 : ', modes_intensity(location(1, ipair), 1, 2), &
            1.d2*modes_intensity(location(1, ipair), 1, 2) / Raman_intensity, '%'
            IF (location(2, ipair) == 2 .and. location(3, ipair) == 1) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' 0 + ', location(1, ipair), ' : ', modes_intensity(location(1, ipair), 2, 1), &
            1.d2*modes_intensity(location(1, ipair), 2, 1) / Raman_intensity, '%'
            IF (location(2, ipair) == 2 .and. location(3, ipair) == 2) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' 0 - ', location(1, ipair), ' : ', modes_intensity(location(1, ipair), 2, 2), &
            1.d2*modes_intensity(location(1, ipair), 2, 2) / Raman_intensity, '%'
          ENDDO
          WRITE(*, *)
          WRITE(*, *)
          !
          ! write the q-resolved Raman intensity of the nphonon_modes branch
          WRITE(fmodes, '(a, i2.2)') 'Raman_mode', iRaman_modes
          OPEN(UNIT = 30, FILE = fmodes, STATUS = 'replace')
          WRITE(30, '(a, 3(6x, a), 10x, a, i2, a)') '#  iq', 'qx', 'qy', 'qz', &
                                                    'Raman intensity of top ', nphonon_modes, ' branches'
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          DO iq = 1, nq
            READ(10, '(9x, 3f10.6)') q(:)
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq0(mu), modes_intensity(mu, 1, :)
            ENDDO
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq0(mu), modes_intensity(mu, 2, :)
            ENDDO
            !
            WRITE(30, '(i5.5, 3f10.6)', ADVANCE='no') iq, q(:)
            DO ipair = 1, nphonon_modes
              IF (location(3, ipair) == 1) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq0(location(1, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 1) * Lorentz_func
              ELSEIF (location(3, ipair) == 2) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq0(location(1, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 2) * Lorentz_func
              ENDIF
            ENDDO ! ipair
            !
            WRITE(30, *)
            !
          ENDDO ! iq
          !
          CLOSE(10)
          CLOSE(30)
          !
          DEALLOCATE(modes_intensity)
          DEALLOCATE(modes_intensity0)
          DEALLOCATE(location)
          DEALLOCATE(mask)
          DEALLOCATE(frq0)
          !
          IF (iRaman_modes == nRaman_modes) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("total time:", f18.4, " secs")') time
            WRITE(*, *)
          ENDIF
          !
        ELSEIF (my_id <= residual_value) THEN
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
            CONTINUE
          ELSE
            DO iq = 1, (my_id - 1) * (q_per_proc + 1)
              READ(10, *)
              DO mu = 1, nmodes
                READ(10, *)
              ENDDO
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 1, 1) = modes_intensity0(mu, 1, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 1, 2) = modes_intensity0(mu, 1, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 2, 1) = modes_intensity0(mu, 2, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 2, 2) = modes_intensity0(mu, 2, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ELSE
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
            CONTINUE
          ELSE
            DO iq = 1, residual_value + (my_id - 1) * q_per_proc
              READ(10, *)
              DO mu = 1, nmodes
                READ(10, *)
              ENDDO
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 1, 1) = modes_intensity0(mu, 1, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 1, 2) = modes_intensity0(mu, 1, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 2, 1) = modes_intensity0(mu, 2, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 2, 2) = modes_intensity0(mu, 2, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ENDIF ! my_id
        !
      ENDDO ! iRaman_modes
      !
    ENDIF ! lRaman_modes
    !
  ! consider hole-phonon scattering
  ! *********************************************************************
  ELSEIF (lhoph) THEN
  ! *********************************************************************
    IF (.NOT. lRaman_modes) THEN
      IF (my_id == 0) THEN
        !
        CALL CPU_TIME(time)
        WRITE(*, '("begin time:",f18.4," secs")') time
        WRITE(*, *)
        !
        ALLOCATE(Raman_tot(nstep))
        ALLOCATE(Raman_ee(nstep))
        ALLOCATE(Raman_hh(nstep))
        ALLOCATE(Raman_eh(nstep))
        ALLOCATE(Raman_he(nstep))
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot(:) = 0.d0
        Raman_ee(:) = 0.d0
        Raman_hh(:) = 0.d0
        Raman_eh(:) = 0.d0
        Raman_he(:) = 0.d0
        !
        ! output Raman data
        OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'replace')
        WRITE(20, '(a, 10x, a, 17x, a, 18x, a, 18x, a, 18x, a)') '#Raman_shift/cm-1', 'I_tot', 'I_ee', 'I_hh', 'I_eh', 'I_he'
        DO proc_index = 1, num_procs-1
          CALL MPI_RECV(Raman_tot0, nstep, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_ee0, nstep, MPI_REAL8, proc_index, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_hh0, nstep, MPI_REAL8, proc_index, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_eh0, nstep, MPI_REAL8, proc_index, 400, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          CALL MPI_RECV(Raman_he0, nstep, MPI_REAL8, proc_index, 500, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          Raman_tot(:) = Raman_tot(:) + Raman_tot0(:)
          Raman_ee(:) = Raman_ee(:) + Raman_ee0(:)
          Raman_hh(:) = Raman_hh(:) + Raman_hh0(:)
          Raman_eh(:) = Raman_eh(:) + Raman_eh0(:)
          Raman_he(:) = Raman_he(:) + Raman_he0(:)
        ENDDO
        Raman_tot(:) = Raman_tot(:) / nq
        Raman_ee(:) = Raman_ee(:) / nq
        Raman_hh(:) = Raman_hh(:) / nq
        Raman_eh(:) = Raman_eh(:) / nq
        Raman_he(:) = Raman_he(:) / nq
        !
        DO istep = 1, nstep
          WRITE(20, '(f8.1, 5x, 5es22.12e3)') Rs_min + (istep - 1) * Rs_inc, Raman_tot(istep), Raman_ee(istep), &
                                                             Raman_hh(istep), Raman_eh(istep), Raman_he(istep)
        ENDDO
        CLOSE(20)
        !
        DEALLOCATE(Raman_tot)
        DEALLOCATE(Raman_ee)
        DEALLOCATE(Raman_hh)
        DEALLOCATE(Raman_eh)
        DEALLOCATE(Raman_he)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
        CALL CPU_TIME(time)
        WRITE(*, '("total time:", f18.4, " secs")') time
        WRITE(*, *)
        !
      ELSEIF (my_id <= residual_value) THEN
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        Raman_ee0(:) = 0.d0
        Raman_hh0(:) = 0.d0
        Raman_eh0(:) = 0.d0
        Raman_he0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
          CONTINUE
        ELSE
          DO iq = 1, (my_id - 1) * (q_per_proc + 1)
            READ(10, *)
            DO mu = 1, nmodes
              READ(10, *)
            ENDDO
            DO nu = 1, nmodes
              READ(10, *)
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 5(1x,2(2x, es22.12e3)))') frq, Raman_p, Raman_m, &
                                        Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                                        Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m
            ! calculate Raman data
            DO istep = 1, nstep
              ! tot
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
              ! ee
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_m * Lorentz_func
              ! hh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_m * Lorentz_func
              ! eh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_m * Lorentz_func
              ! he
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_m * Lorentz_func
            ENDDO
          ENDDO ! mu
          !
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 5(1x,2(2x, es22.12e3)))') frq, Raman_p, Raman_m, &
                                        Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                                        Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m
            ! calculate Raman data
            DO istep = 1, nstep
              ! tot
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
              ! ee
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_m * Lorentz_func
              ! hh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_m * Lorentz_func
              ! eh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_m * Lorentz_func
              ! he
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_m * Lorentz_func
            ENDDO
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_ee0, nstep, MPI_REAL8, 0, 200, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_hh0, nstep, MPI_REAL8, 0, 300, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_eh0, nstep, MPI_REAL8, 0, 400, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_he0, nstep, MPI_REAL8, 0, 500, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
      ELSE
        ALLOCATE(Raman_tot0(nstep))
        ALLOCATE(Raman_ee0(nstep))
        ALLOCATE(Raman_hh0(nstep))
        ALLOCATE(Raman_eh0(nstep))
        ALLOCATE(Raman_he0(nstep))
        ! Raman data initialization
        Raman_tot0(:) = 0.d0
        Raman_ee0(:) = 0.d0
        Raman_hh0(:) = 0.d0
        Raman_eh0(:) = 0.d0
        Raman_he0(:) = 0.d0
        !
        OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
        READ(10, *)
        IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
          CONTINUE
        ELSE
          DO iq = 1, residual_value + (my_id - 1) * q_per_proc
            READ(10, *)
            DO mu = 1, nmodes
              READ(10, *)
            ENDDO
            DO nu = 1, nmodes
              READ(10, *)
            ENDDO
          ENDDO
        ENDIF
        !
        DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
          READ(10, '(4x, i5)') iiq
          IF (iiq /= iq) THEN
            WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
            STOP
          ENDIF
          !
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 5(1x,2(2x, es22.12e3)))') frq, Raman_p, Raman_m, &
                                        Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                                        Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m
            ! calculate Raman data
            DO istep = 1, nstep
              ! tot
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
              ! ee
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_m * Lorentz_func
              ! hh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_m * Lorentz_func
              ! eh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_m * Lorentz_func
              ! he
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_m * Lorentz_func
            ENDDO
          ENDDO ! mu
          !
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            READ(10, '(10x, 2x, f10.4, 5(1x,2(2x, es22.12e3)))') frq, Raman_p, Raman_m, &
                                        Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                                        Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m
            ! calculate Raman data
            DO istep = 1, nstep
              ! tot
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_tot0(istep) = Raman_tot0(istep) + Raman_m * Lorentz_func
              ! ee
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_ee0(istep) = Raman_ee0(istep) + Raman_ee_m * Lorentz_func
              ! hh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_hh0(istep) = Raman_hh0(istep) + Raman_hh_m * Lorentz_func
              ! eh
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_eh0(istep) = Raman_eh0(istep) + Raman_eh_m * Lorentz_func
              ! he
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc - frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Rs_min + (istep - 1) * Rs_inc + frq)**2)
              Raman_he0(istep) = Raman_he0(istep) + Raman_he_m * Lorentz_func
            ENDDO
          ENDDO ! mu
          !
        ENDDO ! iq
        CLOSE(10)
        !
        CALL MPI_SEND(Raman_tot0, nstep, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_ee0, nstep, MPI_REAL8, 0, 200, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_hh0, nstep, MPI_REAL8, 0, 300, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_eh0, nstep, MPI_REAL8, 0, 400, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(Raman_he0, nstep, MPI_REAL8, 0, 500, MPI_COMM_WORLD, ierr)
        DEALLOCATE(Raman_tot0)
        DEALLOCATE(Raman_ee0)
        DEALLOCATE(Raman_hh0)
        DEALLOCATE(Raman_eh0)
        DEALLOCATE(Raman_he0)
        !
      ENDIF ! my_id
      !
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ELSEIF (lRaman_modes) THEN
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO iRaman_modes = 1, nRaman_modes
        IF (my_id == 0) THEN
          !
          IF (iRaman_modes == 1) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("begin time:",f18.4," secs")') time
            WRITE(*, *)
          ENDIF
          !
          ALLOCATE(modes_intensity(nmodes, 2, 2))
          ! nmodes, phonon(1) or defect(2) scattering first, p(1) or m(2)
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          ALLOCATE(location(3,nphonon_modes))
          ALLOCATE(mask(nmodes, 2, 2))
          ALLOCATE(frq0(nmodes))
          modes_intensity(:, :, :) = 0.d0
          mask(:, :, :) = .true.
          !
          ! read the specific total Raman intensiy 
          OPEN(UNIT = 20, FILE = 'Raman.dat', STATUS = 'old')
          READ(20, *)
          DO istep = 1, NINT((Raman_modes(iRaman_modes) - Rs_min) / Rs_inc)
            READ(20, *)
          ENDDO
          ! only read the total Raman intensity
          READ(20, '(13x, es22.12e3)') Raman_intensity
          CLOSE(20)
          !
          DO proc_index = 1, num_procs-1
            CALL MPI_RECV(modes_intensity0, 4*nmodes, MPI_REAL8, proc_index, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            modes_intensity(:, :, :) = modes_intensity(:, :, :) + modes_intensity0(:, :, :)
          ENDDO
          modes_intensity(:, :, :) = modes_intensity(:, :, :) / nq
          !
          ! extract the first nphonon_modes major contribution
          DO ipair = 1, nphonon_modes
            location(:, ipair) = MAXLOC(modes_intensity(:, :, :), mask(:, :, :))
            mask(location(1, ipair), location(2, ipair), location(3, ipair)) = .false.
          ENDDO
          !
          ! record the first nphonon_modes major contribution in the standard output file
          WRITE(*, '(a, i2.2, a, f8.1, a)') 'Raman mode ', iRaman_modes, ' : ', Raman_modes(iRaman_modes), ' cm-1'
          WRITE(*, '(a, i2, a)') 'the top ', nphonon_modes, ' significant branch'
          DO ipair = 1, nphonon_modes
            IF (location(2, ipair) == 1 .and. location(3, ipair) == 1) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' + ', location(1, ipair), ' + 0 : ', modes_intensity(location(1, ipair), 1, 1), &
            1.d2*modes_intensity(location(1, ipair), 1, 1) / Raman_intensity, '%'
            IF (location(2, ipair) == 1 .and. location(3, ipair) == 2) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' - ', location(1, ipair), ' + 0 : ', modes_intensity(location(1, ipair), 1, 2), &
            1.d2*modes_intensity(location(1, ipair), 1, 2) / Raman_intensity, '%'
            IF (location(2, ipair) == 2 .and. location(3, ipair) == 1) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' 0 + ', location(1, ipair), ' : ', modes_intensity(location(1, ipair), 2, 1), &
            1.d2*modes_intensity(location(1, ipair), 2, 1) / Raman_intensity, '%'
            IF (location(2, ipair) == 2 .and. location(3, ipair) == 2) WRITE(*, '(a, i3.3, a, es22.12e3, 5x, f5.2, a)') &
                                ' 0 - ', location(1, ipair), ' : ', modes_intensity(location(1, ipair), 2, 2), &
            1.d2*modes_intensity(location(1, ipair), 2, 2) / Raman_intensity, '%'
          ENDDO
          WRITE(*, *)
          WRITE(*, *)
          !
          ! write the q-resolved Raman intensity of the nphonon_modes branch
          WRITE(fmodes, '(a, i2.2)') 'Raman_mode', iRaman_modes
          OPEN(UNIT = 30, FILE = fmodes, STATUS = 'replace')
          WRITE(30, '(a, 3(6x, a), 10x, a, i2, a)') '#  iq', 'qx', 'qy', 'qz', &
                                                    'Raman intensity of top ', nphonon_modes, ' branches'
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          DO iq = 1, nq
            READ(10, '(9x, 3f10.6)') q(:)
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq0(mu), modes_intensity(mu, 1, :)
            ENDDO
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq0(mu), modes_intensity(mu, 2, :)
            ENDDO
            !
            WRITE(30, '(i5.5, 3f10.6)', ADVANCE='no') iq, q(:)
            DO ipair = 1, nphonon_modes
              IF (location(3, ipair) == 1) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq0(location(1, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 1) * Lorentz_func
              ELSEIF (location(3, ipair) == 2) THEN
                Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq0(location(1, ipair)))**2)
                WRITE(30, '(2x, es22.12e3)', ADVANCE='no') modes_intensity(location(1,ipair), location(2,ipair), 2) * Lorentz_func
              ENDIF
            ENDDO ! ipair
            !
            WRITE(30, *)
            !
          ENDDO ! iq
          !
          CLOSE(10)
          CLOSE(30)
          !
          DEALLOCATE(modes_intensity)
          DEALLOCATE(modes_intensity0)
          DEALLOCATE(location)
          DEALLOCATE(mask)
          DEALLOCATE(frq0)
          !
          IF (iRaman_modes == nRaman_modes) THEN
            CALL CPU_TIME(time)
            WRITE(*, '("total time:", f18.4, " secs")') time
            WRITE(*, *)
          ENDIF
          !
        ELSEIF (my_id <= residual_value) THEN
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > (my_id - 1) * (q_per_proc + 1)) THEN
            CONTINUE
          ELSE
            DO iq = 1, (my_id - 1) * (q_per_proc + 1)
              READ(10, *)
              DO mu = 1, nmodes
                READ(10, *)
              ENDDO
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + (my_id - 1) * (q_per_proc + 1), 1 + my_id * (q_per_proc + 1) - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              ! only read the total Raman intensity
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 1, 1) = modes_intensity0(mu, 1, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 1, 2) = modes_intensity0(mu, 1, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              ! only read the total Raman intensity
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 2, 1) = modes_intensity0(mu, 2, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 2, 2) = modes_intensity0(mu, 2, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ELSE
          ALLOCATE(modes_intensity0(nmodes, 2, 2))
          modes_intensity0(:, :, :) = 0.d0
          !
          OPEN(UNIT = 10, FILE = TRIM(dir_raman), STATUS = 'old')
          READ(10, *)
          IF (1 > residual_value + (my_id - 1) * q_per_proc) THEN
            CONTINUE
          ELSE
            DO iq = 1, residual_value + (my_id - 1) * q_per_proc
              READ(10, *)
              DO mu = 1, nmodes
                READ(10, *)
              ENDDO
              DO nu = 1, nmodes
                READ(10, *)
              ENDDO
            ENDDO
          ENDIF
          !
          DO iq = 1 + residual_value + (my_id - 1) * q_per_proc, 1 + residual_value + my_id * q_per_proc - 1
            READ(10, '(4x, i5)') iiq
            IF (iiq /= iq) THEN
              WRITE(*, '("errors in reading the", i5, " q-point of dir_raman")') iq
              STOP
            ENDIF
            ! first phonon scattering then defect scattering
            DO mu = 1, nmodes
              ! only read the total Raman intensity
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 1, 1) = modes_intensity0(mu, 1, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 1, 2) = modes_intensity0(mu, 1, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            ! first defect scattering then phonon scattering
            DO mu = 1, nmodes
              ! only read the total Raman intensity
              READ(10, '(10x, 2x, f10.4, 1x, 2(2x, es22.12e3))') frq, Raman_p, Raman_m
              ! contribution from respective branch combination
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) - frq)**2)
              modes_intensity0(mu, 2, 1) = modes_intensity0(mu, 2, 1) + Raman_p * Lorentz_func
              Lorentz_func = Lgamma / twopi / ((5.d-1 * Lgamma)**2 + (Raman_modes(iRaman_modes) + frq)**2)
              modes_intensity0(mu, 2, 2) = modes_intensity0(mu, 2, 2) + Raman_m * Lorentz_func
              !
            ENDDO ! mu
            !
          ENDDO ! iq
          CLOSE(10)
          !
          CALL MPI_SEND(modes_intensity0, 4*nmodes, MPI_REAL8, 0, 100, MPI_COMM_WORLD, ierr)
          DEALLOCATE(modes_intensity0)
          !
        ENDIF ! my_id
        !
      ENDDO ! iRaman_modes
      !
    ENDIF ! lRaman_modes
    !
  ENDIF ! lhoph
  !
  !
  ENDSELECT ! select
  !
  CALL MPI_FINALIZE(ierr)
  !
  STOP
END PROGRAM