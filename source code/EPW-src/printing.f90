  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE printing
  !----------------------------------------------------------------------
  !!
  !! This module contains various printing routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_dipole()
    !-----------------------------------------------------------------------
    !!
    !! Added for Raman calculations, by Jianqi Huang
    !! Print the complex dipole matrix elements on the fine interpolated k mesh
    !! with the involved electronic bands located in the Wannier window (nbndsub).
    !!
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE pwcom,            ONLY : ef
    USE epwcom,           ONLY : nbndsub
    USE elph2,            ONLY : dmef, etf, nkqf, xqf, xkf, nkqtotf, nktotf
    USE constants_epw,    ONLY : ryd2ev, zero, czero
    USE mp_global,        ONLY : inter_pool_comm
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : ionode_id
    USE poolgathering,    ONLY : poolgather2, poolgatherc4
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! K-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: nbnd_occ
    !! the number of occupied bands
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: ek
    !! Eigenenergies at k
    REAL(KIND = DP) :: ekk
    !! Eigenenergies at k
    REAL(KIND = DP) :: xkf_all(3, nkqtotf)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP) :: etf_all(nbndsub, nkqtotf)
    !! Collect eigenenergies from all pools in parallel case
    COMPLEX(KIND = DP), ALLOCATABLE :: dmef_all(:, :, :, :)
    !! Collect dipole matrix from all pools in parallel case
    !
    ALLOCATE(dmef_all(3, nbndsub, nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('print_dipole', 'Error allocating dmef_all', 1)
    !
    !
    ! We need quantity from all the pools
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    dmef_all(:, :, :, :) = czero
    !
#if defined(__MPI)
    !
    ! Note that poolgather2 and poolgatherc4 works with the doubled grid (k and k+q)
    !
    CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL poolgatherc4(3, nbndsub, nbndsub, nkqtotf, nkqf, dmef, dmef_all)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    dmef_all = dmef
    !
#endif
    !
    ! Only master writes
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = 10, FILE = 'dipole.dat', STATUS = 'replace')
      WRITE (10, '(" &p_mat nbndsub=", i4, ", nktotf=", i5, " /")') nbndsub, nktotf
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        !
        ! Compute the number of occupied bands at this k point
        DO ibnd = 1, nbndsub
          IF (etf_all(ibnd, ikk) <= ef) nbnd_occ = ibnd
        ENDDO
        ! IF (nbnd_occ == nbndsub) CALL errore('print_dipole', 'No empty band at k point', 1)
        !
        WRITE (10, '(10x, 3f10.6, i7)') xkf_all(:, ikk), nbnd_occ
        !
        ! valence bands in fsthick window
        DO ibnd = 1, nbnd_occ
          ek = etf_all(ibnd, ikk)
          ! conduction bands in fsthick window
          DO jbnd = nbnd_occ + 1, nbndsub
            ekk = etf_all(jbnd, ikk)
            WRITE (10, '(2i6, 2f11.5, $)') ibnd, jbnd, ryd2ev * ek, ryd2ev * ekk
            WRITE (10, '(5x, "(", es20.12e3, ",", es20.12e3, ")", $)') &
            REAL(dmef_all(1, jbnd, ibnd, ikk)), AIMAG(dmef_all(1, jbnd, ibnd, ikk))
            WRITE (10, '(5x, "(", es20.12e3, ",", es20.12e3, ")", $)') &
            REAL(dmef_all(2, jbnd, ibnd, ikk)), AIMAG(dmef_all(2, jbnd, ibnd, ikk))
            WRITE (10, '(5x, "(", es20.12e3, ",", es20.12e3, ")")') &
            REAL(dmef_all(3, jbnd, ibnd, ikk)), AIMAG(dmef_all(3, jbnd, ibnd, ikk))
          ENDDO
        ENDDO
        !
      ENDDO ! ik
      !
      CLOSE(10)
      !
    ENDIF ! master node
    !
    DEALLOCATE(dmef_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('print_dipole', 'Error deallocating dmef_all', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_dipole
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_Raman(iq)
    !-----------------------------------------------------------------------
    !!
    !! Added for Raman calculations, by Jianqi Huang
    !! Print the Raman spectra involving one- and two-phonon cases calculated
    !! with the help of the dipole matrix (dmef) and the el-ph matrix (epf17),
    !! and output the mode-resolved Raman intensity and Raman tensor optionally.
    !!
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE pwcom,            ONLY : ef
    USE modes,            ONLY : nmodes
    USE epwcom,           ONLY : nbndsub, eps_acustic, Raman_type, lhoph, reson_lim, reson_thr, Elaser, Egamma, &
                                 ltensor, qtensor_start, qtensor_end, lBE, temphon, polar, filraman, ei, es, Cq
    USE elph2,            ONLY : etf, ibndmin, ibndmax, nkf, nktotf, dmef, epf17, wf
    USE constants_epw,    ONLY : rydcm1, half, one, cone, two, zero, czero, ci
    USE mp,               ONLY : mp_barrier, mp_sum
    USE mp_global,        ONLY : inter_pool_comm
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    !
    ! Local variables
    CHARACTER(LEN = 13) :: ftensor
    !! file name of q-resolved Raman tensor
    INTEGER :: mu, nu
    !! mode index
    INTEGER :: ik, ikk
    !! k-point index
    INTEGER :: ikq
    !! k+q-point index
    INTEGER :: iiq
    !! q-point index
    INTEGER :: ibnd, jbnd, kbnd, mbnd
    !! band index
    INTEGER :: irow
    !! index of row
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: branch(:)
    !! ordering phonon branch
    REAL(KIND = DP) :: ql
    !! length of q-wavevector
    REAL(KIND = DP) :: frq, frq2
    !! temporary phonon freq. 1 and freq. 2
    REAL(KIND = DP) :: inv_frq, inv_frq2
    !! inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND = DP) :: e_ini, e_exc, e_int, e_int2
    !! initial, excited, intermediate and second intermediate eigenenergy
    !
    ! temporary Raman intensity for Raman_type = 'single'
    REAL(KIND = DP) :: Raman_p, Raman_m, Raman_e_p, Raman_e_m, Raman_h_p, Raman_h_m
    ! temporary Raman intensity for Raman_type = 'double'
    REAL(KIND = DP) :: Raman_pp, Raman_mm, Raman_pm, Raman_mp, &
                       Raman_ee_pp, Raman_ee_mm, Raman_ee_pm, Raman_ee_mp, &
                       Raman_hh_pp, Raman_hh_mm, Raman_hh_pm, Raman_hh_mp, &
                       Raman_eh_pp, Raman_eh_mm, Raman_eh_pm, Raman_eh_mp, &
                       Raman_he_pp, Raman_he_mm, Raman_he_pm, Raman_he_mp
    ! temporary Raman intensity for Raman_type = 'defect'
    REAL(KIND = DP) :: Raman_ee_p, Raman_ee_m, Raman_hh_p, Raman_hh_m, &
                       Raman_eh_p, Raman_eh_m, Raman_he_p, Raman_he_m
    !
    REAL(KIND = DP) :: q(3)
    !! q-coordinates in Cartesian framework
    REAL(KIND = DP), ALLOCATABLE :: bewph_as(:), bewph_s(:)
    !! Bose-Einstein distribution at frequency wph of anti-stokes and stokes phonon
    COMPLEX(KIND = DP) :: dipole_inc(3), dipole_sca(3)
    !! dipole vector in incident and scattered light process
    COMPLEX(KIND = DP) :: elphme, elphme2, Ide
    !! the first, second el-ph matrix element, and the intensity of defect-electron scattering
    COMPLEX(KIND = DP) :: factor
    !! denominator factor for optical transition term
    COMPLEX(KIND = DP) :: factor_ep, factor_ep2
    !! denominator factor for first and second el-ph scattering term
    COMPLEX(KIND = DP) :: tensor_tem(3, 3)
    !! temporary Raman tensor
    !
    ! for Raman_type = 'single'
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_p(:, :, :), tensor_m(:, :, :)
    !! total Raman tensor for one-phonon emission and absorption
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_e_p(:, :, :), tensor_e_m(:, :, :)
    !! Raman tensor for one-phonon emission and absorption with el-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_h_p(:, :, :), tensor_h_m(:, :, :)
    !! Raman tensor for one-phonon emission and absorption with ho-ph
    !
    ! for Raman_type = 'double'
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_pp(:, :, :, :), tensor_mm(:, :, :, :)
    !! total Raman tensor for two-phonon emission and absorption
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_pm(:, :, :, :), tensor_mp(:, :, :, :)
    !! total Raman tensor for first phonon emission (absorption) and second absorption (emission)
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_ee_pp(:, :, :, :), tensor_ee_mm(:, :, :, :)
    !! Raman tensor for two-phonon emission (absorption) with el-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_ee_pm(:, :, :, :), tensor_ee_mp(:, :, :, :)
    !! Raman tensor for first phonon emission (absorption) and second absorption (emission) with el-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_hh_pp(:, :, :, :), tensor_hh_mm(:, :, :, :)
    !! Raman tensor for two-phonon emission (absorption) with ho-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_hh_pm(:, :, :, :), tensor_hh_mp(:, :, :, :)
    !! Raman tensor for first phonon emission (absorption) and second absorption (emission) with ho-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_eh_pp(:, :, :, :), tensor_eh_mm(:, :, :, :)
    !! Raman tensor for two-phonon emission (absorption) with first el-ph and second ho-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_eh_pm(:, :, :, :), tensor_eh_mp(:, :, :, :)
    !! Raman tensor for first phonon emission (absorption) and second absorption (emission)
    !! with first el-ph and second ho-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_he_pp(:, :, :, :), tensor_he_mm(:, :, :, :)
    !! Raman tensor for two-phonon emission (absorption) with first ho-ph and second el-ph
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_he_pm(:, :, :, :), tensor_he_mp(:, :, :, :)
    !! Raman tensor for first phonon emission (absorption) and second absorption (emission)
    !! with first ho-ph and second el-ph
    !
    ! for Raman_type = 'defect'
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_p1(:, :, :), tensor_p2(:, :, :)
    !! total Raman tensor for one-phonon emission. 
    !! 1: first phonon then defect scattering; 2: first defect then phonon scattering
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_m1(:, :, :), tensor_m2(:, :, :)
    !! total Raman tensor for one-phonon absorption
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_ee_p1(:, :, :), tensor_ee_p2(:, :, :)
    !! Raman tensor for one-phonon emission with el-ph and el-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_ee_m1(:, :, :), tensor_ee_m2(:, :, :)
    !! Raman tensor for one-phonon absorption with el-ph and el-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_hh_p1(:, :, :), tensor_hh_p2(:, :, :)
    !! Raman tensor for one-phonon emission with ho-ph and ho-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_hh_m1(:, :, :), tensor_hh_m2(:, :, :)
    !! Raman tensor for one-phonon absorption with ho-ph and ho-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_eh_p1(:, :, :), tensor_eh_p2(:, :, :)
    !! Raman tensor for one-phonon emission with el-ph and ho-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_eh_m1(:, :, :), tensor_eh_m2(:, :, :)
    !! Raman tensor for one-phonon absorption with el-ph and ho-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_he_p1(:, :, :), tensor_he_p2(:, :, :)
    !! Raman tensor for one-phonon emission with ho-ph and el-defect
    COMPLEX(KIND = DP), ALLOCATABLE :: tensor_he_m1(:, :, :), tensor_he_m2(:, :, :)
    !! Raman tensor for one-phonon absorption with ho-ph and el-defect
    LOGICAL :: alive
    !! if qraman file exist
    !
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    !
    !
    SELECT CASE (Raman_type)
    !
    !################################################################
    CASE ('single')   ! single-resonance Raman 
    !################################################################
    IF (iq == 1) THEN
      ! not consider hole-phonon scattering
      ! *********************************************************************
      IF (.NOT. lhoph) THEN
      ! *********************************************************************
        ! allocate Raman tensor
        ALLOCATE(tensor_p(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p', 1)
        ALLOCATE(tensor_m(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m', 1)
        ! allocate Bose-Einstein occupation
        IF (lBE) ALLOCATE(bewph_as(nmodes))
        IF (lBE) ALLOCATE(bewph_s(nmodes))
        !
        tensor_p(:, :, :) = czero
        tensor_m(:, :, :) = czero
        !
        DO mu = 1, nmodes
          ! phonon frequency [Ry]
          frq =  wf(mu, iq)
          ! Min phonon frequency for el-ph calculations
          IF (frq > eps_acustic) THEN
            ! Define the inverse so that we can efficiently multiply instead of dividing
            inv_frq = one / DSQRT(two * frq)
            ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
            IF (lBE) THEN
              bewph_as(mu) = wgauss(frq / temphon, -99)
              bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
              bewph_s(mu)  = bewph_as(mu) + one
            ENDIF
          ELSE
            inv_frq = zero
            IF (lBE) bewph_as(mu) = zero
            IF (lBE) bewph_s(mu) = zero
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            DO ibnd = ibndmin, ibndmax
              ! the initial state at k is sited in valence band
              IF (etf(ibnd, ikk) > ef) CYCLE
              ! Electron eigenenergy [Ry] at initial state at k
              ! note the bands window (nbndsub) defined in etf
              e_ini = etf(ibnd, ikk)
              DO jbnd = ibndmin, ibndmax
                ! the excited state at k is sited in conduction band
                IF (etf(jbnd, ikk) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at excited state at k
                e_exc = etf(jbnd, ikk)
                ! if .true., limiting the energy level seperation 
                ! for resonance consideration
                IF (reson_lim) THEN
                  IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                ! complex denominator term [Ry]
                factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
                ! dipole vector in incident light process
                ! dmef: interpolated dipole matrix elements in units of 1/Bohr
                ! note the args order and bands window (nbndsub) defined in dmef
                dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
                !
                !! Starting the el-ph Raman process
                !
                !! emit one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at the intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, kbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! total Raman tensor
                  tensor_p(:, :, mu) = tensor_p(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
                !! absorb one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at the intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, kbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! total Raman tensor
                  tensor_m(:, :, mu) = tensor_m(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
              ENDDO ! jbnd
            ENDDO ! ibnd
          ENDDO ! ik
        ENDDO ! mu
        !
        !
#if defined(__MPI)
        !
        ! We need quantity from all the pools
        !
        CALL mp_sum(tensor_p, inter_pool_comm )
        CALL mp_sum(tensor_m, inter_pool_comm )
        CALL mp_barrier(inter_pool_comm)
        !
#endif
        !
        ! Only master writes
        IF (mpime == ionode_id) THEN
          ! output Raman tensor
          IF (ltensor) THEN
            WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
            OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
            WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
            WRITE(10, *)
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              WRITE(10, '(a)') 'total Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' = ', frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' = ', -frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('%', 100)
              WRITE(10, *)
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ENDIF ! ltensor
          !
          ! output mode-resolved Raman intensity
          IF (TRIM(polar) == 'custom') THEN
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ELSEIF (TRIM(polar) == 'all') THEN
            ! XX
            ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! YY
            ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! XY
            ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! YX
            ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! LL
            ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! RR
            ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! LR
            ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! RL
            ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'replace')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, frq, Raman_p/nktotf, Raman_m/nktotf
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ENDIF ! polar
          !
        ENDIF ! master node
        !
        DEALLOCATE(tensor_p, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p', 1)
        DEALLOCATE(tensor_m, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m', 1)
        !
        IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
        IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
        !
        !
      ! consider hole-phonon scattering
      ! *********************************************************************
      ELSEIF (lhoph) THEN
      ! *********************************************************************
        ! allocate Raman tensor
        ALLOCATE(tensor_p(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p', 1)
        ALLOCATE(tensor_m(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m', 1)
        ALLOCATE(tensor_e_p(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_e_p', 1)
        ALLOCATE(tensor_e_m(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_e_m', 1)
        ALLOCATE(tensor_h_p(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_h_p', 1)
        ALLOCATE(tensor_h_m(3, 3, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_h_m', 1)
        ! allocate Bose-Einstein occupation
        IF (lBE) ALLOCATE(bewph_as(nmodes))
        IF (lBE) ALLOCATE(bewph_s(nmodes))
        !
        tensor_p(:, :, :) = czero
        tensor_m(:, :, :) = czero
        tensor_e_p(:, :, :) = czero
        tensor_e_m(:, :, :) = czero
        tensor_h_p(:, :, :) = czero
        tensor_h_m(:, :, :) = czero
        !
        DO mu = 1, nmodes
          ! phonon frequency [Ry]
          frq =  wf(mu, iq)
          ! Min phonon frequency for el-ph calculations
          IF (frq > eps_acustic) THEN
            ! Define the inverse so that we can efficiently multiply instead of dividing
            inv_frq = one / DSQRT(two * frq)
            ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
            IF (lBE) THEN
              bewph_as(mu) = wgauss(frq / temphon, -99)
              bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
              bewph_s(mu)  = bewph_as(mu) + one
            ENDIF
          ELSE
            inv_frq = zero
            IF (lBE) bewph_as(mu) = zero
            IF (lBE) bewph_s(mu) = zero
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            DO ibnd = ibndmin, ibndmax
              ! the initial state at k is sited in valence band
              IF (etf(ibnd, ikk) > ef) CYCLE
              ! Electron eigenenergy [Ry] at initial state at k
              ! note the bands window (nbndsub) defined in etf
              e_ini = etf(ibnd, ikk)
              DO jbnd = ibndmin, ibndmax
                ! the excited state at k is sited in conduction band
                IF (etf(jbnd, ikk) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at excited state at k
                e_exc = etf(jbnd, ikk)
                ! if .true., limiting the energy level seperation 
                ! for resonance consideration
                IF (reson_lim) THEN
                  IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                ! complex denominator term [Ry]
                factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
                ! dipole vector in incident light process
                ! dmef: interpolated dipole matrix elements in units of 1/Bohr
                ! note the args order and bands window (nbndsub) defined in dmef
                dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
                !
                !! Starting the el-ph Raman process
                !
                !! emit one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at the intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, kbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! el-ph Raman tensor
                  tensor_e_p(:, :, mu) = tensor_e_p(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
                !! absorb one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at the intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, kbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! el-ph Raman tensor
                  tensor_e_m(:, :, mu) = tensor_e_m(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
                !! Starting the ho-ph Raman process
                !
                !! emit one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at the intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! ho-ph Raman tensor
                  tensor_h_p(:, :, mu) = tensor_h_p(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
                !! absorb one phonon
                !
                DO kbnd = ibndmin, ibndmax
                  ! the intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at the intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-ph coupling [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * factor * factor_ep
                  ! ho-ph Raman tensor
                  tensor_h_m(:, :, mu) = tensor_h_m(:, :, mu) + tensor_tem(:, :)
                  !
                ENDDO ! kbnd
                !
              ENDDO ! jbnd
            ENDDO ! ibnd
          ENDDO ! ik
        ENDDO ! mu
        !
        ! total Raman tensor
        tensor_p(:, :, :) = tensor_e_p(:, :, :) + tensor_h_p(:, :, :)
        tensor_m(:, :, :) = tensor_e_m(:, :, :) + tensor_h_m(:, :, :)
        !
        !
#if defined(__MPI)
        !
        ! We need quantity from all the pools
        !
        CALL mp_sum(tensor_p, inter_pool_comm )
        CALL mp_sum(tensor_m, inter_pool_comm )
        CALL mp_sum(tensor_e_p, inter_pool_comm )
        CALL mp_sum(tensor_e_m, inter_pool_comm )
        CALL mp_sum(tensor_h_p, inter_pool_comm )
        CALL mp_sum(tensor_h_m, inter_pool_comm )
        CALL mp_barrier(inter_pool_comm)
        !
#endif
        !
        ! Only master writes
        IF (mpime == ionode_id) THEN
          ! output Raman tensor
          IF (ltensor) THEN
            WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
            OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
            WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
            WRITE(10, *)
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              WRITE(10, '(a)') 'total Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' = ', frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' = ', -frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'e Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' = ', frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_e_p(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' = ', -frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_e_m(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'h Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' = ', frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_h_p(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' = ', -frq, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_h_m(irow, :, mu)
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('%', 100)
              WRITE(10, *)
              !
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ENDIF ! ltensor
          !
          ! output mode-resolved Raman intensity
          IF (TRIM(polar) == 'custom') THEN
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ELSEIF (TRIM(polar) == 'all') THEN
            ! XX
            ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! YY
            ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! XY
            ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! YX
            ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! LL
            ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! RR
            ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! LR
            ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
            ! RL
            ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
            es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'replace')
            WRITE(10, '(4(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', 'I_e(p, m)', 'I_h(p, m)'
            WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, 0.d0, 0.d0, 0.d0
            DO mu = 1, nmodes
              frq = wf(mu, iq) * rydcm1
              !
              Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p(:,:,mu), ei))))**2
              Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m(:,:,mu), ei))))**2
              Raman_e_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_p(:,:,mu), ei))))**2
              Raman_e_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_e_m(:,:,mu), ei))))**2
              Raman_h_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_p(:,:,mu), ei))))**2
              Raman_h_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_h_m(:,:,mu), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_p = bewph_s(mu)  * Raman_p
              IF (lBE) Raman_m = bewph_as(mu) * Raman_m
              IF (lBE) Raman_e_p = bewph_s(mu)  * Raman_e_p
              IF (lBE) Raman_e_m = bewph_as(mu) * Raman_e_m
              IF (lBE) Raman_h_p = bewph_s(mu)  * Raman_h_p
              IF (lBE) Raman_h_m = bewph_as(mu) * Raman_h_m
              !
              WRITE(10, '(2x, i3, 2x, f10.4, 3(1x, 2(2x, es22.12e3)))') mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                         Raman_e_p/nktotf, Raman_e_m/nktotf, Raman_h_p/nktotf, Raman_h_m/nktotf
            ENDDO ! mu
            !
            CLOSE(10)
            !
          ENDIF ! polar
          !
        ENDIF ! master node
        !
        DEALLOCATE(tensor_p, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p', 1)
        DEALLOCATE(tensor_m, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m', 1)
        DEALLOCATE(tensor_e_p, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_e_p', 1)
        DEALLOCATE(tensor_e_m, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_e_m', 1)
        DEALLOCATE(tensor_h_p, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_h_p', 1)
        DEALLOCATE(tensor_h_m, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_h_m', 1)
        !
        IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
        IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
        !
        !
      ENDIF ! ho-ph
      !
    ENDIF ! iq == 1
    !
    !################################################################
    CASE ('double')   ! double-resonance Raman 
    !################################################################
    ! not consider hole-phonon scattering
    ! *********************************************************************
    IF (.NOT. lhoph) THEN
    ! *********************************************************************
      ! allocate Raman tensor
      ALLOCATE(tensor_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_pp', 1)
      ALLOCATE(tensor_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_mm', 1)
      ALLOCATE(tensor_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_pm', 1)
      ALLOCATE(tensor_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_mp', 1)
      ! allocate Bose-Einstein occupation
      IF (lBE) ALLOCATE(bewph_as(nmodes))
      IF (lBE) ALLOCATE(bewph_s(nmodes))
      !
      tensor_pp(:, :, :, :) = czero
      tensor_mm(:, :, :, :) = czero
      tensor_pm(:, :, :, :) = czero
      tensor_mp(:, :, :, :) = czero
      !
      DO mu = 1, nmodes
        ! First phonon frequency [Ry]
        frq =  wf(mu, iq)
        ! Min phonon frequency for el-ph calculations
        IF (frq > eps_acustic) THEN
          ! Define the inverse so that we can efficiently multiply instead of dividing
          inv_frq = one / DSQRT(two * frq)
          ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
          IF (lBE) THEN
            bewph_as(mu) = wgauss(frq / temphon, -99)
            bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
            bewph_s(mu)  = bewph_as(mu) + one
          ENDIF
        ELSE
          inv_frq = zero
          IF (lBE) bewph_as(mu) = zero
          IF (lBE) bewph_s(mu) = zero
        ENDIF
        DO nu = 1, nmodes
          ! Second phonon frequency [Ry]
          frq2 = wf(nu, iq)
          IF (frq2 > eps_acustic) THEN
            inv_frq2 = one / DSQRT(two * frq2)
          ELSE
            inv_frq2 = zero
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            DO ibnd = ibndmin, ibndmax
              ! the initial state at k is sited in valence band
              IF (etf(ibnd, ikk) > ef) CYCLE
              ! Electron eigenenergy [Ry] at initial state at k
              ! note the bands window (nbndsub) defined in etf
              e_ini = etf(ibnd, ikk)
              DO jbnd = ibndmin, ibndmax
                ! the excited state at k is sited in conduction band
                IF (etf(jbnd, ikk) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at excited state at k
                e_exc = etf(jbnd, ikk)
                ! if .true., limiting the energy level seperation 
                ! for resonance consideration
                IF (reson_lim) THEN
                  IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                ! complex denominator term [Ry]
                factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
                ! dipole vector in incident light process
                ! dmef: interpolated dipole matrix elements in units of 1/Bohr
                ! note the args order and bands window (nbndsub) defined in dmef
                dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
                !
                !! Starting the double el-ph Raman process
                !
                !! emit two phonons (definitely stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 >= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser + frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_pp(:, :, mu, nu) = tensor_pp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb two phonons (definitely anti-stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 <= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser - frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_mm(:, :, mu, nu) = tensor_mm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! emit the 1st phonon and absorb the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 <= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser + frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_pm(:, :, mu, nu) = tensor_pm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb the 1st phonon and emit the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 >= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser - frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_mp(:, :, mu, nu) = tensor_mp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! End the double el-ph Raman process
                !
              ENDDO  ! jbnd
            ENDDO  ! ibnd
          ENDDO  ! ik
          !
        ENDDO ! nu
      ENDDO ! mu
      !
      !
#if defined(__MPI)
      !
      ! We need quantity from all the pools
      !
      CALL mp_sum(tensor_pp, inter_pool_comm )
      CALL mp_sum(tensor_mm, inter_pool_comm )
      CALL mp_sum(tensor_pm, inter_pool_comm )
      CALL mp_sum(tensor_mp, inter_pool_comm )
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      ! Only master writes
      IF (mpime == ionode_id) THEN
        ALLOCATE(branch(nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating branch', 1)
        ! output Raman tensor
        IF (ltensor .AND. qtensor_start <= iq .AND. iq <= qtensor_end) THEN
          !
          OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
          READ(9, *)
          DO iiq = 1, 3 * iq - 2
            READ(9, *)
          ENDDO
          READ(9, *) branch(:)
          CLOSE(9)
          !
          WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
          OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
          WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
          WRITE(10, *)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              WRITE(10, '(a)') 'total Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('%', 100)
              WRITE(10, *)
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! ltensor
        !
        ! output mode-resolved Raman intensity
        OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
        READ(9, *)
        DO iiq = 1, 3 * iq - 2
          READ(9, *)
        ENDDO
        BACKSPACE(9)
        READ(9, '(i5, 3f10.6)') iiq, q(:)
        IF (iiq /= iq) THEN
          WRITE(*, '("errors in reading the", i5, " q-point coordinates in phonon_sort file")') iq
          STOP
        ENDIF
        READ(9, *) branch(:)
        CLOSE(9)
        !
        IF (TRIM(polar) == 'custom') THEN
          INQUIRE (FILE = TRIM(filraman), EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ELSEIF (TRIM(polar) == 'all') THEN
          ! XX
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YY
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! XY
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YX
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LL
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-ll', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RR
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LR
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-lr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RL
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rl', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 1x, 4(2x, es22.12e3))') mu, nu, frq, frq2, &
                          Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! polar
        !
        DEALLOCATE(branch, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating branch', 1)
        !
      ENDIF ! master node
      !
      DEALLOCATE(tensor_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_pp', 1)
      DEALLOCATE(tensor_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_mm', 1)
      DEALLOCATE(tensor_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_pm', 1)
      DEALLOCATE(tensor_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_mp', 1)
      !
      IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
      IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
      !
      !
      ! consider the hole-phonon scattering
    ! *********************************************************************
    ELSEIF (lhoph) THEN
    ! *********************************************************************
      ! allocate Raman tensor
      ALLOCATE(tensor_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_pp', 1)
      ALLOCATE(tensor_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_mm', 1)
      ALLOCATE(tensor_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_pm', 1)
      ALLOCATE(tensor_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_mp', 1)
      !
      ALLOCATE(tensor_ee_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_pp', 1)
      ALLOCATE(tensor_ee_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_mm', 1)
      ALLOCATE(tensor_ee_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_pm', 1)
      ALLOCATE(tensor_ee_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_mp', 1)
      !
      ALLOCATE(tensor_hh_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_pp', 1)
      ALLOCATE(tensor_hh_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_mm', 1)
      ALLOCATE(tensor_hh_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_pm', 1)
      ALLOCATE(tensor_hh_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_mp', 1)
      !
      ALLOCATE(tensor_eh_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_pp', 1)
      ALLOCATE(tensor_eh_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_mm', 1)
      ALLOCATE(tensor_eh_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_pm', 1)
      ALLOCATE(tensor_eh_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_mp', 1)
      !
      ALLOCATE(tensor_he_pp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_pp', 1)
      ALLOCATE(tensor_he_mm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_mm', 1)
      ALLOCATE(tensor_he_pm(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_pm', 1)
      ALLOCATE(tensor_he_mp(3, 3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_mp', 1)
      !
      ! allocate Bose-Einstein occupation
      IF (lBE) ALLOCATE(bewph_as(nmodes))
      IF (lBE) ALLOCATE(bewph_s(nmodes))
      !
      tensor_pp(:, :, :, :) = czero
      tensor_mm(:, :, :, :) = czero
      tensor_pm(:, :, :, :) = czero
      tensor_mp(:, :, :, :) = czero
      tensor_ee_pp(:, :, :, :) = czero
      tensor_ee_mm(:, :, :, :) = czero
      tensor_ee_pm(:, :, :, :) = czero
      tensor_ee_mp(:, :, :, :) = czero
      tensor_hh_pp(:, :, :, :) = czero
      tensor_hh_mm(:, :, :, :) = czero
      tensor_hh_pm(:, :, :, :) = czero
      tensor_hh_mp(:, :, :, :) = czero
      tensor_eh_pp(:, :, :, :) = czero
      tensor_eh_mm(:, :, :, :) = czero
      tensor_eh_pm(:, :, :, :) = czero
      tensor_eh_mp(:, :, :, :) = czero
      tensor_he_pp(:, :, :, :) = czero
      tensor_he_mm(:, :, :, :) = czero
      tensor_he_pm(:, :, :, :) = czero
      tensor_he_mp(:, :, :, :) = czero
      !
      DO mu = 1, nmodes
        ! First phonon frequency [Ry]
        frq =  wf(mu, iq)
        ! Min phonon frequency for e-p calculations
        IF (frq > eps_acustic) THEN
          ! Define the inverse so that we can efficiently multiply instead of dividing
          inv_frq = one / DSQRT(two * frq)
          ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
          IF (lBE) THEN
            bewph_as(mu) = wgauss(frq / temphon, -99)
            bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
            bewph_s(mu)  = bewph_as(mu) + one
          ENDIF
        ELSE
          inv_frq = zero
          IF (lBE) bewph_as(mu) = zero
          IF (lBE) bewph_s(mu) = zero
        ENDIF
        DO nu = 1, nmodes
          ! Second phonon frequency [Ry]
          frq2 = wf(nu, iq)
          IF (frq2 > eps_acustic) THEN
            inv_frq2 = one / DSQRT(two * frq2)
          ELSE
            inv_frq2 = zero
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            DO ibnd = ibndmin, ibndmax
              ! the initial state at k is sited in valence band
              IF (etf(ibnd, ikk) > ef) CYCLE
              ! Electron eigenenergy [Ry] at initial state at k
              ! note the bands window (nbndsub) defined in etf
              e_ini = etf(ibnd, ikk)
              DO jbnd = ibndmin, ibndmax
                ! the excited state at k is sited in conduction band
                IF (etf(jbnd, ikk) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at excited state at k
                e_exc = etf(jbnd, ikk)
                ! if .true., limiting the energy level seperation for resonance consideration
                IF (reson_lim) THEN
                  IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                ! complex denominator term [Ry]
                factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
                ! dipole vector in incident light process
                ! dmef: interpolated dipole matrix elements in units of 1/Bohr
                ! note the args order and bands window (nbndsub) defined in dmef
                dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
                !
                !! Starting the double el-ph Raman process
                !
                !! emit two phonons (definitely stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 >= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser + frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_ee_pp(:, :, mu, nu) = tensor_ee_pp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb two phonons (definitely anti-stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 <= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser - frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_ee_mm(:, :, mu, nu) = tensor_ee_mm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! emit the 1st phonon and absorb the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int >= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 <= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser + frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_ee_pm(:, :, mu, nu) = tensor_ee_pm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb the 1st phonon and emit the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  ! IF (e_int <= e_exc) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in conduction band
                    IF (etf(mbnd, ikk) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = etf(mbnd, ikk)
                    ! IF (e_int2 >= e_int) CYCLE
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 - e_ini - Elaser - frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double el-ph Raman tensor
                    tensor_ee_mp(:, :, mu, nu) = tensor_ee_mp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! End the double el-ph Raman process
                !
                !
                !! Starting the double ho-ph Raman process
                !
                !! emit two phonons (definitely stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in valence band
                    IF (etf(mbnd, ikk) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = -etf(mbnd, ikk)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_exc - Elaser + frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_exc - Elaser + frq + frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double ho-ph Raman tensor
                    tensor_hh_pp(:, :, mu, nu) = tensor_hh_pp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb two phonons (definitely anti-stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in valence band
                    IF (etf(mbnd, ikk) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = -etf(mbnd, ikk)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_exc - Elaser - frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_exc - Elaser - frq - frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double ho-ph Raman tensor
                    tensor_hh_mm(:, :, mu, nu) = tensor_hh_mm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! emit the 1st phonon and absorb the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in valence band
                    IF (etf(mbnd, ikk) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = -etf(mbnd, ikk)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_exc - Elaser + frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_exc - Elaser + frq - frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double ho-ph Raman tensor
                    tensor_hh_pm(:, :, mu, nu) = tensor_hh_pm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb the 1st phonon and emit the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state back to k is sited in valence band
                    IF (etf(mbnd, ikk) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                    e_int2 = -etf(mbnd, ikk)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_exc - Elaser - frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_exc - Elaser - frq + frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! double ho-ph Raman tensor
                    tensor_hh_mp(:, :, mu, nu) = tensor_hh_mp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! End the double ho-ph Raman process
                !
                !
                !! Starting the el-ph & ho-ph Raman process
                !
                !! emit two phonons (definitely stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in valence band
                    IF (etf(mbnd, ikq) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = -etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser + frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser + frq + frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! el-ph & ho-ph Raman tensor
                    tensor_eh_pp(:, :, mu, nu) = tensor_eh_pp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb two phonons (definitely anti-stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in valence band
                    IF (etf(mbnd, ikq) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = -etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser - frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! el-ph & ho-ph Raman tensor
                    tensor_eh_mm(:, :, mu, nu) = tensor_eh_mm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! emit the 1st phonon and absorb the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in valence band
                    IF (etf(mbnd, ikq) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = -etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser + frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! el-ph & ho-ph Raman tensor
                    tensor_eh_pm(:, :, mu, nu) = tensor_eh_pm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb the 1st phonon and emit the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in conduction band
                  IF (etf(kbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                  ! el-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in valence band
                    IF (etf(mbnd, ikq) > ef) CYCLE
                    ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = -etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser - frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser - frq + frq2 - ci*Egamma)
                    ! ho-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, nu, ik) * inv_frq2
                    elphme2 = CONJG(elphme2)
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! el-ph & ho-ph Raman tensor
                    tensor_eh_mp(:, :, mu, nu) = tensor_eh_mp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! End the el-ph & ho-ph Raman process
                !
                !
                !! Starting the ho-ph & el-ph Raman process
                !
                !! emit two phonons (definitely stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry] 
                  ! note the args order and bands window (nbndfst) defined in epf17 
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in conduction band
                    IF (etf(mbnd, ikq) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser + frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser + frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! ho-ph & el-ph Raman tensor
                    tensor_he_pp(:, :, mu, nu) = tensor_he_pp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb two phonons (definitely anti-stokes Raman)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in conduction band
                    IF (etf(mbnd, ikq) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser - frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! ho-ph & el-ph Raman tensor
                    tensor_he_mm(:, :, mu, nu) = tensor_he_mm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! emit the 1st phonon and absorb the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in conduction band
                    IF (etf(mbnd, ikq) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser + frq - frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! ho-ph & el-ph Raman tensor
                    tensor_he_pm(:, :, mu, nu) = tensor_he_pm(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! absorb the 1st phonon and emit the 2nd phonon (may or may not stokes/anti-stokes)
                !
                DO kbnd = ibndmin, ibndmax
                  ! the 1st intermediate state at k+q is sited in valence band
                  IF (etf(kbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                  e_int = -etf(kbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-ph coupling with first phonon [Ry]
                  elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme = CONJG(elphme)
                  DO mbnd = ibndmin, ibndmax
                    ! the 2nd intermediate state at k+q is sited in conduction band
                    IF (etf(mbnd, ikq) <= ef) CYCLE
                    ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                    e_int2 = etf(mbnd, ikq)
                    IF (reson_lim) THEN
                      IF (ABS(e_int2 + e_int - Elaser - frq + frq2) > reson_thr) CYCLE
                    ENDIF
                    factor_ep2 = cone / (e_int2 + e_int - Elaser - frq + frq2 - ci*Egamma)
                    ! el-ph coupling with second phonon [Ry]
                    elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, nu, ik) * inv_frq2
                    ! dipole vector in scattered light process
                    dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                    !
                    ! temporary Raman tensor
                    tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                                 * elphme * elphme2 * factor * factor_ep * factor_ep2
                    ! ho-ph & el-ph Raman tensor
                    tensor_he_mp(:, :, mu, nu) = tensor_he_mp(:, :, mu, nu) + tensor_tem(:, :)
                  ENDDO ! mbnd
                ENDDO ! kbnd
                !
                !! End the ho-ph & el-ph Raman process
                !
                !
              ENDDO  ! jbnd
            ENDDO  ! ibnd
          ENDDO  ! ik
          !
        ENDDO ! nu
      ENDDO ! mu
      !
      tensor_pp(:, :, :, :) = tensor_ee_pp(:, :, :, :) + tensor_hh_pp(:, :, :, :) + &
                              tensor_eh_pp(:, :, :, :) + tensor_he_pp(:, :, :, :)
      tensor_mm(:, :, :, :) = tensor_ee_mm(:, :, :, :) + tensor_hh_mm(:, :, :, :) + &
                              tensor_eh_mm(:, :, :, :) + tensor_he_mm(:, :, :, :)
      tensor_pm(:, :, :, :) = tensor_ee_pm(:, :, :, :) + tensor_hh_pm(:, :, :, :) + &
                              tensor_eh_pm(:, :, :, :) + tensor_he_pm(:, :, :, :)
      tensor_mp(:, :, :, :) = tensor_ee_mp(:, :, :, :) + tensor_hh_mp(:, :, :, :) + &
                              tensor_eh_mp(:, :, :, :) + tensor_he_mp(:, :, :, :)
      !
      !
#if defined(__MPI)
      !
      ! We need quantity from all the pools
      !
      CALL mp_sum(tensor_pp, inter_pool_comm )
      CALL mp_sum(tensor_mm, inter_pool_comm )
      CALL mp_sum(tensor_pm, inter_pool_comm )
      CALL mp_sum(tensor_mp, inter_pool_comm )
      CALL mp_sum(tensor_ee_pp, inter_pool_comm )
      CALL mp_sum(tensor_ee_mm, inter_pool_comm )
      CALL mp_sum(tensor_ee_pm, inter_pool_comm )
      CALL mp_sum(tensor_ee_mp, inter_pool_comm )
      CALL mp_sum(tensor_hh_pp, inter_pool_comm )
      CALL mp_sum(tensor_hh_mm, inter_pool_comm )
      CALL mp_sum(tensor_hh_pm, inter_pool_comm )
      CALL mp_sum(tensor_hh_mp, inter_pool_comm )
      CALL mp_sum(tensor_eh_pp, inter_pool_comm )
      CALL mp_sum(tensor_eh_mm, inter_pool_comm )
      CALL mp_sum(tensor_eh_pm, inter_pool_comm )
      CALL mp_sum(tensor_eh_mp, inter_pool_comm )
      CALL mp_sum(tensor_he_pp, inter_pool_comm )
      CALL mp_sum(tensor_he_mm, inter_pool_comm )
      CALL mp_sum(tensor_he_pm, inter_pool_comm )
      CALL mp_sum(tensor_he_mp, inter_pool_comm )
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      ! Only master writes
      IF (mpime == ionode_id) THEN
        ALLOCATE(branch(nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating branch', 1)
        ! output Raman tensor
        IF (ltensor .AND. qtensor_start <= iq .AND. iq <= qtensor_end) THEN
          !
          OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
          READ(9, *)
          DO iiq = 1, 3 * iq - 2
            READ(9, *)
          ENDDO
          READ(9, *) branch(:)
          CLOSE(9)
          !
          WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
          OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
          WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
          WRITE(10, *)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              WRITE(10, '(a)') 'total Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'ee Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'hh Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'eh Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('*', 100)
              !
              WRITE(10, '(a)') 'he Raman tensor:'
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+', nu, ' = ', frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_pp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '-', nu, ' = ', -frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_mm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band +', mu, '-', nu, ' = ', frq-frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_pm(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a, i3.3, a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+', nu, ' = ', -frq+frq2, ' cm-1'
              DO irow = 1, 3
                WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_mp(irow, :, branch(mu), branch(nu))
              ENDDO
              WRITE(10, *)
              WRITE(10, '(a/)') REPEAT('%', 100)
              WRITE(10, *)
              !
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! ltensor
        !
        ! output mode-resolved Raman intensity
        OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
        READ(9, *)
        DO iiq = 1, 3 * iq - 2
          READ(9, *)
        ENDDO
        BACKSPACE(9)
        READ(9, '(i5, 3f10.6)') iiq, q(:)
        IF (iiq /= iq) THEN
          WRITE(*, '("errors in reading the", i5, " q-point coordinates in phonon_sort file")') iq
          STOP
        ENDIF
        READ(9, *) branch(:)
        CLOSE(9)
        !
        IF (TRIM(polar) == 'custom') THEN
          INQUIRE (FILE = TRIM(filraman), EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ELSEIF (TRIM(polar) == 'all') THEN
          ! XX
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YY
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! XY
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YX
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LL
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-ll', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RR
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LR
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-lr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RL
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rl', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(pp, mm, pm, mp)', &
                                        'I_ee(pp, mm, pm, mp)', 'I_hh(pp, mm, pm, mp)', &
                                        'I_eh(pp, mm, pm, mp)', 'I_he(pp, mm, pm, mp)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            DO nu = 1, nmodes
              frq2 = wf(branch(nu), iq) * rydcm1
              !
              Raman_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_ee_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_hh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_eh_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_mp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pp(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_pm = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_pm(:,:,branch(mu),branch(nu)), ei))))**2
              Raman_he_mp = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_mp(:,:,branch(mu),branch(nu)), ei))))**2
              !
              ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
              IF (lBE) Raman_pp = bewph_s(branch(mu))  * Raman_pp * bewph_s(branch(nu))
              IF (lBE) Raman_mm = bewph_as(branch(mu)) * Raman_mm * bewph_as(branch(nu))
              IF (lBE) Raman_pm = bewph_s(branch(mu))  * Raman_pm * bewph_as(branch(nu))
              IF (lBE) Raman_mp = bewph_as(branch(mu)) * Raman_mp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_pp = bewph_s(branch(mu))  * Raman_ee_pp * bewph_s(branch(nu))
              IF (lBE) Raman_ee_mm = bewph_as(branch(mu)) * Raman_ee_mm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_pm = bewph_s(branch(mu))  * Raman_ee_pm * bewph_as(branch(nu))
              IF (lBE) Raman_ee_mp = bewph_as(branch(mu)) * Raman_ee_mp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_pp = bewph_s(branch(mu))  * Raman_hh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_hh_mm = bewph_as(branch(mu)) * Raman_hh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_pm = bewph_s(branch(mu))  * Raman_hh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_hh_mp = bewph_as(branch(mu)) * Raman_hh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_pp = bewph_s(branch(mu))  * Raman_eh_pp * bewph_s(branch(nu))
              IF (lBE) Raman_eh_mm = bewph_as(branch(mu)) * Raman_eh_mm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_pm = bewph_s(branch(mu))  * Raman_eh_pm * bewph_as(branch(nu))
              IF (lBE) Raman_eh_mp = bewph_as(branch(mu)) * Raman_eh_mp * bewph_s(branch(nu))
              IF (lBE) Raman_he_pp = bewph_s(branch(mu))  * Raman_he_pp * bewph_s(branch(nu))
              IF (lBE) Raman_he_mm = bewph_as(branch(mu)) * Raman_he_mm * bewph_as(branch(nu))
              IF (lBE) Raman_he_pm = bewph_s(branch(mu))  * Raman_he_pm * bewph_as(branch(nu))
              IF (lBE) Raman_he_mp = bewph_as(branch(mu)) * Raman_he_mp * bewph_s(branch(nu))
              !
              WRITE(10, '(2(2x, i3), 2(2x, f10.4), 5(1x,4(2x, es22.12e3)))') mu, nu, frq, frq2, &
                            Raman_pp/nktotf, Raman_mm/nktotf, Raman_pm/nktotf, Raman_mp/nktotf, &
                Raman_ee_pp/nktotf, Raman_ee_mm/nktotf, Raman_ee_pm/nktotf, Raman_ee_mp/nktotf, &
                Raman_hh_pp/nktotf, Raman_hh_mm/nktotf, Raman_hh_pm/nktotf, Raman_hh_mp/nktotf, &
                Raman_eh_pp/nktotf, Raman_eh_mm/nktotf, Raman_eh_pm/nktotf, Raman_eh_mp/nktotf, &
                Raman_he_pp/nktotf, Raman_he_mm/nktotf, Raman_he_pm/nktotf, Raman_he_mp/nktotf
            ENDDO ! nu
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! polar
        !
        DEALLOCATE(branch, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating branch', 1)
        !
      ENDIF ! master node
      !
      DEALLOCATE(tensor_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_pp', 1)
      DEALLOCATE(tensor_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_mm', 1)
      DEALLOCATE(tensor_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_pm', 1)
      DEALLOCATE(tensor_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_mp', 1)
      DEALLOCATE(tensor_ee_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_pp', 1)
      DEALLOCATE(tensor_ee_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_mm', 1)
      DEALLOCATE(tensor_ee_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_pm', 1)
      DEALLOCATE(tensor_ee_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_mp', 1)
      DEALLOCATE(tensor_hh_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_pp', 1)
      DEALLOCATE(tensor_hh_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_mm', 1)
      DEALLOCATE(tensor_hh_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_pm', 1)
      DEALLOCATE(tensor_hh_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_mp', 1)
      DEALLOCATE(tensor_eh_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_pp', 1)
      DEALLOCATE(tensor_eh_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_mm', 1)
      DEALLOCATE(tensor_eh_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_pm', 1)
      DEALLOCATE(tensor_eh_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_mp', 1)
      DEALLOCATE(tensor_he_pp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_pp', 1)
      DEALLOCATE(tensor_he_mm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_mm', 1)
      DEALLOCATE(tensor_he_pm, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_pm', 1)
      DEALLOCATE(tensor_he_mp, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_mp', 1)
      !
      IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
      IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
      !
      !
    ENDIF ! ho-ph
    !
    !################################################################
    CASE ('defect')   ! defect-induced double-resonance Raman
    !################################################################
    ! dependence of defect-electron coupling on the q-wavevector
    OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
    READ(9, *)
    DO iiq = 1, 3 * iq - 2
      READ(9, *)
    ENDDO
    BACKSPACE(9)
    READ(9, '(i5, 3f10.6)') iiq, q(:)
    IF (iiq /= iq) THEN
      WRITE(*, '("errors in reading the", i5, " q-point coordinates in phonon_sort file")') iq
      STOP
    ENDIF
    CLOSE(9)
    ql = SQRT(q(1)**2 + q(2)**2 + q(3)**2)
    IF (ql < 1.d-6) ql = 1.d-6
    Ide = (Cq(1)*ql**(-3.d0) + Cq(2)*ql**(-2.d0) + Cq(3)*ql**(-1.d0) + &
           Cq(4) + Cq(5)*ql + Cq(6)*ql**(2.d0) + Cq(7)*ql**(3.d0))*cone
    !
    ! not consider hole-phonon scattering
    ! *********************************************************************
    IF (.NOT. lhoph) THEN
    ! *********************************************************************
      ! allocate Raman tensor
      ALLOCATE(tensor_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p1', 1)
      ALLOCATE(tensor_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p2', 1)
      ALLOCATE(tensor_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m1', 1)
      ALLOCATE(tensor_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m2', 1)
      ! allocate Bose-Einstein occupation
      IF (lBE) ALLOCATE(bewph_as(nmodes))
      IF (lBE) ALLOCATE(bewph_s(nmodes))
      !
      tensor_p1(:, :, :) = czero
      tensor_p2(:, :, :) = czero
      tensor_m1(:, :, :) = czero
      tensor_m2(:, :, :) = czero
      !
      DO mu = 1, nmodes
        ! First phonon frequency [Ry]
        frq =  wf(mu, iq)
        ! Min phonon frequency for el-ph calculations
        IF (frq > eps_acustic) THEN
          ! Define the inverse so that we can efficiently multiply instead of dividing
          inv_frq = one / DSQRT(two * frq)
          ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
          IF (lBE) THEN
            bewph_as(mu) = wgauss(frq / temphon, -99)
            bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
            bewph_s(mu)  = bewph_as(mu) + one
          ENDIF
        ELSE
          inv_frq = zero
          IF (lBE) bewph_as(mu) = zero
          IF (lBE) bewph_s(mu) = zero
        ENDIF
        !
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          DO ibnd = ibndmin, ibndmax
            ! the initial state at k is sited in valence band
            IF (etf(ibnd, ikk) > ef) CYCLE
            ! Electron eigenenergy [Ry] at initial state at k
            ! note the bands window (nbndsub) defined in etf
            e_ini = etf(ibnd, ikk)
            DO jbnd = ibndmin, ibndmax
              ! the excited state at k is sited in conduction band
              IF (etf(jbnd, ikk) <= ef) CYCLE
              ! Electron eigenenergy [Ry] at excited state at k
              e_exc = etf(jbnd, ikk)
              ! if .true., limiting the energy level seperation 
              ! for resonance consideration
              IF (reson_lim) THEN
                IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
              ENDIF
              ! complex denominator term [Ry]
              factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
              ! dipole vector in incident light process
              ! dmef: interpolated dipole matrix elements in units of 1/Bohr
              ! note the args order and bands window (nbndsub) defined in dmef
              dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
              !
              !! Starting the defect-induced double Raman process
              !
              !! emit one phonon (definitely stokes Raman)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int >= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 >= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! defect-induced double resonance Raman tensor
                  tensor_p1(:, :, mu) = tensor_p1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! emit one phonon (definitely stokes Raman)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int >= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-defect coupling in the first scattering [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 >= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling with phonon [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! defect-induced double resonance Raman tensor
                  tensor_p2(:, :, mu) = tensor_p2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (definitely anti-stokes Raman)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int <= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 <= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! defect-induced double resonance Raman tensor
                  tensor_m1(:, :, mu) = tensor_m1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (definitely anti-stokes Raman)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int <= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 <= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! defect-induced double resonance Raman tensor
                  tensor_m2(:, :, mu) = tensor_m2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! End the double resonance Raman process
              !
            ENDDO  ! jbnd
          ENDDO  ! ibnd
        ENDDO  ! ik
        !
      ENDDO ! mu
      !
      !
#if defined(__MPI)
      !
      ! We need quantity from all the pools
      !
      CALL mp_sum(tensor_p1, inter_pool_comm )
      CALL mp_sum(tensor_p2, inter_pool_comm )
      CALL mp_sum(tensor_m1, inter_pool_comm )
      CALL mp_sum(tensor_m2, inter_pool_comm )
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      ! Only master writes
      IF (mpime == ionode_id) THEN
        ALLOCATE(branch(nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating branch', 1)
        ! output Raman tensor
        IF (ltensor .AND. qtensor_start <= iq .AND. iq <= qtensor_end) THEN
          !
          OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
          READ(9, *)
          DO iiq = 1, 3 * iq - 2
            READ(9, *)
          ENDDO
          READ(9, *) branch(:)
          CLOSE(9)
          !
          WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
          OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
          WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
          WRITE(10, *)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            WRITE(10, '(a)') 'total Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, '+ 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, '+ 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('%', 100)
            WRITE(10, *)
            !
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! ltensor
        !
        ! output mode-resolved Raman intensity
        OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
        READ(9, *)
        DO iiq = 1, 3 * iq - 2
          READ(9, *)
        ENDDO
        BACKSPACE(9)
        READ(9, '(i5, 3f10.6)') iiq, q(:)
        IF (iiq /= iq) THEN
          WRITE(*, '("errors in reading the", i5, " q-point coordinates in phonon_sort file")') iq
          STOP
        ENDIF
        READ(9, *) branch(:)
        CLOSE(9)
        !
        IF (TRIM(polar) == 'custom') THEN
          INQUIRE (FILE = TRIM(filraman), EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ELSEIF (TRIM(polar) == 'all') THEN
          ! XX
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YY
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! XY
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YX
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LL
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-ll', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RR
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LR
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-lr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RL
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rl', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'new')
            WRITE(10, '(2(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          ! first phonon scattering then defect scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          ! first defect scattering then phonon scattering
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 1x, 2(2x, es22.12e3))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! polar
        !
        DEALLOCATE(branch, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating branch', 1)
        !
      ENDIF ! master node
      !
      DEALLOCATE(tensor_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p1', 1)
      DEALLOCATE(tensor_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p2', 1)
      DEALLOCATE(tensor_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m1', 1)
      DEALLOCATE(tensor_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m2', 1)
      !
      IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
      IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
      !
      !
      ! consider the hole-phonon scattering
    ! *********************************************************************
    ELSEIF (lhoph) THEN
    ! *********************************************************************
      ! allocate Raman tensor
      ALLOCATE(tensor_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p1', 1)
      ALLOCATE(tensor_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_p2', 1)
      ALLOCATE(tensor_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m1', 1)
      ALLOCATE(tensor_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_m2', 1)
      !
      ALLOCATE(tensor_ee_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_p1', 1)
      ALLOCATE(tensor_ee_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_p2', 1)
      ALLOCATE(tensor_ee_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_m1', 1)
      ALLOCATE(tensor_ee_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_ee_m2', 1)
      !
      ALLOCATE(tensor_hh_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_p1', 1)
      ALLOCATE(tensor_hh_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_p2', 1)
      ALLOCATE(tensor_hh_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_m1', 1)
      ALLOCATE(tensor_hh_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_hh_m2', 1)
      !
      ALLOCATE(tensor_eh_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_p1', 1)
      ALLOCATE(tensor_eh_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_p2', 1)
      ALLOCATE(tensor_eh_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_m1', 1)
      ALLOCATE(tensor_eh_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_eh_m2', 1)
      !
      ALLOCATE(tensor_he_p1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_p1', 1)
      ALLOCATE(tensor_he_p2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_p2', 1)
      ALLOCATE(tensor_he_m1(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_m1', 1)
      ALLOCATE(tensor_he_m2(3, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating tensor_he_m2', 1)
      !
      ! allocate Bose-Einstein occupation
      IF (lBE) ALLOCATE(bewph_as(nmodes))
      IF (lBE) ALLOCATE(bewph_s(nmodes))
      !
      tensor_p1(:, :, :) = czero
      tensor_p2(:, :, :) = czero
      tensor_m1(:, :, :) = czero
      tensor_m2(:, :, :) = czero
      tensor_ee_p1(:, :, :) = czero
      tensor_ee_p2(:, :, :) = czero
      tensor_ee_m1(:, :, :) = czero
      tensor_ee_m2(:, :, :) = czero
      tensor_hh_p1(:, :, :) = czero
      tensor_hh_p2(:, :, :) = czero
      tensor_hh_m1(:, :, :) = czero
      tensor_hh_m2(:, :, :) = czero
      tensor_eh_p1(:, :, :) = czero
      tensor_eh_p2(:, :, :) = czero
      tensor_eh_m1(:, :, :) = czero
      tensor_eh_m2(:, :, :) = czero
      tensor_he_p1(:, :, :) = czero
      tensor_he_p2(:, :, :) = czero
      tensor_he_m1(:, :, :) = czero
      tensor_he_m2(:, :, :) = czero
      !
      DO mu = 1, nmodes
        ! First phonon frequency [Ry]
        frq =  wf(mu, iq)
        ! Min phonon frequency for e-p calculations
        IF (frq > eps_acustic) THEN
          ! Define the inverse so that we can efficiently multiply instead of dividing
          inv_frq = one / DSQRT(two * frq)
          ! Set Bose-Einstein occupation of stokes/anti-stokes phonons
          IF (lBE) THEN
            bewph_as(mu) = wgauss(frq / temphon, -99)
            bewph_as(mu) = half / (two * bewph_as(mu) - one) - half
            bewph_s(mu)  = bewph_as(mu) + one
          ENDIF
        ELSE
          inv_frq = zero
          IF (lBE) bewph_as(mu) = zero
          IF (lBE) bewph_s(mu) = zero
        ENDIF
          !
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          DO ibnd = ibndmin, ibndmax
            ! the initial state at k is sited in valence band
            IF (etf(ibnd, ikk) > ef) CYCLE
            ! Electron eigenenergy [Ry] at initial state at k
            ! note the bands window (nbndsub) defined in etf
            e_ini = etf(ibnd, ikk)
            DO jbnd = ibndmin, ibndmax
              ! the excited state at k is sited in conduction band
              IF (etf(jbnd, ikk) <= ef) CYCLE
              ! Electron eigenenergy [Ry] at excited state at k
              e_exc = etf(jbnd, ikk)
              ! if .true., limiting the energy level seperation for resonance consideration
              IF (reson_lim) THEN
                IF (ABS(e_exc - e_ini - Elaser) > reson_thr) CYCLE
              ENDIF
              ! complex denominator term [Ry]
              factor = cone / (e_exc - e_ini - Elaser - ci*Egamma)
              ! dipole vector in incident light process
              ! dmef: interpolated dipole matrix elements in units of 1/Bohr
              ! note the args order and bands window (nbndsub) defined in dmef
              dipole_inc(:) = dmef(:, jbnd, ibnd, ikk)
              !
              !! Starting the ee Raman process
              !
              !! emit one phonon (first el-ph and then el-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int >= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 >= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph Raman tensor
                  tensor_ee_p1(:, :, mu) = tensor_ee_p1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! emit one phonon (first el-defect and then el-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int >= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-defect coupling in the first scattering [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 >= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser + frq - ci*Egamma)
                  ! el-ph coupling in the second scattering [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph Raman tensor
                  tensor_ee_p2(:, :, mu) = tensor_ee_p2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first el-ph and then el-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int <= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 <= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph and el-defect Raman tensor
                  tensor_ee_m1(:, :, mu) = tensor_ee_m1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first el-defect and then el-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                ! IF (e_int <= e_exc) CYCLE
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in conduction band
                  IF (etf(mbnd, ikk) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = etf(mbnd, ikk)
                  ! IF (e_int2 <= e_int) CYCLE
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 - e_ini - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 - e_ini - Elaser - frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, ibnd, mbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph and el-defect Raman tensor
                  tensor_ee_m2(:, :, mu) = tensor_ee_m2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! Starting the hh Raman process
              !
              !! emit one phonon (first ho-ph and then ho-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                elphme = CONJG(elphme)
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in valence band
                  IF (etf(mbnd, ikk) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = -etf(mbnd, ikk)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph and ho-defect Raman tensor
                  tensor_hh_p1(:, :, mu) = tensor_hh_p1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! emit one phonon (first ho-defect and then ho-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in valence band
                  IF (etf(mbnd, ikk) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = -etf(mbnd, ikk)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_exc - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_exc - Elaser + frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph and ho-defect Raman tensor
                  tensor_hh_p2(:, :, mu) = tensor_hh_p2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first ho-ph and then ho-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry]
                elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                elphme = CONJG(elphme)
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in valence band
                  IF (etf(mbnd, ikk) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = -etf(mbnd, ikk)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph and ho-defect Raman tensor
                  tensor_hh_m1(:, :, mu) = tensor_hh_m1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first ho-defect and then ho-phonon)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry]
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state back to k is sited in valence band
                  IF (etf(mbnd, ikk) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state back to k
                  e_int2 = -etf(mbnd, ikk)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_exc - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_exc - Elaser - frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(kbnd - ibndmin + 1, mbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, jbnd, ikk)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph and ho-defect Raman tensor
                  tensor_hh_m2(:, :, mu) = tensor_hh_m2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! Starting the eh Raman process
              !
              !! emit one phonon (first el-ph and then ho-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser + frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser + frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in valence band
                  IF (etf(mbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = -etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph & ho-defect Raman tensor
                  tensor_eh_p1(:, :, mu) = tensor_eh_p1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! emit one phonon (first el-defect and then ho-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-defect coupling in the first scattering [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in valence band
                  IF (etf(mbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = -etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - ci*Egamma)
                  ! ho-ph coupling in the second scattering [Ry]
                  elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph & ho-defect Raman tensor
                  tensor_eh_p2(:, :, mu) = tensor_eh_p2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first el-ph and then ho-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser - frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - frq - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = epf17(kbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in valence band
                  IF (etf(mbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = -etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph & ho-defect Raman tensor
                  tensor_eh_m1(:, :, mu) = tensor_eh_m1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first el-defect and then ho-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in conduction band
                IF (etf(kbnd, ikq) <= ef) CYCLE
                ! Electron eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int - e_ini - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int - e_ini - Elaser - ci*Egamma)
                ! el-ph coupling with first phonon [Ry]
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in valence band
                  IF (etf(mbnd, ikq) > ef) CYCLE
                  ! Hole eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = -etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - ci*Egamma)
                  ! ho-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(mbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                  elphme2 = CONJG(elphme2)
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, mbnd, kbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! el-ph & ho-defect Raman tensor
                  tensor_eh_m2(:, :, mu) = tensor_eh_m2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! Starting the he Raman process
              !
              !! emit one phonon (first ho-ph and then el-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser + frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser + frq - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                elphme = CONJG(elphme)
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in conduction band
                  IF (etf(mbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph & el-defect Raman tensor
                  tensor_he_p1(:, :, mu) = tensor_he_p1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! emit one phonon (first ho-defect and then el-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry] 
                ! note the args order and bands window (nbndfst) defined in epf17 
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in conduction band
                  IF (etf(mbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser + frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser + frq - ci*Egamma)
                  ! el-defect coupling in the second scattering [Ry]
                  elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph & el-defect Raman tensor
                  tensor_he_p2(:, :, mu) = tensor_he_p2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first ho-ph and then el-defect)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser - frq) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - frq - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry]
                elphme = epf17(kbnd - ibndmin + 1, ibnd - ibndmin + 1, mu, ik) * inv_frq
                elphme = CONJG(elphme)
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in conduction band
                  IF (etf(mbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - ci*Egamma)
                  ! el-defect couplingin the second scattering [Ry]
                  elphme2 = Ide
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph & el-defect Raman tensor
                  tensor_he_m1(:, :, mu) = tensor_he_m1(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
              !! absorb one phonon (first ho-defect and then el-ph)
              !
              DO kbnd = ibndmin, ibndmax
                ! the 1st intermediate state at k+q is sited in valence band
                IF (etf(kbnd, ikq) > ef) CYCLE
                ! Hole eigenenergy [Ry] at 1st intermediate state at k+q
                e_int = -etf(kbnd, ikq)
                IF (reson_lim) THEN
                  IF (ABS(e_int + e_exc - Elaser) > reson_thr) CYCLE
                ENDIF
                factor_ep = cone / (e_int + e_exc - Elaser - ci*Egamma)
                ! ho-ph coupling with first phonon [Ry]
                elphme = Ide
                DO mbnd = ibndmin, ibndmax
                  ! the 2nd intermediate state at k+q is sited in conduction band
                  IF (etf(mbnd, ikq) <= ef) CYCLE
                  ! Electron eigenenergy [Ry] at 2nd intermediate state at k+q
                  e_int2 = etf(mbnd, ikq)
                  IF (reson_lim) THEN
                    IF (ABS(e_int2 + e_int - Elaser - frq) > reson_thr) CYCLE
                  ENDIF
                  factor_ep2 = cone / (e_int2 + e_int - Elaser - frq - ci*Egamma)
                  ! el-defect couplingin the second scattering [Ry]
                  elphme2 = epf17(mbnd - ibndmin + 1, jbnd - ibndmin + 1, mu, ik) * inv_frq
                  ! dipole vector in scattered light process
                  dipole_sca(:) = dmef(:, kbnd, mbnd, ikq)
                  !
                  ! temporary Raman tensor
                  tensor_tem = MATMUL(RESHAPE(dipole_sca,(/3,1/)), RESHAPE(dipole_inc,(/1,3/))) &
                               * elphme * elphme2 * factor * factor_ep * factor_ep2
                  ! ho-ph & el-defect Raman tensor
                  tensor_he_m2(:, :, mu) = tensor_he_m2(:, :, mu) + tensor_tem(:, :)
                ENDDO ! mbnd
              ENDDO ! kbnd
              !
            ENDDO  ! jbnd
          ENDDO  ! ibnd
        ENDDO  ! ik
        !
      ENDDO ! mu
      !
      tensor_p1(:, :, :) = tensor_ee_p1(:, :, :) + tensor_hh_p1(:, :, :) + &
                              tensor_eh_p1(:, :, :) + tensor_he_p1(:, :, :)
      tensor_p2(:, :, :) = tensor_ee_p2(:, :, :) + tensor_hh_p2(:, :, :) + &
                              tensor_eh_p2(:, :, :) + tensor_he_p2(:, :, :)
      tensor_m1(:, :, :) = tensor_ee_m1(:, :, :) + tensor_hh_m1(:, :, :) + &
                              tensor_eh_m1(:, :, :) + tensor_he_m1(:, :, :)
      tensor_m2(:, :, :) = tensor_ee_m2(:, :, :) + tensor_hh_m2(:, :, :) + &
                              tensor_eh_m2(:, :, :) + tensor_he_m2(:, :, :)
      !
      !
#if defined(__MPI)
      !
      ! We need quantity from all the pools
      !
      CALL mp_sum(tensor_p1, inter_pool_comm )
      CALL mp_sum(tensor_p2, inter_pool_comm )
      CALL mp_sum(tensor_m1, inter_pool_comm )
      CALL mp_sum(tensor_m2, inter_pool_comm )
      CALL mp_sum(tensor_ee_p1, inter_pool_comm )
      CALL mp_sum(tensor_ee_p2, inter_pool_comm )
      CALL mp_sum(tensor_ee_m1, inter_pool_comm )
      CALL mp_sum(tensor_ee_m2, inter_pool_comm )
      CALL mp_sum(tensor_hh_p1, inter_pool_comm )
      CALL mp_sum(tensor_hh_p2, inter_pool_comm )
      CALL mp_sum(tensor_hh_m1, inter_pool_comm )
      CALL mp_sum(tensor_hh_m2, inter_pool_comm )
      CALL mp_sum(tensor_eh_p1, inter_pool_comm )
      CALL mp_sum(tensor_eh_p2, inter_pool_comm )
      CALL mp_sum(tensor_eh_m1, inter_pool_comm )
      CALL mp_sum(tensor_eh_m2, inter_pool_comm )
      CALL mp_sum(tensor_he_p1, inter_pool_comm )
      CALL mp_sum(tensor_he_p2, inter_pool_comm )
      CALL mp_sum(tensor_he_m1, inter_pool_comm )
      CALL mp_sum(tensor_he_m2, inter_pool_comm )
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      ! Only master writes
      IF (mpime == ionode_id) THEN
        ALLOCATE(branch(nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error allocating branch', 1)
        ! output Raman tensor
        IF (ltensor .AND. qtensor_start <= iq .AND. iq <= qtensor_end) THEN
          !
          OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
          READ(9, *)
          DO iiq = 1, 3 * iq - 2
            READ(9, *)
          ENDDO
          READ(9, *) branch(:)
          CLOSE(9)
          !
          WRITE(ftensor, '(a, i5.5)') 'tensor_q', iq
          OPEN(UNIT = 10, FILE = ftensor, STATUS = 'replace')
          WRITE(10, '(a)') 'The Raman tensors are listed in a branch-resolved way'
          WRITE(10, *)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            WRITE(10, '(a)') 'total Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' + 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' + 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('*', 100)
            !
            WRITE(10, '(a)') 'ee Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' + 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' + 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_ee_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('*', 100)
            !
            WRITE(10, '(a)') 'hh Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' + 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' + 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_hh_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('*', 100)
            !
            WRITE(10, '(a)') 'eh Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' + 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' + 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_eh_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('*', 100)
            !
            WRITE(10, '(a)') 'he Raman tensor:'
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band +', mu, ' + 0 = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_p1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 +', mu, ' = ', frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_p2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band -', mu, ' + 0 = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_m1(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a, i3.3, a, f10.4, a)') 'Raman band 0 -', mu, ' = ', -frq, ' cm-1'
            DO irow = 1, 3
              WRITE(10, '(3("(", es14.4e3, ",", es14.4e3, ")", 5x))') tensor_he_m2(irow, :, branch(mu))
            ENDDO
            WRITE(10, *)
            WRITE(10, '(a/)') REPEAT('%', 100)
            WRITE(10, *)
            !
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! ltensor
        !
        ! output mode-resolved Raman intensity
        OPEN(UNIT = 9, FILE = 'phonon_sort', STATUS = 'old')
        READ(9, *)
        DO iiq = 1, 3 * iq - 2
          READ(9, *)
        ENDDO
        BACKSPACE(9)
        READ(9, '(i5, 3f10.6)') iiq, q(:)
        IF (iiq /= iq) THEN
          WRITE(*, '("errors in reading the", i5, " q-point coordinates in phonon_sort file")') iq
          STOP
        ENDIF
        READ(9, *) branch(:)
        CLOSE(9)
        !
        IF (TRIM(polar) == 'custom') THEN
          INQUIRE (FILE = TRIM(filraman), EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = TRIM(filraman), STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ELSEIF (TRIM(polar) == 'all') THEN
          ! XX
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YY
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! XY
          ei = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-xy', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-xy', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! YX
          ei = (/ (0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-yx', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-yx', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LL
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-ll', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-ll', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RR
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! LR
          ei = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-lr', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-lr', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
          ! RL
          ei = (/ (1.d0,0.d0),(0.d0,-1.d0),(0.d0,0.d0) /)
          es = (/ (1.d0,0.d0),(0.d0,1.d0),(0.d0,0.d0) /)
          INQUIRE (FILE = 'qraman-rl', EXIST = alive)
          IF (alive) THEN
            CONTINUE
          ELSE
            OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'new')
            WRITE(10, '(6(a, 10x), a)') '#phonon_branch', 'branch_frequency', 'I_tot(p, m)', &
                                     'I_ee(p, m)', 'I_hh(p, m)', 'I_eh(p, m)', 'I_he(p, m)'
            CLOSE(10)
          ENDIF
          !
          OPEN(UNIT = 10, FILE = 'qraman-rl', STATUS = 'old', POSITION = 'append')
          WRITE(10, '(a, i5, 3f10.6)') '#q  ', iq, q(:)
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p1(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m1(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p1(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m1(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p1(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m1(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p1(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m1(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p1(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m1(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') mu, 0, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          DO mu = 1, nmodes
            frq = wf(branch(mu), iq) * rydcm1
            !
            Raman_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_p2(:,:,branch(mu)), ei))))**2
            Raman_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_m2(:,:,branch(mu)), ei))))**2
            Raman_ee_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_p2(:,:,branch(mu)), ei))))**2
            Raman_ee_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_ee_m2(:,:,branch(mu)), ei))))**2
            Raman_hh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_p2(:,:,branch(mu)), ei))))**2
            Raman_hh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_hh_m2(:,:,branch(mu)), ei))))**2
            Raman_eh_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_p2(:,:,branch(mu)), ei))))**2
            Raman_eh_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_eh_m2(:,:,branch(mu)), ei))))**2
            Raman_he_p = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_p2(:,:,branch(mu)), ei))))**2
            Raman_he_m = (ABS(DOT_PRODUCT(es, MATMUL(tensor_he_m2(:,:,branch(mu)), ei))))**2
            !
            ! consider the effect of Bose-Einstein occupation of phonons on the Raman intensity
            IF (lBE) Raman_p = bewph_s(branch(mu))  * Raman_p
            IF (lBE) Raman_m = bewph_as(branch(mu)) * Raman_m
            IF (lBE) Raman_ee_p = bewph_s(branch(mu))  * Raman_ee_p
            IF (lBE) Raman_ee_m = bewph_as(branch(mu)) * Raman_ee_m
            IF (lBE) Raman_hh_p = bewph_s(branch(mu))  * Raman_hh_p
            IF (lBE) Raman_hh_m = bewph_as(branch(mu)) * Raman_hh_m
            IF (lBE) Raman_eh_p = bewph_s(branch(mu))  * Raman_eh_p
            IF (lBE) Raman_eh_m = bewph_as(branch(mu)) * Raman_eh_m
            IF (lBE) Raman_he_p = bewph_s(branch(mu))  * Raman_he_p
            IF (lBE) Raman_he_m = bewph_as(branch(mu)) * Raman_he_m
            !
            WRITE(10, '(2(2x, i3), 2x, f10.4, 5(1x,2(2x, es22.12e3)))') 0, mu, frq, Raman_p/nktotf, Raman_m/nktotf, &
                                        Raman_ee_p/nktotf, Raman_ee_m/nktotf, Raman_hh_p/nktotf, Raman_hh_m/nktotf, &
                                        Raman_eh_p/nktotf, Raman_eh_m/nktotf, Raman_he_p/nktotf, Raman_he_m/nktotf
          ENDDO ! mu
          !
          CLOSE(10)
          !
        ENDIF ! polar
        !
        DEALLOCATE(branch, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating branch', 1)
        !
      ENDIF ! master node
      !
      DEALLOCATE(tensor_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p1', 1)
      DEALLOCATE(tensor_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_p2', 1)
      DEALLOCATE(tensor_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m1', 1)
      DEALLOCATE(tensor_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_m2', 1)
      DEALLOCATE(tensor_ee_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_p1', 1)
      DEALLOCATE(tensor_ee_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_p2', 1)
      DEALLOCATE(tensor_ee_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_m1', 1)
      DEALLOCATE(tensor_ee_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_ee_m2', 1)
      DEALLOCATE(tensor_hh_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_p1', 1)
      DEALLOCATE(tensor_hh_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_p2', 1)
      DEALLOCATE(tensor_hh_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_m1', 1)
      DEALLOCATE(tensor_hh_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_hh_m2', 1)
      DEALLOCATE(tensor_eh_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_p1', 1)
      DEALLOCATE(tensor_eh_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_p2', 1)
      DEALLOCATE(tensor_eh_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_m1', 1)
      DEALLOCATE(tensor_eh_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_eh_m2', 1)
      DEALLOCATE(tensor_he_p1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_p1', 1)
      DEALLOCATE(tensor_he_p2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_p2', 1)
      DEALLOCATE(tensor_he_m1, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_m1', 1)
      DEALLOCATE(tensor_he_m2, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_Raman', 'Error deallocating tensor_he_m2', 1)
      !
      IF (ALLOCATED(bewph_as)) DEALLOCATE(bewph_as)
      IF (ALLOCATED(bewph_s)) DEALLOCATE(bewph_s)
      !
      !
    ENDIF ! ho-ph
    !
    ENDSELECT ! select
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_Raman
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_gkk(iq)
    !-----------------------------------------------------------------------
    !!
    !! Print the |g| vertex for all n,n' and modes in meV and do average
    !! on degenerate states.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub
    USE elph2,         ONLY : etf, ibndmin, nkqf, xqf, nbndfst,    &
                              nkf, epf17, xkf, nkqtotf, wf, nktotf
    USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE division,      ONLY : fkbounds
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    !
    ! Local variables
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! K-point index
    INTEGER :: ikq
    !! K+q-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: pbnd
    !! Band index
    INTEGER :: nu
    !! Mode index
    INTEGER :: mu
    !! Mode index
    INTEGER :: n
    !! Number of modes
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xkf_all(3, nkqtotf)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP) :: etf_all(nbndsub, nkqtotf)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP) :: wq
    !! Phonon frequency
    REAL(KIND = DP) :: w_1
    !! Temporary phonon freq. 1
    REAL(KIND = DP) :: w_2
    !! Temporary phonon freq. 2
    REAL(KIND = DP) :: gamma
    !! Temporary electron-phonon matrix element
    REAL(KIND = DP) :: ekk
    !! Eigenenergies at k
    REAL(KIND = DP) :: ekq
    !! Eigenenergies at k+q
    REAL(KIND = DP) :: g2
    !! Temporary electron-phonon matrix element square
    REAL(KIND = DP), ALLOCATABLE :: epc(:, :, :, :)
    !! g vectex accross all pools
    REAL(KIND = DP), ALLOCATABLE :: epc_sym(:, :, :)
    !! Temporary g-vertex for each pool
    !
    ! find the bounds of k-dependent arrays in the parallel case in each pool
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    ALLOCATE(epc(nbndfst, nbndfst, nmodes, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('print_gkk', 'Error allocating epc', 1)
    ALLOCATE(epc_sym(nbndfst, nbndfst, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('print_gkk', 'Error allocating epc_sym', 1)
    !
    epc(:, :, :, :)  = zero
    epc_sym(:, :, :) = zero
    !
    ! First do the average over bands and modes for each pool
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      DO nu = 1, nmodes
        wq = wf(nu, iq)
        DO ibnd = 1, nbndfst
          DO jbnd = 1, nbndfst
            gamma = (ABS(epf17(jbnd, ibnd, nu, ik)))**two
            IF (wq > 0.d0) THEN
              gamma = gamma / (two * wq)
            ELSE
              gamma = 0.d0
            ENDIF
            gamma = DSQRT(gamma)
            ! gamma = |g| [Ry]
            epc(ibnd, jbnd, nu, ik + lower_bnd - 1) = gamma
          ENDDO ! jbnd
        ENDDO   ! ibnd
      ENDDO ! loop on modes
      !
      ! Here we "SYMMETRIZE": actually we simply take the averages over
      ! degenerate states, it is only a convention because g is gauge-dependent!
      !
      ! first the phonons
      DO ibnd = 1, nbndfst
        DO jbnd = 1, nbndfst
          DO nu = 1, nmodes
            w_1 = wf(nu, iq)
            g2 = 0.d0
            n  = 0
            DO mu = 1, nmodes
              w_2 = wf(mu, iq)
              IF (ABS(w_2 - w_1) < 0.01/ryd2mev) THEN
                n = n + 1
                g2 = g2 + epc(ibnd, jbnd, mu, ik + lower_bnd - 1) * epc(ibnd, jbnd, mu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      ! Then the k electrons
      DO nu = 1, nmodes
        DO jbnd = 1, nbndfst
          DO ibnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + ibnd, ikk)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikk)
              IF (ABS(w_2 - w_1) < 0.01/ryd2mev) THEN
                n = n + 1
                g2 = g2 + epc(pbnd, jbnd, nu, ik + lower_bnd - 1) * epc(pbnd, jbnd, nu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      !
      ! and finally the k+q electrons
      DO nu = 1, nmodes
        DO ibnd = 1, nbndfst
          DO jbnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + jbnd, ikq)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikq)
              IF (ABS(w_2 - w_1) < 0.01/ryd2mev) then
                n = n + 1
                g2 = g2 + epc(ibnd, pbnd, nu, ik + lower_bnd - 1) * epc(ibnd, pbnd, nu, ik + lower_bnd - 1)
              ENDIF
            ENDDO
            g2 = g2 / FLOAT(n)
            epc_sym(ibnd, jbnd, nu) = DSQRT(g2)
          ENDDO
        ENDDO
      ENDDO
      epc(:, :, :, ik + lower_bnd - 1) = epc_sym
      !
    ENDDO ! k-points
    !
    ! We need quantity from all the pools
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    !
#if defined(__MPI)
    !
    ! Note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL mp_sum(epc, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif
    !
    ! Only master writes
    IF (mpime == ionode_id) THEN
      !
      WRITE(stdout, '(5x, a)') ' Electron-phonon vertex |g| (meV)'
      !
      WRITE(stdout, '(/5x, "iq = ", i7, " coord.: ", 3f12.7)') iq, xqf(:, iq)
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        WRITE(stdout, '(5x, "ik = ", i7, " coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
        WRITE(stdout, '(5x, a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
        WRITE(stdout, '(5x, a)') REPEAT('-', 78)
        !
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk)
          DO jbnd = 1, nbndfst
            ekq = etf_all(ibndmin - 1 + jbnd, ikq)
            DO nu = 1, nmodes
              WRITE(stdout, '(3i9, 2f12.4, 1f20.10, 1e20.10)') ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, &
                   nu, ryd2ev * ekk, ryd2ev * ekq, ryd2mev * wf(nu, iq), ryd2mev * epc(ibnd, jbnd, nu, ik)
            ENDDO
          ENDDO
          !
        ENDDO
        WRITE(stdout, '(5x, a/)') REPEAT('-', 78)
        !
      ENDDO
    ENDIF ! master node
    !
    DEALLOCATE(epc, STAT = ierr)
    IF (ierr /= 0) CALL errore('print_gkk', 'Error deallocating epc', 1)
    DEALLOCATE(epc_sym, STAT = ierr)
    IF (ierr /= 0) CALL errore('print_gkk', 'Error deallocating epc_sym', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_gkk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob_sym(f_out, bztoibz_mat, vkk_all, etf_all, wkf_all, &
                             ef0, sigma, max_mob, xkf_all)
    !-----------------------------------------------------------------------
    !!
    !! This routine prints the mobility using k-point symmetry.
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : ncarrier, nstemp, assume_metal
    USE elph2,         ONLY : nbndfst, gtemp, nktotf, nkqtotf
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    USE symm_base,     ONLY : nsym
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(KIND = DP), INTENT(in) :: f_out(3, nbndfst, nktotf, nstemp)
    !! occupations factor produced by SERTA or IBTE
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level
    REAL(KIND = DP), INTENT(inout) :: sigma(3, 3, nstemp)
    !! Conductivity
    REAL(KIND = DP), INTENT(inout), OPTIONAL :: max_mob(nstemp)
    !! Maximum mobility
    REAL(KIND = DP), INTENT(in), OPTIONAL :: xkf_all(3, nktotf)
    !! K-point coordinate. If passed it will compute state resolved mobility
    !
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma_tmp(3, 3)
    !! Electrical conductivity
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: fi_check(3)
    !! Sum rule on population
    !
    IF (PRESENT(max_mob)) THEN
      max_mob(:) = zero
    ENDIF
    IF (PRESENT(xkf_all)) THEN
      WRITE(stdout, '(a)') ' '
      WRITE(stdout, '(5x, a)') 'Printing the state-resolved mobility to mobility_nk.fmt file'
    ELSE
      CALL prtheader_mob()
    ENDIF
    ! compute conductivity
    DO itemp = 1, nstemp
      carrier_density = 0.0
      etemp = gtemp(itemp)
      sigma_tmp(:, :) = zero
      fi_check(:) = zero
      DO ik = 1,  nktotf
        DO ibnd = 1, nbndfst
          !  energy at k (relative to Ef)
          ekk = etf_all(ibnd, ik) - ef0(itemp)
          fnk = wgauss(-ekk / etemp, -99)
          IF (etf_all(ibnd, ik) < ef0(itemp) .AND. ncarrier < -1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma_sym(f_out(:, ibnd, ik, itemp), bztoibz_mat(:, ik), &
                                   vkk_all(:, ibnd, ik), sigma_tmp, fi_check)
            ! The wkf(ikk) already include a factor 2
            carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk)
          ELSE IF (etf_all(ibnd, ik) > ef0(itemp) .AND. ncarrier > 1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma_sym(f_out(:, ibnd, ik, itemp), bztoibz_mat(:, ik), &
                                   vkk_all(:, ibnd, ik), sigma_tmp, fi_check)
            carrier_density = carrier_density + wkf_all(ik) * fnk
          ELSE IF (assume_metal) THEN
            ! just sum on all bands for metals
            CALL compute_sigma_sym(f_out(:, ibnd, ik, itemp), bztoibz_mat(:, ik), &
                                   vkk_all(:, ibnd, ik), sigma_tmp, fi_check)
            carrier_density = carrier_density + wkf_all(ik) * fnk
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      ! Print the state-resolved mobility to file
      IF (PRESENT(xkf_all)) THEN
        CALL prtmob_nk(itemp, f_out(:, :, :, itemp), vkk_all, bztoibz_mat, carrier_density, ef0(itemp), xkf_all, etemp, etf_all)
      ELSE
        ! Print the resulting mobility
        IF (PRESENT(max_mob)) THEN
          CALL prtmob(itemp, sigma_tmp, carrier_density, fi_check, ef0(itemp), etemp, max_mob(itemp))
        ELSE
          CALL prtmob(itemp, sigma_tmp, carrier_density, fi_check, ef0(itemp), etemp)
        ENDIF
      ENDIF
      sigma(:, :, itemp) = sigma_tmp(:, :)
    ENDDO ! temp
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_mob_sym
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_sigma_sym(f_out, bztoibz_mat, vkk_all, sigma, fi_check)
    !-----------------------------------------------------------------------
    !!
    !!  Computes one element of the sigma tensor when using symmetries.
    !!
    !----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE cell_base,        ONLY : at, bg
    USE epwcom,           ONLY : nkf1, nkf2, nkf3
    USE elph2,            ONLY : s_bztoibz
    USE symm_base,        ONLY : s, nsym
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: bztoibz_mat(nsym)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(KIND = DP), INTENT(in) :: f_out(3)
    !! Occupation function produced by SERTA or IBTE
    REAL(KIND = DP), INTENT(in) :: vkk_all(3)
    !! Velocity
    REAL(KIND = DP), INTENT(inout) :: sigma(3, 3)
    !! Electrical conductivity
    REAL(KIND = DP), INTENT(inout) :: fi_check(3)
    !! Sum rule on population
    !
    ! Local variables
    !
    INTEGER :: nb
    !! Number of points equivalent by sym from BZ to IBTZ
    INTEGER :: ikbz
    !! Index of full BZ points
    INTEGER :: i,j
    !! Cartesian index
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: fi_cart(3)
    !! Cartesian Fi_all
    REAL(KIND = DP) :: fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: sfac
    !! Spin factor
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    vk_cart = vkk_all
    fi_cart = f_out
    !
    ! Loop on the point equivalent by symmetry in the full BZ
    DO nb = 1, nsym
      IF (bztoibz_mat(nb) > 0) THEN
        ikbz = bztoibz_mat(nb)
        !
        ! Transform the symmetry matrix from Crystal to cartesian
        sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
        sb       = MATMUL(bg, sa)
        sr(:, :) = MATMUL(at, TRANSPOSE(sb))
        sr       = TRANSPOSE(sr)
        !
        CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vk_cart, 1 ,0.d0 , v_rot, 1)
        CALL DGEMV('n', 3, 3, 1.d0, sr, 3, fi_cart, 1 ,0.d0 , fi_rot, 1)
        !
        DO j = 1, 3
          DO i = 1, 3
            ! The factor two in the weight at the end is to
            ! account for spin
            sigma(i, j) = sigma(i, j) - (v_rot(j) * fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)
          ENDDO
        ENDDO
        !
        fi_check(:) = fi_check(:) + fi_rot(:) * sfac / (nkf1 * nkf2 * nkf3)
      ENDIF ! BZ
    ENDDO ! ikb
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_sigma_sym
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob(f_out, vkk_all, etf_all, wkf_all, ef0, sigma, max_mob)
    !-----------------------------------------------------------------------
    !!
    !! This routine prints the mobility without k-point symmetry
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : ncarrier, nstemp, assume_metal
    USE elph2,         ONLY : nbndfst, gtemp, nktotf
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: f_out(3, nbndfst, nktotf, nstemp)
    !! Occupation function produced by SERTA or IBTE
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level
    REAL(KIND = DP), INTENT(inout) :: sigma(3, 3, nstemp)
    !! Conductivity
    REAL(KIND = DP), INTENT(inout), OPTIONAL :: max_mob(nstemp)
    !! The maximum mobility computed by thr prtmob routine.
    !
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    !! Cartesian index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma_tmp(3, 3)
    !! Electrical conductivity
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: Fi_check(3)
    !! Sum rule on population
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !
    IF (PRESENT(max_mob)) THEN
      max_mob(:) = zero
    ENDIF
    CALL prtheader_mob()
    DO itemp = 1, nstemp
      carrier_density = 0.0
      etemp = gtemp(itemp)
      sigma_tmp(:, :) = zero
      fi_check(:) = zero
      DO ik = 1,  nktotf
        DO ibnd = 1, nbndfst
          !  energy at k (relative to Ef)
          ekk = etf_all(ibnd, ik) - ef0(itemp)
          fnk = wgauss(-ekk / etemp, -99)
          IF (etf_all(ibnd, ik) < ef0(itemp) .AND. ncarrier < -1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma(f_out(:, ibnd, ik, itemp), vkk_all(:, ibnd, ik), wkf_all(ik), sigma_tmp, fi_check)
            ! The wkf(ikk) already include a factor 2
            carrier_density = carrier_density + wkf_all(ik) * (1.0d0 - fnk)
          ELSE IF (etf_all(ibnd, ik) > ef0(itemp) .AND. ncarrier > 1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma(f_out(:, ibnd, ik, itemp), vkk_all(:, ibnd, ik), wkf_all(ik), sigma_tmp, fi_check)
            ! The wkf(ikk) already include a factor 2
            carrier_density = carrier_density + wkf_all(ik) * fnk
          ELSE IF (assume_metal) THEN
            ! sum on all bands for metals
            CALL compute_sigma(f_out(:, ibnd, ik, itemp), vkk_all(:, ibnd, ik), wkf_all(ik), sigma_tmp, fi_check)
            carrier_density = carrier_density + wkf_all(ik) * fnk
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      IF (PRESENT(max_mob)) THEN
        CALL prtmob(itemp, sigma_tmp, carrier_density, fi_check, ef0(itemp), etemp, max_mob(itemp))
      ELSE
        CALL prtmob(itemp, sigma_tmp, carrier_density, fi_check, ef0(itemp), etemp)
      ENDIF
      sigma(:, :, itemp) = sigma_tmp(:, :)
    ENDDO ! itemp
    !-----------------------------------------------------------------------
    END SUBROUTINE print_mob
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_sigma(f_out, vkk_all, wkf_all, sigma, fi_check)
    !-----------------------------------------------------------------------
    !!
    !! Computes the conductivity tensor without using symetries
    !! With or without B-field
    !!
    !----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE epwcom,           ONLY : nkf1, nkf2, nkf3
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: f_out(3)
    !! Occupation function produced by SERTA or IBTE for a given k-point and band index
    REAL(KIND = DP), INTENT(INOUT) :: sigma(3, 3)
    !! Electrical conductivity
    REAL(KIND = DP), INTENT(INOUT) :: fi_check(3)
    !! Sum rule on population
    REAL(KIND = DP), INTENT(in) :: vkk_all(3)
    !! Velocity for a given k-point and band
    REAL(KIND = DP), INTENT(in) :: wkf_all
    !! Weight of k
    !
    ! Local variables
    !
    INTEGER :: i, j
    !! Dimension loop indices
    REAL(KIND = DP) :: sfac
    !! Spin factor
    IF (noncolin) THEN
      sfac = 1.0
    ELSE
      sfac = 2.0
    ENDIF
    !
    DO j = 1, 3
      DO i = 1, 3
        sigma(i, j) = sigma(i, j) - vkk_all(i) * f_out(j) * wkf_all
      ENDDO
    ENDDO
    fi_check(:) = fi_check(:) + f_out(:) * sfac / (nkf1 * nkf2 * nkf3)
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_sigma
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE prtmob(itemp, sigma, carrier_density, fi_check, ef0, etemp, max_mob)
    !-----------------------------------------------------------------------
    !!
    !! This routine print the mobility (or conducrtivity for metals) in a
    !! nice format and in proper units.
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : assume_metal, system_2d
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : omega, at, alat
    USE elph2,         ONLY : dos
    USE constants_epw, ONLY : zero, kelvin2eV, ryd2ev, eps80, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: sigma(3, 3)
    !! Conductivity tensor
    REAL(KIND = DP), INTENT(in) :: carrier_density
    !! Carrier density in a.u.
    REAL(KIND = DP), INTENT(in) :: fi_check(3)
    !! Integrated population vector
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi-level
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(inout), OPTIONAL :: max_mob
    !! Maximum error for all temperature
    !
    ! Local variables
    REAL(KIND = DP) :: mobility(3, 3)
    !! mobility
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: nden
    !! Carrier density in cm^-3
    !
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
    ELSE
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
    ENDIF
    ! carrier_density in cm^-1
    IF (system_2d == 'no') THEN
      nden = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3.0d0)
    ELSE
      nden = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-2.0d0)
    ENDIF
    mobility(:, :) = (sigma(:, :) * electron_si ** 2 * inv_cell) / (hbarJ * bohr2ang * ang2cm)
    IF (.NOT. assume_metal) THEN
      ! for insulators print mobility so just divide by carrier density
      IF (ABS(nden) < eps80) CALL errore('prtmob', 'The carrier density is 0', 1)
      mobility(:, :) = (mobility(:, :) / (electron_si * carrier_density * inv_cell)) * (bohr2ang * ang2cm) ** 3
      WRITE(stdout, '(5x, 1f8.3, 1f9.4, 1E14.5, 1E14.5, 3E16.6)') etemp * ryd2ev / kelvin2eV, ef0 * ryd2ev, &
           nden, fi_check(1), mobility(1, 1), mobility(1, 2), mobility(1, 3)
    ELSE
      WRITE(stdout, '(5x, 1f8.3, 1f9.4, 1E14.5, 1E14.5, 3E16.6)') etemp * ryd2ev / kelvin2eV, ef0 * ryd2ev, &
           dos(itemp), fi_check(1), mobility(1, 1), mobility(1, 2), mobility(1, 3)
    ENDIF
    WRITE(stdout, '(36x, 1E14.5, 3E16.6)') fi_check(2), mobility(2, 1), mobility(2, 2), mobility(2, 3)
    WRITE(stdout, '(36x, 1E14.5, 3E16.6)') fi_check(3), mobility(3, 1), mobility(3, 2), mobility(3, 3)
    IF (PRESENT(max_mob)) THEN
      max_mob = MAXVAL(mobility(:,:))
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE prtmob
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE prtmob_nk(itemp, f_out, vkk_all, bztoibz_mat, carrier_density, ef0, xkf_all, etemp, etf_all)
    !-----------------------------------------------------------------------
    !!
    !! This routine print the state resolved (nk) mobility in a file called
    !! nice format and in proper units.
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : assume_metal, system_2d, ncarrier, nkf1, nkf2, nkf3, nstemp
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : omega, at, alat, bg
    USE elph2,         ONLY : dos, nkqtotf, s_bztoibz, nbndfst, nktotf
    USE io_var,        ONLY : iufilmu_nk
    USE symm_base,     ONLY : s, nsym
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE constants_epw, ONLY : zero, kelvin2eV, ryd2ev, eps80, eps10, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    INTEGER, INTENT(in) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point from the IBZ, given the index of all k from full BZ
    REAL(KIND = DP), INTENT(in) :: f_out(3, nbndfst, nktotf)
    !! occupations factor produced by SERTA or IBTE
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: carrier_density
    !! Carrier density in a.u.
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi-level
    REAL(KIND = DP), INTENT(in) :: xkf_all(3, nktotf)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies of k
    !
    ! Local variables
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: nb
    !! Number of points equivalent by sym from BZ to IBTZ
    INTEGER :: ikbz
    !! Index of full BZ points
    INTEGER :: i,j
    !! Cartesian index
    REAL(KIND = DP) :: xkk_cart(3)
    !! k-point in cartesian coordinate
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: fi_cart(3)
    !! Cartesian Fi_all
    REAL(KIND = DP) :: xkk_rot(3)
    !! Rotated k-point by the symmetry operation
    REAL(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: fi_rot(3)
    !! Rotated Fi_all by the symmetry operation
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: sfac
    !! Spin factor
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: sigma(3, 3)
    !! Conductivity
    REAL(KIND = DP) :: mobility(3, 3)
    !! Mobility in cm^2/Vs
    REAL(KIND = DP) :: mobility_nk(3, 3)
    !! State-resovled (nk) mobility in cm^s/Vs
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3} or Bohr^{-2}]
    !
    ! Write to file
    IF (mpime == ionode_id) THEN
      IF (itemp == 1) THEN
        OPEN(UNIT = iufilmu_nk, FILE = 'mobility_nk.fmt')
        WRITE(iufilmu_nk, '(a)') '# State-resolved carrier mobility'
      ENDIF
      WRITE(iufilmu_nk, '(a, 1f12.6, a)') '# Temperature', etemp * ryd2ev / kelvin2eV, ' K'
      WRITE(iufilmu_nk, '(a)') &
      '# k-point band isym          k-point (Cartesian)      Energy - ef (eV)             mu_nk(alpha, beta) (cm^2/Vs)'
      !
      ! Spin factor
      IF (noncolin) THEN
        sfac = 1.0
      ELSE
        sfac = 2.0
      ENDIF
      !
      IF (system_2d == 'no') THEN
        inv_cell = 1.0d0 / omega
      ELSE
        ! for 2d system need to divide by area (vacuum in z-direction)
        inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
      ENDIF
      !
      mobility(:, :) = zero
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibnd, ik) - ef0
          vk_cart = vkk_all(:, ibnd, ik)
          fi_cart = f_out(:, ibnd, ik)
          xkk_cart = xkf_all(:, ik)
          CALL cryst_to_cart(1, xkk_cart, bg, 1)
          DO nb = 1, nsym
            IF (bztoibz_mat(nb, ik) > 0) THEN
              sigma(:, :) = zero
              mobility_nk(:, :) = zero
              ikbz = bztoibz_mat(nb, ik)
              ! Transform the symmetry matrix from Crystal to cartesian
              sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
              sb       = MATMUL(bg, sa)
              sr(:, :) = MATMUL(at, TRANSPOSE(sb))
              sr       = TRANSPOSE(sr)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xkk_cart, 1 ,0.d0 , xkk_rot, 1)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vk_cart, 1 ,0.d0 , v_rot, 1)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, fi_cart, 1 ,0.d0 , fi_rot, 1)
              DO j = 1, 3
                DO i = 1, 3
                  ! The factor two in the weight at the end is to
                  ! account for spin
                  sigma(i, j) = - (v_rot(j) * fi_rot(i)) * sfac / (nkf1 * nkf2 * nkf3)
                ENDDO
              ENDDO
              ! nk-resolved mobility $\mu_{nk}^{\alpha\beta}$
              mobility_nk(:, :) = sigma(:, :) * electron_si * (bohr2ang * ang2cm) ** 2 / (hbarJ * carrier_density)
              !
              IF (ekk < eps10 .AND. ncarrier < -1E5 .AND. (mobility_nk(1,1) + mobility_nk(2,2) + mobility_nk(3,3)) > eps80) THEN
                WRITE(iufilmu_nk, '(i9, 2i5, 4f12.8, 3E18.8)') ik, ibnd, nb, xkk_rot(:), ekk * ryd2ev, &
                                                                       mobility_nk(1, 1), mobility_nk(1, 2), mobility_nk(1, 3)
                WRITE(iufilmu_nk, '(67x, 3E18.8)') mobility_nk(2, 1), mobility_nk(2, 2), mobility_nk(2, 3)
                WRITE(iufilmu_nk, '(67x, 3E18.8)') mobility_nk(3, 1), mobility_nk(3, 2), mobility_nk(3, 3)
                mobility(:, :) =  mobility(:, :) + mobility_nk(:, :)
              ELSE IF (ekk > eps10 .AND. ncarrier > 1E5 .AND. (mobility_nk(1,1) + mobility_nk(2,2) + mobility_nk(3,3)) > eps80) THEN
                WRITE(iufilmu_nk, '(i9, 2i5, 4f12.8, 3E18.8)') ik, ibnd, nb, xkk_rot(:), ekk * ryd2ev, &
                                                                       mobility_nk(1, 1), mobility_nk(1, 2), mobility_nk(1, 3)
                WRITE(iufilmu_nk, '(67x, 3E18.8)') mobility_nk(2, 1), mobility_nk(2, 2), mobility_nk(2, 3)
                WRITE(iufilmu_nk, '(67x, 3E18.8)') mobility_nk(3, 1), mobility_nk(3, 2), mobility_nk(3, 3)
                mobility(:, :) =  mobility(:, :) + mobility_nk(:, :)
              ENDIF
            ENDIF ! bztoibz_mat
          ENDDO ! nb
        ENDDO ! ibnd
      ENDDO ! ik
      !
      ! Integrated mobility
      WRITE(iufilmu_nk, '(a)') ' '
      WRITE(iufilmu_nk, '(a)') '# Integrated mobility (cm^/Vs) - Should be the same as in the output. '
      WRITE(iufilmu_nk, '(a, 3E16.6)') 'mu(alpha, beta) = ', mobility(1, 1), mobility(1, 2), mobility(1, 3)
      WRITE(iufilmu_nk, '(a, 3E16.6)') '                  ', mobility(2, 1), mobility(2, 2), mobility(2, 3)
      WRITE(iufilmu_nk, '(a, 3E16.6)') '                  ', mobility(3, 1), mobility(3, 2), mobility(3, 3)
      WRITE(iufilmu_nk, '(a)') ' '
      IF (itemp == nstemp) CLOSE(iufilmu_nk)
    ENDIF ! master core
    !-----------------------------------------------------------------------
    END SUBROUTINE prtmob_nk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE prtheader_mob()
    !-----------------------------------------------------------------------
    !!
    !! This routine print a header for mobility calculation
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : assume_metal, ncarrier, system_2d
    !
    IMPLICIT NONE
    !
    WRITE(stdout, '(/5x, a)') REPEAT('=', 93)
    IF (.NOT. assume_metal) THEN
      IF (ncarrier < -1E5) THEN
        WRITE(stdout, '(5x, "  Temp     Fermi   Hole density  Population SR            Drift Hole mobility ")')
        IF (system_2d == 'no') THEN
          WRITE(stdout, '(5x, "   [K]      [eV]     [cm^-3]      [h per cell]                    [cm^2/Vs]")')
        ELSE
          WRITE(stdout, '(5x, "   [K]      [eV]     [cm^-2]      [h per cell]                    [cm^2/Vs]")')
        ENDIF
      ELSE
        WRITE(stdout, '(5x, "  Temp     Fermi   Elec density  Population SR            Drift Elec mobility ")')
        IF (system_2d == 'no') THEN
          WRITE(stdout, '(5x, "   [K]      [eV]     [cm^-3]      [e per cell]                    [cm^2/Vs]")')
        ELSE
          WRITE(stdout, '(5x, "   [K]      [eV]     [cm^-2]      [e per cell]                    [cm^2/Vs]")')
        ENDIF
      ENDIF
    ELSE
      WRITE(stdout, '(5x, "  Temp     Fermi        DOS        Population SR                 Conductivity ")')
      WRITE(stdout, '(5x, "   [K]      [eV]    [states/Ry] [carriers per cell]               [Ohm.cm]^-1 ")')
    ENDIF
    WRITE(stdout, '(5x, a/)') REPEAT('=', 93)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE prtheader_mob
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE prtheader_supercond(itemp, cal_type)
    !-----------------------------------------------------------------------
    !!
    !! This routine print a header for superconductivity calculation
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : liso, laniso, lreal, imag_read, wscut, &
                              fbw, broyden_beta
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, nsw
    USE constants_epw, ONLY : kelvin2eV, zero
    USE constants,     ONLY : pi
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    INTEGER, INTENT(in) :: cal_type
    !! 1 = limag, 2 = lpade, 3 = lacon, 4 = lreal
    !
    IF (cal_type == 1) THEN
      WRITE(stdout, '(a)') '    '
      WRITE(stdout, '(5x, a, i3, a, f12.5, a, a, i3, a)') 'temp(', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K'
      WRITE(stdout, '(a)') '    '
      IF (liso .AND. .NOT. fbw) &
        WRITE(stdout, '(5x, a)') 'Solve isotropic Eliashberg equations on imaginary-axis'
      IF (liso .AND. fbw) &
        WRITE(stdout, '(5x, a)') 'Solve full-bandwidth isotropic Eliashberg equations on imaginary-axis'
      IF (laniso .AND. .NOT. fbw .AND. .NOT. imag_read) &
        WRITE(stdout, '(5x, a)') 'Solve anisotropic Eliashberg equations on imaginary-axis'
      IF (laniso .AND. fbw .AND. .NOT. imag_read) &
        WRITE(stdout, '(5x, a)') 'Solve full-bandwidth anisotropic Eliashberg equations on imaginary-axis'
      IF (laniso .AND. .NOT. fbw .AND. imag_read .AND. itemp == 1) &
        WRITE(stdout, '(5x, a)') 'Read from file delta and znorm on imaginary-axis'
      IF (laniso .AND. fbw .AND. imag_read .AND. itemp == 1) &
        WRITE(stdout, '(5x, a)') 'Read from file delta and znorm and shift on imaginary-axis'
      WRITE(stdout, '(a)') '    '
      WRITE(stdout, '(5x, a, i6, a, i6)') 'Total number of frequency points nsiw(', itemp, ') = ', nsiw(itemp)
      WRITE(stdout, '(5x, a, f10.4)') 'Cutoff frequency wscut = ', (2.d0 * nsiw(itemp) + 1) * pi * gtemp(itemp)
      IF (broyden_beta < zero ) THEN
        WRITE(stdout, '(5x, a, f12.5)') 'linear mixing factor = ', DABS(broyden_beta)
      ELSE
        WRITE(stdout, '(5x, a, f12.5)') 'broyden mixing factor = ', broyden_beta
      ENDIF
      IF (broyden_beta  < -0.2d0) THEN
        WRITE(stdout, '(5x, a)') 'mixing factor = 0.2 is used for the first three iterations.'
      ENDIF
      WRITE(stdout, '(a)') '    '
    ENDIF
    !
    IF (cal_type == 2) THEN
      WRITE(stdout, '(a)') '    '
      IF (liso .AND. .NOT. fbw) WRITE(stdout, '(5x, a)') &
        'Pade approximant of isotropic Eliashberg equations from imaginary-axis to real-axis'
      IF (laniso .AND. .NOT. fbw) WRITE(stdout, '(5x, a)') &
        'Pade approximant of anisotropic Eliashberg equations from imaginary-axis to real-axis'
      IF (liso .AND. fbw) WRITE(stdout, '(5x, a)') &
        'Pade approximant of full-bandwidth isotropic Eliashberg equations from imaginary-axis to real-axis'
      IF (laniso .AND. fbw) WRITE(stdout, '(5x, a)') &
        'Pade approximant of full-bandwidth anisotropic Eliashberg equations from imaginary-axis to real-axis'
      WRITE(stdout, '(5x, a, f10.4)') 'Cutoff frequency wscut = ', wscut
      WRITE(stdout, '(a)') '    '
    ENDIF
    !
    IF (cal_type == 3) THEN
      WRITE(stdout, '(a)') '    '
      IF (liso) &
        WRITE(stdout, '(5x, a)') 'Analytic continuation of isotropic Eliashberg equations from imaginary-axis to real-axis'
      IF (laniso) &
        WRITE(stdout, '(5x, a)') 'Analytic continuation of anisotropic Eliashberg equations from imaginary-axis to real-axis'
      WRITE(stdout, '(a)') '    '
      WRITE(stdout, '(5x, a, i6)') 'Total number of frequency points nsw = ', nsw
      WRITE(stdout, '(5x, a, f10.4)') 'Cutoff frequency wscut = ', wscut
      WRITE(stdout, '(a)') '    '
    ENDIF
    !
    IF (cal_type == 4) THEN
      WRITE(stdout, '(a)') '    '
      WRITE(stdout, '(5x, a, i3, a, f12.5, a, a, i3, a)') 'temp(', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K'
      WRITE(stdout, '(a)') '    '
      IF (liso .AND. lreal) &
        WRITE(stdout, '(5x, a)') 'Solve isotropic Eliashberg equations on real-axis'
      WRITE(stdout, '(a)') '    '
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE prtheader_supercond
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_clock_epw
    !-----------------------------------------------------------------------
    !!
    !! Adapted from the code PH/print_clock_ph - Quantum-ESPRESSO group
    !!
    USE io_global,     ONLY : stdout
    USE uspp,          ONLY : nlcc_any
    !
    IMPLICIT NONE
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Unfolding on the coarse grid'
    CALL print_clock('dvanqq2')
    CALL print_clock('elphon_wrap')
    WRITE(stdout, '(5x)')
    CALL print_clock('ELPHWAN')
    WRITE(stdout, '(5x, a)') 'INITIALIZATION: '
    CALL print_clock('epq_init')
    WRITE(stdout, '(5x)')
    CALL print_clock('epq_init')
    IF (nlcc_any) CALL print_clock('set_drhoc')
    CALL print_clock('init_vloc')
    CALL print_clock('init_us_1')
    CALL print_clock('newdq2')
    CALL print_clock('dvanqq2')
    CALL print_clock('drho')
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x)')
    !
    ! Electron-Phonon interpolation
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Electron-Phonon interpolation'
    CALL print_clock('ephwann')
    CALL print_clock('ep-interp')
    CALL print_clock('ep-int-ahc')
    CALL print_clock('PH SELF-ENERGY')
    CALL print_clock('ABS SPECTRA')
    CALL print_clock('crys_cart')
    WRITE(stdout, '(5x)')
    CALL print_clock('load data')
    CALL print_clock('Ham: step 1')
    CALL print_clock('Ham: step 2')
    CALL print_clock('Ham: step 3')
    CALL print_clock('Ham: step 4')
    CALL print_clock('ep: step 1')
    CALL print_clock('ep: step 2')
    CALL print_clock('ep: step 3')
    CALL print_clock('ep: step 4')
    CALL print_clock('unfold_sthmat')
    CALL print_clock('collect_sthmat')
    CALL print_clock('dg: step 1')
    CALL print_clock('dg: step 2')
    CALL print_clock('sth: step 1')
    CALL print_clock('sth: step 2')
    CALL print_clock('DynW2B')
    CALL print_clock('HamW2B')
    CALL print_clock('ephW2Bp')
    CALL print_clock('ephW2B')
    CALL print_clock('sthW2Bp')
    CALL print_clock('dgW2B')
    CALL print_clock('dwW2B')
    CALL print_clock('print_ibte')
    CALL print_clock('vmewan2bloch')
    CALL print_clock('vmewan2blochp')
    CALL print_clock('ephW2Bp1')
    CALL print_clock('ephW2Bp2')
    CALL print_clock('kpoint_paral')
    !
    ! Eliashberg equations
    WRITE(stdout, '(5x)')
    CALL print_clock('ELIASHBERG')
    !
    WRITE(stdout, '(5x)')
    WRITE(stdout, '(5x, a)') 'Total program execution'
    CALL print_clock('EPW')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_clock_epw
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE plot_band()
    !-----------------------------------------------------------------------
    !!
    !! This routine writes output files for phonon dispersion and band structure
    !! SP : Modified so that it works with the current plotband.x of QE 5
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, filqf, filkf
    USE elph2,         ONLY : etf, nkf, nqtotf, wf, xkf, xqf, nkqtotf, nktotf
    USE constants_epw, ONLY : ryd2mev, ryd2ev, zero
    USE io_var,        ONLY : iufilfreq, iufileig
    USE elph2,         ONLY : nkqf
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! Global k-point index
    INTEGER :: ikk
    !! Index for the k-point
    INTEGER :: ikq
    !! Index for the q-point
    INTEGER :: ibnd
    !! Band index
    INTEGER :: imode
    !! Mode index
    INTEGER :: iq
    !! Global q-point index
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: dist
    !! Distance from G-point
    REAL(KIND = DP) :: dprev
    !! Previous distance
    REAL(KIND = DP) :: dcurr
    !! Current distance
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! K-points on the full k grid (all pools)
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigenenergies on the full k grid (all pools)
    !
    IF (filqf /= ' ') THEN
      !
      IF (my_pool_id == ionode_id) THEN
        !
        OPEN(iufilfreq, FILE = "phband.freq", FORM = 'formatted')
        WRITE(iufilfreq, '(" &plot nbnd=", i4, ", nks=", i6, " /")') nmodes, nqtotf
        !
        ! crystal to cartesian coordinates
        CALL cryst_to_cart(nqtotf, xqf, bg, 1)
        !
        dist  = zero
        dprev = zero
        dcurr = zero
        DO iq = 1, nqtotf
          !
          IF (iq /= 1) THEN
            dist = DSQRT((xqf(1, iq) - xqf(1, iq - 1)) * (xqf(1, iq) - xqf(1, iq - 1)) &
                       + (xqf(2, iq) - xqf(2, iq - 1)) * (xqf(2, iq) - xqf(2, iq - 1)) &
                       + (xqf(3, iq) - xqf(3, iq - 1)) * (xqf(3, iq) - xqf(3, iq - 1)))
          ELSE
            dist = zero
          ENDIF
          dcurr = dprev + dist
          dprev = dcurr
          WRITE(iufilfreq, '(10x, 3f10.6)') xqf(:, iq)
          WRITE(iufilfreq, '(1000f14.4)') (wf(imode, iq) * ryd2mev, imode = 1, nmodes)
          !
        ENDDO
        CLOSE(iufilfreq)
        !
        ! back from cartesian to crystal coordinates
        CALL cryst_to_cart(nqtotf, xqf, at, -1)
        !
      ENDIF
    ENDIF ! filqf
    !
    IF (filkf /= ' ') THEN
      !
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
      ENDDO
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error allocating etf_all', 1)
      !
#if defined(__MPI)
      CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_barrier(inter_pool_comm)
#else
      !
      xkf_all = xkf
      etf_all = etf
#endif
      !
      IF (my_pool_id == ionode_id) THEN
        !
        OPEN(iufileig, FILE = "band.eig", FORM = 'formatted')
        WRITE(iufileig, '(" &plot nbnd=", i4, ", nks=", i6, " /")') nbndsub, nktotf
        !
        ! crystal to cartesian coordinates
        CALL cryst_to_cart(nkqtotf, xkf_all, bg, 1)
        !
        dist  = zero
        dprev = zero
        dcurr = zero
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF (ikk /= 1) THEN
            dist = DSQRT((xkf_all(1, ikk) - xkf_all(1, ikk - 2)) * (xkf_all(1, ikk) - xkf_all(1, ikk - 2)) &
                       + (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) * (xkf_all(2, ikk) - xkf_all(2, ikk - 2)) &
                       + (xkf_all(3, ikk) - xkf_all(3, ikk - 2)) * (xkf_all(3, ikk) - xkf_all(3, ikk - 2)))
          ELSE
            dist = 0.d0
          ENDIF
          dcurr = dprev + dist
          dprev = dcurr
          WRITE(iufileig, '(10x, 3f10.6)') xkf_all(:, ikk)
          WRITE(iufileig, '(1000f20.12)') (etf_all(ibnd, ikk) * ryd2ev, ibnd = 1, nbndsub)
          !
        ENDDO
        CLOSE(iufileig)
        !
        ! back from cartesian to crystal coordinates
        CALL cryst_to_cart(nkqtotf, xkf_all, at, -1)
        !
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('plot_band', 'Error deallocating etf_all', 1)
      !
    ENDIF ! filkf
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE plot_band
    !----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE plot_fermisurface()
    !-----------------------------------------------------------------------
    !!
    !! This routine writes output files (in .cube) to plot Fermi surfaces
    !! HP : 9/23/2020
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE pwcom,         ONLY : ef
    USE epwcom,        ONLY : mp_mesh_k, nbndsub, nkf1, nkf2, nkf3
    USE elph2,         ONLY : bztoibz, etf, xkf, nkqtotf, nktotf
    USE constants_epw, ONLY : ryd2ev, zero
    USE io_var,        ONLY : iufilFS
    USE elph2,         ONLY : nkqf, ibndmin, ibndmax
    USE io_global,     ONLY : ionode_id, stdout
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE poolgathering, ONLY : poolgather2
    USE io_files,      ONLY : prefix
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! file name
    INTEGER :: ios
    !! IO error message
    INTEGER :: ik
    !! Global k-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: ierr
    !! Error status
    INTEGER :: i
    !! Counter
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! K-points on the full k grid (all pools)
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigenenergies on the full k grid (all pools)
    !
    WRITE(stdout,'(5x, a38)') 'Fermi surface calculation on fine mesh'
    WRITE(stdout, '(5x, a32, f10.6)') 'Fermi level (eV) = ', ef * ryd2ev
    WRITE(stdout,'(5x, i7, a32/)') ibndmax - ibndmin + 1, ' bands within the Fermi window'
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_fermisurface', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_fermisurface', 'Error allocating etf_all', 1)
    !
#if defined(__MPI)
    CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL mp_barrier(inter_pool_comm)
#else
    !
    xkf_all = xkf
    etf_all = etf
#endif
    !
    IF (my_pool_id == ionode_id) THEN
      !
      ! .cube format
      DO ibnd = ibndmin, ibndmax
        !
        IF (ibnd - ibndmin + 1 < 10) THEN
          WRITE(name1, '(a, a4, i1, a5)') TRIM(prefix), '.fs_', ibnd - ibndmin + 1, '.cube'
        ELSEIF (ibnd - ibndmin + 1 < 100) THEN
          WRITE(name1, '(a, a4, i2, a5)') TRIM(prefix), '.fs_', ibnd - ibndmin + 1, '.cube'
        ELSE
          CALL errore( 'plot_fermisurface', 'Too many bands ',1)
        ENDIF
        !
        OPEN(iufilFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('plot_fermisurface', 'error opening file ' // name1, iufilFS)
        WRITE(iufilFS, *) 'Cube file created from EPW calculation'
        WRITE(iufilFS, '(a20, f10.6)') 'Fermi level (eV) = ', ef * ryd2ev
        WRITE(iufilFS, '(i5, 3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilFS, '(i5, 3f12.6)') nkf1, (bg(i, 1) / DBLE(nkf1), i = 1, 3)
        WRITE(iufilFS, '(i5, 3f12.6)') nkf2, (bg(i, 2) / DBLE(nkf2), i = 1, 3)
        WRITE(iufilFS, '(i5, 3f12.6)') nkf3, (bg(i, 3) / DBLE(nkf3), i = 1, 3)
        WRITE(iufilFS, '(i5, 4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
        IF (mp_mesh_k) THEN
          WRITE(iufilFS, '(6f12.6)') ((etf_all(ibnd, 2 * bztoibz(ik) - 1) - ef) * ryd2ev, ik = 1, nkf1 * nkf2 * nkf3)
        ELSE
          WRITE(iufilFS, '(6f12.6)') ((etf_all(ibnd, ik * 2 - 1) - ef) * ryd2ev, ik = 1, nktotf)
        ENDIF
        CLOSE(iufilFS)
        !
      ENDDO
      !
      ! HP: Write in .frmsf format compatible with fermisurfer program
      WRITE(name1, '(a, a3, a6)') TRIM(prefix), '.fs', '.frmsf'
      OPEN(iufilFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('plot_fermisurface', 'error opening file ' // name1, iufilFS)
      WRITE(iufilFS, '(3i5)') nkf1, nkf2, nkf3
      WRITE(iufilFS, '(i5)') 1
      WRITE(iufilFS, '(i5)') ibndmax - ibndmin + 1
      WRITE(iufilFS, '(3f12.6)') (bg(i, 1), i = 1, 3)
      WRITE(iufilFS, '(3f12.6)') (bg(i, 2), i = 1, 3)
      WRITE(iufilFS, '(3f12.6)') (bg(i, 3), i = 1, 3)
      IF (mp_mesh_k) THEN
        WRITE(iufilFS, '(6f12.6)') (((etf_all(ibnd, 2 * bztoibz(ik) - 1) - ef) * ryd2ev, &
                ik = 1, nkf1 * nkf2 * nkf3), ibnd = ibndmin, ibndmax)
      ELSE
        WRITE(iufilFS, '(6f12.6)') (((etf_all(ibnd, ik * 2 - 1) - ef) * ryd2ev, &
                ik = 1, nktotf), ibnd = ibndmin, ibndmax)
      ENDIF
      CLOSE(iufilFS)
      !
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_fermisurface', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plot_fermisurface', 'Error deallocating etf_all', 1)
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE plot_fermisurface
    !----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_mob_b(f_out_b, vkk_all_b, etf_all_b, wkf_all_b, ef0, &
                           carrier_density, sigma, max_mob)
    !-----------------------------------------------------------------------
    !!
    !! This routine prints the mobility without k-point symmetry with finite B-field
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : ncarrier, nstemp, assume_metal
    USE elph2,         ONLY : nbndfst, nkpt_bztau_max, gtemp
    USE constants_epw, ONLY : zero, two, pi, kelvin2eV, ryd2ev, eps10, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: f_out_b(3, nbndfst, nkpt_bztau_max, nstemp)
    !! Occupation function produced by SERTA or IBTE
    REAL(KIND = DP), INTENT(in) :: vkk_all_b(3, nbndfst, nkpt_bztau_max)
    !! Velocity
    REAL(KIND = DP), INTENT(in) :: etf_all_b(nbndfst, nkpt_bztau_max)
    !! Eigenenergies of k
    REAL(KIND = DP), INTENT(in) :: wkf_all_b(nkpt_bztau_max)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level
    REAL(KIND = DP), INTENT(inout) :: carrier_density(nstemp)
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP), INTENT(inout) :: sigma(3, 3, nstemp)
    !! Conductivity
    REAL(KIND = DP), INTENT(inout), OPTIONAL :: max_mob(nstemp)
    !! The maximum mobility computed by thr prtmob routine.
    !
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: sigma_tmp(3, 3)
    !! Electrical conductivity
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP) :: Fi_check(3)
    !! Sum rule on population
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !
    IF (PRESENT(max_mob)) THEN
      max_mob(:) = zero
    ENDIF
    CALL prtheader_mob()
    DO itemp = 1, nstemp
      carrier_density(itemp) = 0.0
      etemp = gtemp(itemp)
      sigma_tmp(:, :) = zero
      fi_check(:) = zero
      DO ik = 1,  nkpt_bztau_max
        DO ibnd = 1, nbndfst
          !  energy at k (relative to Ef)
          ekk = etf_all_b(ibnd, ik) - ef0(itemp)
          fnk = wgauss(-ekk / etemp, -99)
          IF (etf_all_b(ibnd, ik) < ef0(itemp) .AND. ncarrier < -1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma(f_out_b(:, ibnd, ik, itemp), vkk_all_b(:, ibnd, ik), wkf_all_b(ik), sigma_tmp, fi_check)
            ! The wkf(ikk) already include a factor 2
            carrier_density(itemp) = carrier_density(itemp) + wkf_all_b(ik) * (1.0d0 - fnk)
          ELSEIF (etf_all_b(ibnd, ik) > ef0(itemp) .AND. ncarrier > 1E5 .AND. .NOT. assume_metal) THEN
            CALL compute_sigma(f_out_b(:, ibnd, ik, itemp), vkk_all_b(:, ibnd, ik), wkf_all_b(ik), sigma_tmp, fi_check)
            ! The wkf(ikk) already include a factor 2
            carrier_density(itemp) = carrier_density(itemp) + wkf_all_b(ik) * fnk
          ELSEIF (assume_metal) THEN
            ! sum on all bands for metals
            CALL compute_sigma(f_out_b(:, ibnd, ik, itemp), vkk_all_b(:, ibnd, ik), wkf_all_b(ik), sigma_tmp, fi_check)
            carrier_density(itemp) = carrier_density(itemp) + wkf_all_b(ik) * fnk
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      IF (PRESENT(max_mob)) THEN
        CALL prtmob(itemp, sigma_tmp, carrier_density(itemp), fi_check, ef0(itemp), etemp, max_mob(itemp))
      ELSE
        CALL prtmob(itemp, sigma_tmp, carrier_density(itemp), fi_check, ef0(itemp), etemp)
      ENDIF
      sigma(:, :, itemp) = sigma_tmp(:, :)
    ENDDO ! itemp
    !-----------------------------------------------------------------------
    END SUBROUTINE print_mob_b
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_hall(carrier_density, sigma_serta, sigma_bte, sigmab_serta, sigmab_bte)
    !-----------------------------------------------------------------------
    !!
    !! This routine print the final conductivity, mobility and Hall factor
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : system_2d, nstemp, bfieldx, bfieldy, bfieldz
    USE low_lvl,       ONLY : matinv3
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : omega, at, alat
    USE elph2,         ONLY : gtemp
    USE constants_epw, ONLY : zero, kelvin2eV, ryd2ev, eps80, cm2m, &
                              bohr2ang, ang2cm, hbarJ
    USE constants,     ONLY : electron_si
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: carrier_density(nstemp)
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP), INTENT(in) :: sigma_serta(3, 3, nstemp)
    !! SERTA conductivity tensor without Magnetic field.
    REAL(KIND = DP), INTENT(in) :: sigma_bte(3, 3, nstemp)
    !! BTE conductivity tensor without Magnetic field.
    REAL(KIND = DP), INTENT(in) :: sigmab_serta(3, 3, nstemp)
    !! SERTA conductivity tensor with Magnetic field.
    REAL(KIND = DP), INTENT(in) :: sigmab_bte(3, 3, nstemp)
    !! BTE conductivity tensor with Magnetic field.
    !
    ! Local variables
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: carrier_density_cm(nstemp)
    !! Carrier density (cm^(-d))
    REAL(KIND = DP) :: etemp
    !! Temperature value in Ry
    REAL(KIND = DP) :: b_norm
    !! Norm of the B-field
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume or Area
    REAL(KIND = DP) :: mob_serta(3, 3, nstemp)
    !! SERTA mobility tensor without Magnetic field.
    REAL(KIND = DP) :: mob_bte(3, 3, nstemp)
    !! BTE mobility tensor without Magnetic field.
    REAL(KIND = DP) :: mobb_serta(3, 3, nstemp)
    !! SERTA mobility tensor with Magnetic field.
    REAL(KIND = DP) :: mobb_bte(3, 3, nstemp)
    !! BTE mobility tensor with Magnetic field.
    REAL(KIND = DP) :: sigma_serta_si(3, 3, nstemp)
    !! SERTA conductivity tensor without Magnetic field in SI unit.
    REAL(KIND = DP) :: sigma_bte_si(3, 3, nstemp)
    !! BTE conductivity tensor without Magnetic field in SI unit.
    REAL(KIND = DP) :: sigmab_serta_si(3, 3, nstemp)
    !! SERTA conductivity tensor with Magnetic field in SI unit.
    REAL(KIND = DP) :: sigmab_bte_si(3, 3, nstemp)
    !! BTE conductivity tensor with Magnetic field in SI unit.
    REAL(KIND = DP) :: hall_serta(3, 3, nstemp)
    !! Hall factor
    REAL(KIND = DP) :: hall(3, 3, nstemp)
    !! Hall factor
    REAL(KIND = DP) :: mob_hall(3, 3)
    !! Hall mobility
    REAL(KIND = DP) :: sigma_inv(3, 3, nstemp)
    !! Inverse conductivity tensor
    REAL(KIND = DP) :: mob_inv(3, 3, nstemp)
    !! Inverse mobility tensor
    !
    ! Convert carrier number in true carrier density
    carrier_density_cm(:) = zero
    DO itemp = 1, nstemp
      IF (system_2d == 'no') THEN
        inv_cell = 1.0d0 / omega
        carrier_density_cm(itemp) = carrier_density(itemp) * (inv_cell * (bohr2ang * ang2cm)**(-3.0d0))
      ELSE
        inv_cell = (at(3, 3) * alat) / omega
        carrier_density_cm(itemp) = carrier_density(itemp) * (inv_cell * (bohr2ang * ang2cm)**(-2.0d0))
      ENDIF
      !
      mob_serta(:, :, itemp)  = (sigma_serta(:, :, itemp) * electron_si * (bohr2ang * ang2cm)**2) &
                              / (carrier_density(itemp) * hbarJ)
      mob_bte(:, :, itemp)    = (sigma_bte(:, :, itemp) * electron_si * (bohr2ang * ang2cm)**2) &
                              / (carrier_density(itemp) * hbarJ)
      mobb_serta(:, :, itemp) = (sigmab_serta(:, :, itemp) * electron_si * (bohr2ang * ang2cm)**2) &
                              / (carrier_density(itemp) * hbarJ)
      mobb_bte(:, :, itemp)   = (sigmab_bte(:, :, itemp) * electron_si * (bohr2ang * ang2cm)**2) &
                              / (carrier_density(itemp) * hbarJ)
      !
      ! Convert conductivity tensor in SI units [Siemens m^-1=Coulomb s^-1 V^-1 m^-d ]
      ! in 3d: cm^2 s^-1 V^-1 * (cm ^-2  cmtom^-1 C) = Coulomb s^-1 V^-1
      IF (system_2d == 'no') THEN
        sigma_serta_si(:, :, itemp)  = mob_serta(:, :, itemp) * (electron_si * carrier_density_cm(itemp) * cm2m**(-1))
        sigma_bte_si(:, :, itemp)    = mob_bte(:, :, itemp) * (electron_si * carrier_density_cm(itemp) * cm2m**(-1))
        sigmab_serta_si(:, :, itemp) = mobb_serta(:, :, itemp) * (electron_si * carrier_density_cm(itemp) * cm2m**(-1))
        sigmab_bte_si(:, :, itemp)   = mobb_bte(:, :, itemp) * (electron_si * carrier_density_cm(itemp) * cm2m**(-1))
      ELSE
        sigma_serta_si(:, :, itemp)  = mob_serta(:, :, itemp) * (electron_si * carrier_density_cm(itemp))
        sigma_bte_si(:, :, itemp)    = mob_bte(:, :, itemp) * (electron_si * carrier_density_cm(itemp))
        sigmab_serta_si(:, :, itemp) = mobb_serta(:, :, itemp) * (electron_si * carrier_density_cm(itemp))
        sigmab_bte_si(:, :, itemp)   = mobb_bte(:, :, itemp) * (electron_si * carrier_density_cm(itemp))
      ENDIF
      !
      ! To make the diagonal of mobb zero.
      mobb_serta(:, :, itemp) = mobb_serta(:, :, itemp) - mob_serta(:, :, itemp)
      mobb_bte(:, :, itemp) = mobb_bte(:, :, itemp) - mob_bte(:, :, itemp)
      !

    ENDDO ! itemp
    !
    b_norm = SQRT(bfieldx**2 + bfieldy**2 + bfieldz**2)
    WRITE(stdout, '(/5x, a)') REPEAT('=',93)
    WRITE(stdout, '(5x, a)') 'BTE in the self-energy relaxation time approximation (SERTA)'
    WRITE(stdout, '(5x, a)') REPEAT('=',93)
    DO itemp = 1, nstemp
      etemp = gtemp(itemp)
      !
      IF (system_2d == 'no') THEN ! We suppose vacuum is in the z direction
        mob_inv(:, :, itemp) = matinv3(mob_serta(:, :, itemp))
      ELSE
        mob_serta(3, 3, :) = 1d0
        mob_inv(:, :, itemp) = matinv3(mob_serta(:, :, itemp))
        mob_inv(3, 3, :) = 0d0
      ENDIF
      hall_serta(:, :, itemp) = MATMUL(MATMUL(mob_inv(:, :, itemp), mobb_serta(:, :, itemp)), &
                          mob_inv(:, :, itemp)) / (b_norm * hbarJ ) * electron_si * (bohr2ang * ang2cm)**2
      !
      mob_hall(:, :)   = MATMUL(hall_serta(:, :, itemp),mob_serta(:, :, itemp))
      !
      !
      WRITE(stdout, '(5x,a)') ' '
      WRITE(stdout, '(5x,a,1f10.4,a)') 'Temperature: ', etemp * ryd2ev / kelvin2eV, ' K'
      IF (system_2d == 'no') THEN
        WRITE(stdout, '(5x,a,1E18.6)') 'Conductivity tensor without magnetic field | with magnetic field [Siemens/m]'
      ELSE
        WRITE(stdout, '(5x,a,1E18.6)') 'Conductivity tensor without magnetic field | with magnetic field [Siemens]'
      ENDIF
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_serta_si(1, :, itemp), '  |', sigmab_serta_si(1, :, itemp)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_serta_si(2, :, itemp), '  |', sigmab_serta_si(2, :, itemp)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_serta_si(3, :, itemp), '  |', sigmab_serta_si(3, :, itemp)
      !
      WRITE(stdout, '(5x,a)') 'Mobility tensor without magnetic field     | Hall mobility [cm^2/Vs]'
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_serta(1, :, itemp), '  |', mob_hall(1, :)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_serta(2, :, itemp), '  |', mob_hall(2, :)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_serta(3, :, itemp), '  |', mob_hall(3, :)
      !
      ! bfield is energy*sec/lenght**2, mobility is in cm**2 V**-1 sec**-1.
      ! To convert bfield to the same units of the mobility I do the same conversion as passing from the sigma
      ! tensor to the mobility: the energy is converted with electron_SI/hbarJ and the lenght**2 with (bohr2ang * ang2cm)**2
      WRITE(stdout, '(5x, a, 1E18.6)') 'Hall factor'
      WRITE(stdout, '(3x, 3E16.6)') hall_serta(1, :, itemp)
      WRITE(stdout, '(3x, 3E16.6)') hall_serta(2, :, itemp)
      WRITE(stdout, '(3x, 3E16.6)') hall_serta(3, :, itemp)
    ENDDO
    !
    WRITE(stdout, '(/5x, a)') REPEAT('=',93)
    WRITE(stdout, '(5x, a)') 'BTE'
    WRITE(stdout, '(5x, a)') REPEAT('=',93)
    DO itemp = 1, nstemp
      etemp = gtemp(itemp)
      !
      IF (system_2d == 'no') THEN
        mob_inv(:, :, itemp) = matinv3(mob_bte(:, :, itemp))
      ELSE ! We suppose vacuum is in the z direction
        mob_bte(3, 3, :) = 1d0
        mob_inv(:, :, itemp) = matinv3(mob_bte(:, :, itemp))
        mob_inv(3, 3, :) = 0d0
      ENDIF
      hall(:, :, itemp) = MATMUL(MATMUL(mob_inv(:, :, itemp), mobb_bte(:, :, itemp)), &
                          mob_inv(:, :, itemp)) / (b_norm * hbarJ ) * electron_si * (bohr2ang * ang2cm)**2
      mob_hall(:, :)   = MATMUL(hall(:, :, itemp),mob_bte(:, :, itemp))
      !
      WRITE(stdout, '(5x,a)') ' '
      WRITE(stdout, '(5x,a,1f10.4,a)') 'Temperature: ', etemp * ryd2ev / kelvin2eV, ' K'
      IF (system_2d == 'no') THEN
        WRITE(stdout, '(5x,a,1E18.6)') 'Conductivity tensor without magnetic field | with magnetic field [Siemens/m]'
      ELSE
        WRITE(stdout, '(5x,a,1E18.6)') 'Conductivity tensor without magnetic field | with magnetic field [Siemens]'
      ENDIF
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_bte_si(1, :, itemp), '  |', sigmab_bte_si(1, :, itemp)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_bte_si(2, :, itemp), '  |', sigmab_bte_si(2, :, itemp)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') sigma_bte_si(3, :, itemp), '  |', sigmab_bte_si(3, :, itemp)
      !
      WRITE(stdout, '(5x,a)') 'Mobility tensor without magnetic field     | Hall mobility [cm^2/Vs]'
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_bte(1, :, itemp), '  |', mob_hall(1, :)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_bte(2, :, itemp), '  |', mob_hall(2, :)
      WRITE(stdout, '(4x,3E14.5,a,3E14.5)') mob_bte(3, :, itemp), '  |', mob_hall(3, :)
      !
      ! bfield is energy*sec/lenght**2, mobility is in cm**2 V**-1 sec**-1.
      ! To convert bfield to the same units of the mobility I do the same conversion as passing from the sigma
      ! tensor to the mobility: the energy is converted with electron_SI/hbarJ and the lenght**2 with (bohr2ang * ang2cm)**2
      WRITE(stdout, '(5x, a, 1E18.6)') 'Hall factor'
      WRITE(stdout, '(3x, 3E16.6)') hall(1, :, itemp)
      WRITE(stdout, '(3x, 3E16.6)') hall(2, :, itemp)
      WRITE(stdout, '(3x, 3E16.6)') hall(3, :, itemp)
    ENDDO
    !
    ! SP: To be removed when adding to the main QE repo.
    !WRITE(stdout, '(5x, a)')'---------------------------------------'
    !WRITE(stdout, '(5x, a)')' '
    !WRITE(stdout, '(4x, 2f12.2, 2f12.3)') mob_serta(1, 1, 1), mob_bte(1, 1, 1), hall_serta(1, 2, 1), hall(1, 2, 1)
    !WRITE(stdout, '(4x, 2f12.2, 2f12.3)') mob_serta(1, 1, 2), mob_bte(1, 1, 2), hall_serta(1, 2, 2), hall(1, 2, 2)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE print_hall
    !-----------------------------------------------------------------------
  !-------------------------------------------------------------------------
  END MODULE printing
  !-------------------------------------------------------------------------
