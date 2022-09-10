MODULE analytic_module
USE precisions, ONLY : ki4, kr8
USE globals, ONLY : print_log, pi
USE errors_module, ONLY : fatal_error
USE string_module, ONLY : str
USE xs_types, ONLY : macro_assm_xs_type
IMPLICIT NONE

PRIVATE
PUBLIC :: anal

CONTAINS

  SUBROUTINE anal (xkeff, xflux, anal_ref, assm_map, assm_pitch, assm_xs, &
                   core_x_size, core_y_size, num_eg)
    REAL(kr8), INTENT(IN) :: xkeff
    REAL(kr8), INTENT(IN) :: xflux(:,:,:)
    CHARACTER(*), INTENT(IN) :: anal_ref
    INTEGER(ki4), INTENT(IN) :: assm_map(:,:)
    REAL(kr8), INTENT(IN) :: assm_pitch
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER(ki4), INTENT(IN) :: core_x_size, core_y_size, num_eg

    CALL print_log('')
    CALL print_log('=== ANALYSIS ===')

    ! we can analyze keff, flux, and flux ratio
    ! but not all reference solutions may have all options available
    ! so we do some switching here
    SELECT CASE (anal_ref)
      CASE ('2d2g')
        CALL anal_frat (xflux, anal_ref, assm_map, assm_pitch, assm_xs, &
                        core_x_size, core_y_size)
        CALL anal_keff (xkeff, anal_ref, assm_map, assm_pitch, assm_xs, &
                        core_x_size, core_y_size)
        CALL anal_flux (xflux, anal_ref, assm_pitch, core_x_size, core_y_size, &
                        num_eg)
      CASE DEFAULT
        CALL fatal_error (&
          'analysis not implemented for reference problem: ' // anal_ref)
      ENDSELECT

    RETURN
  ENDSUBROUTINE  anal

  REAL(kr8) FUNCTION calc_keff (anal_ref, assm_map, assm_pitch, assm_xs, &
                                core_x_size, core_y_size)
    CHARACTER(*), INTENT(IN) :: anal_ref
    INTEGER(ki4), INTENT(IN) :: assm_map(:,:)
    REAL(kr8), INTENT(IN) :: assm_pitch
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER, INTENT(IN) :: core_x_size, core_y_size

    REAL(kr8)    :: Lx, Ly
    REAL(kr8)    :: bsq_geo, frat
    INTEGER(ki4) :: mat_num

    SELECT CASE (anal_ref)
      CASE ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        bsq_geo = (pi/Lx)**2 + (pi/Ly)**2
        mat_num = assm_map(1,1)
        frat = calc_frat(anal_ref, assm_map, assm_pitch, assm_xs, &
                         core_x_size, core_y_size)
        calc_keff = (assm_xs(mat_num)%nusigma_f(1) + &
          assm_xs(mat_num)%nusigma_f(2)*frat) / &
          (assm_xs(mat_num)%D(1)*bsq_geo + assm_xs(mat_num)%sigma_r(1))
        RETURN
      CASE DEFAULT
        CALL fatal_error (&
          'analytic keff calculation not implemented for reference problem: ' &
          // anal_ref)
    ENDSELECT

  ENDFUNCTION calc_keff

  SUBROUTINE anal_keff (xkeff, anal_ref, assm_map, assm_pitch, assm_xs, &
                        core_x_size, core_y_size)
    REAL(kr8), INTENT(IN) :: xkeff
    CHARACTER(*), INTENT(IN) :: anal_ref
    INTEGER(ki4), INTENT(IN) :: assm_map(:,:)
    REAL(kr8), INTENT(IN) :: assm_pitch
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER(ki4), INTENT(IN) :: core_x_size, core_y_size

    REAL(kr8) :: anal, calc, err

    SELECT CASE (anal_ref)
      CASE ('2d2g')
        anal = calc_keff(anal_ref, assm_map, assm_pitch, assm_xs, &
                         core_x_size, core_y_size)
        calc = xkeff
      CASE DEFAULT
        CALL fatal_error (&
          'analytic keff calculation not implemented for reference problem: ' &
          // anal_ref)
    ENDSELECT

    err = (anal - calc)*1d5

    CALL print_log('XKEFF_anal = ' // str(anal, 10, 'f'))
    CALL print_log('XKEFF_calc = ' // str(calc, 10, 'f'))
    CALL print_log('XKEFF_err  = ' // TRIM(ADJUSTL(str(err , 6, 'f'))) // ' [pcm]')

    RETURN
  ENDSUBROUTINE anal_keff

  REAL(kr8) FUNCTION calc_frat (anal_ref, assm_map, assm_pitch, assm_xs, &
                                core_x_size, core_y_size)
    CHARACTER(*), INTENT(IN) :: anal_ref
    INTEGER(ki4), INTENT(IN) :: assm_map(:,:)
    REAL(kr8), INTENT(IN) :: assm_pitch
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER, INTENT(IN) :: core_x_size, core_y_size

    REAL(kr8) :: Lx, Ly
    REAL(kr8) :: bsq_geo
    INTEGER   :: mat_num

    SELECT CASE (anal_ref)
      CASE ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        bsq_geo = (pi/Lx)**2 + (pi/Ly)**2
        mat_num = assm_map(1,1)
        calc_frat = assm_xs(mat_num)%sigma_scat(2,1) / &
          (assm_xs(mat_num)%D(2)*bsq_geo + assm_xs(mat_num)%sigma_r(2))
        RETURN
      CASE DEFAULT
        CALL fatal_error (&
          'analytic ratio calculation not implemented for reference problem: ' &
          // anal_ref)
    ENDSELECT

  ENDFUNCTION calc_frat

  SUBROUTINE anal_frat (xflux, anal_ref, assm_map, assm_pitch, assm_xs, &
                        core_x_size, core_y_size)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:)
    CHARACTER(*), INTENT(IN) :: anal_ref
    INTEGER(ki4), INTENT(IN) :: assm_map(:,:)
    REAL(kr8), INTENT(IN) :: assm_pitch
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER(ki4), INTENT(IN) :: core_x_size, core_y_size

    REAL(kr8) :: anal, calc, err

    SELECT CASE (anal_ref)
      CASE ('2d2g')
        anal = calc_frat(anal_ref, assm_map, assm_pitch, assm_xs, core_x_size, core_y_size)
        calc = maxval(abs(xflux(:,:,2))) / maxval(abs(xflux(:,:,1))) ! TODO rewrite with inf norm
      CASE DEFAULT
        CALL fatal_error (&
          'analytic ratio calculation not implemented for reference problem: ' &
          // anal_ref)
        ENDSELECT

    err = anal - calc

    CALL print_log('RATIO_anal = ' // str(anal, 6, 'es'))
    CALL print_log('RATIO_calc = ' // str(calc, 6, 'es'))
    CALL print_log('RATIO_err  = ' // str(err , 6, 'es'))

    RETURN
  ENDSUBROUTINE anal_frat

  REAL(kr8) FUNCTION calc_flux(x, y, z, anal_ref, assm_pitch, core_x_size, &
                               core_y_size)
    REAL(kr8), INTENT(IN) :: x, y, z
    CHARACTER(*), INTENT(IN) :: anal_ref
    REAL(kr8), INTENT(IN) :: assm_pitch
    INTEGER(ki4), INTENT(IN) :: core_x_size, core_y_size

    REAL(kr8) :: Lx, Ly

    SELECT CASE (anal_ref)
      CASE ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        ! note: ONLY fast group (i.e. g=1 with normalization constant k1=1)
        calc_flux = 1d0*SIN(pi/Lx*x)*SIN(pi/Ly*y)
        RETURN
      CASE DEFAULT
        CALL fatal_error (&
          'analytic flux calculation not implemented for reference problem: ' &
          // anal_ref)
    ENDSELECT

  ENDFUNCTION calc_flux

  SUBROUTINE anal_flux (xflux, anal_ref, assm_pitch, core_x_size, core_y_size, &
                        num_eg)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:)
    CHARACTER(*), INTENT(IN) :: anal_ref
    REAL(kr8), INTENT(IN) :: assm_pitch
    INTEGER(ki4), INTENT(IN) :: core_x_size, core_y_size, num_eg

    INTEGER   :: i, j
    REAL(kr8) :: dx, dy, x, y
    REAL(kr8) :: diff, diff_two, diff_inf, diff_rms
    REAL(kr8), ALLOCATABLE :: cpy_flux(:,:,:)

    ! necessary so that other routines can renormalize to the referenced
    ! analytic solution
    ALLOCATE(cpy_flux(core_x_size, core_y_size, num_eg))
    cpy_flux = xflux
    cpy_flux = cpy_flux / MAXVAL(cpy_flux)

    diff_two = 0d0
    diff_inf = 0d0

    dx = assm_pitch
    dy = assm_pitch

    DO j = 1,core_y_size
      y = dy*(j-0.5d0)
      DO i = 1,core_x_size
        x = dx*(i-0.5d0)
        diff = calc_flux(x, y, 0d0, anal_ref, assm_pitch, core_x_size, core_y_size) - cpy_flux(i,j,1)
        diff_two = diff_two + diff**2
        diff_inf = MAX(diff_inf, ABS(diff))
      ENDDO
    ENDDO

    diff_rms = diff_two / SQRT(DBLE(core_x_size*core_y_size))

    CALL print_log('XFLUX_rms = ' // str(diff_rms))
    CALL print_log('XFLUX_two = ' // str(diff_two))
    CALL print_log('XFLUX_inf = ' // str(diff_inf))

    DEALLOCATE(cpy_flux)

    RETURN
  ENDSUBROUTINE anal_flux

ENDMODULE analytic_module
