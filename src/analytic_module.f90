MODULE analytic_module
USE precisions, ONLY : ki4, kr8
IMPLICIT NONE

PRIVATE
PUBLIC :: anal

! TODO move this to some sort of "constants" module
! (probably with Avogadro's number)

! Note: this is considered the best way of calculating pi as it is guaranteed to
! have the maximal precision. Alternatively, just type a lot of digits.
real(kr8), parameter :: pi = 4d0*atan(1d0)

CONTAINS

  subroutine anal (xkeff, xflux, anal_ref)
    use globals , only : print_log
    use errors_module, only : fatal_error
    IMPLICIT NONE
    real(kr8), intent(in) :: xkeff
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref


    call print_log('')
    call print_log('=== ANALYSIS ===')

    ! we can analyze keff, flux, and flux ratio
    ! but not all reference solutions may have all options available
    ! so we do some switching here
    select case (anal_ref)
      case ('2d2g')
        call anal_frat (xflux, anal_ref)
        call anal_keff (xkeff, anal_ref)
        call anal_flux (xflux, anal_ref)
      case default
        call fatal_error (&
          'analysis not implemented for reference problem: ' // anal_ref)
    endselect

    return
  endsubroutine anal

  real(kr8) function calc_keff (anal_ref)
    use globals, only : assm_xs, assm_map, assm_pitch, core_x_size, core_y_size
    use errors_module, only : fatal_error
    IMPLICIT NONE
    character(*), intent(in) :: anal_ref

    real(kr8) :: Lx, Ly
    real(kr8) :: bsq_geo
    integer   :: mat_num
    
    select case (anal_ref)
      case ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        bsq_geo = (pi/Lx)**2 + (pi/Ly)**2
        mat_num = assm_map(1,1)
        calc_keff = (assm_xs(mat_num)%nusigma_f(1) + &
          assm_xs(mat_num)%nusigma_f(2)*calc_frat(anal_ref)) / &
          (assm_xs(mat_num)%D(1)*bsq_geo + &
          (assm_xs(mat_num)%sigma_t(1) - assm_xs(mat_num)%sigma_scat(1,1)))
        return
      case default
        call fatal_error (&
          'analytic keff calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

  endfunction calc_keff

  subroutine anal_keff (xkeff, anal_ref)
    use errors_module, only : fatal_error
    use globals, only : print_log
    use string_module, only : str
    IMPLICIT NONE
    real(kr8), intent(in) :: xkeff
    character(*), intent(in) :: anal_ref

    real(kr8) :: anal, calc, err

    select case (anal_ref)
      case ('2d2g')
        anal = calc_keff(anal_ref)
        calc = xkeff
      case default
        call fatal_error (&
          'analytic keff calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

    err = (anal - calc)*1d5

    call print_log('XKEFF_anal = ' // str(anal, 10, 'f'))
    call print_log('XKEFF_calc = ' // str(calc, 10, 'f'))
    call print_log('XKEFF_err  = ' // TRIM(ADJUSTL(str(err , 6, 'f'))) // ' [pcm]')

    return
  endsubroutine anal_keff

  real(kr8) function calc_frat (anal_ref)
    use globals, only : assm_xs, assm_map, assm_pitch, core_x_size, core_y_size
    use errors_module, only : fatal_error
    IMPLICIT NONE
    character(*), intent(in) :: anal_ref

    real(kr8) :: Lx, Ly
    real(kr8) :: bsq_geo
    integer   :: mat_num

    select case (anal_ref)
      case ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        bsq_geo = (pi/Lx)**2 + (pi/Ly)**2
        mat_num = assm_map(1,1)
        calc_frat = assm_xs(mat_num)%sigma_scat(2,1) / &
          (assm_xs(mat_num)%D(2)*bsq_geo + &
          (assm_xs(mat_num)%sigma_t(2)-assm_xs(mat_num)%sigma_scat(2,2))) ! TODO rewrite with sigma_r
        return
      case default
        call fatal_error (&
          'analytic ratio calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

  endfunction calc_frat

  subroutine anal_frat (xflux, anal_ref)
    use errors_module, only : fatal_error
    use globals, only : print_log
    use string_module, only : str
    IMPLICIT NONE
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref

    real(kr8) :: anal, calc, err

    select case (anal_ref)
      case ('2d2g')
        anal = calc_frat(anal_ref)
        calc = maxval(abs(xflux(:,:,2))) / maxval(abs(xflux(:,:,1))) ! TODO rewrite with inf norm
      case default
        call fatal_error (&
          'analytic ratio calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

    err = anal - calc

    call print_log('RATIO_anal = ' // str(anal, 6, 'es'))
    call print_log('RATIO_calc = ' // str(calc, 6, 'es'))
    call print_log('RATIO_err  = ' // str(err , 6, 'es'))

    return
  endsubroutine anal_frat

  REAL(kr8) FUNCTION calc_flux(x, y, z, anal_ref)
    use globals, only : assm_pitch, core_x_size, core_y_size
    use errors_module, only : fatal_error
    REAL(kr8), INTENT(in) :: x, y, z
    CHARACTER(*), INTENT(in) :: anal_ref

    REAL(kr8) :: Lx, Ly

    select case (anal_ref)
      case ('2d2g')
        Lx = assm_pitch*core_x_size
        Ly = assm_pitch*core_y_size
        ! note: only fast group (i.e. g=1 with normalization constant k1=1)
        calc_flux = 1d0*SIN(pi/Lx*x)*SIN(pi/Ly*y)
        RETURN
      case default
        call fatal_error (&
          'analytic flux calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

  ENDFUNCTION calc_flux

  subroutine anal_flux (xflux, anal_ref)
    use globals , only : core_x_size, core_y_size, num_eg, print_log, assm_pitch
    use string_module, only : str
    IMPLICIT NONE
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref

    INTEGER   :: i, j
    REAL(kr8) :: dx, dy, x, y
    REAL(kr8) :: diff, diff_two, diff_inf, diff_rms
    REAL(kr8), ALLOCATABLE :: cpy_flux(:,:,:)

    ! necessary so that other routines can renormalize to the referenced
    ! analytic solution
    allocate(cpy_flux(core_x_size, core_y_size, num_eg))
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
        diff = calc_flux(x, y, 0d0, anal_ref) - cpy_flux(i,j,1)
        diff_two = diff_two + diff**2
        diff_inf = MAX(diff_inf, ABS(diff))
      ENDDO
    ENDDO

    diff_rms = diff_two / SQRT(DBLE(core_x_size*core_y_size))

    CALL print_log('XFLUX_rms = ' // str(diff_rms))
    CALL print_log('XFLUX_two = ' // str(diff_two))
    CALL print_log('XFLUX_inf = ' // str(diff_inf))

    DEALLOCATE(cpy_flux)

    return
  endsubroutine anal_flux

ENDMODULE analytic_module
