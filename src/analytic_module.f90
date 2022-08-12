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
    use globals , only : core_x_size, core_y_size, num_eg, print_log
    use errors_module, only : fatal_error
    IMPLICIT NONE
    real(kr8), intent(in) :: xkeff
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref


    real(kr8), allocatable :: cpy_flux(:,:,:)

    ! necessary so that other routines can renormalize to the referenced
    ! analytic solution
    allocate(cpy_flux(core_x_size, core_y_size, num_eg))
    cpy_flux = xflux

    call print_log('')
    call print_log('=== ANALYSIS ===')


    ! we can analyze keff, flux, and flux ratio
    ! but not all reference solutions may have all options available
    ! so we do some switching here
    select case (anal_ref)
      case ('2d2g')
        call anal_keff (xkeff, anal_ref)
        call anal_flux (xflux, anal_ref)
        call anal_frat (xflux, anal_ref)
      case default
        call fatal_error (&
          'analytic ratio calculation not implemented for reference problem: ' &
          // anal_ref)
    endselect

    deallocate(cpy_flux)

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
    call print_log('XKEFF_err  = ' // str(err , 6, 'f'))

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

  subroutine anal_flux (xflux, anal_ref)
    IMPLICIT NONE
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref
    return
  endsubroutine anal_flux

ENDMODULE analytic_module
