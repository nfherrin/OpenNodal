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

  subroutine anal_keff (xkeff, anal_ref)
    IMPLICIT NONE
    real(kr8), intent(in) :: xkeff
    character(*), intent(in) :: anal_ref
    return
  endsubroutine anal_keff

  subroutine anal_frat (xflux, anal_ref)
    use errors_module, only : fatal_error
    use globals, only : assm_xs, print_log
    use string_module, only : str
    IMPLICIT NONE
    real(kr8), intent(in) :: xflux(:,:,:)
    character(*), intent(in) :: anal_ref

    real(kr8) :: Lx, Ly
    real(kr8) :: bsq_geo
    real(kr8) :: anal, calc, err

    select case (anal_ref)
      case ('2d2g')
        Lx = 1d2
        Ly = 1d2
        bsq_geo = (pi/Lx)**2 + (pi/Ly)**2
        anal = assm_xs(2)%sigma_scat(2,1) / &
          (assm_xs(2)%D(2)*bsq_geo + &
          (assm_xs(2)%sigma_t(2)-assm_xs(2)%sigma_scat(2,2))) ! TODO rewrite with sigma_r
        write(0,*) 'sigma_s(2,1)=', assm_xs(2)%sigma_scat(2,1)
        write(0,*) 'D(2)        =', assm_xs(2)%D(2)
        write(0,*) 'sigma_r(2)  =', assm_xs(2)%sigma_t(2)-assm_xs(2)%sigma_scat(2,2)
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
