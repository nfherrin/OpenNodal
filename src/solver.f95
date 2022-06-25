! TODO finish implementation
subroutine sor(aa, b, x, rank, omega, tol_inner_x, tol_inner_maxit)
  use precisions, only : ki4, kr8
  IMPLICIT NONE
  ! designed for a 5-stripe matrix with the diagonal in the fifth position
  real(kr8), intent(in) :: aa(:,:) ! (5, rank)
  real(kr8), intent(in) :: b(:)    ! (rank)
  real(kr8), intent(inout) :: x(:)    ! (rank)
  integer(ki4), intent(in) :: rank
  real(kr8), intent(in) :: omega
  real(kr8), intent(in) :: tol_inner_x
  integer(ki4), intent(in) :: tol_inner_maxit

  integer(ki4) :: i, iter

  do iter = 1,tol_inner_maxit
  enddo

endsubroutine sor

subroutine solver()
  use globals, only : stdout_unit, & ! TODO output to log as well
    core_x_size, core_y_size, &
    xkeff, xflux, tol_max_iter, tol_xkeff, tol_xflux
  use precisions, only : ki4, kr8
  USE globals, ONLY : print_log
  USE string_module, ONLY : str
  IMPLICIT NONE

  integer(ki4) :: iter
  real(kr8) :: conv_xflux, conv_xkeff

  real(kr8), allocatable :: bvec(:) ! (core_x_size*core_y_size*2)
  real(kr8), allocatable :: amat(:,:) ! (5, core_x_size*core_y_size*2)

  ! TODO build amat
  allocate(amat(5, core_x_size*core_y_size*2)) ! matrix is 5-stripe, 2 energy groups

  ! TODO create a routine to build bvec
  allocate(bvec(core_x_size*core_y_size*2)) ! 2 energy groups

  CALL print_log('Iter | Keff     | Conv_Keff | Conv_Flux')

  conv_xflux = 1d2*tol_xflux + 1d0
  conv_xkeff = 1d2*tol_xkeff + 1d0

  do iter = 1,tol_max_iter

    CALL print_log(TRIM(str(iter,4))//'   '//TRIM(str(xkeff,6,'F'))//'   '//TRIM(str(conv_xkeff,2)) &
      //'   '//TRIM(str(conv_xflux,2)))

    if ((conv_xflux < tol_xflux) .and. (conv_xkeff < tol_xkeff)) exit

  enddo

  CALL print_log('ITERATIONS FINISHED')
  CALL print_log('XKEFF = '//str(xkeff,16,'F'))

endsubroutine solver
