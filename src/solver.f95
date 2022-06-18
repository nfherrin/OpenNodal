! TODO sor solver

! TODO proper outputs... (log & stdout)

subroutine solver()
  use globals, only : stdout_unit, &
    core_x_size, core_y_size, &
    xkeff, xflux, tol_max_iter, tol_xkeff, tol_xflux
  use precisions, only : ki4, kr8
  IMPLICIT NONE

  integer(ki4) :: iter
  real(kr8) :: conv_xflux, conv_xkeff

  real(kr8), allocatable :: bvec(:) ! (core_x_size*core_y_size*2)
  real(kr8), allocatable :: amat(:,:) ! (5, core_x_size*core_y_size*2)

  ! TODO build amat
  allocate(amat(5, core_x_size*core_y_size*2)) ! matrix is 5-stripe, 2 energy groups

  ! TODO write a routine to build bvec
  allocate(bvec(core_x_size*core_y_size*2)) ! 2 energy groups

  write(stdout_unit, '(a)') 'Iter Keff Conv_Keff Conv_Flux'

  conv_xflux = 1d2*tol_xflux + 1d0
  conv_xkeff = 1d2*tol_xkeff + 1d0

  do iter = 1,tol_max_iter

    write(stdout_unit, '(i3, x, f8.6, x, es8.2, x, es8.2)') &
      iter, xkeff, conv_xkeff, conv_xflux

    if ((conv_xflux < tol_xflux) .and. (conv_xkeff < tol_xkeff)) exit

  enddo

  write(stdout_unit, '(a)') 'ITERATIONS FINISHED'
  write(stdout_unit, '(a,f18.16)') 'XKEFF = ', xkeff

endsubroutine solver
