! TODO finish implementation
!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Nodal solvers module.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE solvers_module
  USE globals
  USE errors_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solver
CONTAINS
!
!---------------------------------------------------------------------------------------------------
!> @brief This subroutine uses SOR to solve an Ax=b problem for x given A and b
!> @param aa - A matrix
!> @param b - RHS b vector
!> @param x - x solution vector
!> @param omega -
!> @param tol_inner_x -
!> @param tol_inner_maxit -
!>
  subroutine sor(aa, b, x, omega, tol_inner_x, tol_inner_maxit)
    use precisions, only : ki4, kr8
    IMPLICIT NONE
    ! designed for a 5-stripe matrix with the diagonal in the fifth position
    real(kr8), intent(in) :: aa(:,:) ! (5, rank)
    real(kr8), intent(in) :: b(:)    ! (rank)
    real(kr8), intent(inout) :: x(:)    ! (rank)
    real(kr8), intent(in) :: omega
    real(kr8), intent(in) :: tol_inner_x
    integer(ki4), intent(in) :: tol_inner_maxit

    integer(ki4) :: i, iter, rank

    rank=SIZE(b)

    do iter = 1,tol_inner_maxit
    enddo

  endsubroutine sor

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine solves the nodal problem
!>
  subroutine solver()
    use globals, only : stdout_unit, core_x_size, core_y_size, &
      xkeff, xflux, tol_max_iter, tol_xkeff, tol_xflux
    use precisions, only : ki4, kr8
    USE globals, ONLY : print_log
    USE string_module, ONLY : str
    IMPLICIT NONE

    integer(ki4) :: iter
    real(kr8) :: conv_xflux, conv_xkeff

    real(kr8), allocatable :: bvec(:) ! (core_x_size*core_y_size*2)
    real(kr8), allocatable :: amat(:,:) ! (core_x_size*core_y_size*num_eg, core_x_size*core_y_size*num_eg)

    !allocate dtildes
    ALLOCATE(dtilde_x(num_eg,core_x_size+1,core_y_size))
    ALLOCATE(dtilde_y(num_eg,core_x_size,core_y_size+1))
    dtilde_x=0.0
    dtilde_y=0.0

    ! TODO build amat
    ALLOCATE(amat(core_x_size*core_y_size*num_eg, core_x_size*core_y_size*num_eg)) ! eg by problem size
    CALL build_amatrix(amat)

    ! TODO create a routine to build bvec
    ALLOCATE(bvec(core_x_size*core_y_size*num_eg)) ! 2 energy groups

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

  !amatrix builder for direct solve
  SUBROUTINE build_amatrix(amatrix)
    REAL(kr8), INTENT(INOUT) :: amatrix(:,:)
    INTEGER(ki4) :: g,j,i,cell_idx,neigh_idx

    amatrix=0.0D0
    DO g=1,num_eg
      DO j=1,core_y_size
        DO i=1,core_x_size
          cell_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i
          !total xs term for the base of the A matrix diagonal
          amatrix(cell_idx,cell_idx)=assm_xs(assm_map(i,j))%sigma_t(g)*assm_pitch
          !left term
          IF(i .NE. 1)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i-1
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(g,i,j)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(g,i,j)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_x(g,i,j)
          ENDIF
          !right term
          IF(i .NE. core_x_size)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i+1
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(g,i+1,j)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(g,i+1,j)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_x(g,i+1,j)
          ENDIF
          !above term
          IF(j .NE. 1)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j-2)*core_y_size+i
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(g,i,j)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(g,i,j)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_y(g,i,j)
          ENDIF
          !below term
          IF(j .NE. core_y_size)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j)*core_y_size+i
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(g,i,j+1)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(g,i,j+1)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_y(g,i,j+1)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    amatrix=amatrix/assm_pitch
  ENDSUBROUTINE build_amatrix

ENDMODULE solvers_module