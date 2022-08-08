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
    real(kr8) :: conv_xflux, conv_xkeff,keff_old

    real(kr8), allocatable :: amat(:,:),flux_old(:,:,:) ! (core_x_size*core_y_size*num_eg, core_x_size*core_y_size*num_eg)

    ! build amatrix. This will change with dtilde
    ALLOCATE(amat(core_x_size*core_y_size*num_eg, core_x_size*core_y_size*num_eg)) ! eg by problem size
    ALLOCATE(flux_old(core_x_size,core_y_size,num_eg))
    CALL build_amatrix(amat)

    CALL print_log('Iter | Keff     | Conv_Keff | Conv_Flux')

    conv_xflux = 1d2*tol_xflux + 1d0
    conv_xkeff = 1d2*tol_xkeff + 1d0
    flux_old=xflux
    keff_old=xkeff

    do iter = 1,tol_max_iter

      !solve the cmfd problem for the given dtilde values
      CALL solve_cmfd(amat)
      conv_xkeff=ABS(xkeff-keff_old)/xkeff
      conv_xflux=calc_2norm(RESHAPE(xflux-flux_old,(/1/)))
      keff_old=xkeff
      flux_old=xflux
      CALL print_log(TRIM(str(iter,4))//'   '//TRIM(str(xkeff,6,'F'))//'   '//TRIM(str(conv_xkeff,2)) &
        //'   '//TRIM(str(conv_xflux,2)))

      if ((conv_xflux < tol_xflux) .and. (conv_xkeff < tol_xkeff)) exit

      CALL build_amatrix(amat)

    enddo

    CALL print_log('ITERATIONS FINISHED')
    CALL print_log('XKEFF = '//str(xkeff,16,'F'))

  endsubroutine solver

  !amatrix builder for direct solve for the CMFD portion
  SUBROUTINE build_amatrix(amatrix)
    REAL(kr8), INTENT(OUT) :: amatrix(:,:)
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
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(i,j,g)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(i,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_x(i,j,g)
          ENDIF
          !right term
          IF(i .NE. core_x_size)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i+1
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(i+1,j,g)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(i+1,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_x(i+1,j,g)
          ENDIF
          !above term
          IF(j .NE. 1)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j-2)*core_y_size+i
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(i,j,g)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(i,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_y(i,j,g)
          ENDIF
          !below term
          IF(j .NE. core_y_size)THEN
            neigh_idx=(g-1)*core_y_size*core_x_size+(j)*core_y_size+i
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(i,j+1,g)
            amatrix(cell_idx,neigh_idx)=amatrix(cell_idx,neigh_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(i,j+1,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(cell_idx,cell_idx)=amatrix(cell_idx,cell_idx)+dtilde_y(i,j+1,g)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    amatrix=amatrix/assm_pitch
  ENDSUBROUTINE build_amatrix

  !solves the cmfd problem based on the current amatrix
  SUBROUTINE solve_cmfd(amatrix)
    REAL(kr8), INTENT(IN) :: amatrix(:,:)
    REAL(kr8), ALLOCATABLE :: bvec(:),atemp(:,:),fiss_src_sum(:),flux_old(:,:,:)
    REAL(kr8) :: flux_err,keff_err,keff_old
    INTEGER :: i,j,g,prob_size,err1,cell_idx,ii
    INTEGER,ALLOCATABLE :: ipiv(:)

    prob_size=core_x_size*core_y_size*num_eg
    ! build bvec. This will change each inner iteration since we're doing power iterations
    ALLOCATE(bvec(prob_size))
    ALLOCATE(atemp(prob_size,prob_size))
    ALLOCATE(ipiv(prob_size))
    ALLOCATE(fiss_src_sum(2))
    ALLOCATE(flux_old(core_x_size,core_y_size,num_eg))
    flux_old=1.0D0
    fiss_src_sum=1.0D0
    xflux=xflux/SUM(xflux)
    keff_old=xkeff

    flux_err = 1d2*tol_xflux + 1d0
    keff_err = 1d2*tol_xkeff + 1d0

    DO WHILE(keff_err .GE. tol_xkeff .OR. flux_err .GE. tol_xflux)
      !build the bvec based on current keff and flux
      CALL build_bvec(bvec)

      atemp=amatrix
      CALL DGESV(prob_size,1,atemp,prob_size,ipiv,bvec,prob_size,err1)
      IF(err1 .NE. 0)CALL fatal_error('DGESV failed')

      CALL calc_fiss_src_sum(fiss_src_sum(1))
      !assign xflux to the new bvec
      DO j=1,core_y_size
        DO i=1,core_x_size
          DO g=1,num_eg
            cell_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i
            xflux(i,j,g)=bvec(cell_idx)
          ENDDO
        ENDDO
      ENDDO

      CALL calc_fiss_src_sum(fiss_src_sum(2))
      !calculate the eigenvalue
      keff_old=xkeff
      xkeff=xkeff*fiss_src_sum(2)/fiss_src_sum(1)
      keff_err=ABS(xkeff-keff_old)/xkeff
      !normalize xflux to 1
      xflux=xflux/SUM(xflux)
      flux_err=calc_2norm(RESHAPE(xflux-flux_old,(/1/)))
      flux_old=xflux
    ENDDO
  ENDSUBROUTINE solve_cmfd

  !bvector builder for current flux
  SUBROUTINE build_bvec(bvec)
    REAL(kr8), INTENT(OUT) :: bvec(:)
    INTEGER(ki4) :: i,j,g,gp,cell_idx,loc_id
    REAL(kr8) :: fiss_src,scat_src
    bvec=0.0D0
    DO j=1,core_y_size
      DO i=1,core_x_size
        !compute fission source for this cell
        loc_id=assm_map(i,j)
        fiss_src=0.0D0
        DO gp=1,num_eg
          fiss_src=fiss_src+xflux(i,j,gp)*assm_xs(loc_id)%nusigma_f(gp)
        ENDDO
        DO g=1,num_eg
          cell_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_y_size+i
          !set bvec to fission source
          bvec(cell_idx)=fiss_src*assm_xs(loc_id)%chi(g)/xkeff
          !and add scattering source
          DO gp=1,num_eg
            bvec(cell_idx)=bvec(cell_idx)+xflux(i,j,gp)*assm_xs(loc_id)%sigma_scat(g,gp)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE build_bvec

  !fission source calculator
  SUBROUTINE calc_fiss_src_sum(fiss_src)
    REAL(kr8), INTENT(OUT) :: fiss_src
    INTEGER(ki4) :: i,j,gp,loc_id
    fiss_src=0.0D0
    DO i=1,core_x_size
      DO j=1,core_y_size
        loc_id=assm_map(i,j)
        DO gp=1,num_eg
          fiss_src=fiss_src+xflux(i,j,gp)*assm_xs(loc_id)%nusigma_f(gp)
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE calc_fiss_src_sum

  !calculate a 2 norm
  FUNCTION calc_2norm(vec)
    REAL(kr8) :: calc_2norm
    REAL(kr8),INTENT(IN) :: vec(:)
    INTEGER :: i
    calc_2norm=0.0D0
    DO i=1,SIZE(vec)
      calc_2norm=vec(i)**2
    ENDDO
    calc_2norm=SQRT(calc_2norm)
  ENDFUNCTION calc_2norm
ENDMODULE solvers_module