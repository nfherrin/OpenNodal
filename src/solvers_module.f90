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
  PUBLIC :: solver,solver_init

  CHARACTER(*), PARAMETER :: inner_solve_method = 'dgesv'

CONTAINS

  INTEGER PURE FUNCTION calc_idx(i, j, g)
    USE globals, ONLY : core_x_size, core_y_size
    INTEGER, INTENT(IN) :: i, j, g
    calc_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_x_size+i
    RETURN
  ENDFUNCTION calc_idx

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine solves the nodal problem
!>
  SUBROUTINE solver_init()
    INTEGER :: i,ii,j,jj,g,ip,jp
    INTEGER, ALLOCATABLE :: temp_map(:,:)

    IF(nsplit .GT. 1)THEN
      !in this case we are splitting the system and need to remake the assemply map
      !and recompute the core_x_size and core_y_size
      ALLOCATE(temp_map(core_x_size,core_y_size))
      temp_map=assm_map
      DEALLOCATE(assm_map)
      ALLOCATE(assm_map(nsplit*core_x_size,nsplit*core_y_size))
      assm_map=0
      jp=0
      DO j=1,core_y_size
        DO jj=1,nsplit
          jp=jp+1
          ip=0
          DO i=1,core_x_size
            DO ii=1,nsplit
              ip=ip+1
              assm_map(ip,jp)=temp_map(i,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      core_x_size=nsplit*core_x_size
      core_y_size=nsplit*core_y_size
      assm_pitch=assm_pitch/(nsplit*1.0D0)
      DEALLOCATE(temp_map)
    ENDIF
    xkeff = 1d0
    allocate(xflux(core_x_size,core_y_size,num_eg))
    xflux = 1d0
    !allocate dtildes
    ALLOCATE(dtilde_x(core_x_size+1,core_y_size,num_eg))
    ALLOCATE(dtilde_y(core_x_size,core_y_size+1,num_eg))
    dtilde_x=0.0
    dtilde_y=0.0
    !set initial CMFD boundary conditions for dtilde
    SELECT CASE(prob_sym)
      CASE('full') !vacuum on all sides
        dtilde_x(1,:,:)=0.5d0
        dtilde_x(core_x_size+1,:,:)=0.5d0
      CASE('qtr','half') !vacuum on right and bot sides or on right bot and top sides
        dtilde_x(core_x_size+1,:,:)=0.5d0
      CASE DEFAULT
        CALL fatal_error('Invalid symmetry option: '//TRIM(ADJUSTL(prob_sym)))
    ENDSELECT
    SELECT CASE(prob_sym)
      CASE('full','half') !vacuum on all sides or right and top/bot sides
        dtilde_y(:,1,:)=0.5d0
        dtilde_y(:,core_y_size+1,:)=0.5d0
      CASE('qtr')!vacuum on right and bot sides
        dtilde_y(:,core_y_size+1,:)=0.5d0
      CASE DEFAULT
        CALL fatal_error('Invalid symmetry option: '//TRIM(ADJUSTL(prob_sym)))
    ENDSELECT
  ENDSUBROUTINE solver_init
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
  SUBROUTINE sor(aa, b, flux, rank, omega, tol_inner_x, tol_inner_maxit)
    USE precisions, ONLY : ki4, kr8
    USE globals, ONLY : core_x_size, core_y_size, num_eg, print_log
    IMPLICIT NONE
    ! designed for a 5-stripe matrix with the diagonal in the fifth position
    REAL(kr8),    INTENT(in) :: aa(:,:) ! (5, rank)
    REAL(kr8),    INTENT(in) :: b(:)    ! (rank)
    REAL(kr8),    INTENT(inout) :: flux(:,:,:)    ! (rank)
    INTEGER(ki4), INTENT(IN) :: rank
    REAL(kr8),    INTENT(in) :: omega
    REAL(kr8),    INTENT(in) :: tol_inner_x
    INTEGER(ki4), INTENT(in) :: tol_inner_maxit

    INTEGER(ki4) :: i, j, g, iter
    INTEGER(ki4) :: cell_idx
    REAL(kr8)    :: zum, xdif, xmax
    REAL(kr8)    :: fold, fnew
    REAL(kr8)    :: converge

    ! TODO implement gg array so i dont have to track boundaries everywhere

    DO iter = 1,tol_inner_maxit
      xdif = 0d0
      xmax = 0d0
      DO g = 1,num_eg
        DO j = 1,core_y_size
          DO i = 1,core_x_size

            zum = 0d0
            cell_idx=calc_idx(i,j,g)

            ! west
            IF (i .NE. 1)           zum = zum + aa(2,cell_idx)*flux(i-1,j,g)
            ! east
            IF (i .NE. core_x_size) zum = zum + aa(3,cell_idx)*flux(i+1,j,g)
            ! south
            IF (j .NE. 1)           zum = zum + aa(4,cell_idx)*flux(i,j-1,g)
            ! north
            IF (j .NE. core_y_size) zum = zum + aa(5,cell_idx)*flux(i,j+1,g)

            zum  = (b(cell_idx) - zum) / aa(1,cell_idx)
            fold = flux(i,j,g)
            fnew = fold + omega*(zum - fold)
            xdif = MAX(xdif, ABS(fnew-fold))
            xmax = MAX(xmax, ABS(fnew))

            flux(i,j,g) = fnew

          ENDDO ! i = 1,core_x_size
        ENDDO ! j = 1,core_y_size
      ENDDO ! g = 1,num_eg
      converge = xdif/xmax
      IF (converge < tol_inner_x) RETURN
    ENDDO ! iter = 1,toL_inner_maxit

    CALL print_log('WARNING: SOR not converged')
  ENDSUBROUTINE sor

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine solves the nodal problem
!>
  subroutine solver()
    USE globals, ONLY : core_x_size, core_y_size, &
      xkeff, xflux, tol_max_iter, tol_xkeff, tol_xflux
    USE precisions, ONLY : ki4, kr8
    USE globals, ONLY : print_log
    USE string_module, ONLY : str
    USE linalg_module, ONLY : norm
    IMPLICIT NONE

    INTEGER(ki4) :: i, j, gp
    INTEGER(ki4) :: iter
    REAL(kr8) :: conv_xflux, conv_xkeff,keff_old

    REAL(kr8), ALLOCATABLE :: amat(:,:),flux_old(:,:,:) ! (core_x_size*core_y_size*num_eg, core_x_size*core_y_size*num_eg)
    REAL(kr8), ALLOCATABLE :: bvec(:)

    REAL(kr8) :: fiss_src_sum(2)

    INTEGER :: prob_size

    prob_size=core_x_size*core_y_size*num_eg

    ! build amatrix. This will change with dtilde
    ALLOCATE(amat(5, core_x_size*core_y_size*num_eg)) ! eg by problem size
    ALLOCATE(flux_old(core_x_size,core_y_size,num_eg))
    ALLOCATE(bvec(prob_size))

    ! TODO this will need to be done on every iteration
    CALL build_amatrix(amat)

    CALL print_log('Iter | Keff     | Conv_Keff | Conv_Flux')

    xflux = 1d0

    conv_xflux = 1d2*tol_xflux + 1d0
    conv_xkeff = 1d2*tol_xkeff + 1d0
    flux_old=xflux
    keff_old=xkeff

    DO iter = 1,tol_max_iter

      ! build the bvec based on current keff and flux
      CALL build_bvec(bvec)

      ! TODO implement inner tolerances
      CALL inner_solve(inner_solve_method, prob_size, 1d-3*MIN(tol_xflux,tol_xkeff), 10000, &
                       amat, bvec, xflux)

      CALL calc_fiss_src_sum(fiss_src_sum(2))

      IF (iter > 1) xkeff=xkeff*fiss_src_sum(2)/fiss_src_sum(1)
      conv_xkeff=ABS(xkeff-keff_old) ! absolute

      ! maximal relative change in flux
      ! TODO should be measured in power/fission source
      conv_xflux = 0d0
      DO gp = 1,num_eg
        DO j = 1,core_y_size
          DO i = 1,core_x_size
            conv_xflux = MAX(conv_xflux, &
              ABS(xflux(i,j,gp)-flux_old(i,j,gp))/xflux(i,j,gp))
          ENDDO ! i = 1,core_x_size
        ENDDO ! j = 1,core_y_size
      ENDDO ! gp = 1,num_eq

      ! store old data
      keff_old=xkeff
      flux_old=xflux
      fiss_src_sum(1) = fiss_src_sum(2)

      CALL print_log(TRIM(str(iter,4))//'   '//TRIM(str(xkeff,6,'F'))//'   '//TRIM(str(conv_xkeff,2)) &
        //'   '//TRIM(str(conv_xflux,2)))
      !CALL comp_dtilde()

      IF ((conv_xflux < tol_xflux) .AND. (conv_xkeff < tol_xkeff)) EXIT

    ENDDO ! iter = 1,tol_max_iter

    DO j=1,core_y_size
      !write(*,'(10000ES16.8)')xflux(:,j,1)
    ENDDO

    DEALLOCATE(amat, flux_old, bvec)

    CALL print_log('ITERATIONS FINISHED')
    CALL print_log('XKEFF = '//str(xkeff,16,'F'))

  ENDSUBROUTINE solver

  !amatrix builder
  SUBROUTINE build_amatrix(amatrix)
    REAL(kr8), INTENT(OUT) :: amatrix(:,:)
    INTEGER(ki4) :: g,j,i,cell_idx

    ! (1   , 2   , 3    , 4     , 5   )
    ! (diag, left, right, below, above)
    ! (diag, left, right, down , up   )
    ! (diag, west, east , south, north)

    amatrix=0.0D0
    DO g=1,num_eg
      DO j=1,core_y_size
        DO i=1,core_x_size
          cell_idx=calc_idx(i,j,g)
          !total xs term for the base of the A matrix diagonal
          amatrix(1,cell_idx)=assm_xs(assm_map(i,j))%sigma_t(g)*assm_pitch
          !left term
          IF(i .NE. 1)THEN
            amatrix(1,cell_idx)=amatrix(1,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(i,j,g)
            amatrix(2,cell_idx)=amatrix(2,cell_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i-1,j))%D(g))**(-1)-dtilde_x(i,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx)=amatrix(1,cell_idx)+dtilde_x(i,j,g)
          ENDIF
          !right term
          IF(i .NE. core_x_size)THEN
            amatrix(1,cell_idx)=amatrix(1,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(i+1,j,g)
            amatrix(3,cell_idx)=amatrix(3,cell_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i+1,j))%D(g))**(-1)+dtilde_x(i+1,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx)=amatrix(1,cell_idx)+dtilde_x(i+1,j,g)
          ENDIF
          !below term
          IF(j .NE. 1)THEN
            amatrix(1,cell_idx)=amatrix(1,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(i,j,g)
            amatrix(4,cell_idx)=amatrix(4,cell_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j-1))%D(g))**(-1)-dtilde_y(i,j,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx)=amatrix(1,cell_idx)+dtilde_y(i,j,g)
          ENDIF
          !above term
          IF(j .NE. core_y_size)THEN
            amatrix(1,cell_idx)=amatrix(1,cell_idx)&
              +2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(i,j+1,g)
            amatrix(5,cell_idx)=amatrix(5,cell_idx)&
              -2.0D0*(assm_pitch/assm_xs(assm_map(i,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j+1))%D(g))**(-1)+dtilde_y(i,j+1,g)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx)=amatrix(1,cell_idx)+dtilde_y(i,j+1,g)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    amatrix=amatrix/assm_pitch
  ENDSUBROUTINE build_amatrix

  SUBROUTINE stripe_to_dense(stripe, dense)
    USE globals, ONLY : num_eg, core_x_size, core_y_size
    REAL(kr8), INTENT(IN)  :: stripe(:,:)
    REAL(kr8), INTENT(OUT) :: dense(:,:)

    INTEGER :: i, j, g
    INTEGER :: cell_idx, neigh_idx

    dense = 0d0

    DO g = 1,num_eg
      DO j = 1,core_y_size
        DO i = 1,core_x_size

          cell_idx=calc_idx(i,j,g)

          ! self
          dense(cell_idx,cell_idx) = stripe(1,cell_idx)

          ! west
          IF (i .NE. 1) THEN
            neigh_idx=calc_idx(i-1,j,g)
            dense(cell_idx,neigh_idx) = stripe(2,cell_idx)
          ENDIF

          ! east
          IF (i .NE. core_x_size) THEN
            neigh_idx=calc_idx(i+1,j,g)
            dense(cell_idx,neigh_idx) = stripe(3,cell_idx)
          ENDIF

          ! south
          IF (j .NE. 1) THEN
            neigh_idx=calc_idx(i,j-1,g)
            dense(cell_idx,neigh_idx) = stripe(4,cell_idx)
          ENDIF

          ! north
          IF (j .NE. core_y_size) THEN
            neigh_idx=calc_idx(i,j+1,g)
            dense(cell_idx,neigh_idx) = stripe(5,cell_idx)
          ENDIF

        ENDDO ! i = 1,core_x_size
      ENDDO ! j = 1,core_y_size
    ENDDO ! g = 1,num_eq

    RETURN
  ENDSUBROUTINE stripe_to_dense

  ! solve inner linear system of equations (allows switching solver)
  SUBROUTINE inner_solve(method, rank, tol_inner_x, tol_inner_maxit, &
                         amat, bvec, xflux)
    USE globals, ONLY : num_eg, core_x_size, core_y_size
    CHARACTER(*), INTENT(IN)    :: method
    INTEGER,      INTENT(IN)    :: rank
    REAL(kr8),    INTENT(IN)    :: tol_inner_x
    INTEGER,      INTENT(IN)    :: tol_inner_maxit
    REAL(kr8),    INTENT(IN)    :: amat(:,:)
    REAL(kr8),    INTENT(IN)    :: bvec(:)
    REAL(kr8),    INTENT(INOUT) :: xflux(:,:,:)

    ! for dgesv
    REAL(kr8), ALLOCATABLE :: atemp(:,:)
    INTEGER(kr4), ALLOCATABLE :: ipiv(:)
    INTEGER :: ierr
    INTEGER :: cell_idx
    INTEGER :: i, j, g

    SELECT CASE (method)
      CASE ('dgesv')
        ALLOCATE(atemp(rank,rank))
        ALLOCATE(ipiv(rank))
        call stripe_to_dense(amat, atemp)
        CALL DGESV(rank,1,atemp,rank,ipiv,bvec,rank,ierr)
        IF(ierr .NE. 0)CALL fatal_error('DGESV failed')
        !assign xflux to the new bvec
        DO g=1,num_eg
          DO j=1,core_y_size
            DO i=1,core_x_size
              cell_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_x_size+i
              xflux(i,j,g)=MAX(bvec(cell_idx),0d0)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(atemp, ipiv)
      CASE ('sor')
        ! TODO set omega better than this
        CALL sor(amat, bvec, xflux, rank, 1.2d0, tol_inner_x, tol_inner_maxit)
      CASE DEFAULT
        call fatal_error('selected inner_solve method not implemented')
    ENDSELECT

    RETURN
  ENDSUBROUTINE inner_solve

  !bvector builder for current flux
  SUBROUTINE build_bvec(bvec)
    REAL(kr8), INTENT(OUT) :: bvec(:)
    INTEGER(ki4) :: i,j,g,gp,cell_idx,loc_id
    REAL(kr8) :: fiss_src
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
          cell_idx=(g-1)*core_y_size*core_x_size+(j-1)*core_x_size+i
          !set bvec to fission source
          bvec(cell_idx)=fiss_src*assm_xs(loc_id)%chi(g)/xkeff
          ! TODO scattering should be in LHS
          ! TODO scattering transposed...
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
    DO gp=1,num_eg
      DO j=1,core_y_size
        DO i=1,core_x_size
          loc_id=assm_map(i,j)
          fiss_src=fiss_src+xflux(i,j,gp)*assm_xs(loc_id)%nusigma_f(gp)
        ENDDO ! i = 1,core_x_size
      ENDDO ! j = 1,core_y_size
    ENDDO ! gp = 1,num_eg
  ENDSUBROUTINE calc_fiss_src_sum

  SUBROUTINE comp_dtilde()
    REAL(kr8) :: s_bar_x(core_x_size,core_y_size,num_eg),s_bar_y(core_x_size,core_y_size,num_eg)
    CALL comp_s(s_bar_x,s_bar_y)
    STOP 'comp_dtilde not yet complete'
  ENDSUBROUTINE comp_dtilde

  SUBROUTINE comp_s(s_bar_x,s_bar_y)
    REAL(kr8),INTENT(OUT) :: s_bar_x(:,:,:),s_bar_y(:,:,:)
    REAL(kr8) :: j_x(core_x_size+1,core_y_size,num_eg),j_y(core_x_size,core_y_size+1,num_eg)
    REAL(kr8) :: l_bar_x(core_x_size,core_y_size,num_eg),l_bar_y(core_x_size,core_y_size,num_eg)
    INTEGER :: i,j,g

    !compute x direction currents at each face
    j_x=0.0D0
    DO i=1,core_x_size+1
      DO j=1,core_y_size
        DO g=1,num_eg
          IF(i .EQ. 1)THEN
            j_x(i,j,g)=-dtilde_x(i,j,g)*xflux(i,j,g)
          ELSEIF(i .EQ. core_x_size+1)THEN
            j_x(i,j,g)=dtilde_x(i,j,g)*xflux(i-1,j,g)
          ELSE !not a boundary
            j_x(i,j,g)=2.0D0*(assm_pitch/assm_xs(assm_map(i-1,j))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j))%D(g))*(xflux(i-1,j,g)-xflux(i,j,g))&
              +dtilde_x(i,j,g)*(xflux(i-1,j,g)+xflux(i,j,g))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !compute y direction currents at each face
    j_y=0.0D0
    DO i=1,core_x_size
      DO j=1,core_y_size+1
        DO g=1,num_eg
          IF(j .EQ. 1)THEN
            j_y(i,j,g)=-dtilde_y(i,j,g)*xflux(i,j,g)
          ELSEIF(j .EQ. core_y_size+1)THEN
            j_y(i,j,g)=dtilde_y(i,j,g)*xflux(i,j-1,g)
          ELSE !not a boundary
            j_y(i,j,g)=2.0D0*(assm_pitch/assm_xs(assm_map(i,j-1))%D(g)&
              +assm_pitch/assm_xs(assm_map(i,j))%D(g))*(xflux(i,j-1,g)-xflux(i,j,g))&
              +dtilde_y(i,j,g)*(xflux(i,j-1,g)+xflux(i,j,g))
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !!x direction currents
    !write(*,*)'j_x'
    !DO j=1,core_y_size
    !  WRITE(*,'(10000ES16.8)')j_x(:,j,1)
    !ENDDO
    !!y direction currents
    !write(*,*)'j_y'
    !DO j=1,core_y_size+1
    !  WRITE(*,'(10000ES16.8)')j_y(:,j,1)
    !ENDDO

    !compute l_bar_x values
    DO i=1,core_x_size
      DO j=1,core_y_size
        DO g=1,num_eg
          l_bar_x(i,j,g)=j_x(i+1,j,g)-j_x(i,j,g)
        ENDDO
      ENDDO
    ENDDO
    !compute l_bar_y values
    DO i=1,core_x_size
      DO j=1,core_y_size
        DO g=1,num_eg
          l_bar_y(i,j,g)=j_y(i,j+1,g)-j_y(i,j,g)
        ENDDO
      ENDDO
    ENDDO

    !!x direction l
    !write(*,*)'l_x'
    !DO j=1,core_y_size
    ! WRITE(*,'(10000ES16.8)')l_bar_x(:,j,1)
    !ENDDO
    !!y direction l
    !write(*,*)'l_y'
    !DO j=1,core_y_size
    ! WRITE(*,'(10000ES16.8)')l_bar_y(:,j,1)
    !ENDDO

    !compute s_bar_x values
    DO i=1,core_x_size
      DO j=1,core_y_size
        DO g=1,num_eg
          s_bar_x(i,j,g)=l_bar_y(i,j,g)/assm_pitch
        ENDDO
      ENDDO
    ENDDO
    !compute s_bar_y values
    DO i=1,core_x_size
      DO j=1,core_y_size
        DO g=1,num_eg
          s_bar_y(i,j,g)=l_bar_x(i,j,g)/assm_pitch
        ENDDO
      ENDDO
    ENDDO

    !x direction l
    write(*,*)'s_x'
    DO j=1,core_y_size
      WRITE(*,'(10000ES16.8)')s_bar_x(:,j,1)
    ENDDO
    !y direction l
    write(*,*)'s_y'
    DO j=1,core_y_size
      WRITE(*,'(10000ES16.8)')s_bar_y(:,j,1)
    ENDDO
  ENDSUBROUTINE comp_s
ENDMODULE solvers_module
