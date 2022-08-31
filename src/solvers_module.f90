!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Nodal solvers module.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE solvers_module
  USE globals
  USE errors_module
  USE xs_types
  USE string_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solver,solver_init

  CHARACTER(*), PARAMETER :: inner_solve_method = 'sor'

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine computes the cell id based on the 2D map
!> @param i - x location index
!> @param j - y location index
!> @param core_x_size - core size in the x direction
!>
  INTEGER(ki4) PURE FUNCTION calc_idx(i, j, core_x_size)
    INTEGER, INTENT(IN) :: i, j, core_x_size
    calc_idx=(j-1)*core_x_size+i
    RETURN
  ENDFUNCTION calc_idx

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine initializes values for the nodal solver (i.e. BC/splitting stuff)
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_map - assembly map
!> @param refl_mat - reflector material
!> @param num_assm_reg - number of unique assemblies
!> @param assm_xs - assembly level cross sections
!> @param ax_buckle - axial buckling for 2D problems
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!> @param nsplit - nsplit value, for decomposing nodes into nsplit number of sub-nodes
!> @param assm_pitch - assembly pitch
!> @param xkeff - eigenvalue
!> @param xflux - scalar flux
!> @param dtilde_x - current correction cross section used in the x direction
!> @param dtilde_y - current correction cross section used in the y direction
!> @param bc_opt - boundary condition option
!> @param albedos - albedo boundary conditions
!> @param prob_sym - problem symmetry
!>
  SUBROUTINE solver_init(core_x_size,core_y_size,num_eg,assm_map,refl_mat,num_assm_reg,assm_xs, &
                          ax_buckle,h_x,h_y,nsplit,assm_pitch,xkeff,xflux,dtilde_x,dtilde_y,bc_opt, &
                          albedos,prob_sym)
    !input/output variables
    INTEGER(ki4), INTENT(INOUT) :: core_x_size,core_y_size
    INTEGER(ki4), INTENT(INOUT), ALLOCATABLE :: assm_map(:,:)
    INTEGER(ki4), INTENT(IN) :: num_eg,refl_mat,num_assm_reg,nsplit
    TYPE(macro_assm_xs_type), INTENT(INOUT) :: assm_xs(:)
    REAL(kr8), INTENT(IN) :: ax_buckle,albedos(:)
    REAL(kr8), INTENT(INOUT) :: assm_pitch
    REAL(kr8), INTENT(OUT) :: xkeff
    REAL(kr8), INTENT(INOUT), ALLOCATABLE :: h_x(:),h_y(:)
    REAL(kr8), INTENT(OUT), ALLOCATABLE :: xflux(:,:,:),dtilde_x(:,:,:),dtilde_y(:,:,:)
    CHARACTER(*), INTENT(IN) :: bc_opt,prob_sym
    !local variables
    INTEGER(ki4) :: i,ii,j,jj,g,ip,jp,temp_map(core_x_size,core_y_size)
    INTEGER(ki4) :: temp_hx(core_x_size),temp_hy(core_y_size)
    REAL(kr8) :: gamma(num_eg)

    !set gaps in ragged core equal to the reflector material
    IF(MINVAL(assm_map) .LE. 0 .AND. refl_mat .EQ. 0)&
      CALL fatal_error('Ragged core specified, but no reflector material given')
    DO j=1,core_y_size
      DO i=1,core_x_size
        IF(assm_map(i,j) .LE. 0)assm_map(i,j)=refl_mat
      ENDDO
    ENDDO

    !create the removal cross section by adding buckling and subtracting self-scatter
    DO i=1,num_assm_reg
      DO g=1,num_eg
        assm_xs(i)%sigma_r(g)=assm_xs(i)%sigma_t(g)+assm_xs(i)%D(g)*ax_buckle &
          -assm_xs(i)%sigma_scat(g,g)
      ENDDO
    ENDDO

    IF(nsplit .GT. 1)THEN
      !in this case we are splitting the system and need to remake the assembly map
      !and recompute the core_x_size and core_y_size
      temp_map=assm_map
      temp_hx=h_x
      temp_hy=h_y
      DEALLOCATE(assm_map,h_x,h_y)
      ALLOCATE(assm_map(nsplit*core_x_size,nsplit*core_y_size))
      ALLOCATE(h_x(nsplit*core_x_size),h_y(nsplit*core_y_size))
      assm_map=0
      jp=0
      DO j=1,core_y_size
        DO jj=1,nsplit
          jp=jp+1
          h_y(jp)=temp_hy(j)/(nsplit*1.0D0)
          ip=0
          DO i=1,core_x_size
            DO ii=1,nsplit
              ip=ip+1
              h_x(ip)=temp_hx(i)/(nsplit*1.0D0)
              assm_map(ip,jp)=temp_map(i,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      core_x_size=nsplit*core_x_size
      core_y_size=nsplit*core_y_size
      assm_pitch=assm_pitch/(nsplit*1.0D0)
    ENDIF
    xkeff = 1d0
    ALLOCATE(xflux(core_x_size,core_y_size,num_eg))
    xflux = 1d0
    !allocate dtildes
    ALLOCATE(dtilde_x(core_x_size+1,core_y_size,num_eg))
    ALLOCATE(dtilde_y(core_x_size,core_y_size+1,num_eg))
    dtilde_x=0.0D0
    dtilde_y=0.0D0

    !set the gamma value
    SELECTCASE(bc_opt)
      CASE('vac','vacuum') !vacuum gamma is 2
        gamma=2.0D0
      CASE('zero') !zero bc gama is 0
        gamma=0.0D0
      CASE('reflective') !nothing to do for now, we'll set all dtildes to zero later...
      CASE('reflector')
        !finite difference reflector, does not account for out-scattering
        DO g=1,num_eg
          gamma(g)=SQRT(assm_xs(refl_mat)%D(g)/assm_xs(refl_mat)%sigma_t(g))/assm_xs(refl_mat)%D(g)
        ENDDO
      CASE('albedo') !will eventually be supported so give a debugging stop, not a fatal error
        DO g=1,num_eg
          gamma(g)=2.0D0+4.0D0/(albedos(g)**(-1)-1)
        ENDDO
    ENDSELECT

    !set CMFD boundary conditions for dtilde based on symmetry
    DO g=1,num_eg
      DO j=1,core_y_size
        SELECT CASE(prob_sym)
          CASE('full') !vacuum on all sides
                dtilde_x(1,j,g)=(0.5D0*h_x(1)/assm_xs(assm_map(1,j))%D(g)+gamma(g))**(-1)
                dtilde_x(core_x_size+1,j,g)=(0.5D0*h_x(core_x_size)/assm_xs(assm_map(core_x_size,j))%D(g)+gamma(g))**(-1)
          CASE('qtr','half') !vacuum on right and bot sides or on right bot and top sides
            dtilde_x(core_x_size+1,j,g)=(0.5D0*h_x(core_x_size)/assm_xs(assm_map(core_x_size,j))%D(g)+gamma(g))**(-1)
          CASE DEFAULT
            CALL fatal_error('Invalid symmetry option: '//TRIM(ADJUSTL(prob_sym)))
        ENDSELECT
      ENDDO
      DO i=1,core_x_size
        SELECT CASE(prob_sym)
          CASE('full','half') !vacuum on all sides or right and top/bot sides
            dtilde_y(i,1,g)=(0.5D0*h_y(1)/assm_xs(assm_map(i,1))%D(g)+gamma(g))**(-1)
            dtilde_y(i,core_y_size+1,g)=(0.5D0*h_y(core_y_size)/assm_xs(assm_map(i,core_y_size))%D(g)+gamma(g))**(-1)
          CASE('qtr')!vacuum on right and bot sides
            dtilde_y(i,core_y_size+1,g)=(0.5D0*h_y(core_y_size)/assm_xs(assm_map(i,core_y_size))%D(g)+gamma(g))**(-1)
          CASE DEFAULT
            CALL fatal_error('Invalid symmetry option: '//TRIM(ADJUSTL(prob_sym)))
        ENDSELECT
      ENDDO
    ENDDO

    !if we have reflective just set back to 0
    SELECTCASE(bc_opt)
      CASE('vac','vacuum','albedo','reflector') !nothing to do, already taken care of in previous block
      CASE('reflective') !all dtildes should start as zero since that is for reflective at boundary
        dtilde_x=0.0D0
        dtilde_y=0.0D0
    ENDSELECT
  ENDSUBROUTINE solver_init
!
!---------------------------------------------------------------------------------------------------
!> @brief This subroutine uses SOR to solve an Ax=b problem for x given A and b
!> @param aa - A matrix
!> @param b - RHS b vector
!> @param xflux - scalar flux
!> @param rank - system size
!> @param omega - Over-relaxation factor
!> @param tol_inner_x - Inner tolerance
!> @param tol_inner_maxit - Number of max iterations
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!>
  SUBROUTINE sor(aa, b, flux, rank, omega, tol_inner_x, tol_inner_maxit,core_x_size,core_y_size)
    IMPLICIT NONE
    ! designed for a 5-stripe matrix with the diagonal in the fifth position
    REAL(kr8),    INTENT(IN) :: aa(:,:) ! (5, rank)
    REAL(kr8),    INTENT(IN) :: b(:)    ! (rank)
    REAL(kr8),    INTENT(INOUT) :: flux(:,:)    ! (rank)
    INTEGER(ki4), INTENT(IN) :: rank,core_x_size,core_y_size
    REAL(kr8),    INTENT(IN) :: omega
    REAL(kr8),    INTENT(IN) :: tol_inner_x
    INTEGER(ki4), INTENT(IN) :: tol_inner_maxit

    INTEGER(ki4) :: i, j, iter
    INTEGER(ki4) :: cell_idx
    REAL(kr8)    :: zum, xdif, xmax
    REAL(kr8)    :: fold, fnew
    REAL(kr8)    :: converge

    ! TODO implement gg array so i dont have to track boundaries everywhere

    DO iter = 1,tol_inner_maxit
      xdif = 0d0
      xmax = 0d0
      DO j = 1,core_y_size
        DO i = 1,core_x_size

          zum = 0d0
          cell_idx=calc_idx(i,j,core_x_size)

          ! west
          IF (i .NE. 1)           zum = zum + aa(2,cell_idx)*flux(i-1,j)
          ! east
          IF (i .NE. core_x_size) zum = zum + aa(3,cell_idx)*flux(i+1,j)
          ! south
          IF (j .NE. 1)           zum = zum + aa(4,cell_idx)*flux(i,j-1)
          ! north
          IF (j .NE. core_y_size) zum = zum + aa(5,cell_idx)*flux(i,j+1)

          zum  = (b(cell_idx) - zum) / aa(1,cell_idx)
          fold = flux(i,j)
          fnew = fold + omega*(zum - fold)
          xdif = MAX(xdif, ABS(fnew-fold))
          xmax = MAX(xmax, ABS(fnew))

          flux(i,j) = fnew

        ENDDO ! i = 1,core_x_size
      ENDDO ! j = 1,core_y_size
      converge = xdif/xmax
      IF (converge < tol_inner_x) RETURN
    ENDDO ! iter = 1,toL_inner_maxit

    CALL print_log('WARNING: SOR not converged')
  ENDSUBROUTINE sor

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine solves the nodal problem
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param tol_xflux - flux convergence tolerance
!> @param tol_xkeff - keff convergence tolerance
!> @param xflux - scalar flux
!> @param xkeff - eigenvalue
!> @param tol_max_iter - maximum number of iterations
!> @param nodal_method - nodal method option
!> @param assm_map - assembly map
!> @param assm_xs - assembly level cross sections
!> @param dtilde_x - current correction cross section used in the x direction
!> @param dtilde_y - current correction cross section used in the y direction
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!>
  subroutine solver(core_x_size,core_y_size,num_eg,tol_xflux,tol_xkeff,xflux,xkeff,tol_max_iter, &
                    nodal_method,assm_map,assm_xs,dtilde_x,dtilde_y,h_x,h_y)
    INTEGER, INTENT(IN) :: core_x_size,core_y_size,num_eg,tol_max_iter,assm_map(:,:)
    REAL(kr8), INTENT(IN) :: tol_xflux,tol_xkeff,h_x(:),h_y(:)
    REAL(kr8), INTENT(INOUT) :: xflux(:,:,:),xkeff,dtilde_x(:,:,:),dtilde_y(:,:,:)
    CHARACTER(*), INTENT(IN) :: nodal_method
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: i, j, gp, g
    INTEGER(ki4) :: iter
    REAL(kr8) :: conv_xflux, conv_xkeff,ckeff_old
    REAL(kr8), ALLOCATABLE :: amat(:,:,:),flux_old(:,:,:),amat_base(:,:,:)
    REAL(kr8), ALLOCATABLE :: bvec(:,:)
    REAL(kr8) :: fiss_src_sum(2)
    !computed keff taking account for the weiland shift
    REAL(kr8) :: ckeff
    !lambda prime for weilandt shift
    REAL(kr8) :: lambda_p
    INTEGER :: prob_size

    prob_size=core_x_size*core_y_size

    ! build amatrix. This will change with dtilde
    ALLOCATE(amat(5, core_x_size*core_y_size,num_eg)) ! eg by problem size
    ALLOCATE(flux_old(core_x_size,core_y_size,num_eg))
    ALLOCATE(bvec(prob_size,num_eg))

    !build that amatrix
    CALL build_amatrix(amat,core_x_size,core_y_size,num_eg,assm_xs,assm_map,h_x,h_y,dtilde_x, &
                        dtilde_y)
    !set the base amatrix if we're doing weilandt shift
    IF(dl_weilandt .GT. 0.0D0)THEN
      ALLOCATE(amat_base(5, core_x_size*core_y_size,num_eg))
      amat_base=amat
    ENDIF

    CALL print_log(' Iter | Keff     | Conv_Keff | Conv_Flux')

    conv_xflux = 1d2*tol_xflux + 1d0
    conv_xkeff = 1d2*tol_xkeff + 1d0
    flux_old=xflux
    ckeff=xkeff
    ckeff_old=ckeff


    DO iter = 1,tol_max_iter

      ! build the bvec based on current keff and flux
      CALL build_bvec(bvec,core_x_size,core_y_size,num_eg,assm_xs,xkeff,xflux,assm_map)
      DO g=1,num_eg
        CALL add_downscatter(bvec,g,core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux)
      ENDDO

      IF(dl_weilandt .GT. 0.0D0 .AND. iter > 10)THEN
        lambda_p=ckeff_old+dl_weilandt
        amat=weilandt_shift_amatrix(amat_base,lambda_p,assm_xs,core_x_size,core_y_size,num_eg, &
                                  assm_map,xflux)
      ENDIF
      ! TODO implement inner tolerances
      DO g=1,num_eg
        CALL inner_solve(inner_solve_method, prob_size, 1d-1*MIN(tol_xflux,tol_xkeff), 10000, &
                        amat(:,:,g), bvec(:,g), xflux(:,:,g),core_x_size,core_y_size)
      ENDDO

      fiss_src_sum(2)=calc_fiss_src_sum(core_x_size,core_y_size,num_eg,assm_map,xflux,assm_xs)

      IF (iter > 1) THEN
        xkeff=xkeff*fiss_src_sum(2)/fiss_src_sum(1)
        ckeff=xkeff
        !if it's weilandt shift, update appropriately
        IF(dl_weilandt .GT. 0.0D0 .AND. iter > 10)THEN
          write(*,*)lambda_p,xkeff
          ckeff=(1.0D0/lambda_p+1.0D0/xkeff)**(-1)
          write(*,*)ckeff
          IF(ckeff .LE. 0)STOP 'negative ckeff'
        ENDIF
      ENDIF
      conv_xkeff=ABS(ckeff-ckeff_old) ! absolute

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
      ckeff_old=ckeff
      flux_old=xflux
      fiss_src_sum(1) = fiss_src_sum(2)

      CALL print_log(TRIM(str(iter,5))//'   '//TRIM(str(ckeff,6,'F'))//'   '//TRIM(str(conv_xkeff,2)) &
        //'   '//TRIM(str(conv_xflux,2)))
      !update the dtilde factors if we are using the poly expansion method
      IF(nodal_method .EQ. 'poly')THEN
        CALL comp_dtilde(core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux,dtilde_x, &
                        dtilde_y,h_x,h_y)
        CALL build_amatrix(amat,core_x_size,core_y_size,num_eg,assm_xs,assm_map,h_x,h_y,dtilde_x, &
                          dtilde_y)
        IF(dl_weilandt .GT. 0.0D0)THEN
          amat_base=amat
        ENDIF
      ENDIF

      IF ((conv_xflux < tol_xflux) .AND. (conv_xkeff < tol_xkeff)) EXIT

    ENDDO ! iter = 1,tol_max_iter

    xkeff=ckeff

    DEALLOCATE(amat, flux_old, bvec)

    CALL print_log('ITERATIONS FINISHED')
    CALL print_log('XKEFF = '//str(xkeff,16,'F'))

  ENDSUBROUTINE solver

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine builds the amatrix for the CMFD problem
!> @param amatrix - amatrix to build
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_xs - assembly level cross sections
!> @param assm_map - assembly map
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!> @param dtilde_x - current correction cross section used in the x direction
!> @param dtilde_y - current correction cross section used in the y direction
!>
  SUBROUTINE build_amatrix(amatrix,core_x_size,core_y_size,num_eg,assm_xs,assm_map,h_x,h_y,dtilde_x, &
                          dtilde_y)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    REAL(kr8), INTENT(IN) :: h_x(:),h_y(:),dtilde_x(:,:,:),dtilde_y(:,:,:)
    REAL(kr8), INTENT(OUT) :: amatrix(:,:,:)
    !loca variables
    INTEGER(ki4) :: g,j,i,cell_idx

    ! (1   , 2   , 3    , 4     , 5   )
    ! (diag, left, right, below, above)
    ! (diag, left, right, down , up   )
    ! (diag, west, east , south, north)

    amatrix=0.0D0
    DO j=1,core_y_size
      DO i=1,core_x_size
        cell_idx=calc_idx(i,j,core_x_size)
        DO g=1,num_eg
          !total xs term for the base of the A matrix diagonal
          amatrix(1,cell_idx,g)=assm_xs(assm_map(i,j))%sigma_r(g)
          !left term
          IF(i .NE. 1)THEN
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)&
              +2.0D0*(h_x(i)/assm_xs(assm_map(i,j))%D(g)&
              +h_x(i-1)/assm_xs(assm_map(i-1,j))%D(g))**(-1)/h_x(i)-dtilde_x(i,j,g)/h_x(i)
            amatrix(2,cell_idx,g)=amatrix(2,cell_idx,g)&
              -2.0D0*(h_x(i)/assm_xs(assm_map(i,j))%D(g)&
              +h_x(i-1)/assm_xs(assm_map(i-1,j))%D(g))**(-1)/h_x(i)-dtilde_x(i,j,g)/h_x(i)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)+dtilde_x(i,j,g)/h_x(i)
          ENDIF
          !right term
          IF(i .NE. core_x_size)THEN
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)&
              +2.0D0*(h_x(i)/assm_xs(assm_map(i,j))%D(g)&
              +h_x(i+1)/assm_xs(assm_map(i+1,j))%D(g))**(-1)/h_x(i)+dtilde_x(i+1,j,g)/h_x(i)
            amatrix(3,cell_idx,g)=amatrix(3,cell_idx,g)&
              -2.0D0*(h_x(i)/assm_xs(assm_map(i,j))%D(g)&
              +h_x(i+1)/assm_xs(assm_map(i+1,j))%D(g))**(-1)/h_x(i)+dtilde_x(i+1,j,g)/h_x(i)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)+dtilde_x(i+1,j,g)/h_x(i)
          ENDIF
          !below term
          IF(j .NE. 1)THEN
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)&
              +2.0D0*(h_y(j)/assm_xs(assm_map(i,j))%D(g)&
              +h_y(j-1)/assm_xs(assm_map(i,j-1))%D(g))**(-1)/h_y(j)-dtilde_y(i,j,g)/h_y(j)
            amatrix(4,cell_idx,g)=amatrix(4,cell_idx,g)&
              -2.0D0*(h_y(j)/assm_xs(assm_map(i,j))%D(g)&
              +h_y(j-1)/assm_xs(assm_map(i,j-1))%D(g))**(-1)/h_y(j)-dtilde_y(i,j,g)/h_y(j)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)+dtilde_y(i,j,g)/h_y(j)
          ENDIF
          !above term
          IF(j .NE. core_y_size)THEN
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)&
              +2.0D0*(h_y(j)/assm_xs(assm_map(i,j))%D(g)&
              +h_y(j+1)/assm_xs(assm_map(i,j+1))%D(g))**(-1)/h_y(j)+dtilde_y(i,j+1,g)/h_y(j)
            amatrix(5,cell_idx,g)=amatrix(5,cell_idx,g)&
              -2.0D0*(h_y(j)/assm_xs(assm_map(i,j))%D(g)&
              +h_y(j+1)/assm_xs(assm_map(i,j+1))%D(g))**(-1)/h_y(j)+dtilde_y(i,j+1,g)/h_y(j)
          ELSE
            !boundary condition specified by the dtilde factor computation in the polynomial portion
            amatrix(1,cell_idx,g)=amatrix(1,cell_idx,g)+dtilde_y(i,j+1,g)/h_y(j)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE build_amatrix

!---------------------------------------------------------------------------------------------------
!> @brief This converts a stripe represented amatrix to a dense amatrix
!> @param stripe - stripe matrix to convert
!> @param dense - dense matrix that is converted to
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!>
  SUBROUTINE stripe_to_dense(stripe, dense,core_x_size,core_y_size)
    REAL(kr8), INTENT(IN)  :: stripe(:,:)
    INTEGER(ki4) :: core_x_size,core_y_size
    REAL(kr8), INTENT(OUT) :: dense(:,:)

    INTEGER(ki4) :: i, j
    INTEGER(ki4) :: cell_idx, neigh_idx

    dense = 0d0

    DO j = 1,core_y_size
      DO i = 1,core_x_size

        cell_idx=calc_idx(i,j,core_x_size)

        ! self
        dense(cell_idx,cell_idx) = stripe(1,cell_idx)

        ! west
        IF (i .NE. 1) THEN
          neigh_idx=calc_idx(i-1,j,core_x_size)
          dense(cell_idx,neigh_idx) = stripe(2,cell_idx)
        ENDIF

        ! east
        IF (i .NE. core_x_size) THEN
          neigh_idx=calc_idx(i+1,j,core_x_size)
          dense(cell_idx,neigh_idx) = stripe(3,cell_idx)
        ENDIF

        ! south
        IF (j .NE. 1) THEN
          neigh_idx=calc_idx(i,j-1,core_x_size)
          dense(cell_idx,neigh_idx) = stripe(4,cell_idx)
        ENDIF

        ! north
        IF (j .NE. core_y_size) THEN
          neigh_idx=calc_idx(i,j+1,core_x_size)
          dense(cell_idx,neigh_idx) = stripe(5,cell_idx)
        ENDIF

      ENDDO ! i = 1,core_x_size
    ENDDO ! j = 1,core_y_size

    RETURN
  ENDSUBROUTINE stripe_to_dense

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine solves inner linear system of equations (allows switching solvers)
!> @param method - linear solver method to use
!> @param rank - system size
!> @param tol_inner_x - Inner tolerance
!> @param tol_inner_maxit - Number of max iterations
!> @param amat - A matrix
!> @param bvec - RHS b vector
!> @param flux - scalar flux
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!>
  SUBROUTINE inner_solve(method, rank, tol_inner_x, tol_inner_maxit, &
                         amat, bvec, flux,core_x_size,core_y_size)
    CHARACTER(*), INTENT(IN)    :: method
    INTEGER,      INTENT(IN)    :: rank,core_x_size,core_y_size
    REAL(kr8),    INTENT(IN)    :: tol_inner_x
    INTEGER,      INTENT(IN)    :: tol_inner_maxit
    REAL(kr8),    INTENT(IN)    :: amat(:,:)
    REAL(kr8),    INTENT(IN)    :: bvec(:)
    REAL(kr8),    INTENT(INOUT) :: flux(:,:)

    ! for dgesv
    REAL(kr8), ALLOCATABLE :: atemp(:,:)
    INTEGER(ki4), ALLOCATABLE :: ipiv(:)
    INTEGER(ki4) :: ierr
    INTEGER(ki4) :: cell_idx
    INTEGER(ki4) :: i, j

    SELECT CASE (method)
      CASE ('dgesv')
        ALLOCATE(atemp(rank,rank))
        ALLOCATE(ipiv(rank))
        call stripe_to_dense(amat, atemp,core_x_size,core_y_size)
        CALL DGESV(rank,1,atemp,rank,ipiv,bvec,rank,ierr)
        IF(ierr .NE. 0)CALL fatal_error('DGESV failed')
        !assign xflux to the new bvec
        DO j=1,core_y_size
          DO i=1,core_x_size
            cell_idx=calc_idx(i,j,core_x_size)
            flux(i,j)=MAX(bvec(cell_idx),0d0)
          ENDDO
        ENDDO
        DEALLOCATE(atemp, ipiv)
      CASE ('sor')
        ! TODO set omega better than this
        CALL sor(amat, bvec, flux, rank, 1.2d0, tol_inner_x, tol_inner_maxit,core_x_size,core_y_size)
      CASE DEFAULT
        call fatal_error('selected inner_solve method not implemented')
    ENDSELECT

    RETURN
  ENDSUBROUTINE inner_solve

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine builds the bvector for current flux (sans downscattering)
!> @param bvec - RHS b vector
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_xs - assembly level cross sections
!> @param xkeff - eigenvalue
!> @param xflux - scalar flux
!> @param assm_map - assembly map
!>
  SUBROUTINE build_bvec(bvec,core_x_size,core_y_size,num_eg,assm_xs,xkeff,xflux,assm_map)
    REAL(kr8), INTENT(OUT) :: bvec(:,:)
    REAL(kr8), INTENT(IN) :: xkeff,xflux(:,:,:)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: i,j,g,gp,cell_idx,loc_id
    REAL(kr8) :: fiss_src
    bvec=0.0D0
    DO j=1,core_y_size
      DO i=1,core_x_size
        cell_idx=calc_idx(i,j,core_x_size)
        !compute fission source for this cell
        loc_id=assm_map(i,j)
        fiss_src=0.0D0
        DO gp=1,num_eg
          fiss_src=fiss_src+xflux(i,j,gp)*assm_xs(loc_id)%nusigma_f(gp)
        ENDDO
        DO g=1,num_eg
          !set bvec to fission source
          bvec(cell_idx,g)=fiss_src*assm_xs(loc_id)%chi(g)/xkeff
          !and add up-scattering source
          DO gp=g+1,num_eg
            bvec(cell_idx,g)=bvec(cell_idx,g)+xflux(i,j,gp)*assm_xs(loc_id)%sigma_scat(g,gp)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE build_bvec

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine adds downscatter from the just solved group to the bvector
!> @param bvec - RHS b vector
!> @param gin - group we are scattering from to get the added downscattering
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_xs - assembly level cross sections
!> @param assm_map - assembly map
!> @param xflux - scalar flux
!>
  SUBROUTINE add_downscatter(bvec,gin,core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux)
    REAL(kr8), INTENT(INOUT) :: bvec(:,:)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:)
    INTEGER(ki4), INTENT(IN) :: gin,core_x_size,core_y_size,num_eg,assm_map(:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: i,j,g,cell_idx,loc_id
    DO j=1,core_y_size
      DO i=1,core_x_size
        cell_idx=calc_idx(i,j,core_x_size)
        loc_id=assm_map(i,j)
        DO g=gin+1,num_eg
          !add down-scattering source
          bvec(cell_idx,g)=bvec(cell_idx,g)+xflux(i,j,gin)*assm_xs(loc_id)%sigma_scat(g,gin)
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE add_downscatter

!---------------------------------------------------------------------------------------------------
!> @brief This function calculates the fission source
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_map - assembly map
!> @param xflux - scalar flux
!> @param assm_xs - assembly level cross sections
!>
  REAL(kr8) FUNCTION calc_fiss_src_sum(core_x_size,core_y_size,num_eg,assm_map,xflux,assm_xs)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: i,j,gp,loc_id
    calc_fiss_src_sum=0.0D0
    DO gp=1,num_eg
      DO j=1,core_y_size
        DO i=1,core_x_size
          loc_id=assm_map(i,j)
          calc_fiss_src_sum=calc_fiss_src_sum+xflux(i,j,gp)*assm_xs(loc_id)%nusigma_f(gp)
        ENDDO ! i = 1,core_x_size
      ENDDO ! j = 1,core_y_size
    ENDDO ! gp = 1,num_eg
  ENDFUNCTION calc_fiss_src_sum

  !shift the amatrix using the weilandt shift
  FUNCTION weilandt_shift_amatrix(amat_base,lambda_p,assm_xs,core_x_size,core_y_size,num_eg, &
                                  assm_map,xflux)
    REAL(kr8),INTENT(IN) :: amat_base(:,:,:),lambda_p,xflux(:,:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    !local variables
    REAL(kr8) :: weilandt_shift_amatrix(5,core_x_size*core_y_size,num_eg)
    INTEGER(ki4) :: g,i,j,cell_idx,loc_id
    REAL(kr8) :: fiss_src
    weilandt_shift_amatrix=amat_base
    DO i=1,core_x_size
      DO j=1,core_y_size
        loc_id=assm_map(i,j)
        cell_idx=calc_idx(i,j,core_x_size)
        fiss_src=0.0D0
        DO g=1,num_eg
          fiss_src=fiss_src+xflux(i,j,g)*assm_xs(loc_id)%nusigma_f(g)
        ENDDO
        DO g=1,num_eg
          weilandt_shift_amatrix(1,cell_idx,g)=weilandt_shift_amatrix(1,cell_idx,g)&
            -(fiss_src*assm_xs(loc_id)%chi(g))/(lambda_p*xflux(i,j,g))
        ENDDO
      ENDDO
    ENDDO
  ENDFUNCTION weilandt_shift_amatrix

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine calculates dtilde
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_xs - assembly level cross sections
!> @param assm_map - assembly map
!> @param xflux - scalar flux
!> @param dtilde_x - current correction cross section used in the x direction
!> @param dtilde_y - current correction cross section used in the y direction
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!>
  SUBROUTINE comp_dtilde(core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux,dtilde_x, &
                          dtilde_y,h_x,h_y)
    INTEGER, INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:),h_x(:),h_y(:)
    REAL(kr8), INTENT(INOUT) :: dtilde_x(:,:,:),dtilde_y(:,:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    REAL(kr8) :: s_bar_x(core_x_size,core_y_size,num_eg),s_bar_y(core_x_size,core_y_size,num_eg)
    CALL comp_s(s_bar_x,s_bar_y,core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux,dtilde_x, &
                dtilde_y,h_x,h_y)
    STOP 'comp_dtilde not yet complete'
  ENDSUBROUTINE comp_dtilde

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine the sbar values
!> @param s_bar_x - sbar in the x direction
!> @param s_bar_y - sbar in the y direction
!> @param core_x_size - core size in the x direction
!> @param core_y_size - core size in the y direction
!> @param num_eg - number of energy groups
!> @param assm_xs - assembly level cross sections
!> @param assm_map - assembly map
!> @param xflux - scalar flux
!> @param dtilde_x - current correction cross section used in the x direction
!> @param dtilde_y - current correction cross section used in the y direction
!> @param h_x - node widths in the x direction
!> @param h_y - node widths in the y direction
!>
  SUBROUTINE comp_s(s_bar_x,s_bar_y,core_x_size,core_y_size,num_eg,assm_xs,assm_map,xflux,dtilde_x, &
                    dtilde_y,h_x,h_y)
    INTEGER, INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:),h_x(:),h_y(:)
    REAL(kr8), INTENT(INOUT) :: dtilde_x(:,:,:),dtilde_y(:,:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    REAL(kr8), INTENT(OUT) :: s_bar_x(:,:,:),s_bar_y(:,:,:)
    REAL(kr8) :: j_x(core_x_size+1,core_y_size,num_eg),j_y(core_x_size,core_y_size+1,num_eg)
    REAL(kr8) :: l_bar_x(core_x_size,core_y_size,num_eg),l_bar_y(core_x_size,core_y_size,num_eg)
    INTEGER(ki4) :: i,j,g

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
            j_x(i,j,g)=2.0D0*(h_x(i-1)/assm_xs(assm_map(i-1,j))%D(g)&
              +h_x(i)/assm_xs(assm_map(i,j))%D(g))*(xflux(i-1,j,g)-xflux(i,j,g))&
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
            j_y(i,j,g)=2.0D0*(h_y(j-1)/assm_xs(assm_map(i,j-1))%D(g)&
              +h_y(j)/assm_xs(assm_map(i,j))%D(g))*(xflux(i,j-1,g)-xflux(i,j,g))&
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
          s_bar_x(i,j,g)=l_bar_y(i,j,g)/h_y(j)
        ENDDO
      ENDDO
    ENDDO
    !compute s_bar_y values
    DO i=1,core_x_size
      DO j=1,core_y_size
        DO g=1,num_eg
          s_bar_y(i,j,g)=l_bar_x(i,j,g)/h_x(i)
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
