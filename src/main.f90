!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Program driver. A Fortran based nodal solver
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM opennodal
  USE analytic_module
  USE errors_module
  USE globals
  USE input_module
  USE output_module
  USE solvers_module
  USE xs_types
  USE string_module
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  IMPLICIT NONE

  !> problem dimension
  INTEGER(ki4) :: prob_dim=2
  !> core size, assumed square
  INTEGER(ki4) :: core_x_size=0,core_y_size=0
  !> assembly pitch
  REAL(kr8) :: assm_pitch=0
  !> problem symmetry
  CHARACTER(100) :: prob_sym='full'
  !> assembly map
  INTEGER(ki4), ALLOCATABLE :: assm_map(:,:)
  !> node widths
  REAL(kr8), ALLOCATABLE :: h_x(:),h_y(:)
  !> boundary condition option
  CHARACTER(100) :: bc_opt='vacuum'
  !> albedo boundary conditions
  REAL(kr8), ALLOCATABLE :: albedos(:)
  !> number energy groups
  INTEGER(ki4) :: num_eg
  !> nsplit value, for decomposing nodes into nsplit number of sub-nodes
  INTEGER(ki4) :: nsplit=1
  !>  keff convergence tolerance
  REAL(kr8) :: tol_xkeff = 1d-6
  !>  flux convergence tolerance
  REAL(kr8) :: tol_xflux = 1d-5
  !>  maximum number of iterations
  INTEGER(ki4) :: tol_max_iter = 100
  !> nodal method option
  CHARACTER(100) :: nodal_method='fd'
  !> number of unique assemblies
  INTEGER(ki4) :: num_assm_reg
  !> assembly level cross sections for the problem
  TYPE(macro_assm_xs_type), ALLOCATABLE :: assm_xs(:)
  !> reflector material
  INTEGER(ki4) :: refl_mat=0
  !> axial buckling for 2D problems
  REAL(kr8) :: ax_buckle=0.0D0
  !>  eigenvalue
  REAL(kr8) :: xkeff=1d0 ! TODO implement an initial user guess input AND output (restart files)
  !>  scalar flux
  REAL(kr8), ALLOCATABLE :: xflux(:,:,:) ! (nx,ny,ng
  !>D-tilde correction factors for each surface
  REAL(kr8), ALLOCATABLE  :: dtilde_x(:,:,:),dtilde_y(:,:,:)
  !> Weilandt shift value
  REAL(kr8) :: dl_wielandt=-1.0D0
  !> Analytic solution title
  CHARACTER(16) :: anal_ref = 'NONE'

  !variables to open log file
  CHARACTER(1024) :: t_char
  CHARACTER(1024) :: words(100)
  INTEGER(ki4) :: nwords,i

  !standard output
  stdout_unit=OUTPUT_UNIT

  !read command line arguments (just input file for now
  CALL read_cmd_args()

  !open the log file
  t_char=base_in
  CALL parse(t_char,'.',words,nwords)
  t_char=TRIM(ADJUSTL(words(1)))
  DO i=2,nwords
    IF((words(i) .EQ. 'inp' .OR. words(i) .EQ. 'input' .OR. words(i) .EQ. 'in') .AND. i .EQ. nwords)EXIT
    t_char=TRIM(ADJUSTL(t_char))//'.'//TRIM(ADJUSTL(words(i)))
  ENDDO
  OPEN(UNIT=log_unit, FILE=TRIM(ADJUSTL(t_char))//'.log', STATUS="REPLACE", ACTION="WRITE")

  !print the header
  CALL print_log('***************************** OpenNodal - Version 1 *****************************')
  CALL print_log('*********************************************************************************')
  CALL print_log('*********************************************************************************')
  CALL print_log('*********************************************************************************')
  CALL print_log('**                ____                   _   __          __      __            **')
  CALL print_log('**               / __ \____  ___  ____  / | / /___  ____/ /___ _/ /            **')
  CALL print_log('**              / / / / __ \/ _ \/ __ \/  |/ / __ \/ __  / __ `/ /             **')
  CALL print_log('**             / /_/ / /_/ /  __/ / / / /|  / /_/ / /_/ / /_/ / /              **')
  CALL print_log('**             \____/ .___/\___/_/ /_/_/ |_/\____/\__,_/\__,_/_/               **')
  CALL print_log('**                 /_/                                                         **')
  CALL print_log('**                                                                             **')
  CALL print_log('**   By N.F. Herring and W.C. Dawn                                             **')
  CALL print_log('*********************************************************************************')
  CALL print_log('*********************************************************************************')
  CALL print_log('*********************************************************************************')
  CALL print_log('***************************** OpenNodal - Version 1 *****************************')

  !read the base input file data and cross sections
  CALL read_files(prob_dim,core_x_size,core_y_size,assm_pitch,prob_sym,assm_map,h_x,h_y,bc_opt, &
                  albedos,num_eg,nsplit,tol_xkeff,tol_xflux,tol_max_iter,nodal_method,num_assm_reg, &
                  assm_xs,refl_mat,ax_buckle,dl_wielandt,anal_ref)

  ! initialize memory and the solver
  ! this is done here to allow for multi-state runs in the future rather than
  ! initializing inside of the solver
  CALL solver_init(core_x_size,core_y_size,num_eg,assm_map,refl_mat,num_assm_reg,assm_xs,ax_buckle, &
                  h_x,h_y,nsplit,assm_pitch,xkeff,xflux,dtilde_x,dtilde_y,bc_opt,albedos,prob_sym)

  !call the solver
  CALL solver(core_x_size,core_y_size,num_eg,tol_xflux,tol_xkeff,xflux,xkeff,tol_max_iter, &
              nodal_method,assm_map,assm_xs,dtilde_x,dtilde_y,h_x,h_y,dl_wielandt)

  if (anal_ref /= 'NONE') call anal (xkeff, xflux, anal_ref, assm_map, assm_pitch, assm_xs, &
                                     core_x_size, core_y_size, num_eg)

  !output the results
  CALL output_results(xflux,xkeff,core_x_size,core_y_size,nsplit,num_eg,assm_xs,h_x,h_y,assm_map)

ENDPROGRAM opennodal
