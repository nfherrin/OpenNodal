!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for globals. Both variables and routines
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE globals
  USE precisions
  USE xs_types
  IMPLICIT NONE

  !> input filename, assume no longer than 100 characters
  CHARACTER(100) :: base_in=''

  !> problem title
  CHARACTER(100) :: prob_title=''

  !> problem dimensions
  INTEGER(ki4) :: prob_dim=0

  !> core size, assumed square
  INTEGER(ki4) :: core_size=0,core_x_size=0,core_y_size=0

  !> assembly pitch
  REAL(kr8) :: assm_pitch=0

  !> assembly map
  INTEGER(ki4),ALLOCATABLE :: assm_map(:,:)

  !> number of unique assemblies
  INTEGER(ki4) :: num_assm_reg

  !> assembly level cross sections for the problem
  TYPE(macro_assm_xs_type),ALLOCATABLE :: assm_xs(:)

  !>D-tilde correction factors for each surface
  REAL(kr8),ALLOCATABLE  :: dtilde_x(:,:,:),dtilde_y(:,:,:)

  !> number energy groups
  INTEGER(ki4) :: num_eg

  !> problem symmetry
  CHARACTER(100) :: prob_sym='full'

  !> standard output unit
  INTEGER(ki4) :: stdout_unit=0

  !> log input unit
  INTEGER(ki4), PARAMETER :: log_unit=9999

  !>  eigenvalue
  REAL(kr8) :: xkeff=1d0 ! TODO implement an initial user guess
  !>  scalar flux
  REAL(kr8), ALLOCATABLE :: xflux(:,:,:) ! (nx,ny)

  !>  maximum number of iterations
  INTEGER(ki4) :: tol_max_iter = 100
  !>  keff convergence tolerance
  REAL(kr8) :: tol_xkeff = 1d-6
  !>  flux convergence tolerance
  REAL(kr8) :: tol_xflux = 1d-5

  !> nsplit value, for decomposing nodes into nsplit number of sub-nodes
  INTEGER(ki4) :: nsplit=1

  !> nodal method option
  CHARACTER(100) :: nodal_method='fd'

  !> boundary condition option
  CHARACTER(100) :: bc_opt='vacuum'

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine prints to both the log file and the screen
!> @param log_msg - string to print to log and screen
!> @param advancing - logical to indicate if the printed message should advance the line
!>
  SUBROUTINE print_log(log_msg,advancing)
    CHARACTER(*),INTENT(IN) :: log_msg
    LOGICAL, OPTIONAL, INTENT(IN) :: advancing
    LOGICAL :: advopt

    advopt=.TRUE.
    IF(PRESENT(advancing))advopt=advancing

    IF(advopt)THEN
      WRITE(stdout_unit,'(A)')TRIM(log_msg)
      WRITE(log_unit,'(A)')TRIM(log_msg)
    ELSE
      WRITE(stdout_unit,'(A)',ADVANCE='NO')TRIM(log_msg)
      WRITE(log_unit,'(A)',ADVANCE='NO')TRIM(log_msg)
    ENDIF
  ENDSUBROUTINE print_log
ENDMODULE globals
