!globals module. Both variables and variables
MODULE globals
  USE precisions
  USE xs_types
  IMPLICIT NONE

  !input filename, assume no longer than 100 characters
  CHARACTER(100) :: base_in=''

  !problem title
  CHARACTER(100) :: prob_title=''

  !problem dimensions
  INTEGER(ki4) :: prob_dim=0

  !core size, assumed square
  INTEGER(ki4) :: core_size=0,core_x_size=0,core_y_size=0

  !assembly pitch
  REAL(kr8) :: assm_pitch=0

  !assembly map
  INTEGER(ki4),ALLOCATABLE :: assm_map(:,:)

  !number of unique assemblies
  INTEGER(ki4) :: num_assm_reg

  !assembly level cross sections for the problem
  TYPE(macro_assm_xs_type),ALLOCATABLE :: assm_xs(:)

  !number energy groups
  INTEGER(ki4) :: num_eg

  !problem symmetry
  CHARACTER(100) :: prob_sym='full'

  !standard output unit
  INTEGER(ki4) :: stdout_unit=0

  !log input unit
  INTEGER(ki4), PARAMETER :: log_unit=9999
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE printlog(log_msg,advancing)
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
  ENDSUBROUTINE printlog
ENDMODULE globals
