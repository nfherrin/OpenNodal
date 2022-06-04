!globals module. Both variables and variables
MODULE globals
  IMPLICIT NONE

  !input filename, assume no longer than 100 characters
  CHARACTER(100) :: base_in

  !standard output unit
  INTEGER :: stdout_unit

  !log input unit
  INTEGER, PARAMETER :: log_unit=9999
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
