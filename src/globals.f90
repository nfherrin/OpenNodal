!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for globals. Both variables and routines
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE globals
  USE precisions
  IMPLICIT NONE

  !> input filename, assume no longer than 100 characters
  CHARACTER(100) :: base_in=''
  !> problem title
  CHARACTER(100) :: prob_title=''

  !> standard output unit
  INTEGER(ki4) :: stdout_unit=0

  !> log input unit
  INTEGER(ki4), PARAMETER :: log_unit=9999

  !> it's pi
  REAL(kr8),PARAMETER :: pi=4*ATAN(1.d0)

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine prints to both the log file and the screen
!> @param log_msg - string to print to log and screen
!> @param advancing - logical to indicate if the printed message should advance the line
!>
  SUBROUTINE print_log(log_msg,advancing)
    CHARACTER(*), INTENT(IN) :: log_msg
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

!---------------------------------------------------------------------------------------------------
!> @brief Checks if software is installed
!> @param software_name - Name of software to check for
!>
  LOGICAL FUNCTION check_for(software_name)
    CHARACTER(*), INTENT(IN) :: software_name
    INTEGER(ki4) :: iexit
    CALL EXECUTE_COMMAND_LINE('which '//TRIM(software_name)//' 2> /dev/null > /dev/null', exitstat=iexit)
    check_for = (iexit == 0)
  ENDFUNCTION
ENDMODULE globals
