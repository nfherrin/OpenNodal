!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module error and warning reporting.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE errors_module
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fatal_error, raise_warning
CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine raises a fatal error and kills the program
!> @param message - error message to print
!>
  SUBROUTINE fatal_error(message)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: message

    ! Print message
    CALL print_log("*****************************************************************************")
    CALL print_log("*****************************************************************************")
    CALL print_log("*****************************************************************************")
    IF (PRESENT(message)) THEN
      CALL print_log('FATAL ERROR!')
      CALL print_log('ERROR: '//TRIM(ADJUSTL(message)))
    ELSE
      CALL print_log('FATAL ERROR!')
      CALL print_log('...')
      CALL print_log('No error message given')
    ENDIF
    CALL print_log('>> OpenNodal encountered a fatal error!')
    CALL print_log('>> Execution of OpenNodal terminated UNsuccessfully!')
    CALL print_log("*****************************************************************************")
    CALL print_log("*****************************************************************************")
    CALL print_log("*****************************************************************************")

    STOP
  ENDSUBROUTINE fatal_error

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine raises a warning before continuing the program
!> @param message - warning message to print
!>
  SUBROUTINE raise_warning(message)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: message

    IF (PRESENT(message)) THEN
      CALL print_log(' *** WARNING: '//TRIM(ADJUSTL(message)))
    ELSE
      CALL print_log(' *** WARNING: No warning message given')
    ENDIF
  ENDSUBROUTINE raise_warning
ENDMODULE errors_module
