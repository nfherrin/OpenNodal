!input functions
MODULE errors_module
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fatal_error
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fatal_error(message)
    CHARACTER(*), OPTIONAL,INTENT(IN) :: message

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
ENDMODULE errors_module
