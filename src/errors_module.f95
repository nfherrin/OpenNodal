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
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")
    IF (PRESENT(message)) THEN
      CALL printlog('FATAL ERROR!')
      CALL printlog('ERROR: '//TRIM(ADJUSTL(message)))
    ELSE
      CALL printlog('FATAL ERROR!')
      CALL printlog('...')
      CALL printlog('No error message given')
    ENDIF
    CALL printlog('>> OpenNodal encountered a fatal error!')
    CALL printlog('>> Execution of OpenNodal terminated UNsuccessfully!')
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")

    STOP
  ENDSUBROUTINE fatal_error
ENDMODULE errors_module
