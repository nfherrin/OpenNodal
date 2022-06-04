!A Fortran based neural network accelerated simulated annealing software.
PROGRAM fortNNASA
  USE globals
  USE input_module
  USE output_module
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

    !standard output
    stdout_unit=OUTPUT_UNIT

    !read command line arguments (just input file for now
    CALL read_cmd_args()

    OPEN(UNIT=log_unit, FILE=TRIM(ADJUSTL(base_in))//'.log', STATUS="REPLACE", ACTION="WRITE")

    CALL printlog('OpenNodal - Version 1')

    !read the base input file
    CALL read_base_input()

ENDPROGRAM fortNNASA
