!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Program driver. A Fortran based nodal solver
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM opennodal
  use errors_module
  USE globals
  USE input_module
  USE output_module
  USE solvers_module
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  IMPLICIT NONE

    !standard output
    stdout_unit=OUTPUT_UNIT

    !read command line arguments (just input file for now
    CALL read_cmd_args()

    OPEN(UNIT=log_unit, FILE=TRIM(ADJUSTL(base_in))//'.log', STATUS="REPLACE", ACTION="WRITE")

    CALL print_log('*************************** OpenNodal - Version 1 ***************************')

    !read the base input file
    CALL read_files()

    ! initialize memory
    ! this is done here to allow for multi-state runs in the future rather than
    ! initializing inside of the solver
    if (num_eg /= 2) call fatal_error('only supporting 2 energy groups')

    CALL solver_init()

    CALL solver()

    CALL output_results()

    deallocate(xflux)

ENDPROGRAM opennodal
