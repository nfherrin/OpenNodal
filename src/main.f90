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

    !read the base input file
    CALL read_files()

    ! initialize memory
    ! this is done here to allow for multi-state runs in the future rather than
    ! initializing inside of the solver

    CALL solver_init()

    CALL solver()

    CALL output_results()

    deallocate(xflux)

ENDPROGRAM opennodal
