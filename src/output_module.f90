!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for outputting results.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE output_module
  USE globals
  USE errors_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: output_results

  INTEGER(ki4),PARAMETER :: out_unit=31
CONTAINS

  SUBROUTINE output_results()
    INTEGER(ki4) :: t_int,g,j
    CHARACTER(100) :: t_char

    OPEN(UNIT=out_unit, FILE=TRIM(base_in)//'_results.out', STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      CALL fatal_error(t_char)
    ENDIF

    WRITE(out_unit,'(A,F24.16)')'Final Eigenvalue = ',xkeff

    DO g=1,num_eg
      WRITE(out_unit,'(A,I0)')'Flux G = ',g
      DO j=1,core_y_size
        WRITE(out_unit,'(10000ES24.16)')xflux(:,j,g)
      ENDDO
    ENDDO
  ENDSUBROUTINE output_results
ENDMODULE output_module
