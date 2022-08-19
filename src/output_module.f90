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

  !output results
  SUBROUTINE output_results()
    INTEGER(ki4) :: t_int,g,j,i
    CHARACTER(100) :: t_char

    OPEN(UNIT=out_unit, FILE=TRIM(base_in)//'_results.out', STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      CALL fatal_error(t_char)
    ENDIF

    WRITE(out_unit,'(A,F24.16)')'Final Eigenvalue = ',xkeff
    IF(num_eg .EQ. 2)WRITE(out_unit,'(A,F21.16)')'Final 2G Flux Ratio = ',SUM(xflux(:,:,2))/SUM(xflux(:,:,1))

    DO g=1,num_eg
      WRITE(out_unit,'(A,I0)')'Node Averaged Flux G = ',g
      DO j=1,core_y_size,nsplit
        DO i=1,core_x_size,nsplit
          WRITE(out_unit,'(ES24.16)',ADVANCE='NO')SUM(xflux(i:i+nsplit-1,j:j+nsplit-1,g))/(1.0D0*nsplit**2)
        ENDDO
        WRITE(out_unit,*)
      ENDDO
    ENDDO

    CLOSE(out_unit)

    CALL create_flux_csv()

    !only plot if gnuplot is found
    IF(check_for('gnuplot'))CALL plot_flux()
  ENDSUBROUTINE output_results

  !creates the flux csv output
  SUBROUTINE create_flux_csv()
    INTEGER(ki4) :: t_int,g,j,i
    CHARACTER(100) :: t_char

    DO g=1,num_eg
      WRITE(t_char,'(A,I0,A)')TRIM(base_in)//'_flux_g',g,'.csv'
      OPEN(UNIT=out_unit, FILE=t_char, STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
      IF(t_int .NE. 0)THEN
        CALL fatal_error(t_char)
      ENDIF

      !print out the CSV flux data
      !this commented stuff can help avoid blank regions of the plot
      ! WRITE(out_unit,'(3ES16.8)')0.0D0,0.0D0,xflux(1,1,g)
      ! DO j=1,core_y_size
        ! WRITE(out_unit,'(3ES16.8)')0.0D0,SUM(h_y(1:j))-0.5D0*h_y(j),xflux(1,j,g)
      ! ENDDO
      ! WRITE(out_unit,'(3ES16.8)')0.0D0,SUM(h_y),xflux(1,core_y_size,g)
      ! WRITE(out_unit,*)
      DO i=1,core_x_size
        ! WRITE(out_unit,'(3ES16.8)')SUM(h_x(1:i))-0.5D0*h_x(i),0.0D0,xflux(i,1,g)
        DO j=1,core_y_size
          WRITE(out_unit,'(3ES16.8)')SUM(h_x(1:i))-0.5D0*h_x(i),SUM(h_y(1:j))-0.5D0*h_y(j),xflux(i,j,g)
        ENDDO
        ! WRITE(out_unit,'(3ES16.8)')SUM(h_x(1:i))-0.5D0*h_x(i),SUM(h_y),xflux(i,core_y_size,g)
        WRITE(out_unit,*)
      ENDDO
      ! WRITE(out_unit,'(3ES16.8)')SUM(h_x),0.0D0,xflux(core_x_size,1,g)
      ! DO j=1,core_y_size
        ! WRITE(out_unit,'(3ES16.8)')SUM(h_x),SUM(h_y(1:j))-0.5D0*h_y(j),xflux(core_x_size,j,g)
      ! ENDDO
      ! WRITE(out_unit,'(3ES16.8)')SUM(h_x),SUM(h_y),xflux(core_x_size,core_y_size,g)
      ! WRITE(out_unit,*)

      CLOSE(out_unit)
    ENDDO

  ENDSUBROUTINE create_flux_csv

  !plot the flux in the CSV files using gnuplot, like a real engineer...
  SUBROUTINE plot_flux()
    INTEGER(ki4) :: t_int,g
    CHARACTER(100) :: t_char
    INTEGER(ki4) :: out_unit_temp=41
    DO g=1,num_eg
      OPEN(UNIT=out_unit_temp, FILE='temp.plotcommands.temp', STATUS='REPLACE', ACTION = "WRITE", &
          IOSTAT=t_int, IOMSG=t_char)
      IF(t_int .NE. 0)THEN
        CALL fatal_error(t_char)
      ENDIF

      !output the plot commands
      WRITE(out_unit_temp,'(A)')'# plot.plt'
      WRITE(out_unit_temp,'(A)')'set term png'
      WRITE(out_unit_temp,'(A,I0,A)')'set output "flux_g',g,'.png"'
      WRITE(out_unit_temp,'(A,I0,A)')'set title "Flux Group = ',g,'"'
      WRITE(out_unit_temp,'(A)')'set grid'
      WRITE(out_unit_temp,'(A)')'set xlabel "x [cm]"'
      WRITE(out_unit_temp,'(A)')'set ylabel "y [cm]"'
      WRITE(out_unit_temp,'(A,F16.8)')'set size ratio ',SUM(h_y(:))/SUM(h_x(:))
      WRITE(out_unit_temp,'(A,ES16.8,A)')'set xrange [0:',SUM(h_x(:)),']'
      WRITE(out_unit_temp,'(A,ES16.8,A)')'set yrange [',SUM(h_y(:)),':0]'
      WRITE(out_unit_temp,'(A,I0,A)')'plot "'//TRIM(base_in)//'_flux_g',g,'.csv" with image'

      CALL EXECUTE_COMMAND_LINE('gnuplot -c temp.plotcommands.temp')

      CLOSE(out_unit_temp)
    ENDDO

    CALL EXECUTE_COMMAND_LINE('rm temp.plotcommands.temp')
  ENDSUBROUTINE plot_flux
ENDMODULE output_module
