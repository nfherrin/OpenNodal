!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for outputting results.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE output_module
  USE globals
  USE errors_module
  USE string_module
  USE xs_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: output_results

  INTEGER(ki4),PARAMETER :: out_unit=31
CONTAINS

  !output results
  SUBROUTINE output_results(xflux,xkeff,core_x_size,core_y_size,nsplit,num_eg,assm_xs,h_x,h_y,assm_map)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,nsplit,num_eg,assm_map(:,:)
    REAL(kr8), INTENT(IN) :: xflux(:,:,:),xkeff,h_x(:),h_y(:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: t_int,g,j,i,nwords
    CHARACTER(100) :: t_char,words(100)

    CALL parse(base_in,'.',words,nwords)
    prob_title=TRIM(ADJUSTL(words(1)))
    DO i=2,nwords
      IF(words(i) .EQ. 'inp')EXIT
      prob_title=TRIM(ADJUSTL(prob_title))//'.'//TRIM(ADJUSTL(words(i)))
    ENDDO

    OPEN(UNIT=out_unit, FILE=TRIM(prob_title)//'_results.out', STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
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

    CALL create_flux_csv(core_x_size,core_y_size,num_eg,xflux,h_x,h_y,assm_xs,assm_map)

    !only plot if gnuplot is found
    IF(check_for('gnuplot'))CALL plot_flux(num_eg,h_x,h_y)
  ENDSUBROUTINE output_results

  !creates the flux csv output
  SUBROUTINE create_flux_csv(core_x_size,core_y_size,num_eg,xflux,h_x,h_y,assm_xs,assm_map)
    INTEGER(ki4), INTENT(IN) :: core_x_size,core_y_size,num_eg,assm_map(:,:)
    REAL(kr8), INTENT(IN) :: h_x(:),h_y(:),xflux(:,:,:)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    !local variables
    INTEGER(ki4) :: t_int,g,j,i
    CHARACTER(100) :: t_char
    REAL(kr8),ALLOCATABLE :: fiss_src(:,:)

    ALLOCATE(fiss_src(core_x_size,core_y_size))
    fiss_src=0.0D0

    DO g=1,num_eg
      WRITE(t_char,'(A,I0,A)')TRIM(prob_title)//'_flux_g',g,'.csv'
      OPEN(UNIT=out_unit, FILE=t_char, STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
      IF(t_int .NE. 0)THEN
        CALL fatal_error(t_char)
      ENDIF

      !print out the CSV flux data
      DO i=1,core_x_size
        DO j=1,core_y_size
          WRITE(out_unit,'(ES16.8,A,ES16.8,A,ES16.8)')SUM(h_x(1:i))-0.5D0*h_x(i),', '&
            ,SUM(h_y(1:j))-0.5D0*h_y(j),', ',xflux(i,j,g)
          fiss_src(i,j)=fiss_src(i,j)+xflux(i,j,g)*assm_xs(assm_map(i,j))%nusigma_f(g)/assm_xs(assm_map(i,j))%nu(g)
        ENDDO
        WRITE(out_unit,*)
      ENDDO

      CLOSE(out_unit)
    ENDDO

    WRITE(t_char,'(A)')TRIM(prob_title)//'_fiss_src.csv'
    OPEN(UNIT=out_unit, FILE=t_char, STATUS='REPLACE', ACTION = "WRITE", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      CALL fatal_error(t_char)
    ENDIF

    !print out the CSV flux data
    DO i=1,core_x_size
      DO j=1,core_y_size
        WRITE(out_unit,'(ES16.8,A,ES16.8,A,ES16.8)')SUM(h_x(1:i))-0.5D0*h_x(i),', '&
          ,SUM(h_y(1:j))-0.5D0*h_y(j),', ',fiss_src(i,j)
      ENDDO
      WRITE(out_unit,*)
    ENDDO

    CLOSE(out_unit)

    DEALLOCATE(fiss_src)

  ENDSUBROUTINE create_flux_csv

  !plot the flux in the CSV files using gnuplot, like a real engineer...
  SUBROUTINE plot_flux(num_eg,h_x,h_y)
    INTEGER(ki4), INTENT(IN) :: num_eg
    REAL(kr8), INTENT(IN) :: h_x(:),h_y(:)
    !local variables
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
      WRITE(out_unit_temp,'(A,I0,A)')'set output "'//TRIM(prob_title)//'_flux_g',g,'.png"'
      WRITE(out_unit_temp,'(A,I0,A)')'set title "Flux Group = ',g,'"'
      WRITE(out_unit_temp,'(A)')'set grid'
      WRITE(out_unit_temp,'(A)')'set xlabel "x [cm]"'
      WRITE(out_unit_temp,'(A)')'set ylabel "y [cm]"'
      WRITE(out_unit_temp,'(A,F16.8)')'set size ratio ',SUM(h_y(:))/SUM(h_x(:))
      WRITE(out_unit_temp,'(A,ES16.8,A)')'set xrange [0:',SUM(h_x(:)),']'
      WRITE(out_unit_temp,'(A,ES16.8,A)')'set yrange [',SUM(h_y(:)),':0]'
      WRITE(out_unit_temp,'(A,I0,A)')'plot "'//TRIM(prob_title)//'_flux_g',g,'.csv" with image'

      CALL EXECUTE_COMMAND_LINE('gnuplot -c temp.plotcommands.temp')

      CLOSE(out_unit_temp)
    ENDDO

    OPEN(UNIT=out_unit_temp, FILE='temp.plotcommands.temp', STATUS='REPLACE', ACTION = "WRITE", &
        IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      CALL fatal_error(t_char)
    ENDIF

    !output the plot commands
    WRITE(out_unit_temp,'(A)')'# plot.plt'
    WRITE(out_unit_temp,'(A)')'set term png'
    WRITE(out_unit_temp,'(A)')'set output "'//TRIM(prob_title)//'_fiss_src.png"'
    WRITE(out_unit_temp,'(A)')'set title "Fission Source"'
    WRITE(out_unit_temp,'(A)')'set grid'
    WRITE(out_unit_temp,'(A)')'set xlabel "x [cm]"'
    WRITE(out_unit_temp,'(A)')'set ylabel "y [cm]"'
    WRITE(out_unit_temp,'(A,F16.8)')'set size ratio ',SUM(h_y(:))/SUM(h_x(:))
    WRITE(out_unit_temp,'(A,ES16.8,A)')'set xrange [0:',SUM(h_x(:)),']'
    WRITE(out_unit_temp,'(A,ES16.8,A)')'set yrange [',SUM(h_y(:)),':0]'
    WRITE(out_unit_temp,'(A,A)')'plot "'//TRIM(prob_title)//'_fiss_src.csv" with image'

    CALL EXECUTE_COMMAND_LINE('gnuplot -c temp.plotcommands.temp')

    CLOSE(out_unit_temp)

    CALL EXECUTE_COMMAND_LINE('rm temp.plotcommands.temp')
  ENDSUBROUTINE plot_flux
ENDMODULE output_module
