!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for performing edits
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE edits_module
IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine edits the cross sections, calculating the kinf for each assembly
!>
  SUBROUTINE edit_xs()
    USE string_module, ONLY : str
    USE xs_types, ONLY : macro_assm_xs_type
    USE precisions, ONLY : ki4, kr8
    USE globals, ONLY : stdout_unit, num_assm_reg, assm_xs, num_eg, print_log
    USE errors_module, ONLY : fatal_error
    IMPLICIT NONE

    INTEGER(ki4) :: i
    REAL(kr8) :: flux_ratio, kinf

    DO i = 1,num_assm_reg
      IF (assm_xs(i)%fissile .AND. num_eg == 2) THEN
        flux_ratio = assm_xs(i)%sigma_scat(2,1)/assm_xs(i)%sigma_a(2)
        kinf = (assm_xs(i)%nusigma_f(1) + assm_xs(i)%nusigma_f(2)*flux_ratio) / &
          (assm_xs(i)%sigma_a(1) + assm_xs(i)%sigma_scat(2,1))
        CALL print_log('kinf '//TRIM(ADJUSTL(assm_xs(i)%mat_id))//' '//str(kinf,6,'F'))
      ENDIF
    ENDDO

  ENDSUBROUTINE edit_xs

ENDMODULE edits_module
