!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for performing edits
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE edits_module
  USE globals
  USE xs_types
  USE string_module
IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine edits the cross sections, calculating the kinf for each assembly
!> @param assm_xs - assembly level cross sections
!> @param num_assm_reg - number of unique assemblies (or at least, unique assembly IDs)
!> @param num_eg - number of energy groups
!>
  SUBROUTINE edit_xs(assm_xs,num_assm_reg,num_eg)
    TYPE(macro_assm_xs_type), INTENT(IN) :: assm_xs(:)
    INTEGER, INTENT(IN) :: num_assm_reg,num_eg
    !local variables
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
