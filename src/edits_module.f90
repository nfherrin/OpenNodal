module edits_module
IMPLICIT NONE

contains

  subroutine edit_xs()
    USE string_module, ONLY : str
    use xs_types, only : macro_assm_xs_type
    use precisions, only : ki4, kr8
    use globals, only : stdout_unit, & ! TODO output to log as well
      num_assm_reg, assm_xs, num_eg, print_log
    use errors_module, only : fatal_error
    IMPLICIT NONE

    integer(ki4) :: i
    real(kr8) :: flux_ratio, kinf

    if (num_eg /= 2) call fatal_error('only supporting 2 energy groups')

    do i = 1,num_assm_reg
      if (assm_xs(i)%fissile) then
        flux_ratio = assm_xs(i)%sigma_scat(2,1)/assm_xs(i)%sigma_a(2)
        kinf = (assm_xs(i)%nusigma_f(1) + assm_xs(i)%nusigma_f(2)*flux_ratio) / &
          (assm_xs(i)%sigma_a(1) + assm_xs(i)%sigma_scat(2,1))
        CALL print_log('kinf '//trim(adjustl(assm_xs(i)%mat_id))//' '//str(kinf,6,'F'))
      endif
    enddo

  endsubroutine edit_xs

endmodule edits_module
