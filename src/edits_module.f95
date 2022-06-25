module edits_module
IMPLICIT NONE

contains

  subroutine edit_xs()
    use xs_types, only : macro_assm_xs_type
    use precisions, only : ki4, kr8
    use globals, only : stdout_unit, & ! TODO output to log as well
      num_assm_reg, assm_xs, num_eg
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
        write(stdout_unit, '(a, x, a16, x, f8.6)') & ! TODO name capped at 16
          'kinf', trim(adjustl(assm_xs(i)%mat_id)), kinf
      endif
    enddo

  endsubroutine edit_xs

endmodule edits_module
