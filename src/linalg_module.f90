!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for general linear algebra routines.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE linalg_module
  USE precisions
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: norm

CONTAINS

!---------------------------------------------------------------------------------------------------
!> @brief This subroutine takes the norm of a vector
!> @param vec - vector to take the norm of
!> @param ell_in - norm order to take
!>
  REAL(kr8) PURE FUNCTION norm(vec, ell_in)
    REAL(kr8), INTENT(IN) :: vec(:)
    INTEGER(ki4), INTENT(IN), OPTIONAL :: ell_in
    !local variables
    INTEGER(ki4) :: ell, i

    IF (PRESENT(ell_in)) THEN
      ell = ell_in
    ELSE
      ell = 2 ! default to l2 norm
    ENDIF

    SELECT CASE (ell)
      CASE (-1) ! infinite/maximal norm
        norm = MAXVAL(ABS(vec))
      CASE (1)
        norm = SUM(ABS(vec))
      CASE (2)
        norm = SQRT(SUM(vec**2))
      CASE DEFAULT
        norm = 0d0
        DO i = 1,SIZE(vec)
          norm = norm + ABS(vec(i))**ell
        ENDDO ! j = 1,size(vec)
        norm = norm**(1d0/DBLE(ell))
    ENDSELECT ! ell

    RETURN
  ENDFUNCTION norm

ENDMODULE linalg_module
