!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module for precision parameters.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE precisions

  ! Real kinds

  !> single precision real
  INTEGER, PARAMETER :: kr4 = SELECTED_REAL_KIND(6,37)
  !>  precision real
  INTEGER, PARAMETER :: kr8 = SELECTED_REAL_KIND(15,307)

  ! Integer kinds

  !> single precision integer
  INTEGER, PARAMETER :: ki4 = SELECTED_INT_KIND(9)
  !> double precision real
  INTEGER, PARAMETER :: ki8 = SELECTED_INT_KIND(18)

  !Complex kinds

  !> single precision complex
  INTEGER, PARAMETER :: kc4 = kr4
  !> double precision complex
  INTEGER, PARAMETER :: kc8 = kr8

END MODULE precisions
