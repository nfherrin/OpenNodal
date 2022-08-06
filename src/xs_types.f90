!OpenNodal is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module defines the cross section types for the problem.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE xs_types
  USE precisions

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: macro_assm_xs_type

  !> This data type stores cross section material information
  TYPE :: macro_assm_xs_type
    !> material name
    CHARACTER(64) :: mat_id
    !> indicator if material is fissile
    LOGICAL :: fissile=.TRUE.
    !> Diffusion coefficient
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: D
    !> fission spectrum
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: chi
    !> fission production cross section nuSigma_f
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: nusigma_f
    !> fission production nu
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: nu
    !> absorption cross section Sigma_a
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: sigma_a
    !> total cross section Sigma_t
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: sigma_t
    !> scattering matrix
    REAL(kr8), DIMENSION(:,:), ALLOCATABLE :: sigma_scat
  ENDTYPE macro_assm_xs_type

END MODULE xs_types