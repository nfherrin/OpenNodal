!cross section types
MODULE xs_types
  USE precisions

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: macro_assm_xs_type

  !cross section material type
  TYPE :: macro_assm_xs_type
    !material name
    CHARACTER(64) :: mat_id
    !indicator if material is fissile
    LOGICAL :: fissile=.TRUE.
    !Diffusion coefficient
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: D
    !fission spectrum
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: chi
    !fission production cross section nuSigma_f
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: nusigma_f
    !fission production nu
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: nu
    !absorption cross section Sigma_a
    REAL(kr8), DIMENSION(:), ALLOCATABLE :: sigma_a
    !scattering matrix
    REAL(kr8), DIMENSION(:,:), ALLOCATABLE :: sigma_scat
  ENDTYPE macro_assm_xs_type

END MODULE xs_types