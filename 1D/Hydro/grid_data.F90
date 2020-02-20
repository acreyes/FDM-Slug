module grid_data
  implicit none
  real, allocatable, dimension(:), save :: gr_xCoord
  real, save :: gr_xbeg,gr_xend,gr_dx,gr_K
  integer, save :: gr_i0,gr_ibeg,gr_iend,gr_imax,gr_ngc,gr_nx,gr_radius

  real, allocatable, dimension(:,:) :: gr_U ! conservative vars
  real, allocatable, dimension(:,:) :: gr_V ! primitive vars
  real, allocatable, dimension(:,:) :: gr_W ! characteristic vars

  real, allocatable, dimension(:,:) :: gr_vR       ! Plus char fluxes
  real, allocatable, dimension(:,:) :: gr_vL       ! Minus char fluxes
  real, allocatable, dimension(:,:) :: gr_flux     ! fluxes at interfaces
!  real, allocatable, dimension(:,:) :: gr_CntrFlux !fluxes at cell centers
  real, allocatable, dimension(:,:) :: gr_vP       ! Plus char fluxes
  real, allocatable, dimension(:,:) :: gr_vM       ! Minus char fluxes

  real, allocatable, dimension(:,:)   :: gr_eigval ! eigenvalues
  real, allocatable, dimension(:,:,:) :: gr_reigvc ! right eigenvectors
  real, allocatable, dimension(:,:,:) :: gr_leigvc ! left  eigenvectors
  real, allocatable, dimension(:)     :: gr_maxAlphas ! max char speeds throughout domain

  !GP vars
  real, allocatable, dimension(:  ) :: gr_GPv
  real, allocatable, dimension(:,:) :: gr_GPZ

end module grid_data
