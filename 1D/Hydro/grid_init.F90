subroutine grid_init()

#include "definition.h"

  use grid_data
  use read_initFile

  implicit none

!!$  call read_initFileInt ('slug.init','gr_nx',   gr_nx)
!!$  call read_initFileInt ('slug.init','gr_radius', gr_radius)
!!$  call read_initFileInt ('slug.init','gr_ngc',  gr_ngc)
!!$  call read_initFileReal('slug.init','gr_xbeg', gr_xbeg)
!!$  call read_initFileReal('slug.init','gr_xend', gr_xend)
!!$  call read_initFileReal('slug.init','gr_K', gr_K)

  ! allocate cell coordinates
  allocate(gr_xCoord(gr_nx+2*gr_ngc)); gr_xCoord = 0.0

  ! the first and the last interior cell index
  gr_ibeg = gr_ngc + 1
  gr_iend = gr_ngc + gr_nx

  ! the min and max of the entire cell index
  gr_i0   = 1
  gr_imax = 2*gr_ngc + gr_nx

  ! grid delta
  gr_dx = (gr_xend - gr_xbeg)/gr_nx

  ! allocate grid variables
  allocate(gr_U(NSYS_VAR,gr_imax)); gr_U = 0.
  allocate(gr_V(NUMB_VAR,gr_imax)); gr_V = 0.
  allocate(gr_W(NSYS_VAR,gr_imax)); gr_W = 0.

  ! allocate grid Riemann states
  allocate(gr_vR(NUMB_VAR,gr_imax)); gr_vR = 0.
  allocate(gr_vL(NUMB_VAR,gr_imax)); gr_vL = 0.
  !allocate(gr_vP(NUMB_VAR,gr_imax)); gr_vP = 0.
  !allocate(gr_vM(NUMB_VAR,gr_imax)); gr_vM = 0.

  ! allocate grid fluxes
  allocate(gr_flux(NSYS_VAR,gr_imax)); gr_flux = 0.
  !allocate(gr_CntrFlux(NSYS_VAR,gr_imax)); gr_CntrFlux = 0.

  ! allocate grid eigensystem
  allocate(gr_eigval(NUMB_WAVE,gr_imax)); gr_eigval = 0.
  allocate(gr_leigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_leigvc = 0.
  allocate(gr_reigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_reigvc = 0.
  allocate(gr_maxAlphas(NSYS_VAR)); gr_maxAlphas = 0.


  ! allocate GP variables
  allocate(gr_GPv(2*gr_radius+1   )); gr_GPv = 0.
  allocate(gr_GPZ(2, 2*gr_radius+1)); gr_GPZ = 0.
  return
end subroutine grid_init
