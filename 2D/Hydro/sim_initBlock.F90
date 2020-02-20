subroutine sim_initBlock()
 
#include "definition.h"
  
  use sim_data
  use grid_data, only : gr_V,gr_U,gr_ngc,gr_xbeg &
                       ,gr_i0,gr_imax,gr_xCoord,gr_dx, gr_yCoord, gr_dy
  use block_data, only   : bl_i, bl_j, bl_nBlock  
  use primconsflux, only : prim2cons
  use bc
  
  implicit none

  integer :: i,j,in,jn
  real :: ekin, eint, x, y, ranx, rany, small, r2, beta, T, dr2, E, gamm
  real :: rt3, xmin, x0, ymin, sig, thkns, lngth, loc

  small = 0.01
  sig = 0.05/SQRT(2.)
  call RANDOM_SEED()
  
  ! generate x-coordinate
!!$  do i = gr_i0(XDIM),gr_imax(XDIM)
!!$     gr_xCoord(i) = (real(i-gr_ngc + (bl_i-1)*bl_nBlock(XDIM)))*gr_dx + gr_xbeg(XDIM)
!!$  end do
!!$
!!$  ! generate y-coordinates
!!$  do i = gr_i0(YDIM),gr_imax(YDIM)
!!$     gr_yCoord(i) = (real(i-gr_ngc+(bl_j-1)*bl_nBlock(YDIM)))*gr_dy + gr_xbeg(YDIM)
!!$  end do
!!$
  do i = gr_i0(XDIM),gr_imax(XDIM)
     gr_xCoord(i) = (real(i-gr_ngc + (bl_i-1)*bl_nBlock(XDIM))-0.5)*gr_dx + gr_xbeg(XDIM)
  end do

  ! generate y-coordinates
  do i = gr_i0(YDIM),gr_imax(YDIM)
     gr_yCoord(i) = (real(i-gr_ngc+(bl_j-1)*bl_nBlock(YDIM))-0.5)*gr_dy + gr_xbeg(YDIM)
  end do

  !some handy constants
  gamm = sim_gamma
  !these are used for sedov
  dr2 = (3.5*MIN(gr_dx, gr_dy))**2
  E = 1.

  rt3 = sqrt(3.)
  xmin = 1./6.
  x0 = xmin + 1./rt3
  
  do i = gr_i0(XDIM),gr_imax(XDIM)
     do j = gr_i0(YDIM),gr_imax(YDIM)
#ifdef BDRY_VAR
        gr_V(BDRY_VAR,i,j) = -1.0
#endif        
        if (sim_icType == 'gaussx') then
           !advect gaussian in the x-direction
           x = gr_xCoord(i)-0.5
           y = gr_yCoord(j)-0.5
           gr_V(DENS_VAR,i,j) = 1. + EXP(-100.*(y**2))
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 1.
           gr_V(PRES_VAR,i,j) = 1./sim_gamma

        elseif (sim_icType == 'gaussxy') then
           x = gr_xCoord(i)-0.5
           y = gr_yCoord(j)-0.5
           
           gr_V(DENS_VAR,i,j) = 1. + EXP(-(x**2+y**2)/(0.1)**2)
           gr_V(VELX_VAR,i,j) = 1.
           gr_V(VELY_VAR,i,j) = 1.
           gr_V(PRES_VAR,i,j) = 1./sim_gamma

        elseif (sim_icType == 'DMR') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           ymin = (x-xmin)*rt3
           if (y > ymin) then
              !in the shock region
              gr_V(DENS_VAR,i,j) = 8.
              gr_V(VELX_VAR,i,j) = 7.1447096
              gr_V(VELY_VAR,i,j) = -4.125
              gr_V(PRES_VAR,i,j) = 116.5
           else
              !outside shock
              gr_V(DENS_VAR,i,j) = 1.4
              gr_V(VELX_VAR,i,j) = 0.
              gr_V(VELY_VAR,i,j) = 0.
              gr_V(PRES_VAR,i,j) = 1.
           end if
           
        elseif (sim_icType == '2DRP') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           if (x > sim_x0 .and. y .ge. sim_y0) then
              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q1(:)
           elseif (x .le. sim_x0 .and. y .ge. sim_y0) then
              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q2(:)
           elseif ( x .le. sim_x0 .and. y < sim_y0) then
              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q3(:)
           elseif (x > sim_x0 .and. y < sim_y0) then
              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q4(:)
           end if
           
        elseif (sim_icType == 'shocky') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           if (y < 0.5) then
              gr_V(DENS_VAR,i,j) = 1.
              gr_V(PRES_VAR,i,j) = 1.
           else
              gr_V(DENS_VAR,i,j) = 0.125
              gr_V(PRES_VAR,i,j) = 0.1
           end if
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.
        elseif (sim_icType == 'shockx') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           if (x < 0.5) then
              gr_V(DENS_VAR,i,j) = 1.
              gr_V(PRES_VAR,i,j) = 1.
           else
              gr_V(DENS_VAR,i,j) = 0.125
              gr_V(PRES_VAR,i,j) = 0.1
           end if
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.
        elseif (sim_icType == 'shockxy') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           r2 = x**2 + y**2
           if (x < y) then
              gr_V(DENS_VAR,i,j) = 1.
              gr_V(PRES_VAR,i,j) = 1.
           else
              gr_V(DENS_VAR,i,j) = 1.
              gr_V(PRES_VAR,i,j) = 1.
           end if
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.
        elseif(sim_icType == 'shock_periodic') then
           x = gr_xCoord(i)+.5*sqrt(2.)*gr_dx
           y = gr_yCoord(j)+.5*sqrt(2.)*gr_dx
           r2 = 1./sqrt(2.)*(x+y)

           if ((r2 .le. -1.5) .or. (r2 > -0.5 .and. r2 .le. 0.5) .or. (r2 > 1.5)) then
              !initial "left" state
              gr_V(DENS_VAR,i,j) = 1.
              gr_V(PRES_VAR,i,j) = 1.
           else
              !right state
              gr_V(DENS_VAR,i,j) = 0.125
              gr_V(PRES_VAR,i,j) = 0.1
           end if
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.

        elseif (sim_icType == 'vortex') then
           beta = 5.
           x = gr_xCoord(i) - 10.
           y = gr_yCoord(j) - 10.
           r2 = x**2+y**2
           T = 1. - (sim_gamma-1.)*beta*beta*EXP(1.-r2)/(8.*sim_gamma*PI*PI)
           
           gr_V(DENS_VAR,i,j) = (T)**(1./(sim_gamma-1.))
           gr_V(VELX_VAR,i,j) = 1. - y*beta/(2.*PI)*EXP(0.5*(1.-r2))
           gr_V(VELY_VAR,i,j) = 1. + x*beta/(2.*PI)*EXP(0.5*(1.-r2))
           gr_V(PRES_VAR,i,j) = gr_V(DENS_VAR,i,j)*T

        elseif (sim_icType == 'sedov') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           r2 = x**2 + y**2
           if (sqrt(r2) < sqrt(dr2)) then
              gr_V(PRES_VAR,i,j) = (gamm-1.)*E/(PI*dr2)
           else
              gr_V(PRES_VAR,i,j) = 1.e-5
           end if

           gr_V(DENS_VAR,i,j) = 1.
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.

        elseif (sim_icType == 'sedovChamber') then
           loc = .25
           thkns = 0.025
           lngth = .15
           x = gr_xCoord(i)
           y = gr_yCoord(j)
           r2 = x**2 + y**2
           if (sqrt(r2) < sqrt(dr2)) then
              gr_V(PRES_VAR,i,j) = (gamm-1.)*E/(PI*dr2)
           else
              gr_V(PRES_VAR,i,j) = 1.e-5
           end if

           gr_V(DENS_VAR,i,j) = 1.
           gr_V(VELX_VAR,i,j) = 0.
           gr_V(VELY_VAR,i,j) = 0.

#ifdef BDRY_VAR
           if ( (abs(x-loc) .le. thkns) .and. (abs(y) .le. lngth) ) then
              gr_V(BDRY_VAR,i,j) = 1.
              gr_V(DENS_VAR,i,j) = 0.
           elseif ( (abs(x+loc) .le. thkns) .and. (abs(y) .le. lngth) ) then
              gr_V(BDRY_VAR,i,j) = 1.
              gr_V(DENS_VAR,i,j) = 0.
           elseif ( (abs(y+loc) .le. thkns) .and. (abs(x) .le. lngth) ) then
              gr_V(BDRY_VAR,i,j) = 1.
              gr_V(DENS_VAR,i,j) = 0.
           elseif ( (abs(y-loc) .le. thkns) .and. (abs(x) .le. lngth) ) then
              gr_V(BDRY_VAR,i,j) = 1.
              gr_V(DENS_VAR,i,j) = 0.
           else
              gr_V(BDRY_VAR,i,j) = -1.0
           endif
#endif

        elseif (sim_icType == 'KH') then
           x = gr_xCoord(i)
           y = gr_yCoord(j)

           if (abs(y) > 0.25) then
              gr_V(VELX_VAR,i,j) = -0.5
              gr_V(DENS_VAR,i,j) = 1.
           else
              gr_V(VELX_VAR,i,j) = 0.5
              gr_V(DENS_VAR,i,j) = 2.
           end if
           gr_V(PRES_VAR,i,j) = 2.5
           gr_V(VELX_VAR,i,j) = gr_V(VELX_VAR,i,j) + sim_boost !+ small*ranx
           gr_V(VELY_VAR,i,j) = sim_boost +  0.1*SIN(4.*PI*x)*(EXP(-(y+0.25)**2/(0.05**2)) + EXP(-(y-0.25)**2/(0.05**2)))!small*rany


           elseif (sim_icType == 'windtunnel') then
              x = gr_xCoord(i)
              y = gr_yCoord(j)
              ! set BDRY_VAR = 1 for solid cells
              ! The step is located at
              ! 0.6 <= x <= 3 & 0 <= y <= 0.2
              ! the entire domain is [0, 3] x [0, 1]
              if ((x .ge. 0.6) .and. (y .le. 0.2) ) then
                 gr_V(DENS_VAR,i,j) = 10.4
                 gr_V(PRES_VAR,i,j) = 1.0
                 gr_V(VELX_VAR,i,j) = 0.0
                 gr_V(VELY_VAR,i,j) = 0.0
#ifdef BDRY_VAR                 
                 gr_V(BDRY_VAR,i,j) = 1.
                 gr_V(DENS_VAR,i,j) = 0.
#endif
                 
              else
                 ! these are fluid cells
#ifdef BDRY_VAR   
                 gr_V(BDRY_VAR,i,j) = -1.0
#endif
                 gr_V(DENS_VAR,i,j) = 1.4
                 gr_V(PRES_VAR,i,j) = 1.0
                 gr_V(VELX_VAR,i,j) = 3.0
                 gr_V(VELY_VAR,i,j) = 0.0
              endif
           

           elseif (sim_icType == 'highmach') then
              gr_V(DENS_VAR,i,j) = 5.0
              gr_V(VELX_VAR,i,j) = 0.0
              gr_V(VELY_VAR,i,j) = 0.0
              gr_V(PRES_VAR,i,j) = 0.4217

           elseif (sim_icType == 'shockdefraction') then
              x = gr_xCoord(i)
              y = gr_yCoord(j)


           end if

        gr_V(GAMC_VAR,i,j) = sim_gamma
        gr_V(GAME_VAR,i,j) = sim_gamma
        gr_V(EINT_VAR,i,j) = gr_V(PRES_VAR,i,j)/(gr_V(GAME_VAR,i,j)-1.)/gr_V(DENS_VAR,i,j)
        
     end do
  end do
  !initialize cons vars
  do i = gr_i0(XDIM),gr_imax(XDIM)
     do j = gr_i0(YDIM),gr_imax(YDIM)
        call prim2cons(gr_V(:,i,j), gr_U(DENS_VAR:ENER_VAR,i,j))
     end do
  end do

  call bc_apply(gr_V)
  
  
!!$  do i = gr_i0,gr_imax
!!$     if (sim_icType == 'shock') then
!!$        if (gr_xCoord(i) < sim_shockLoc) then
!!$           gr_V(DENS_VAR,i) = sim_densL
!!$           gr_V(VELX_VAR,i) = sim_velxL
!!$           gr_V(PRES_VAR,i) = sim_presL
!!$        else
!!$           gr_V(DENS_VAR,i) = sim_densR
!!$           gr_V(VELX_VAR,i) = sim_velxR
!!$           gr_V(PRES_VAR,i) = sim_presR
!!$        end if
!!$
!!$     elseif (sim_icType == 'blast') then
!!$        !do blas IC
!!$        if (gr_xCoord(i) <= 0.9 .AND. gr_xCoord(i) > .1) then
!!$           !middle state
!!$           gr_V(DENS_VAR,i) = 1.
!!$           gr_V(VELX_VAR,i) = 0.
!!$           gr_V(PRES_VAR,i) = .01
!!$        elseif (gr_xCoord(i) <= .1) then
!!$           gr_V(DENS_VAR,i) = 1.
!!$           gr_V(VELX_VAR,i) = 0.
!!$           gr_V(PRES_VAR,i) = 1000.
!!$        else
!!$           gr_V(DENS_VAR,i) = 1.
!!$           gr_V(VELX_VAR,i) = 0.
!!$           gr_V(PRES_VAR,i) = 100.
!!$        endif
!!$           
!!$     elseif (sim_icType == 'shu') then
!!$        !do IC for shu-osher problem
!!$        !transform the domain [0,1] onto the one give for the problem:[-4.5,4.5]
!!$        x = gr_xCoord(i) - 4.5
!!$        if (x < -4.) then
!!$           !left state
!!$           gr_V(DENS_VAR,i) = 3.857143
!!$           gr_V(VELX_VAR,i) = 2.629369
!!$           gr_V(PRES_VAR,i) = 10.33333
!!$        else
!!$           gr_V(DENS_VAR,i) = 1 + .2*SIN(5.*x)
!!$           gr_V(VELX_VAR,i) = 0.
!!$           gr_V(PRES_VAR,i) = 1.
!!$        end if
!!$
!!$     elseif (sim_icType == 'gauss') then
!!$        x = gr_xCoord(i) - .5
!!$        gr_V(DENS_VAR,i) = 1. + EXP(-100.*x**2)
!!$        gr_V(VELX_VAR,i) = -1.
!!$        gr_V(PRES_VAR,i) = 1./sim_gamma
!!$
!!$     elseif (sim_icType == 'sine') then
!!$        x = gr_xCoord(i)
!!$        gr_V(DENS_VAR,i) = 1.5 + SIN(2.*PI*x)
!!$        gr_V(VELX_VAR,i) = -1.
!!$        gr_V(PRES_VAR,i) = 1./sim_gamma
!!$        
!!$     end if
!!$     gr_V(GAMC_VAR,i) = sim_gamma
!!$     gr_V(GAME_VAR,i) = sim_gamma
!!$     gr_V(EINT_VAR,i) = gr_V(PRES_VAR,i)/(gr_V(GAME_VAR,i)-1.)/gr_V(DENS_VAR,i)
!!$     
!!$  end do

!!$  ! also initialize conservative vars
!!$  do i = gr_i0,gr_imax
!!$     call prim2cons(gr_V(:,i), gr_U(DENS_VAR:ENER_VAR,i))
!!$  end do

  
end subroutine sim_initBlock
