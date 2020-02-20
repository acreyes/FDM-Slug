module bc

#include "definition.h"

  use DMR
  use grid_data
  use sim_data, only : sim_xBC, sim_yBC, sim_bcTypex, sim_bcTypey, sim_reconmultiD
  implicit none

contains

  subroutine bc_init(bc_char, bc_int)
    implicit none

    character(len=MAX_STRING_LENGTH), intent(IN ) :: bc_char
    integer,                          intent(OUT) :: bc_int

    if (bc_char == 'periodic') then
       bc_int = PERIODIC
    elseif (bc_char == 'outflow') then
       bc_int = OUTFLOW
    elseif (bc_char == 'DMR') then
       bc_int = DBLMCH
    elseif (bc_char == 'outflow45') then
       bc_int = OUTFLOW45
    elseif (bc_char == 'inflow') then
       bc_int = INFLOW
    elseif (bc_char == 'reflect') then
       bc_int = REFLECT
    elseif (bc_char == 'jet') then
       bc_int = REFLECT
    end if

  end subroutine bc_init
!!$  subroutine bc_apply(V)
!!$    implicit none
!!$    real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: V
!!$    return
!!$  end subroutine bc_apply

  subroutine bc_apply(V)
    use block_data
    implicit none
    real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: V

    integer :: dest, gc, i, j

    real, dimension(NUMB_VAR,gr_ngc, gr_ny) :: loc_buffL, loc_buffR
    real, dimension(NUMB_VAR,gr_nx, gr_ngc) :: loc_buffT, loc_buffB

    
    
!!$    !load halo data into image buffers
    bl_buffL(:, :, :) = &
         V(:, gr_ibeg(XDIM):gr_ibeg(XDIM)+gr_ngc-1, gr_ibeg(YDIM):gr_iend(YDIM))
    bl_buffR(:, :, :) = &
         V(:, gr_iend(XDIM)-gr_ngc+1:gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))
    bl_buffB(:, :, :) = &
         V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_ibeg(YDIM)+gr_ngc-1)
    bl_buffT(:, :, :) = &
         V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_iend(YDIM)-gr_ngc+1:gr_iend(YDIM))

    Sync All

    !do halo exchanges along non-domain boundaries
    dest = bl_BC(1)
    if (dest > 0) then
       loc_buffL(:,:,:) = bl_buffR(:,:,:)[dest]
    end if
    
    dest = bl_BC(2)
    if (dest > 0) then
       loc_buffT(:,:,:) = bl_buffB(:,:,:)[dest]
    end if
    
    dest = bl_BC(3)
    if (dest > 0) then
       loc_buffR(:,:,:) = bl_buffL(:,:,:)[dest]
    end if
    
    dest = bl_BC(4)
    if (dest > 0) then
       loc_buffB(:,:,:) = bl_buffT(:,:,:)[dest]
    end if

    !now we apply domain boundary conditions
    if (sim_bcTypex == 'periodic') then
       !for periodic we can just handle exchange here
       if     (bl_i == 1        ) then
          dest = bl_grid(bl_iProcs, bl_j)
          loc_buffL(:,:,:) = bl_buffR(:,:,:)[dest]
       end if
       if (bl_i == bl_iProcs) then
          dest = bl_grid(1        , bl_j)
          loc_buffR(:,:,:) = bl_buffL(:,:,:)[dest]
       end if
    elseif (sim_bcTypex == 'outflow') then
       if (bl_i == 1) then
          !loc_buffL(:,:,:) = V(:,gr_i0(XDIM)+1:gr_ibeg(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))
          do i = gr_i0(XDIM),gr_ibeg(XDIM)-1
             loc_buffL(:,i,:) = V(:,gr_ibeg(XDIM),gr_ibeg(YDIM):gr_iend(YDIM))
          end do
       end if
       if (bl_i == bl_iProcs) then
          !loc_buffR(:,:,:) = V(:, gr_iend(XDIM):gr_imax(XDIM)-1, gr_ibeg(YDIM):gr_iend(YDIM))
          do i = gr_iend(XDIM)+1,gr_imax(XDIM)
             loc_buffR(:,i-gr_iend(XDIM),:) = V(:,gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM))
          end do
       end if
    elseif (sim_bcTypex == 'outflow45') then
       if (bl_i == 1) then
          loc_buffL(:,:,:) = V(:,gr_i0(XDIM)+1:gr_ibeg(XDIM), gr_ibeg(YDIM)+1:gr_iend(YDIM)+1)
       end if
       if (bl_i == bl_iProcs) then
          loc_buffR(:,:,:) = V(:, gr_iend(XDIM):gr_imax(XDIM)-1, gr_ibeg(YDIM)-1:gr_iend(YDIM)-1)
       end if
    elseif (sim_bcTypex == 'DMR') then
       if (bl_i == 1) then
          loc_buffL(:,:,:) = V(:,gr_i0(XDIM):gr_ibeg(XDIM)-1, gr_ibeg(YDIM):gr_iend(YDIM))
          call DMR_IN(loc_buffL(:,:,:))
       end if
       if (bl_i == bl_iProcs) then
          !do outflow
          do i = gr_iend(XDIM)+1,gr_imax(XDIM)
             loc_buffR(:,i-gr_iend(XDIM),:) = V(:,gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM))
          end do
       end if
    elseif (sim_bcTypex == 'reflect') then
       if (bl_i == 1) then
          loc_buffL(:,:,:) = V(:,gr_ibeg(XDIM)+gr_ngc:gr_ibeg(XDIM):-1,gr_ibeg(YDIM):gr_iend(YDIM))
          loc_buffL(VELX_VAR,:,:) = -loc_buffL(VELX_VAR,:,:)
       end if

       if (bl_i == bl_iProcs) then
          loc_buffR(:,:,:) = V(:,gr_iend(XDIM)-gr_ngc:gr_iend(XDIM):-1,gr_ibeg(YDIM):gr_iend(YDIM))
          loc_buffR(VELX_VAR,:,:) = -loc_buffR(VELX_VAR,:,:)
       end if
    elseif (sim_bcTypex == 'inflow') then
       if (bl_i == 1) then
          !let's leave as initial conditions on the inflow side
          loc_buffL(:,:,:) = V(:,gr_i0(XDIM):gr_ibeg(XDIM)-1,gr_ibeg(YDIM):gr_iend(YDIM))
       end if

       if (bl_i == bl_iProcs) then
          !outflow on the other side
          do i = gr_iend(XDIM)+1,gr_imax(XDIM)
             loc_buffR(:,i-gr_iend(XDIM),:) = V(:,gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM))
          end do
       end if
    elseif (sim_bcTypex == 'jet') then
          if (bl_i == 1) then
            loc_buffL(:,:,:) = V(:,gr_i0(XDIM):gr_ibeg(XDIM)-1,gr_ibeg(YDIM):gr_iend(YDIM)-1)
          end if

          if (bl_i == bl_iProcs) then
            do i = gr_iend(XDIM)+1,gr_imax(XDIM)
             loc_buffR(:,i-gr_iend(XDIM),:) = V(:,gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM))
          end do
       end if

    end if

    if (sim_bcTypey == 'periodic') then
       !again for periodic we can handle exchange here
       if (bl_j == 1) then
          dest = bl_grid(bl_i, bl_jProcs)
          loc_buffB(:,:,:) = bl_buffT(:,:,:)[dest]
       end if
       if (bl_j == bl_jProcs) then
          dest = bl_grid(bl_i, 1)
          loc_buffT(:,:,:) = bl_buffB(:,:,:)[dest]
       end if
    elseif (sim_bcTypey == 'outflow' ) then
       if (bl_j == 1) then
          !loc_buffB(:,:,:) = V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_i0(YDIM)+1:gr_ibeg(YDIM))
          do i = gr_i0(YDIM), gr_ibeg(YDIM)-1
             loc_buffB(:,:,i) = V(:,gr_ibeg(XDIM):gr_iend(XDIM),gr_ibeg(YDIM))
          end do
       end if
       if (bl_j == bl_jProcs) then
          !loc_buffT(:,:,:) = V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_iend(YDIM):gr_imax(YDIM)-1)
          do i = gr_iend(YDIM)+1, gr_imax(YDIM)
             loc_buffT(:,:,i-gr_iend(YDIM)) = V(:,gr_ibeg(XDIM):gr_iend(XDIM),gr_iend(YDIM))
          end do
       end if
    elseif (sim_bcTypey == 'outflow45' ) then
       if (bl_j == 1) then
          loc_buffB(:,:,:) = V(:, gr_ibeg(XDIM)+1:gr_iend(XDIM)+1, gr_i0(YDIM)+1:gr_ibeg(YDIM))
       end if
       if (bl_j == bl_jProcs) then
          loc_buffT(:,:,:) = V(:, gr_ibeg(XDIM)-1:gr_iend(XDIM)-1, gr_iend(YDIM):gr_imax(YDIM)-1)
       end if
    elseif (sim_bcTypey == 'DMR') then
       if (bl_j == 1) then
          loc_buffB(:,:,:) = V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_ibeg(YDIM)+gr_ngc)
          call DMR_bot(loc_buffB(:,:,:))
       end if
       if (bl_j == bl_jProcs) then
          loc_buffT(:,:,:) = V(:,gr_ibeg(XDIM):gr_iend(XDIM), gr_iend(YDIM)+1:gr_imax(YDIM))
          call DMR_top(loc_buffT(:,:,:))
       end if

    elseif(sim_bcTypey == 'reflect') then
       if (bl_j == 1) then
          do i = gr_ibeg(XDIM), gr_iend(XDIM)
             loc_buffB(:,i-gr_ngc,:) = V(:,i,gr_ibeg(YDIM)+gr_ngc:gr_ibeg(YDIM):-1)
             loc_buffB(VELY_VAR,i-gr_ngc,:) = -loc_buffB(VELY_VAR,i-gr_ngc,:)
          end do
       end if

       if (bl_j == bl_jProcs) then
          loc_buffT(:,:,:) = V(:,gr_ibeg(XDIM):gr_iend(XDIM),gr_iend(YDIM):gr_iend(YDIM)-gr_ngc:-1)
          loc_buffT(VELY_VAR,:,:) = -loc_buffT(VELY_VAR,:,:)
       end if
    end if
       
       

    !fill halo data
    !do this guard cell by guard cell so that it fills int he same order
    V(:, gr_i0(XDIM):gr_ibeg(XDIM)-1, gr_ibeg(YDIM):gr_iend(YDIM)) = &
         loc_buffL(:,:,:)
    V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_iend(YDIM)+1:gr_imax(YDIM)) = &
         loc_buffT(:,:,:)
    V(:, gr_iend(XDIM)+1:gr_imax(XDIM), gr_ibeg(YDIM):gr_iend(YDIM)) = &
         loc_buffR(:,:,:)
    V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_i0(YDIM):gr_ibeg(YDIM)-1)  = &
         loc_buffB(:,:,:)

    if (sim_reconMultiD) call bc_corners(V)

!!$#ifdef BDRY_VAR
!!$    if (sim_entFix) then
!!$       call wind_entFix(V)
!!$    end if
!!$#endif  
    sync all
  end subroutine bc_apply

  subroutine bc_outflowX(V, loc_buffL, loc_buffR)
    implicit none
    real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: V
    real, dimension(NUMB_VAR, gr_ngc, gr_ny), intent(INOUT) :: loc_buffL, loc_buffR
    real, dimension(NUMB_VAR) :: Vbeg, Vend
    integer :: i,j

    !x-region
!!$    loc_buffL(:,:,:) = V(:, gr_iend(XDIM):gr_imax(XDIM)-1)
!!$    do i = 1, gr_ngc
!!$       loc_buffL
!!$    end do
!!$    do j = gr_ibeg(YDIM),gr_iend(YDIM)!1,gr_ny
!!$       do i = 1,gr_ngc
!!$          Vbeg(1:NUMB_VAR) = V(1:NUMB_VAR, gr_ibeg(XDIM)+i-1,j)
!!$          Vend(1:NUMB_VAR) = V(1:NUMB_VAR, gr_iend(XDIM)-i+1,j)
!!$
!!$          V(1:NUMB_VAR, gr_iend(XDIM)+i,j) = Vbeg(1:NUMB_VAR)
!!$          V(1:NUMB_VAR, gr_ibeg(XDIM)-i,j) = Vend(1:NUMB_VAR)
!!$       end do
!!$    end do

       return

     end subroutine bc_outflowX

  
  subroutine bc_periodic(V)
    !I need to think about how to deal with the corner cases only affects GP
    implicit none
    real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: V
    real, dimension(NUMB_VAR) :: Vbeg, Vend, Vne, Vnw, Vse, Vsw
    integer :: i,j

    !first do periodic in x
    do j = gr_ibeg(YDIM),gr_iend(YDIM)!1,gr_ny
       do i = 1,gr_ngc
          Vbeg(1:NUMB_VAR) = V(1:NUMB_VAR, gr_ibeg(XDIM)+i-1,j)
          Vend(1:NUMB_VAR) = V(1:NUMB_VAR, gr_iend(XDIM)-i+1,j)

          V(1:NUMB_VAR, gr_iend(XDIM)+i,j) = Vbeg(1:NUMB_VAR)
          V(1:NUMB_VAR, gr_ibeg(XDIM)-i,j) = Vend(1:NUMB_VAR)
       end do
    end do

    !now periodic in y
    do i = gr_ibeg(XDIM),gr_iend(XDIM)!1,gr_nx
       do j = 1,gr_ngc
          Vbeg(1:NUMB_VAR) = V(1:NUMB_VAR, i, gr_ibeg(YDIM)+j-1)
          Vend(1:NUMB_VAR) = V(1:NUMB_VAR, i, gr_iend(YDIM) - j+1)

          V(1:NUMB_VAR, i,gr_iend(YDIM)+j) = Vbeg(1:NUMB_VAR)
          V(1:NUMB_VAR, i,gr_ibeg(YDIM)-j) = Vend(1:NUMB_VAR)
       end do
    end do

    !deal with the corners
    do i = 1, gr_ngc
       do j = 1, gr_ngc
          Vnw(:) = V(1:NUMB_VAR, gr_ibeg(XDIM) + i-1, gr_iend(YDIM) - j+1)
          Vne(:) = V(1:NUMB_VAR, gr_iend(XDIM) - i+1, gr_iend(YDIM) -j+1)
          Vsw(:) = V(1:NUMB_VAR, gr_ibeg(XDIM) + i-1, gr_ibeg(YDIM) + j-1)
          Vse(:) = V(1:NUMB_VAR, gr_iend(XDIM) - i+1, gr_ibeg(YDIM) + j-1)

          V(1:NUMB_VAR, gr_ibeg(XDIM) -i, gr_iend(YDIM) + j) = Vse !nw->se
          V(1:NUMB_VAR, gr_iend(XDIM) +i, gr_iend(YDIM) + j) = Vsw !ne->sw
          V(1:NUMB_VAR, gr_ibeg(XDIM) -i, gr_ibeg(YDIM) - j) = Vne !sw->ne
          V(1:NUMB_VAR, gr_iend(XDIM) +i, gr_ibeg(YDIM) - j) = Vnw !se->nw
       end do
    end do
    
  end subroutine bc_periodic

!!$  subroutine bc_user(V)
!!$    implicit none
!!$    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
!!$    !BC for shu-osher problem.
!!$    !don't do anything!
!!$  end subroutine bc_user
!!$  
!!$  subroutine bc_outflow(V)
!!$    implicit none
!!$    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
!!$    integer :: i
!!$ 
!!$    do i = 1, gr_ngc
!!$       ! on the left GC
!!$       V(1:NUMB_VAR,i) = V(1:NUMB_VAR,gr_ibeg) 
!!$
!!$       ! on the right GC
!!$       V(1:NUMB_VAR,gr_iend+i) = V(1:NUMB_VAR,gr_iend)
!!$    end do
!!$
!!$    return
!!$  end subroutine bc_outflow
!!$
!!$  subroutine bc_reflect(V)
!!$    implicit none
!!$    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
!!$    integer :: i,k0,k1
!!$
!!$    do i = 1, gr_ngc
!!$       k0 = 2*gr_ngc+1
!!$       k1 = gr_iend-gr_ngc
!!$
!!$       ! on the left GC
!!$       V(       :,i) = V(       :,k0-i)
!!$       V(VELX_VAR,i) =-V(VELX_VAR,k0-i)
!!$
!!$       ! on the right GC
!!$       V(       :,k1+k0-i) = V(         :,k1+i)
!!$       V(VELX_VAR,k1+k0-i) =-V(  VELX_VAR,k1+i)
!!$    end do
!!$
!!$    return
!!$  end subroutine bc_reflect


end module bc
