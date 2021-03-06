module io
 
#include "definition.h"
  
  use grid_data, only : gr_xCoord, gr_ibeg, gr_iend &
       ,gr_yCoord, gr_V, gr_nx, gr_ny &
       , gr_glb_nx, gr_glb_ny, gr_xbeg &
       , gr_dx, gr_dy, gr_betas
  
  use sim_data, only : sim_name, sim_hdf5
  use block_data
  implicit none

  
  integer, save :: nCounter
  
contains

  subroutine io_writeOutput(t,nstep,ioCounter,eTime)

    implicit none
    real, intent(IN) :: t, eTime
    integer, intent(IN) :: nstep, ioCounter
    if (sim_hdf5) then
       call io_writeHDF5(t,nstep,ioCounter,eTime)
    else
       call io_writeASCII(t,nstep,ioCounter)
    end if
  end subroutine io_writeOutput

  subroutine io_writeHDF5(t,nstep,ioCounter,eTime)

    use HDF5
    implicit none

    real, intent(IN) :: t, eTime
    integer, intent(IN) :: nstep,ioCounter

    

    integer :: i,j,nVar,nCell,dest
    character(len=50) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend
    real, dimension(gr_glb_nx) :: xCoord, xread
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny) :: V, V_read
    integer, dimension(NDIM)   :: img_size
    
    integer, allocatable :: nBlock(:)[:]
    real,    allocatable :: img_V(:,:,:)[:], loc_V(:,:,:)

    character(len=50) :: dset_prim, dset_x, dset_y
    integer(HID_T)    :: file_id, dspace_id, dset_id
    integer           :: error, rank_V, rank_XYZ, rankT
    integer(HSIZE_T), dimension(1) :: dimT
    integer(HSIZE_T),allocatable, dimension(:) :: dims_V, dims_XYZ, dims_B

    !make image buffers for grid
    allocate(img_V(NUMB_VAR, gr_nx, gr_ny)[*])
    allocate(nBlock(NDIM)[*])
    
    nBlock(XDIM) = gr_nx
    nBlock(YDIM) = gr_ny
    img_V(:, :, :) = gr_V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))

    Sync All

    if (bl_ID == 1) then
       !I am root
       dset_prim = "prim_vars"
       dset_x    = "xCoord"
       dset_y    = "yCoord"
       do i = 1, gr_glb_nx
          xCoord(i) = ( real(i)-0.5 )*gr_dx + gr_xbeg(XDIM)
       end do
       do j = 1, gr_glb_ny
          yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
       end do
       do i = 1, bl_iProcs
          do j = 1, bl_jProcs
             dest = bl_grid(i,j)
             img_size(:) = nBlock(:)[dest]

             ibeg = (i-1)*bl_nBlock(XDIM) + 1
             jbeg = (j-1)*bl_nBlock(YDIM) + 1
             iend = ibeg -1 + img_size(XDIM)
             jend = jbeg -1 + img_size(YDIM)
             
             allocate(loc_V(NUMB_VAR,img_size(XDIM),img_size(YDIM)))
             loc_V(:,:,:) = img_V(:,:,:)[dest]
             V(:, ibeg:iend, jbeg:jend) = loc_V(:,:,:)
             deallocate(loc_V)
          end do
       end do
       
       ! convert conter number to character
       write(cCounter,910) ioCounter + 10000

       ! file name for ascii output
       !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
       ofile = trim(sim_name)//'_'//cCounter//'.slug'
       allocate(dims_V(3))
       allocate(dims_B(4))
       allocate(dims_XYZ(1))
       rank_XYZ = 1
       dims_V = (/ NSYS_VAR,  gr_glb_nx, gr_glb_ny/)
       rank_V = 3
       
       !hdf5 fortran interface
       call h5open_f(error)
       !create data file
       call h5fcreate_f(ofile, H5F_ACC_TRUNC_F, file_id, error)

       !datatspace for primitive vars
       call h5screate_simple_f(rank_V, dims_V, dspace_id, error)
       !dataset for primitive vars
       call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write prim vars
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, V(DENS_VAR:NSYS_VAR,:,:), dims_V, error)
       
       !close dataset
       call h5dclose_f(dset_id, error)
       !close dataspace
       call h5sclose_f(dspace_id, error)

       
       !gp betas
       rank_V = 4
       dims_B = (/ NSYS_VAR, gr_glb_nx, gr_glb_ny, NDIM /)
       dset_prim = 'betas'
       call h5screate_simple_f(rank_V, dims_B, dspace_id, error)
       call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gr_betas(:,gr_ibeg(XDIM):gr_iend(XDIM),gr_ibeg(YDIM):gr_iend(YDIM),:), dims_B, error)

        !close dataset
       call h5dclose_f(dset_id, error)
       !close dataspace
       call h5sclose_f(dspace_id, error)

       !dataspace for x-coordinates
       dims_XYZ = (/gr_glb_nx/)
       call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
       call h5dcreate_f(file_id, dset_x, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xCoord, dims_XYZ, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)

       !dataspace for x-coordinates
       dims_XYZ = gr_glb_ny
       call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
       call h5dcreate_f(file_id, dset_y, h5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yCoord, dims_XYZ, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)

       !dataspace for elapsed cpu time
       rankT = 1
       dimT  = 1
       call h5screate_simple_f(rankT,dimT,dspace_id,error)
       call h5dcreate_f(file_id, 'eTime', h5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eTime, dimT, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)
      



       !close hdf5 file and interface
       call h5fclose_f(file_id,error)
       call h5close_f(error)

       deallocate(dims_V)
       deallocate(dims_XYZ)

       !call h5screatet_simple_f(rank, dims1)
       
!!$       open(unit=20,file=ofile,status='unknown')
!!$       do i=1, gr_glb_nx
!!$          do j = 1, gr_glb_ny
!!$             write(20,920)xCoord(i),ycoord(j),(V(nVar,i,j),nVar=1,EINT_VAR)
!!$          end do
!!$          write(20,100)
!!$       end do

    end if
    Sync All
    deallocate(img_V)
    deallocate(nBlock)

100 format()
910 format(i5)
920 format(1x,f16.8,1x,f16.8,1x,NUMB_VAR f32.16)
    
    close(20)

  end subroutine io_writeHDF5

  subroutine io_writeASCII(t,nstep,ioCounter)
    implicit none

    real, intent(IN) :: t
    integer, intent(IN) :: nstep,ioCounter
    
    integer :: i,j,nVar,nCell,dest
    character(len=50) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend
    real, dimension(gr_glb_nx) :: xCoord
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny) :: V
    integer, dimension(NDIM)   :: img_size
    
    integer, allocatable :: nBlock(:)[:]
    real,    allocatable :: img_V(:,:,:)[:], loc_V(:,:,:)

    !make image buffers for grid
    allocate(img_V(NUMB_VAR, gr_nx, gr_ny)[*])
    allocate(nBlock(NDIM)[*])
    
    nBlock(XDIM) = gr_nx
    nBlock(YDIM) = gr_ny
    img_V(:, :, :) = gr_V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))

    Sync All

    if (bl_ID == 1) then
       !I am root
       do i = 1, gr_glb_nx
          xCoord(i) = ( real(i)-0.5 )*gr_dx + gr_xbeg(XDIM)
       end do
       do j = 1, gr_glb_ny
          yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
       end do
       do i = 1, bl_iProcs
          do j = 1, bl_jProcs
             dest = bl_grid(i,j)
             img_size(:) = nBlock(:)[dest]

             ibeg = (i-1)*bl_nBlock(XDIM) + 1
             jbeg = (j-1)*bl_nBlock(YDIM) + 1
             iend = ibeg -1 + img_size(XDIM)
             jend = jbeg -1 + img_size(YDIM)
             
             allocate(loc_V(NUMB_VAR,img_size(XDIM),img_size(YDIM)))
             loc_V(:,:,:) = img_V(:,:,:)[dest]
             V(:, ibeg:iend, jbeg:jend) = loc_V(:,:,:)
             deallocate(loc_V)
          end do
       end do
       
       ! convert conter number to character
       write(cCounter,910) ioCounter + 10000

       ! file name for ascii output
       !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
       ofile = trim(sim_name)//'_'//cCounter//'.dat'
       
       open(unit=20,file=ofile,status='unknown')
       do i=1, gr_glb_nx
          do j = 1, gr_glb_ny
             write(20,920)xCoord(i),ycoord(j),(V(nVar,i,j),nVar=1,EINT_VAR)
          end do
          write(20,100)
       end do

    end if
    Sync All
    deallocate(img_V)
    deallocate(nBlock)

100 format()
910 format(i5)
920 format(1x,f16.8,1x,f16.8,1x,NUMB_VAR f32.16)
    
    close(20)
  end subroutine io_writeASCII
  
end module io
