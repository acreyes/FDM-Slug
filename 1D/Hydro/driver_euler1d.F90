program driver_euler1d

#include "definition.h"
  
  use sim_data
  use grid_data
  use io, only : io_writeOutput
  use bc
  use eos, only : eos_all

  implicit none

  real :: t,dt
  integer :: nStep,ioCounter,ioTimeFreqCounter
  real :: ioCheckTime
  
  t = 0.
  nStep = 0
  ioCounter = 0
  ioTimeFreqCounter = 0
  
  
  ! grid_init should be called first before sim_init
  call read_pars()
  call grid_init()
  call sim_init()

  write(*,*)''
  write(*,*)'================================================='
  write(*,*)'     AMS 260 - CFD, 1D Euler FVM Code            '
  write(*,*)'      Written by Prof. Dongwook Lee              '
  write(*,*)'           Winter Quarter, 2015                  '
  write(*,*)'================================================='
  write(*,*)''

  
  ! write the initial condition
  write(*,*)''
  write(*,*)'       Initial condition was written!            '
  write(*,*)'================================================='
  write(*,*)'   Steps      Time              dt               '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter)

  if (sim_order == 10) then
     !we're doing GP and need to initialize
     call sim_GPinit
     call gp_Fluxinit
     call gp_eigens
  end if

  do while ( (t < sim_tmax))
     
     if (sim_fixDt) then
        dt = sim_dt
     else
        call cfl(dt)
     end if

     if (dt .le. 0.) then
        exit
     end if
     !check to see if there is a reason to stop
     if (  sim_nlim .and. (nStep .ge. sim_nstep)) then
        exit
     elseif ( abs(t - sim_tmax) .le. dt ) then
        dt = abs(t - sim_tmax)
        !exit
     end if
     call soln_ReconEvolveAvg(dt)
     call soln_update(dt)


     ! call BC on primitive vars
     call bc_apply(gr_V)
     
     ! call eos to make sure all through GC regions
     !call eos_all
     !lets remove this because dongwook say its better
     
     ! write outputs every ioNfreq cycle or ioTfreq cycle
     
     ioCheckTime = sim_ioTfreq*real(ioTimeFreqCounter+1)
     if (t-dt < ioCheckTime .and. t>ioCheckTime) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        ioTimeFreqCounter = ioTimeFreqCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif

     if (sim_ioNfreq > 0) then
     if (mod(nStep, sim_ioNfreq) == 0) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif
     endif
     
     ! update your time and step count
     t = t + dt
     nStep = nStep + 1

     write(*,900)nstep,t,dt
     !pm = pm*-1.
  enddo


  !! Let's write the final result before exiting
  write(*,*)''
  write(*,*)' Final output no.',ioCounter+1, 'has been written'
  write(*,*)'================================================='
  write(*,*)'        The final tmax has reached, bye!         '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter+1)

  !! finalize and deallocate memories
  call grid_finalize()

900 format(1x,i5,f16.8,1x,f16.8)
  
  return
end program driver_euler1d
