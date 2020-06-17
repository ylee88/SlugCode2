program SlugCode

#include "definition.h"

  use sim_data
  use grid_data
  use block_data, only : bl_ID
  use io,         only : io_writeOutput
  use bc,         only : bc_apply
  use mpi,        only : MPI_Wtime

  implicit none

  real :: t, dt, dt_small, dt_init
  integer :: nStep, ioCounter, ioTimeFreqCounter
  real :: ioCheckTime, start, finish
  character(len=MAX_STRING_LENGTH) :: init_file_name

  t = 0.
  nStep = 0
  ioCounter = 0
  ioTimeFreqCounter = 0
  dt_init = 2.e-10

  ! for DMR bc, we need current time
  sim_t = t

  ! read init file as an argument
  call GETARG(1, init_file_name)
  if (init_file_name == '') then
    init_file_name = './slug.init'
  end if

  ! grid_init should be called first before sim_init
  call read_pars(init_file_name)

  call block_init()
  call grid_init()
  call bc_init()
  call sim_init()

  if (bl_id == 0) then
    write(*,*)''
    write(*,*)'================================================='
    write(*,*)'                  SlugCode2                      '
    write(*,*)'      Applied Math Dept., UC Santa Cruz          '
    write(*,*)'               Prof. Dongwook Lee                '
    write(*,*)'================================================='
    write(*,*)''
    write(*,*)'   Read parameters from ', '[', trim(init_file_name), ']   '
    write(*,*)''

    ! write the initial condition
    write(*,*)''
    write(*,*)'    Initial condition was written!               '
    write(*,*)'================================================='
    write(*,*)'   Steps      Time              dt               '
    write(*,*)'================================================='
    write(*,*)''
  end if

  ! write initial condition
  call io_writeOutput(t, nStep, ioCounter, 0.)

  ! call cpu_time(start)
  start = MPI_Wtime()
  ! looping over time
  do while ( t < sim_tmax )
    ! choose dt
    if (sim_fixDt) then
      dt = sim_dt
    else
      call cfl(dt)
    end if

    ! slow starter for handling discontinuities in IC
    dt_small = dt_init*2.**nstep
    if (dt_small < dt) dt = dt_small

    !check to see if there is a reason to stop
    if (sim_nlim .and. (nStep >= sim_nStep)) then
      exit
    elseif ( abs(t - sim_tmax) <= dt ) then
      dt = abs(t - sim_tmax)
    end if

    ! save current time to global variable,
    ! just for DMR
    sim_t = t

    ! do high order intp/recon/temporal
    call soln_numeric(dt)
    ! update solution
    call soln_update(dt)
    ! call BC on primitive vars
    call bc_apply(gr_V)

    ! write outputs every ioNfreq cycle or ioTfreq cycle
    ioCheckTime = sim_ioTfreq*real(ioTimeFreqCounter+1)
    if (t-dt < ioCheckTime .and. t>ioCheckTime) then
      if (bl_ID == 0) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
      end if
      ioCounter = ioCounter + 1
      ioTimeFreqCounter = ioTimeFreqCounter + 1
      finish = MPI_Wtime()
      call io_writeOutput(t, nStep, ioCounter, finish-start)
    endif

    if (sim_ioNfreq > 0) then
      if (mod(nStep, sim_ioNfreq) == 0) then
        if (bl_ID == 0) then
          write(*,*)''
          write(*,*)' Output no.',ioCounter+1, 'has been written      '
          write(*,*)'================================================='
          write(*,*)'   Steps      Time              dt               '
          write(*,*)'================================================='
          write(*,*)''
        end if
        ioCounter = ioCounter + 1
        finish = MPI_Wtime()
        call io_writeOutput(t, nStep, ioCounter, finish-start)
      endif
    endif

    ! update your time and step count
    t = t + dt
    nStep = nStep + 1

    ! print current infos to STDOUT
    if (bl_ID == 0) write(*,900) nstep, t, dt

    ! do we really need this?
    if (dt .le. 0.) then
      exit
    end if

  end do


  !! Let's write the final result before exiting
  if (bl_ID == 0) then
    write(*,*)''
    write(*,*)' Final output no.',ioCounter+1, 'has been written'
    write(*,*)'================================================='
    write(*,*)'        The final tmax has reached, bye!         '
    write(*,*)'================================================='
    write(*,*)''
  end if
  finish = MPI_Wtime()
  call io_writeOutput(t, nStep,ioCounter+1, finish-start)

  !! finalize and deallocate memories
  call grid_finalize()
  call block_finalize()


900 format(1x,i5,f16.8,1x,f16.8,f16.3)

end program SlugCode

