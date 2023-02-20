!     ==============================================================
!                 WRAPPER AROUND DISKEVOL_MODULE.F90
!     ==============================================================

!     --------------------------------------------------------------
!                             MAIN ROUTINE
!     --------------------------------------------------------------
program diskevol
  use natconst_module
  use diskevol_module
  use dust_module
  implicit none
  doubleprecision :: time
  doubleprecision :: rin,rout,sig0,plsig0,alpha0,alphadust,dt,dtmin,dtrel,timesave,dtold,tinfall
  doubleprecision :: tsav0,tend,rdisk0,dtdust,dusttime
  integer :: nr,ir,itsave=0,it,maxstep,nsave,itemp
  logical :: fex,dump,simstop
  character(len=100)               :: buffer, label
  integer                          :: pos
  integer                          :: ios = 0
  integer                          :: line = 0
  integer, parameter               :: fh = 1
  !
  ! Read the Rosseland mean opacity table for dust
  !
  inquire(file='rossmean_dust.inp',exist=fex)
  if(fex) then
     open(unit=1,file='rossmean_dust.inp')
     read(1,*) ross_ntemp
     allocate(ross_temp(ross_ntemp),ross_kappa(ross_ntemp,2))   ! (:,1) = dust, (:,2) = gas
     do itemp=1,ross_ntemp
        read(1,*) ross_temp(itemp),ross_kappa(itemp,1)
     enddo
     close(1)
  else
     stop 1
  endif
  !
  ! Set important parameters of the model
  !
  ! to be read from parameters file: default values set here
  !
  nr                           = 1000
  diskevol_gtd                 = 100.
  diskevol_flang               = 0.05
  diskevol_tempbg              = 10.d0
  diskevol_cloud_omega         = 5.d-15
  alpha0                       = 1.d-3
  alphadust                    = 1.d-3     !alphadust <= alpha0
  !
  ! other parameters:
  !
  diskevol_mstar               = 1.d-4*MS
  diskevol_rstar               = RS
  diskevol_tstar               = TS
  diskevol_mugas               = 2.3
  diskevol_sigma_min           = 1d-20
  diskevol_isw_qvisc           = 1
  diskevol_isw_qirrad          = 1
  diskevol_isw_growstar        = 1
  diskevol_isw_incmdisk        = 0
  diskevol_isw_lumbl           = 0
  diskevol_isw_selfirr         = 0
  diskevol_visheat_mode        = 0
  diskevol_stlumgrow           = 0.d0
  diskevol_gravinst            = 1
  diskevol_ginst_plaw          = 7.d0
  diskevol_ginst_qcrit         = 2.0
  diskevol_cloud_mass          = MS
  diskevol_cloud_ismooth       = 0
  diskevol_cloud_smoothparam1  = 5d-1
  diskevol_cloud_smoothparam2  = 100.*AU
  diskevol_cloud_idistr        = 1
  rin                          = 60.*RS
  rout                         = 10*pc
  sig0                         = 1.d-10
  plsig0                       = -1.0d0
  rdisk0                       = 100.*au
  maxstep                      = 100000000
  dtmin                        = 1*year
  dtrel                        = 0.00001
  tinfall                      = 1.3d5*year  ! starting from this time, outputs are dumped more often
  tsav0                        = 1d2*year  ! the first output time
  tend                         = 5d6*year  ! end time of the simulation
  nsave                        = 200
  !
  ! Read parameters values from the control file
  !
  open(fh,file='params.par',action='read')

      ! ios is negative if an end of record condition is encountered or if
      ! an endfile condition was detected.  It is positive if an error was
      ! detected.  ios is zero otherwise.

      do while (ios == 0)
         read(fh, '(A)', iostat=ios) buffer
         if (ios == 0) then
            line = line + 1

            pos = scan(buffer, ' 	')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label)
            case ('number_of_radial_zones')
               read(buffer, *, iostat=ios) nr
               print *, 'Number of radial zones: ', nr
            case ('gas_to_dust')
               read(buffer, *, iostat=ios) diskevol_gtd
               print *, 'Gas-to-dust ratio: ', diskevol_gtd
            case ('flaring_angle')
               read(buffer, *, iostat=ios) diskevol_flang
               print *, 'Flaring angle: ', diskevol_flang
            case('background_temperature')
               read(buffer, *, iostat=ios) diskevol_tempbg
               print *, 'Background temperature: ', diskevol_tempbg
            case('cloud_omega')
               read(buffer, *, iostat=ios) diskevol_cloud_omega
               print *, 'Cloud omega: ', diskevol_cloud_omega
            case('viscous_alpha')
               read(buffer, *, iostat=ios) alpha0
               print *, 'Viscous alpha: ', alpha0
            case('dust_alpha')
               read(buffer, *, iostat=ios) alphadust
               print *, 'Dust alpha: ', alphadust
            case default
               print *, 'Skipping invalid label at line', line
            end select
         end if
      end do
  !
  ! Initialize the diskevol_module
  !
  call diskevol_init(nr)
  !
  ! Set the grid
  !
  do ir=1,nr
     diskevol_grid_r(ir) = rin * (rout/rin)**((ir-1.d0)/(nr-1.d0))
  enddo
  !
  ! Set the initial surface density
  !
  do ir=1,nr
     diskevol_tprsigma(ir) = 2*pi*diskevol_grid_r(ir)*diskevol_sigma_min
     if(diskevol_grid_r(ir).le.rdisk0) then
        diskevol_tprsigma(ir) = diskevol_tprsigma(ir) +                      &
             2*pi*diskevol_grid_r(ir)*sig0*(diskevol_grid_r(ir)/au)**plsig0
     endif
  enddo
  !
  ! Init dust
  !
  call dust_init(diskevol_grid_r,diskevol_tprsigma,rdisk0,diskevol_gtd)
  !
  ! Set the alpha
  !
  do ir=1,nr
     diskevol_alpha(ir) = alpha0
  enddo
  !
  ! Reset some stuff related to data saving
  !
  time     = 0.d0
  dusttime = 0.d0
  itsave   = 0
  call diskevol_reset_savefiles()
  !
  ! First save data of the initial state (itsave=0)
  !
  call save_data(time,itsave)
  itsave   = 1
  timesave = tsav0
  dt = dtmin
  dtdust = 1.d-4*year
  !
  ! Now do the time loop
  !
  do it=1,maxstep
     !
     ! Reset flags
     !
     dump=.false.
     simstop=.false.
     !
     ! Check if we are about to cross the time-limit
     !
     if(time+dt.gt.tend) then
        dt = tend - time
        simstop = .true.
        dump    = .true.
     else
        !
        ! Check if we are about to cross the save time
        !
        if(time+dt.ge.timesave) then
           dt   = timesave - time
           dump = .true.
        endif
     endif
     !
     ! Message
     !
     write(*,*) 'Time step ',it,', Time = ',time/year,' year', ' dt = ',dt/year
     !
     ! Call the time step integrator
     !
     call diskevol_one_timestep(time,dt)
     !
     ! Update time
     !
     time = time + dt
     !
     ! Dust evolution: substepping to reach the gas time
     !
     do while (dusttime < time)
        if (dusttime + dtdust > time) dtdust = time - dusttime
        call dust_evol(dusttime,dtdust,diskevol_grid_r,diskevol_tprsigma,diskevol_temp,diskevol_nu, &
                    diskevol_mstar,diskevol_ml_sigdot,diskevol_vr,alphadust,diskevol_gtd)
        dusttime = dusttime + dtdust
     enddo
     !
     ! Save, if necessary
     !
     if(dump) then
        !
        ! Save the data
        !
        call save_data(time,itsave)
        !
        ! Compute the new timesave
        !
        itsave   = itsave + 1
        if (time < tinfall) then
           timesave = timesave * 10.
           if (timesave > tinfall) timesave = tinfall
        else
           timesave = tinfall * (tend/tinfall)**((itsave-1.d0)/(nsave-1.d0))
        endif
     endif
     !
     ! End simulation, if necessary
     !
     if(simstop) then
        goto 20
     endif
     !
     ! Set new gas timestep
     !
     if (time < 5.d5 * year) then
        dt = dtmin
     else
        dt = max(dtmin, time * dtrel)
     endif
     !
  enddo
  !
  ! Ended without reaching final time
  !
  write(*,*) 'Error: reached max nr of time steps before end of time.'
  stop
  !
  ! Ended correctly
  !
20 continue
  !
  ! Cleanup the diskevol_module
  !
  call diskevol_cleanup()
  call dust_cleanup()
  !
end program diskevol
