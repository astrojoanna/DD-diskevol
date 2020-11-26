!===============================================================================================
! DISK VISCOUS EVOLUTION MODEL
! by C.P. Dullemond
!
! This is a module for evolving a viscous disk model with infall, using an implicit
! integration scheme.
!
! See papers:
!
!   Dullemond, Natta & Testi (2006) ApJ 645, L69
!   Dullemond, Apai & Walch (2006) ApJ 640, L67
!   Hueso & Guillot (2005) A&A 442, 703
!
! This module is designed to be linked in together with other code.
!
! Before use you must call:
!
!   diskevol_init(nr)
!
! where nr is the number of radial grid points for the disk model. Next you must
! set up the grid by filling the array
!
!   diskevol_grid_r(:)
!
! Next you must define the initial condition for the gas surface density by filling
! the array
!
!   diskevol_tprsigma(:)
!
! The values in this array are not Sigma_gas, but 2*pi*r*Sigma_gas, where r is the
! radius (see diskevol_grid_r). This is because this is the conserved quantity.
! Next you must also define the value of the viscous alpha (the MRI one, i.e. not
! including grav inst):
!
!   diskevol_alpha(:)
!
! You should also fill the Rosseland mean opacity arrays. This requires you to
! allocate and fill the arrays ross_temp(:) and ross_kappa(:,:) (the second
! index of ross_kappa can be have dimension 1, i.e. ross_kappa(ntemp,1)).
!
! And finally you may want to modify the various scalar parameters and switches
! (which all have their default values which should be ok) to the values you wish.
!
! Now you you can make a time step of the surface density Sigma(r) by calling the
! subroutine:
!
!   diskevol_one_timestep(time,dt)
!
! This will update the diskevol_tprsigma(:) array to the new values. The reason
! why you must specify also "time" is because the infalling cloud model is
! dependent on the time after onset of collapse. The infall model is very
! simple: Shu combined with Ulrich.
!
! Note that at the same time the diskevol_tprsigma(:) is update, also various
! other arrays are updated, some of which might be useful for futher modeling:
!
!  diskevol_vr         The radial velocity of the gas in cm/s (interface centered!)
!  diskevol_temp       The gas midplane temperature in K
!  diskevol_nu         The viscosity in cm^2/s
!  diskevol_mdot       The accretion rate in g/s (interface centered!)
!  diskevol_teff       The disk effective surface temperature in K
!  diskevol_ml_sigdot  The mass loading rate in g/cm^2/s
!  diskevol_dcoef      The diffusion coefficient in cm^2/s
!  diskevol_alpnew     (if you have diskevol_isw_incmdisk=1): the actual alpha used
!  diskevol_mstar      (if you have diskevol_isw_growstar=1): The stellar mass
!  diskevol_photoev_sigdot  (if photoevaporation is active; currently not)
!
! After you finish using this module, you can free all allocated arrays by calling:
!
!   diskevol_cleanup()
!
!===============================================================================================
module diskevol_module
  use natconst_module
  use nrecip_module
  use dust_module
  !
  !-----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  !-----------------------------------------------------------------------------------
  !
  ! The global disk evolution parameters with their default values
  !
  doubleprecision :: diskevol_mstar      = 1d-4*MS    ! Default = Msun*1e-4
  doubleprecision :: diskevol_rstar      = RS         ! Default = Rsun
  doubleprecision :: diskevol_tstar      = TS         ! Default = Tsun
  doubleprecision :: diskevol_gtd        = 100.       ! The default gas-to-dust ratio (not important!)
  doubleprecision :: diskevol_flang      = 0.05       ! Fixed flaring angle
  doubleprecision :: diskevol_mugas      = 2.3        ! Fixed value for mugas
  doubleprecision :: diskevol_sigma_min  = 1d-20      ! Lower floor for surface density
  doubleprecision :: diskevol_tempbg     = 10.d0      ! Lower limit to the temperature of the disk
  !
  ! Switches for the different physics components of the disk evolution
  !
  integer ::         diskevol_isw_qvisc         = 1   ! Switch on the viscous heating
  integer ::         diskevol_isw_qirrad        = 1   ! Switch on the irradiative heating
  integer ::         diskevol_isw_growstar      = 1   ! Switch on the star mass growth
  integer ::         diskevol_isw_incmdisk      = 0   ! Switch on the inclusion of mass of disk in mstar
  integer ::         diskevol_isw_lumbl         = 0   ! Default: no irradiation by boundary layer
  integer ::         diskevol_isw_selfirr       = 0   ! Default: no self-irradiation
  doubleprecision :: diskevol_visheat_mode=0          ! How to compute T_c from T_eff for active accr.
  doubleprecision :: diskevol_stlumgrow  = 0.d0       ! If >0, then grow Lstar with Mstar
  integer ::         diskevol_gravinst   = 0          ! Include enhanced viscosity by gravitational instab?
  doubleprecision :: diskevol_ginst_plaw = 7.d0       ! The default p value for the grav-inst (1-Q)^p
  doubleprecision :: diskevol_ginst_qcrit= 2.0        ! The critical Toomre parameter
  !
  ! Defaults for infalling envelope
  !
  doubleprecision :: diskevol_cloud_mass        = MS                 ! Mass of infalling cloud
  doubleprecision :: diskevol_cloud_omega       = 2.d-14             ! Rotation rate of cloud (Ulrich model)
  integer ::         diskevol_cloud_ismooth     = 0                  ! =1 ---> Smooth the end of the infall phase
  doubleprecision :: diskevol_cloud_smoothparam1= 5d-1               ! Smoothing parameter
  doubleprecision :: diskevol_cloud_smoothparam2= 100.*1.4960000d+13 ! Smoothing parameter
  integer ::         diskevol_cloud_idistr      = 1                  ! Default: Hueso & Guillot way of distrib
                                                                     !   infalling matter onto the disk
  !
  ! Controls for the EUV photoevaporation
  !
  integer ::         diskevol_isw_euv_evap      = 0       ! Switch for photoevaporation by EUV photons
  doubleprecision :: euvp_phi                   = 0.d0    ! EUV luminosity (in number of photons)
  !
  ! Technical switches (methods)
  !
  integer ::         ml_implicit       = 1          ! Switch for mass loading in implicit way
  integer ::         ml_mode           = 1          ! How to include mass loading
  !
  ! The epsilon value for the Crank-Nicholson routine
  !    eps=0 -> exp, eps=1 -> imp, eps=0.5 -> CN
  !
  doubleprecision :: cranknich_eps     = 1.d0       ! Fully implicit
  !
  ! The arrays for the 2*pi*Sigma_gas(r) and alpha(r).
  ! The tprsigma is the current value of 2*pi*Sigma_gas(r) and will be evolved over time.
  ! An initial value has to be set by the user.
  ! The alpha has to be set by the user.
  !
  doubleprecision, allocatable :: diskevol_tprsigma(:)       ! INPUT AND OUTPUT: 2*pi*Sigma_gas(r) in g/cm^2
  doubleprecision, allocatable :: diskevol_alpha(:)          ! INPUT: alpha_mri(r)
  !
  ! The grid array
  !
  integer :: diskevol_grid_nr
  doubleprecision, allocatable :: diskevol_grid_r(:)
  !
  !-----------------------------------------------------------------------------------
  ! INTERNAL VARIABLES AND ARRAYS
  !-----------------------------------------------------------------------------------
  !
  ! Variables that are calculated on-the-fly
  !
  doubleprecision :: diskevol_cloud_asound      = 0.d0   ! Sound speed of infalling cloud (Shu model) (must be calculated)
  doubleprecision :: ml_mdot                    = 0.d0
  doubleprecision :: ml_mdotcap                 = 0.d0
  doubleprecision :: ml_rcentr                  = 0.d0
  !
  ! Data for the Rosseland mean opacity tables
  !
  doubleprecision, allocatable :: ross_temp(:)
  doubleprecision, allocatable :: ross_kappa(:,:)
  integer :: ross_ntemp
  !
  ! Data for the temperature solving
  !
  doubleprecision :: st_fact_qpl       = 0.d0
  doubleprecision :: st_fact_toomre    = 0.d0
  doubleprecision :: st_sigmakap       = 0.d0
  doubleprecision :: st_alpha_grav     = 0.d0
  doubleprecision :: st_alpha_mri      = 0.d0
  doubleprecision :: st_alpha          = 0.d0
  doubleprecision :: st_qplus          = 0.d0
  doubleprecision :: st_tempeff        = 0.d0
  !
  ! The disk variables (i.e. the actual model!)
  !
  doubleprecision, allocatable :: diskevol_vr(:)             ! v_r(r) = radial velo in cm/s,  Note: interface, so nr+1 elements
  doubleprecision, allocatable :: diskevol_temp(:)           ! T(r) = midplane temperature in K
  doubleprecision, allocatable :: diskevol_nu(:)             ! nu(r) = viscosity
  doubleprecision, allocatable :: diskevol_v0(:)             ! A kind of velocity (internal use)
  doubleprecision, allocatable :: diskevol_iv0(:)            ! A kind of velocity (internal use)
  doubleprecision, allocatable :: diskevol_mdot(:)           ! Mdot(r) = accr rate in g/s,  Note: interface, so nr+1 elements
  doubleprecision, allocatable :: diskevol_mgrav(:)          ! Mgrav(r) = the gravitational mass of the disk inside of r
  doubleprecision, allocatable :: diskevol_teff(:)           ! T_eff(r) = effective surf temp in K
  doubleprecision, allocatable :: diskevol_photoev_sigdot(:) ! The photoevaporation rate in g/cm^2/s
  doubleprecision, allocatable :: diskevol_ml_sigdot(:)      ! The mass loading rate of cloud onto disk in g/cm^2/s
  doubleprecision, allocatable :: diskevol_alpnew(:)         ! The alpha resulting from viscosity + grav inst
  doubleprecision, allocatable :: diskevol_dcoef(:)          ! The diffusion coefficient resulting from viscosity + grav inst
  !
  ! Helper arrays
  !
  integer :: diskevol_tmp_nrcap = 30                        ! Nr of grid points for infall inside of grid (rcap)
  doubleprecision, allocatable :: diskevol_tmp_a(:)
  doubleprecision, allocatable :: diskevol_tmp_b(:)
  doubleprecision, allocatable :: diskevol_tmp_c(:)
  doubleprecision, allocatable :: diskevol_tmp_q(:)
  doubleprecision, allocatable :: diskevol_tmp_null(:)
  doubleprecision, allocatable :: diskevol_tmp_onee(:)
  doubleprecision, allocatable :: diskevol_tmp_photev(:)
  doubleprecision, allocatable :: diskevol_tmp_infall(:)
  doubleprecision, allocatable :: diskevol_tmp_backup(:)
  doubleprecision, allocatable :: diskevol_tmp_src(:)
  doubleprecision, allocatable :: diskevol_tmp_rcap(:)
  doubleprecision, allocatable :: diskevol_tmp_rdum(:)

contains


  !--------------------------------------------------------------
  !                INITIALIZATION OF THIS MODULE
  ! Nr is the number of cells in the radial grid.
  ! Note that for cell interface quantities the arrays will be
  ! nr+1 elements.
  !--------------------------------------------------------------
  subroutine diskevol_init(nr)
    implicit none
    integer :: nr
    diskevol_grid_nr = nr
    allocate(diskevol_grid_r(nr)   )
    allocate(diskevol_tprsigma(nr) )
    allocate(diskevol_alpha(nr)    )
    allocate(diskevol_vr(nr+1)     )
    allocate(diskevol_temp(nr)     )
    allocate(diskevol_nu(nr)       )
    allocate(diskevol_v0(nr)       )
    allocate(diskevol_iv0(nr)      )
    allocate(diskevol_mdot(nr+1)   )
    allocate(diskevol_mgrav(nr)    )
    allocate(diskevol_teff(nr)     )
    allocate(diskevol_photoev_sigdot(nr))
    allocate(diskevol_ml_sigdot(nr))
    allocate(diskevol_alpnew(nr))
    allocate(diskevol_dcoef(nr))
    allocate(diskevol_tmp_a(nr)       )
    allocate(diskevol_tmp_b(nr)       )
    allocate(diskevol_tmp_c(nr)       )
    allocate(diskevol_tmp_q(nr)       )
    allocate(diskevol_tmp_null(nr)    )
    allocate(diskevol_tmp_onee(nr)    )
    allocate(diskevol_tmp_photev(nr)  )
    allocate(diskevol_tmp_infall(nr+1))
    allocate(diskevol_tmp_backup(nr)  )
    allocate(diskevol_tmp_src(nr)     )
    allocate(diskevol_tmp_rcap(diskevol_tmp_nrcap))
    allocate(diskevol_tmp_rdum(diskevol_tmp_nrcap))
  end subroutine diskevol_init


  !--------------------------------------------------------------
  !               CLEANUP OF THIS MODULE
  !--------------------------------------------------------------
  subroutine diskevol_cleanup()
    implicit none
    if(allocated(diskevol_grid_r       )) deallocate(diskevol_grid_r   )
    if(allocated(diskevol_tprsigma     )) deallocate(diskevol_tprsigma )
    if(allocated(diskevol_alpha        )) deallocate(diskevol_alpha    )
    if(allocated(diskevol_vr           )) deallocate(diskevol_vr       )
    if(allocated(diskevol_temp         )) deallocate(diskevol_temp     )
    if(allocated(diskevol_nu           )) deallocate(diskevol_nu       )
    if(allocated(diskevol_v0           )) deallocate(diskevol_v0       )
    if(allocated(diskevol_iv0          )) deallocate(diskevol_iv0      )
    if(allocated(diskevol_mdot         )) deallocate(diskevol_mdot     )
    if(allocated(diskevol_mgrav        )) deallocate(diskevol_mgrav    )
    if(allocated(diskevol_teff         )) deallocate(diskevol_teff     )
    if(allocated(diskevol_photoev_sigdot)) deallocate(diskevol_photoev_sigdot)
    if(allocated(diskevol_ml_sigdot    )) deallocate(diskevol_ml_sigdot)
    if(allocated(diskevol_alpnew       )) deallocate(diskevol_alpnew)
    if(allocated(diskevol_dcoef        )) deallocate(diskevol_dcoef)
    if(allocated(diskevol_tmp_a        )) deallocate(diskevol_tmp_a      )
    if(allocated(diskevol_tmp_b        )) deallocate(diskevol_tmp_b      )
    if(allocated(diskevol_tmp_c        )) deallocate(diskevol_tmp_c      )
    if(allocated(diskevol_tmp_q        )) deallocate(diskevol_tmp_q      )
    if(allocated(diskevol_tmp_null     )) deallocate(diskevol_tmp_null   )
    if(allocated(diskevol_tmp_onee     )) deallocate(diskevol_tmp_onee   )
    if(allocated(diskevol_tmp_photev   )) deallocate(diskevol_tmp_photev )
    if(allocated(diskevol_tmp_infall   )) deallocate(diskevol_tmp_infall )
    if(allocated(diskevol_tmp_backup   )) deallocate(diskevol_tmp_backup )
    if(allocated(diskevol_tmp_src      )) deallocate(diskevol_tmp_src    )
    if(allocated(diskevol_tmp_rcap     )) deallocate(diskevol_tmp_rcap   )
    if(allocated(diskevol_tmp_rdum     )) deallocate(diskevol_tmp_rdum   )
  end subroutine diskevol_cleanup



  !--------------------------------------------------------------
  !          FUNCTION: FIND ROSSELAND MEAN OPACITY
  !--------------------------------------------------------------
  function rossmean(temp,icmp)
    implicit none
    doubleprecision :: rossmean,temp
    integer :: icmp
    doubleprecision :: eps
    integer :: jlo
    !
    ! Check if the ross_* arrays are allocated
    !
    if(.not.allocated(ross_kappa)) then
       write(*,*) 'Error: Rosseland mean opacity arrays not set'
       stop 5
    endif
    !
    ! For the moment we only allow icmp=1
    !
    if(icmp.ne.1) stop 4
    !
    ! Find the position of temp in the rosseland mean opac table
    !
    call hunt(ross_temp,ross_ntemp,temp,jlo)
    !
    ! Check if out of bound
    ! If so, then return the closest kappa
    !
    if(jlo.le.0) then
       rossmean = ross_kappa(1,icmp)
       return
    endif
    if(jlo.ge.ross_ntemp) then
       rossmean = ross_kappa(ross_ntemp,icmp)
       return
    endif
    !
    ! Return the linear interpol of the kappa table
    !
    eps    = (temp-ross_temp(jlo))/(ross_temp(jlo+1)-ross_temp(jlo))
    if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 83
    rossmean  = (1.d0-eps)*ross_kappa(jlo,icmp)+eps*ross_kappa(jlo+1,icmp)
    !
    return
  end function rossmean

  !--------------------------------------------------------------
  !          HELPER FUNCTION FOR SOLVE TEMPERATURE
  !               FOR PURE ACCRETIONAL HEATING
  !--------------------------------------------------------------
  function helptemp(temp)
    implicit none
    doubleprecision :: helptemp,temp
    doubleprecision :: kappa,qmin,qtoomre
    !
    ! Find the rosseland mean opacity for this temperature
    !
    kappa    = rossmean(temp,1)
    !
    ! Find the Toomre parameter
    !
    if(st_fact_toomre.ne.0.d0) then
       qtoomre       = st_fact_toomre * sqrt(temp)
       if(qtoomre.lt.diskevol_ginst_qcrit) then
          if(qtoomre.lt.0.d0) stop 13
          st_alpha_grav = (1.d0-(qtoomre/diskevol_ginst_qcrit))**diskevol_ginst_plaw
       else
          st_alpha_grav = 0.d0
       endif
    else
       st_alpha_grav = 0.d0
    endif
    !
    ! Find the qplus for the accretion
    !
    st_qplus = st_fact_qpl * temp * ( st_alpha_mri + st_alpha_grav )
    !
    ! Find the qmin (the cooling rate)
    !
    if(diskevol_visheat_mode.eq.0) then
       !
       ! The usual purely radiative cooling through an optically thick disk
       !
       qmin     = 8.*ss*temp**4 / (3.*st_sigmakap*kappa)
       !
    elseif(diskevol_visheat_mode.eq.-1) then
       !
       ! Absolute extreme case: T_c = T_eff. So ignore the blanket-effect.
       ! NOTE: Just for testing.
       !
       qmin     = 2*ss*temp**4
    else
       stop 3354
    endif
    !
    ! Return the differences of the logs
    !
    helptemp = log(st_qplus) - log(qmin)
    return
  end function helptemp


  !--------------------------------------------------------------
  !              FUNCTION: COMPUTE SELF-RADIATION
  !
  ! A very simple self-irradiation recipe. We assume that the disk
  ! is optically thick.
  !--------------------------------------------------------------
  function find_temp_selfirr(ir)
    implicit none
    doubleprecision :: flux,find_temp_selfirr
    integer :: ir,ir2
    !
    flux = 0.d0
    if(ir.gt.2) then
       do ir2=2,ir-1
          flux = flux + ss *                                          &
               (diskevol_flang*(1.d0-(diskevol_grid_r(ir2)/diskevol_grid_r(ir))))**2 *  &
               0.5 * ( diskevol_teff(ir2-1)**4 + diskevol_teff(ir2)**4 ) *    &
               ( diskevol_grid_r(ir2)**2 - diskevol_grid_r(ir2-1)**2 ) /                &
               diskevol_grid_r(ir)**2
       enddo
    endif
    find_temp_selfirr = (flux/ss)**0.25
    !
  end function find_temp_selfirr


  !--------------------------------------------------------------
  !                 FUNCTION: SOLVE TEMPERATURE
  !
  ! This routine solves the local disk temperature including the
  ! heating via viscous dissipation as well as the heating by
  ! irradiation by the central star. As an important by-producy
  ! this routine also computes the viscosity coefficient alpha
  ! which may have been modified by gravitational instability.
  !--------------------------------------------------------------
  function diskevol_solvetemp(ir,r,tprsigma,temp0)
    implicit none
    doubleprecision :: diskevol_solvetemp,tprsigma,temp0,mu,r
    integer :: ir
    doubleprecision :: sigma,temp_visc,lumbl,rdx
    doubleprecision :: temp_irrad,temp_selfirr,temp_blirr
    doubleprecision :: eps,dummy,dm1,dm2
    !
    mu             = diskevol_mugas
    !
    ! First check if the ir indeed belongs to the r
    ! Just an internal consistency check.
    !
    if((ir.lt.1).or.(ir.gt.diskevol_grid_nr)) stop 8832
    if(ir.lt.diskevol_grid_nr) then
       if((r.lt.diskevol_grid_r(ir)).or.(r.ge.diskevol_grid_r(ir+1))) stop 8832
    endif
    !
    ! Compute the basic alpha at location r by linear interpolation
    ! Note: the real alpha (including grav inst effects) is computed
    ! on the fly during the temperature solving.
    !
    if(ir.ne.diskevol_grid_nr) then
       eps    = (r-diskevol_grid_r(ir))/(diskevol_grid_r(ir+1)-diskevol_grid_r(ir))
       st_alpha_mri = (1.d0-eps)*diskevol_alpha(ir)+eps*diskevol_alpha(ir+1)
    else
       st_alpha_mri = diskevol_alpha(diskevol_grid_nr)
    endif
    !
    ! First the temperature for accretion only
    !
    if(diskevol_isw_qvisc.ne.0) then
       !
       ! ...prepare the helper function for zbrent
       !
       sigma          = (tprsigma/(2*pi*r))
       st_fact_qpl    = (9./4.)*sigma*(kk/(mu*mp))*sqrt(GG*diskevol_mgrav(ir)/r**3)
       if(diskevol_gravinst.ne.0) then
          st_fact_toomre = sqrt(kk/(mu*mp))*sqrt(GG*diskevol_mstar/r**3)/(pi*GG*sigma)
       else
          st_fact_toomre = 0.d0
       endif
       !
       ! Compute the sigma of the dust (only used for the opacity!)
       !
       st_sigmakap     = sigma/diskevol_gtd
       !
       ! ...call zbrent to find the temperature
       !
       !
       if(helptemp(1d0).gt.0.d0) then
          dm1 = helptemp(1d0)
          dm2 = helptemp(1d10)
          if(dm2.gt.0.d0) then
             ! temp_visc = 1d10
             write(*,*) 'ERROR: No convergence in temperature at ir=',ir
             write(*,*) '  helptemp(1) = ',dm1
             write(*,*) '  helptemp(2) = ',dm2
             stop
          else
             temp_visc = zbrent(helptemp,1d0,1d10,1d-2)
          endif
       else
          temp_visc = 0.d0
       endif
       !
       ! One more call to the temperature routine, to make sure that
       ! we get the current values, both for the self-irradiation
       ! (the st_qplus) and for the gravitational instability alpha.
       ! (24-11-05)
       !
       if(temp_visc.gt.0.d0) then
          dummy      = helptemp(temp_visc)
       else
          st_qplus   = 0.d0
       endif
       !
    else
       temp_visc = 0.d0
       st_qplus  = 0.d0
    endif
    !
    ! Now the temperature from the irradiation
    !
    ! ...First the irradiation by the star
    !
    if(diskevol_isw_qirrad.ne.0) then
       !
       ! Check if we wish the star luminosity to grow with star mass
       !
       if(diskevol_stlumgrow.gt.0.d0) then
          if(ml_mode.eq.1) then
             rdx = (diskevol_mstar/diskevol_cloud_mass)**diskevol_stlumgrow
             if(rdx.gt.1.d0) rdx=1.d0
          else
             rdx = 1.d0
          endif
       else
          rdx = 1.d0
       endif
       !
       ! Compute the irradiation temperature
       !
       temp_irrad = diskevol_tstar*(0.5*diskevol_flang*rdx*(diskevol_rstar/r)**2)**0.25
    else
       temp_irrad = 0.d0
    endif
    !
    ! ...Then the irradiation by the disk itself (self-irradiation)
    !    NEWSTUFF 24-11-05
    !
    if(diskevol_isw_selfirr.ne.0) then
       temp_selfirr = find_temp_selfirr(ir)
    else
       temp_selfirr = 0.d0
    endif
    !
    ! ...Then the irradiation by the boundary layer of the disk
    !    (or in other words the `accretion shock' on the star)
    !    NEWSTUFF 24-11-05
    !
    if(diskevol_isw_lumbl.ne.0) then
       if(diskevol_isw_lumbl.eq.1) then
          !
          ! Use Calvet & Gullbring magnetospheric accretion model
          !
          lumbl      = (1.d0-diskevol_rstar/diskevol_grid_r(1)) *      &
                       ( GG * abs(diskevol_mdot(1)) * diskevol_mstar / &
                       diskevol_rstar )
          temp_blirr = (0.5*diskevol_flang*lumbl/(4*pi*r*r*ss))**0.25
       else
          !
          ! Simple boundary layer (not yet done)
          !
          stop 8264
       endif
    else
       temp_blirr = 0.d0
    endif
    !
    ! Return the 'sum' of the two
    ! BUGFIX 13-07-05: Add the background temperature as well
    !
    diskevol_solvetemp = (temp_visc**4+temp_irrad**4+diskevol_tempbg**4+   &
                          temp_blirr**4+temp_selfirr**4)**0.25
    !
    ! Compute (for use in self-irradiation) the effective temperature
    ! assuming a blackbody disk.
    ! NOTE: I omit the background and the self-irradiation here
    ! NEWSTUFF 24-11-05
    !
    st_tempeff     = ( (st_qplus/(2*ss)) + temp_irrad**4 )**0.25
    !
    ! If gravitational instability mode active, then compute new alpha
    !
    if(diskevol_gravinst.ne.0) then
       dummy      = helptemp(diskevol_solvetemp)
       st_alpha   = st_alpha_mri + st_alpha_grav
    else
       st_alpha   = st_alpha_mri
    endif
    !
    return
  end function diskevol_solvetemp


  !---------------------------------------------------------------
  !                FUNCTION: SIGMA * NU * SQRT(R)
  !
  ! This function contains all the microphysics needed for the
  ! Shakura-Sunyaev disk. It computes the Sigma * nu * sqrt(R)
  ! value, given the radius R (r) and the 2*pi*Sigma*R (tprsigma). The
  ! arguments temp and nu are meant for returning the values
  ! of temperature and viscosity, in case these are needed.
  !    A crucial property of this function is that the resulting
  ! value only depends on r and tprsigma! This is because the
  ! function will be called by the functions that compute the
  ! double-log derivative coefficients for the linearized
  ! equations. The 'ir' argument is only there to make it
  ! easier to locate where we are on the grid.
  !
  ! ARGUMENTS:
  !    ir                Index
  !    r                 Radius in cm
  !    tprsigma          Sigma*2*pi*R in g/cm
  !    diskevol_temp         Initial guess of temperature (usually
  !                      this is not necessary, but may become
  !                      important for disk instabilities due to
  !                      the multiple temperature solutions).
  !
  ! RETURNS:
  !    sigma_nu_sqr      The value of Sigma*nu*sqrt(R)
  !    temp              The midplane temperature of the disk
  !    nu                The viscosity coefficient
  !    diskevol_alpnew     The new value of the alpha with grav inst
  !---------------------------------------------------------------
  function sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
    implicit none
    integer :: ir
    doubleprecision :: sigma_nu_sqrtr,r,tprsigma,temp,nu
    doubleprecision :: mu
    !
    mu    = diskevol_mugas
    !
    ! Compute the midplane temperature, and as a by-product compute
    ! the modified alpha (modified by gravitational instability)
    !
    temp  = diskevol_solvetemp(ir,r,tprsigma,diskevol_temp(ir))
    !
    ! The viscosity coefficient, using the alpha value returned by
    ! the diskevol_solvetemp() routine (=st_alpha)
    ! NOTE: The constant = k_B / (m_p*sqrt(G))
    !
    nu    = 3.1957972d11 * st_alpha * r**1.5 * temp /     &
                    ( sqrt(diskevol_mgrav(ir)) * mu )
    !
    ! Do a check
    !
    if(nu.le.0.d0) then
       write(*,*) 'ERROR: nu zero or negative...'
       stop
    endif
    if(number_invalid(nu).ne.0) then
       write(*,*) 'ERROR: Invalid nu'
       write(*,*) '     ir    = ',ir
       write(*,*) '     r     = ',r
       write(*,*) '     nu    = ',nu
       write(*,*) '     temp  = ',temp
       write(*,*) '     alpha = ',st_alpha
       write(*,*) '     mgrav = ',diskevol_mgrav(ir)
       write(*,*) '     mstar = ',diskevol_mstar
       stop
    endif
    !
    ! Return the final value
    !
    sigma_nu_sqrtr = (tprsigma/(2*pi*r)) * nu * sqrt(r)
    !
    return
  end function sigma_nu_sqrtr



  !--------------------------------------------------------------
  !      DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(SIGMA*R)
  !--------------------------------------------------------------
  function deriv_signusr_sigr(ir,r,tprsigma)
    implicit none
    integer :: ir
    doubleprecision :: deriv_signusr_sigr,r,tprsigma
    !
    doubleprecision :: eps,y1,y2,dum1,dum2
    parameter(eps=1.d-4)
    !
    y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2)
    y2 = sigma_nu_sqrtr(ir,r,tprsigma*(1.d0+eps),dum1,dum2)
    deriv_signusr_sigr = 2.d0*(y2-y1)/(eps*(y2+y1))
    !
    return
  end function deriv_signusr_sigr


  !--------------------------------------------------------------
  !        DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(R)
  !--------------------------------------------------------------
  function deriv_signusr_r(ir,r,tprsigma)
    implicit none
    integer :: ir
    doubleprecision :: deriv_signusr_r,r,tprsigma
    !
    doubleprecision :: eps,y1,y2,dum1,dum2
    parameter(eps=1.d-4)
    !
    y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2)
    y2 = sigma_nu_sqrtr(ir,r*(1.d0+eps),tprsigma,dum1,dum2)
    deriv_signusr_r = 2.d0*(y2-y1)/(eps*(y2+y1))
    !
    return
  end function deriv_signusr_r


  !--------------------------------------------------------------
  !        COMPUTE RADIAL VELOCITY AND DIFFUSION COEFFICIENT
  !
  ! This routine computes the radial velocity vr according to the
  ! physics of Shakura-Sunyaev disks. It also computes the
  ! effective radial velocity v0 and the diffusion coefficient
  ! D that you get when you write the continuity equation out
  ! into a diffusion-advection type equation. Because the pure
  ! advection form with vr is numerically unstable, we use the
  ! diffusion-advection type equation for the model.
  !--------------------------------------------------------------
  subroutine compute_vr_v0_difcoef()
    implicit none
    doubleprecision :: sigma,mu,tillum,taues,temp_illum,r,tprsigma
    doubleprecision :: signusr,pd_r,pd_sigr,nu,temp
    integer :: ir
    parameter(tillum=3.481278d9)
    !
    mu    = diskevol_mugas
    !
    ! First check something
    !
    if(diskevol_sigma_min.eq.0.d0) then
       write(*,*) 'Error: Must define diskevol_sigma_min'
       write(*,*) '       in order to keep log derivs intact'
       stop 13
    endif
    !
    ! For the numerical model we use v0 and D. These are found as
    ! partial derivatives of Sigma * nu * sqrt(R) to (Sigma*R)_{const R}
    ! and to (R)_{const Sigma*R}. This is simply using the chain rule.
    !
    !                 / dlg(Sigma*nu*sqrt(R)) \           nu
    !       v_0 = -3  | --------------------- |           --
    !                 \        dlg(R)         /{Sigma*R}  R
    !
    !                 / dlg(Sigma*nu*sqrt(R)) \                        .
    !       D   =  3  | --------------------- |     nu
    !                 \      dlg(Sigma*R)     /{R}
    !
    do ir=1,diskevol_grid_nr
       r                 = diskevol_grid_r(ir)
       tprsigma          = diskevol_tprsigma(ir) + 2*pi*r*diskevol_sigma_min
       signusr           = sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
       if(temp.le.0.d0) then
          write(*,*) 'ERROR: Temperature <=0 found!'
          stop
       endif
       diskevol_temp(ir)     = temp
       diskevol_nu(ir)       = nu
       diskevol_teff(ir)     = st_tempeff
       diskevol_alpnew(ir)   = st_alpha
       pd_r                  = deriv_signusr_r(ir,r,tprsigma)
       pd_sigr               = deriv_signusr_sigr(ir,r,tprsigma)
       diskevol_v0(ir)       = -3.d0*pd_r*nu/r
       diskevol_dcoef(ir)    = 3.d0*pd_sigr*nu
    enddo
    !
    ! Compute the intervace values of of diskevol_v0
    !
    ! **** NOTE: **** 27-08-04   [20-02-05: I put this note here
    !                             later. Hope I interpreted it OK.]
    !      I changed the indexing of vint, so that it is
    !      consistent with the indexing of 'flux' in the
    !      aposteriori_flux_and_source() routine!!!
    !
    do ir=2,diskevol_grid_nr
       diskevol_iv0(ir) = 0.5d0 * ( diskevol_v0(ir) + diskevol_v0(ir-1) )
    enddo
    diskevol_iv0(1) = 0.d0
    !
    !     Done
    !
  end subroutine compute_vr_v0_difcoef




  ! --------------------------------------------------------------
  !       THE CRANK-NICHOLSON SCHEME FOR DIFFUSION-ADVECTION
  !
  !     Perform one time step for the following PDE:
  !
  !        du    d /     \    d /           d  / u  \ \
  !        -- + -- | u v | - -- | h(x) D(x) -- |----| | = K + L u
  !        dt   dx \     /   dx \           dx \h(x)/ /
  !
  !     with boundary conditions
  !
  !         du/h |            |
  !       p ---- |      + q u |       = r
  !          dx  |x=xbc       |x=xbc
  !
  !     ARGUMENTS:
  !       n             Nr of grid points
  !       x             The grid
  !       u             The current values of u(x)
  !       vint          The values for v(x) [NOTE: vint defined at
  !                                                cell interface!]
  !       d             The values for D(x)
  !       h             The values for h(x)
  !       k             The values for K(x)
  !       l             The values for L(x)
  !       dt            The time step
  !       pl,pr         The p value for the lhs/rhs of the BC equation
  !       ql,qr         The q value for the lhs/rhs of the BC equation
  !       rl,rr         The r value for the lhs/rhs of the BC equation
  !       a,b,c,q       Dummy arrays used by the cranknich() routine
  !
  !     NOTE: This is the new version, in which the bugs are fixed,
  !           and in which a K and L term are added on the right-
  !           hand-side (01-11-03).
  !
  !     NOTE: An even newer version (changed name to cranknich_general)
  !           which also allows `relative diffusion': the h(x) factor
  !           in the equation (27-01-04).
  !
  !     *** IMPORTANT NOTE: ***
  !           This is an adapted version of the cranknich_general
  !           routine in which the velocity v is given at the cell
  !           boundaries!!! Specially designed for the radial mixing
  !           problem...
  !           **** NOTE: **** 27-08-04
  !                I changed the indexing of vint, so that it is
  !                consistent with the indexing of 'flux' in the
  !                aposteriori_flux_and_source() routine!!!
  !
  !--------------------------------------------------------------
  subroutine cranknich_general(n,x,u,vint,d,h,k,l,dt,pl,pr,ql,qr,rl,rr,a,b,c,q)
    implicit none
    integer :: n
    doubleprecision :: pl,pr,ql,qr,rl,rr,dt
    doubleprecision :: x(n),u(n),vint(n+1),d(n),a(n),b(n),c(n),q(n)
    doubleprecision :: k(n),l(n),h(n)
    integer :: i
    doubleprecision :: eps,eps1,dum,dpl,dmn,vpl,vmn
    !
    ! The epsilon
    !
    eps  = cranknich_eps
    eps1 = 1.d0-eps
    !
    if(vint(1).ne.0.d0) then
       write(*,*) 'ERROR: This is a special version of '
       write(*,*) '   cranknich_general: v defined at interface!'
       write(*,*) '   Found vint(1).ne.0... Should not be!'
       stop
    endif
    !
    ! The main body of the equation
    !
    do i=2,n-1
       dum  = 2.0*dt/(x(i+1)-x(i-1))
       dpl  = dum*0.25*(h(i+1)+h(i))*(d(i+1)+d(i))/(x(i+1)-x(i))
       dmn  = dum*0.25*(h(i)+h(i-1))*(d(i)+d(i-1))/(x(i)-x(i-1))
       vpl  = dt*vint(i+1)/(x(i+1)-x(i-1))
       vmn  = dt*vint(i)/(x(i+1)-x(i-1))
       b(i) = 1.d0+eps*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))
       a(i) = -eps*dmn/h(i-1)-eps*vmn
       c(i) = -eps*dpl/h(i+1)+eps*vpl
       q(i) = u(i)+k(i)*dt                                      &
                 -eps1*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))*u(i)    &
                 +eps1*dmn*u(i-1)/h(i-1)+eps1*dpl*u(i+1)/h(i+1) &
                 +eps1*vmn*u(i-1)-eps1*vpl*u(i+1)
    enddo
    !
    ! The left (i=1) boundary condition
    !
    b(1) = ql-pl/(h(1)*(x(2)-x(1)))
    c(1) = pl/(h(2)*(x(2)-x(1)))
    q(1) = rl
    !
    ! The right (i=n) boundary condition
    !
    b(n) = qr+pr/(h(n)*(x(n)-x(n-1)))
    a(n) = -pr/(h(n-1)*(x(n)-x(n-1)))
    q(n) = rr
    !
    ! Solve the implicit differencing system
    !
    call tridag(a,b,c,q,u,n)
    !
    !     Check whether the solution is indeed a solution
    !
    !      do i=2,n-1
    !         tmp(i) = a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1)
    !      enddo
    !      tmp(1) = b(1)*u(1) + c(1)*u(2)
    !      tmp(n) = a(n)*u(n-1) + b(n)*u(n)
    !      do i=1,n
    !         error = abs(tmp(i)-q(i))/(abs(tmp(i))+abs(q(i))+   &
    !                    abs(b(i)*uorig(i)))
    !         if(error.gt.1d-3) then
    !            write(*,*) 'ERROR in matrix inversion. Error=',error
    !            write(*,*) '   at i=',i,', tmp,q,u=',tmp(i),q(i),u(i), &
    !                 uorig(i)
    !            stop
    !         endif
    !      enddo
    !
  end subroutine cranknich_general


  ! --------------------------------------------------------------
  !        A-POSTERIORI COMPUTATION OF FLUXES AND SOURCES
  !
  ! This routine computes to machine precision what the fluxes
  ! at the cell interfaces were over the last time step and what
  ! the source terms were. For the cell interfaces between two
  ! cell centers this is straight-forward. For the cell interfaces
  ! at the edges this is not straightforward, but the fluxes can
  ! be computed indirectly. For this to work one has to specify
  ! where the boundaries of the outermost cells are. This can be
  ! directly at the location of the grid point (possibility 1) or
  ! a linear extrapolation (possibility 2):
  !
  !      *---*-----*-----*
  !      | |    |     |      (possibility 1)
  !
  !      *---*-----*-----*
  !    |   |    |     |      (possibility 2)
  !
  ! The resulting flux and source are such that one can verify:
  !
  !                        2*dt     /             \
  !   unew_i - uold_i = ----------- | F_i+1 - F_i | + dt * Src_i
  !                     x_i+1-x_i-1 \             /
  !
  ! *** IMPORTANT NOTE: ***
  !       This is an adapted version of the cranknich_general
  !       routine in which the velocity v is given at the cell
  !       boundaries!!! Specially designed for the radial mixing
  !       problem...
  !       **** NOTE: **** 27-08-04
  !            I changed the indexing of vint, so that it is
  !            consistent with the indexing of 'flux' in the
  !            aposteriori_flux_and_source() routine!!!
  !
  !--------------------------------------------------------------
  subroutine aposteriori_flux_and_source(n,x,uold,unew,vint,d,h,k,l,dt,ibnd,flux,src)
    implicit none
    integer :: n,ibnd
    doubleprecision :: dt,flux(n+1),src(n)
    doubleprecision :: x(n),uold(n),unew(n)
    doubleprecision :: vint(n+1),d(n),k(n),l(n),h(n)
    integer :: i
    doubleprecision :: eps,eps1,flold,flnew,v_av,d_av,uold_av,unew_av,dx
    doubleprecision :: dxl,dxr,h_av
    !
    eps  = cranknich_eps
    eps1 = 1.d0-eps
    !
    if(vint(1).ne.0.d0) then
       write(*,*) 'ERROR: This is a special version of '
       write(*,*) '   cranknich_general: v defined at interface!'
       write(*,*) '   Found vint(1).ne.0... Should not be!'
       stop
    endif
    !
    ! Compute the fluxes for i=2,n
    !
    do i=2,n
       v_av    = vint(i)
       d_av    = 0.5*(d(i-1)+d(i))
       h_av    = 0.5*(h(i-1)+h(i))
       uold_av = 0.5*(uold(i-1)+uold(i))
       unew_av = 0.5*(unew(i-1)+unew(i))
       dx      = x(i)-x(i-1)
       flold   = v_av*uold_av - h_av*d_av*(uold(i)/h(i)-uold(i-1)/h(i-1))/dx
       flnew   = v_av*unew_av - h_av*d_av*(unew(i)/h(i)-unew(i-1)/h(i-1))/dx
       flux(i) = eps1*flold + eps*flnew
       src(i)  = k(i) + l(i) * ( eps1*uold(i) + eps*unew(i) )
    enddo
    !
    ! For i=1 and i=n+1 (the boundaries) we know what the source
    ! should be, and we know what the boundary conditions have
    ! produced (unew(1) and unew(n+1)). So we can compute implicitly
    ! what the flux must have been.
    !
    ! First compute the source terms:
    !
    src(1)   = k(1) + l(1) *  ( eps1*uold(1) + eps*unew(1) )
    !
    ! Now specify the dx of the boundary cells according to possibility
    ! 1 or 2
    !
    if(ibnd.eq.1) then
       dxl = 0.5*(x(2)-x(1))
       dxr = 0.5*(x(n)-x(n-1))
    elseif(ibnd.eq.2) then
       dxl = x(2)-x(1)
       dxr = x(n)-x(n-1)
    else
       write(*,*) 'ERROR in aposteriori_flux_and_source():'
       write(*,*) '   Dont know ibnd=',ibnd
       stop
    endif
    !
    ! Now find the flux at the left and right boundaries
    !
    flux(1) = (dxl/dt)*(unew(1)-uold(1))-dxl*src(1)+flux(2)
    flux(n+1) = -(dxr/dt)*(unew(n)-uold(n))+dxr*src(n)+flux(n)
    !
    ! Done...
    !
  end subroutine aposteriori_flux_and_source


  !--------------------------------------------------------------
  !                     INTEGRATION ROUTINE
  !--------------------------------------------------------------
  function integrate(n,x,f)
    implicit none
    integer :: n,i
    doubleprecision :: x(n),f(n),integrate,int
    int=0.d0
    do i=2,n
       int=int+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
    enddo
    if(x(n).gt.x(1)) then
       integrate = int
    else
       integrate = -int
    endif
    return
  end function integrate


  !--------------------------------------------------------------
  !              UPDATE THE SIGMA FROM DIFFUSION
  !
  ! This is the main routine for the viscous evolution of the
  ! disk.
  !--------------------------------------------------------------
  subroutine update_sigma_diff(dt)
    implicit none
    doubleprecision :: dt,r,tprsigma,signusr,temp,nu
    doubleprecision :: pl,pr,ql,qr,rl,rr
    integer :: ir,iflag
    !
    ! Set null to null
    !
    do ir=1,diskevol_grid_nr
       diskevol_tmp_null(ir) = 0.d0
       diskevol_tmp_onee(ir) = 1.d0
    enddo
    !
    ! Backup current tprsigma
    !
    do ir=1,diskevol_grid_nr
       diskevol_tmp_backup(ir) = diskevol_tprsigma(ir)
    enddo
    !
    ! Set boundary conditions at inner edge
    !
    pl      = 0.d0
    ql      = 2.d0
    rl      = 2.d0*diskevol_sigma_min*2*pi*diskevol_grid_r(1)
    !
    ! Set boundary conditions at outer edge
    !
    r       = diskevol_grid_r(diskevol_grid_nr)
    tprsigma= diskevol_tprsigma(diskevol_grid_nr)
    signusr = sigma_nu_sqrtr(diskevol_grid_nr,r,tprsigma,temp,nu)
    pr      = deriv_signusr_sigr(diskevol_grid_nr,r,tprsigma) * nu
    qr      = deriv_signusr_r(diskevol_grid_nr,r,tprsigma) * nu / r
    rr      = 0.d0
    !
    ! First reset source term
    !
    do ir=1,diskevol_grid_nr
       diskevol_tmp_infall(ir) = 0.d0
    enddo
    !
    ! Include mass loading elements
    !
    if((ml_mode.ne.0).and.(ml_implicit.eq.1)) then
       do ir=1,diskevol_grid_nr
          diskevol_tmp_infall(ir) = diskevol_tmp_infall(ir) + 2*pi*diskevol_grid_r(ir)*diskevol_ml_sigdot(ir)
       enddo
    endif
    !
    ! Include photoevap sink
    !
    if(diskevol_isw_euv_evap.ne.0) then
       do ir=1,diskevol_grid_nr
          diskevol_tmp_photev(ir) = -2*pi*diskevol_grid_r(ir)*diskevol_photoev_sigdot(ir)/diskevol_tprsigma(ir)
       enddo
    else
       do ir=1,diskevol_grid_nr
          diskevol_tmp_photev(ir) = 0.d0
       enddo
    endif
    !
    ! Perform Crank-Nicholson step and the sources
    !
    call cranknich_general(diskevol_grid_nr,diskevol_grid_r,diskevol_tprsigma,   &
         diskevol_iv0,diskevol_dcoef,diskevol_tmp_onee,                          &
         diskevol_tmp_infall,diskevol_tmp_photev,dt,pl,pr,ql,qr,rl,rr,           &
         diskevol_tmp_a,diskevol_tmp_b,diskevol_tmp_c,diskevol_tmp_q)
    !
    ! First check if the results are reasonable
    !
    do ir=1,diskevol_grid_nr
       if(number_invalid(diskevol_tprsigma(ir)).ne.0) then
          write(*,*) 'ERROR: Invalid diskevol_tprsigma()'
          write(*,*) '     ir    = ',ir
          write(*,*) '     nu    = ',nu
          write(*,*) '     alpha = ',diskevol_alpha(ir)
          write(*,*) '     T     = ',diskevol_temp(ir)
          write(*,*) '     tprs  = ',diskevol_tprsigma(ir)
          write(*,*) pl,ql,rl,pr,qr,rr
          stop
       endif
    enddo
    !
    ! Check for positiveness
    !
    iflag=0
    do ir=1,diskevol_grid_nr
       if(diskevol_tprsigma(ir).lt.0.d0) then
          iflag=1
       endif
    enddo
    if(iflag.ne.0) then
       write(*,*) 'Negative Sigma detected'
       stop
    endif
    !
    ! Implement floor value
    !
    do ir=1,diskevol_grid_nr
       if(diskevol_tprsigma(ir).lt.2.d0*pi*diskevol_grid_r(ir)*diskevol_sigma_min) then
          diskevol_tprsigma(ir) = 2.d0*pi*diskevol_grid_r(ir)*diskevol_sigma_min
       endif
    enddo
    !
    ! Now a-posteriori exact computation of mass flux at interfaces
    !
    call aposteriori_flux_and_source(diskevol_grid_nr,diskevol_grid_r,             &
                diskevol_tmp_backup,diskevol_tprsigma,diskevol_iv0,diskevol_dcoef, &
                diskevol_tmp_onee,diskevol_tmp_null,diskevol_tmp_null,dt,2,        &
                diskevol_mdot,diskevol_tmp_src)
    !
    ! Compute the exact v_R from flux/(2*pi*r*sigma)
    !
    do ir=2,diskevol_grid_nr
       diskevol_vr(ir) = diskevol_mdot(ir) / ( 0.5*(diskevol_tprsigma(ir)+     &
                                         diskevol_tprsigma(ir-1)) )
    enddo
    diskevol_vr(1) = 0.d0
    diskevol_vr(diskevol_grid_nr+1) = 0.d0
    !
    ! Done
    !
  end subroutine update_sigma_diff


  !--------------------------------------------------------------
  !           CALCULATE THE INFALL ONTO THE DISK SURFACE
  !
  ! Use Shu model coupled to Ulrich model to calculate the
  ! accretion rate onto the disk surface.
  !--------------------------------------------------------------
  subroutine shuulrich_infall(nr,r,mcloud,cs,omega,time,psi,mdot,  &
       mdotcap,rcentr,ismooth,smoothparam1,smoothparam2,idistr)
    implicit none
    integer :: nr,icen,ir,ismooth,idistr
    doubleprecision :: r(nr),psi(nr),dum(nr)
    doubleprecision :: fact,dm
    doubleprecision :: smoothparam1,smoothparam2
    doubleprecision :: mcloud,cs,omega,time,mcentr,psicap,eps,jdot
    doubleprecision, parameter :: m0 = 0.975
    doubleprecision :: mdot,rcloud,tcloud,tcloudm,mfree,mcoll
    doubleprecision :: rcoll,rcoll0,mdotcap,rcentr,mmu0,rhocmueq,vz
    integer :: it
    !
    ! Switch between Shu collapse and the externally specified collapse
    !
    if(ml_mode.eq.1) then
       !
       ! Shu model
       !
       ! Calculate the global numbers for the Shu infall model with
       ! rotation.
       !
       mdot   = m0*cs**3/GG
       rcloud = GG * mcloud / ( 2*cs**2 )
       tcloud = rcloud/cs
       tcloudm= tcloud*(2./m0)
       if(time.lt.tcloudm) then
          mcentr = (m0/2.)*mcloud*(time/tcloud)
          if(time.lt.tcloud) then
             mcoll  = mcloud*(time/tcloud)
          else
             mcoll  = mcloud
          endif
       else
          mcentr = mcloud
          mcoll  = mcloud
       endif
       mfree  = mcoll-mcentr
       rcoll  = cs*time
       rcoll0 = rcoll * (m0/2.d0)
       if(time.gt.tcloud) then
          rcoll = rcloud
       endif
       if(time.gt.tcloudm) then
          rcoll0 = rcloud
       endif
       rcentr = omega**2*rcoll0**4/(GG*(mcentr+1d0))
       !
       ! Now include, if asked for, the smoothing-off of the infall
       !
       if(ismooth.eq.0) then
          !
          ! Original Shu type model: an abrupt end to the infall
          !
          if(time.gt.tcloudm) then
             mdot   = 0.d0
          endif
       elseif(ismooth.eq.1) then
          !
          ! A linearly smoothed-off version
          !
          if(smoothparam1.gt.1.d0) stop 7290
          if(smoothparam1.lt.0.d0) stop 7291
          if(time/tcloudm.gt.1.d0-smoothparam1) then
             dm = (time/tcloudm-1.d0+smoothparam1)/(2*smoothparam1)
             if(dm.lt.0.d0) stop 7292
             mdot = (1.d0-dm)*mdot
             if(mdot.lt.0.d0) mdot=0.d0
          endif
       elseif(ismooth.eq.11) then
          !
          ! A constant tail of infall
          !
          if(time.gt.tcloudm) then
             mdot   = smoothparam1*mdot
             rcentr = smoothparam2
          endif
       elseif(ismooth.eq.21) then
          !
          ! A constant tail of infall
          !
          if(time.gt.tcloudm) then
             mdot   = smoothparam1*mdot*(tcloudm/time)
             rcentr = smoothparam2
          endif
       endif
    else
       stop 3097
    endif
    !
    ! Compute the accretion rate onto the disk
    !
    if(mdot.gt.0.d0) then
       !
       ! Yes, the mass-loading is active at present...
       !
       icen = 0
       do ir=1,nr
          if(r(ir).lt.rcentr) icen=ir
       enddo
       psi(1) = 0.d0
       if(idistr.eq.0) then
          !
          ! Distribute the matter over the disk in the Ulrich way
          !
          do ir=1,icen-1
             mmu0      = sqrt(1-r(ir)/rcentr)
             rhocmueq  = mdot*(r(ir)/rcentr)/((4*pi*sqrt(GG*mcentr*r(ir)**3))*(2*mmu0**2))
             vz        = sqrt(GG*mcentr/r(ir))*mmu0
             psi(ir)   = 2*rhocmueq*vz
          enddo
       elseif(idistr.eq.1) then
          !
          ! Distribute the matter over the disk in such a way that
          ! there is no friction of the infalling matter with the disk
          ! (see Hueso & Guillot 2005)
          !
          do ir=1,icen-1
             psi(ir)   = (mdot/(8*pi*r(ir)*rcentr)) * ((r(ir)/rcentr)* &
                         (1.d0-sqrt(r(ir)/rcentr)))**(-0.5)
          enddo
       else
          stop 8109
       endif
       do ir=max(icen,1),nr
          psi(ir) = 0.d0
       enddo
       !
       ! Compute rate of accretion inward of inner R coordinate
       ! This mass is assumed to be directly loaded onto the star
       !
       if(rcentr.le.r(1)) then
          !
          ! All matter falls within R_in, which is by definition the
          ! capture radius of the central star.
          !
          if(icen.ne.0) stop 19009
          mdotcap = mdot
       else
          if(idistr.eq.0) then
             !
             ! Ulrich way
             !
             do ir=1,diskevol_tmp_nrcap
                diskevol_tmp_rcap(ir) = r(1)*ir/(1.d0*diskevol_tmp_nrcap)
                mmu0                  = sqrt(1-diskevol_tmp_rcap(ir)/rcentr)
                rhocmueq   = mdot*(diskevol_tmp_rcap(ir)/rcentr)/             &
                           ((4*pi*sqrt(GG*mcentr*diskevol_tmp_rcap(ir)**3))*  &
                           (2*mmu0**2))
                vz         = sqrt(GG*mcentr/diskevol_tmp_rcap(ir))*mmu0
                psicap     = 2*rhocmueq*vz
                diskevol_tmp_rdum(ir) = 2*pi*diskevol_tmp_rcap(ir)*psicap
             enddo
          elseif(idistr.eq.1) then
             !
             ! Hueso & Guillot way
             !
             do ir=1,diskevol_tmp_nrcap
                diskevol_tmp_rcap(ir) = r(1)*ir/(1.d0*diskevol_tmp_nrcap)
                psicap     = (mdot/(8*pi*diskevol_tmp_rcap(ir)*rcentr)) *       &
                              ((diskevol_tmp_rcap(ir)/rcentr)*                  &
                              (1.d0-sqrt(diskevol_tmp_rcap(ir)/rcentr)))**(-0.5)
                diskevol_tmp_rdum(ir) = 2*pi*diskevol_tmp_rcap(ir)*psicap
             enddo
          else
             stop 7204
          endif
          mdotcap = integrate(diskevol_tmp_nrcap,diskevol_tmp_rcap,diskevol_tmp_rdum)
          if(mdotcap.gt.mdot) mdotcap = mdot !stop 90909
          if(mdotcap.lt.0.d0) then
             write(*,*) (diskevol_tmp_rcap(ir),ir=1,diskevol_tmp_nrcap)
             write(*,*) (diskevol_tmp_rdum(ir),ir=1,diskevol_tmp_nrcap)
             stop
          endif
       endif
       !
       ! Normalize this to make sure the mass loading is perfect
       !
       do ir=1,nr
          dum(ir) = 2*pi*r(ir)*psi(ir)
       enddo
       fact = mdot / ( mdotcap + integrate(nr,r,dum) )
       if(number_invalid(fact).ne.0) then
          write(*,*) 'fact = ',fact
          stop
       endif
       if((fact.gt.1.5).or.(fact.lt.1.d0/1.5)) then
          write(*,*) 'ERROR: Something wrong with normalization'
          write(*,*) '       of the mass loading...'
          write(*,*) '       Factor wrong is: ',fact
          write(*,*) mdot,mdotcap,integrate(nr,r,dum)
          stop
       endif
       do ir=1,nr
          psi(ir) = psi(ir) * fact
       enddo
       mdotcap = mdotcap * fact
    else
       !
       ! Apparently there is no mass loading (anymore)
       !
       mdotcap = 0.d0
       do ir=1,nr
          psi(ir) = 0.d0
       enddo
    endif
    !
    ! Done...
    !
  end subroutine shuulrich_infall


  !-------------------------------------------------------------
  ! Reset the save files
  !-------------------------------------------------------------
  subroutine diskevol_reset_savefiles()
    implicit none
    integer :: ndust=0,ir
    open(unit=1,file='sigma.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='sigmad.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='sigmaice.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='sigmavap.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='sigmaplts.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='velo.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='vdust.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='mflux.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='vdrift.dat',status='unknown')
    write(1,*) diskevol_grid_nr,ndust
    close(1)
    open(unit=1,file='temperature.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='visc.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='alpha.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='mdot.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='mstar.dat',status='unknown')
    write(1,*)
    close(1)
    open(unit=1,file='infall.dat',status='unknown')
    write(1,*)
    close(1)
    open(unit=1,file='infall2.dat',status='unknown')
    write(1,*)
    close(1)
    open(unit=1,file='time.dat',status='unknown')
    write(1,*)
    close(1)
    open(unit=1,file='dustsize.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='stokesnr.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='stfrag.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='stdf.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='stdrift.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    open(unit=1,file='etamid.dat',status='unknown')
    write(1,*) diskevol_grid_nr
    close(1)
    !
    open(unit=1,file='grid.info',status='unknown')
    write(1,*) diskevol_grid_nr
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,*) diskevol_grid_r(ir)
    enddo
    close(1)
  end subroutine diskevol_reset_savefiles


  !--------------------------------------------------------------
  !                 SAVE THE DATA BY APPENDING TO FILE
  !--------------------------------------------------------------
  subroutine save_data(time,isave)
    implicit none
    integer :: ir,isave
    doubleprecision :: time
    !
    write(*,*) 'Saving at time ',time/year,' year'
    !
    open(unit=1,file='sigma.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_tprsigma(ir)/(2*pi*diskevol_grid_r(ir))+1.d-98
10     format(E13.6)
11     format(2(E13.6,1X))
    enddo
    close(1)
    !
    open(unit=1,file='sigmad.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) sigmad(ir)
    enddo
    close(1)
    !
    open(unit=1,file='sigmaice.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 5.d-29
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) sice(ir)
    enddo
    close(1)
    !
    open(unit=1,file='sigmavap.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 5.d-29
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) svap(ir)
    enddo
    close(1)
    !
    open(unit=1,file='sigmaplts.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) sigmaplts(ir)
    enddo
    close(1)
    !
    open(unit=1,file='dustsize.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) a1(ir)
    enddo
    close(1)
    !
    open(unit=1,file='stokesnr.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) st1(ir)
    enddo
    close(1)
    !
    open(unit=1,file='stfrag.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) stfrag(ir)
    enddo
    close(1)
    !
    open(unit=1,file='stdf.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) stdf(ir)
    enddo
    close(1)
    !
    open(unit=1,file='stdrift.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) stdrift(ir)
    enddo
    close(1)
    !
    open(unit=1,file='etamid.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 1.d-28
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) etamid(ir)
    enddo
    close(1)
    !
    open(unit=1,file='vdust.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 0.d0
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) vdust(ir)
    enddo
    close(1)
    !
    open(unit=1,file='mflux.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 0.d0
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) mflux(ir)
    enddo
    close(1)
    !
    open(unit=1,file='vdrift.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,nshift
       write(1,10) 0.d0
    enddo
    do ir=1,diskevol_grid_nr-nshift
       write(1,10) vdrift(ir)
    enddo
    close(1)
    !
    open(unit=1,file='velo.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_vr(ir)
    enddo
    close(1)
    !
    open(unit=1,file='temperature.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_temp(ir)
    enddo
    close(1)
    !
    open(unit=1,file='visc.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_nu(ir)
    enddo
    close(1)
    !
    open(unit=1,file='alpha.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_alpnew(ir)
    enddo
    close(1)
    !
    open(unit=1,file='mdot.dat',status='unknown',access='append')
    write(1,*)
    do ir=1,diskevol_grid_nr
       write(1,10) diskevol_mdot(ir)
    enddo
    close(1)
    !
    open(unit=1,file='mstar.dat',status='unknown',access='append')
    write(1,10) diskevol_mstar
    close(1)
    !
    open(unit=1,file='infall.dat',status='unknown',access='append')
    write(1,11) ml_mdotcap,ml_rcentr
    close(1)
    !
    open(unit=1,file='infall2.dat',status='unknown',access='append')
    write(1,11) ml_mdot
    close(1)
    !
    open(unit=1,file='time.dat',status='unknown',access='append')
    write(1,10) time
    close(1)
    !
    open(unit=1,file='time.info',status='unknown')
    write(1,*) isave
    close(1)
    !
  end subroutine save_data



  !--------------------------------------------------------------
  !                       DO ONE TIME STEP
  !
  ! ARGUMENTS:
  !
  !  time           Time in seconds after onset of cloud collapse
  !  dt             Time step in seconds
  !
  ! OTHER INPUT:
  !
  !  Please make sure to set all the parameters in this module:
  !  the diskevol_** parameters at the top of this file.
  !--------------------------------------------------------------
  subroutine diskevol_one_timestep(time,dt)
    implicit none
    integer :: it,itsave,ir,idust,irs
    doubleprecision :: dt,dummy,diskmass,time
    !
    ! Check
    !
    if(dt.le.0.d0) then
       write(*,*) 'ERROR: Cannot have dt<=0'
       stop
    endif
    !
    ! Prepare mass loading sources
    !
    if(ml_mode.ne.0) then
       !
       ! Calculate the sound speed from the background temperature
       !
       diskevol_cloud_asound = 5991.*sqrt(diskevol_tempbg)
       !
       ! Compute the rotating Shu/Ulrich infall model, and compute the
       ! mass accretion rate onto the surface of the disk (psi(R)) and
       ! the mass accretion rate inward of the capture radius (the
       ! `effective radius' of the central star).
       !
       call shuulrich_infall(diskevol_grid_nr,diskevol_grid_r,diskevol_cloud_mass,   &
                             diskevol_cloud_asound,diskevol_cloud_omega,time,        &
                             diskevol_ml_sigdot,ml_mdot,ml_mdotcap,ml_rcentr,        &
                             diskevol_cloud_ismooth,diskevol_cloud_smoothparam1,     &
                             diskevol_cloud_smoothparam2,diskevol_cloud_idistr)
!       if (sum(diskevol_ml_sigdot(:)).eq.0.d0) write(*,*) 'in module: zero'
    endif
    !
    ! Photoevaporation
    !
    !if(diskevol_isw_euv_evap.ne.0) then
    !   call photoevap(diskevol_grid_nr,diskevol_mstar,diskevol_grid_r,                     &
    !        diskevol_tprsigma,diskevol_temp,euvp_phi,diskevol_photoev_sigdot,photoev_izev)
    !endif
    !
    diskevol_photoev_sigdot(:) = 0.d0    ! Zero for now...
    !
    ! (Re)compute the gravitational mass
    !
    if(diskevol_isw_incmdisk.ne.0) then
       !
       ! Include the mass of the disk inward of R to the total
       ! mass inducing the gravity.
       !
       diskevol_mgrav(1)  = diskevol_mstar
       do ir=2,diskevol_grid_nr
          diskevol_mgrav(ir) = diskevol_mgrav(ir-1) + 0.5 *                   &
                       ( diskevol_tprsigma(ir) + diskevol_tprsigma(ir-1) ) *  &
                       ( diskevol_grid_r(ir) - diskevol_grid_r(ir-1) )
       enddo
    else
       !
       ! Only use the star as the source of gravity
       !
       do ir=1,diskevol_grid_nr
          diskevol_mgrav(ir) = diskevol_mstar
       enddo
    endif
    !
    ! Compute radial velocity
    !
    call compute_vr_v0_difcoef()
    !
    ! Now do the diffusion step for the gas+dust
    !
    call update_sigma_diff(dt)
    !
    ! If matter infalling onto disk:
    ! Update the mass of the central star by matter falling directly
    ! from the envelope onto the star (i.e. into the capture radius
    ! of the star).
    !
    ! BUGFIX 2017.07.04: If diskevol_isw_growstar was 0, it still loaded mass
    !                    onto the star directly from the envelope. Fixed.
    !
    if((ml_mode.ne.0).and.(diskevol_isw_growstar.ne.0)) then
       diskevol_mstar = diskevol_mstar + dt * ml_mdotcap
    endif
    !
    ! Update the surface density of the disk by the matter falling
    ! onto the disk from the envelope (explicit method).
    !
    if((ml_mode.ne.0).and.(ml_implicit.eq.0)) then
       do ir=1,diskevol_grid_nr
          diskevol_tprsigma(ir) = diskevol_tprsigma(ir) +                                 &
                                  dt*2*pi*diskevol_grid_r(ir)*diskevol_ml_sigdot(ir)
       enddo
    endif
    !
    ! Update mass of the star (if active)
    !
    if(diskevol_isw_growstar.ne.0) then
       diskevol_mstar = diskevol_mstar + dt * abs(diskevol_mdot(1))
    endif
    !
    ! Done
    !
  end subroutine diskevol_one_timestep


end module diskevol_module
