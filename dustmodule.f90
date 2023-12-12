!===============================================================================================
! DUST EVOLUTION AND PLANETESIMAL FORMATION MODULE
! by J. Drazkowska
!
! For details, see:
!   Drazkowska & Alibert (2017) A&A 608, A92
!   Drazkowska & Dullemond (2018) A&A 614, A62
!
! Solids consist of ice and refractory component. Ice sublimation and re-condensation,
! and water vapour evolution are considered here. Dust sizes are regulated by growth,
! fragmentation, and radial drift. A simple prescription similar to proposed by
! Birnstiel, Klahr & Ercolano (2012) is used for dust sizes. Collective drift effect is
! included.
!===============================================================================================
module dust_module
   use natconst_module
   implicit none
   public
   private :: i
!----------------parameters------------
   doubleprecision, parameter :: a0 = 1.d-4            ! monomer size in cm
   doubleprecision, parameter :: rhoice = 1.d0         ! ice material density in g / cm3
   doubleprecision, parameter :: rhosil = 3.d0         ! rock material density in g / cm3
   integer, parameter         :: nshift = 4            ! number of cells we skip at the inner edge (because of timestep issues)
!-----------variables------------------
   integer                                    :: i
   doubleprecision, dimension(:), allocatable :: sigmad
   doubleprecision, dimension(:), allocatable :: sigmaplts
   doubleprecision, dimension(:), allocatable :: gammad
   doubleprecision, dimension(:), allocatable :: etadrift
   doubleprecision, dimension(:), allocatable :: vdrift
   doubleprecision, dimension(:), allocatable :: lmfp
   doubleprecision, dimension(:), allocatable :: st0
   doubleprecision, dimension(:), allocatable :: st1
   doubleprecision, dimension(:), allocatable :: stfrag
   doubleprecision, dimension(:), allocatable :: stdrift
   doubleprecision, dimension(:), allocatable :: stdf
   doubleprecision, dimension(:), allocatable :: a1
   doubleprecision, dimension(:), allocatable :: dsigmaplts
   doubleprecision, dimension(:), allocatable :: aini, ainiprev
   doubleprecision, dimension(:), allocatable :: stini
   doubleprecision, dimension(:), allocatable :: stbar
   doubleprecision, dimension(:), allocatable :: etamid
   doubleprecision, dimension(:), allocatable :: vdust
   doubleprecision, dimension(:), allocatable :: dr
   doubleprecision, dimension(:), allocatable :: srock
   doubleprecision, dimension(:), allocatable :: sice
   doubleprecision, dimension(:), allocatable :: svap
   doubleprecision, dimension(:), allocatable :: rhop
   doubleprecision, dimension(:), allocatable :: mflux
   contains

!--------------------------------------------------------------
! This subroutine initializes solids at a constant fraction of gas density
! and at a constant monomer size.
!--------------------------------------------------------------
   subroutine dust_init(rgrid,tprsigmag,rout,gtd)
      implicit none
      integer :: NN
      doubleprecision                            :: rout
      doubleprecision                            :: gtd
      doubleprecision, dimension(:), allocatable :: tprsigmag
      doubleprecision, dimension(:), allocatable :: rgrid
      doubleprecision                            :: eta

      eta = 1. / gtd

      NN = size(rgrid(:))-nshift
      allocate(sigmad(NN), sigmaplts(NN), srock(NN), sice(NN), svap(NN))
      allocate(gammad(NN), etadrift(NN), vdrift(NN), lmfp(NN), st0(NN), st1(NN), stfrag(NN), stdrift(NN), &
               stdf(NN), a1(NN), dsigmaplts(NN), aini(NN), stini(NN), stbar(NN), &
               etamid(NN), vdust(NN), dr(NN), rhop(NN), &
               ainiprev(NN), mflux(NN))

      do i = 1, NN
         if (rgrid(i+nshift) < rout) then
            sice(i) = 0.5 * eta * tprsigmag(i+nshift) / (2.*pi*rgrid(i+nshift))
            srock(i) = 0.5 * eta * tprsigmag(i+nshift) / (2.*pi*rgrid(i+nshift))
            svap(i) = 5.e-29
         else
            sice(i) = 5.d-29
            srock(i) = 5.d-29
            svap(i) = 5.e-29
         endif
      enddo
      sigmad(:) = sice(:) + srock(:)
      sigmaplts(:) = 1.d-28
      ainiprev(:) = a0

      do i = 1, NN-1
         dr(i) = rgrid(i+1+nshift) - rgrid(i+nshift)
      enddo
      dr(NN) = dr(NN-1)

   end subroutine dust_init

!--------------------------------------------------------------
! This subroutine de-allocates arrays associated to solids.
!--------------------------------------------------------------
   subroutine dust_cleanup
      implicit none

      deallocate(sigmad, sigmaplts, srock, sice, svap, gammad, etadrift, vdrift, &
               lmfp, st0, st1, stfrag, stdrift, stdf, a1, dsigmaplts, aini, &
               stini, stbar, etamid, vdust, dr, rhop, ainiprev, mflux)

   end subroutine dust_cleanup
   
!--------------------------------------------------------------
! This subroutine calculates total flux and performs advection of given 
! component based on velocity and denisty arrays and the diffusion 
! coefficient provided as input.
!--------------------------------------------------------------

   subroutine advection(dt, rgrid, dr, sigma, sigmag, velocity, Diff, massf)
      doubleprecision, dimension(:), allocatable :: rgrid, sigmag, dr
      doubleprecision, dimension(:) :: Diff, velocity
      doubleprecision :: dt
      doubleprecision, dimension(:), allocatable :: massf
      doubleprecision, dimension(:), allocatable :: sigma
      doubleprecision, dimension(:), allocatable :: rsigma, ddtgdr, diffl, r12, limit, F12
      integer :: NN
      
      NN = size(rgrid(:))-nshift
      
      allocate(rsigma(NN), ddtgdr(NN), diffl(NN), r12(NN), limit(NN), F12(NN))
   
      rsigma(1:NN) = sigma(1:NN) * rgrid(1+nshift:NN+nshift)
      
      ! calculate the diffusion flux
      ddtgdr(2:NN) = (sigma(2:NN)/sigmag(2:NN)-sigma(1:NN-1)/sigmag(1:NN-1)) / &
                     (rgrid(2+nshift:NN+nshift)-rgrid(1+nshift:NN-1+nshift))
      ddtgdr(1) = max(0.d0,ddtgdr(2) - (ddtgdr(3) - ddtgdr(2)) * dr(1)/dr(2))
      diffl(1:NN) = rgrid(1+nshift:NN+nshift) * Diff(1:NN) * sigmag(1:NN) * ddtgdr(1:NN)
   
      ! calculate the flux limiter
      do i = 3, NN-1
         if (velocity(i) < 0) then
            if (abs(rsigma(i) - rsigma(i-1)) > 1.d-28) then
               r12(i) = (rsigma(i+1) - rsigma(i)) / (rsigma(i) - rsigma(i-1))
            else
               r12(i) = 0.0
            endif
         else
            if (abs(rsigma(i) - rsigma(i-1)) > 1.d-28) then
               r12(i) = (rsigma(i-1) - rsigma(i-2)) / (rsigma(i) - rsigma(i-1))
            else
               r12(i) = 0.0
            endif
         endif
      enddo
         if (velocity(NN) < 0) then
            if (abs(rsigma(NN) - rsigma(NN-1)) > 1.d-28) then
               r12(NN) = ( -1.d0* rsigma(NN)) / (rsigma(NN) - rsigma(NN-1))
            else
               r12(NN) = 0.0
            endif
         else
            if (abs(rsigma(NN) - rsigma(NN-1)) > 1.d-28) then
               r12(NN) = (rsigma(NN-1) - rsigma(NN-2)) / (rsigma(NN) - rsigma(NN-1))
            else
               r12(NN) = 0.0
            endif
         endif
         if (velocity(2) < 0) then
            if (abs(rsigma(2) - rsigma(1)) > 1.d-28) then
               r12(2) = (rsigma(3) - rsigma(2)) / (rsigma(2) - rsigma(1))
            else
               r12(2) = 0.0
            endif
         else
            if (abs(rsigma(2) - rsigma(1)) > 1.d-28) then
               r12(2) = rsigma(1) / (rsigma(2) - rsigma(1))
            else
               r12(2) = 0.0
            endif
         endif
         r12(1) = 0.d0

      limit(:) = (r12(:) + abs(r12(:))) / (1.d0 + abs(r12(:)))
      
      ! calculate the total flux
      F12(2:NN) = 0.5d0*velocity(2:NN)*((1.d0+sign(1.d0,velocity(2:NN)))*rsigma(1:NN-1)+(1.d0-sign(1.d0,velocity(2:NN)))  &
                *rsigma(2:NN)) + 0.5d0 * dabs(velocity(2:NN)) * (1.d0-dabs(velocity(2:NN)*dt/dr(2:NN)))*limit(2:NN)  &
                *(rsigma(2:NN)-rsigma(1:NN-1)) - diffl(2:NN)
      F12(1) = min(0.d0, velocity(1)*rsigma(1)) - diffl(1)
      
      massf(1:NN) = 2.*pi*(F12(1:NN)) ! mass flux for the output
      
      ! do the advection - modify the surface density
      rsigma(1:NN-1) = max(rsigma(1:NN-1) - (F12(2:NN) - F12(1:NN-1)) * dt / dr(1:NN-1), 5.d-29*rgrid(1+nshift:NN-1+nshift))
      rsigma(NN) = 5.d-29 * rgrid(NN+nshift)
      sigma(1:NN) = rsigma(1:NN) / rgrid(1+nshift:NN+nshift)
      
      deallocate(rsigma, ddtgdr, diffl, r12, limit, F12)
   
   end subroutine advection

!--------------------------------------------------------------
! This subroutine performs one time-step of dust evolution and planetesimal
! formation.
!--------------------------------------------------------------
   subroutine dust_evol(time,dt,rgrid,tprsigmag,Tmid,Dgas,Mstar,dsigmag,vgas,alphadust,gtd)
      implicit none
!-------input from diskevol---------------------
      doubleprecision                            :: time
      doubleprecision, dimension(:), allocatable :: rgrid
      doubleprecision, dimension(:), allocatable :: tprsigmag
      doubleprecision, dimension(:), allocatable :: vgas
      doubleprecision, dimension(:), allocatable :: Tmid
      doubleprecision, dimension(:), allocatable :: Dgas
      doubleprecision                            :: alphadust
      doubleprecision, dimension(:), allocatable :: dsigmag
      doubleprecision                            :: Mstar
      doubleprecision                            :: gtd
!-----dt regulated by dust----------------------
      doubleprecision                            :: dt
!---------local variables-----------------------
      integer  :: NN
      doubleprecision, dimension(:), allocatable :: sigmag
      doubleprecision, dimension(:), allocatable :: cs
      doubleprecision, dimension(:), allocatable :: vK
      doubleprecision, dimension(:), allocatable :: omegaK
      doubleprecision, dimension(:), allocatable :: Pg
      doubleprecision, dimension(:), allocatable :: Hg
      doubleprecision, dimension(:), allocatable :: Peq, Pvap
      doubleprecision, dimension(:), allocatable :: Rcond, Revap, difference
      doubleprecision, dimension(:), allocatable :: vfrag
      doubleprecision, dimension(:), allocatable :: Ddust
      doubleprecision, dimension(:), allocatable :: mrock, mice, mvap
      doubleprecision                            :: eta, largefrac

      eta = 1./gtd

      NN = size(rgrid(:))-nshift

      allocate(sigmag(NN), cs(NN), vK(NN), omegaK(NN), Pg(NN), Hg(NN), Peq(NN), Pvap(NN), vfrag(NN), &
               Ddust(NN), difference(NN), Rcond(NN), Revap(NN), mrock(NN), mice(NN), mvap(NN))

   ! retrive gas variables from the input
      sigmag(1:NN) = tprsigmag(1+nshift:NN+nshift) / (2.*pi*rgrid(1+nshift:NN+nshift))
      cs(1:NN) = dsqrt(kk * Tmid(1+nshift:NN+nshift) / mH2)
      omegaK(1:NN) = dsqrt(GG*MStar/rgrid(1+nshift:NN+nshift)**3.)
      vK(1:NN) = omegaK(1:NN) * rgrid(1+nshift:NN+nshift)
      Pg(1:NN) = sigmag(1:NN) * cs(1:NN) * omegaK(1:NN) / sqrt(2.*pi)  ! ???
      Hg(1:NN) = cs(1:NN) / omegaK(1:NN)   ! ???

   ! infall
      srock(1:NN) = srock(1:NN) + 0.5 * eta * dsigmag(1+nshift:NN+nshift) * dt
      sice(1:NN) = sice(1:NN) + 0.5 * eta * dsigmag(1+nshift:NN+nshift) * dt
      
   ! evaporation and recondensation
      Peq(1:NN) = 1.14d+13*exp(-1.d0*6062./Tmid(1+nshift:NN+nshift))
      Rcond(1:NN) = 6.*dsqrt(kk*Tmid(1+nshift:NN+nshift)/mH2O)/Hg(1:NN)/(pi*a1(1:NN)*rhop(1:NN)) ! condensation rate
      Revap(1:NN) = 6.*dsqrt(2.*pi*mH2O/(kk*Tmid(1+nshift:NN+nshift)))*Peq(1:NN)/(pi*a1(1:NN)*rhop(1:NN)) ! evaporation rate
      difference(1:NN) = (Rcond(1:NN)*svap(1:NN) - Revap(1:NN))*sice(1:NN)*dt
      where (difference(1:NN) >0.)      ! condensation
         sice(1:NN) = sice(1:NN) + min(difference(1:NN),svap(1:NN)-5.d-29)
         svap(1:NN) = svap(1:NN) - min(difference(1:NN),svap(1:NN)-5.d-29)
      elsewhere    ! evaporation
         svap(1:NN) = svap(1:NN) + min(-1.*difference(1:NN), sice(1:NN)-5.d-29)
         sice(1:NN) = sice(1:NN) - min(-1.*difference(1:NN),sice(1:NN) - 5.d-29)
      end where
      sigmad(1:NN) = srock(1:NN) + sice(1:NN) 

   ! internal density of grains
      rhop(1:NN) = rhoice*rhosil*sigmad(1:NN)/(rhoice*srock(1:NN)+rhosil*sice(1:NN))

   ! fragmentation velocity
      vfrag(1:NN) = 1000./(1. + 9.*dexp(-2.*200.*(sice(1:NN)/sigmad(1:NN)-0.01)))
      where(vfrag < 100.) vfrag = 100.

   ! pressure gradient and maximum drift speed
      do i = 1, NN-1
         gammad(i) = rgrid(i+nshift) * (Pg(i+1)-Pg(i)) / (Pg(i) * (rgrid(i+1+nshift)-rgrid(i+nshift)))
      enddo
      gammad(NN) = gammad(NN-1) + (gammad(NN-1) - gammad(NN-2)) / (rgrid(NN-1+nshift) - &
                rgrid(NN-2+nshift)) * (rgrid(NN+nshift) - rgrid(NN-1+nshift))
      etadrift(:) = cs(:)**2. * gammad(:) / (2. * vK(:)**2.)
      vdrift(:) = etadrift(:) * vK(:)       ! maximum drift speed

   ! mean free path in gas
      lmfp(:) = dsqrt(2.*pi) * Hg(:) * mh2 / (sigmag(:) * ah2)

   ! Stokes number of minimum size / monomer
      do i = 1, NN
         if (a0 > 2.25 * lmfp(i)) then ! Stokes regime
            st0(i) = sqrt(2.*pi) * rhop(i) * ah2 * a0**2. * omegaK(i) / (9. * mh2 * cs(i) )
         else                          ! Epstein regime
            st0(i) = pi * a0 * rhop(i) /  (2. * sigmag(i))
         endif
      enddo

   ! maximum Stokes number
      do i = 1, NN
         if (sigmad(i) > 1.d-6) then
            stfrag(i) = 0.37 * vfrag(i)**2. / (3. * alphadust * cs(i)**2.)
            stdrift(i) = 0.055 * sigmad(i) / (sigmag(i) * dabs(etadrift(i)))
            stdf(i) = 0.37 * vfrag(i) / (2. *dabs(vdrift(i)))
            st1(i) = min(stfrag(i), stdrift(i), stdf(i))
            a1(i) = 2. * sigmag(i) * st1(i) / (pi * rhop(i))
            if (a1(i) > 2.25 * lmfp(i)) then
               a1(i) = dsqrt(9.*cs(i)*mh2*st1(i)/(dsqrt(2.*pi)*rhop(i)*omegaK(i)*ah2))
            endif
            aini(i) = ainiprev(i) * dexp(0.8*dt*omegaK(i)*sigmad(i)/sigmag(i))
            if (aini(i) > 2.25 * lmfp(i)) then ! Stokes regime
               stini(i) = sqrt(2.*pi) * rhop(i) * ah2 * aini(i)**2. * omegaK(i) / (9. * mh2 * cs(i) )
            else                        ! Epstein regime
               stini(i) = pi * aini(i) * rhop(i) /  (2. * sigmag(i))
            endif
            st1(i) = min(st1(i),stini(i))
            a1(i) = min(a1(i),aini(i))
            st1(i) = max(st1(i),st0(i))
            a1(i) = max(a1(i),a0)
         else ! if sigmad < 1.d-6 then all the dust stays at monomer size
            stfrag(i) = st0(i)
            stdrift(i) = st0(i)
            stdf(i) = st0(i)
            stini(i) = st0(i)
            aini(i) = a0
            st1(i) = st0(i)
            a1(i) = a0
         endif
         ainiprev(i) = a1(i)
      enddo

   ! mass-weighted average Stokes number and midplane dust-to-gas ratio
      do i = 1, NN ! decide for each cell which regime is that and thus what is the size distribution
        if (st1(i) > st0(i)) then
           if ( (stfrag(i) < min(stdrift(i),stini(i),stdf(i))) .or. (stdf(i) < min(stfrag(i),stdrift(i),stini(i))) ) then ! fragmentation regime
              stbar(i) = 0.75 * st1(i) + 0.25 * st0(i)
              etamid(i) = 0.75 * sigmad(i) / sigmag(i) * dsqrt((alphadust + st1(i))/alphadust) +  &
                          0.25 * sigmad(i) / sigmag(i) * dsqrt((alphadust + st0(i))/alphadust)
           else ! drift / initial growth regime
              stbar(i) = 0.97 * st1(i) + 0.03 * st0(i)
              etamid(i) = 0.97 * sigmad(i) / sigmag(i) * dsqrt((alphadust + st1(i))/alphadust) +  &
                          0.03 * sigmad(i) / sigmag(i) * dsqrt((alphadust + st0(i))/alphadust)
           endif
        else
           stbar(i) = st0(i)
           etamid(i) = sigmad(i) / sigmag(i) * dsqrt((alphadust + st0(i))/alphadust)
        endif
      enddo

   ! drift velocity
      vdust(1:NN) = (2. * vdrift(1:NN) * stbar(1:NN) + vgas(1+nshift:NN+nshift) * (1. + etamid(1:NN))) / &
                    (stbar(1:NN)**2. + (1 + etamid(1:NN))**2.)

   ! dust diffusion coefficient
      Ddust(1:NN) = alphadust * cs(1:NN) * Hg(1:NN) / (1 + stbar(1:NN)**2.)
        
   ! do advection for each component separately
      call advection(dt, rgrid, dr, srock, sigmag, vdust, Ddust, mrock)
      call advection(dt, rgrid, dr, sice,  sigmag, vdust, Ddust, mice )
      call advection(dt, rgrid, dr, svap,  sigmag, vgas(1+nshift:NN+nshift), Dgas(1+nshift:NN+nshift), mvap)
      sigmad(:) = srock(:) + sice(:)
      mflux(:) = mrock(:) + mice(:) ! this is the pebble flux, just for the output 
      
   ! planetesimal formation
      do i = 1, NN
         if ((sigmad(i) > 1.d-6) .and. (st1(i) > 1.d-2) .and. (etamid(i) > 1.d0)) then
            dsigmaplts(i) = min(1.d-5 * sigmad(i) / (2. * pi / omegaK(i)) * dt, sigmad(i) - 1.d-28)
            sigmaplts(i) = sigmaplts(i) + dsigmaplts(i)
            sice(i) = sice(i) - min(sice(i) / sigmad(i), 1.d0) * dsigmaplts(i)
            srock(i) = srock(i) - min(srock(i) / sigmad(i), 1.d0) * dsigmaplts(i)
            sigmad(i) = sigmad(i) - dsigmaplts(i)
         endif
      enddo

    ! next timestep limitation
      dt = 1.5 * dt
      dt = min(dt, 0.01*minval(dr(:) / dabs(vdust(:))))
      dt = min(dt, 0.01*minval(dr(:)**2./Dgas(1+nshift:NN+nshift)))

      deallocate(sigmag, cs, vK, omegaK, Pg, Hg, vfrag, Pvap, Peq, Rcond, Revap, difference, mrock, mice, mvap)

      return
   end subroutine dust_evol

end module
