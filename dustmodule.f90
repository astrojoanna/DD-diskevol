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
   doubleprecision, dimension(:), allocatable :: rsigmad
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
   doubleprecision, dimension(:), allocatable :: F12, Fice12, Fvap12
   doubleprecision, dimension(:), allocatable :: limit, limitice, limitvap
   doubleprecision, dimension(:), allocatable :: r12, rice12, rvap12
   doubleprecision, dimension(:), allocatable :: diffl, difflice, diffvap
   doubleprecision, dimension(:), allocatable :: aini, ainiprev
   doubleprecision, dimension(:), allocatable :: stini
   doubleprecision, dimension(:), allocatable :: stbar
   doubleprecision, dimension(:), allocatable :: etamid
   doubleprecision, dimension(:), allocatable :: vdust
   doubleprecision, dimension(:), allocatable :: ddtgdr, ditgdr, dvtgdr
   doubleprecision, dimension(:), allocatable :: dr
   doubleprecision, dimension(:), allocatable :: time0
   doubleprecision, dimension(:), allocatable :: sice, rsice
   doubleprecision, dimension(:), allocatable :: svap, rsvap
   doubleprecision, dimension(:), allocatable :: rhop
   doubleprecision, dimension(:), allocatable :: srock
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
      allocate(sigmad(NN), sigmaplts(NN), srock(NN), sice(NN), rsice(NN), svap(NN), rsvap(NN))
      allocate(rsigmad(NN), gammad(NN), etadrift(NN), vdrift(NN), lmfp(NN), st0(NN), st1(NN), stfrag(NN), stdrift(NN), &
               stdf(NN), a1(NN), dsigmaplts(NN), F12(NN), limit(NN), r12(NN), diffl(NN), aini(NN), stini(NN), stbar(NN), &
               etamid(NN), vdust(NN), ddtgdr(NN), dr(NN), time0(NN), rhop(NN), rice12(NN), limitice(NN), &
               Fice12(NN), difflice(NN), ditgdr(NN), dvtgdr(NN), diffvap(NN), rvap12(NN), limitvap(NN), Fvap12(NN), &
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
      time0(:) = 0.d0
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

      deallocate(rsigmad, gammad, etadrift, vdrift, lmfp, st0, st1, stfrag, stdrift, dvtgdr, diffvap, rvap12, limitvap, &
               stdf, a1, dsigmaplts, F12, Fice12, limit, limitice, r12, rice12, diffl, aini, stini, stbar, &
               etamid, vdust, ddtgdr, ditgdr, dr, sigmad, sigmaplts, sice, rsice, svap, rsvap, srock, rhop, difflice, &
               Fvap12, ainiprev)

   end subroutine dust_cleanup

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
      doubleprecision, dimension(:), allocatable :: sigmaprev
      doubleprecision, dimension(:), allocatable :: Peq, Pvap
      doubleprecision                            :: difference
      doubleprecision, dimension(:), allocatable :: vfrag
      doubleprecision, dimension(:), allocatable :: Ddust
      doubleprecision                            :: eta, largefrac

      eta = 1./gtd

      NN = size(rgrid(:))-nshift

      allocate(sigmag(NN), cs(NN), vK(NN), omegaK(NN), Pg(NN), Hg(NN), sigmaprev(NN), Peq(NN), Pvap(NN), vfrag(NN), &
               Ddust(NN))

      sigmaprev(:) = sigmad(:)

   ! retrive gas variables from the input
      sigmag(1:NN) = tprsigmag(1+nshift:NN+nshift) / (2.*pi*rgrid(1+nshift:NN+nshift))
      cs(1:NN) = dsqrt(kk * Tmid(1+nshift:NN+nshift) / mH2)
      omegaK(1:NN) = dsqrt(GG*MStar/rgrid(1+nshift:NN+nshift)**3.)
      vK(1:NN) = omegaK(1:NN) * rgrid(1+nshift:NN+nshift)
      Pg(1:NN) = sigmag(1:NN) * cs(1:NN) * omegaK(1:NN) / sqrt(2.*pi)  ! ???
      Hg(1:NN) = cs(1:NN) / omegaK(1:NN)   ! ???

   ! infall
      sigmad(1:NN) = sigmad(1:NN) + eta * dsigmag(1+nshift:NN+nshift) * dt
      sice(1:NN) = sice(1:NN) + 0.5 * eta * dsigmag(1+nshift:NN+nshift) * dt
!      write(*,*) 'in dust mod:', eta, sum(dsigmag(1+nshift:NN+nshift)), dt

      do i = 1, NN
         if ( (sigmaprev(i) < 1.d-6) .and. (sigmad(i) > 1.d-6) ) then
            time0(i) = time   ! this is to delay the start of initial growth until sigmad reaches a reasonable value
         endif
      enddo

   ! internal density of grains
      rhop(1:NN) = rhoice*rhosil*sigmad(1:NN)/(rhoice*(sigmad(1:NN)-sice(1:NN))+rhosil*sice(1:NN))

   ! fragmentation velocity
      do i = 1, NN
         if (sice(i) > 0.01 * sigmad(i)) then
            vfrag(i) = 1.d+3
         else
            vfrag(i) = 1.d+2
         endif
      enddo

   ! start dust and vapour: r*Sigma will be evolved
      rsigmad(1:NN) = sigmad(1:NN) * rgrid(1+nshift:NN+nshift)
      rsice(1:NN) = sice(1:NN) * rgrid(1+nshift:NN+nshift)
      rsvap(1:NN) = svap(1:NN) * rgrid(1+nshift:NN+nshift)

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

   ! dust diffusion coeficient
      Ddust(1:NN) = alphadust * cs(1:NN) * Hg(1:NN) / (1 + stbar(1:NN)**2.)

   ! fluxes calculation - diffusion fluxes
      ddtgdr(2:NN) = (sigmad(2:NN)/sigmag(2:NN)-sigmad(1:NN-1)/sigmag(1:NN-1)) / &
                     (rgrid(2+nshift:NN+nshift)-rgrid(1+nshift:NN-1+nshift))
      ddtgdr(1) = max(0.d0,ddtgdr(2) - (ddtgdr(3) - ddtgdr(2)) * dr(1)/dr(2))
      diffl(1:NN) = rgrid(1+nshift:NN+nshift) * Ddust(1:NN) * sigmag(1:NN) * ddtgdr(1:NN)

      ditgdr(2:NN) = (sice(2:NN)/sigmag(2:NN)-sice(1:NN-1)/sigmag(1:NN-1)) / &
                     (rgrid(2+nshift:NN+nshift)-rgrid(1+nshift:NN-1+nshift))
      ditgdr(1) = max(0.,ditgdr(2) - (ditgdr(3) - ditgdr(2)) * dr(1)/dr(2))
      difflice(1:NN) = rgrid(1+nshift:NN+nshift) * Ddust(1:NN) * sigmag(1:NN) * ditgdr(1:NN)

      dvtgdr(2:NN) = (svap(2:NN)/sigmag(2:NN)-svap(1:NN-1)/sigmag(1:NN-1)) / &
                     (rgrid(2+nshift:NN+nshift)-rgrid(1+nshift:NN-1+nshift))
      dvtgdr(1) = max(0.,dvtgdr(2) - (dvtgdr(3) - dvtgdr(2)) * dr(1) / dr(2))
      diffvap(1:NN) = rgrid(1+nshift:NN+nshift) * Dgas(1+nshift:NN+nshift) * sigmag(1:NN) * dvtgdr(1:NN)

   ! flux limiters
      do i = 3, NN-1
         if (vdust(i) < 0) then
            if (abs(rsigmad(i) - rsigmad(i-1)) > 1.d-28) then
               r12(i) = (rsigmad(i+1) - rsigmad(i)) / (rsigmad(i) - rsigmad(i-1))
               rice12(i) = (rsice(i+1) - rsice(i)) / (rsice(i) - rsice(i-1))
            else
               r12(i) = 0.0
               rice12(i) = 0.0
            endif
         else
            if (abs(rsigmad(i) - rsigmad(i-1)) > 1.d-28) then
               r12(i) = (rsigmad(i-1) - rsigmad(i-2)) / (rsigmad(i) - rsigmad(i-1))
               rice12(i) = (rsice(i-1) - rsice(i-2)) / (rsice(i) - rsice(i-1))
            else
               r12(i) = 0.0
               rice12(i) = 0.0
            endif
         endif
      enddo
         if (vdust(NN) < 0) then
            if (abs(rsigmad(NN) - rsigmad(NN-1)) > 1.d-28) then
               r12(NN) = ( -1.d0* rsigmad(NN)) / (rsigmad(NN) - rsigmad(NN-1))
               rice12(NN) = ( -1.d0* rsice(NN)) / (rsice(NN) - rsice(NN-1))
            else
               r12(NN) = 0.0
               rice12(NN) = 0.0
            endif
         else
            if (abs(rsigmad(NN) - rsigmad(NN-1)) > 1.d-28) then
               r12(NN) = (rsigmad(NN-1) - rsigmad(NN-2)) / (rsigmad(NN) - rsigmad(NN-1))
               rice12(NN) = (rsice(NN-1) - rsice(NN-2)) / (rsice(NN) - rsice(NN-1))
            else
               r12(NN) = 0.0
               rice12(NN) = 0.0
            endif
         endif
         if (vdust(2) < 0) then
            if (abs(rsigmad(2) - rsigmad(1)) > 1.d-28) then
               r12(2) = (rsigmad(3) - rsigmad(2)) / (rsigmad(2) - rsigmad(1))
               rice12(2) = (rsice(3) - rsice(2)) / (rsice(2) - rsice(1))
            else
               r12(2) = 0.0
               rice12(2) = 0.0
            endif
         else
            if (abs(rsigmad(2) - rsigmad(1)) > 1.d-28) then
               r12(2) = rsigmad(1) / (rsigmad(2) - rsigmad(1))
               rice12(2) = rsice(1) / (rsice(2) - rsice(1))
            else
               r12(2) = 0.0
               rice12(2) = 0.0
            endif
         endif
         r12(1) = 0.d0
         rice12(1) = 0.d0

      limit(:) = (r12(:) + abs(r12(:))) / (1.d0 + abs(r12(:)))
      limitice(:) = (rice12(:) + abs(rice12(:))) / (1.d0 + abs(rice12(:)))

      do i = 3, NN-1
         if (rsvap(i) /= rsvap(i-1)) then
            if (vgas(i+nshift) < 0) then
               rvap12(i) = (rsvap(i+1) - rsvap(i)) / (rsvap(i) - rsvap(i-1))
            else
               rvap12(i) = (rsvap(i-1) - rsvap(i-2)) / (rsvap(i) - rsvap(i-1))
            endif
         else
            rvap12(i) = 0.0
         endif
      enddo
      rvap12(1:2) = 0.0
      rvap12(NN) = 0.0

      limitvap(1:NN) = (rvap12(1:NN) + abs(rvap12(1:NN))) / (1. + abs(rvap12(1:NN)))

    ! total fluxes
      F12(2:NN) = 0.5d0*vdust(2:NN)*((1.d0+sign(1.d0,vdust(2:NN)))*rsigmad(1:NN-1)+    &
                 (1.d0-sign(1.d0,vdust(2:NN)))*rsigmad(2:NN)) + 0.5d0 * dabs(vdust(2:NN)) * &
                 (1.d0-dabs(vdust(2:NN)*dt/dr(2:NN)))*limit(2:NN)*(rsigmad(2:NN)-rsigmad(1:NN-1))  &
                - diffl(2:NN)       ! total flux
      F12(1) = min(0.d0, vdust(1)*rsigmad(1)) - diffl(1)

      Fice12(2:NN) = 0.5d0*vdust(2:NN)*((1.d0+sign(1.d0,vdust(2:NN)))*rsice(1:NN-1)+(1.d0-sign(1.d0,vdust(2:NN)))  &
                *rsice(2:NN)) + 0.5d0 * dabs(vdust(2:NN)) * (1.d0-dabs(vdust(2:NN)*dt/dr(2:NN)))*limitice(2:NN)  &
                *(rsice(2:NN)-rsice(1:NN-1)) - difflice(2:NN)
      Fice12(1) = min(0.d0, vdust(1)*rsice(1)) - difflice(1)

      Fvap12(2:NN) = 0.5d0*vgas(2+nshift:NN+nshift)*((1.d0+sign(1.d0,vgas(2+nshift:NN+nshift)))*rsvap(1:NN-1) &
                 +(1.d0-sign(1.d0,vgas(2+nshift:NN+nshift)))*rsvap(2:NN)) + 0.5d0*dabs(vgas(2+nshift:NN+nshift)) &
                 *(1.d0-dabs(vgas(2+nshift:NN+nshift)*dt/dr(2:NN)))*limitvap(2:NN)*(rsvap(2:NN)-rsvap(1:NN-1)) &
                 - diffvap(2:NN)
      Fvap12(1) = min(0.d0, vgas(1+nshift) * rsvap(1)) - diffvap(1)

   ! advection
      rsigmad(1:NN-1) = max(rsigmad(1:NN-1) - (F12(2:NN) - F12(1:NN-1)) * dt / dr(1:NN-1), 1.d-28*rgrid(1+nshift:NN-1+nshift))
      mflux(1:NN) = 2.*pi*F12(1:NN)
      rsigmad(NN) = 1.d-28 *rgrid(NN+nshift)
      sigmad(1:NN) = rsigmad(1:NN) / rgrid(1+nshift:NN+nshift)

      rsice(1:NN-1) = max(rsice(1:NN-1) - (Fice12(2:NN) - Fice12(1:NN-1)) * dt / dr(1:NN-1), 5.e-29*rgrid(1+nshift:NN+nshift))
      rsice(NN) = 5.e-29 * rgrid(NN+nshift)
      sice(1:NN) = rsice(1:NN) / rgrid(1+nshift:NN+nshift)

      rsvap(1:NN-1) = max(rsvap(1:NN-1) - (Fvap12(2:NN) - Fvap12(1:NN-1)) * dt / dr(1:NN-1), rgrid(1+nshift:NN+nshift)*5.e-29)
      rsvap(NN) = rgrid(NN+nshift) * 5.e-29
      svap(1:NN) = rsvap(1:NN) / rgrid(1+nshift:NN+nshift)

    ! evaporation
      Peq(1:NN) = 1.14d+13*exp(-1.d0*6062./Tmid(1+nshift:NN+nshift))
      Pvap(1:NN) = svap(1:NN) * kk * Tmid(1+nshift:NN+nshift) * omegaK(1:NN) / (2.d0 * cs(1:NN) * mH2O)
      do i = 1, NN
         if ((log(Peq(i)/Pvap(i))>0.d0).and.(sice(i)>5.d-29)) then
            difference = min((Peq(i)-Pvap(i))*2.*cs(i)*mh2/(omegaK(i)*kk*Tmid(i+nshift)),sice(i)-5.d-29)
            sigmad(i) = max(sigmad(i) - difference,1.e-28)
            svap(i) = svap(i) + difference
            sice(i) = max(sice(i) - difference,5.e-29)
         endif
      enddo

    ! recondensation
      do i = 1, NN
         if ((log(Peq(i)/Pvap(i))<0.d0).and.(svap(i)>1.d-7).and.(sigmad(i)>1.e-6)) then
            difference = min((Pvap(i)-Peq(i))*2.*cs(i)*mh2/(omegaK(i)*kk*Tmid(i+nshift)),svap(i)-5.e-29)
            sigmad(i) = sigmad(i) + difference
            sice(i) = sice(i) + difference
            svap(i) = max(svap(i) - difference,5.e-29)
         endif
      enddo

   ! planetesimal formation
      do i = 1, NN
         if ((sigmad(i) > 1.e-6) .and. (st1(i) > 1.d-2) .and. (etamid(i) > 1.d0)) then
            dsigmaplts(i) = 1.d-5 * sigmad(i) / (2. * pi / omegaK(i)) * dt
            sigmaplts(i) = sigmaplts(i) + dsigmaplts(i)
            sice(i) = sice(i) - min(sice(i) / sigmad(i),1.) * dsigmaplts(i)
            sigmad(i) = sigmad(i) - dsigmaplts(i)
            rsigmad(i) = sigmad(i) * rgrid(i+nshift)
            rsice(i) = sice(i) * rgrid(i+nshift)
         endif
      enddo

   ! next timestep limitation
      dt = 1.5 * dt
      dt = min(dt, 0.2*minval(dr(:) / dabs(vdust(:))))
      dt = min(dt, 0.2*minval(dr(:)**2./Dgas(1+nshift:NN+nshift)))

      deallocate(sigmag, cs, vK, omegaK, Pg, Hg, vfrag, Pvap, Peq)

      return
   end subroutine dust_evol

end module
