module natconst_module
  !
  ! Natural constants in CGS
  !
  doubleprecision GG,mp,me,hh,kk,ee,cc,st,ss,aa
  parameter(GG  = 6.673d-8)      ! Gravitational constant
  parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
  parameter(me  = 9.1095d-28)    ! Mass of electron        [g]
  parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
  parameter(hh  = 6.6262d-27)    ! Planck's constant       [erg.s]
  parameter(ee  = 4.8032d-10)    ! Unit charge             
  parameter(cc  = 2.997924580000000d10) ! Light speed             [cm/s]
  parameter(st  = 6.6524d-25)    ! Thompson cross-section  [cm^2]
  parameter(ss  = 5.6703d-5)     ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
  parameter(aa  = 7.5657d-15)    ! 4 ss / cc               [erg/cm^3/K^4]
  !
  ! Gas constants                                            
  !                                                              
  doubleprecision muh2, mh2, mH2O, ah2                                    
  parameter(muh2= 2.3000d0)      ! Mean molec weight H2+He 
  parameter(mh2= muh2*mp)        ! Mean weight H2+He gas atom
  parameter(mH2O= 18.*mp)         ! weight of H2O molecule
  parameter(ah2=2.d-15)          ! surface area of gas atom
  !                                                              
  ! Alternative units                                        
  !                                                              
  doubleprecision ev,kev,micr,km,angs
  parameter(ev  = 1.6022d-12)    ! Electronvolt            [erg]
  parameter(kev = 1.6022d-9)     ! Kilo electronvolt       [erg]
  parameter(micr= 1.d-4)         ! Micron                  [cm]
  parameter(km  = 1.d5)          ! Kilometer               [cm]
  parameter(angs= 1.d-8)         ! Angstroem               [cm]
  !                                                              
  ! Astronomy constants                                      
  !                                                              
  doubleprecision LS,RS,MS,AU,pc                           
  parameter(LS  = 3.8525d33)     ! Solar luminosity        [erg/s]
  parameter(RS  = 6.96d10)       ! Solar radius            [cm]
  parameter(MS  = 1.98892d33)    ! Solar mass              [g]
  parameter(TS  = 5.78d3)        ! Solar temperature       [K]
  parameter(AU  = 1.49598d13)    ! Astronomical Unit       [cm]
  parameter(pc  = 3.08572d18)    ! Parsec                  [cm]
  !
  ! Time units                                               
  !
  doubleprecision year,hour,day                            
  parameter(year= 3.15576d7)      ! Year                    [s]
  parameter(hour= 3.6000d3)      ! Hour                    [s]
  parameter(day = 8.64d4)        ! Day                     [s]
  !
  ! Math constants
  !
  doubleprecision pi
  parameter(pi  = 3.1415926535897932385d0)
end module natconst_module
