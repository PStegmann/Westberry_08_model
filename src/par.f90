module par
!
! par.f90
!
! Description:
!  module to compute Photosynthetically active radiation
!  (PAR) at a given depth z.
!  The algorithm is based on:
!
!  Morel, A., S. Maritorena (2001): Bio-optical properties
!  of oceanic waters: A reappraisal. J. Geoph. Res. 
!  106, 7163-80.
!
!  This file is part of wb08.f90.
!
! Copyright 2018 Patrick Stegmann
! 

use interp, only: lin_interp
implicit none

integer, parameter :: resolution = 150

integer :: ii 

real, parameter :: pi = 4.0*atan(1.0)

real, parameter, dimension(resolution) :: depth = &
    (/(real(ii)-1.0, ii=1,resolution) /) ! units: m

real, dimension(resolution) :: Chl = &
    (/(0., ii=1,resolution) /) ! units: m

!real, parameter :: E_d0   ! spectral irradiance 
                           ! at the surface.
                           ! units: W*m^-2

real, parameter :: c = 299792458.0             ! speed of light in [m/s] 
real, parameter :: h = 6.626070040*10.**(-34) ! Planck's constant in [J*s]
real, parameter :: k = 1.38064852*10.**(-23)  ! Boltzmann constant in []

real, parameter, dimension(71) :: lambda = &
(/(real(ii)*10.**(-9), ii=350,700,5) /) ! units: nm

real, parameter, dimension(71) :: K_w = &
(/ 0.02710, 0.02380, 0.02160, 0.01880, 0.01770, &
   0.01595, 0.01510, 0.01376, 0.01271, 0.01208, &
   0.01042, 0.00890, 0.00812, 0.00765, 0.00758, &
   0.00768, 0.00770, 0.00792, 0.00885, 0.00990, &
   0.01148, 0.01182, 0.01188, 0.01211, 0.01251, &
   0.01320, 0.01444, 0.01526, 0.01660, 0.01885, &
   0.02188, 0.02701, 0.03385, 0.04090, 0.04214, &
   0.04287, 0.04454, 0.04630, 0.04846, 0.05212, &
   0.05746, 0.06053, 0.06280, 0.06507, 0.07034, &
   0.07801, 0.09038, 0.11076, 0.13584, 0.16792, &
   0.22310, 0.25838, 0.26506, 0.26843, 0.27612, &
   0.28400, 0.29218, 0.30176, 0.31134, 0.32553, &
   0.34052, 0.37150, 0.41048, 0.42947, 0.43946, &
   0.44844, 0.46543, 0.48642, 0.51640, 0.55939, &
   0.62438 /) ! units: m^-1

real, parameter, dimension(71) :: eee = &
(/ 0.77800, 0.76700, 0.75600, 0.73700, 0.72000, &
   0.70000, 0.68500, 0.67300, 0.67000, 0.66000, &
   0.64358, 0.64776, 0.65175, 0.65555, 0.65917, &
   0.66259, 0.66583, 0.66889, 0.67175, 0.67443, &
   0.67692, 0.67923, 0.68134, 0.68327, 0.68501, &
   0.68657, 0.68794, 0.68903, 0.68955, 0.68947, &
   0.68880, 0.68753, 0.68567, 0.68320, 0.68015, &
   0.67649, 0.67224, 0.66739, 0.66195, 0.65591, &
   0.64927, 0.64204, 0.64000, 0.63000, 0.62300, &
   0.61500, 0.61000, 0.61400, 0.61800, 0.62200, &
   0.62600, 0.63000, 0.63400, 0.63800, 0.64200, &
   0.64700, 0.65300, 0.65800, 0.66300, 0.66700, &
   0.67200, 0.67700, 0.68200, 0.68700, 0.69500, &
   0.69700, 0.69300, 0.66500, 0.64000, 0.62000, &
   0.60000 /) ! units: [-]

real, parameter, dimension(71) :: chi = &
(/ 0.15300, 0.14900, 0.14400, 0.14000, 0.13600, &
   0.13100, 0.12700, 0.12300, 0.11900, 0.11800, &
   0.11748, 0.12066, 0.12259, 0.12326, 0.12269, &
   0.12086, 0.11779, 0.11372, 0.10963, 0.10560, &
   0.10165, 0.09776, 0.09393, 0.09018, 0.08649, &
   0.08287, 0.07932, 0.07584, 0.07242, 0.06907, &
   0.06579, 0.06257, 0.05943, 0.05635, 0.05341, &
   0.05072, 0.04829, 0.04611, 0.04419, 0.04253, &
   0.04111, 0.03996, 0.03900, 0.03750, 0.03600, &
   0.03400, 0.03300, 0.03280, 0.03250, 0.03300, &
   0.03400, 0.03500, 0.03600, 0.03750, 0.03850, &
   0.04000, 0.04200, 0.04300, 0.04400, 0.04450, &
   0.04500, 0.04600, 0.04750, 0.04900, 0.05150, &
   0.05200, 0.05050, 0.04400, 0.03900, 0.03400, &
   0.0300 /) ! units: m^2*mg^(-1)

contains

real function PRAD(z)
    !
    ! Photosynthetically active radiation
    ! at depth z.
    !
    ! units: mol*m^(-2)*d^(-1)
    !
    implicit none
    real, intent(in) :: z ! units: [m]
    integer :: jj
    integer :: kk 
    integer :: di  
    real :: tmp 
    real :: ext 
    tmp = 0.
    di = 1
    idx_loop: do while (z .gt. depth(di))
      di = di + 1
    end do idx_loop
   
    spectral_loop: do jj=10,70-1
      ext = 0.
      !depth_loop: do while (z < depth(kk) .and. kk < resolution)
      depth_loop: do kk = 1, di-1 
        ext = ext + 0.5 &
          *( K_pl(depth(kk+1), lambda(jj+1), Chl(kk+1)) &
          + K_pl(depth(kk), lambda(jj), Chl(kk)) ) &
          *( depth(kk+1) - depth(kk) )
      end do depth_loop
      tmp = tmp + planck_src(lambda(jj),5800.0)*exp(-1.*ext) 
    end do spectral_loop
    PRAD = tmp * ( 2.5*10.**(-9))
    return 
end function PRAD

real function K_pl(z, lbda, Chl)
    !
    ! function to compute the extinction 
    ! coefficient of sea water based on 
    ! Morel et al., 2001
    ! 
    ! units: [m^-1]
    !
    implicit none
    real, intent(in) :: z      ! units: [m]
    real, intent(in) :: lbda   ! units: [nm]
    real, intent(in) :: Chl    ! units: [mg/m^3]
    logical, parameter :: debug = .false.
    integer :: jj
    jj = find_index(lbda)
    !if (debug .and. lin_interp(lambda(jj:jj+1),K_w(jj:jj+1),lbda) .lt. 0.) then
    !  stop "K_w negative!"
    !endif
    if (debug .and. K_bio(z,lbda,Chl) .lt. 0.) then
      stop "K_bio negative!"
    end if 
    if (jj .eq. 1) then
      K_pl = K_w(1) + K_bio(z,lambda(1),Chl)
    else if (jj .ge. 70)  then
      K_pl = K_w(70) + K_bio(z,lambda(69),Chl)
    else
      K_pl = lin_interp(lambda(jj:jj+1),K_w(jj:jj+1),lbda) &
           +  K_bio(z, lbda, Chl)
    end if
    return 
end function K_pl

real function K_bio(z, lbda, Chl)
    !
    ! function to compute the extinction 
    ! coefficient of plankton based on 
    ! Morel et al., 2001
    ! 
    ! units: [m^-1]
    !
    implicit none
    real, intent(in) :: z      ! units: [m]
    real, intent(in) :: lbda ! units: [nm]
    real, intent(in) :: Chl    ! units: [mg/m^3]
    real :: chi_interp, ee_interp ! interpolation tempvar
    integer :: jj
    ! find wavelength bracket:
    jj = find_index(lbda)
    ! interpolate chi and exponent:
    chi_interp = lin_interp(lambda(jj:jj+1), &
                            chi(jj:jj+1), &
                            lbda)
    ee_interp = lin_interp(lambda(jj:jj+1), &
                           eee(jj:jj+1), &
                           lbda)
    ! final result:
    K_bio = chi_interp*Chl**ee_interp
    return 
end function K_bio

integer function find_index(val) result(idx)
    ! 
    ! Utility function to find lambda
    ! array index.
    !
    implicit none
    real, intent(in) :: val
    !character(3), intent(in) :: key
    integer :: jj
    ! real, pointer :: zeiger(70)
    ! select case(key)
    ! case ('K_w')
    !   zeiger => K_w(:)
    ! case ('eee')
    !   zeiger => eee(:) 
    ! case ('chi')
    !   zeiger => chi(:)
    ! case default
    !   stop "Wrong array key. in par::find_index."
    ! end select
    jj = 1
    idx_loop: do while (val .ge. lambda(jj) )
      jj = jj + 1
    end do idx_loop
    !nullify(zeiger)
    idx = jj
end function find_index

real(kind=8) function planck_src(lbda, T)
    ! 
    ! Returns the spectral specific emissivity
    ! as a function of wavelength in units of 
    ! [W m^-2 m^-1]
    ! 
    implicit none
    real, intent(in) :: lbda ! wavelength in [nm]
    real, intent(in) :: T    ! temperature in [K]
    planck_src = 10000.*(2.0*pi*h*c**2.)/(lbda**5.) &
           * 1.0/(exp(h*c/(lbda*k*T)) - 1.0 ) &
           / (h*c/lbda) &
           / (6.022e28) &
           * 0.8695 !... Transmission coefficient for water (n=1.333)
    return 
end function planck_src

subroutine print_parameters() 
    implicit none
    do ii=1,71
      write(*,*) lambda(ii), K_w(ii), eee(ii), chi(ii)
    end do 
end subroutine print_parameters

subroutine print_PAR()
    implicit none
    do ii=1,SIZE(depth)
      write(*,*) depth(ii), PRAD(depth(ii))
    end do 
end subroutine print_PAR

subroutine print_planck()
    implicit none
    integer jj
    do jj = 1, size(lambda)
      write(*,*) lambda(jj), planck_src(lambda(jj), 5800.0)
    end do 
end subroutine

subroutine print_Kpl()
    implicit none
    integer jj,kk 
    !jj = 1
    do jj = 1, size(depth)
      do kk = 1,70-1
        write(*,*) lambda(kk:kk+1), depth(jj), K_pl(depth(jj), lambda(kk), 0.1), K_w(kk:kk+1)
      end do  
    end do 
end subroutine print_Kpl

end module par 