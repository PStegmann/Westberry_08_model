!-------------------------------------------------------
!
! global.f90
!
! Description:
! ============
! 	This file is part of wb08.f90
!
! Comments:
! =========
! 	Requested by reviewer 1 for ocean lidar 
!	manuscript.
!
!
! Copyright © 2018 Patrick Stegmann
!
! This file is part of Westberry_08_model.
!
!-------------------------------------------------------

module GLOBAL

    use par, only: depth, Chl, PRAD
    use interp, only: lin_interp

    implicit none

    real, parameter :: R = 0.1 ! background losses
                               ! units: d^-1
                               
contains

subroutine retrieve_carbon(nn,output)
    implicit none
    integer, intent(in) :: nn ! depth index
    real, intent(out) :: output
    real :: tmp1, tmp2
    real, parameter :: eps = 1.e-6
    integer, parameter :: maxiter = 1.e6
    integer :: mm 
    tmp1 = 2.*eps
    tmp2 = 0.
    mm = 0
    iteration: do while ( abs(tmp1-tmp2) > eps )
      tmp2 = tmp1
      tmp1 = Carbon(depth(nn))
      Chl(nn) = Chlorophyl(depth(nn))
      mm = mm + 1
      ! write(*,*) "iteration: ", mm
      if (mm .gt. maxiter) then
        stop "iteration loop failed to converge."
      endif
    end do iteration
    output = tmp1 
    return
end subroutine retrieve_carbon

real function Carbon(z)
    !
    ! Phytoplankton Carbon
    ! units: "mg C"*m^-3 
    !
    implicit none
    real, intent(in) :: z
    real, parameter :: C_0 = 10.5 ! units: [mg m^-3]
    real, parameter :: z_MLD = 83. ! [m]
    if (z .le. z_MLD) then
      Carbon = C_0
      return
    else
      if ( R .lt. mu(z) ) then
        Carbon = C_0
        return
      else if ( mu(z) .lt. R ) then
        Carbon = C_0*(mu(z)/R)
        return
      else
        STOP "Error in the computation of Carbon(z)."
      end if
    end if
end function Carbon

real function mu(z)
    ! Phytoplankton growth rate
    ! units: 1/d 
    implicit none
    real :: z
    real, parameter :: mu_max = 2. ! d^-1
    real, parameter :: Chl2Cmu0 = 0.0003 ! [mg Chl (mg C)^-1 ]
    mu = mu_max &
       * (Chl2C(z) - Chl2Cmu0)/(Chl2CNTmax(z) - Chl2Cmu0) &
       * (1.0 - exp(-5.*PRAD(z)))
    ! write(*,*) mu, PRAD(z), z, exp(-5.*PRAD(z))
    return 
end function mu

real function Chl2C(z)
    ! Chlorophyl to Carbon ratio as a function
    ! of depth z.
    ! 
    ! units: (mg Chl) (mg C)^-1
    !
    ! Dependencies: PAR 
    implicit none
    real, intent(in) :: z
    real, parameter :: nutrient_stress = 0.
    real, parameter :: depthNO3 = 67 ! [m]
    Chl2C = ( 0.022 + (0.045 - 0.022)*EXP(-3.0*PRAD(z)) ) &
          - (nutrient_stress*(1.0 - EXP(-0.075*(z-depthNO3))))
    return
end function Chl2C

real function Chl2CNTmax(z)
    ! Chlorophyl to Carbon ratio as a function
    ! of depth z.
    ! 
    ! units: (mg Chl) (mg C)^-1
    !
    ! Dependencies: PAR 
    implicit none
    real, intent(in) :: z
    Chl2CNTmax = ( 0.022 + (0.045 - 0.022)*EXP(-3.0*PRAD(z)) )
    return 
end function Chl2CNTmax

real function Chlorophyl(z)
    implicit none
    real, intent(in) :: z
    Chlorophyl = Chl2C(z)*Carbon(z)
    return
end function Chlorophyl

end module GLOBAL 