module interp
!
! interp.f90
!
! Description:
!  Module containing simple interpolation
!  routines.
!  This file is part of wb08.f90.
!
! Copyright 2018 Patrick Stegmann
! 
implicit none

contains
    
pure real function lin_interp(a,b,t) result(x)
    implicit none
    real, dimension(2), intent(in) :: a,b
    real, intent(in) :: t
    x = b(1) + (t - a(1))*(b(2) - b(1))/(a(2)-a(1))
end function lin_interp

end module interp
