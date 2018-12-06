!
! wb08.f90
!
! Description:
! ============
! 	Implementation of the phytoplankton NPP
!	model of Westberry et al. 2008.
!
! Comments:
! =========
! 	Requested by reviewer 1 for ocean lidar 
!	manuscript.
!
! Copyright 2018 Patrick Stegmann
!
! Westberry_08_model is free software:
! you can redistribute it and/or modify it under the terms of
! the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful,
! but WITHOUT ANY WARRANTY; without even the implied 
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
! PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------

program wb08
  use global
  use par, only: depth 
  implicit none 
  integer hh
  real, dimension(size(depth)) :: Kohlenstoff = &
  (/(0.1, hh=1,size(depth)) /) ! units: mg C m^-3
  ! Actual computations:
  write(*,*) "Calculations started."
  depth_calc: do hh = 1, size(depth)
    call retrieve_carbon(hh,Kohlenstoff(hh))
  end do depth_calc
  output: do hh = 1, size(depth)
    if (Kohlenstoff(hh) .lt. 0.) then
      stop "Negative Carbon density."
    endif
    write(*,*) depth(hh), Kohlenstoff(hh)
  end do output
  !write(*,*) "Program finished."
end program wb08