! SNAP: Servere Nuclear Accident Programme
! Copyright (C) 1992-2017   Norwegian Meteorological Institute

! This file is part of SNAP. SNAP is free software: you can
! redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module checkDomainML
  implicit none
  private

  public :: constrain_to_domain

  contains

!> Purpose:
!>   check if particle is inside domain (sets outside_domain = .true. if outside)
!>   move particle to lowest (highest) level if below lowest (highest) level
subroutine constrain_to_domain(part, outside_domain)
  USE snapdimML, only: nx,ny,nk
  use snapgrdML, only: vlevel
  use particleML, only: Particle

  type(Particle), intent(inout) :: part
  logical, intent(out) :: outside_domain

  real :: vmin, vmax

  outside_domain = .false.

  if (part%x < 1. .OR. part%y < 1. .OR. &
      part%x > nx .OR. part%y > ny) then
      outside_domain = .true.
      return
   endif

  vmin = vlevel(nk)
  vmax = vlevel( 1)
  part%z = min(part%z, vmax)
  part%z = max(part%z, vmin)
end subroutine constrain_to_domain

end module checkDomainML
