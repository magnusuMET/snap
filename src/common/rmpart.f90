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

module rmpartML
  implicit none
  private
    integer, save :: has_reduction = -1
  public rmpart

  contains

!> Remove particles which are inactive
!> or have lost (almost) all mass, in the last case the
!> remaining mass is transferred to to the other particles
!> in the same plume (or to the next plume if none left).
subroutine rmpart(rmlimit)
  USE ISO_FORTRAN_ENV, only: real64
  USE particleML, only: pdata, numeric_limit_rad
  USE snapparML, only: ncomp, run_comp, iparnum, def_comp
  USE releaseML, only: iplume, plume_release, nplume, npart

  !> if particles content drop to rmlimit*initial-content
  !> the particle will be removed and it's content redistributed
  real, intent(in) :: rmlimit

  integer :: m,n,npl,i,i1,i2,iredist, mm

  integer, allocatable, save :: npkeep(:)
  real, allocatable, save :: pbqmin(:), pbqtotal(:), pbqlost(:)
  real, allocatable, save :: pbqdist(:)

  if (.not.allocated(npkeep)) allocate(npkeep(ncomp))
  if (.not.allocated(pbqmin)) allocate(pbqmin(ncomp))
  if (.not.allocated(pbqtotal)) allocate(pbqtotal(ncomp))
  if (.not.allocated(pbqlost)) allocate(pbqlost(ncomp))
  if (.not.allocated(pbqdist)) allocate(pbqdist(ncomp))

  !..initialize
  if (has_reduction .eq. -1) then 
    has_reduction = 0
    do n=1,ncomp
      m = run_comp(n)%to_defined
      if(def_comp(m)%kdrydep == 1 .or. def_comp(m)%kwetdep == 1 &
          .or. def_comp(m)%kdecay == 1) then
      has_reduction = 1
     end if
    end do
  end if

  pbqlost = 0.

  n=0

  do npl=1,nplume

    i1 = iplume(npl)%start
    i2 = iplume(npl)%end

  !..redistribution of lost mass (within one plume)
    if(has_reduction == 1 .AND. i1 > 0) then
      pbqtotal = 0.0
      npkeep = 0
      do i=i1,i2
        m=pdata(i)%icomp
        mm = def_comp(m)%to_running
        if((pdata(i)%rad > numeric_limit_rad) .and. &
           (pdata(i)%rad > (plume_release(npl, mm)*rmlimit))) then
          pbqtotal(mm)=pbqtotal(mm)+pdata(i)%rad
          npkeep(mm) = npkeep(mm) + 1
        else
          pbqlost(mm) = pbqlost(mm) + pdata(i)%rad
          pdata(i)%rad=0.
          pdata(i)%active = .FALSE.
          pdata(i)%x=0.
          pdata(i)%y=0.
        end if
      end do
      iredist=0
      do m=1,ncomp
        pbqdist(m)=0.
        if(pbqlost(m) > 0. .AND. npkeep(m) > 0) then
          pbqdist(m)=pbqlost(m)/float(npkeep(m))
          pbqlost(m) = 0.0
          iredist=1
        end if
      end do
      if(iredist == 1) then
        do i=i1,i2
          if(pdata(i)%rad > 0.0) then
            m=pdata(i)%icomp
            pdata(i)%rad= pdata(i)%rad+pbqdist(m)
          end if
        end do
      end if
    end if

  ! removal of particles outside of the domain
  ! by reordering of plumes
    iplume(npl)%start = n + 1
    do i=i1,i2
    ! reorder all particles, only keep active
      if(pdata(i)%active) then
        n=n+1
        if(n /= i) then
        !             moving paricle to new position in pdata (n < i)
          pdata(n) = pdata(i)
          iparnum(n)=iparnum(i)
        end if
      endif
    end do

  ! updating plume-particle relation, or making plume empty
  ! (plumes are never removed!)
    iplume(npl)%end = n
    if(iplume(npl)%start > iplume(npl)%end) then
      iplume(npl)%start = 0
      iplume(npl)%end = -1
    end if

  end do

! updating the used number of particles
  npart=n

!..note: if pmlost>0 we lost mass inside the grid area
!..      (no later plumes to take the mass).


end subroutine rmpart
end module rmpartML
