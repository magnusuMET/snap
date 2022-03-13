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

module readfield_ncML
  USE ftestML, only: ftest
  USE om2edotML, only: om2edot
  USE milibML, only: mapfield
  USE snaptabML, only: t2thetafac
  USE netcdf
  USE snapdebug, only: iulog, idebug

  implicit none

  private
  public readfield_nc, check, nfcheckload, calc_2d_start_length, find_index, &
         compute_vertical_coordinates_at_half

  real, parameter, public :: mean_surface_air_pressure = 1013.26

  interface nfcheckload
    module procedure nfcheckload1d, nfcheckload2d, nfcheckload3d
  end interface

  contains

!..search in list of available timesteps with model level data
!> Next index to read from, ntav == 0 if no available meteo
integer function find_index(first, backward, itimei, ihr1, ihr2) result(ntav)
  USE snapfilML, only: kavail, iavail
  USE datetime, only: datetime_t, duration_t

  !> Whether this is the first input (which ignores ihr1)
  logical, intent(in) :: first
  !> Calculate time backwards
  logical, value :: backward
  !> Last time of reading
  type(datetime_t), intent(in) :: itimei
  !> bounds in the age of the next time
  integer, value :: ihr1, ihr2

  integer :: current_index
  type(datetime_t) :: test_date, itime(2)

  if (first) then
    ihr1 = 0
    ! Search opposite direction to find suitable initial fields
    backward = .not.backward
  endif
  ntav = 0

!..search in list of available timesteps with model level data
  if(.not.backward) then
    current_index = kavail(1)
    itime(1) = itimei + duration_t(ihr1)
    itime(2) = itimei + duration_t(ihr2)
  else
  !..using the backward list
    current_index = kavail(2)
    itime(1) = itimei - duration_t(ihr1)
    itime(2) = itimei - duration_t(ihr2)
  end if

  if (.not.backward) then
    write(iulog,*) '*READFIELD* Requested time: ', itimei
    write(iulog,*) '                Time limit: ', itime(1), itime(2)
  else
    write(iulog,*) '*READFIELD* Requested time: ', itimei
    write(iulog,*) '                Time limit: ', itime(2), itime(1)
  endif

  do while (current_index > 0)
    test_date = datetime_t(iavail(current_index)%aYear, &
                           iavail(current_index)%aMonth, &
                           iavail(current_index)%aDay, &
                           iavail(current_index)%aHour) + &
                duration_t(iavail(current_index)%fcHour)

    !..pointer to next timestep (possibly same time)
    if (.not.backward) then
      if (test_date >= itime(1) .and. test_date <= itime(2)) then
        ntav = current_index
        exit
      endif
      current_index = iavail(current_index)%nAvail
    else
      if (test_date <= itime(1) .and. test_date >= itime(2)) then
        ntav = current_index
        exit
      endif
      current_index = iavail(current_index)%pAvail
    end if
  end do
end function

!> Read fields from NetCDF files
subroutine readfield_nc(istep, backward, itimei, ihr1, ihr2, &
    itimefi,ierror)
  USE iso_fortran_env, only: error_unit
  USE snapfilML, only: nctype, iavail, filef
  USE snapfldML, only: &
      xm, ym, u1, u2, v1, v2, w1, w2, t1, t2, ps1, ps2, pmsl1, pmsl2, &
      hbl1, hbl2, hlayer1, hlayer2, garea, dgarea, hlevel1, hlevel2, &
      hlayer1, hlayer2, bl1, bl2, enspos, precip
  USE snapgrdML, only: alevel, blevel, vlevel, ahalf, bhalf, vhalf, &
      gparam, kadd, klevel, ivlevel, imslp, igtype, ivlayer
  USE snapmetML, only: met_params
  USE snapdimML, only: nx, ny, nk
  USE datetime, only: datetime_t, duration_t
!> current timestep (always positive), negative istep means reset
  integer, intent(in) :: istep
!> whether meteorology should be read backwards
  logical, intent(in) :: backward
!> minimal time-offset after itimei
  integer, intent(in) :: ihr1
!> maximal time-offset after itimei
  integer, intent(in) :: ihr2
!> time since last file input
  type(datetime_t), intent(in) :: itimei
!> final time (output)
  type(datetime_t), intent(out) :: itimefi
!> error (output)
  integer, intent(out) :: ierror

! local variables
  integer, save :: ncid = 0
  integer, save :: ntav1, ntav2 = 0
  character(len=1024), save :: file_name = ""
  logical, save :: first_time_read = .true.

  integer :: i, k, ilevel, i1, i2
  integer :: nhdiff
  real :: dxgrid, dygrid
  integer :: kk
  real :: dred, red, p, px

  integer :: timepos, timeposm1
  integer :: start3d(7), start4d(7), count3d(7), count4d(7)

  ierror = 0

  if (istep < 0) then
  ! set 'save' variables to default values,
  ! ncid not needed, will close automatically
    ntav1 = 0
    ntav2 = 0
  end if

!..get time offset in hours (as iavail(n)%oHour)
  ntav1 = ntav2
  ntav2 = find_index(istep < 0, backward, itimei, ihr1, ihr2)

  if(ntav2 == 0) then
    write(iulog,*) '*READFIELD* No model level data available'
    write(error_unit,*) '*READFIELD* No model level data available'
    ierror=1
    return
  end if


  if(idebug == 1) then
    write(iulog,*) 'MODEL LEVEL SEARCH LIST.   ntav2=',ntav2
    write(iulog,*) 'nx,ny,nk: ',nx,ny,nk
    write(iulog,*) 'istep: ',istep
    write(iulog,*) 'itimei(5), ihr1, ihr2:',itimei,ihr1,ihr2
    write(iulog,fmt='(7(1x,i4),1x,i6,2i5)') (iavail(ntav2))
    flush(iulog)
  end if

! time between two inputs
! open the correct file, if required
  if (file_name /= filef(iavail(ntav2)%fileNo)) then
    if (ncid /= 0) then
      call check(nf90_close(ncid), "close ncid")
    end if
    file_name = filef(iavail(ntav2)%fileNo)
    call check(nf90_open(file_name, NF90_NOWRITE, ncid), file_name)
  end if

!     set timepos and nhdiff
  nhdiff = 3
  if (ntav1 /= 0) &
    nhdiff = abs(iavail(ntav2)%oHour - iavail(ntav1)%oHour)

  timepos = iavail(ntav2)%timePos
  timeposm1 = timepos ! Default: No deaccumulation possible
  if (iavail(ntav2)%pavail_same_file /= 0) then
    ! previous timestep in same file for deaccumulation, even if not in list
    timeposm1 = iavail(iavail(ntav2)%pavail_same_file)%timePos
  endif
  itimefi = datetime_t(iavail(ntav2)%aYear, &
                       iavail(ntav2)%aMonth, &
                       iavail(ntav2)%aDay, &
                       iavail(ntav2)%aHour)
  itimefi = itimefi + duration_t(iavail(ntav2)%fcHour)

  if(idebug == 1) then
    write(iulog,*) 'READING DATA FROM file=',trim(file_name)
    write(iulog,*) 'READING DATA FROM position=',timepos, ' for ', &
      itimefi, ', prev. position=',timeposm1,', hours:',nhdiff
  end if




  if( .TRUE. ) then
  !..move data from input time step 2 to 1

    u1(:,:,:) = u2
    v1(:,:,:) = v2
    w1(:,:,:) = w2
    t1(:,:,:) = t2
    hlevel1(:,:,:) = hlevel2
    hlayer1(:,:,:) = hlayer2

    ps1(:,:) = ps2
    bl1(:,:) = bl2
    hbl1(:,:) = hbl2

    if(imslp /= 0) then
      pmsl1(:,:) = pmsl2
    end if

  end if

  do k=nk-kadd,2,-1

  !..input model level no.
    ilevel=klevel(k)

    ! dummy-dim only on 2d
    call calc_2d_start_length(start4d, count4d, nx, ny, ilevel, &
        enspos, timepos, has_2d_dummy_height=.false.)

  !..u
  !     Get the varid of the data variable, based on its name.
    call nfcheckload(ncid, met_params%xwindv, start4d, count4d, u2(:,:,k))

  !..v
    call nfcheckload(ncid, met_params%ywindv, start4d, count4d, v2(:,:,k))
  ! bug in chernobyl borders from destaggering
    where (v2 >= 1e+30)
      v2 = 0.0
    end where

  !..pot.temp. or abs.temp.
    call nfcheckload(ncid, met_params%pottempv, start4d, count4d, t2(:,:,k))

  !..sigma_dot/eta_dot (0 at surface)
  !..eta: eta_dot (or omega) stored in the same levels as u,v,th.
    if (met_params%sigmadotv == '') then
      w2 = 0
    else
      call nfcheckload(ncid, met_params%sigmadotv, &
          start4d, count4d, w2(:,:,k))
    end if

  end do ! k=nk-kadd,2,-1


!..surface pressure, 10m wind and possibly mean sea level pressure,
!..precipitation
  call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
      enspos, timepos, met_params%has_dummy_dim)


! ps
  call nfcheckload(ncid, met_params%psv, start3d, count3d, ps2(:,:))
!  input ps, must be hPa, otherwise:
  if (nctype == 'arome' .OR. nctype == 'dmi_eps' .OR. &
  nctype == 'ec_det' .OR. nctype == 'h12_grib' .OR. &
  nctype == "ec_n1s" .OR. nctype == "SLIM" .OR. &
  nctype == 'gfs_grib_filter_fimex') then
    ps2 = ps2*0.01
  endif


! u10m
! v10m
  if (.not.met_params%use_model_wind_for_10m) then
    call nfcheckload(ncid, met_params%xwind10mv, start3d, count3d, u2(:,:,1))
    call nfcheckload(ncid, met_params%ywind10mv, start3d, count3d, v2(:,:,1))
  else
    if (enspos >= 0) then
      call nfcheckload(ncid, met_params%xwindv, [1, 1, enspos+1, nk, timepos], &
          [nx, ny, 1, 1, 1], u2(:,:,1))
      call nfcheckload(ncid, met_params%ywindv, [1, 1, enspos+1, nk, timepos], &
          [nx, ny, 1, 1, 1], v2(:,:,1))
    else
      call nfcheckload(ncid, met_params%xwindv, [1, 1, nk, timepos], &
          [nx, ny, 1, 1], u2(:,:,1))
      call nfcheckload(ncid, met_params%ywindv, [1, 1, nk, timepos], &
          [nx, ny, 1, 1], v2(:,:,1))
    endif
  endif

!..mean sea level pressure, not used in computations,
!..(only for output to results file)
  if(imslp /= 0) then
    if ( .NOT. met_params%mslpv == '') then
      write(iulog,*) 'Mslp not found. Not important.'
      imslp=0
    else
      call nfcheckload(ncid, met_params%mslpv, start3d, count3d, pmsl2(:,:))
    end if
  end if

  if (met_params%need_precipitation) then
    call read_precipitation(ncid, nhdiff, timepos, timeposm1)
  else
    precip = 0.0
  endif

! first time initialized data
  if (first_time_read) then
    first_time_read = .false.

    call read_vertical_coordinates(ncid, alevel, blevel, vlevel)
    call compute_vertical_coordinates_at_half(alevel, blevel, vlevel, ahalf, bhalf, vhalf)

  !..compute map ratio
    call mapfield(1,0,igtype,gparam,nx,ny,xm,ym,&
        xm, & ! Ignored when icori = 0
        dxgrid,dygrid,ierror)
    if(ierror /= 0) then
      write(iulog,*) 'MAPFIELD ERROR. ierror= ',ierror
      write(error_unit,*) 'MAPFIELD ERROR. ierror= ',ierror
      error stop 255
    end if
    gparam(7)=dxgrid
    gparam(8)=dygrid
  !..size of each grid square (m**2)
    garea(:,:) = abs((dxgrid/xm) * (dygrid/ym))
    dgarea(:,:) = garea

  ! end initialization
  end if

  if (met_params%temp_is_abs) then
  !..abs.temp. -> pot.temp.
    do k=2,nk-kadd
      associate(p => alevel(k) + blevel(k)*ps2(:,:))
        t2(:,:,k) = t2(:,:,k)*t2thetafac(p)
      end associate
    end do
  end if

  if (met_params%sigmadot_is_omega) then
  !..omega -> etadot, or rather etadot derived from continuity-equation (mean of both)
    call om2edot
  else if (met_params%sigmadotv == '') then
  !..omega -> etadot, or rather etadot derived from continuity-equation
    call om2edot
  ! om2edot take means of omega (=0) and continuity-equation, -> use only continuity equation
    w2 = 2.0*w2
  end if

!..sigma_dot/eta_dot 0 at surface
  w2(:,:,1) = 0.0

!..no temperature at or near surface (not used, yet)
  t2(:,:,1) = -999.0
  if(kadd > 0) then
  !..levels added at the top
    dred=0.5/float(kadd)
    red=1.
    kk=nk-kadd
    do k=nk-kadd+1,nk
      red=red-dred
      u2(:,:,k) = u2(:,:,kk)
      v2(:,:,k) = v2(:,:,kk)
      w2(:,:,k) = w2(:,:,kk)*red
      t2(:,:,k) = t2(:,:,kk)
    end do
  end if

  if(backward) then
  ! backward-calculation, switch sign of winds
    u2 = -u2
    v2 = -v2
    w2 = -w2
  end if


! test---------------------------------------------------------------
  write(iulog,*) 'k,k_model,alevel,blevel,vlevel,p,dp:'
  px=alevel(nk)+blevel(nk)*1000.
  do k=nk,1,-1
    p=alevel(k)+blevel(k)*1000.
    write(iulog,fmt='(1x,2i5,f9.2,2f9.5,f8.0,f6.0)') &
        k,klevel(k),alevel(k),blevel(k),vlevel(k),p,p-px
    px=p
  end do

! test---------------------------------------------------------------

  if(idebug == 1) then
    call ftest('u  ', u2, reverse_third_dim=.true.)
    call ftest('v  ', v2, reverse_third_dim=.true.)
    call ftest('w  ', w2, reverse_third_dim=.true.)
    call ftest('t  ', t2, reverse_third_dim=.true.)
    call ftest('ps ', ps2)
    if (istep > 0) &
      call ftest('pre', precip(:,:))
  end if


  if (istep == 0) then
  ! test---------------------------------------------------------------
    write(iulog,*) 'k,ahalf,bhalf,vhalf,p,dp:'
    px=ahalf(nk)+bhalf(nk)*1000.
    do k=nk,1,-1
      p=ahalf(k)+bhalf(k)*1000.
      write(iulog,fmt='(1x,i5,f9.2,2f9.5,f8.0,f6.0)') &
          k,ahalf(k),bhalf(k),vhalf(k),p,p-px
      px=p
    end do
  ! test---------------------------------------------------------------

  !..level table for (vertical) interpolation
  !..(remember that fields are stored bottom to top
  !.. and that all parameters now are in the same levels)
    write(iulog,*) 'ivlevel:'
    write(iulog,*) 'k,i1,i2,vlevel(k+1),vlevel(k)'
    i2=-1
    do k=nk-1,1,-1
      i1=i2+1
      i2=vlevel(k)*10000.
      if(k == 1) i2=10000
      do i=i1,i2
        ivlevel(i)=k
      end do
      write(iulog,*) k,i1,i2,vlevel(k+1),vlevel(k)
    end do

  !..level table for concentration in each sigma/eta layer
  !..(layers here as in the input model, no '10m' layer,
  !.. but ordering bottom to top, reorder at time of output)
    write(iulog,*) 'ivlayer:'
    write(iulog,*) 'k,i1,i2,vhalf(k+1),vhalf(k)'
    i2=-1
    do k=nk-1,1,-1
      i1=i2+1
      i2=nint(vhalf(k)*10000.)
      if(k == 1) i2=10000
      do i=i1,i2
        ivlayer(i)=k
      end do
      write(iulog,*) k,i1,i2,vhalf(k+1),vhalf(k)
    end do
  end if

  return
end subroutine readfield_nc

!> read precipitation
subroutine read_precipitation(ncid, nhdiff, timepos, timeposm1)
  use iso_fortran_env, only: error_unit
  use snapdebug, only: iulog
  use snapmetML, only: met_params
  use snapfldML, only: field1, field2, field3, field4, precip, &
      enspos
  use snapdimML, only: nx, ny
  USE snapfilML, only: nctype

!> open netcdf file
  integer, intent(in) :: ncid
!> time difference in hours between two precip fields
  integer, intent(in) :: nhdiff
!> timestep in file
  integer, intent(in) :: timepos
!> previous timestep
  integer, intent(in) :: timeposm1


  integer :: start3d(7), count3d(7)
  real :: unitScale
  real :: totalprec

  if (met_params%precaccumv /= '') then
  !..precipitation between input time 't1' and 't2'
    if (timepos /= 1) then
      call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
          enspos, timeposm1, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%precaccumv, &
          start3d, count3d, field1(:,:))

      call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
          enspos, timepos, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%precaccumv, &
          start3d, count3d, field2(:,:))

    !..the difference below may get negative due to different scaling
      precip(:, :) = (field2 - field1)/nhdiff
      where (precip < 0.0)
        precip = 0.0
      end where
    end if
  else if (met_params%precstratiaccumv /= '') then
  ! accumulated stratiform and convective precipitation
  !..precipitation between input time 't1' and 't2'
    if (timepos /= 1) then
      call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
          enspos, timeposm1, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%precstratiaccumv, &
          start3d, count3d, field1)
      call nfcheckload(ncid, met_params%precconaccumv, &
          start3d, count3d, field2)
      call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
          enspos, timepos, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%precstratiaccumv, &
          start3d, count3d, field3)
      call nfcheckload(ncid, met_params%precconaccumv, &
          start3d, count3d, field4)
    !..the difference below may get negative due to different scaling
      unitScale = 1.
      if (nctype == 'ec_det' .or. nctype == "ec_n1s") unitScale = 1000.
      precip(:,:) = (field3 + field4) - (field1 + field2)
      precip(:,:) = precip/nhdiff*unitScale
      where (precip < 0.0)
        precip = 0.0
      end where
    else
    ! timepos eq 1, check if precipitation already present / assume dummy step 0
      call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
          enspos, timepos, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%precstratiaccumv, &
          start3d, count3d, field3(:,:))
      call nfcheckload(ncid, met_params%precconaccumv, &
          start3d, count3d, field4(:,:))

      field1 = 0.0
      field2 = 0.0
      totalprec = sum(field3) + sum(field4)

      if (totalprec > 1e-5) then
      !..the difference below may get negative due to different scaling
        write(iulog,*) "found precip in first timestep, assuming ", &
            "empty 0 timestep to deaccumulate precip"
        unitScale = 1.
        if (nctype == 'ec_det' .or. nctype == "ec_n1s") unitScale = 1000.
        precip(:,:) = (field3 + field4) - (field1 + field2)
        precip(:,:) = precip/nhdiff*unitScale
        where (precip < 0.0)
          precip = 0.0
        end where
      endif
    end if
  else if (met_params%total_column_rain /= '') then
      call calc_2d_start_length(start3d, count3d, nx, ny, 1, &
          enspos, timepos, met_params%has_dummy_dim)
      call nfcheckload(ncid, met_params%total_column_rain, &
          start3d, count3d, field3)
      precip(:,:) = field3
      write(error_unit, *) "Check precipation correctness"
  else
  !..non-accumulated emissions in stratiform an convective
    call calc_2d_start_length(start3d, count3d, nx, ny, -1, &
              enspos, timepos, met_params%has_dummy_dim)
    call nfcheckload(ncid, met_params%precstrativrt, &
        start3d, count3d, field1(:,:))
    if (met_params%precconvrt /= '') then
      call nfcheckload(ncid, met_params%precconvrt, &
          start3d, count3d, field2(:,:))
    else
      field2 = 0.
    endif

    unitScale = 3600.*1000.0 ! m/s -> mm/h
    if (nctype == 'gfs_grib_filter_fimex') unitScale = 3600. !kg/m3/s -> mm/h
    precip(:,:) = (field1 + field2)*unitScale
    where (precip < 0.0)
      precip = 0.0
    end where

  end if
end subroutine

subroutine parse_formula_terms_hybrid_sigma_pressure(formula, a, is_ap, b, ps, p0)
  character(len=*), intent(in) :: formula

  character(len=32), intent(out) :: a, b, ps, p0
  logical, intent(out) :: is_ap

  character(len=32) :: name, var

  integer :: istart, iend

  a = ""
  is_ap = .false.
  b = ""
  ps = ""
  p0 = ""

  ! Split after <SPACE><NAME>: <VAR>

  istart = 1
  do while (istart < len(formula))
    do iend=istart,len(formula)
      select case (formula(iend:iend))
        case (" ")
          istart = istart + 1
          cycle
        case (":")
          exit
      end select
    end do
    name = formula(istart:iend-1)
    istart = iend+1

    do iend=istart,len(formula)
      if (formula(iend:iend) /= " ") exit
      istart = istart + 1
    enddo

    do iend=istart,len(formula)
      if (formula(iend:iend) == " ") then
        exit
      endif
    enddo
    var = formula(istart:iend)
    istart = iend+1

    select case (name)
      case("a")
        a = var
      case("ap")
        a = var
        is_ap = .true.
      case("b")
        b = var
      case("ps")
        ps = var
      case("p0")
        p0 = var
      case default
        error stop "Unknown coordinate name"
    end select
  enddo
end subroutine

!> Read vertical coordinates for transformation from k-level
!> to pressure
!>
!> Pressure can be computed with p = a + b*ps, where [p] = hPa, [ps] = hPa
subroutine read_vertical_coordinates(iunit, alevel, blevel, vlevel)
  USE iso_fortran_env, only: error_unit
  USE snapmetML, only: met_params
  USE snapdimML, only: nk
  USE snapgrdML, only: klevel, kadd
  USE netcdf

  integer, intent(in) :: iunit
  real, intent(out) :: alevel(:)
  real, intent(out) :: blevel(:)
  real, intent(out) :: vlevel(:)

  integer :: varid_lev
  integer :: nf_str_len
  character(len=:), allocatable :: standard_name
  character(len=:), allocatable :: formula_terms

  integer :: k, stat

  stat = nf90_inq_varid(iunit, "lev", varid_lev)
  if (stat /= NF90_NOERR) then
    call check(nf90_inq_varid(iunit, "hybrid", varid_lev), "lev")
  endif

  call check(nf90_inquire_attribute(iunit, varid_lev, "standard_name", len=nf_str_len), &
    "standard_name")
  allocate(character(len=nf_str_len) :: standard_name)
  call check(nf90_get_att(iunit, varid_lev, "standard_name", standard_name), "standard_name")

  select case(standard_name)
    case("atmosphere_hybrid_sigma_pressure_coordinate")
      block
      character(len=32) :: name_a, name_b, name_ps, name_p0
      integer :: varid_p0
      logical :: a_is_ap
      real scaling_factor_pressure

      call check(nf90_inquire_attribute(iunit, varid_lev, &
        "formula_terms", len=nf_str_len), &
        "formula_terms")
      allocate(character(len=nf_str_len) :: formula_terms)
      call check(nf90_get_att(iunit, varid_lev, "formula_terms", formula_terms), "formula_terms")
      call parse_formula_terms_hybrid_sigma_pressure(formula_terms, name_a, a_is_ap, name_b, name_ps, name_p0)

      if (name_ps /= met_params%psv) then
        error stop "Mismatch between surface pressure names"
      endif

      do k=nk-kadd,2,-1
        call nfcheckload(iunit, name_a, klevel(k:k), [1], alevel(k:k))
        call nfcheckload(iunit, name_b, klevel(k:k), [1], blevel(k:k))
      enddo

      !..surface
      alevel(1)=0.
      blevel(1)=1.

      ! Normalize alevel by surface pressure
      if (.not. a_is_ap) then
        block
        real :: p0array(1)
        character(len=:), allocatable :: p0_units
        call check(nf90_inq_varid(iunit, name_p0, varid_p0), "Retrieve p0")
        call check(nf90_inquire_attribute(iunit, varid_p0, "units", len=nf_str_len), "get units of p0")
        allocate(character(len=nf_str_len) :: p0_units)
        call check(nf90_get_att(iunit, varid_p0, "units", p0_units), "get units of p0")
        call nfcheckload(iunit, name_p0, [1], [1], p0array)

        scaling_factor_pressure = 1.0
        select case(p0_units)
          case ("Pa")
            scaling_factor_pressure = 1.0/100.0
          case ("hPa")
            scaling_factor_pressure = 1.0
          case default
            error stop "Unknown pressure unit"
        end select

        alevel(:) = alevel(:) * p0array(1) * scaling_factor_pressure
        blevel(:) = blevel(:)

        vlevel(:) = alevel/mean_surface_air_pressure + blevel
        end block
      else
        if (len_trim(name_p0) == 0) then
          error stop "Formula terms contains ap, but no p0 specified"
        endif
        block
        character(len=:), allocatable :: ap_units
        integer :: varid_ap
        call check(nf90_inq_varid(iunit, name_a, varid_ap), "Retrieve ap")
        call check(nf90_inquire_attribute(iunit, varid_ap, "units", len=nf_str_len), "get units of p0")
        allocate(character(len=nf_str_len) :: ap_units)
        call check(nf90_get_att(iunit, varid_ap, "units", ap_units), "get units of p0")
        scaling_factor_pressure = 1.0
        select case(ap_units)
          case ("Pa")
            scaling_factor_pressure = 1.0/100.0
          case ("hPa")
            scaling_factor_pressure = 1.0
          case default
            error stop "Unknown pressure unit"
        end select

        alevel(:) = alevel(:) * scaling_factor_pressure
        blevel(:) = blevel(:)

        vlevel(:) = alevel/mean_surface_air_pressure + blevel
        end block
      endif
      end block
    case("atmosphere_ln_pressure_coordinate","atmosphere_sigma_coordinate")
      error stop "Parametric vertical coordinate not yet supported"
    case default
      write(error_unit, *) "Parametric vertical coordinate: ", standard_name
      error stop "Unknown parametric vertical coordinate"
  end select
end subroutine

subroutine compute_vertical_coordinates_at_half(alevel, blevel, vlevel, ahalf, bhalf, vhalf)
  USE snapmetML, only: met_params
  USE snapdimML, only: nk
  USE snapgrdML, only: klevel
  real, intent(in) :: alevel(:)
  real, intent(in) :: blevel(:)
  real, intent(in) :: vlevel(:)

  real, intent(out) :: ahalf(:)
  real, intent(out) :: bhalf(:)
  real, intent(out) :: vhalf(:)

  integer :: k

  !..half levels where height is found,
  !..alevel and blevel are in the middle of each layer
    ahalf(1)=alevel(1)
    bhalf(1)=blevel(1)
    vhalf(1)=vlevel(1)
  !..check if subselection of levels
    do k=2,nk-1
      if (klevel(k+1) /= klevel(k)-1) then
        met_params%manual_level_selection = .TRUE.
      endif
    end do
    do k=2,nk-1
      if ( .NOT. met_params%manual_level_selection) then
        ahalf(k)=alevel(k)+(alevel(k)-ahalf(k-1))
        bhalf(k)=blevel(k)+(blevel(k)-bhalf(k-1))
        vhalf(k)=ahalf(k)/mean_surface_air_pressure+bhalf(k)
      else
        ahalf(k)=(alevel(k)+alevel(k+1))*0.5
        bhalf(k)=(blevel(k)+blevel(k+1))*0.5
        vhalf(k)=ahalf(k)/mean_surface_air_pressure+bhalf(k)
      end if
    end do
    ahalf(nk)=alevel(nk)
    bhalf(nk)=blevel(nk)
    vhalf(nk)=vlevel(nk)
end subroutine

!> calculate the start and length paramters for slicing
!> a 2d field from a 3-5d dataset
subroutine calc_2d_start_length(start, length, nx, ny, zpos, &
    enspos, tpos, has_2d_dummy_height)
  !> Offset into field
  integer, intent (out) :: start(7)
  !> Length of each field
  integer, intent(out) :: length(7)
  !> size (x) of field
  integer, intent (in) :: nx
  !> size (y) of field
  integer, intent(in) :: ny
  !> hybrid position
  integer, intent(in) :: zpos
  !> ensemble member, should be given in C notation,
  !> with 0 being the first member. Setting to a negative
  !> value disables this.
  integer, intent(in) :: enspos
  !> time on which to compute
  integer, intent(in) :: tpos
  !> if set to true, an empty dimension (size 1) has been used
  !> for height. Ignored if \p zpos >= 1
  logical, intent (in) :: has_2d_dummy_height

  integer :: pos

  start = 1
  length = 1

  length(1) = nx
  length(2) = ny

  pos = 2
! ensemble enspos given in C notation, e.g. first member = 0
  if (enspos >= 0) then
    pos = pos + 1
    start(pos) = enspos+1
  end if
! z
  if (zpos >= 1) then
    pos = pos + 1
    start(pos) = zpos
  elseif (has_2d_dummy_height) then
    pos = pos + 1
  end if
! time
  pos = pos + 1
  start(pos) = tpos

end subroutine calc_2d_start_length


subroutine check(status, errmsg)
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use snapdebug, only: iulog, idebug
  integer, intent ( in) :: status
  character(len=*), intent(in), optional :: errmsg

  if(status /= NF90_NOERR) then
    if (present(errmsg)) then
      write(OUTPUT_UNIT,*) trim(nf90_strerror(status)), ": ", trim(errmsg)
      write(ERROR_UNIT,*) trim(nf90_strerror(status)), ": ", trim(errmsg)
      if (idebug == 1) write(iulog,*) trim(nf90_strerror(status)), ": ", trim(errmsg)
    else
      write(OUTPUT_UNIT, *) trim(nf90_strerror(status))
      write(ERROR_UNIT, *) trim(nf90_strerror(status))
      if (idebug == 1) write(iulog, *) trim(nf90_strerror(status))
    endif
    error stop 1
  endif
end subroutine check

subroutine fillscaleoffset(ncid, varid, fillvalue, scalefactor, offset, status)
  use iso_fortran_env, only: real32

  integer, intent(in) :: ncid, varid
  real(kind=real32), intent(out) :: fillvalue, scalefactor, offset
  integer, intent(out) :: status

  status = nf90_get_att(ncid, varid, "_FillValue", fillvalue)
  if (status == NF90_ENOTATT) then
    fillvalue = NF90_FILL_FLOAT
  else if (status /= NF90_NOERR) then
    return
  endif

  status = nf90_get_att(ncid, varid, "scale_factor", scalefactor)
  if (status == NF90_ENOTATT) then
    scalefactor = 1
  else if (status /= NF90_NOERR) then
    return
  endif

  status = nf90_get_att(ncid, varid, "add_offset", offset)
  if (status == NF90_ENOTATT) then
    offset = 0
    status = NF90_NOERR
  else if (status /= NF90_NOERR) then
    return
  endif
end subroutine fillscaleoffset

subroutine nfcheckload1d(ncid, varname, start, length, field, return_status)
  use ieee_arithmetic, only: ieee_value, IEEE_QUIET_NAN
  use iso_fortran_env, only: real32

  integer, intent(in) :: ncid, start(:), length(:)
  character(len=*), intent(in) :: varname
  real(real32), intent(out) :: field(:)
  !> Return status instead of panic
  integer, intent(out), optional :: return_status

  real(real32) :: factor, offset, fillvalue
  integer :: varid, status

  if (present(return_status)) return_status = NF90_NOERR

  status = nf90_inq_varid(ncid, varname, varid)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  write (iulog,*) "reading "//trim(varname)//", dims: ", "start(1):",start, " size:",length
  status = nf90_get_var(ncid, varid, field, start=start, count=length)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  call fillscaleoffset(ncid, varid, fillvalue, factor, offset, status)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status)

  where (field == fillvalue)
    field = IEEE_VALUE(fillvalue, IEEE_QUIET_NAN)
  end where

  if (factor /= 1. .OR. offset /= 0.) then
    field = field*factor + offset
  end if
end subroutine nfcheckload1d

subroutine nfcheckload2d(ncid, varname, start, length, field, return_status)
  use ieee_arithmetic, only: ieee_value, IEEE_QUIET_NAN
  use iso_fortran_env, only: real32

  integer, intent(in) :: ncid, start(:), length(:)
  character(len=*), intent(in) :: varname
  real(real32), intent(out) :: field(:,:)
  !> Return status instead of panic
  integer, intent(out), optional :: return_status

  real(real32) :: factor, offset, fillvalue
  integer :: varid, status

  if (present(return_status)) return_status = NF90_NOERR

  status = nf90_inq_varid(ncid, varname, varid)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  write (iulog,*) "reading "//trim(varname)//", dims: ", "start(1):",start, " size:",length
  status = nf90_get_var(ncid, varid, field, start=start, count=length)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  call fillscaleoffset(ncid, varid, fillvalue, factor, offset, status)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status)

  where (field == fillvalue)
    field = IEEE_VALUE(fillvalue, IEEE_QUIET_NAN)
  end where

  if (factor /= 1. .OR. offset /= 0.) then
    field = field*factor + offset
  end if
end subroutine nfcheckload2d

subroutine nfcheckload3d(ncid, varname, start, length, field, return_status)
  use ieee_arithmetic, only: ieee_value, IEEE_QUIET_NAN
  use iso_fortran_env, only: real32

  integer, intent(in) :: ncid, start(:), length(:)
  character(len=*), intent(in) :: varname
  real(real32), intent(out) :: field(:,:,:)
  !> Return status instead of panic
  integer, intent(out), optional :: return_status

  real(real32) :: factor, offset, fillvalue
  integer :: varid, status

  if (present(return_status)) return_status = NF90_NOERR

  status = nf90_inq_varid(ncid, varname, varid)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  write (iulog,*) "reading "//trim(varname)//", dims: ", "start(1):",start, " size:",length
  status = nf90_get_var(ncid, varid, field, start=start, count=length)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status, varname)

  call fillscaleoffset(ncid, varid, fillvalue, factor, offset, status)
  if (status /= NF90_NOERR .and. present(return_status)) then
    return_status = status
    return
  endif
  call check(status)

  where (field == fillvalue)
    field = IEEE_VALUE(fillvalue, IEEE_QUIET_NAN)
  end where

  if (factor /= 1. .OR. offset /= 0.) then
    field = field*factor + offset
  end if
end subroutine nfcheckload3d
end module readfield_ncML
