      subroutine fldout_nc(iwrite,iunit,filnam,itime,tf1,tf2,tnow,tstep,
     *                  istep,nsteph,ierror)
! netcdf
c  Purpose:  Accumulation for average fields (iwrite=0,1).
c            Make and write output fields (iwrite=1).
c	     Initialization of output accumulation arrays (iwrite=-1).
c            Fields written to a sequential 'model output' file,
c            not opened here.
c
c---------------------------------------------------------------------
c  Field parameter numbers used here:
c     *   8 - surface pressure (if model level output) (hPa)
c     *  17 - precipitation accummulated from start of run (mm)
c     *  58 - mean sea level pressure, mslp (if switched on) (hPa)
c     * 500 - instant height of boundary layer (m)
c     * 501 - average height of boundary layer (m)
c     * 502 - precipitation accummulated between field output (mm)
c	      (better scaling than param. 17 for long runs)
c     * 510 - instant concentration in boundary layer (Bq/m3)
c     * 511 - average concentration in boundary layer (Bq/m3)
c     * 512 - dry deposition (for one time interval)  (Bq/m2)
c     * 513 - wet deposition (for one time interval)  (Bq/m2)
c     * 514 - dry deposition (accumulated from start) (Bq/m2)
c     * 515 - wet deposition (accumulated from start) (Bq/m2)
c     * 516 - instant part of Bq in the boundary layer  (%)
c     * 517 - average part of Bq in the boundary layer  (%)
c     * 518 - accumulated concentration in the lower layer (Bq*hr/m3)
c     * 519 - instant     concentration in the lower layer (Bq/m3)
c     * 521 - BOMB dry deposition (for one time interval),  % of release
c     * 522 - BOMB wet deposition (for one time interval),  % of release
c     * 523 - BOMB dry deposition (accumulated from start), % of release
c     * 524 - BOMB wet deposition (accumulated from start), % of release
c       540-569 - instant concentration in each layer (Bq/m3)
c       570-599 - average concentration in each layer (Bq/m3)
c     * 901 - geographic latitude  (degrees)
c     * 902 - geographic longitude (degrees)
c     * 903 - grid square area (m2)
c     ^------- current output
c---------------------------------------------------------------------
c  Notes:
c    -  data type (field ident. no. 3) is set to:
c              3 if forecast length (itime(5)) is 0 (valid time, +0)
c              2 if forecast length (itime(5)) is greater than 0
c              4 if geographic latitude, longitude, grid square area
c    -  for parameter 510-517 and 521-524 the component is identified
c	in field ident. no. 7 (usually 'level').
c    -  for parameter 540-569 : 540=total, 540+idcomp for each type
c       for parameter 570-599 : 570=total, 570+idcomp for each type
c	(the idcomp maximum is then 29, i.e. max 29 components)
c    -  parameters 500 - 509 for fields not dependant on component
c	parameters 510 - 539 for fields dependant on component
c	(level = idcomp to identify components,  0 = total)
c    -  possible use of grid square area, garea  (param. 903):
c	    Bq in boundary layer
c		= concentration/(garea*hbl)
c	    Bq above boundary layer
c		= (1.0-part_of_bq_in_bl)*concentration/(garea*hbl)
c---------------------------------------------------------------------
c
#if defined(DRHOOK)
      USE PARKIND1  ,ONLY : JPIM     ,JPRB
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
#endif
      implicit none
#if defined(DRHOOK)
      REAL(KIND=JPRB) :: ZHOOK_HANDLE ! Stack variable i.e. do not use SAVE
#endif
c
      include 'snapdim.inc'
      include 'snapfil.inc'
      include 'snapgrd.inc'
      include 'snapfld.inc'
      include 'snappar.inc'
      include 'snaptab.inc'
      include 'snapargos.inc'
      include 'snapdebug.inc'
      include 'netcdf.inc'
c
c     *   ps 8 - surface pressure (if model level output) (hPa)
c     *  accum_prc 17 - precipitation accummulated from start of run (mm)
c     *  mslp 58 - mean sea level pressure, mslp (if switched on) (hPa)
c     * ihbl 500 - instant height of boundary layer (m)
c     * ahbl 501 - average height of boundary layer (m)
c     * 502 - precipitation accummulated between field output (mm)
c       (better scaling than param. 17 for long runs)
c     * icbl 510 - instant concentration in boundary layer (Bq/m3)
c     * acbl 511 - average concentration in boundary layer (Bq/m3)
c     * idd 512 - dry deposition (for one time interval)  (Bq/m2)
c     * iwd 513 - wet deposition (for one time interval)  (Bq/m2)
c     * accdd 514 - dry deposition (accumulated from start) (Bq/m2)
c     * accwd 515 - wet deposition (accumulated from start) (Bq/m2)
c     * 516 - instant part of Bq in the boundary layer  (%)
c     * 517 - average part of Bq in the boundary layer  (%)
c     * acsurf 518 - accumulated concentration in the surface layer (Bq*hr/m3)
c     * icsurf 519 - instant     concentration in the surface layer (Bq/m3)
c     * 521 - BOMB dry deposition (for one time interval),  % of release
c     * 522 - BOMB wet deposition (for one time interval),  % of release
c     * 523 - BOMB dry deposition (accumulated from start), % of release
c     * 524 - BOMB wet deposition (accumulated from start), % of release
c       icml 540-569 - instant concentration in each layer (Bq/m3)
c       acml 570-599 - average concentration in each layer (Bq/m3)
c     * 901 - geographic latitude  (degrees)
c     * 902 - geographic longitude (degrees)
c     * 903 - grid square area (m2)
c     ^------- current output
      logical do_accum_prc, do_ihbl, do_ahbl,
     +   do_icbl, do_acbl, do_idd, do_iwd, do_accdd,
     +   do_accwd, do_acsurf, do_ic_surf, do_icml, do_acml

      integer, save :: ps_varid, accum_prc_varid, mslp_varid,
     +    ihbl_varid(mcomp), ahbl_varid(mcomp),
     +    icbl_varid(mcomp),
     +    acbl_varid(mcomp), idd_varid(mcomp), iwd_varid(mcomp),
     +    accdd_varid(mcomp), accwd_varid(mcomp), acsurf_varid(mcomp),
     +    ic_surf_varid(mcomp), icml_varid(mcomp), acml_varid(mcomp)

      integer   iwrite,iunit,istep,nsteph,ierror
      integer   itime(5)
      real      tf1,tf2,tnow,tstep
      character*(*) filnam
c
      integer   igeofield,naverage,initacc,initargos,ihr
      real      geoparam(6)
c

      integer dimids2d(2),dimids3d(3),dimids4d(4), ipos(4), isize(4)
      integer, save :: x_dimid, y_dimid, k_dimid, t_dimid
      integer, save :: varid, t_varid
      integer, save  :: iftime(5), ihrs, ihrs_pos


      integer          nptot1,nptot2
      double precision bqtot1,bqtot2
      double precision dblscale
      double precision dblfield(nx,ny)
c
      integer iyear,month,iday,ihour,minute,isecond
      integer i,j,k,m,mm,n,id03,id04,ivlvl,idextr,idry,iwet,loop,iparx
      integer itab,ko,lvla,lvlb,ipar,ierr,numfields,ios,iu,i1,i2
      real    rt1,rt2,scale,undef,average,averinv,cscale,dscale,hbl
      real    avg,hrstep,dh,splon,splat
c
      integer   itypef,ltimef,itimef(5),icodef,lspecf,loptf,ioptf
      integer   itimeargos(5)
      integer*2 ispecf(3)
c
      character*256 filename
c
      character*256 string
      integer lenstr
c
      real field1print(nx*ny),field2print(nx*ny)

      data do_accum_prc, do_ihbl, do_ahbl
     +     /.false.,.false.,.false./
      data do_icbl, do_acbl, do_idd, do_iwd, do_accdd
     +     /.false.,.false.,.false.,.false.,.true./
      data do_accwd, do_acsurf, do_ic_surf, do_icml, do_acml
     +     /.true.,.true.,.true.,.true.,.true./

c
c      equivalence (field1(1,1),field1print(1))
c      equivalence (field2(1,1),field2print(1))
c
c..used in xyconvert (x,y -> longitude,latitude)
      data geoparam/1.,1.,1.,1.,0.,0./
c
      data initacc,initargos,igeofield,naverage/0,0,0,0/
      data numfields/0/
c
      data itimeargos/5*0/
c
#if defined(DRHOOK)
      ! Before the very first statement
      IF (LHOOK) CALL DR_HOOK('FLDOUT_NC',0,ZHOOK_HANDLE)
#endif
c
      ierror=0
c
c..initialization
c
      if(imodlevel.eq.1 .and.(nxmc.ne.nx .or. nymc.ne.ny)) imodlevel=0
c
      if(initacc.eq.0) then
        do m=1,ncomp
          do j=1,ny
            do i=1,nx
              depdry(i,j,m)=0.0d0
              depwet(i,j,m)=0.0d0
              accdry(i,j,m)=0.0d0
              accwet(i,j,m)=0.0d0
              concen(i,j,m)=0.0d0
              concacc(i,j,m)=0.0d0
           end do
          end do
        end do
       do j=1,ny
         do i=1,nx
           accprec(i,j)=0.0d0
         end do
       end do
       initacc=1
      end if
c
      if(iwrite.eq.-1) then
c..initialization of routine (not of file)
c..iwrite=-1: istep is no. of field outputs
       n=ncomp
       if(ncomp.gt.1 .and. itotcomp.eq.1) n=n+1
       numfields= 2+n*9
       if(itprof.eq.2) numfields= numfields + ncomp*4
       if(inprecip .gt.0) numfields=numfields+2
       if(imslp    .gt.0) numfields=numfields+1
       if(imodlevel.gt.0) numfields=numfields+n*nk*2+nk+1
        numfields= numfields*istep + 4
       if(numfields.gt.32767) numfields=32767
       do i=1,5
         itimeargos(i)=itime(i)
       end do
#if defined(DRHOOK)
c     before the return statement
      IF (LHOOK) CALL DR_HOOK('FLDOUT_NC',1,ZHOOK_HANDLE)
#endif
       return
      end if
c
c..accumulation for average fields......................................
c
      if(naverage.eq.0) then
c
        do j=1,ny
          do i=1,nx
             avghbl(i,j)=0.0d0
            avgprec(i,j)=0.0d0
         end do
        end do
c
       do m=1,ncomp
          do j=1,ny
            do i=1,nx
              avgbq1(i,j,m)=0.0d0
              avgbq2(i,j,m)=0.0d0
           end do
          end do
        end do
c
c..note: model level output on if nxmc=nx, nymc=ny and imodlevel=1
       if(imodlevel.eq.1) then
         do m=1,ncomp
            do k=1,nk-1
              do j=1,nymc
                do i=1,nxmc
                  avgbq(i,j,k,m)=0.0d0
               end do
              end do
            end do
          end do
       end if
c
      end if
c
      naverage=naverage+1
c
c..for time interpolation
      rt1=(tf2-tnow)/(tf2-tf1)
      rt2=(tnow-tf1)/(tf2-tf1)
      hrstep=1./float(nsteph)
c
c..height of boundary layer
      do j=1,ny
        do i=1,nx
          avghbl(i,j)=avghbl(i,j)+(rt1*hbl1(i,j)+rt2*hbl2(i,j))
        end do
      end do
c
c..precipitation (no time interpolation, but hourly time intervals)
      scale=tstep/3600.
      do j=1,ny
        do i=1,nx
          avgprec(i,j)=avgprec(i,j)+scale*precip(i,j,iprecip)
        end do
      end do
c
      do n=1,npart
        i=nint(pdata(1,n))
        j=nint(pdata(2,n))
ccc     ivlvl=pdata(3,n)*10000.
ccc     k=ivlevel(ivlvl)
       m=iruncomp(icomp(n))
        if(pdata(3,n).ge.pdata(4,n)) then
c..in boundary layer
          avgbq1(i,j,m)=avgbq1(i,j,m)+pdata(9,n)
        else
c..above boundary layer
          avgbq2(i,j,m)=avgbq2(i,j,m)+pdata(9,n)
        end if
      end do
c
c..accumulated/integrated concentration
c
      do m=1,ncomp
       do j=1,ny
         do i=1,nx
           concen(i,j,m)=0.0d0
         end do
       end do
      end do
c
      do n=1,npart
        ivlvl=pdata(3,n)*10000.
        k=ivlayer(ivlvl)
        if(k.eq.1) then
          i=nint(pdata(1,n))
          j=nint(pdata(2,n))
         m=iruncomp(icomp(n))
          concen(i,j,m)= concen(i,j,m)+dble(pdata(9,n))
       end if
      end do
c
      do m=1,ncomp
       do j=1,ny
         do i=1,nx
           if(concen(i,j,m).gt.0.0d0) then
             dh= rt1*hlayer1(i,j,1)+rt2*hlayer2(i,j,1)
             concen(i,j,m)= concen(i,j,m)/(dh*dgarea(i,j))
             concacc(i,j,m)= concacc(i,j,m) + concen(i,j,m)*hrstep
           end if
         end do
       end do
      end do
c
      if(imodlevel.eq.1) then
c
        do n=1,npart
          i=nint(pdata(1,n))
          j=nint(pdata(2,n))
          ivlvl=pdata(3,n)*10000.
          k=ivlayer(ivlvl)
          m=iruncomp(icomp(n))
c..in each sigma/eta (input model) layer
          avgbq(i,j,k,m)=avgbq(i,j,k,m)+pdata(9,n)
       end do
c
      end if
c
      if(iwrite.eq.0) then
#if defined(DRHOOK)
c     before the return statement
      IF (LHOOK) CALL DR_HOOK('FLDOUT_NC',1,ZHOOK_HANDLE)
#endif
        return
      end if
c
c
c..output...............................................................
c
      write(9,*) '*FLDOUT_NC*', numfields
c
      if(numfields.gt.0) then
c..initialization of file
c..remove an existing file and create a completely new one
       if (iunit.ne.30)
         call check(nf_close(iunit))
       numfields=0
       call rmfile(filnam,0,ierror)
       write(9,*) 'creating fldout_nc: ',filnam
       ihrs_pos = 0
       call check(nf_create(filnam, NF_NETCDF4, iunit), filnam)
       call check(nf_def_dim(iunit, "time", NF_UNLIMITED, t_dimid),
     +             "t-dim")
       call check(nf_def_dim(iunit, "x", nx, x_dimid), "x-dim")
       call check(nf_def_dim(iunit, "y", ny, y_dimid), "y-dim")
       call check(nf_def_dim(iunit, "k", nk, k_dimid), "k-dim")

       call check(nf_def_var(iunit, "time",NF_FLOAT,1,t_dimid,t_varid))
       write(string,'(A12,I4,A1,I0.2,A1,I0.2,A1,I0.2,A12)')
     + "hours since ",itime(1),"-",itime(2),"-",itime(3)," ",
     + itime(4),":00:00 +0000"
       call check(nf_put_att_text(iunit, t_varid, "units",
     +      len_trim(string), string))
c
c..store the files base-time
      do i=1,4
        iftime(i) = itime(i)
      end do
      iftime(5) = 0

       dimids2d(1) = x_dimid
       dimids3d(1) = x_dimid
       dimids4d(1) = x_dimid
       dimids2d(2) = y_dimid
       dimids3d(2) = y_dimid
       dimids4d(2) = y_dimid
       dimids3d(3) = t_dimid
       dimids4d(3) = k_dimid
       dimids4d(4) = t_dimid

       if (imodlevel.eq.1) then
         call check(nf_def_var(iunit,"surface_air_pressure",
     +           NF_FLOAT,3,dimids3d,ps_varid))
         call check(NF_DEF_VAR_DEFLATE(iunit, ps_varid, 1,1,1))
         call check(nf_put_att_text(iunit, ps_varid, "units",
     +      3, "hPa"))
         call check(nf_put_att_text(iunit, ps_varid, "standard_name",
     +      20, "surface_air_pressure"))
       end if
       if (imslp.eq.1) then
         call check(nf_def_var(iunit,"air_pressure_at_sea_level",
     +           NF_FLOAT,3,dimids3d,mslp_varid))
         call check(NF_DEF_VAR_DEFLATE(iunit,mslp_varid, 1,1,1))
         call check(nf_put_att_text(iunit, mslp_varid, "units",
     +      3, "hPa"))
         call check(nf_put_att_text(iunit, mslp_varid, "standard_name",
     +      25, "air_pressure_at_sea_level"))
       end if
       if (do_accdd) then
         do m=1,ncomp
           mm=idefcomp(m)
           call check(nf_def_var(iunit,"acc_dry_dep_"//compnamemc(mm),
     +             NF_FLOAT,3,dimids3d,accdd_varid(m)))
           call check(NF_DEF_VAR_DEFLATE(iunit,accdd_varid(m), 1,1,1))
           call check(nf_put_att_text(iunit,accdd_varid(m),"units",
     +        5, "Bq/m2"))
           call check(nf_put_att_text(iunit,accdd_varid(m),"metno_name",
     +        20+len_trim(compnamemc(mm)),
     +        "accumulated_dry_dep_"//TRIM(compnamemc(mm))))
         end do
       end if
       filename=filnam(1:lenstr(filnam,1))//'_level_names'
       open (90,file=filename,access='sequential',form='formatted')
       write(90,1090) 0,'Total'
       do m=1,ncomp
         mm=idefcomp(m)
         k=lenstr(compnamemc(mm),1)
         write(90,1090) idcomp(mm),compnamemc(mm)(1:k)
       end do
 1090   format(1x,i5,1x,'"',a,'"')
       close(90)
       call check(nf_enddef(iunit))
      end if

c set the runtime
      ihrs_pos = ihrs_pos + 1
      call hrdiff(0,0,iftime,itime,ihrs,ierror,ierror)
      call check(NF_PUT_VAR1_REAL(iunit,t_varid,ihrs_pos,FLOAT(ihrs)))
      ipos(1) = 1
      ipos(2) = 1
      ipos(3) = ihrs_pos
      isize(1) = nx
      isize(2) = ny
      isize(3) = 1

c
c..open output felt (field) file
c      call mwfelt(1,filnam,iunit,1,nx*ny,field1,1.0,
c     +            ldata,idata,ierror)
c      if(ierror.ne.0) goto 920
c      call check(nf_open(filnam, NF_WRITE, iunit))
c
c..common field identification.............
c
      do i=1,20
        idata(i)=0
      end do
      idata( 1)=iprodr
      idata( 2)=igridr
      idata( 3)=2
      if(itime(5).eq.0) idata(3)=3
      idata( 4)=itime(5)
c.... idata( 5)= .......... vertical coordinate
c.... idata( 6)= .......... parameter no.
c.... idata( 7)= .......... level or level no. or component id
      idata( 8)=0
c.... idata( 9)= .......... set by gridpar
c.... idata(10)= .......... set by gridpar
c.... idata(11)= .......... set by gridpar
      idata(12)=itime(1)
      idata(13)=itime(2)*100+itime(3)
      idata(14)=itime(4)*100
c.... idata(15)= .......... set by gridpar
c.... idata(16)= .......... set by gridpar
c.... idata(17)= .......... set by gridpar
c.... idata(18)= .......... set by gridpar
c.... idata(19)= .......... 0 or sigma*10000 value
c.... idata(20)= .......... field scaling (automatic the best possible)
c
c..put grid parameters into field identification
c..(into the first 20 words and possibly also after space for data)
      call gridpar(-1,ldata,idata,igtype,nx,ny,gparam,ierror)
c
c..geographic coordinates etc.
      if(igeofield.eq.0) then
        do j=1,ny
          do i=1,nx
            field1(i,j)=float(i)
            field2(i,j)=float(j)
          end do
        end do
        id03=idata(3)
        id04=idata(4)
        idata( 3)=4
        idata( 4)=0
        idata( 5)=0
        idata( 7)=0
        idata( 8)=0
        idata(19)=0
       call xyconvert(nx*ny,field1,field2,
     +		       igtype,gparam,2,geoparam,ierror)
        if(idebug.eq.1) call ftest('lat',1,1,nx,ny,1,field2,0)
        idata( 6)=901
        idata(20)=-32767
c      call mwfelt(2,filnam,iunit,1,nx*ny,field2,1.0,
c    +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
        if(idebug.eq.1) call ftest('long',1,1,nx,ny,1,field1,0)
        idata( 6)=902
        idata(20)=-32767
c       call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
        if(idebug.eq.1) call ftest('area',1,1,nx,ny,1,garea,0)
        idata( 6)=903
        idata(20)=-32767
c       call mwfelt(2,filnam,iunit,1,nx*ny,garea,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
        idata(3)=id03
        idata(4)=id04
        igeofield=1
      end if
c
      idata( 5)=2
      idata( 7)=1000
      idata( 8)=0
      idata(19)=0
c
      undef=+1.e+35
c
      average=float(naverage)
      averinv=1./float(naverage)
      naverage=0
c
c..fixed base scaling for concentrations (unit 10**-12 g/m3 = 1 picog/m3)
ccc   cscale=10.**12
c
c..fixed base scaling for depositions (unit 10**-9 g/m2 = 1 nanog/m3)
ccc   dscale=10.**9
c
      cscale= 1.0
      dscale= 1.0
c
c..for linear interpolation in time
      rt1=(tf2-tnow)/(tf2-tf1)
      rt2=(tnow-tf1)/(tf2-tf1)
c
c..surface pressure (if model level output, for vertical crossections)
      if(imodlevel.eq.1) then
       do j=1,ny
         do i=1,nx
           field1(i,j)=rt1*ps1(i,j)+rt2*ps2(i,j)
         end do
       end do
        if(idebug.eq.1) call ftest('ps',1,1,nx,ny,1,field1,0)
        idata( 6)=8
        idata(20)=-32767
        call check(NF_PUT_VARA_REAL(iunit, ps_varid,ipos,isize,field1))
c       call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
      end if
c
c..total accumulated precipitation from start of run
      if(inprecip.eq.1) then
       do j=1,ny
         do i=1,nx
           accprec(i,j)=accprec(i,j)+avgprec(i,j)
           field1(i,j)=accprec(i,j)
         end do
       end do
       idextr=nint(float(istep)/float(nsteph))
        if(idebug.eq.1) call ftest('accprec',1,1,nx,ny,1,field1,0)
        idata( 6)=17
        idata(19)=idextr
        idata(20)=-32767
c       call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
        idata(19)=0
      end if
c
c..mslp (if switched on)
      if(imslp.eq.1) then
       do j=1,ny
         do i=1,nx
           field1(i,j)=rt1*pmsl1(i,j)+rt2*pmsl2(i,j)
         end do
       end do
        if(idebug.eq.1) call ftest('mslp',1,1,nx,ny,1,field1,0)
c        idata( 6)=58
c        idata(20)=-32767
c       call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        call check(NF_PUT_VARA_REAL(iunit, mslp_varid, ipos, isize,
     +               field1))
      end if
c
c..instant height of boundary layer
      do j=1,ny
        do i=1,nx
          field4(i,j)=rt1*hbl1(i,j)+rt2*hbl2(i,j)
        end do
      end do
      if(idebug.eq.1) call ftest('hbl',1,1,nx,ny,1,field4,0)
      idata( 6)=500
      idata(20)=-32767
c      call mwfelt(2,filnam,iunit,1,nx*ny,field4,1.0,
c     +            ldata,idata,ierror)
      if(ierror.ne.0) goto 900
c
c..average height of boundary layer
      do j=1,ny
        do i=1,nx
          field1(i,j)=avghbl(i,j)*averinv
        end do
      end do
      if(idebug.eq.1) call ftest('avghbl',1,1,nx,ny,1,field1,0)
      idata( 6)=501
      idata(20)=-32767
c      call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +            ldata,idata,ierror)
      if(ierror.ne.0) goto 900
c
c..precipitation accummulated between field output
      if(inprecip.eq.1) then
       do j=1,ny
         do i=1,nx
           field1(i,j)=avgprec(i,j)
         end do
       end do
       idextr=nint(average*tstep/3600.)
        if(idebug.eq.1) call ftest('prec',1,1,nx,ny,1,field1,0)
        idata( 6)=502
        idata(19)=idextr
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
        idata(19)=0
      end if
c
c..parameters for each component......................................
c
ccc   idata( 5)=3
      idata( 5)=0
      idata( 8)=0
      idata(19)=0
c
      do m=1,ncomp
c
       mm=idefcomp(m)
c
c..using the field level identifier to identify the component
        idata(7)=idcomp(mm)
c
c..instant Bq in and above boundary layer
        do j=1,ny
          do i=1,nx
            field1(i,j)=0.
            field2(i,j)=0.
          end do
        end do
       bqtot1=0.0d0
       bqtot2=0.0d0
       nptot1=0
       nptot2=0
c
        do n=1,npart
         if(icomp(n).eq.mm) then
            i=nint(pdata(1,n))
            j=nint(pdata(2,n))
            if(pdata(3,n).ge.pdata(4,n)) then
              field1(i,j)=field1(i,j)+pdata(9,n)
             bqtot1=bqtot1+dble(pdata(9,n))
             nptot1=nptot1+1
            else
              field2(i,j)=field2(i,j)+pdata(9,n)
             bqtot2=bqtot2+dble(pdata(9,n))
             nptot2=nptot2+1
            end if
          end if
        end do
c
       write(9,*) ' component: ',compname(mm)
       write(9,*) '   Bq,particles in    abl: ',bqtot1,nptot1
       write(9,*) '   Bq,particles above abl: ',bqtot2,nptot2
       write(9,*) '   Bq,particles          : ',bqtot1+bqtot2,
     +						 nptot1+nptot2
c
c..instant part of Bq in boundary layer
       scale=100.
       do j=1,ny
         do i=1,nx
           if(field1(i,j)+field2(i,j).gt.0.) then
             field3(i,j)=scale*field1(i,j)/(field1(i,j)+field2(i,j))
           else
             field3(i,j)=undef
           end if
         end do
       end do
c
c..instant concentration in boundary layer
        do j=1,ny
          do i=1,nx
ccc         hbl=rt1*hbl1(i,j)+rt2*hbl2(i,j)
            hbl=field4(i,j)
            field2(i,j)=cscale*field1(i,j)/(hbl*garea(i,j))
          end do
        end do
        if(idebug.eq.1) call ftest('conc',1,1,nx,ny,1,field2,0)
        idata( 6)=510
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,1,nx*ny,field2,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..average concentration in boundary layer
        do j=1,ny
         do i=1,nx
            field1(i,j)=cscale*avgbq1(i,j,m)
     +			      /(garea(i,j)*avghbl(i,j))
         end do
        end do
        if(idebug.eq.1) call ftest('avgconc',1,1,nx,ny,1,field1,0)
        idata( 6)=511
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..dry deposition
        if(kdrydep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*sngl(depdry(i,j,m))/garea(i,j)
              accdry(i,j,m)=accdry(i,j,m)+depdry(i,j,m)
            end do
          end do
         if(idebug.eq.1) call ftest('dry',1,1,nx,ny,1,field1,0)
          idata( 6)=512
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..wet deposition
        if(kwetdep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*sngl(depwet(i,j,m))/garea(i,j)
              accwet(i,j,m)=accwet(i,j,m)+depwet(i,j,m)
            end do
          end do
          if(idebug.eq.1) call ftest('wet',1,1,nx,ny,1,field1,0)
          idata( 6)=513
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..accumulated dry deposition
        if(kdrydep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*sngl(accdry(i,j,m))/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('adry',1,1,nx,ny,1,field1,0)
          idata( 6)=514
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          call check(NF_PUT_VARA_REAL(iunit, accdd_varid(m),ipos,isize,
     +            field1))

          if(ierror.ne.0) goto 900
        end if
c
c..accumulated wet deposition
        if(kwetdep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*sngl(accwet(i,j,m))/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('awet',1,1,nx,ny,1,field1,0)
          idata( 6)=515
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..instant part of Bq in boundary layer
        if(idebug.eq.1) call ftest('pbq',1,1,nx,ny,1,field3,1)
        idata( 6)=516
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..average part of Bq in boundary layer
       scale=100.
       do j=1,ny
         do i=1,nx
           if(avgbq1(i,j,m)+avgbq2(i,j,m).gt.0.) then
             field3(i,j)=scale*avgbq1(i,j,m)
     +			       /(avgbq1(i,j,m)+avgbq2(i,j,m))
           else
             field3(i,j)=undef
           end if
         end do
       end do
        if(idebug.eq.1) call ftest('apbq',1,1,nx,ny,1,field3,1)
        idata( 6)=517
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..accumulated/integrated concentration
       do j=1,ny
         do i=1,nx
           field3(i,j)= sngl(concacc(i,j,m))
         end do
       end do
        if(idebug.eq.1) call ftest('concac',1,1,nx,ny,1,field3,1)
        idata( 6)=518
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c.....end do m=1,ncomp
      end do
c
c
c..total parameters (sum of all components).............................
c
      if(ncomp.gt.1 .and. itotcomp.eq.1) then
c
c..using the field level identifier to identify component, 0=total
        idata(7)=0
c
c..total instant Bq in and above boundary layer
        do j=1,ny
          do i=1,nx
            field1(i,j)=0.
            field2(i,j)=0.
          end do
        end do
c
        do n=1,npart
          i=nint(pdata(1,n))
          j=nint(pdata(2,n))
          if(pdata(3,n).ge.pdata(4,n)) then
            field1(i,j)=field1(i,j)+pdata(9,n)
          else
            field2(i,j)=field2(i,j)+pdata(9,n)
          end if
        end do
c
c..total instant part of Bq in boundary layer
       scale=100.
       do j=1,ny
         do i=1,nx
           if(field1(i,j)+field2(i,j).gt.0.) then
             field3(i,j)=scale*field1(i,j)/(field1(i,j)+field2(i,j))
           else
             field3(i,j)=undef
           end if
         end do
       end do
c
c..total instant concentration in boundary layer
        do j=1,ny
          do i=1,nx
ccc         hbl=rt1*hbl1(i,j)+rt2*hbl2(i,j)
            hbl=field4(i,j)
            field2(i,j)=cscale*field1(i,j)/(hbl*garea(i,j))
          end do
        end do
        if(idebug.eq.1) call ftest('tconc',1,1,nx,ny,1,field2,0)
        idata( 6)=510
        idata(20)=-32767
        call mwfelt(2,filnam,iunit,1,nx*ny,field2,1.0,
     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..total average concentration in boundary layer
        do j=1,ny
         do i=1,nx
            field1(i,j)=0.
         end do
        end do
       do m=1,ncomp
          do j=1,ny
           do i=1,nx
              field1(i,j)=field1(i,j)+avgbq1(i,j,m)
           end do
         end do
        end do
        do j=1,ny
         do i=1,nx
            field1(i,j)=cscale*field1(i,j)
     +			      /(garea(i,j)*avghbl(i,j))
         end do
        end do
        if(idebug.eq.1) call ftest('tavgconc',1,1,nx,ny,1,field1,0)
        idata( 6)=511
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
        idry=0
        iwet=0
       do m=1,ncomp
         mm=idefcomp(m)
         if(kdrydep(mm).eq.1) idry=1
         if(kwetdep(mm).eq.1) iwet=1
       end do
c
c..total dry deposition
       if(idry.eq.1) then
         do j=1,ny
           do i=1,nx
             field1(i,j)=0.
           end do
         end do
         do m=1,ncomp
           mm=idefcomp(m)
            if(kdrydep(mm).eq.1) then
              do j=1,ny
                do i=1,nx
                  field1(i,j)=field1(i,j)+sngl(depdry(i,j,m))
                end do
              end do
           end if
         end do
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*field1(i,j)/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('tdry',1,1,nx,ny,1,field1,0)
          idata( 6)=512
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..total wet deposition
        if(iwet.eq.1) then
         do j=1,ny
           do i=1,nx
             field1(i,j)=0.
           end do
         end do
         do m=1,ncomp
           mm=idefcomp(m)
            if(kwetdep(mm).eq.1) then
              do j=1,ny
                do i=1,nx
                  field1(i,j)=field1(i,j)+sngl(depwet(i,j,m))
                end do
              end do
           end if
         end do
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*field1(i,j)/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('twet',1,1,nx,ny,1,field1,0)
          idata( 6)=513
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..total accumulated dry deposition
       if(idry.eq.1) then
         do j=1,ny
           do i=1,nx
             field1(i,j)=0.
           end do
         end do
         do m=1,ncomp
           mm=idefcomp(m)
            if(kdrydep(mm).eq.1) then
              do j=1,ny
                do i=1,nx
                  field1(i,j)=field1(i,j)+sngl(accdry(i,j,m))
                end do
              end do
           end if
         end do
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*field1(i,j)/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('tadry',1,1,nx,ny,1,field1,0)
          idata( 6)=514
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..total accumulated wet deposition
        if(iwet.eq.1) then
         do j=1,ny
           do i=1,nx
             field1(i,j)=0.
           end do
         end do
         do m=1,ncomp
           mm=idefcomp(m)
            if(kwetdep(mm).eq.1) then
              do j=1,ny
                do i=1,nx
                  field1(i,j)=field1(i,j)+sngl(accwet(i,j,m))
                end do
              end do
           end if
         end do
          do j=1,ny
            do i=1,nx
              field1(i,j)=dscale*field1(i,j)/garea(i,j)
            end do
          end do
          if(idebug.eq.1) call ftest('tawet',1,1,nx,ny,1,field1,0)
          idata( 6)=515
          idata(20)=-32767
c          call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                ldata,idata,ierror)
          if(ierror.ne.0) goto 900
        end if
c
c..total instant part of Bq in boundary layer
        if(idebug.eq.1) call ftest('tpbq',1,1,nx,ny,1,field3,1)
        idata( 6)=516
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..total average part of Bq in boundary layer
       scale=100.
       do j=1,ny
         do i=1,nx
           field1(i,j)=0.
           field2(i,j)=0.
         end do
       end do
       do m=1,ncomp
         do j=1,ny
           do i=1,nx
             field1(i,j)=field1(i,j)+sngl(avgbq1(i,j,m))
             field2(i,j)=field2(i,j)+sngl(avgbq2(i,j,m))
           end do
         end do
       end do
       do j=1,ny
         do i=1,nx
           if(field1(i,j)+field2(i,j).gt.0.) then
             field3(i,j)=scale*field1(i,j)
     +			       /(field1(i,j)+field2(i,j))
           else
             field3(i,j)=undef
           end if
         end do
       end do
        if(idebug.eq.1) call ftest('tapbq',1,1,nx,ny,1,field3,1)
        idata( 6)=517
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c..total accumulated/integrated concentration
       do j=1,ny
         do i=1,nx
           field3(i,j)=0.
         end do
       end do
       do m=1,ncomp
         do j=1,ny
           do i=1,nx
             field3(i,j)= field3(i,j) + sngl(concacc(i,j,m))
           end do
         end do
       end do
        if(idebug.eq.1) call ftest('concac',1,1,nx,ny,1,field3,1)
        idata( 6)=518
        idata(20)=-32767
c        call mwfelt(2,filnam,iunit,2,nx*ny,field3,1.0,
c     +              ldata,idata,ierror)
        if(ierror.ne.0) goto 900
c
c.....end if(ncomp.gt.1 .and. itotcomp.eq.1) then
      end if
c
c
c..BOMB fields..........................................................
c
      if (itprof.eq.2) then
c
c..bomb parameters for each component.........
c
ccc     idata( 5)=3
        idata( 5)=0
        idata( 8)=0
        idata(19)=0
c
        do m=1,ncomp
c
         mm= idefcomp(m)
c
         if(idebug.eq.1) write(9,*) ' component: ',compname(mm)
c
c..using the field level identifier to identify the component
          idata(7)=idcomp(mm)
c
c..scale to % of total released Bq (in a single bomb)
         dblscale= 100.0d0/dble(totalbq(mm))
c
c..dry deposition
          if(kdrydep(mm).eq.1) then
            do j=1,ny
              do i=1,nx
                field1(i,j)=sngl(dblscale*depdry(i,j,m))
              end do
            end do
            if(idebug.eq.1) call ftest('dry%',1,1,nx,ny,1,field1,0)
            idata( 6)=521
            idata(20)=-32767
c            call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
            if(ierror.ne.0) goto 900
          end if
c
c..wet deposition
          if(kwetdep(mm).eq.1) then
            do j=1,ny
              do i=1,nx
                field1(i,j)=sngl(dblscale*depwet(i,j,m))
              end do
            end do
            if(idebug.eq.1) call ftest('wet%',1,1,nx,ny,1,field1,0)
            idata( 6)=522
            idata(20)=-32767
c            call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
            if(ierror.ne.0) goto 900
          end if
c
c..accumulated dry deposition
          if(kdrydep(mm).eq.1) then
            do j=1,ny
              do i=1,nx
                field1(i,j)=sngl(dblscale*accdry(i,j,m))
              end do
            end do
            if(idebug.eq.1) call ftest('adry%',1,1,nx,ny,1,field1,0)
            idata( 6)=523
            idata(20)=-32767
c            call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
            if(ierror.ne.0) goto 900
          end if
c
c..accumulated wet deposition
          if(kwetdep(mm).eq.1) then
            do j=1,ny
              do i=1,nx
                field1(i,j)=sngl(dblscale*accwet(i,j,m))
              end do
            end do
            if(idebug.eq.1) call ftest('awet%',1,1,nx,ny,1,field1,0)
            idata( 6)=524
            idata(20)=-32767
c            call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
            if(ierror.ne.0) goto 900
          end if
c
c.......end do m=1,ncomp
        end do
c
      end if
c
c
c..model level fields...................................................
c
      if(imodlevel.ne.1) goto 800
c
c..concentration in each layer
c..(height only computed at time of output)
c
      idata( 5)=ivcoor
c
c..loop for 1=average and 2=instant concentration
c..(now computing average first, then using the same arrays for instant)
c
      do loop=1,2
c
        if(loop.eq.1) then
c
         avg=average
         iparx=570
c
        else
c
         avg=1.
         iparx=540
c
         do m=1,ncomp
            do k=1,nk-1
              do j=1,nymc
                do i=1,nxmc
                  avgbq(i,j,k,m)=0.0d0
               end do
              end do
            end do
          end do
c
          do n=1,npart
            i=nint(pdata(1,n))
            j=nint(pdata(2,n))
            ivlvl=pdata(3,n)*10000.
            k=ivlayer(ivlvl)
           m=iruncomp(icomp(n))
c..in each sigma/eta (input model) layer
            avgbq(i,j,k,m)=avgbq(i,j,k,m)+pdata(9,n)
         end do

        end if
c
        do k=1,nk-1
         do j=1,ny
           do i=1,nx
              dh=rt1*hlayer1(i,j,k)+rt2*hlayer2(i,j,k)
             field1(i,j)=dh
             field4(i,j)=dh*garea(i,j)*avg
           end do
         end do
c.. write out layer-height
         if (loop .eq. 1) then
           ko=klevel(k+1)
           lvla=nint(alevel(k+1)*10.)
           lvlb=nint(blevel(k+1)*10000.)
           if(ivcoor.eq.2) lvla=0
c use parameter z (1)
           idata( 6)=1
           idata( 7)=ko
           idata( 8)=lvla
           idata(19)=lvlb
           idata(20)=-32767
c           call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +               ldata,idata,ierror)
           if(ierror.ne.0) goto 900
         end if

         do m=1,ncomp
           do j=1,ny
             do i=1,nx
           avgbq(i,j,k,m)=avgbq(i,j,k,m)/field4(i,j)
             end do
           end do
         end do
        end do
c
c..average concentration in each layer for each type
       do m=1,ncomp
          do k=1,nk-1
            do j=1,ny
              do i=1,nx
                field1(i,j)=cscale*sngl(avgbq(i,j,k,m))
              end do
            end do
            if(idebug.eq.1) call ftest('avconcl',1,1,nx,ny,1,field1,0)
           ko=klevel(k+1)
           lvla=nint(alevel(k+1)*10.)
           lvlb=nint(blevel(k+1)*10000.)
           if(ivcoor.eq.2) lvla=0
           ipar=iparx+idcomp(m)
            idata( 6)=ipar
            idata( 7)=ko
            idata( 8)=lvla
            idata(19)=lvlb
            idata(20)=-32767
c don't write average currently, only instant (loop = 2)
c        if (loop .eq. 2)
c     +       call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
c
            if(ierror.ne.0) goto 900
          end do
       end do
c
c..total average concentration in each layer
       if(ncomp.gt.1 .and. itotcomp.eq.1) then
         do m=2,ncomp
            do k=1,nk-1
              do j=1,ny
                do i=1,nx
                  avgbq(i,j,k,1)=avgbq(i,j,k,1)+avgbq(i,j,k,m)
                end do
              end do
           end do
         end do
          do k=1,nk-1
            do j=1,ny
              do i=1,nx
                field1(i,j)=cscale*avgbq(i,j,k,1)
              end do
            end do
            if(idebug.eq.1) call ftest('tavconcl',1,1,nx,ny,1,field1,0)
           ko=klevel(k+1)
           lvla=nint(alevel(k+1)*10.)
           lvlb=nint(blevel(k+1)*10000.)
           if(ivcoor.eq.2) lvla=0
           ipar=iparx+0
            idata( 6)=ipar
            idata( 7)=ko
            idata( 8)=lvla
            idata(19)=lvlb
            idata(20)=-32767
c            call mwfelt(2,filnam,iunit,1,nx*ny,field1,1.0,
c     +                  ldata,idata,ierror)
            if(ierror.ne.0) goto 900
          end do
       end if
c
c.....end do loop=1,2
      end do
c
  800 ierror=0
c
      do m=1,ncomp
       mm=idefcomp(m)
        if(kdrydep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              depdry(i,j,m)=0.0d0
            end do
          end do
       end if
        if(kwetdep(mm).eq.1) then
          do j=1,ny
            do i=1,nx
              depwet(i,j,m)=0.0d0
            end do
          end do
        end if
      end do
c
c..close output felt (field) file
c       call check(nf_close(iunit))
c      call mwfelt(13,filnam,iunit,1,nx*ny,field1,1.0,
c     +            ldata,idata,ierror)
#if defined(DRHOOK)
c     before the return statement
      IF (LHOOK) CALL DR_HOOK('FLDOUT_NC',1,ZHOOK_HANDLE)
#endif
      return
c
  900 ierror=1
c..close output felt (field) file
      call mwfelt(13,filnam,iunit,1,nx*ny,field1,1.0,
     +            ldata,idata,ierr)
  920 write(9,*) '*FLDOUT_NC*  Terminates due to write error.'
c
#if defined(DRHOOK)
c     before the return statement
      IF (LHOOK) CALL DR_HOOK('FLDOUT_NC',1,ZHOOK_HANDLE)
#endif
      return
      end subroutine fldout_nc