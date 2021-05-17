      program basicstate
c     This program prepares the basic state of the zonal averaged 
c     atmosphere (mean zonal wind (U), Brunt-Vaisala frequency (N^2) and
c     Potencial Vorticity (q_y))
c
c     the size of the files and relevant parameters are specified 
c     in wave1.h which is also used by the wave model.
      implicit integer(i-n)
      real wide,Ro,go,Ho,ae,om,top,bottom,cphreal,cphimag
      include 'wave1.h'
      parameter(nzadd=10)
c      parameter(ny=50,nz=50)
      real 
     $     N2(ny,nz),U(ny,nz),qy(ny,nz),y(ny),z(nz),fcor(Ny),beta(Ny),
     $     fz,fy,Uznz,Uyy,Uy,roN2(ny,nz),n2top,n2bot,
     $     T(ny,nz),Tz(nz),
     $     Uobs(nyobs,nzobs),Tobs(nyobs,nzobs),ztemp(nzobs-1+nzadd),
     $     Utemp(nyobs,nzobs-1+nzadd),Ttemp(nyobs,nzobs-1+nzadd),
     $     pobs(nzobs),zobs(nzobs),yobs(nyobs),
     $     Uint1(nzobs-1+nzadd),Uint2(nz),
     $     Tint1(nzobs-1+nzadd),Tint2(nz),
     $     Uread(ndobs,nzobs,nyobs),Tread(ndobs,nzobs,nyobs)

c      character*4 year(24)
      character dirIn*39,dirOut*37,UOfile*52,N2Ofile*53,PVOfile*53
      character*56 UIfile, TIfile
      character thld*5
c
      data pobs
     $     / 1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,150.,200.,250.,300.
     $      ,400.,500.,700.,775.,850.,925.,990./
c     $     / 1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,125.,150.,175.,200.
c     $      ,225.,250.,300.,350.,400.,450.,500/  
c      data pobs	  
c     $     / 1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,125.,150.,175.,200.
c     $      ,225.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.
c     $      ,775.,800.,825.,850.,875.,900.,925.,950.,975.,1000./
c      data pobs
c     $    /1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,700.,
c     $      650.,600.,550.,500.,450.,400.,350.,300.,250.,225.,200.,175.,
c     $      150.,125.,100.,70.,50.,30.,20.,10., 7., 5., 3., 2., 1./
c      data year
c     $     /'1979','1980','1981','1982','1983','1984','1985','1986',
c     $     '1987','1988','1989','1990','1991','1992','1993','1994',
c     $     '1995','1996','1997','1998','1999','2000','2001','2002'/
c     
c     Open the wave3sp.dat file, which contains input information:
      open (11,file='wave3sp.dat',status='unknown')
      read(11,*)ae,om,top,bottom
      read(11,*)s,wide,cphreal,cphimag
      read(11,*)alheight,alwidth,endy,endz
      read(11,*)dum1,dum2,dum3,dum4
      if(dum1.ne.ny.or.dum2.ne.nz) print*,'check wave3.dat, ny,nz'
      read(11,*)alampy,alscaley,amrad,ahrad
      read(11,*)awrad,alamp,alscale,Ho
      read(11,*)idays,iyrs,ny1,ny2
      if(idays.ne.ndobs) print*,'check wave3 nt1- number of days'
      close (11)
      cph=cphreal+(0.,1.)*cphimag
      print*,ndobs,ntobs
c     print*,'enter hemisphere (SH: negative num, NH: positive or zero)'
c      read*,hem
c      thld=2.
c      print*,'enter threshold (for U210) directory : e.g. thld2'
c      read*,thld

      hem=1 !Northern Hemisphere
c      hem=-1 !Southern Hemisphere
c     
c     units are dimensional!
c
      delz=(top-bottom)/float(nz)*Ho
      width=wide*pi/180.
      dely=width/float(ny-1)
c     
c     define y,z 
c
      delyobs=2.*pi/180.
c      delyobs=2.5*pi/180.
      if(hem.ge.0.) then        
         do i=1,ny
            y(i)=dely*float(i-1)
         end do
         do i=1,nyobs
            yobs(i)=(i-1)*delyobs
         end do
      else
         do i=1,ny
            y(i)=dely*float(i-1)-0.5*pi
         end do
         do i=1,nyobs
            yobs(i)=(i-1)*delyobs-pi/2.
         end do
      endif
      do j=1,nz
         z(j)=delz*float(j-1)+bottom*Ho
      end do
      do i=1,ny
         fcor(i)=2*om*sin(y(i))
         beta(i)=2*om*cos(y(i))/ae
      end do

c     create vertical grid for interpolation of observations
      do j=1,nzobs
c     flip direction of height coordinate
         zobs(j)=-log(pobs(nzobs-j+1)/pobs(nzobs))*Ho
c     do not flip direction of height coordinate
c     zobs(j)=-log(pobs(j)/pobs(1))*Ho
      end do
      do j=1,nzobs-1
         ztemp(j)=zobs(j)
      end do
      do j=1,nzadd
         ztemp(nzobs-1+j)=zobs(nzobs)+
     $        float(j)*(top*Ho-zobs(nzobs))/float(nzadd)
      end do
      
      do iyr=1,1
c     open input files
c         print*,'year: ',year(iyr)
c'
         dirIn='./'
         UIfile='ERAUzm19792015NH.dat'
         TIfile='ERATzm19792015NH.dat'
         print*,UIfile,TIfile
         open(unit=41,file=UIfile,form='unformatted',access='direct',
     $        recl=4*nyobs*nzobs*ndobs)             
         open(unit=42,file=TIfile,form='unformatted',access='direct',
     $        recl=4*nyobs*nzobs*ndobs) 
c     
c     call odb_open (idfile,
c     $        '/home/harnik/data/ERA40/U210/UTrunmean/2002.U5drm.nc', 0)
c     call odb_rdvar (idfile, 'U5d', Uread)
c     call odb_close (idfile)
c     
c     open output files
         dirOut='./'
         UOfile ='ERAU19792015NH.dat'
         N2Ofile='ERAN219792015NH.dat'
         PVOfile='ERAPV19792015NH.dat'        
         print*,UOfile,N2Ofile,PVOfile
         open(unit=71,file=UOfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=72,file=N2Ofile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=73,file=PVOfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
c         open(71,file=UOfile)
c         open(72,file=N2Ofile)
c         open(73,file=PVOfile)
      print*, idays

c         do iy=1,idays
c     read basic state
c     flip lat for NH to have equator first
         read(41,rec=1)Uread
         read(42,rec=1)Tread
         print*,Tread(2,10,10)
         do iy=1,idays
            irec=iy
            do j=1,nzobs
               do i=1,nyobs
                  if (hem.lt.0) then
                     Uobs(i,j)=Uread(iy,nzobs-j+1,i)
                     Tobs(i,j)=Tread(iy,nzobs-j+1,i)
                  else
                     Uobs(i,j)=Uread(iy,nzobs-j+1,nyobs-i+1)
                     Tobs(i,j)=Tread(iy,nzobs-j+1,nyobs-i+1)
                  endif
               end do           ! y loop            
            end do              ! z loop
c     
c     interpolate to model grid but first extend obs grid to top of model
c        
            do j=1,nzobs-1
               do i=1,nyobs
                  Utemp(i,j)=Uobs(i,j)
                  Ttemp(i,j)=Tobs(i,j)
               end do
            end do
c     nnU(T)=0: sponge has value of nzobs-1 level, nnU(T)=1: of nzobs level
            nnU=1
            do j=nzobs,nzobs-1+nzadd         
               Utemp(1,j)=Uobs(1,nzobs+nnU-1)
               Utemp(nyobs,j)=Uobs(nyobs,nzobs+nnU-1)
               do i=2,nyobs-1
                  Utemp(i,j)=0.25*Uobs(i-1,nzobs+nnU-1)+
     $                 0.5*Uobs(i,nzobs+nnU-1)+
     $                 0.25*Uobs(i+1,nzobs+nnU-1)
               end do
               Ttemp(1,j)=Tobs(1,nzobs+nnU-1)
               Ttemp(nyobs,j)=Tobs(nyobs,nzobs+nnU-1)
               do i=2,nyobs-1
                  Ttemp(i,j)=0.25*Tobs(i-1,nzobs+nnU-1)+
     $                 0.5*Tobs(i,nzobs+nnU-1)+
     $                 0.25*Tobs(i+1,nzobs+nnU-1)
               end do
            end do
            if (ny.ne.nyobs) then
               call twodint(nyobs,nzobs-1+nzadd,yobs,ztemp,
     $              Utemp,ny,nz,y,z,U,nerr)
               call twodint(nyobs,nzobs-1+nzadd,yobs,ztemp,
     $              Ttemp,ny,nz,y,z,T,nerr)
            else
               do i=1,nyobs        
                  do j=1,nzobs-1+nzadd
                     Uint1(j)=Utemp(i,j)
                     Tint1(j)=Ttemp(i,j)
                  end do
                  call interp(nzobs-1+nzadd,ztemp,Uint1,nz,z,Uint2,nerr)
                  call interp(nzobs-1+nzadd,ztemp,Tint1,nz,z,Tint2,nerr)
                  do j=1,nz
                     U(i,j)=Uint2(j)
                     T(i,j)=Tint2(j)
                  end do
               end do
            end if   
c     smooth T by doing 121 filter on z 10 times, and U 2 times
            do i=1,Ny            
               do jj=1,10
                  do j=2,Nz-1
                     T(i,j)=0.25*T(i,j-1)+0.5*T(i,j)+0.25*T(i,j+1)
                  end do
               end do      
               do jj=1,2
                  do j=2,Nz-1
                     U(i,j)=0.25*U(i,j-1)+0.5*U(i,j)+0.25*U(i,j+1)
                  end do
               end do
c     calculate N2 from temperature
               do j=2,nz-1
                  Tz(j)=(T(i,j+1)-T(i,j-1))/delz/2.
               end do
               Tz(nz)=0.
               Tz(1)=Tz(2)            
               do j=1,nz
                  N2(i,j)=(Tz(j)+T(i,j)*2./7./Ho)*Ro/Ho
               end do
c     smooth N2 by doing 121 filter on z 10 times
               do jj=1,10
                  do j=2,Nz-1
                     N2(i,j)=0.25*N2(i,j-1)+0.5*N2(i,j)+0.25*N2(i,j+1)
                  end do
               end do
               roN2(i,1)=exp(-z(1)/Ho)/N2(i,1)
               roN2(i,nz)=exp(-z(nz)/Ho)/N2(i,nz)
               do j=2,Nz-1
                  roN2(i,j)=exp(-z(j)/Ho)/N2(i,j)
               end do
            end do              ! loop on y direction
c     smooth U 2 times in y direction
            do j=1,Nz
               do jj=1,2
                  do i=2,Ny-1
                     U(i,j)=0.25*U(i-1,j)+0.5*U(i,j)+0.25*U(i+1,j)
                  end do
               end do
            end do
c     
c     calculate PV gradient
c     
c     quality check- set N2=999 for bad data values
            do i=1,ny                 
               do j=1,nz
                  if(U(i,j).lt.-180.0.or.U(i,j).gt.180.0) N2(i,j)=999.
                  if(N2(i,j).lt.0.0) N2(i,j)=999.
               end do
            end do
            do j=1,nz
               call rzderv(U,z,Ho,roN2,i,j,delz,ny,nz,fz,Uznz)
               qy(1,j)=0.
               do i=2,ny
                  call rzderv(U,z,Ho,roN2,i,j,delz,ny,nz,fz,Uznz)
                  call ryderv(U,i,j,dely,ny,nz,Uy,Uyy)            
                  qy(i,j)=beta(i)+1./ae/ae*(
     $                 -Uyy+Uy*tan(y(i))+U(i,j)/cos(y(i))/cos(y(i)))-
     $                 fcor(i)*fcor(i)*Uznz
               end do
            end do
c     
c     write basic state to files
c     
c            do j=1,nz
c               write(71,*)(U(i,j),i=1,ny)
c               write(72,*)(N2(i,j),i=1,ny)
c               write(73,*)(qy(i,j),i=1,ny)
c            end do
            print*,U(1,10),U(10,1)
            write(71,rec=irec)U
c            print*,U
            write(72,rec=irec)N2
            write(73,rec=irec)qy
            print*,'wrote'
c     print*,U(10,10),qy(10,10),N2(10,10)
c     call writeout(U,ny,nz,'outdir/wind.out                          ')
c     call writeout(N2,1,nz,'outdir/nsquared.out                      ')
c     call writeout(qy,ny,nz,'outdir/pvgrad.out                       ')
 777        continue
         end do                 !loop on day
         close(41)
         close(42)
         close(71)
         close(72)
         close(73)
      end do                    !loop on year
      end

c     --------------------------------------------------------------
c
c     subroutine writeout write a variable to a file. input variable,
c     its dimentions and file name
c
      subroutine writeout(var,n1,n2,filename)
c
      real var(n1,n2)
      integer n1,n2
      character*25 filename
      open(77,file=filename)

      do j=1,n2
         write(77,*)(var(i,j),i=1,n1)
      end do
      close(77)
c
      return
      end
c     ---------------------------------------------------------------
c
c     subroutine writeout write a variable to a file. input variable,
c     its dimentions and file name
c
      subroutine rwriteout(var,n1,n2,filename)
c
      complex var(n1,n2)
      integer n1,n2
      character*25 filename
      open(77,file=filename)

      do j=1,n2
         write(77,*)(var(i,j),i=1,n1)
      end do
      close(77)
c
      return
      end
c     ---------------------------------------------------------------
c     subroutine ryderv calculates y derivatives of a real fuction
c
      subroutine ryderv(f,i,j,dely,ny,nz,fy,fyy)
c
      real dely,f(ny,nz),fy,fyy
c
      if (i.eq.1) then
         fy=(f(i+1,j)-f(i,j))/dely
         fyy=(f(i+2,j)-2.*f(i+1,j)+f(i,j))/dely/dely
      else if (i.eq.ny) then
         fy=(f(i,j)-f(i-1,j))/dely
         fyy=(f(i-2,j)-2.*f(i-1,j)+f(i,j))/dely/dely
      else
         fy=(f(i+1,j)-f(i-1,j))/2./dely
         fyy=(f(i+1,j)-2.*f(i,j)+f(i-1,j))/dely/dely   
      end if
c
      return
      end
c     ----------------------------------------------------------
c     subroutine rzderv calculates the z derivatives of a real function
c     taking into account density and N2 effects...
c
      subroutine rzderv(f,z,Ho,roN2,i,j,delz,ny,nz,fz,fznz)
c
      real delz,Ho,roN2(ny,nz),z(nz),f(ny,nz),fz,fznz
c
      if (j.eq.1) then
         fz=(f(i,j+1)-f(i,j))/delz
         fznz=exp(z(j+1)/Ho)/2./delz/delz*(
     $        (roN2(i,j+2)+roN2(i,j+1))*f(i,j+2)-
     $        (roN2(i,j+2)+2.*roN2(i,j+1)+roN2(i,j))*f(i,j+1)+
     $        (roN2(i,j+1)+roN2(i,j))*f(i,j))
      else if (j.eq.nz) then
         fz=(f(i,j)-f(i,j-1))/delz
         fznz=exp(z(j-1)/Ho)/2./delz/delz*(
     $        (roN2(i,j)+roN2(i,j-1))*f(i,j)-
     $        (roN2(i,j)+2.*roN2(i,j-1)+roN2(i,j-2))*f(i,j-1)+
     $        (roN2(i,j-1)+roN2(i,j-2))*f(i,j-2))
      else
         fz=(f(i,j+1)-f(i,j-1))/2./delz
         fznz=exp(z(j)/Ho)/2./delz/delz*(
     $        (roN2(i,j+1)+roN2(i,j))*f(i,j+1)-
     $        (roN2(i,j+1)+2.*roN2(i,j)+roN2(i,j-1))*f(i,j)+
     $        (roN2(i,j)+roN2(i,j-1))*f(i,j-1))
      end if
c      
      return
      end
c     ------------------------------------------------------------
c     subroutine interp- given a matrix f1 of which the SECOND 
c     dimention is m, and the FIRST is defined on a vector z1 that
c     is an n1 size array (i.e. m different profiles...)
c     interpolate onto a new set of z points, denoted by z2(n2),
c     for each of the m profiles. The output is f2.
c     uses a bicubic spline interpolation.
c     input: m,n1,z1(n1),f1(n1,m),n2,z2(m2) 
c     output:f2(n2,m),nerr
c     nerr is an error parameter.
c     call interp(n1,z1,f1,n2,z2,f2,nerr)
c     
      subroutine interp(n1,z1,f1,n2,z2,f2,nerr)
c
c     i goes with the m index. n,j go with z direction.
c
      integer n1,n2,jln,j,nerr
      real 
     $     f1(n1),f2(n2),z1(n1),z2(n2)
      parameter(mmax=100,nmax=100)
      real f2derv(nmax)
c     if the point is out of the range, it uses the last 2*i(j)block
c     grid points to extraplate to the new point.the main program 
c     should then quit the integration.
      nerr=0
      call spline(z1,f1,n1,1.e30,1.e30,f2derv)
      jln=0
      do 32 j=1,n2           
         call splint(z1,f1,f2derv,n1,z2(j),f2(j))
         call hunt(z1,n1,z2(j),jln)
         if(f2(j).eq.min(f1(jln),f1(jln+1),f2(j))
     $        .and.jln.ne.0.and.jln.ne.n1) then
            nerr=nerr+1
c            print*,'error min',jln,z1(jln)
c            print*,f1(jln),f1(jln+1),f2(j)
         else if(f2(j).eq.max(f1(jln),f1(jln+1),f2(j))
     $           .and.jln.ne.0.and.jln.ne.n1) then
            nerr=nerr+1
c            print*,'error max',jln, z1(jln)
c            print*,f1(jln),f1(jln+1),f2(j)
         end if
 32   continue
c      if (nerr.ne.0) print*,'nerr=',nerr
      return
      end
c     -----------------------------------------------------------
c     Given arrays X and Y of length N containinf a tabulated function, 
c     i.e. Yi=f(xi), with x1<x2...<xN, and given values YP1 and YPN for 
c     the first derivative of the interpolating function at points 1 and N,
c     respectively, this routine returns an array Y2 of length N which 
c     contains the second derivatives of the intepolating function at the
c     tabualted points Xi. If YP1 adn/or YPN are equal to 1x10E30 or larger, 
c     the routine is signalled to set the corresponding boundary condition
c     for a natural spline, with zero second derivative on that boundary.
c     
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+
     $        1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     $        u(i-1))/p
 11   continue
      if (ypn.gt..99e30) then
         qn=0.
         un=0.
      else
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 12   continue
      return
      END
c     -----------------------------------------------------------
c     subroutine splint- given arrays XA and YA of length N whcih tabulate
c     a function (with the XAi's is order), and given the array Y2A, which 
c     is the output from SPLINE above, and given a value of X, this routine
c     returns a cubic-spline interpolated value Y.
c
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     $     2)/6.
      return
      END
c     -----------------------------------------------------------
c     subroutine hunt- gien an array XX of length N, and given a value X,
c     returns a value JLO such that X is between XX(JLO) and XX(JLO+1). 
c     XX must be a monotonic, either increasing or decreasing. JLO=0 or JLO=N
c     is returned to indicate that X is out or range. JLO on input is taken 
c     as the initial guess for JLO on output.
c
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1       jhi=jlo+inc
         if(jhi.gt.n)then
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
 2       jlo=jhi-inc
         if(jlo.lt.1)then
            jlo=0
         else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
      END
c     -----------------------------------------------------------

c     subroutine twodint- given a function f1 on a 2-D grid that is defined
c     by y1, z1 which are m1, n1 size arrays, respectively interpolate
c     onto a new grid
c     y2,z2 (m2xn2 matrices). The output is f2.
c     uses a bicubic spline interpolation.
c     input: m1,n1,y1(m1),z1(n1),f1(m1,n1),
c            m2,n2,y2(m2,n2),z2(m2,n2) 
c     output:f2(m2,n2),nerr
c     nerr is an error parameter.
c     call twoint(ny,nz,y,z,f1,ns,nr,ys,zs,f2,nerr)
c     call twodint(ny,nz,y,z,f1,nr,ns,yr,zr,f2,nerr)
c
      subroutine twodint(m1,n1,y1,z1,f1,m2,n2,y2,z2,f2,nerr)
c
c     m,i go with y direction. n,j go with z direction.
c
      integer mm,nn,m1,n1,m2,n2,jlm,jln,i,j,nerr
      real 
     $     f1(m1,n1),f2(m2,n2),y1(m1),z1(n1),
     $     y2(m2),z2(n2)
      parameter(mmax=100,nmax=100)
      real f2derv(mmax,nmax)
c     if the point is out of the range, it uses the last 2*i(j)block
c     grid points to extraplate to the new point.the main program 
c     should then quit the integration.
      nerr=0
      call splie2(y1,z1,f1,m1,n1,f2derv)
      jlm=0
      jln=0
      do 31 i=1,m2
         do 32 j=1,n2           
            call splin2(y1,z1,f1,f2derv,m1,n1,
     $           y2(i),z2(j),f2(i,j))
            call hunt(y1,m1,y2(i),jlm)
            call hunt(z1,n1,z2(j),jln)
            if(f2(i,j).eq.min(f1(jlm,jln+1),f1(jlm+1,jln+1),
     $           f1(jlm,jln),f1(jlm+1,jln),f2(i,j)).and.jlm.ne.0.and.
     $           jlm.ne.m1.and.jln.ne.0.and.jln.ne.n1) then
               nerr=nerr+1
            else if(f2(i,j).eq.max(f1(jlm,jln+1),f1(jlm+1,jln+1),
     $              f1(jlm,jln),f1(jlm+1,jln),f2(i,j)).and.
     $              jlm.ne.0.and.jlm.ne.m1.and.jln.ne.0.and.
     $              jln.ne.n1) then
               nerr=nerr+1
            end if
 32      continue
 31   continue
      return
      end
      
c     ---------------------------------------------------------
c     subroutine splie2- given an M by N tabulated function YA, and 
c     tabulated independent variables X1A (N values), this routine 
c     constructs one dimensioanl natuaral cubic splines of the rows of 
c     YA and returns the second derivatives in the array Y2A
c      
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
C     USES spline
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN)
      do 13 j=1,m
         do 11 k=1,n
            ytmp(k)=ya(j,k)
 11      continue
         call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
         do 12 k=1,n
            y2a(j,k)=y2tmp(k)
 12      continue
 13   continue
      return
      END
c     -----------------------------------------------------------
c     subroutine splin2- given X1A, X2A, YA,M,N as described in SPLIE2,
c     and Y2A as produced by that routine; and given a desired interpolating
c     point X1, X2; this routine returns an interpolated funcion value Y by
c     bicubic spline interpolation.
c
      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
C     USES spline,splint
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
         do 11 k=1,n
            ytmp(k)=ya(j,k)
            y2tmp(k)=y2a(j,k)
 11      continue
         call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
 12   continue
      call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
c     -----------------------------------------------------------

c       if(nysponge.gt.0) then
c            do i=nyobs-nysponge+1,nyobs+10
c               Uobs(i,j)=Uobs(nyobs-nysponge,j)
c            end do
c         end if





