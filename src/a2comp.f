      program anywave
c
c     spherical coordiante model using Plumb's formulation 
c     (streamfunction=geopotential height/coriolis force)
c     poisson eq. solver - 
c     driver for blktoc.f subroutine
c
      implicit integer(i-n)
      include 'wave1.h'
      real
     $     width,Ho,dumry,dumalz,alheight,alamp,
     $     roN2top(Ny),roN2bot(Ny),N2bot(Ny),N2top(Ny),
     $     alscale,alwidth,ae,om,top,bottom,s,phio,
     $     gammac,wide,dely,delz,amrad,ahrad,awrad
      real
     $     ampf(ny,nz), phasef(ny,nz), y(ny), z(nz), realf(ny,nz),
     $     N2(ny,nz),roN2(ny,nz),U(ny,nz),qy(ny,nz),alpha(ny,nz),
     $     rayd(ny,nz),imagf(ny,nz),fcor(Ny),beta(Ny),Fn2(ny,nz),
     $     n2ref(ny,nz),n2refdamp(ny,nz),aqprime(ny,nz),cph(ny,nz),
     $     rqprime(ny,nz),ly2(ny,nz),mz2(ny,nz),iqprime(ny,nz),
     $     rcheck2(ny,nz),icheck2(ny,nz),epfy(ny,nz),epfz(ny,nz),
     $     rtprime(ny,nz),itprime(ny,nz),endy,endz,alampy,alscaley
      complex
     $     phixx,phiyy,phiznz,phiy,phiz,qprime(ny,nz),aitop,aibot,
     $     psizz,check2(ny,nz),bibot(ny),bitop,damp,gamma(ny,nz)
      character*35 filename
      character dirdat*43,Ufile*52,N2file*53,PVfile*53
      character*54 phrfile,phifile,tprfile,tpifile,ly2file,mz2file
      character*54 n2rfile,epyfile,epzfile
      character*5 thld
c     character*4 year(24)
c     matrix solver definitions
      complex a(ny,nz,9),f(ny,nz)
      complex e(ny,ny),v(ny,2),w(ny,ny),ef(ny,ny,nz)
      integer ier,lp

c      data year
c     $     /'1979','1980','1981','1982','1983','1984','1985','1986',
c     $     '1987','1988','1989','1990','1991','1992','1993','1994',
c     $     '1995','1996','1997','1998','1999','2000','2001','2002'/

      open (11,file='wave3sp.dat',status='unknown')
      read(11,*)ae,om,top,bottom
      read(11,*)s,wide,freq,gammac
      read(11,*)alheight,alwidth,endy,endz
      read(11,*)dum1,dum2,dum3,dum4
      if(dum1.ne.ny.or.dum2.ne.nz) print*,'check wave3.dat, ny,nz'
      if(dum3.ne.ns.or.dum4.ne.nr) print*,'check wave3.dat, nr,ns'
      read(11,*)alampy,alscaley,amrad,ahrad
      read(11,*)awrad,alamp,alscale,Ho
      read(11,*)idays,iyrs,ny1,ny2
      if(idays.ne.ndobs) print*,'check wave3 nt1- number of days'
      close (11)
      print*,iyrs
      print*,'s,freq,wide,Ho:' 
      print*,s,freq,wide,Ho
      print*,'ae,om,top,bottom ',ae,om,top,bottom
      print*,'alheight,alwidth,alamp,alscale:'
      print*,alheight,alwidth,alamp,alscale
      print*,'endy, endz ',endy,endz
      print*,'alampy,alscaley',alampy,alscaley
      print*,'amrad,ahrad,awrad',amrad,ahrad,awrad
c     
c      print*,'enter hemisphere (SH: negative num, NH: positive or zero)'
c      read*,hem
c      print*,'enter threshold directory: e.g. thld2'
c      read*,thld
      hem=1   !Northern hemisphere
      
      if(hem.ge.0.) then        
         alscaley=abs(alscaley)*pi/180.
      else
         alscaley=-abs(alscaley)*pi/180.
      endif
      awrad=awrad*pi/180.
      endy=endy*pi/180.
      alwidth=alwidth*pi/180.
      width=wide*pi/180.
      top=top*Ho
      bottom=bottom*Ho
      alheight=alheight*Ho
      alscale=alscale*Ho
      endz=endz*Ho
      ahrad=ahrad*Ho
      alamp=alamp/86400.
      alampy=alampy/86400.
      amrad=amrad/86400.
      
      delz=(top-bottom)/float(nz)
      dely=width/float(ny-1)
c     
c     define y,z 
c
      if(hem.ge.0.) then        
         do i=1,ny
            y(i)=dely*float(i-1)
         end do
      else
         do i=1,ny
            y(i)=dely*float(i-1)-0.5*pi
         end do
      endif
      do j=1,nz
         z(j)=delz*float(j-1)+bottom
      end do
      do i=1,ny
         fcor(i)=2*om*sin(y(i))
         beta(i)=2*om*cos(y(i))/ae
      end do

c     get phase speed in m/sec from frequency, which is given in days^(-1)
      do j=1,nz
         do i=1,ny
            cph(i,j)=freq*ae*cos(y(i))/s/86400.0
         end do
      end do

      
c     bottom forcing
c      open(66,file='outdir/fbottom.out')
c      read(66,*)bibot
c      close (66)     
c     specify 
c     
      do i=2,ny-1
         yi=float(i)
         yny=float(ny)
         bibot(i)=1.0
      end do
      call rwriteout(fbot,1,ny,'outdir/fbottom.out         ')
c
c     define the damping coefficients rayd (linear momentum damping 
c     coefficient) and alpha (newtonian damping coefficient)
c
      do j=1,nz
         do i=1,ny
            rayd(i,j)=0.5*alamp*(tanh((z(j)-alheight)/alscale)-
     $           tanh((z(1)-alheight)/alscale))
            if(hem.lt.0.) then
               alwidth=-abs(alwidth)
               rayd(i,j)=rayd(i,j)+alampy*0.5*
     $              (tanh((y(i)-alwidth)/alscaley)-
     $              tanh((y(1)-alwidth)/alscaley))
            else
               alwidth=abs(alwidth)
               rayd(i,j)=rayd(i,j)+alampy*0.5*
     $              (tanh((y(ny)-alwidth)/alscaley)-
     $              tanh((y(i)-alwidth)/alscaley))
            end if
            alpha(i,j)=rayd(i,j)+amrad*
     $           (0.45*exp(-(z(j)-50000.0)**2/(13000.0)**2)+
     $           0.05)
            gamma(i,j)=(0.,1.)*gammac*ae*cos(y(i))/s/86400.
         end do
      end do
c
      dirdat='./'
c     start loop on years
      do iyr=1,1
         print*,iyr 
c     basic state input files 
         Ufile  ='ERAU19792015NH.dat'
         N2file ='ERAN219792015NH.dat'
         PVfile ='ERAPV19792015NH.dat'        
         open(unit=22,file=Ufile,form='unformatted',access='direct',
     $        recl=4*ny*nz)             
         open(unit=33,file=N2file,form='unformatted',access='direct',
     $        recl=4*ny*nz)             
         open(unit=55,file=PVfile,form='unformatted',access='direct',
     $        recl=4*ny*nz)             
c     output files 
         phrfile='phrcat_obs.dat'
         phifile='phicat_obs.dat'
         tprfile='tprcat_obs.dat'        
         tpifile='tpicat_obs.dat'
         ly2file='ly2cat_obs.dat'
         mz2file='mz2cat_obs.dat'        
         n2rfile='n2rcat_obs.dat'
         epyfile='epycat_obs.dat'
         epzfile='epzcat_obs.dat'        

         open(unit=71,file=phrfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=72,file=phifile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=73,file=tprfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=74,file=tpifile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=75,file=ly2file,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=76,file=mz2file,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=77,file=n2rfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=78,file=epyfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
         open(unit=79,file=epzfile,form='unformatted',access='direct',
     $        status='replace',recl=4*ny*nz)             
c     
c     start loop on days (months/time steps)
         do iy=1,idays
            irec=iy
            read(22,rec=irec)U
            read(33,rec=irec)N2
            read(55,rec=irec)qy
c     print*,U(10,10),qy(10,10),N2(10,10)
c     
c     define basic state functions (N2,density effects)
c     
            do i=1,Ny
               do j=1,nz
                   roN2(i,j)=exp(-z(j)/Ho)/N2(i,j)
               end do
               N2bot(i)=N2(i,1)
               N2top(i)=N2(i,Nz)
               roN2bot(i)=exp(-(z(1)-delz)/Ho)/N2bot(i)
               roN2top(i)=exp(-(z(nz)+delz)/Ho)/N2top(i)
            end do
c     
c     set up coefficients for the matrix solver. 
c     i is "latitude" index, j is "height" index.
c     
            do j=1,nz
c     left boundary condition
               a(1,j,1)=(0.,0.)
               a(1,j,2)=(0.,0.)
               a(1,j,3)=(0.,0.)
               a(1,j,4)=(0.,0.)
               a(1,j,5)=(1.,0.)
               a(1,j,6)=(0.,0.)
               a(1,j,7)=(0.,0.)
               a(1,j,8)=(0.,0.)
               a(1,j,9)=(0.,0.)
c     right boundary condition
               a(ny,j,1)=(0.,0.)
               a(ny,j,2)=(0.,0.)
               a(ny,j,3)=(0.,0.)
               a(ny,j,4)=(0.,0.)
               a(ny,j,5)=(1.,0.)
               a(ny,j,6)=(0.,0.)
               a(ny,j,7)=(0.,0.)
               a(ny,j,8)=(0.,0.)
               a(ny,j,9)=(0.,0.)
            end do
c     
            do i=1,ny
               aibot=(0.,0.)
               aitop=(0.,0.)
c     aitop unescessary if we specify phi at lower boundary, so set to zero.
c     aitop is zero if there is no shear at the top boundary and a lid. 
c     lower boundary condition 
               a(i,1,1)=(0.,0.)
               a(i,1,2)=(0.,0.)
               a(i,1,3)=(0.,0.)
               a(i,1,4)=(0.,0.)
               a(i,1,5)=(1.,0.)
               a(i,1,6)=(0.,0.)
               a(i,1,7)=(0.,0.)
               a(i,1,8)=(0.,0.)
               a(i,1,9)=(0.,0.)
c     
c     top boundary condition- nz is the top boundary
               a(i,nz,1)=(0.,0.)
               a(i,nz,2)=(0.,0.)
               a(i,nz,3)=(0.,0.)
               a(i,nz,4)=(0.,0.)
               a(i,nz,5)=(1.,0.)
               a(i,nz,6)=(0.,0.)
               a(i,nz,7)=(0.,0.)
               a(i,nz,8)=(0.,0.)
               a(i,nz,9)=(0.,0.)
            end do
c     
c     initialize all f(i,j) - forcing - to zero
c     
            do i=1,ny
               do j=1,nz
                  if(i.eq.1.or.i.eq.ny.or.j.eq.1.or.j.eq.nz)then
                     f(i,j)=(0.,0.)
                  else
                     f(i,j)=(0.,0.)
                  end if
               end do
            end do
            
c     bottom forcing- wavemaker- 
c     make sure it is straemfunction, not geopotential heihgt.
            do i=2,ny-1
               f(i,1)=bibot(i)*100.0*go/fcor(i)
            end do

c     set all interior coeffs.
            do i=2,ny-1
               do j=2,nz-1
                  sn=sin(y(i))
                  cs=cos(y(i))
                  ez=exp(z(j)/Ho)
                  damp=(0.,1.)*s*(U(i,j)-cph(i,j)-gamma(i,j))/ae/cs
                  a(i,j,1)=(0.,0.)
                  a(i,j,2)=(roN2(i,j)+roN2(i,j-1))/2./delz/delz*
     $                 fcor(i)*fcor(i)*ez*(damp+alpha(i,j))-
     $                 fcor(i)*fcor(i)/N2(i,j)/2./delz*
     $                 (alpha(i,j+1)-alpha(i,j-1))/2./delz
                  a(i,j,3)=(0.,0.)
                  a(i,j,4)=(damp+rayd(i,j))/ae/ae/dely/dely+
     $                 (damp+rayd(i,j))/ae/ae*sn/cs/2./dely-
     $                 (rayd(i+1,j)-rayd(i-1,j))/2./dely/ae/ae/2./dely
                  a(i,j,5)=(0.,1.)*s/ae/cs*qy(i,j)-
     $                 s*s/ae/ae/cs/cs*(damp+rayd(i,j))-
     $                 (roN2(i,j+1)+2.*roN2(i,j)+roN2(i,j-1))/2./delz/
     $                 delz*
     $                 fcor(i)*fcor(i)*ez*(damp+alpha(i,j))-
     $                 2./dely/dely/ae/ae*(damp+rayd(i,j))
                  a(i,j,6)=(damp+rayd(i,j))/ae/ae/dely/dely-
     $                 (damp+rayd(i,j))/ae/ae*sn/cs/2./dely+
     $                 (rayd(i+1,j)-rayd(i-1,j))/2./dely/ae/ae/2./dely
                  a(i,j,7)=(0.,0.)
                  a(i,j,8)=(roN2(i,j+1)+roN2(i,j))/2./delz/delz*
     $                 fcor(i)*fcor(i)*ez*
     $                 (damp+alpha(i,j))+
     $                 fcor(i)*fcor(i)/N2(i,j)/2./delz*
     $                 (alpha(i,j+1)-alpha(i,j-1))/2./delz
                  a(i,j,9)=(0.,0.)
               end do
            end do
c     call main matrix solver. f is forcing on input and the solution
c     on output
c      print*,'calling matrix solver'
c            call blktoc(a,ny,nz,ny,nz,f,ny,e,v,w,44,ier,lp)
            call blktoc1(a,ny,nz,ny,nz,f,ny,e,ef,v,w,44,ier,lp)
c     **********************************************************
c     calculate diagnostics here- call varous subroutines
c     **********************************************************  
c     
c     Fn2 is the N2 and density contribution to the index of refraction

            call Fn2calc (z,N2,nz,ny,roN2,delz,Ho,roN2top,roN2bot,
     $           N2top,N2bot,Fn2) 
c     

            do i=1,ny
c              print*,'calculating lat:', i
               do j=1,nz
c     
c     calulate the PV perturbation     
c     
                  call fyderv(f,i,j,dely,ny,nz,phiy,phiyy)
                  call fzderv(f,z,Ho,N2,roN2,Fn2,i,j,delz,ny,nz,
     $                 phiz,phiznz,psizz)
                  call raydy(rayd,i,j,ny,nz,dely,dumry)
                  call alphaz(alpha,i,j,ny,nz,delz,dumalz)
c     
                  if (i.eq.1) then
                     qprime(i,j)=0.
                     check2(i,j)=(0.,0.)
                  else
                     phixx=-s*s*f(i,j)/ae/ae/cos(y(i))**2
                     qprime(i,j)=fcor(i)*fcor(i)*phiznz+
     $                    phixx+phiyy/ae/ae-tan(y(i))/ae/ae*phiy
                     damp=(0.,1.)*s/ae/cos(y(i))
                     check2(i,j)=qprime(i,j)*
     $                    (U(i,j)-cph(i,j)-gamma(i,j))*damp+
     $                    damp*qy(i,j)*f(i,j)+
     $                    rayd(i,j)*(phixx+phiyy/ae/ae-
     $                    tan(y(i))/ae/ae*phiy)+
     $                    alpha(i,j)*fcor(i)*fcor(i)*phiznz+
     $                    dumry/ae/ae*phiy+
     $                    fcor(i)*fcor(i)/N2(i,j)*dumalz*phiz
                  end if
                  rcheck2(i,j)=real(check2(i,j))
                  icheck2(i,j)=imag(check2(i,j))
c     
c     calculate temperature perturbation- dimensional
c     
                  rtprime(i,j)=real(phiz)*Ho/Ro*fcor(i)
                  itprime(i,j)=imag(phiz)*Ho/Ro*fcor(i)
c     
c     calculate wave numbers using second derivative of phi (ly2,mz2)
c     
                  call wavn2(f,phiyy,psizz,z,Ho,N2,y,phiy,i,j,ny,nz,
     $                 ly2(i,j),mz2(i,j))
c     
c     calculate EP fluxes
                  epfy(i,j)=exp(-z(j)/Ho)*s*0.5*
     $                 imag(conjg(f(i,j))*phiy)/ae
                  epfz(i,j)=roN2(i,j)*s*0.5*
     $                 imag(conjg(f(i,j))*phiz)*fcor(i)**2
c     
c     calculate index of refraction squared
c     (n2ref, n2refdamp) with and without damping terms
c     
                  if(i.eq.1) then
                     n2ref(i,j)=0.
                  else
                     dumqy=qy(i,j)/(U(i,j)-cph(i,j))*ae*ae
                     dums2=-s*s/cos(y(i))**2
                     dumFn2=Fn2(i,j)*ae*ae*fcor(i)*fcor(i)
                     n2ref(i,j)=dumqy+dums2+dumFn2
c     n2ref(i,j)=n2ref(i,j)/ae/ae/fcor(i)/fcor(i)
                  end if
               end do
            end do 
c     
c     write results to output files
c     convert to geopotential height by dividing by fcor(i)
c     
            do i=1,ny
               do j=1,nz
                  ampf(i,j)=sqrt(real(f(i,j))**2+imag(f(i,j))**2)
                  imagf(i,j)=imag(f(i,j)*fcor(i))
                  realf(i,j)=real(f(i,j)*fcor(i))
                  rqprime(i,j)=real(qprime(i,j))
                  iqprime(i,j)=imag(qprime(i,j))
                  aqprime(i,j)=sqrt(real(qprime(i,j))**2+
     $                 imag(qprime(i,j))**2)
                  if(f(i,j).ne.0.)then
                     phasef(i,j)=atan(imag(f(i,j))/real(f(i,j)))
                  else
                     phasef(i,j)=0.
                  end if
               end do
            end do            

            write(71,rec=irec)realf
            write(72,rec=irec)imagf
            write(73,rec=irec)rtprime
            write(74,rec=irec)itprime
            write(75,rec=irec)ly2
            write(76,rec=irec)mz2
            write(77,rec=irec)n2ref
            write(78,rec=irec)epfy
            write(79,rec=irec)epfz
         end do                 !loop on day
         print*, iy
      end do                    !loop on year
      print*,'put negative for evanescence regions, in ly, mz'
      call writeout(alpha,ny,nz,'outdir/alpha.out                 ')
      call writeout(rayd,ny,nz,'outdir/rayd.out                 ')
c      call writeout(n2refdamp,ny,nz,'outdir/n2refdamp.out')
c
      close(22)
      close(33)
      close(55)
      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      close(79)
      if(ier.ne.0) print*,'An error occurred (ier != 0)'
      stop 
      end

c**************************************************************
c     subroutines to calculate varous diagnostics
c**************************************************************
      subroutine Fn2calc (z,N2,nz,ny,roN2,delz,Ho,roN2top,roN2bot,
     $     N2top,N2bot,Fn2)
c
c     calculates Fn2, the contribution of N2 adn density to the index
c     of refracton.
c      
      real Fn2(ny,nz),roN2(ny,nz),N2(ny,nz),z(nz)
      real Ho,delz,roN2bot(ny),roN2top(ny),N2top(ny),N2bot(ny)
      do i=1,Ny
         do j=1,nz
            if (j.eq.1) then
               dumcp1=exp(z(j+1)/Ho/2.)*sqrt(N2(i,j+1))
               dumc0=exp(z(j)/Ho/2.)*sqrt(N2(i,j))
               dumcm1=exp((z(j)-delz)/Ho/2.)*sqrt(N2bot(i))
               duma0=dumc0/N2(i,j)
               dumbp1=roN2(i,j+1)
               dumb0=roN2(i,j)
               dumbm1=roN2bot(i)
            else if (j.eq.nz) then
               dumcp1=exp((z(j)+delz)/Ho/2.)*sqrt(N2top(i))
               dumc0=exp(z(j)/Ho/2.)*sqrt(N2(i,j))
               dumcm1=exp(z(j-1)/Ho/2.)*sqrt(N2(i,j-1))
               duma0=dumc0/N2(i,j)
               dumbp1=roN2top(i)
               dumb0=roN2(i,j)
               dumbm1=roN2(i,j-1)
            else
               dumcp1=exp(z(j+1)/Ho/2.)*sqrt(N2(i,j+1))
               dumc0=exp(z(j)/Ho/2.)*sqrt(N2(i,j))
               dumcm1=exp(z(j-1)/Ho/2.)*sqrt(N2(i,j-1))
               duma0=dumc0/N2(i,j)
               dumbp1=roN2(i,j+1)
               dumb0=roN2(i,j)
               dumbm1=roN2(i,j-1)
            end if
            Fn2(i,j)=duma0/2./delz/delz*(
     $           (dumbp1+dumb0)*dumcp1-
     $           (dumbp1+2.*dumb0+dumbm1)*dumc0+
     $           (dumb0+dumbm1)*dumcm1)
         end do

         Fn2(i,1)=Fn2(i,2)
      end do

c
      return
      end
c     --------------------------------------------------------
c     subroutine fyderv calculates y derivatives of f (phi)
c
      subroutine fyderv(f,i,j,dely,ny,nz,phiy,phiyy)
c
      real dely
      complex f(ny,nz),phiy,phiyy
c
      if (i.eq.1) then
         phiy=(f(i+1,j)-f(i,j))/dely
         phiyy=(f(i+2,j)-2.*f(i+1,j)+f(i,j))/dely/dely
      else if (i.eq.ny) then
         phiy=(f(i,j)-f(i-1,j))/dely
         phiyy=(f(i-2,j)-2.*f(i-1,j)+f(i,j))/dely/dely
      else
         phiy=(f(i+1,j)-f(i-1,j))/2./dely
         phiyy=(f(i+1,j)-2.*f(i,j)+f(i-1,j))/dely/dely   
      end if
c
      return
      end
c     ----------------------------------------------------------
c     subroutine raydy calculates the y derivative of rayd.
c
      subroutine raydy(rayd,i,j,ny,nz,dely,dumry)
c
      real dely,dumry,rayd(ny,nz)
c      
      if (i.eq.1) then
         dumry=(rayd(i+1,j)-rayd(i,j))/dely
      else if (i.eq.ny) then
         dumry=(rayd(i,j)-rayd(i-1,j))/dely
      else
         dumry=(rayd(i+1,j)-rayd(i-1,j))/2./dely
      end if
c      
      return
      end
c     -----------------------------------------------------------
c     subroutine fzderv calculates the z derivatives of f (phi) 
c
      subroutine fzderv(f,z,Ho,N2,roN2,Fn2,i,j,delz,ny,nz,
     $     phiz,phiznz,psizz)
c
      real delz, Ho
      real N2(ny,nz),Fn2(ny,nz),z(nz),roN2(ny,nz)
      complex f(ny,nz),phiz,phiznz,psizz
c
      if (j.eq.1) then
         phiz=(f(i,j+1)-f(i,j))/delz
         phiznz=exp(z(j+1)/Ho)/2./delz/delz*(
     $        (roN2(i,j+2)+roN2(i,j+1))*f(i,j+2)-
     $        (roN2(i,j+2)+2.*roN2(i,j+1)+roN2(i,j))*f(i,j+1)+
     $        (roN2(i,j+1)+roN2(i,j))*f(i,j))
      else if (j.eq.nz) then
         phiz=(f(i,j)-f(i,j-1))/delz
         phiznz=exp(z(j-1)/Ho)/2./delz/delz*(
     $        (roN2(i,j)+roN2(i,j-1))*f(i,j)-
     $        (roN2(i,j)+2.*roN2(i,j-1)+roN2(i,j-2))*f(i,j-1)+
     $        (roN2(i,j-1)+roN2(i,j-2))*f(i,j-2))
      else
         phiz=(f(i,j+1)-f(i,j-1))/2./delz
         phiznz=exp(z(j)/Ho)/2./delz/delz*(
     $        (roN2(i,j+1)+roN2(i,j))*f(i,j+1)-
     $        (roN2(i,j+1)+2.*roN2(i,j)+roN2(i,j-1))*f(i,j)+
     $        (roN2(i,j)+roN2(i,j-1))*f(i,j-1))
      end if
      psizz=exp(-z(j)/Ho*0.5)*sqrt(N2(i,j))*(phiznz-
     $     Fn2(i,j)*f(i,j))
c      
      return
      end
c     ---------------------------------------------------------
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
c     subroutine alphaz calculates the z derivative of alpha
c
      subroutine alphaz(alpha,i,j,ny,nz,delz,dumalz)
c
      real delz, dumalz
      real alpha(ny,nz)
c
      if (j.eq.1) then
         dumalz=(alpha(i,j+1)-alpha(i,j))/delz
      else if (j.eq.nz) then
         dumalz=(alpha(i,j)-alpha(i,j-1))/delz
      else
         dumalz=(alpha(i,j+1)-alpha(i,j-1))/2./delz               
      end if
c
      return
      end
c     ------------------------------------------------------------
c     subroutine wavn2 calulates teh meridional and vertical wavenumbers
c     ly2,mz2, respectively, using the second derivative of phi.
c
      subroutine wavn2(f,phiyy,psizz,z,Ho,N2,y,phiy,i,j,ny,nz,ly2,mz2)
c
      real N2(ny,nz),z(nz),ly2,mz2,Ho,y(Ny)
      complex f(ny,nz)
      complex phiyy,psizz,phiy
c
      if (f(i,j).eq.0.) then
         ly2=0.
         mz2=0.
      else
         ly2=real((phiyy-tan(y(i))*phiy)/f(i,j))
         mz2=real(psizz/f(i,j))*exp(z(j)/Ho*0.5)*sqrt(N2(i,j))
      end if
      if (ly2.gt.0.) then
c         ly2=0
         ly2=-sqrt(ly2)
      else
         ly2=sqrt(-ly2)
      end if
      if (mz2.gt.0.) then
c         mz2=0
         mz2=-sqrt(mz2)
      else
         mz2=sqrt(-mz2)
      end if
c
      return
      end
c     --------------------------------------------------------------
c     .......... COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI) .........
c
      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)

c      DOUBLE PRECISION AR,AI,BR,BI,CR,CI
c      DOUBLE PRECISION S,ARS,AIS,BRS,BIS
      real AR,AI,BR,BI,CR,CI
      real S,ARS,AIS,BRS,BIS

c      S = DABS(BR) + DABS(BI)
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
c      S = BRS**2.D0 + BIS**2.D0
      S = BRS**2. + BIS**2.
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S


      RETURN
      END
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
c     ---------------------------------------------------------------
         SUBROUTINE BLKTOC1 (A, M, N, IA1, IA2, F, IFF, E,EF, 
     $     V, W,LU,IER, LP)
C ********************************************************************
C		MODIFIED BY NILI HARNIK, DEC 1997, TO USE MEMORY INSTEAD OF
C		WRITING TO HARD DISK AND BACK READING. ADDED EF(M,M,N) VARIABLE
C		TO CALLING SEQUENCE. DON'T NEED TO OPEN FILE OF UNIT LU, BUT
C		THE LINES USING IT ARE STILL COMMENTED OUT. 
C	
C        tridiagonal solvers documentation.
C

C        this routine solves the following system for f(j)
C

C                        b(j) * f(j) + c(j) * f(j+1) = d(j), j = 1,
C        a(j) * f(j-1) + b(j) * f(j) + c(j) * f(j+1) = d(j), j = 2, n-1
C        a(j) * f(j-1) + b(j) * f(j)                       , j = n,
C        where each a(j), b(j) and c(j) is a m by m tridiagonal matrix
C        and each d(j) and f(j) is an m-vector.
C

C        an extension of the thomas algorithm is used.
C        refs. richtmyer and morton, 1957
C              lindzen, r.s. and kuo, h-l. (1969). mon. wea. rev., 97,
C              pp732-4.
C

C        only 3 diagonals of a, b and c are stored, so that
C        a(i,1) is in fact a(i,i-1)
C        a(i,2) is in fact a(i,i  )
C        a(i,3) is in fact a(i,i+1).
C

C        on input, d is stored in array f. on output f contains the
C        solution.
C        e, v and w are working arrays.
C

C        ier is returned as non-zero if an error occurs.

C


C        almost all matrix operations are carried out by subsidiary

C        routines named, for example, mcxtf.

C        the first letter is m for matrix

C        the second letter is c for complex

C        the third letter denotes the operation, and is one of a, s, m,

C           x for +, -, *, and transfer.

C        the fourth and fifth letters denote the form of the first and

C        second matrix operands respectively. the letters can be t, f

C        or v for tridiagonal matrix, full matrix and vector.

C        the sixth character, if present, indicates that the third

C        vector argument is to be identical to the first or second

C        argument, and is therefore omitted.


C        only one of the e(j), j = 1, n-1 is stored in memory at any
C        time. the other e(j) are written sequentially to file lu, and
C        read backwards when required in the back substitution phase.
C

C

C        in this version, the coefficients are stored in the single
C        array a in the form generated by an elliptic fde, i.e.
C

C     a(i,j,7)*phi(i-1,j+1)+a(i,j,8)*phi(i,j+1)+a(i,j,9)*phi(i+1,j+1) +
C     a(i,j,4)*phi(i-1,j  )+a(i,j,5)*phi(i,j  )+a(i,j,6)*phi(i+1,j  ) +
C     a(i,j,1)*phi(i-1,j-1)+a(i,j,2)*phi(i,j-1)+a(i,j,3)*phi(i+1,j-1) =
C     f(i,j), i = 1, m,  j = 1, n.
C
C     EXTERNAL CALLS:   To IMSL routines LEQT1C
C                                        UERTST
C

C *********************************************************************
      COMPLEX A(IA1,IA2,9), F(IFF,N)
      COMPLEX E(M,M),EF(M,M,N),V(M,2), W(M,M)
	IA1IA2 = IA1 * IA2
	IJOB = 0
      J = 1
      N1 = N - 1
      CALL MCXTF (A(1,J,4), W, M, IA1IA2, M)
      CALL MCXTF (A(1,J,7), E, M, IA1IA2, M)
      CALL LEQT1C (W, M, M, E, M, M, IJOB, V, IER)
      IF (IER. NE. 0) GOTO 30
	do iy1=1,M
		do iy2=1,M
			EF(iy1,iy2,J)=E(iy1,iy2)
		end do
	end do
c		WRITE (LU) E
		CALL LEQT1C (W, M, M, F(1,1), 1, IFF, 2, V, IER)
		DO 10 J = 2, N
			CALL MCMTF (A(1,J,1), E, W, M, M, IA1IA2, M, M)
			CALL MCSTF2 (A(1,J,4), W, M, IA1IA2, M)
			IF (J .NE. N)
     $		CALL MCXTF (A(1,J,7), E, M, IA1IA2, M)
			IF (J .EQ. N) IJOB = 1
			CALL LEQT1C (W, M, M, E, M, M, IJOB, V, IER)
			IF (IER .NE. 0) GOTO 30
		IF (J .LT. N1) THEN 
			do iy1=1,M
				do iy2=1,M
					EF(iy1,iy2,J)=E(iy1,iy2)
				end do
			end do
		end if	
c		IF (J .LT. N1) THEN 
c				WRITE (LU) E
c			END IF	
			CALL MCMTV (A(1,J,1), F(1,J-1), V(1,2), M, IA1IA2)
			CALL MCSVV1 (V(1,2), F(1,J), M)
			CALL LEQT1C (W, M, M, F(1,J), 1, IFF, 2, V, IER)
10		CONTINUE
			DO 20 JJ = 1, N1
			J = N - JJ
			CALL MCMFV (E, F(1,J+1), V, M, M, M)
			CALL MCSVV1 (V, F(1,J), M)
			IF (J.EQ.1) GOTO 20 
			do iy1=1,M
				do iy2=1,M
					E(iy1,iy2)=EF(iy1,iy2,J-1)
				end do
			end do
c			BACKSPACE LU
c			read (LU, END=1000) E
c1000			continue
c			BACKSPACE LU
20		CONTINUE
	RETURN
30	CONTINUE
      CALL CPRNT (W, M, M, 1, M, M, LP, 10HL-U DECOMP)
      CALL CPRNT (V, M, 1, 1, M, 1, LP, 10HPIVOT INDS)
c      CALL UERRV1 ('BLKTOC1, 1 1 4 1 2, J, 1101,' ')
      STOP
	END

c---------


         SUBROUTINE MCMFV (A, B, C, M, N, IA)
C
C        this routine multiplies a complex full matrix a by a complex
C        vector b and stores the result in c.
C
      COMPLEX A(IA,N), B(N), C(M)
      DO 20 I = 1,M
         C(I) = 0.0
         DO 10 K = 1, N
            C(I) = C(I) + A(I,K) * B(K)
10       CONTINUE
20    CONTINUE
         RETURN
         END
         SUBROUTINE MCMTF (A, B, C, M, N, IA, IB, IC)
C
C        this routine multiplies a complex tridiagonal matrix a by a
C        complex full matrix b and stores the result in c.
C        only the diagonals of a are stored, so that
C        a(i,1) is in fact a(i,i-1)
C        a(i,2) is in fact a(i,i  )
C        a(i,3) is in fact a(i,i+1).
C
      COMPLEX A(IA,3), B(IB,N), C(IC,N)
      M1 = M - 1
      DO 20 J = 1, N
         I = 1
         C(I,J) = A(I,2) * B(I,J) + A(I,3) * B(I+1,J)
         DO 10 I = 2, M1
           C(I,J) = A(I,1) * B(I-1,J) + A(I,2) * B(I,J) +
     1               A(I,3) * B(I+1,J)
10       CONTINUE
         I = M
         C(I,J) = A(I,1) * B(I-1,J) + A(I,2) * B(I,J)
20    CONTINUE
         RETURN
         END
         SUBROUTINE M C M T V (A, B, C, M, IA)
C
C        this routine multiplies a complex tridiagonal matrix a by a
C        complex vector b and stores the result in c.
C        only the diagonals of a are stored, so that
C        a(i,1) is in fact a(i,i-1)
C        a(i,2) is in fact a(i,i  )
C        a(i,3) is in fact a(i,i+1).
C
      COMPLEX A(IA,3), B(M), C(M)
      C(1) = A(1,2) * B(1) + A(1,3) * B(2)
      M1 = M - 1
      DO 10 I = 2, M1
         C(I) = A(I,1) * B(I-1) + A(I,2) * B(I) + A(I,3) * B(I+1)
10    CONTINUE
      I = M
      C(I) = A(I,1) * B(I-1) + A(I,2) * B(I)
         RETURN
         END
         SUBROUTINE M C X T F (A, B, M, IA, IB)
C
C        this routine transfers a complex tridiagonal matrix a into a
C        complex full matrix b.
C        only the diagonals of a are stored, so that
C        a(i,1) is in fact a(i,i-1)
C        a(i,2) is in fact a(i,i  )
C        a(i,3) is in fact a(i,i+1).
C
      COMPLEX A(IA,3), B(IB,M)
      DO 20 J = 1, M
         DO 10 I = 1, M
            B(I,J) = 0.0
10       CONTINUE
20    CONTINUE
      B(1,1) = A(1,2)
      B(1,2) = A(1,3)
      M1 = M - 1
      DO 30 I = 2, M1
         B(I,I-1) = A(I,1)
         B(I,I  ) = A(I,2)
         B(I,I+1) = A(I,3)
30    CONTINUE
      I = M
      B(I,I-1) = A(I,1)
      B(I,I) = A(I,2)
         RETURN
         END
         SUBROUTINE M C S T F 2 (A, C, M, IA, IC)
C
C        this routine subtracts a complex full matrix c from a complex
C        tridiagonal matrix a and stores the result in c.
C        only the diagonals of a are stored, so that
C        a(i,1) is in fact a(i,i-1)
C        a(i,2) is in fact a(i,i  )
C        a(i,3) is in fact a(i,i+1).
C
      COMPLEX A(IA,3), C(IC,M)
      DO 20 J = 1, M
         DO 10 I = 1, M
            C(I,J) = -C(I,J)
10       CONTINUE
20    CONTINUE
      C(1,1) = C(1,1) + A(1,2)
      C(1,2) = C(1,2) + A(1,3)
      M1 = M - 1
      DO 30 I = 2, M1
         C(I,I-1) = C(I,I-1) + A(I,1)
         C(I,I  ) = C(I,I  ) + A(I,2)
         C(I,I+1) = C(I,I+1) + A(I,3)
30    CONTINUE
      I = M
      C(I,I-1) = C(I,I-1) + A(I,1)
      C(I,I  ) = C(I,I  ) + A(I,2)
         RETURN
         END
         SUBROUTINE M C S V V 1 (B, C, M)
C
C        this routine subtracts a complex vector b from a complex
C        vector c and stores the result in c.
C
      COMPLEX B(M), C(M)
      DO 10 I = 1, M
         C(I) = C(I) - B(I)
10    CONTINUE
         RETURN
         END
         SUBROUTINE CPRNT (A, M, N, KK, MD, ND, LP, LABEL)
C
C        this routine prints a complex array a of up to three dimensions
C
      COMPLEX A(MD,ND,KK)
      DATA MACR / 12 /
      MSECT = (M - 1) / MACR + 1
      WRITE (LP,196) LABEL
196   FORMAT(1H0,A10)
      DO 30 K = 1, KK
         IF (KK .GT. 1) WRITE (LP,199) K
199      FORMAT('0 K = ',I5)
         DO 20 IS = 1, MSECT
            IA = (IS - 1) * MACR + 1
            IB = MIN0 (IS * MACR, M)
            WRITE (LP,198) (I, I = IA, IB)
198         FORMAT(1H0,6X,12(I5,5X))
            DO 10 JJ = 1, N
               J = N + 1 - JJ
               WRITE (LP,197) J, (REAL(A(I,J,K)), I = IA, IB)
               WRITE (LP,197) J, (IMAG(A(I,J,K)), I = IA, IB)
197            FORMAT(1X,I5,1X,1P12G10.3)
10          CONTINUE
20       CONTINUE
30    CONTINUE
         RETURN
         END
      SUBROUTINE UERRV1 (NAME, CODE, VALUE, LENMSS, MESS)
C
C        this routine handles error and message processing.
C
C        for full details of its use, see the documentation common deck
C        errordoc.
C
      CHARACTER NAME*(*), MESS*(*)
      INTEGER CODE, VALUE, LENMSS
C
      CHARACTER STRING*137
      LOGICAL DAYFIL, OUTPUT, TRACE
      INTEGER SEVERE, VALTYP, CODX
      CHARACTER SEVLEV(4)*8
      INTEGER NCHARS, NCHARN, ERRNUM
      INTEGER NSPL, LMSPL(3), LENSPL(3)
      CHARACTER SPLMSS(3)*63, CHVAL*10
      DATA NSPL / 3 /
      DATA LENSPL(1), LENSPL(2), LENSPL(3) /
     1            20,        24,        26 /
      DATA SEVLEV(1), SEVLEV(2), SEVLEV(3), SEVLEV(4) /
     1       'NOTE.','WARNING.',   'STOP.',  'ABORT.' /
      DATA LMSPL(1), LMSPL(2), LMSPL(3) /
     1         1001,     1101,     1102 /
      DATA SPLMSS(1) / 'END OF FILE ON UNIT ' /
      DATA SPLMSS(2) / 'SINGULAR MATRIX FOR J = ' /
      DATA SPLMSS(3) / 'DENOM = 0.0 AT 1000*I+J = ' /
      CODX = ABS (CODE)
C
C        unravel codes.
      SEVERE = MIN0 (4, MOD (CODX / 100, 10))
      IF (SEVERE .EQ. 0) RETURN
      DAYFIL = MOD (CODX / 10000, 10) .EQ. 1
      OUTPUT = MOD (CODX /  1000, 10) .EQ. 1
      TRACE  = MOD (CODX /    10, 10) .EQ. 1
      VALTYP = MOD (CODX        , 10)
      ISPL = 0
      IF (LENMSS .LE. 1000) GO TO 20
C
C        look for special message.
      DO 10 I = 1, NSPL
         ISPL = I
         IF (LENMSS .EQ. LMSPL(I)) GO TO 20
10    CONTINUE
C
C        special message not found - abort.
      ISPL = 0
      SEVERE = 4
      TRACE = .TRUE.
      OUTPUT = .TRUE.
      DAYFIL = .TRUE.
20    CONTINUE
C
C        construct the full message in string.
      IF (.NOT. (DAYFIL .OR. OUTPUT)) GO TO 30
      STRING = ' ' // NAME // ' ' // SEVLEV(SEVERE) // ' '
      NCHARS = 3 + LEN (NAME) + LEN (SEVLEV(SEVERE))
      IF (ISPL .EQ. 0) THEN
         NCHARN = NCHARS + LEN (MESS)
         STRING(NCHARS+1:NCHARN) = MESS
      ELSE
         NCHARN = NCHARS + LENSPL(ISPL)
         STRING(NCHARS+1:NCHARN) = SPLMSS(ISPL)
      ENDIF
      NCHARS = NCHARN
      IF (VALTYP .NE. 0) THEN
         WRITE (CHVAL,'(I10)') VALUE
         STRING(NCHARS+1:NCHARS+10) = CHVAL
         NCHARS = NCHARS + 10
      ENDIF
      IF (DAYFIL) WRITE(*,*) STRING(:NCHARS)
      IF (OUTPUT) STRING(:1) = '0'
30    CONTINUE
      IF (.NOT. OUTPUT) THEN
         NCHARS = 1
         STRING(:1) = ' '
      ENDIF
      ERRNUM = 51
      IF (SEVERE .GE. 4) ERRNUM = 52
      IF (TRACE .OR. SEVERE .GE. 4)
     1   CALL SYSTEM ('echo '//STRING(:NCHARS))
      IF (OUTPUT .AND. .NOT. TRACE) PRINT '(A)', STRING(:NCHARS)
      IF (SEVERE .EQ. 3) STOP
      RETURN
      END

      
c     -----------------------------------------------------------------
      SUBROUTINE LEQT1C (A,N,IA,B,M,IB,IJOB,WA,IER) 
      INTEGER            N,IA,M,IB,IJOB,IER 
      COMPLEX            A(IA,N),B(IB,M)
      REAL               WA(N)
      REAL               P,Q,ZERO,ONE,T(2),RN,BIG 
      COMPLEX            SUM,TEMP 
      INTEGER            I,J,JM1,IM1,K,IMAX,JP1,IW,N1 
      EQUIVALENCE        (SUM,T(1)) 
      DATA               ZERO/0.0/,ONE/1.E0/
      IER = 0 
      IF (IJOB .EQ. 2) GO TO 75 
      RN = N
      DO 10 I=1,N 
         BIG = ZERO 
         DO 5 J=1,N 
            TEMP = A(I,J) 
             P = ABS(TEMP)
            IF (P .GT. BIG) BIG = P 
    5    CONTINUE 
         IF (BIG .EQ. ZERO) GO TO 105 
         WA(I) = ONE/BIG
   10 CONTINUE
      DO 70 J = 1,N 
         JM1 = J-1
         IF (JM1 .LT. 1) GO TO 25 
         DO 20 I=1,JM1
            SUM = A(I,J)
            IM1 = I-1 
            IF (IM1 .LT. 1) GO TO 20
            DO 15 K=1,IM1 
               SUM = SUM-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J) = SUM
   20    CONTINUE 
   25    P = ZERO 
         DO 45 I=J,N
            SUM = A(I,J)
            IF (JM1 .LT. 1) GO TO 40
            DO 35 K=1,JM1 
               SUM = SUM-A(I,K)*A(K,J)
   35       CONTINUE
            A(I,J) = SUM
   40       Q = WA(I)*ABS(SUM) 
            IF (P .GE. Q) GO TO 45
            P = Q 
            IMAX = I
   45    CONTINUE 
         Q = RN+P 
         IF (Q .EQ. RN) GO TO 105 
         IF (J .EQ. IMAX) GO TO 60
         DO 50 K=1,N
            TEMP = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = TEMP 
   50    CONTINUE 
         WA(IMAX) = WA(J) 
   60    WA(J) = IMAX 
         JP1 = J+1
         IF (JP1 .GT. N) GO TO 70 
         TEMP = A(J,J)
         DO 65 I = JP1,N
            A(I,J) = A(I,J)/TEMP
   65    CONTINUE 
   70 CONTINUE
   75 IF (IJOB .EQ. 1) GO TO 9005 
      DO 103 K = 1,M
         IW = 0 
         DO 90 I = 1,N
            IMAX = WA(I)
            SUM = B(IMAX,K) 
            B(IMAX,K) = B(I,K)
            IF (IW .EQ. 0) GO TO 85 
            IM1 = I-1 
            DO 80 J = IW,IM1
               SUM = SUM-A(I,J)*B(J,K)
   80       CONTINUE
            GO TO 88
   85       IF (T(1) .NE. ZERO .OR. T(2) .NE. ZERO) IW = I
   88       B(I,K) = SUM
   90    CONTINUE 
         N1 = N+1 
         DO 100 IW = 1,N
            I = N1-IW 
            JP1 = I+1 
            SUM = B(I,K)
            IF (JP1 .GT. N) GO TO 98
            DO 95 J = JP1,N 
               SUM = SUM-A(I,J)*B(J,K)
   95       CONTINUE
   98       B(I,K) = SUM/A(I,I) 
  100    CONTINUE 
  103 CONTINUE
      GO TO 9005
  105 IER = 129 
 9000 CONTINUE
      CALL UERTST(IER,6HLEQT1C) 
 9005 RETURN
      END 
C=====================================================================
      SUBROUTINE UERTST (IER,NAME)
C
C-uertst----------------library 3--------------------------------------
C
C   function            - error message generation
C   usage               - call uertst(ier,name)
C   parameters   ier    - error parameter. type + n  where
C                           type= 128 implies terminal error
C                                  64 implies warning with fix
C                                  32 implies warning
C                           n   = error code relevant to calling routin
C                name   - input scalar containing the name of the
C                           calling routine as a 6-character literal
C                           string.
C   language            - fortran
C----------------------------------------------------------------------
C   latest revision     - august 1, 1973
C
      CHARACTER*10 ITYPE (2,4)
      INTEGER IBIT(4)
      INTEGER            WARN,WARF,TERM,PRINTR
      EQUIVALENCE        (IBIT(1),WARN),(IBIT(2),WARF),(IBIT(3),TERM)
      DATA     ITYPE    /'WARNING   ','          ',
     *                   'WARNING(WI','TH FIX)   ',
     *                   'TERMINAL  ','          ',
     *                   'NON DEFINE','D         '/
                   DATA IBIT/32,64,128,0/
      IER2=IER
      IF (IER2 .GE. WARN) GO TO 5
C                                  non-defined
      IER1=4
      GO TO 20
   5  IF (IER2 .LT. TERM) GO TO 10
C                                  terminal
      IER1=3
      GO TO 20
  10  IF (IER2 .LT. WARF) GO TO 15
C                                  warning(with fix)
      IER1=2
      GO TO 20
C                                  warning
  15  IER1=1
C                                  extract *n*
  20  IER2=IER2-IBIT(IER1)
C                                  print error message
      WRITE(16,25) (ITYPE(I,IER1),I=1,2),NAME,IER2,IER
   25 FORMAT(26H *** I M S L(UERTST) ***  ,2A10,4X,A6,4X,I2,
     1   8H (IER = ,I3,1H))
      RETURN
      END
C=====================================================================
