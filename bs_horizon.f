      program bs_horizon         
c///////////////////////////////////////////////////////////////////////
c    Appendix
c    Nonogaki S, Masumoto S and Shiono K (2008) Optimal Determination
c    of Geologic Boundary Surface Using Cubic B-Spline. Geoinformatics,
c    vol.19, no.2, pp.61-77.
c
c    Fortran77 program for determination of an optimal geologic surface
c    using both elevation data and strike-dip data.
c
c    Note: Run a command "ulimit -s unlimited" if you get stack size
c          error in running this program.
c
c                                         Producted by Susumu Nonogaki
c                                         Last update February 17, 2009
c    1st update(version2): revise the method of knot calculation
c    2nd update(version3): change output format of optimal surface
c///////////////////////////////////////////////////////////////////////
c    Basic concepts
c    1) A surface f(x,y) is expressed by a tensor product form of Cubic
c       B-spline functions:
c         f(x,y) = double sigma (c(ij)Ni(x)Nj(y))
c                               (i = 1 to mx+3; j = 1 to my+3)
c       where Ni(x) and Nj(y) are cubic B-splines, and c(ij) are
c       coefficients multiplied by each tensor product which may be
c       determined in calculation.
c    2) Location of outcrop (xh(k), yh(k), zh(k)) provides equality and
c       inequality constraints for the elevation of surface;
c         f(xh(k),yh(k)) = zh(k)  if lm(k) = 0
c         f(xh(k),yh(k)) < zh(k)  if lm(k) < 0
c         f(xh(k),yh(k)) > zh(k)  if lm(k) > 0
c       where lm(k) is a indicator of data types.
c    3) Strike-dip information at location (xd(k), yd(k), zd(k)) also
c       provides constraints for the partial derivative of surface;
c         fx(xd(k),yd(k)) + sin(tr(k)) * tan(dp(k)) = 0
c         fy(xd(k),yd(k)) + cos(tr(k)) * tan(dp(k)) = 0
c       where tr(k) is trend of the maximum slope measured clockwise
c       from north and dp(k) is dip angle from horizontal surface.
c    4) The surface determination problem is considered as a following
c       constrained optimization problem:
c         Find an optimal solution f(x,y) that minimizes objective
c         function J(f) which evaluates the smoothness of the surface;
c           J(f) = m1 * J1(f) + m2 * J2(f)
c                = m1 * [double integral {(fx)**2 + (fy)**2} dxdy]
c                + m2 * [double integral {(fxx)**2+2*(fxy)**2+(fyy)**2}dxdy]
c         subject to
c           f(xh(k),yh(k)) = zh(k) ; if lm(k) = 0
c           f(xh(k),yh(k)) < zh(k) ; if lm(k) < 0
c           f(xh(k),yh(k)) > zh(k) ; if lm(k) > 0
c         and
c           f(xd(k),yd(k)) + sin(tr(k))*tan(dp(k)) = 0
c           f(xd(k),yd(k)) + cos(tr(k))*tan(dp(k)) = 0.
c    5) To solve the problem, an augmented objective function:
c         Q(f;alpha) = J(f) + alpha * R(f)
c       is introduced based on the exterior penalty function method. Where
c       R(f), defined in a form of residual mean of squares, evaluates
c       the goodness of fit and alpha is a constant called a penalty.
c       The optimal surface is a solution of linear simultaneous equations:
c         A c = b
c       derived from dQ/dc(ij) = 0. The equations are solved by Choleski's method.
c    6) The optimal surface can be saved as the form of both function and
c       grid data. In the case of function output, the coefficients c(ij),
c       which give the optimal surface, are saved. In grid data output,
c       a set of values f(ig,jg) is estimated by assigning the normalized
c       coordinates of grid points (xg(ig),yg(jg)) (ig = 1 to nxg; jg = 1 to nyg)
c       to the equation of surface. 
c///////////////////////////////////////////////////////////////////////
c    Contents of output function
c      numbera of divisions : mx,my
c      origin of knots      : (xmin,ymin)
c      knot intervals       : hx,hy
c      coefficients         : c(ij)        i = 1 to mx+3; j = 1 to my+3
c    Contents of output grid data 
c      numbers of grid      : nxg,nyg
c      origin of grid       : (xming,yming)
c      grid intervals       : dx,dy
c      elevation            : f(ig,jg)     ig = 1 to nxg, jg = 1 to nyg
c///////////////////////////////////////////////////////////////////////
c    Note the limitation of array dimension
c      mxhdt            = 100000  :maximum number of elevation data
c      mxddt            = 10000   :maximum number of strike-dip data
c      mxgrd            = 1001    :maximum number of grid (nxg & nyg)
c      mxm              = 200     :maximum number of division (mx & my)
c                                  (set bigger value between mx and my)
c      mxr=mxm**2       = 40000   :maximum number of rows    of matrix A
c      mxc=(my+3)*3+4   = 613     :maximum number of columns of matrix A
c///////////////////////////////////////////////////////////////////////
      parameter (mxhdt=100000,mxddt=10000,mxgrd=1001,mxm=200,
     &           mxr=40000,mxc=613)
      implicit real*8 (a-h,o-z)
      dimension xh(mxhdt),yh(mxhdt),zh(mxhdt),lm(mxhdt)
      dimension xd(mxddt),yd(mxddt),tr(mxddt),dp(mxddt)
      dimension dfdx(mxddt),dfdy(mxddt)
      dimension c(mxr)
c
c----------Input equality-inequality elevation data---------------------
  100 call xyzdt(mxhdt,xh,yh,zh,lm,nh,xmind,xmaxd,ymind,ymaxd,
     &           zmind,zmaxd)
c
c----------Input trend & dip data---------------------------------------
      call dipdt(mxddt,xd,yd,tr,dp,dfdx,dfdy,nd,xmind,xmaxd,
     &           ymind,ymaxd)
c
c----------Set a region for surface estimation--------------------------
  200 call setreg(mxm,nh,nd,xmind,xmaxd,ymind,ymaxd,zmind,zmaxd,
     &            xmin,xmax,ymin,ymax,mx,my)
c
c----------Solve the simultaneous equations-----------------------------
  300 call solve(mxhdt,mxddt,mxr,mxc,xh,yh,zh,lm,nh,xd,yd,dfdx,dfdy,nd,
     &           xmin,xmax,ymin,ymax,mx,my,hx,hy,neq,c)
c
c----------Output the result -------------------------------------------
  400 call cfout(mxgrd,mxr,xmin,xmax,ymin,ymax,mx,my,hx,hy,neq,c)
c
c----------Set next work------------------------------------------------
      write(*,*)
      write(*,10)
   10 format("#####Next work#####"/
     &       "Change input data            = 1     "/
     &       "Change calculation region    = 2     "/
     &       "Change calculation parameter = 3     "/
     &       "Output result to files       = 4     "/
     &       "End                          = others"/
     &       "                   next work = ",$)
      read(*,*) ians
      goto (100,200,300,400) ians
c
c-----------------------------------------------------------------------
      stop
      end
c
c//////////////////subroutine xyzdt/////////////////////////////////////
c    Input equality-inequality elevation data.
c       Data format :  id, xh(i), yh(i), zh(i), lm(i)   (i = 1 to nh)
c       End of file :   0,   9e9,   9e9,   9e9,    9
c          id                : id number (do not have a relation to calculation)
c          xh(i),yh(i),zh(i) : location of outcrop
c          lm(i)             : relation to surface
c                               0 <=> surface must pass the point
c                              -1 <=> surface must pass under the point
c                               1 <=> surface must pass above the point
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       mxhdt                : maximum number of elevation data
c       xh(nh),yh(nh),zh(nh) : location of outcrop
c       lm(nh)               : relation to surface
c       nh                   : number of elevation data
c       *mind,*maxd          : minimum & maximun value of *
c///////////////////////////////////////////////////////////////////////
      subroutine xyzdt(mxhdt,xh,yh,zh,lm,nh,xmind,xmaxd,ymind,ymaxd,
     &                 zmind,zmaxd)
      implicit real*8 (a-h,o-z)
      dimension xh(mxhdt),yh(mxhdt),zh(mxhdt),lm(mxhdt)
      character*80 flh
c
      write(*,10)
   10 format("#####Data input#####"/
     &       "File name for elevation data   = ",$)
      read(*,11) flh
   11 format(a80)
c
      lu=50
      open(lu,file=flh,err=999)
      read(lu,*) id,xh(1),yh(1),zh(1),lm(1)
      xmind=xh(1)
      xmaxd=xh(1)
      ymind=yh(1)
      ymaxd=yh(1)
      zmind=zh(1)
      zmaxd=zh(1)
      do 100 i=2,mxhdt
        read(lu,*) id,xh(i),yh(i),zh(i),lm(i)
        if(dabs(xh(i)).gt.1.0D9) goto 110 
        if(xh(i).lt.xmind) xmind=xh(i)
        if(xh(i).gt.xmaxd) xmaxd=xh(i)
        if(yh(i).lt.ymind) ymind=yh(i)
        if(yh(i).gt.ymaxd) ymaxd=yh(i)
        if(zh(i).lt.zmind) zmind=zh(i)
        if(zh(i).gt.zmaxd) zmaxd=zh(i)
  100 continue
  110 nh=i-1
      close(lu)
      return
c
  999 write(*,99) flh
   99 format("Error at open file :",a80)
      stop
c
      end
c
c//////////////////subroutine dipdt/////////////////////////////////////
c    Input trend & dip data.
c       Data format : id, xd(i), yd(i), zd(i), tr(i), dp(i)   (i = 1 to nd)
c       End of file :  0,  9e9,  9e9,  9e9,  9e9,  9e9
c          id                : id number (do not have a relation to calculation)
c          xd(i),yd(i),zd(i) : location of outcrop
c          tr(i),dp(i)       : trend & dip in degree unit
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       mxddt                :  maximum number of strike-dip data
c       xd(nd),yd(nd),zd(nd) :  location of outcrop
c       tr(nd),dp(nd)        :  trend & dip in degree unit
c       dfdx(nd),dfdy(nd)    :  partial derivative with respect to x & y
c       nd                   :  number of trend & dip data
c///////////////////////////////////////////////////////////////////////
      subroutine dipdt(mxddt,xd,yd,tr,dp,dfdx,dfdy,nd,xmind,xmaxd,
     &                 ymind,ymaxd)
      implicit real*8 (a-h,o-z)
      dimension xd(mxddt),yd(mxddt),zd(mxddt),tr(mxddt),dp(mxddt)
      dimension dfdx(mxddt),dfdy(mxddt)
      character*80 fld
c
      write(*,10)
   10 format("Do you have trend & dip data ? "/
     &       "                : 1(yes),0(no) = ",$)
      read(*,*) ians
      if(ians.eq.0) then
        nd=0
        return
      endif
c
      write(*,11)
   11 format("File name for trend & dip data = ",$)
      read(*,12) fld
   12 format(a80)
c
      pi=3.14159265358979/1.80D2
      lu=50
      open(lu,file=fld,err=999)
      do 100 i=1,mxddt
        read(lu,*) id,xd(i),yd(i),zd(i),tr(i),dp(i)
        if(dabs(xd(i)).gt.1.0D9) goto 110
        if(xd(i).lt.xmind) xmind=xd(i)
        if(xd(i).gt.xmaxd) xmaxd=xd(i)
        if(yd(i).lt.ymind) ymind=yd(i)
        if(yd(i).gt.ymaxd) ymaxd=yd(i)
        dfdx(i)=-dsin(pi*tr(i))*dtan(pi*dp(i))
        dfdy(i)=-dcos(pi*tr(i))*dtan(pi*dp(i))
  100 continue
  110 nd=i-1
      close(lu)
      if(nd.eq.0) goto 998
      return
c
  999 write(*,99) fld
   99 format("Error at open file :",a80)
      stop
c
  998 write(*,98)
   98 format("Error !! No trend & dip data in file.")
      stop
c
      end
c
c//////////////////subroutine setreg////////////////////////////////////
c    Set a region and the number of division for surface estimation.
c
c          ----------(150,100)
c         |    |    |    |      In the case of
c         |----|----|----|            xmin,xmax=0,150
c         |    |    |    |            ymin,ymax=0,100
c       (0,0)------------             mx,my    =3,2
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       *min,*max       :  edge of calculation region
c       mx,my           :  number of division of the region
c///////////////////////////////////////////////////////////////////////
      subroutine setreg(mxm,nh,nd,xmind,xmaxd,ymind,ymaxd,zmind,zmaxd,
     &                 xmin,xmax,ymin,ymax,mx,my)
      implicit real*8 (a-h,o-z)
c
      write(*,*)
      write(*,10)
   10 format("#####Data information#####")
c
      write(*,11) nh,nd,xmind,xmaxd,ymind,ymaxd,zmind,zmaxd
   11 format("Number      : NH(Height data)      =",i8/
     &       "            : ND(Trend & dip data) =",i8/
     &       "Data area   : x(min),x(max)        =",2f12.4/
     &       "            : y(min),y(max)        =",2f12.4/
     &       "            : z(min),z(max)        =",2f12.4)
c
      write(*,*)
      write(*,12)
   12 format("#####Calculation region#####")
c
      write(*,13)
   13 format("Calculation region : x(min),x(max) = ",$)
      read(*,*) xmin,xmax
      if(xmax.le.xmin) goto 999
c
      write(*,14)
   14 format("                   : y(min),y(max) = ",$)
      read(*,*) ymin,ymax
      if(ymax.le.ymin) goto 999
c
      write(*,15) mxm
   15 format("Number of divisio    Mx,My must be =<",i3/
     &       "                   : Mx,My         = ",$)
      read(*,*) mx,my
      if((mx.le.0).or.(mx.gt.mxm)) goto 999
      if((my.le.0).or.(my.gt.mxm)) goto 999
      return
c
  999 write(*,99)
   99 format("Error: invalid parameter input")
      stop
c
      end
c
c//////////////////subroutine solve/////////////////////////////////////
c    Solve the simultaneous equations based on the exterior penalty method.
c       Q(f;alpha) = J(smoothness) + alpha * R(residual) ---> min
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       hx,hy      : knot interval
c       c(neq)     : solution (coefficients of approximation)
c    Variables
c       alpmin,alpmax : initial and final penalty in iterative calculation
c       nitr          : number of iteration
c       gam           : weight of trend & dip data to elevation data
c       pm1,pm2       : weight of 1st & 2nd order partial derivative term in J
c       neq           : number of rows    of matrix A (number of equations)
c       nbs           : number of columns of matrix A (for choleski's method)
c       A(neq,nbs)    : coefficient matrix of simultaneous equations
c       b(neq)        : constant matrix of simultaneous equations
c///////////////////////////////////////////////////////////////////////
      subroutine solve(mxhdt,mxddt,mxr,mxc,xh,yh,zh,lm,nh,xd,yd,dfdx,
     &                 dfdy,nd,xmin,xmax,ymin,ymax,mx,my,hx,hy,neq,c)
      implicit real*8 (a-h,o-z)
      dimension xh(mxhdt),yh(mxhdt),zh(mxhdt),lm(mxhdt)
      dimension xd(mxddt),yd(mxddt),dfdx(mxddt),dfdy(mxddt)
      dimension A(mxr,mxc),b(mxr),c(mxr)
c
      write(*,*)
      write(*,10)
   10 format("#####Calculation parameter#####")
c
c----------Penalty : alpha for elevation data---------------------------
      itry=0
      do 100 i=1,nh
        itry=itry+lm(i)**2
  100 continue
      if(itry.gt.0) goto 110
c
      write(*,11)
   11 format("No iteration is required."/
     &       "Q= J+alpha*R       : alpha         = ",$)
      read(*,*) alpmin
      if(alpmin.le.0.0D0) goto 999
      ratio=1.0D0
      nitr=1
      goto 120
c
  110 write(*,12)
   12 format("Some iteration is required."/
     &       "Q= J+alpha*R       : alpha(min)    = ",$)
      read(*,*) alpmin
      if(alpmin.le.0.0D0) goto 999
      write(*,13)
   13 format("                   : alpha(max)    = ",$)
      read(*,*) alpmax
      if(alpmax.lt.alpmin) goto 999
      write(*,14)
   14 format("            number of iteration    = ",$)
      read(*,*) nitr
      if(nitr.le.1) goto 999
      rn=1.0D0/dfloat(nitr-1)
      ratio=(alpmax/alpmin)**rn
c
c----------Penalty : gamma for trend & dip data-------------------------
  120 if(nd.eq.0) then
        gam=0.0D0
      else  
        write(*,15)
   15   format("R=Rh+gamma*Rd      : gamma         = ",$)
        read(*,*) gam
        if(gam.lt.0.0D0) goto 999
      endif
c
c----------Parameter : pm1,pm2------------------------------------------
      write(*,16)
   16 format("J=m1*[(fx)**2+(fy)**2]+"/
     &       "  m2*[(fxx)**2+2(fxy)**2+(fyy)**2]"/
     &       "                   : m1,m2         = ",$)
      read(*,*) pm1,pm2
      if((pm1.lt.0.0D0).or.(pm2.lt.0.0D0)) goto 999
c
c----------Set knot interval--------------------------------------------
      hx=(xmax-xmin)/dfloat(mx)
      hy=(ymax-ymin)/dfloat(my)
c
c----------Calculate the number of rows & columns of matrix A-----------
      neq=(mx+3)*(my+3)
      nbs=(mx+3)*3+4
c
c----------Initialize matrix c------------------------------------------
      do 200 i=1,neq
        c(i)=0.0D0
  200 continue
c
c----------Iterative calculation----------------------------------------
      do 300 itr=1,nitr
        do 310 i=1,neq
          b(i)=0.0D0
          do 320 j=1,nbs
            A(i,j)=0.0D0
  320     continue
  310   continue
        alp=alpmin*(ratio**dfloat(itr-1))
c----------Set matrix b and Ah------------------------------------------
        call setbAh(mxhdt,mxr,mxc,xh,yh,zh,lm,nh,xmin,xmax,ymin,ymax,
     &              mx,my,alp,c,b,A,hx,hy,nhdt,neq,nbs)
c----------Set matrix b and Ad------------------------------------------
        call setbAd(mxddt,mxr,mxc,xd,yd,dfdx,dfdy,nd,xmin,xmax,
     &              ymin,ymax,mx,my,alp,gam,b,A,hx,hy,nddt)
        write(*,*)
        write (*,20) nhdt,nddt
   20   format("#####Number of constraints#####"/
     &         "elevation constraints : NH         =",i7/
     &         "    slope constraints : ND         =",i7)
c----------Set matrix Ax,Ay,Axx,Axy and Ayy-----------------------------
        call setJ(mxr,mxc,mx,my,pm1,pm2,A,hx,hy)
c----------Solve simultaneous equations---------------------------------
        call choles(mxr,mxc,neq,A,b,c,nbs)
c----------Evaluate Rh--------------------------------------------------
        call Rhvalue(mxhdt,mxr,xh,yh,zh,lm,nh,xmin,xmax,ymin,ymax,
     &               mx,my,hx,hy,nhdt,c,QRh,QRhav)
c----------Evaluate Rd--------------------------------------------------
        call Rdvalue(mxddt,mxr,xd,yd,dfdx,dfdy,nd,xmin,xmax,
     &               ymin,ymax,mx,my,hx,hy,nddt,c,QRd,QRdav)
c----------Evaluate J---------------------------------------------------
        call Jvalue(mxr,mx,my,pm1,pm2,hx,hy,c,QJx,QJy,
     &              QJxx,QJxy,QJyy,QJ1,QJ2,QJ)
c----------Evaluate Q---------------------------------------------------
        QR=QRh+gam*QRd
        aR=alp*QR
        Q=QJ+aR
c----------Output the result--------------------------------------------
        write(*,21)
   21   format("#####Estimation result#####")
        write(*,22) itr,alp,gam
   22   format("iteration=",i4,"  alpha=",E12.5,"  gamma=",E12.5/
     &         "------------------------------------------------------")
        write(*,23) nhdt,QRh,QRhav
   23   format("H : data =",i6,"  Rh =",E12.5,"  Error=",E12.5)
        write(*,24) nddt,QRd,QRdav
   24   format("D : data =",i6,"  Rd =",E12.5,"  Error=",E12.5/
     &         "------------------------------------------------------")
        write(*,25) QJx,QJy
   25   format("Jx =",E12.5,"  Jy =",E12.5)
        write(*,26) QJxx,QJxy,QJyy
   26   format("Jxx=",E12.5,"  Jxy=",E12.5,"  Jyy=",E12.5)
        write(*,27) QJ1,QJ2
   27   format("J1 =",E12.5,"  J2 =",E12.5/
     &         "------------------------------------------------------")
        write(*,28) Q,QJ,aR
   28   format("Q  =",E12.5,"  J  =",E12.5,"  aR =",E12.5)
  300 continue
c
      return
c
  999 write(*,99)
   99 format("Error: invalid parameter input")
      stop
c
      end
c
c//////////////////subroutine bspl//////////////////////////////////////
c    Calculate 4 non-zero B-spline at data point.
c       4 basic functions:
c         bn1(r) =       r**3                         / 6
c         bn2(r) = (-3 * r**3 + 3 * r**2 + 3 * r + 1) / 6
c         bn3(r) = ( 3 * r**3 - 6 * r**2         + 4) / 6
c         bn4(r) = (-    r**3 + 3 * r**2 - 3 * r + 1) / 6
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       r             :  normalized coordinate (0.0 to 1.0)
c       h             :  knot interval
c       ider          :  order of derivative
c       bn(4)         :  4 non-zero B-spline
c///////////////////////////////////////////////////////////////////////
      subroutine bspl(r,h,ider,bn)
      implicit real*8 (a-h,o-z)
      dimension bn(4)
c
      if(ider.eq.0) then
        a=1.0D0/6.0D0
        bn(1)=r**3*a
        bn(2)=(3.0D0*(-r**3+r**2+r)+1.0D0)*a
        bn(3)=(3.0D0*r**3-6.0D0*r**2+4.0D0)*a
        bn(4)=(-r**3+3.0D0*(r**2-r)+1.0D0)*a
      else if(ider.eq.1) then
        a=0.5D0/h
        bn(1)=r**2*a
        bn(2)=(-3.0D0*r**2+2.0D0*r+1.0D0)*a
        bn(3)=(3.0D0*r**2-4.0D0*r)*a
        bn(4)=(-r**2+2.0D0*r-1.0D0)*a
      else if(ider.eq.2) then
        a=1.0D0/(h**2)
        bn(1)=r*a
        bn(2)=(-3.0D0*r+1.0D0)*a
        bn(3)=(3.0D0*r-2.0D0)*a
        bn(4)=(-r+1.0D0)*a
      endif
c
      return
      end
c
c//////////////////subroutine setbAh////////////////////////////////////
c    Set the elements of matrix b & A (Ah) derived from elevation data.
c       Ah(ij,i'j') = alpha / nh * sigma (Nx(i)Nx(i')Ny(j)Ny(j'))
c       b(ij)       = alpha / nh * sigma (Nx(i)Ny(j)zh)
c       (i,i' = 1 to mx + 3; j,j' = 1 to my + 3)
c    When knot number on right-hand neighbor is 'i', only 4 B-spline
c    have non-zero values, that is, N(i), N(i+1), N(i+2) and N(i+3).
c    Using 4 non-zero B-spline, N(i+j) can be expressed as follows.
c       N(ip+j) = b(4-j)    (j = 0 to 3)
c-----------------------------------------------------------------------
c    Variables
c       ip* : knot number on right-hand neighbor
c       r*  : normalized coordinate from knot on left-hand neigbor
c///////////////////////////////////////////////////////////////////////
      subroutine setbAh(mxhdt,mxr,mxc,xh,yh,zh,lm,nh,xmin,xmax,ymin,ymax
     &                 ,mx,my,alp,c,b,A,hx,hy,nhdt,neq,nbs)
      implicit real*8 (a-h,o-z)
      dimension xh(mxhdt),yh(mxhdt),zh(mxhdt),lm(mxhdt)
      dimension A(mxr,mxc),b(mxr),c(mxr)
      dimension bx0(4),by0(4)
c
      nhdt=0
      do 100 id=1,nh
        if(xh(id).lt.xmin.or.xh(id).gt.xmax) goto 100
        if(yh(id).lt.ymin.or.yh(id).gt.ymax) goto 100
        nhdt=nhdt+1
  100 continue
c
      wt=alp/dfloat(nhdt)
      do 200 id=1,nh
        if(xh(id).lt.xmin.or.xh(id).gt.xmax) goto 200
        if(yh(id).lt.ymin.or.yh(id).gt.ymax) goto 200
c---Determine knot number on right-hand neighbor------------------------
        ipx=dint((xh(id)-xmin)/hx)+1
        ipy=dint((yh(id)-ymin)/hy)+1
c---Consider the data on final knot as the one in the final section-----
        if(ipx.eq.mx+1) ipx=mx
        if(ipy.eq.my+1) ipy=my
c---Calcurate normalized coordinate from knot on left-hand neigbor------
        rx=(xh(id)-xmin)/hx-dfloat(ipx-1)
        ry=(yh(id)-ymin)/hy-dfloat(ipy-1)
c---Calculate 4 non-zero B-spline with respect to x & y-----------------
        call bspl(rx,hx,0,bx0)
        call bspl(ry,hy,0,by0)
c---Examine whether data is satisfied with constraints or not-----------
        fh=0.0D0
        do 210 j=0,3
          do 220 i=0,3
            ij=ipx+i+(ipy-1+j)*(mx+3)
            fh=fh+c(ij)*bx0(4-i)*by0(4-j)
  220     continue
  210   continue
        hres=fh-zh(id)
        if((lm(id).ne.0).and.(hres*dfloat(lm(id)).gt.0.0D0)) goto 200
c---Set b & Ah using 4 non-zero B-spline--------------------------------
        do 230 j1=0,3
          do 240 i1=0,3
            ia=ipx+i1+(ipy-1+j1)*(mx+3)
            b(ia)=b(ia)+wt*bx0(4-i1)*by0(4-j1)*zh(id)
            do 250 j2=0,3
              do 260 i2=0,3
                ja=i2-i1+1+(j2-j1)*(mx+3)
                if(ja.le.0) goto 260
                A(ia,ja)=A(ia,ja)
     &                   +wt*bx0(4-i1)*by0(4-j1)*bx0(4-i2)*by0(4-j2)
  260         continue
  250       continue
  240      continue
  230   continue
  200 continue
c
      return
      end
c
c//////////////////subroutine setbAd////////////////////////////////////
c     Set the elements of matrix b and A (Ad) derived from trend & dip data.
c        Ad(ij,i'j') =   alpha * gamma / nd 
c          * sigma (Nx'(i)Nx'(i')Ny(j)Ny(j')+Nx(i)Nx(i')Ny'(j)Ny'(j'))
c        b(ij)       = - alpha * gamma / nd
c          * sigma (Nx'(i)Ny(j)sin(tr)tan(dp)+Nx(i)Ny'(j)cos(tr)tan(dp))
c        (i,i' = 1 to mx + 3; j,j' = 1 to my + 3)
c///////////////////////////////////////////////////////////////////////
      subroutine setbAd(mxddt,mxr,mxc,xd,yd,dfdx,dfdy,nd,xmin,xmax,
     &                  ymin,ymax,mx,my,alp,gam,b,A,hx,hy,nddt)
      implicit real*8 (a-h,o-z)
      dimension xd(mxddt),yd(mxddt),dfdx(mxddt),dfdy(mxddt)
      dimension A(mxr,mxc),b(mxr)
      dimension bx0(4),by0(4),bx1(4),by1(4)
c
      nddt=0
      if(nd.eq.0) return
      do 100 id=1,nd
        if(xd(id).lt.xmin.or.xd(id).gt.xmax) goto 100
        if(yd(id).lt.ymin.or.yd(id).gt.ymax) goto 100
        nddt=nddt+1
  100 continue
c
      beta=alp*gam
      wt=beta/dfloat(nddt)
      do 200 id=1,nd
        if(xd(id).lt.xmin.or.xd(id).gt.xmax) goto 200
        if(yd(id).lt.ymin.or.yd(id).gt.ymax) goto 200
c---Determine knot number on right-hand neighbor------------------------
        ipx=dint((xd(id)-xmin)/hx)+1
        ipy=dint((yd(id)-ymin)/hy)+1
c---Consider the data on final knot as the one in the final section-----
        if(ipx.eq.mx+1) ipx=mx
        if(ipy.eq.my+1) ipy=my
c---Calcurate normalized coordinate from knot on left-hand neigbor------
        rx=(xd(id)-xmin)/hx-dfloat(ipx-1)
        ry=(yd(id)-ymin)/hy-dfloat(ipy-1)
c---Calculate 4 non-zero B-spline with respect to x & y-----------------
        call bspl(rx,hx,0,bx0)
        call bspl(rx,hx,1,bx1)
        call bspl(ry,hy,0,by0)
        call bspl(ry,hy,1,by1)
c---Set b & Ad using 4 non-zero B-spline--------------------------------
        do 210 j1=0,3
          do 220 i1=0,3
            ia=ipx+i1+(ipy-1+j1)*(mx+3)
            b(ia)=b(ia)
     &            +wt*bx1(4-i1)*by0(4-j1)*dfdx(id)
     &            +wt*bx0(4-i1)*by1(4-j1)*dfdy(id)
            do 230 j2=0,3
              do 240 i2=0,3
                ja=i2-i1+1+(j2-j1)*(mx+3)
                if(ja.le.0) goto 240
                A(ia,ja)=A(ia,ja)
     &                   +wt*bx1(4-i1)*by0(4-j1)*bx1(4-i2)*by0(4-j2)
     &                   +wt*bx0(4-i1)*by1(4-j1)*bx0(4-i2)*by1(4-j2)
  240         continue
  230       continue
  220      continue
  210   continue
  200 continue
c
      return
      end
c
c//////////////////subroutine setJ//////////////////////////////////////
c     Set the elements of matrix A (Ax, Ay, Axx, Ayy and Axy) derived
c     from the integration of 1st & 2nd order partial derivative squared
c     with respect to x & y.
c        Ax(ij,i'j')  =     m1 * [integral Nx'(i)Nx'(i') dx]
c                                        * [integral Ny(j)Ny(j') dy]
c        Ay(ij,i'j')  =     m1 * [integral Nx(i)Nx(i') dx]
c                                        * [integral Ny'(j)Ny'(j') dy]
c        Axx(ij,i'j') =     m2 * [integral Nx''(i)Nx''(i') dx]
c                                        * [integral Ny(j)Ny(j') dy]
c        Axy(ij,i'j') = 2 * m2 * [integral Nx'(i)Nx'(i') dx]
c                                        * [integral Ny'(j)Ny'(j') dy]
c        Ayy(ij,i'j') =     m2 * [integral Nx(i)Nx(i') dx]
c                                        * [integral Ny''(j)Ny''(j') dy]
c        (i,i' = 1 to mx + 3; j,j' = 1 to my + 3)
c///////////////////////////////////////////////////////////////////////
      subroutine setJ(mxr,mxc,mx,my,pm1,pm2,A,hx,hy)
      implicit real*8 (a-h,o-z)
      dimension ffix0(4),ffix1(4),ffix2(4)
      dimension ffiy0(4),ffiy1(4),ffiy2(4)
      dimension A(mxr,mxc)
c
      wx=pm1*(hy/hx)/dfloat(mx*my)
      wy=pm1*(hx/hy)/dfloat(mx*my)
      wxx=pm2*hy/(hx**3)
      wxy=2.0D0*pm2/(hx*hy)
      wyy=pm2*hx/(hy**3)
      do 100 j1=1,my+3
c---Calculate the integration with respect to y-------------------------
        call integral(my,0,j1,ffiy0)
        call integral(my,1,j1,ffiy1)
        call integral(my,2,j1,ffiy2)
        j22=4
        if(j1.ge.my+1) j22=my+4-j1
        do 110 i1=1,mx+3
c---Calculate the integration with respect to x-------------------------
          call integral(mx,0,i1,ffix0)
          call integral(mx,1,i1,ffix1)
          call integral(mx,2,i1,ffix2)
          ia=i1+(j1-1)*(mx+3)
          i22=4
          if(i1.ge.mx+1) i22=mx+4-i1
c---Set Ax, Ay, Axx, Ayy and Axy----------------------------------------
          call setJtmp(mxr,mxc,mx,A,wx,ia,i22,j22,ffix1,ffiy0)
          call setJtmp(mxr,mxc,mx,A,wy,ia,i22,j22,ffix0,ffiy1)
          call setJtmp(mxr,mxc,mx,A,wxx,ia,i22,j22,ffix2,ffiy0)
          call setJtmp(mxr,mxc,mx,A,wxy,ia,i22,j22,ffix1,ffiy1)
          call setJtmp(mxr,mxc,mx,A,wyy,ia,i22,j22,ffix0,ffiy2)
  110  continue
  100 continue
c
      return
      end
c
c//////////////////subroutine setJtmp///////////////////////////////////
c     Calculate the individual element values of A by multiply weight
c     with two kinds of integration.
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       w         : weight of element 
c       ia        : row number of objective element
c       i22,j22   : indicator wheather element may be copyed or not
c       ffix,ffiy : integration value
c///////////////////////////////////////////////////////////////////////
      subroutine setJtmp(mxr,mxc,mx,A,w,ia,i22,j22,ffix,ffiy)
      implicit real*8 (a-h,o-z)
      dimension ffix(4),ffiy(4)
      dimension A(mxr,mxc)
c
      do 100 j2=1,j22
        do 110 i2=1,i22
          ja=i2+(j2-1)*(mx+3)
          tmp=w*ffix(i2)*ffiy(j2)
          A(ia,ja)=A(ia,ja)+tmp
          if((j2.ge.2).and.(i2.ge.2)) then
            iia=ia+i2-1
            jja=ja-2*(i2-1)
            A(iia,jja)=A(iia,jja)+tmp
          endif
  110   continue
  100 continue
c
      return
      end
c
c//////////////////subroutine integral//////////////////////////////////
c    Calculate integration of 4 basic non-zero B-spline functions
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       m      : number of division of the region
c       ider   : order of derivative
c       ia     : knot number
c       ffi(4) : integration value of 4 basic functions
c///////////////////////////////////////////////////////////////////////
      subroutine integral(m,ider,ia,ffi)
      implicit real*8 (a-h,o-z)
      dimension ffi(4)
c
c---Calculate the integration of basic functions------------------------
      if(ider.eq.0) then
        if(ia.eq.1) then
          ffi(1)=  1.0D0/ 2.52D2 ; ffi(2)=  4.3D1/1.680D3
          ffi(3)=  1.0D0/  8.4D1 ; ffi(4)=  1.0D0/5.040D3
        else if(ia.eq.2) then
          ffi(1)= 1.51D2/ 6.30D2 ; ffi(2)=  5.9D1/ 2.80D2
          ffi(3)=  1.0D0/  4.2D1 ; ffi(4)=  1.0D0/5.040D3
        else if(ia.eq.3) then
          ffi(1)= 5.99D2/1.260D3 ; ffi(2)= 3.97D2/1.680D3
          ffi(3)=  1.0D0/  4.2D1 ; ffi(4)=  1.0D0/5.040D3
        else if(ia.eq.m+1) then
          ffi(1)= 5.99D2/1.260D3 ; ffi(2)=  5.9D1/ 2.80D2
          ffi(3)=  1.0D0/  8.4D1 ; ffi(4)=  0.0D0
        else if(ia.eq.m+2) then
          ffi(1)= 1.51D2/ 6.30D2 ; ffi(2)=  4.3D1/1.680D3
          ffi(3)=  0.0D0         ; ffi(4)=  0.0D0
        else if(ia.eq.m+3) then
          ffi(1)=  1.0D0/ 2.52D2 ; ffi(2)=  0.0D0
          ffi(3)=  0.0D0         ; ffi(4)=  0.0D0
        else
          ffi(1)= 1.51D2/ 3.15D2 ; ffi(2)= 3.97D2/1.680D3
          ffi(3)=  1.0D0/  4.2D1 ; ffi(4)=  1.0D0/5.040D3
        endif
c---Calculate the integration of 1st order
c          partial derivative of basic functions------------------------
      else if(ider.eq.1) then
        if(ia.eq.1) then
          ffi(1)=  1.0D0/ 2.0D1 ; ffi(2)=  7.0D0/1.20D2
          ffi(3)=- 1.0D0/ 1.0D1 ; ffi(4)=- 1.0D0/1.20D2
        else if(ia.eq.2) then
          ffi(1)=  1.0D0/ 3.0D0 ; ffi(2)=- 1.1D1/ 6.0D1
          ffi(3)=- 1.0D0/ 5.0D0 ; ffi(4)=- 1.0D0/1.20D2
        else if(ia.eq.3) then
          ffi(1)=  3.7D1/ 6.0D1 ; ffi(2)=- 1.0D0/ 8.0D0
          ffi(3)=- 1.0D0/ 5.0D0 ; ffi(4)=- 1.0D0/1.20D2
        else if(ia.eq.m+1) then
          ffi(1)=  3.7D1/ 6.0D1 ; ffi(2)=- 1.1D1/ 6.0D1
          ffi(3)=- 1.0D0/ 1.0D1 ; ffi(4)=  0.0D0
        else if(ia.eq.m+2) then
          ffi(1)=  1.0D0/ 3.0D0 ; ffi(2)=  7.0D0/1.20D2
          ffi(3)=  0.0D0        ; ffi(4)=  0.0D0
        else if(ia.eq.m+3) then
          ffi(1)=  1.0D0/ 2.0D1 ; ffi(2)=  0.0D0
          ffi(3)=  0.0D0        ; ffi(4)=  0.0D0
        else
          ffi(1)=  2.0D0/ 3.0D0 ; ffi(2)=- 1.0D0/ 8.0D0
          ffi(3)=- 1.0D0/ 5.0D0 ; ffi(4)=- 1.0D0/1.20D2
        endif
c---Calculate the integration of 2nd order
c          partial derivative of basic functions------------------------
      else if(ider.eq.2) then
        if(ia.eq.1) then
          ffi(1)= 1.0D0/3.0D0 ; ffi(2)=-1.0D0/2.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 1.0D0/6.0D0
        else if(ia.eq.2) then
          ffi(1)= 4.0D0/3.0D0 ; ffi(2)=-1.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 1.0D0/6.0D0
        else if(ia.eq.3) then
          ffi(1)= 7.0D0/3.0D0 ; ffi(2)=-3.0D0/2.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 1.0D0/6.0D0
        else if(ia.eq.m+1) then
          ffi(1)= 7.0D0/3.0D0 ; ffi(2)=-1.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 0.0D0
        else if(ia.eq.m+2) then
          ffi(1)= 4.0D0/3.0D0 ; ffi(2)=-1.0D0/2.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 0.0D0
        else if(ia.eq.m+3) then
          ffi(1)= 1.0D0/3.0D0 ; ffi(2)= 0.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 0.0D0
        else
          ffi(1)= 8.0D0/3.0D0 ; ffi(2)=-3.0D0/2.0D0
          ffi(3)= 0.0D0       ; ffi(4)= 1.0D0/6.0D0
        endif
      endif
c
      return
      end
c
c//////////////////subroutine choles////////////////////////////////////
c     Solving a normal equation with choleski's method
c     See Ichida & Yoshimoto (1979) for original program.
c///////////////////////////////////////////////////////////////////////
      subroutine choles(mxr,mxc,neq,A,b,c,nbs)
      implicit real*8 (a-h,o-z)
      dimension A(mxr,mxc),b(mxr),c(mxr)
c
      A(1,1)=dsqrt(A(1,1))
      b(1)=b(1)/A(1,1)
      do 10 jt=2,nbs
        A(1,jt)=A(1,jt)/A(1,1)
   10 continue
      do 15 ia=2,neq
        iam1=ia-1
        if(ia.ge.nbs) goto 20
        ias=1
        goto 25
   20   ias=ia-nbs+1
   25   tki=0.0D0
        do 30 ka=ias,iam1
          iat=ia-ka+1
          tki=tki+A(ka,iat)**2
   30   continue
        A(ia,1)=dsqrt(A(ia,1)-tki)
        if(neq-ia.ge.nbs) goto 31
        jtl=neq-ia+1
        if(jtl.eq.1) goto 38
        goto 32
   31   jtl=nbs
   32   do 35 jt=2,jtl
          if(ia.ge.nbs) goto 36
          if(jt.gt.nbs-ia) goto 34
          iasj=1
          goto 37
   34     iasj=jt-(nbs-ia)
          goto 37
   36     iasj=ias+jt-1
   37     tki=0.0D0
          if(iasj.gt.iam1) goto 41
          do 40 ka=iasj,iam1
            iat=ia-ka+1
            jtt=iat+jt-1
            tki=tki+A(ka,iat)*A(ka,jtt)
   40     continue
   41     A(ia,jt)=(A(ia,jt)-tki)/A(ia,1)
   35   continue
   38   bki=0.0D0
        do 45 ka=ias,iam1
          iat=ia-ka+1
          bki=bki+A(ka,iat)*b(ka)
   45   continue
        b(ia)=(b(ia)-bki)/A(ia,1)
   15 continue
      c(neq)=b(neq)/A(neq,1)
      neqm1=neq-1
      do 50 iat=1,neqm1
        ia=neq-iat
        if(iat.ge.nbs) goto 55
        net=iat+1
        goto 60
   55   net=nbs
   60   tx=0.0D0
        do 65 ka=2,net
          kat=ia+ka-1
          tx=tx+A(ia,ka)*c(kat)
   65   continue
        c(ia)=(b(ia)-tx)/A(ia,1)
   50 continue
c
      return
      end
c
c//////////////////subroutine Rhvalue///////////////////////////////////
c    Evaluate the residual mean of squares between elevation data and
c    final solution. 
c        Rh(f) = sigma [{double sigma (c(ij)Nx(i)Ny(j)) - zh}**2] / nh
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       QRh   : residual mean of squares
c       QRhav : square root of QRh
c///////////////////////////////////////////////////////////////////////
      subroutine Rhvalue(mxhdt,mxr,xh,yh,zh,lm,nh,xmin,xmax,ymin,ymax,
     &                   mx,my,hx,hy,nhdt,c,QRh,QRhav)
      implicit real*8 (a-h,o-z)
      dimension xh(mxhdt),yh(mxhdt),zh(mxhdt),lm(mxhdt)
      dimension c(mxr)
      dimension bx0(4),by0(4)
c
      nhdt=0
      Rh=0.0D0
      do 100 id=1,nh
        if(xh(id).lt.xmin.or.xh(id).gt.xmax) goto 100
        if(yh(id).lt.ymin.or.yh(id).gt.ymax) goto 100
        ipx=dint((xh(id)-xmin)/hx)+1
        ipy=dint((yh(id)-ymin)/hy)+1
        if(ipx.eq.mx+1) ipx=mx
        if(ipy.eq.my+1) ipy=my
        rx=(xh(id)-xmin)/hx-dfloat(ipx-1)
        ry=(yh(id)-ymin)/hy-dfloat(ipy-1)
        call bspl(rx,hx,0,bx0)
        call bspl(ry,hy,0,by0)
        fh=0.0D0
        do 110 j=0,3
          do 120 i=0,3
            ij=ipx+i+(ipy-1+j)*(mx+3)
            fh=fh+c(ij)*bx0(4-i)*by0(4-j)
  120     continue
  110   continue
        hres=fh-zh(id)
        if((lm(id).ne.0).and.(hres*dfloat(lm(id)).gt.0.0D0)) goto 100
        Rh=Rh+hres**2
        nhdt=nhdt+1
  100 continue
      QRh=Rh/dfloat(nhdt)
      QRhav=dsqrt(QRh)
c
      return
      end
c
c//////////////////subroutine Rdvalue///////////////////////////////////
c    Evaluate the residual mean of squares between trend & dip data and
c    final solution.
c        Rd(f) = sigma [
c                {double sigma (c(ij)Nx'(i)Ny(j))+sin(tr)tan(dp)}**2
c                {double sigma (c(ij)Nx(i)Ny'(j))+cos(tr)tan(dp)}**2
c                ] / nd
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       QRd   : residual mean of squares
c       QRdav : square root of QRd
c///////////////////////////////////////////////////////////////////////
      subroutine Rdvalue(mxddt,mxr,xd,yd,dfdx,dfdy,nd,xmin,xmax,
     &                   ymin,ymax,mx,my,hx,hy,nddt,c,QRd,QRdav)
      implicit real*8 (a-h,o-z)
      dimension xd(mxddt),yd(mxddt),dfdx(mxddt),dfdy(mxddt)
      dimension c(mxr)
      dimension bx0(4),by0(4),bx1(4),by1(4)
c
      if(nd.eq.0) then
       QRd=0.0D0
       QRdav=0.0D0
       return
      endif
c
      Rdx=0.0D0
      Rdy=0.0D0
      do 100 id=1,nd
        if(xd(id).lt.xmin.or.xd(id).gt.xmax) goto 100
        if(yd(id).lt.ymin.or.yd(id).gt.ymax) goto 100
        ipx=dint((xd(id)-xmin)/hx)+1
        ipy=dint((yd(id)-ymin)/hy)+1
        if(ipx.eq.mx+1) ipx=mx
        if(ipy.eq.my+1) ipy=my
        rx=(xd(id)-xmin)/hx-dfloat(ipx-1)
        ry=(yd(id)-ymin)/hy-dfloat(ipy-1)
        call bspl(rx,hx,0,bx0)
        call bspl(rx,hx,1,bx1)
        call bspl(ry,hy,0,by0)
        call bspl(ry,hy,1,by1)
        fdx=0.0D0
        fdy=0.0D0
        do 110 j=0,3
          do 120 i=0,3
            ij=ipx+i+(ipy-1+j)*(mx+3)
            fdx=fdx+c(ij)*bx1(4-i)*by0(4-j)
            fdy=fdy+c(ij)*bx0(4-i)*by1(4-j)
  120     continue
  110   continue
        dxres=fdx-dfdx(id)
        dyres=fdy-dfdy(id)
        Rdx=Rdx+dxres**2
        Rdy=Rdy+dyres**2
  100 continue
      QRd=(Rdx+Rdy)/dfloat(nddt)
      QRdav=dsqrt(QRd)
c
      return
      end
c
c//////////////////subroutine Jvalue////////////////////////////////////
c    Evaluate the smoothness of final solution.
c       J(f) = m1 * J1(f) + m2 * J2(f)
c            = m1 * [double integral {(fx)**2 + (fy)**2} dxdy]
c            + m2 * [double integral {(fxx)**2+2*(fxy)**2+(fyy)**2}dxdy]
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       QJx   : double integration of (fx)**2
c       QJy   : double integration of (fy)**2
c       QJxx  : double integration of (fxx)**2
c       QJxy  : double integration of (fxy)**2
c       QJyy  : double integration of (fyy)**2
c       QJ1   : value of J1(f)
c       QJ2   : value of J2(f)
c       QJ    : value of J(f)
c///////////////////////////////////////////////////////////////////////
      subroutine Jvalue(mxr,mx,my,pm1,pm2,hx,hy,c,QJx,QJy,
     &                  QJxx,QJxy,QJyy,QJ1,QJ2,QJ)
      implicit real*8 (a-h,o-z)
      dimension ffix0(4),ffix1(4),ffix2(4)
      dimension ffiy0(4),ffiy1(4),ffiy2(4)
      dimension c(mxr)
c
      wx=hy/hx/dfloat(mx*my)
      wy=hx/hy/dfloat(mx*my)
      wxx=hy/(hx**3)
      wxy=1.0D0/(hx*hy)
      wyy=hx/(hy**3)
      QJx=0.0D0
      QJy=0.0D0
      QJxx=0.0D0
      QJxy=0.0D0
      QJyy=0.0D0
      do 100 j1=1,my+3
        call integral(my,0,j1,ffiy0)
        call integral(my,1,j1,ffiy1)
        call integral(my,2,j1,ffiy2)
        j22=4
        if(j1.ge.my+1) j22=my+4-j1
        do 110 i1=1,mx+3
          call integral(mx,0,i1,ffix0)
          call integral(mx,1,i1,ffix1)
          call integral(mx,2,i1,ffix2)
          i22=4
          if(i1.ge.mx+1) i22=mx+4-i1
          call Jtmp(mxr,mx,c,wx,i1,j1,i22,j22,ffix1,ffiy0,QJx)
          call Jtmp(mxr,mx,c,wy,i1,j1,i22,j22,ffix0,ffiy1,QJy)
          call Jtmp(mxr,mx,c,wxx,i1,j1,i22,j22,ffix2,ffiy0,QJxx)
          call Jtmp(mxr,mx,c,wxy,i1,j1,i22,j22,ffix1,ffiy1,QJxy)
          call Jtmp(mxr,mx,c,wyy,i1,j1,i22,j22,ffix0,ffiy2,QJyy)
  110    continue
  100 continue
      QJ1=(QJx+QJy)
      QJ2=QJxx+2.0D0*QJxy+QJyy
      QJ=pm1*QJ1+pm2*QJ2
c
      return
      end
c
c//////////////////subroutine Jtmp//////////////////////////////////////
c     Calculate the double integration of partial derivative squared
c     using final coefficients matrix c.
c        Jx  = ct * Ax  * c
c        Jy  = ct * Ay  * c
c        Jxx = ct * Axx * c
c        Jxy = ct * Axy * c * 2
c        Jyy = ct * Ayy * c
c     where ct is an inverse matrix of c.
c-----------------------------------------------------------------------
c    Arguments used in subroutine
c       w       : weight of element 
c       i1,j1   : element position indicator
c       i22,j22 : indicator where is the limit of copy
c       QJ      : double integration partial derivative squared
c///////////////////////////////////////////////////////////////////////
      subroutine Jtmp(mxr,mx,c,w,i1,j1,i22,j22,ffix,ffiy,QJ)
      implicit real*8 (a-h,o-z)
      dimension ffix(4),ffiy(4)
      dimension c(mxr)
c
c      ja=i2+(ia-1)+(j2-1)*(mx+3)
c      i = i1
c      i'= i1+i2-1
c      j = j1
c      j'= j1+j2-1
      do 100 j2=1,j22
        do 110 i2=1,i22
        tmp=w*ffix(i2)*ffiy(j2)
        ia1=i1+(j1-1)*(mx+3)
        ja1=(i1+i2-1)+((j1+j2-1)-1)*(mx+3)
        if(ia1.eq.ja1) then
          QJ=QJ+c(ia1)*tmp*c(ja1)
        else
          QJ=QJ+c(ia1)*tmp*c(ja1)*2.0D0
          if((j2.ne.1).and.(i2.ne.1)) then
            ia2=(i1+i2-1)+(j1-1)*(mx+3)
            ja2=i1+((j1+j2-1)-1)*(mx+3)
            QJ=QJ+c(ia2)*tmp*c(ja2)*2.0D0
          endif
        endif
  110   continue
  100 continue
c
      return
      end
c
c//////////////////subroutine cfout/////////////////////////////////////
c    Output the DEM for optimal surface f(ig,jg) and/or the final
c    coefficients matrix c to defined files. To generate the DEM,
c    set the output region and coordinates of grid points (xg(ig),yg(jg)).
c    Output parameters are as below.
c
c          ----------(150,100)
c         |    |    |    |      In the case of
c         |----|----|----|            xmin,xmax=0,150
c         |    |    |    |            ymin,ymax=0,100
c       (0,0)------------             mx,my    =3,2
c-----------------------------------------------------------------------
c    Valiables
c       *ming,*maxg     :  edge of the output region
c       nxg,nyg         :  number of grid points
c       xg(nxg),yg(nyg) :  grid points
c       dx,dy           :  interval of grid
c       f               :  final solution
c///////////////////////////////////////////////////////////////////////
      subroutine cfout(mxgrd,mxr,xmin,xmax,ymin,ymax,mx,my,hx,hy,neq,c)
      implicit real*8 (a-h,o-z)
      dimension xg(mxgrd),yg(mxgrd),f(mxgrd,mxgrd)
      dimension c(mxr),bx0(4),by0(4)
      character*80 flg,flo
c
      write(*,*)
      write(*,10)
   10 format("#####Output result to files#####"/
     &       "File name for DEM or no (skip)     = ",$)
      read(*,11) flg
   11 format(a80)
      if((flg.eq."NO").or.(flg.eq."no").or.(flg.eq."No")) goto 400
c
      write(*,12)
   12 format("Output region      : x(min),x(max) = ",$)
      read(*,*) xming,xmaxg
      if((xming.lt.xmin).or.(xmaxg.gt.xmax)) goto 999
      write(*,13)
   13 format("                   : y(min),y(max) = ",$)
      read(*,*) yming,ymaxg
      if((yming.lt.ymin).or.(ymaxg.gt.ymax)) goto 999
      write(*,14) mxgrd
   14 format("Number of grid       Nx,Ny must be =< ",i4/
     &       "                   : Nx,Ny         = ",$)
      read(*,*) nxg,nyg
      if((nxg.lt.2).or.(nxg.gt.mxgrd)) goto 999
      if((nyg.lt.2).or.(nyg.gt.mxgrd)) goto 999
c
c----------Calculate the interval of grid-------------------------------
      dx=(xmaxg-xming)/dfloat(nxg-1)
      dy=(ymaxg-yming)/dfloat(nyg-1)
c
c----------Calculate the coodinates of grid points----------------------
      do 100 i=1,nxg
        xg(i)=xming+dfloat(i-1)*dx
  100 continue
      do 110 i=1,nyg
        yg(i)=yming+dfloat(i-1)*dy
  110 continue
c
c----------Calculate the height value f(ig,jg)--------------------------
      do 200 ig=1,nxg
        do 210 jg=1,nyg
          ipx=dint((xg(ig)-xmin)/hx)+1
          ipy=dint((yg(jg)-ymin)/hy)+1
          if(ipx.eq.mx+1) ipx=mx
          if(ipy.eq.my+1) ipy=my
          rx=(xg(ig)-xmin)/hx-dfloat(ipx-1)
          ry=(yg(jg)-ymin)/hy-dfloat(ipy-1)
          call bspl(rx,hx,0,bx0)
          call bspl(ry,hy,0,by0)
          f(ig,jg)=0.0D0
          do 220 j=0,3
            do 230 i=0,3
              ij=ipx+i+(ipy-1+j)*(mx+3)
              f(ig,jg)=f(ig,jg)+c(ij)*bx0(4-i)*by0(4-j)
  230       continue
  220     continue
  210   continue
  200 continue
c
c----------Output DEM---------------------------------------------------
      lu=60
      open(lu,file=flg)
      write(lu,*) nxg,nyg,xming,yming,dx,dy
      do 300 i=1,nxg
        write(lu,*) (f(i,j),j=1,nyg)
  300 continue
      close(lu)
c
  400 write(*,15)
   15 format("File name for optimal surface"/
     &       "                      or no (skip) = ",$)
      read(*,11) flo
      if((flo.eq."NO").or.(flo.eq."no").or.(flo.eq."No")) return
c
c----------Output the optimal surface-----------------------------------
      lu=61
      open(lu,file=flo)
      write(lu,*) mx,my,xmin,ymin,hx,hy
      do 410 i=1,neq
        write(lu,16) c(i)
   16   format(E16.9)
  410 continue
      close(lu)
      return
c
  999 write(*,99) 
   99 format("Invalid parameter input.")
      stop
c
      end
