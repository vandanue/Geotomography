c  ************************  elecfem3d.f  ****************************
c  BACKGROUND

c  This program solves Laplace's equation in a random conducting
c  material using the finite element method.  Each pixel in the 3-D digital
c  image is a cubic tri-linear finite element,  having its own conductivity.
c  Periodic boundary conditions are maintained.  In the comments below,
c  (USER) means that this is a section of code that the user might 
c  have to change for his particular problem.  Therefore the user is
c  encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the energy
c  1/2 uAu + b u + C, where A is the Hessian matrix composed of the 
c  stiffness matrices (dk) for each pixel/element, b is a constant vector
c  and C is a constant that are determined by the applied field and 
c  the periodic boundary conditions, and u is a vector of all the voltages.
c  The method used is the conjugate gradient relaxation algorithm.
c  Other variables are:  gb is the gradient = Au+b, h and Ah are 
c  auxiliary variables used in the conjugate gradient algorithm (in dembx),
c  dk(n,i,j) is the stiffness matrix of the n'th phase, sigma(n,i,j) is
c  the conductivity tensor of the n'th phase, pix is a vector that gives
c  the phase label of each pixel, ib is a matrix that gives the labels of
c  the 27 (counting itself) neighbors of a given node, prob is the volume 
c  fractions of the various phases, and currx, curry, currz are the 
c  volume averaged total currents in the x, y, and z directions.

c DIMENSIONS

c  The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
c  ns=nx*ny*nz, where the digital image of the microstructure considered
c  is a rectangular parallelipiped ( nx x ny x nz) in size.
c  The arrays pix and ib are also dimensioned to the system size.
c  The array ib has 27 components, for the 27 neighbors of a node.
c  Note that the program is set up at present to have at most 100 
c  different phases.  This can easily be changed, simply by changing 
c  the dimensions of dk, prob, and sigma. Nphase gives the number of
c  phases being considered.
c  All arrays are passed between subroutines using simple common statements.

c  STRONGLY SUGGESTED:  READ THE MANUAL BEFORE USING PROGRAM!!

c  (USER) Change these dimensions and in other subroutines at the same time.
c  For example, search and replace all occurrences throughout the program
c  of "(8000000" by "(64000", to go from a 20 x 20 x 20 system to a 
c  40 x 40 x 40 system.
	real u(8000000),gb(8000000),b(8000000),dk(100,8,8)
        real h(8000000),Ah(8000000)
	real sigma(100,3,3),prob(100)
	integer in(27),jn(27),kn(27)
	integer*4 ib(8000000,27)
        integer*2 pix(8000000)
	character*40 poreFileName
	integer iwriteFlag

	common/list1/currx,curry,currz
	common/list2/ex,ey,ez
	common/list3/ib
      common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb
	common/list8/sigma
	common/list9/h,Ah

c  (USER) Unit 9 is the microstructure input file, unit 7 is
c  the results output file.
	  open(unit=8, file='ec3d.pam')
c        open(9,file='microstructure.dat')
        open(7,file='outputfile.out')

c  (USER) nx,ny,nz give the size of the lattice
c      nx=20
c      ny=20
c      nz=20
c ns=total number of sites
		
	read(8,'(a40)') poreFileName
	print*, 'Pore File name: ',poreFileName

	read(8,*) nx, ny, nz
	print*, 'NX, NY, NZ    : ',nx, ny, nz

	read(8,*) gtest
	print*, 'Conv. Criteria: ',gtest

	read(8,*) kmax
	print*, 'Max Inversin  : ',kmax

	read(8,*) iwriteFlag
	print*, 'Write Local Current: ',iwriteFlag

	read(8,*) nphase
	print*, 'Number of Phase   : ',nphase

	
	do 6060 n=1,nphase
		read(8,*) sigma(n,1,1)
		print*, sigma(n,1,1)
6060  continue
	
	read(8,*) ex,ey,ez
	print*, 'Electric Filed Applied   : ',ex,ey,ez

	close(8)
	open (unit=9,file=poreFileName)



      ns=nx*ny*nz
      write(7,9010) nx,ny,nz,ns
9010  format('nx= ',i4,' ny= ',i4,' nz= ',i4,' ns = ',i8)

c  (USER) nphase is the number of phases being considered in the problem.
c  The values of pix(m) will run from 1 to nphase.
c      nphase=2

c  (USER) gtest is the stopping criterion, compared to gg=gb*gb.  
c  gtest=abc*ns, so that when gg < gtest, that average value per pixel
c  of gb is less than sqrt(abc).
      gtest=gtest*ns

c  Construct the neighbor table, ib(m,n)

c  First construct 27 neighbor table in terms of delta i, delta j, delta k
c  (See Table 3 in manual)
      in(1)=0
      in(2)=1
      in(3)=1
      in(4)=1
      in(5)=0
      in(6)=-1
      in(7)=-1
      in(8)=-1

      jn(1)=1
      jn(2)=1
      jn(3)=0
      jn(4)=-1
      jn(5)=-1
      jn(6)=-1
      jn(7)=0
      jn(8)=1

      do 555 n=1,8
      kn(n)=0
      kn(n+8)=-1
      kn(n+16)=1
      in(n+8)=in(n)
      in(n+16)=in(n)
      jn(n+8)=jn(n)
      jn(n+16)=jn(n)
555   continue
      in(25)=0
      in(26)=0
      in(27)=0
      jn(25)=0
      jn(26)=0
      jn(27)=0
      kn(25)=-1
      kn(26)=1
      kn(27)=0

c  Now construct neighbor table according to 1-d labels
c  Matrix ib(m,n) gives the 1-d label of the n'th neighbor (n=1,27) of
c  the node labelled m.
      nxy=nx*ny
      do 1020 k=1,nz
      do 1020 j=1,ny
      do 1020 i=1,nx
      m=nxy*(k-1)+nx*(j-1)+i
      do 1004 n=1,27
      i1=i+in(n)
      j1=j+jn(n)
      k1=k+kn(n)
      if(i1.lt.1) i1=i1+nx
      if(i1.gt.nx) i1=i1-nx
      if(j1.lt.1) j1=j1+ny
      if(j1.gt.ny) j1=j1-ny
      if(k1.lt.1) k1=k1+nz
      if(k1.gt.nz) k1=k1-nz
      m1=nxy*(k1-1)+nx*(j1-1)+i1
      ib(m,n)=m1
1004  continue
1020  continue 

c  Compute the electrical conductivity of each microstructure.
c  (USER) npoints is the number of microstructures to use.
        npoints=1
        do 8000 micro=1,npoints

c  Read in a microstructure in subroutine ppixel, and set up pix(m) 
c  with the appropriate phase assignments.
        call ppixel(nx,ny,nz,ns,nphase)
c Count and output the volume fractions of the different phases
        call assig(ns,nphase,prob)
	do 805 i=1,nphase
	write(7,*) 'Volume fraction of phase ',i,' = ',prob(i)
805     continue

c  (USER) sigma(100,3,3) is the electrical conductivity tensor of each phase
c  The user can make the value of sigma to be different for each 
c  phase of the microstructure if so desired (up to 100 phases as currently
c  dimensioned).

	do 9909 i=1,nphase

      sigma(i,1,2)=0.0
      sigma(i,1,3)=0.0
      sigma(i,2,2)=sigma(i,1,1)
      sigma(i,2,3)=0.0
      sigma(i,3,3)=sigma(i,1,1)
      sigma(i,2,1)=sigma(i,1,2)
      sigma(i,3,1)=sigma(i,1,3)
      sigma(i,3,2)=sigma(i,2,3)



9909  continue

c  write out the phase electrical conductivity tensors
      do 11 i=1,nphase 
      write(7,*) 'Phase ',i,' conductivity tensor is:'
      write(7,*) sigma(i,1,1),sigma(i,1,2),sigma(i,1,3)
      write(7,*) sigma(i,2,1),sigma(i,2,2),sigma(i,2,3)
      write(7,*) sigma(i,3,1),sigma(i,3,2),sigma(i,3,3)
11    continue

c  (USER) Set applied electric field.
      write(7,*) 'Applied field components:'
      write(7,*) 'ex = ',ex,' ey = ',ey,' ez = ',ez

c  Set up the finite element "stiffness" matrices and the Constant and
c  vector required for the energy

	call femat(nx,ny,nz,ns,nphase)

c  Apply a homogeneous macroscopic electric field as the initial condition
	do 1050 k=1,nz
	do 1050 j=1,ny
        do 1050 i=1,nx
		m=nxy*(k-1)+nx*(j-1)+i
		x=float(i-1)
		y=float(j-1)
		z=float(k-1)
		u(m)=-x*ex-y*ey-z*ez
1050	continue
	

c  Relaxation Loop
c  (USER) kmax is the maximum number of times dembx will be called, with
c  ldemb conjugate gradient steps done during each call. The total
c  number of conjugate gradient cycles allowed for a given conductivity
c  computation is kmax*ldemb.
c        kmax=40
        ldemb=50
        ltot=0

c  Call energy to get initial energy and initial gradient
        call energy(nx,ny,nz,ns,utot)
c  gg is the norm squared of the gradient (gg=gb*gb)
 	gg=0.0
        do 100 m=1,ns
        gg=gg+gb(m)*gb(m)
100     continue
	write(7,*) 'Initial energy = ',utot,'gg = ',gg
	print*, 'Initial energy = ',utot,' gg = ',gg
c       call flush(7)

        do 5000 kkk=1,kmax
c  Call dembx to go into conjugate gradient solver
        call dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
        ltot=ltot+Lstep
c  Call energy to compute energy after dembx call. If gg < gtest, this
c  will be the final energy. If gg is still larger than gtest, then this
c  will give an intermediate energy with which to check how the relaxation 
c  process is coming along.
        call energy(nx,ny,nz,ns,utot)
	write(7,*) 'Energy = ',utot,'gg = ',gg
	write(7,*) ltot, ' conj. grad. steps'
	print*, 'Energy = ',utot,' gg = ',gg
	print*, 'Number of conjugate steps  = ',ltot
        if(gg.lt.gtest) goto 444

c  If relaxation process will continue, compute and output currents as an
c  additional aid to judge how the relaxation procedure if progressing.
       call current(nx,ny,nz,ns,0)
c Output intermediate currents
        write(7,*)
	write(7, *) ' Current in x direction = ',currx
	write(7, *) ' Current in y direction = ',curry
	write(7, *) ' Current in z direction = ',currz
	
c        call flush(7)

5000   continue

444    call current(nx,ny,nz,ns,iwriteFlag)

c Output final currents
        write(7,*)
	write(7, *) ' Current in x direction = ',currx
	write(7, *) ' Current in y direction = ',curry
	write(7, *) ' Current in z direction = ',currz
	print*, ' Current in x direction = ',currx
	print*, ' Current in y direction = ',curry
	print*, ' Current in z direction = ',currz
c        call flush(7)

8000    continue

        end

c  Subroutine that sets up the stiffness matrices, linear term in 
c  voltages, and constant term C that appear in the total energy due 
c  to the periodic boundary conditions.

      subroutine femat(nx,ny,nz,ns,nphase)
      real dk(100,8,8),xn(8),b(8000000),C
      real dndx(8),dndy(8),dndz(8)
      real g(3,3,3),sigma(100,3,3)
      real es(3,8)
      integer is(8)
      integer*4 ib(8000000,27)
      integer*2 pix(8000000)

	common/list2/ex,ey,ez
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list8/sigma

      nxy=nx*ny

c  initialize stiffness matrices
      do 40 m=1,nphase
      do 40 j=1,8
      do 40 i=1,8
      dk(m,i,j)=0.0
40    continue

c  set up Simpson's rule integration weight vector
      do 30 k=1,3
      do 30 j=1,3
      do 30 i=1,3
      nm=0
      if(i.eq.2) nm=nm+1
      if(j.eq.2) nm=nm+1
      if(k.eq.2) nm=nm+1
      g(i,j,k)=4.0**nm
30    continue

c  loop over the nphase kinds of pixels and Simpson's rule quadrature
c  points in order to compute the stiffness matrices.  Stiffness matrices
c  of trilinear finite elements are quadratic in x, y, and z, so that
c  Simpson's rule quadrature is exact.
      do 4000 ijk=1,nphase
      do 3000 k=1,3
      do 3000 j=1,3
      do 3000 i=1,3
      x=float(i-1)/2.0
      y=float(j-1)/2.0
      z=float(k-1)/2.0
c  dndx means the negative derivative with respect to x of the shape 
c  matrix N (see manual, Sec. 2.2), dndy, dndz are similar.
      dndx(1)=-(1.0-y)*(1.0-z)
      dndx(2)=(1.0-y)*(1.0-z)
      dndx(3)=y*(1.0-z)
      dndx(4)=-y*(1.0-z)
      dndx(5)=-(1.0-y)*z
      dndx(6)=(1.0-y)*z
      dndx(7)=y*z
      dndx(8)=-y*z
      dndy(1)=-(1.0-x)*(1.0-z)
      dndy(2)=-x*(1.0-z)
      dndy(3)=x*(1.0-z)
      dndy(4)=(1.0-x)*(1.0-z)
      dndy(5)=-(1.0-x)*z
      dndy(6)=-x*z
      dndy(7)=x*z
      dndy(8)=(1.0-x)*z
      dndz(1)=-(1.0-x)*(1.0-y)
      dndz(2)=-x*(1.0-y)
      dndz(3)=-x*y
      dndz(4)=-(1.0-x)*y
      dndz(5)=(1.0-x)*(1.0-y)
      dndz(6)=x*(1.0-y)
      dndz(7)=x*y
      dndz(8)=(1.0-x)*y
c  now build electric field matrix
      do 2799 n1=1,3
      do 2799 n2=1,8
      es(n1,n2)=0.0
2799  continue
      do 2797 n=1,8
      es(1,n)=dndx(n)
      es(2,n)=dndy(n)
      es(3,n)=dndz(n)
2797  continue
c  now do matrix multiply to determine value at (x,y,z), multiply by
c  proper weight, and sum into dk, the stiffness matrix
      do 900 ii=1,8
      do 900 jj=1,8
c  Define sum over field matrices and conductivity tensor that defines
c  the stiffness matrix.
      sum=0.0
      do 890 kk=1,3
      do 890 ll=1,3
      sum=sum+es(kk,ii)*sigma(ijk,kk,ll)*es(ll,jj)
890   continue
      dk(ijk,ii,jj)=dk(ijk,ii,jj)+g(i,j,k)*sum/216.
900   continue

3000  continue
4000  continue

c  Set up vector for linear term, b, and constant term, C, 
c  in the electrical energy.  This is done using the stiffness matrices,
c  and the periodic terms in the applied field that come in at the boundary
c  pixels via the periodic boundary conditions and the condition that
c  an applied macroscopic field exists (see Sec. 2.2 in manual).

      do 5000 m=1,ns
      b(m)=0.0
5000  continue

c  For all cases, correspondence between 1-8 finite element node labels
c  and 1-27 neighbor labels is:  1:ib(m,27),2:ib(m,3),3:ib(m,2),
c  4:ib(m,1),5:ib(m,26),6:ib(m,19),7:ib(m,18),8:ib(m,17) 
c  (see Table 4 in manual)
      is(1)=27
      is(2)=3
      is(3)=2
      is(4)=1
      is(5)=26
      is(6)=19
      is(7)=18
      is(8)=17

      C=0.0
c  x=nx face
      i=nx
      do 2001 i8=1,8
      xn(i8)=0.0
      if(i8.eq.2.or.i8.eq.3.or.i8.eq.6.or.i8.eq.7) then
      xn(i8)=-ex*nx
      end if
2001  continue
      do 2000 j=1,ny-1
      do 2000 k=1,nz-1
      m=nxy*(k-1)+j*nx
      do 1900 mm=1,8
      sum=0.0
      do 1899 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
1899  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1900  continue
2000  continue
c  y=ny face
      j=ny
      do 2011 i8=1,8
      xn(i8)=0.0
      if(i8.eq.3.or.i8.eq.4.or.i8.eq.7.or.i8.eq.8) then
      xn(i8)=-ey*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      do 2010 k=1,nz-1
      m=nxy*(k-1)+nx*(ny-1)+i
      do 1901 mm=1,8
      sum=0.0
      do 2099 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2099  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1901  continue
2010  continue
c  z=nz face
      k=nz
      do 2021 i8=1,8
      xn(i8)=0.0
      if(i8.eq.5.or.i8.eq.6.or.i8.eq.7.or.i8.eq.8) then
      xn(i8)=-ez*nz
      end if
2021  continue
      do 2020 i=1,nx-1
      do 2020 j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+i
      do 1902 mm=1,8
      sum=0.0
      do 2019 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2019  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1902  continue
2020  continue
c  x=nx y=ny edge 
      i=nx
      y=ny
      do 2031 i8=1,8
      xn(i8)=0.0
      if(i8.eq.2.or.i8.eq.6) then
      xn(i8)=-ex*nx
      end if
      if(i8.eq.4.or.i8.eq.8) then
      xn(i8)=-ey*ny
      end if
      if(i8.eq.3.or.i8.eq.7) then
      xn(i8)=-ey*ny-ex*nx
      end if
2031  continue
      do 2030 k=1,nz-1
      m=nxy*k
      do 1903 mm=1,8
      sum=0.0
      do 2029 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2029  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1903  continue
2030  continue
c  x=nx z=nz edge 
      i=nx
      k=nz
      do 2041 i8=1,8
      xn(i8)=0.0
      if(i8.eq.2.or.i8.eq.3) then
      xn(i8)=-ex*nx
      end if
      if(i8.eq.5.or.i8.eq.8) then
      xn(i8)=-ez*nz
      end if
      if(i8.eq.6.or.i8.eq.7) then
      xn(i8)=-ez*nz-ex*nx
      end if
2041  continue
      do 2040 j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+nx
      do 1904 mm=1,8
      sum=0.0
      do 2039 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2039  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1904  continue
2040  continue
c  y=ny z=nz edge 
      j=ny
      k=nz
      do 2051 i8=1,8
      xn(i8)=0.0
      if(i8.eq.5.or.i8.eq.6) then
      xn(i8)=-ez*nz
      end if
      if(i8.eq.3.or.i8.eq.4) then
      xn(i8)=-ey*ny
      end if
      if(i8.eq.7.or.i8.eq.8) then
      xn(i8)=-ey*ny-ez*nz
      end if
2051  continue
      do 2050 i=1,nx-1
      m=nxy*(nz-1)+nx*(ny-1)+i
      do 1905 mm=1,8
      sum=0.0
      do 2049 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2049  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1905  continue
2050  continue
c  x=nx y=ny z=nz corner 
      i=nx
      j=ny
      k=nz
      do 2061 i8=1,8
      xn(i8)=0.0
      if(i8.eq.2) then
      xn(i8)=-ex*nx
      end if
      if(i8.eq.4) then
      xn(i8)=-ey*ny
      end if
      if(i8.eq.5) then
      xn(i8)=-ez*nz
      end if
      if(i8.eq.8) then
      xn(i8)=-ey*ny-ez*nz
      end if
      if(i8.eq.6) then
      xn(i8)=-ex*nx-ez*nz
      end if
      if(i8.eq.3) then
      xn(i8)=-ex*nx-ey*ny
      end if
      if(i8.eq.7) then
      xn(i8)=-ex*nx-ey*ny-ez*nz
      end if
2061  continue
      m=nx*ny*nz
      do 1906 mm=1,8
      sum=0.0
      do 2059 m8=1,8
      sum=sum+xn(m8)*dk(pix(m),m8,mm)
      C=C+0.5*xn(m8)*dk(pix(m),m8,mm)*xn(mm)
2059  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1906  continue

      return
      end

c  Subroutine computes the total energy, utot, and gradient, gb

      subroutine energy(nx,ny,nz,ns,utot)
	real u(8000000),gb(8000000)
	real b(8000000),C
	real dk(100,8,8)
	real utot
 	integer*4 ib(8000000,27)
        integer*2 pix(8000000)
	
	common/list2/ex,ey,ez
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb

	do 2090 m=1,ns
	gb(m)=0.0
2090	continue

c Energy loop. Do global matrix multiply via small stiffness matrices,
c  gb=Au + b.  The long statement below correctly brings in all the
c  terms from the global matrix A using only the small stiffness matrices.

      do 3000 m=1,ns
      gb(m)=gb(m)+u(ib(m,1))*( dk(pix(ib(m,27)),1,4)+
     &dk(pix(ib(m,7)),2,3)+
     &dk(pix(ib(m,25)),5,8)+dk(pix(ib(m,15)),6,7) )+
     &u(ib(m,2))*( dk(pix(ib(m,27)),1,3)+dk(pix(ib(m,25)),5,7) )+
     &u(ib(m,3))*( dk(pix(ib(m,27)),1,2)+dk(pix(ib(m,5)),4,3)+
     &dk(pix(ib(m,13)),8,7)+dk(pix(ib(m,25)),5,6) )+
     &u(ib(m,4))*( dk(pix(ib(m,5)),4,2)+dk(pix(ib(m,13)),8,6) )+
     &u(ib(m,5))*( dk(pix(ib(m,6)),3,2)+dk(pix(ib(m,5)),4,1)+
     &dk(pix(ib(m,14)),6,7)+dk(pix(ib(m,13)),8,5) )+
     &u(ib(m,6))*( dk(pix(ib(m,6)),3,1)+dk(pix(ib(m,14)),7,5) )+
     &u(ib(m,7))*( dk(pix(ib(m,6)),3,4)+dk(pix(ib(m,7)),2,1)+
     &dk(pix(ib(m,14)),7,8)+dk(pix(ib(m,15)),6,5) )+
     &u(ib(m,8))*( dk(pix(ib(m,7)),2,4)+dk(pix(ib(m,15)),6,8) )+
     &u(ib(m,9))*( dk(pix(ib(m,25)),5,4)+dk(pix(ib(m,15)),6,3) )+
     &u(ib(m,10))*( dk(pix(ib(m,25)),5,3) )+
     &u(ib(m,11))*( dk(pix(ib(m,13)),8,3)+dk(pix(ib(m,25)),5,2) )+
     &u(ib(m,12))*( dk(pix(ib(m,13)),8,2) )+
     &u(ib(m,13))*( dk(pix(ib(m,13)),8,1)+dk(pix(ib(m,14)),7,2) )+
     &u(ib(m,14))*( dk(pix(ib(m,14)),7,1) )+
     &u(ib(m,15))*( dk(pix(ib(m,14)),7,4)+dk(pix(ib(m,15)),6,1) )+
     &u(ib(m,16))*( dk(pix(ib(m,15)),6,4) )+
     &u(ib(m,17))*( dk(pix(ib(m,27)),1,8)+dk(pix(ib(m,7)),2,7) )+
     &u(ib(m,18))*( dk(pix(ib(m,27)),1,7) )+
     &u(ib(m,19))*( dk(pix(ib(m,27)),1,6)+dk(pix(ib(m,5)),4,7) )+
     &u(ib(m,20))*( dk(pix(ib(m,5)),4,6) )+
     &u(ib(m,21))*( dk(pix(ib(m,5)),4,5)+dk(pix(ib(m,6)),3,6) )+
     &u(ib(m,22))*( dk(pix(ib(m,6)),3,5) )+
     &u(ib(m,23))*( dk(pix(ib(m,6)),3,8)+dk(pix(ib(m,7)),2,5) )+
     &u(ib(m,24))*( dk(pix(ib(m,7)),2,8) )+
     &u(ib(m,25))*( dk(pix(ib(m,14)),7,3)+dk(pix(ib(m,13)),8,4)+
     &dk(pix(ib(m,15)),6,2)+dk(pix(ib(m,25)),5,1) )+
     &u(ib(m,26))*( dk(pix(ib(m,6)),3,7)+dk(pix(ib(m,5)),4,8)+
     &dk(pix(ib(m,27)),1,5)+dk(pix(ib(m,7)),2,6) )+
     &u(ib(m,27))*( dk(pix(ib(m,27)),1,1)+dk(pix(ib(m,7)),2,2)+
     &dk(pix(ib(m,6)),3,3)+dk(pix(ib(m,5)),4,4)+dk(pix(ib(m,25)),5,5)+
     &dk(pix(ib(m,15)),6,6)+dk(pix(ib(m,14)),7,7)+
     &dk(pix(ib(m,13)),8,8) )
3000  continue

	utot=0.0
	do 3100 m=1,ns
	utot=utot+0.5*u(m)*gb(m)+b(m)*u(m)
	gb(m)=gb(m)+b(m)
3100	continue

	utot=utot+C

        return
        end           	       

c    Subroutine that carries out the conjugate gradient relaxation process

      subroutine dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
      real u(8000000),gb(8000000),dk(100,8,8)    
      real Ah(8000000),h(8000000),B,lambda,gamma
      integer*4 ib(8000000,27)
      integer*2 pix(8000000)

	common/list3/ib
        common/list4/pix
	common/list6/u
	common/list7/gb
	common/list9/h,Ah
      
c  Initialize the conjugate direction vector on first call to dembx only.
c  For calls to dembx after the first, we want to continue using the
c  value fo h determined in the previous call.  Of course, if npooints
c  is greater than 1, then this initialization step will be run every
c  a new microstructure is used, as kkk will be reset to 1 every time 
c  the counter micro is increased.
      if(kkk.eq.1) then
      do 50 m=1,ns
      h(m)=gb(m) 
50    continue                                                          
      end if
c  Lstep counts the number of conjugate gradient steps taken in each call
c  to dembx.  
      Lstep=0                                                              

c     Conjugate gradient loop

      do 800 ijk=1,ldemb
      Lstep=Lstep+1
          
      do 290 m=1,ns
      Ah(m)=0.0
290   continue

c  Do global matrix multiply via small stiffness matrices, Ah = A * h.
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices.

      do 400 m=1,ns
      Ah(m)=Ah(m)+h(ib(m,1))*( dk(pix(ib(m,27)),1,4)+
     &dk(pix(ib(m,7)),2,3)+
     &dk(pix(ib(m,25)),5,8)+dk(pix(ib(m,15)),6,7) )+
     &h(ib(m,2))*( dk(pix(ib(m,27)),1,3)+dk(pix(ib(m,25)),5,7) )+
     &h(ib(m,3))*( dk(pix(ib(m,27)),1,2)+dk(pix(ib(m,5)),4,3)+
     &dk(pix(ib(m,13)),8,7)+dk(pix(ib(m,25)),5,6) )+
     &h(ib(m,4))*( dk(pix(ib(m,5)),4,2)+dk(pix(ib(m,13)),8,6) )+
     &h(ib(m,5))*( dk(pix(ib(m,6)),3,2)+dk(pix(ib(m,5)),4,1)+
     &dk(pix(ib(m,14)),6,7)+dk(pix(ib(m,13)),8,5) )+
     &h(ib(m,6))*( dk(pix(ib(m,6)),3,1)+dk(pix(ib(m,14)),7,5) )+
     &h(ib(m,7))*( dk(pix(ib(m,6)),3,4)+dk(pix(ib(m,7)),2,1)+
     &dk(pix(ib(m,14)),7,8)+dk(pix(ib(m,15)),6,5) )+
     &h(ib(m,8))*( dk(pix(ib(m,7)),2,4)+dk(pix(ib(m,15)),6,8) )+
     &h(ib(m,9))*( dk(pix(ib(m,25)),5,4)+dk(pix(ib(m,15)),6,3) )+
     &h(ib(m,10))*( dk(pix(ib(m,25)),5,3) )+
     &h(ib(m,11))*( dk(pix(ib(m,13)),8,3)+dk(pix(ib(m,25)),5,2) )+
     &h(ib(m,12))*( dk(pix(ib(m,13)),8,2) )+
     &h(ib(m,13))*( dk(pix(ib(m,13)),8,1)+dk(pix(ib(m,14)),7,2) )+
     &h(ib(m,14))*( dk(pix(ib(m,14)),7,1) )+
     &h(ib(m,15))*( dk(pix(ib(m,14)),7,4)+dk(pix(ib(m,15)),6,1) )+
     &h(ib(m,16))*( dk(pix(ib(m,15)),6,4) )+
     &h(ib(m,17))*( dk(pix(ib(m,27)),1,8)+dk(pix(ib(m,7)),2,7) )+
     &h(ib(m,18))*( dk(pix(ib(m,27)),1,7) )+
     &h(ib(m,19))*( dk(pix(ib(m,27)),1,6)+dk(pix(ib(m,5)),4,7) )+
     &h(ib(m,20))*( dk(pix(ib(m,5)),4,6) )+
     &h(ib(m,21))*( dk(pix(ib(m,5)),4,5)+dk(pix(ib(m,6)),3,6) )+
     &h(ib(m,22))*( dk(pix(ib(m,6)),3,5) )+
     &h(ib(m,23))*( dk(pix(ib(m,6)),3,8)+dk(pix(ib(m,7)),2,5) )+
     &h(ib(m,24))*( dk(pix(ib(m,7)),2,8) )+
     &h(ib(m,25))*( dk(pix(ib(m,14)),7,3)+dk(pix(ib(m,13)),8,4)+
     &dk(pix(ib(m,15)),6,2)+dk(pix(ib(m,25)),5,1) )+
     &h(ib(m,26))*( dk(pix(ib(m,6)),3,7)+dk(pix(ib(m,5)),4,8)+
     &dk(pix(ib(m,27)),1,5)+dk(pix(ib(m,7)),2,6) )+
     &h(ib(m,27))*( dk(pix(ib(m,27)),1,1)+dk(pix(ib(m,7)),2,2)+
     &dk(pix(ib(m,6)),3,3)+dk(pix(ib(m,5)),4,4)+dk(pix(ib(m,25)),5,5)+
     &dk(pix(ib(m,15)),6,6)+dk(pix(ib(m,14)),7,7)+
     &dk(pix(ib(m,13)),8,8) )
 
400   continue

      hAh=0.0                                                            
      do 530 m=1,ns                                                    
      hAh=hAh+h(m)*Ah(m)                                                 
530   continue                                                          

      lambda=gg/hAh
      do 540 m=1,ns                                                    
      u(m)=u(m)-lambda*h(m)                                        
      gb(m)=gb(m)-lambda*Ah(m)                                        
540   continue
                                                                
      gglast=gg
      gg=0.0
      do 550 m=1,ns                                                     
      gg=gg+gb(m)*gb(m)       
550   continue
      if(gg.le.gtest) goto 1000

      gamma=gg/gglast
      do 570 m=1,ns                                                   
      h(m)=gb(m)+gamma*h(m)                                       
570   continue                                                        

800   continue

1000  continue                                                        

      return                                                           
      end        

c Subroutine that computes average current in three directions

      subroutine current(nx,ny,nz,ns,iwriteFlag)

      real af(3,8)
      real u(8000000),uu(8)
      real sigma(100,3,3)
      integer*4 ib(8000000,27)
      integer*2 pix(8000000)
	integer*4 nnss;

	common/list1/currx,curry,currz
	common/list2/ex,ey,ez
	common/list3/ib
        common/list4/pix
	common/list6/u
	common/list8/sigma

      nxy=nx*ny
	nnss = (nx-2)*(ny-2)*(nz-2);
c  af is the average field matrix, average field in a pixel is af*u(pixel).
c  The matrix af relates the nodal voltages to the average field in the pixel.

c Set up single element average field matrix

      af(1,1)=0.25
      af(1,2)=-0.25
      af(1,3)=-0.25
      af(1,4)=0.25
      af(1,5)=0.25
      af(1,6)=-0.25
      af(1,7)=-0.25
      af(1,8)=0.25
      af(2,1)=0.25
      af(2,2)=0.25
      af(2,3)=-0.25
      af(2,4)=-0.25
      af(2,5)=0.25
      af(2,6)=0.25
      af(2,7)=-0.25
      af(2,8)=-0.25
      af(3,1)=0.25
      af(3,2)=0.25
      af(3,3)=0.25
      af(3,4)=0.25
      af(3,5)=-0.25
      af(3,6)=-0.25
      af(3,7)=-0.25
      af(3,8)=-0.25

c  now compute current in each pixel
      currx=0.0
      curry=0.0
      currz=0.0
c  compute average field in each pixel

	if(iwriteFlag .eq. 1) open(unit=3, file='currField.dat')

      do 470 k=2,nz-1
      do 470 j=2,ny-1
      do 470 i=2,nx-1
      m=(k-1)*nxy+(j-1)*nx+i
c  load in elements of 8-vector using pd. bd. conds.
      uu(1)=u(m) 
      uu(2)=u(ib(m,3))
      uu(3)=u(ib(m,2))
      uu(4)=u(ib(m,1))
      uu(5)=u(ib(m,26))
      uu(6)=u(ib(m,19))
      uu(7)=u(ib(m,18))
      uu(8)=u(ib(m,17))
c  Correct for periodic boundary conditions, some voltages are wrong
c  for a pixel on a periodic boundary. Since they come from an opposite
c  face, need to put in applied fields to correct them.
      if(i.eq.nx) then 
      uu(2)=uu(2)-ex*nx
      uu(3)=uu(3)-ex*nx
      uu(6)=uu(6)-ex*nx
      uu(7)=uu(7)-ex*nx
      end if
      if(j.eq.ny) then
      uu(3)=uu(3)-ey*ny
      uu(4)=uu(4)-ey*ny
      uu(7)=uu(7)-ey*ny
      uu(8)=uu(8)-ey*ny
      end if
      if(k.eq.nz) then
      uu(5)=uu(5)-ez*nz
      uu(6)=uu(6)-ez*nz
      uu(7)=uu(7)-ez*nz
      uu(8)=uu(8)-ez*nz
      end if
c  cur1, cur2, cur3 are the local currents averaged over the pixel
      cur1=0.0
      cur2=0.0
      cur3=0.0
      do 465 n=1,8
      do 465 nn=1,3
      cur1=cur1+sigma(pix(m),1,nn)*af(nn,n)*uu(n)
      cur2=cur2+sigma(pix(m),2,nn)*af(nn,n)*uu(n)
      cur3=cur3+sigma(pix(m),3,nn)*af(nn,n)*uu(n)
465   continue
c  sum into the global average currents
      currx=currx+cur1
      curry=curry+cur2
      currz=currz+cur3

	if(iwriteFlag .eq. 1) write(3,6070) cur1, cur2, cur3
6070	format(3(e12.5, 2x))

470   continue

c Volume average currents
      currx=currx/float(nnss)
      curry=curry/float(nnss)
      currz=currz/float(nnss)

	if(iwriteFlag .eq. 1) close(3)

      return
      end

c  Subroutine that counts phase volume fractions

      subroutine assig(ns,nphase,prob)

      integer ns,nphase
      integer*2 pix(8000000)
      real prob(100)
      common/list4/pix

	do 90 i=1,nphase
	prob(i)=0.0
90      continue

      do 100 m=1,ns
	do 100 i=1,nphase
        if(pix(m).eq.i) then
	prob(i)=prob(i)+1
	endif
100   continue

	do 110 i=1,nphase
	prob(i)=prob(i)/float(ns)
110     continue

      return
      end

c  Subroutine that sets up microstructural image

      subroutine ppixel(nx,ny,nz,ns,nphase)
      integer*2 pix(8000000)
      common/list4/pix
c	integer*2 itmp

c  (USER) If you want to set up a test image inside the program, instead
c  of reading it in from a file, this should be done inside this subroutine.

      do 100 k=1,nz
      do 100 j=1,ny
      do 100 i=1,nx 
      m=nx*ny*(k-1)+nx*(j-1)+i
      read(9,*) pix(m)
	pix(m)=pix(m)+1;
100   continue

c  Check for wrong phase labels--less than 1 or greater than nphase
       do 500 m=1,ns
       if(pix(m).lt.1) then
        write(7,*) 'Phase label in pix < 1--error at ',m
       end if
       if(pix(m).gt.nphase) then
	  write(*,*) pix(m)
        write(7,*) 'Phase label in pix > nphase--error at ',m
       end if
500    continue

      return
      end
