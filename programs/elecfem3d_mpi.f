c  *********************** elecfem3d_mpi.f **************************

c  This is the new MPI version of the elecfem3d.f code from
c  Section 9.3.1 of NISTIR 6269.

c  The main differences with this code compared to the serial
c  version are:

c  1. Removal of ib array.
c  2. Change of dimensionality on pix from pix(m) to pix(i,j,k)
c  Maximum value of m = nx*ny*nz (nx,ny,nz are the array dims).
c  3. All important arrays (pix,vox,gb,b,u,h,Ah) are dynamically allocated.
c
c  IN THIS VERSION:

c  The USER needs the following input:
c  (Search for occurences of "USER" in the code).

c  1. A 3-D pixel value data file with input & output names.
c  2. The values of the 3 dimensions: (nx,ny,nz)
c  3. The number of phases in the mixture: nphase
c  4. A convergence value: gtest
c  5. Applied electric field: ex,ey,ez
c  6. Values for DEMBX_MPI and how long it will run: kmax & ldemb

c  7. Flag for printing timing info for all data
c  passing MPI routines ( FEMAT_MPI, ENERGY_MPI, DEMBX_MPI)
c  from MAIN is called: pflag
c  pflag Values = 0,1 0=no timing info; 1=print timing info

c  pflag is a common value.
c
c  Timing info for the RELAXATION loop is not
c  influenced by the pflag and will always be printed.

c  User may edit the code to supress the printing.

c  8. Timing info stored in arrays namex X_time(i)
c  Where X=n,f,e ie.
c  n_time is in MAIN
c  f_time is in FEMAT_MPI
c  e_time is in ENERGY_MPI

c  NB: One also needs to insure that the values for
c  phasemod(i,j) are initialized correctly in
c  SUBROUTINE phasemod_init.


c  END of NEW comments.

c  BEGIN ORIGINAL comments.

c  BACKGROUND


c  This program solves Laplace's equation in a random conducting
c  material using the finite element method. Each pixel in the 3-D digital
c  image is a cubic tri-linear finite element, having its own conductivity.
c  Periodic boundary conditions are maintained. In the comments below,
c  (USER) means that this is a section of code that the user might
c  have to change for his particular problem. Therefore the user is
c  encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the energy
c  1/2 uAu + b u + C, where A is the Hessian matrix composed of the
c  stiffness matrices (dk) for each pixel/element, b is a constant vector
c and C is a constant that are determined by the applied field and
c  the periodic boundary conditions, and u is a vector of all the voltages.
c  The method used is the conjugate gradient relaxation algorithm.
c  Other variables are: gb is the gradient = Au+b, h and Ah are
c  auxiliary variables used in the conjugate gradient algorithm (in dembx),
c  dk(n,i,j) is the stiffness matrix of the n'th phase, sigma(n,i,j) is
c  the conductivity tensor of the n'th phase, pix is a vector that gives
c  the phase label of each pixel, ib is a matrix that gives the labels of
c  the 27 (counting itself) neighbors of a given node, prob is the volume
c  fractions of the various phases, and currx, curry, currz are the
c  volume averaged total currents in the x, y, and z directions.

c  DIMENSIONS

c  The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
c  ns=nx*ny*nz, where the digital image of the microstructure considered
c  is a rectangular parallelipiped ( nx x ny x nz) in size.
c  The arrays pix and ib are also dimensioned to the system size.
c  The array ib has 27 components, for the 27 neighbors of a node.
c  Note that the program is set up at present to have at most 100
c  different phases. This can easily be changed, simply by changing
c  the dimensions of dk, prob, and sigma. Nphase gives the number of
c  phases being considered.
c  All arrays are passed between subroutines using simple common statements.

c  STRONGLY SUGGESTED: READ THE MANUAL BEFORE USING PROGRAM!!
      implicit none
      include 'mpif.h'
	  
c  (USER) Change the nx,ny,nz dimensions at the beginning.
c  All important arrays are dynamically allocated.

      integer*2, allocatable :: dat(:,:,:), datn(:,:,:)
      integer*2, allocatable :: pix(:,:,:), pixn(:,:,:)
      integer*2, allocatable :: vox(:,:,:)

      integer, allocatable :: d1s(:),d2s(:)

      double precision, allocatable :: b(:,:,:,:)
      double precision, allocatable :: gb(:,:,:,:)
      double precision, allocatable :: u(:,:,:,:)
      double precision, allocatable :: h(:,:,:,:)

      double precision, allocatable :: sigma(:,:,:), prob(:)
      double precision, allocatable :: dk(:,:,:)

      double precision dgg,gg,utot,gtest,C
      double precision ex,ey,ez
      double precision x,y,z,saves

      double precision cuxxp,cuyyp,cuzzp
      double precision currx,curry,currz

      integer d1,d2,ns,sxip,kkk,mxy
      integer i,j,k,nx,ny,nz,nxy,nphase
      integer count,rem
      integer sz,sized
      integer npoints,micro,m
      integer kmax,ldemb,ltot,lstep
      integer pflag

      integer myrank,ierr,nprocs,irank
      integer status(MPI_STATUS_SIZE)

      double precision starttime,endtime, start_npoint, end_npoint
      double precision kkk_start,kkk_end
      double precision elapsed_time,stress_loop
      double precision n_time(24)

      common/list1/pflag,nphase
      common/list2/ex,ey,ez
      common/list3/currx,curry,currz
      common/list4/cuxxp,cuyyp,cuzzp

      call MPI_INIT(ierr)

      starttime = MPI_Wtime(ierr)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if (myrank.eq.0) then
      write(*,*) "There are ",nprocs," processors running this job."
      end if

c  USER: Change nx,ny,nz,nphase values to match your data.

      nx=512
      ny=512
      nz=512
      nphase=2

      nxy=nx*ny
      ns=nx*ny*nz
      sz=nz/nprocs
      mxy= 3*nx*ny

      gtest=1.d-10*ns
	  
c  pflag=0 for no timing info printed.
c  pflag=1 for timing info printed.
      pflag = 0

c  End this USER section.

      utot =0.0d0

      allocate( sigma(nphase,3,3) )
      allocate( dk(nphase,8,8) )
      allocate( prob(nphase) )

c  (USER) sigma(nphase,3,3) is the electrical conductivity tensor of each phase
c  The user can make the value of sigma to be different for each
c  phase of the microstructure if so desired.

c  Only diagonal elements need values

      sigma=0.0d0

      sigma(1,1,1)=1.0d0
      sigma(1,2,2)=1.0d0
      sigma(1,3,3)=1.0d0

      sigma(2,1,1)=200.0d0
      sigma(2,2,2)=200.0d0
      sigma(2,3,3)=200.0d0

c  (USER) Set applied electric field.
      ex=1.0d0
      ey=1.0d0
      ez=1.0d0

c  Calculate d1s(n) & d2s(n); These hold the d1 and d2
c  values for processor n.

      if (myrank.eq.0) then

      allocate (d1s(0:nprocs-1))
      allocate (d2s(0:nprocs-1))

      do irank=0,nprocs-1
            d1s(irank)=irank*sz+1
            d2s(irank)=(irank+1)*sz
      end do
		
      rem = nz - nprocs*sz
		
      if (rem.ne.0) then
            do j=1,rem
                  irank=nprocs-rem+j-1
                  d1s(irank)=d1s(irank)+ j-1
                  d2s(irank)=d2s(irank)+ j
            end do 
      end if 

c  Send all d1s(i) and d2s(i) from ROOT
c  to NODE i & store into d1 & d2

      do i=0,nprocs-1
      call MPI_SEND(d1s(i),1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(d2s(i),1,MPI_INTEGER,i,1,MPI_COMM_WORLD,ierr)
      end do

      end if 

      call MPI_RECV(d1,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(d2,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierr)
      write(*,*) "Rank#",myrank,"d1= ",d1," d2= ",d2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c  Allocate other arrays which need d1&d2 values.

      allocate (gb(nx,ny,d1-1:d2+1))
      gb=0.0d0
      allocate(b(nx,ny,d1-1:d2+1))
      b = 0.0d0

      allocate (u(nx,ny,d1-1:d2+1))
      allocate (h(nx,ny,d1-1:d2+1))

c  Want the ability to calculate on a series
c  of input files based on a value & some if statements.

c  Compute the average stress and strain in each microstructure.
c  (USER) npoints is the number of microstructures to use.

      npoints=1

c  (USER) Unit 9 is the microstructure input file,
c  unit 7 is the results output file.

      n_time(1) = MPI_Wtime(ierr)

      do micro=1,npoints

c  Allocate pix, so root can read it.

      if (myrank.eq.0) then
            allocate (pix(nx,ny,nz))
      end if

      start_npoint=MPI_Wtime(ierr)
      n_time(2) = MPI_Wtime(ierr)

      if (myrank.eq.0) then

c  Get pix from input file (unit=9).

c  (USER) Unit 9 is the microstructure input file, unit 7 is
c  the results output file.

      open(9,file='512_6.txt')
      open(7,file='512_elecmpi.out')

      write(*,*) "MICRO = ", micro
      write(7,*) "MICRO = ", micro

c  Finally... read in pix

      write(*,*) "call dpixel"
      Call dpixel(nx,ny,nz,ns,pix)
      write(*,*) "back from dpixel"

c  ns=total number of sites
      write(7,9010) nx,ny,nz,ns,nprocs
9010  format('nx= ',i4,' ny= ',i4,' nz= ',i4,' ns= 'i8,' nprocs= ',i4)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c  Now that the nodes are set up correctly,
c  one can pass the data from the root node (myrank=0)
c  to all the rest.

      allocate(dat(nx,ny,d1:d2))      
      sized = SIZE(dat)
      dat=0

      n_time(3)=MPI_Wtime(ierr)

      if (nprocs.eq.1) then
            dat=pix
            write(*,*) "dat=pix"
      end if

      if (nprocs.gt.1) then

            if (myrank.eq.0) then
                  dat(:,:,d1:d2)=pix(:,:,d1:d2)
                  do i=1,nprocs-1
                        allocate (pixn(nx,ny,d1s(i):d2s(i)))
                        pixn = pix(:,:,d1s(i):d2s(i))
                        sxip = SIZE(pixn)
                        call MPI_SEND(pixn,2*sxip,MPI_BYTE, i,7,
     $                     MPI_COMM_WORLD,status,ierr)
                        deallocate(pixn)
                  end do
            else
                  allocate(datn(nx,ny,d1:d2))
                  call MPI_RECV(datn,2*sized,MPI_BYTE,0,7,
     $               MPI_COMM_WORLD,status,ierr)
                  dat(:,:,d1:d2) = datn
                  deallocate(datn)
            end if
      end if

      n_time(4)=MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank, " time to get original data= ",
     $   n_time(4)-n_time(3)
      endif

      allocate(vox(nx,ny,d1-1:d2+1))


c  Make the copy

      do k=d1,d2
            vox(:,:,k) = dat(:,:,k)
      end do
      deallocate(dat)

c  Call z_ghost_int to make Z ghost layers of INTEGER*2 values (aka vox).

      call z_ghost_int(vox,nx,ny,d1,d2)
	  
77    format(3(a5,i5,2x))
78    format(a,3(i5,2x))

      if (myrank.eq.0) then
            call dassig(nx,ny,nz,prob,pix)

      do i=1,nphase
            write(7,*) 'Volume fraction of phase ',i,' = ',prob(i)
      end do

      call flush(7)
      deallocate(pix)

      end if
 
   	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (myrank.eq.0) then

c  write out the phase electrical conductivity tensors
      do 11 i=1,nphase
      write(7,*) 'Phase ',i,' conductivity tensor is:'
      write(7,*) sigma(i,1,1),sigma(i,1,2),sigma(i,1,3)
      write(7,*) sigma(i,2,1),sigma(i,2,2),sigma(i,2,3)
      write(7,*) sigma(i,3,1),sigma(i,3,2),sigma(i,3,3)
11 continue

      write(7,*) 'Applied field components:'
      write(7,*) 'ex = ',ex,' ey = ',ey,' ez = ',ez

      call flush(7)
      end if
		
c  Set up the finite element "stiffness" matrices and the Constant and
c  vector required for the energy

      count=0

      n_time(9)=MPI_Wtime(ierr)
      call femat_mpi(nx,ny,nz,d1,d2,vox,sigma,b,dk,C)
      n_time(10)=MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank," femat_mpi time=",n_time(10)-n_time(9)
      endif

      do k=d1,d2
      do j=1,ny
      do i=1,nx
            x=dfloat(i-1)
            y=dfloat(j-1)
            z=dfloat(k-1)
            u(i,j,k)=-x*ex-y*ey-z*ez
      end do; end do; end do
	  

c  Call z_ghost_dp to make Z ghost layers of DOUBLE PRECISION values (aka u).

      call z_ghost_dp(u,nx,ny,d1,d2)
	  
c  RELAXATION LOOP
c  (USER) kmax is the maximum number of times dembx_mpi will be called, with
c  ldemb conjugate gradient steps performed during each call. The total
c  number of conjugate gradient steps allowed for a given elastic
c  computation is kmax*ldemb.
      kmax=100
      ldemb=100
      ltot=0
c  Call energy to get initial energy and initial gradient

      n_time(15)=MPI_Wtime(ierr)

      call energy_mpi(u,dk,b,C,nx,ny,nz,d1,d2,gb,utot,vox)

      n_time(16)=MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank,"Initial energy_mpi time=",n_time(16)-n_time(15)
      endif

c  gg is the norm squared of the gradient (gg=gb*gb)
      dgg= 0.0d0
      gg = 0.0d0
      dgg = SUM(gb(:,:,d1:d2)*gb(:,:,d1:d2))
      call MPI_ALLREDUCE(dgg,gg,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)

      n_time(17)=MPI_Wtime(ierr)

      if (myrank.eq.0) then
            write(*,*) " Initial Energy = ",utot, " gg = ",gg
            write(7,*) " Initial Energy = ",utot, " gg = ",gg
            call flush(7)
      end if

      elapsed_time=0.0d0

      n_time(18)=MPI_Wtime(ierr)

      kkk=0
      kkk=kkk+1
      do kkk=1,kmax
         kkk_start = MPI_Wtime(ierr)

c  call dembx_mpi to go into the conjugate gradient solver

      call dembx_mpi(nx,ny,nz,d1,d2,Lstep,gb,u,vox,h, gg,dk,gtest,ldemb,
     $   kkk)
	 
      ltot=ltot+Lstep
      call energy_mpi(u,dk,b,C,nx,ny,nz,d1,d2,gb,utot,vox)

      kkk_end = MPI_Wtime(ierr)
      elapsed_time=elapsed_time+(kkk_end-kkk_start)

      if (myrank.eq.0) then
            write(7,*) "Energy = ",utot," gg = ",gg
            write(7,*) "Number of conjugate steps = ",ltot
            write(7,*) "Root took ",kkk_end-kkk_start," s for ",
     $         ltot, "conjugate steps."
            write(7,*) "Elapsed time=",elapsed_time," s for ",
     $         ltot, "conjugate steps."
        
            write(*,*) "Energy = ",utot," gg = ",gg
            write(*,*) "Number of conjugate steps = ",ltot
            write(*,*) "Root took ",kkk_end-kkk_start," s for ",
     $         ltot, "conjugate steps."
            write(*,*) "Elapsed time= ",elapsed_time," s for ",
     $         ltot, "conjugate steps."
            call flush(7)
      end if

c  Call energy_mpi to compute energy after dembx_mpi call. If gg < gtest,
c  this will be the final energy. If gg is still larger than gtest,
c  then this will give an intermediate energy with which to check how the
c  relaxation process is coming along.

c  If relaxation process is finished, jump out of loop

      if(gg.le.gtest) goto 444

c  If relaxation process will continue, compute and output stresses
c  and strains as an additional aid to judge how the
c  relaxation procedure is progressing.

      call current_mpi(nx,ny,nz,ns,sigma,vox,u,d1,d2)	

      if (myrank.eq.0) then
c  Output intermediate currents
            write(7,*)
            write(7,*) ' Current in x direction = ',cuxxp
            write(7,*) ' Current in y direction = ',cuyyp
            write(7,*) ' Current in z direction = ',cuzzp
            call flush(7)
      end if
	  
      end do
	  
444   call current_mpi(nx,ny,nz,ns,sigma,vox,u,d1,d2)

      if (myrank.eq.0) then
c  Output final currents
            write(7,*)
            write(7,*) ' Current in x direction = ',cuxxp
            write(7,*) ' Current in y direction = ',cuyyp
            write(7,*) ' Current in z direction = ',cuzzp
            write(7,*)
            call flush(7)
		 
c  Output final currents
            write(*,*)
            write(*,*) ' Current in x direction = ',cuxxp
            write(*,*) ' Current in y direction = ',cuyyp
            write(*,*) ' Current in z direction = ',cuzzp

            close(unit=9)
            close(unit=7)

      end if

      deallocate(vox)

      end do


c  Do another calculation using loop var: npoints


      deallocate(u)
      deallocate(b)
      deallocate(gb)
      deallocate(h)
	  
      n_time(24) = MPI_Wtime(ierr)

      CALL MPI_FINALIZE(ierr)

      end

c  **********************************************************

      subroutine femat_mpi(nx,ny,nz,d1,d2,vox,sigma,b,dk,C)

      implicit none

      include 'mpif.h'

      integer i,ierr,nx,j,ny,nz
      integer d1,d2,myrank,nprocs,mxy
      integer ipx,ipy,ipz
      integer nxy,k,nm,ijk,mm,ii,jj,kk,ll
      integer i8,dn,m,m8
      integer pflag,nphase

      integer status(MPI_STATUS_SIZE)
      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision sum_num,cterm,cpos,cneg
      double precision c,c3,x,y,z
      double precision f_time(24)

      double precision dk(nphase,8,8), sigma(nphase,3,3)
      double precision dndx(8),dndy(8),dndz(8)
      double precision g(3,3,3)
      double precision es(3,8),xn(8)
      double precision b(nx,ny,d1-1:d2+1)
      double precision, allocatable :: ab(:,:), ba(:,:)
      double precision ex,ey,ez

      common/list1/pflag,nphase
      common/list2/ex,ey,ez

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      f_time(1) = MPI_Wtime(ierr)
      nxy=nx*ny
      mxy=3*nxy

      allocate (ab(nx,ny))
      allocate (ba(nx,ny))

c  initialize stiffness matrices

      dk=0.0d0
	  
c  set up Simpson's integration rule weight vector
      do k=1,3
      do j=1,3
      do i=1,3
      nm=0
      if(i.eq.2) nm=nm+1
      if(j.eq.2) nm=nm+1
      if(k.eq.2) nm=nm+1
      g(i,j,k)=4.0d0**nm
      end do
      end do
      end do

c  loop over the nphase kinds of pixels and Simpson's rule quadrature
c  points in order to compute the stiffness matrices. Stiffness matrices
c  of trilinear finite elements are quadratic in x, y, and z, so that
c  Simpson's rule quadrature gives exact results.
      do ijk=1,nphase
      do k=1,3
      do j=1,3
      do i=1,3
      x=dfloat(i-1)/2.0d0
      y=dfloat(j-1)/2.0d0
      z=dfloat(k-1)/2.0d0
c  dndx means the negative derivative, with respect to x, of the shape
c  matrix N (see manual, Sec. 2.2), dndy, and dndz are similar.
      dndx(1)=-(1.0d0-y)*(1.0d0-z)
      dndx(2)=(1.0d0-y)*(1.0d0-z)
      dndx(3)=y*(1.0d0-z)
      dndx(4)=-y*(1.0d0-z)
      dndx(5)=-(1.0d0-y)*z
      dndx(6)=(1.0d0-y)*z
      dndx(7)=y*z
      dndx(8)=-y*z
      dndy(1)=-(1.0d0-x)*(1.0d0-z)
      dndy(2)=-x*(1.0d0-z)
      dndy(3)=x*(1.0d0-z)
      dndy(4)=(1.0d0-x)*(1.0d0-z)
      dndy(5)=-(1.0d0-x)*z
      dndy(6)=-x*z
      dndy(7)=x*z
      dndy(8)=(1.0d0-x)*z
      dndz(1)=-(1.0d0-x)*(1.0d0-y)
      dndz(2)=-x*(1.0d0-y)
      dndz(3)=-x*y
      dndz(4)=-(1.0d0-x)*y
      dndz(5)=(1.0d0-x)*(1.0d0-y)
      dndz(6)=x*(1.0d0-y)
      dndz(7)=x*y
      dndz(8)=(1.0d0-x)*y

c  now build electric field matrix

      es=0.0d0

      es(1,:)=dndx
      es(2,:)=dndy
      es(3,:)=dndz

c  Matrix multiply to determine value at (x,y,z), multiply by
c  proper weight, and sum_num into dk, the stiffness matrix

      f_time(2) = MPI_Wtime(ierr)

      do ii=1,8
      do jj=1,8
c  Define sum over strain matrices and elastic moduli matrix for
c  stiffness matrix
      sum_num=0.0d0
      do kk=1,3
      do ll=1,3
      sum_num=sum_num+es(kk,ii)*sigma(ijk,kk,ll)*es(ll,jj)
      end do; end do
      dk(ijk,ii,jj)=dk(ijk,ii,jj)+g(i,j,k)*sum_num/216.

      end do; end do
      end do; end do; end do; end do

      f_time(3) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank, "time to calculate dk = ",f_time(3)-f_time(2)
      endif

c  Initialize b and C

      b=0.0d0
      C=0.0d0
      c3=0.0d0

999   format(4(i4,1x,),3(f9.6,1x))


c  x=nx face

      i=nx
      do i8=1,8
         xn(i8) = 0.0d0
         if(i8.eq.2.or.i8.eq.3.or.i8.eq.6.or.i8.eq.7) then
            xn(i8)=-ex*nx
         end if
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      dn=d2
      if (dn.eq.nz) then
            dn = nz-1
      end if

      cpos=0.0d0; cneg=0.0d0
      cterm=0.0d0

      do k=d1,dn
      do j=1,ny-1

      m=nxy*(k-1)+j*nx

      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0

      do m8=1,8

      cterm =0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
      cpos = cpos + cterm
      else
      cneg = cneg + cterm
      end if

      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)

      end do

c  Assign b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num

      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num
                                                                                                                                         
      end do
      end do; end do

c  y=ny face

      j=ny
      do i8=1,8
         xn(i8)=0.0d0
         if(i8.eq.3.or.i8.eq.4.or.i8.eq.7.or.i8.eq.8) then
            xn(i8)=-ey*ny
         end if
      end do

      do i=1,nx-1
      do k=d1,dn
      m=nxy*(k-1)+nx*(ny-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m8=1,8

      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)

      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
            cpos = cpos + cterm
      else
         cneg = cneg + cterm
      end if

      end do

      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num

      end do
      end do; end do



c  Zface calcs

c  Only the last node does these series of calculations since
c  it contains all the necessary data therefore no data transfer
c  occurs.

      if (myrank.eq.nprocs-1) then

      k = nz
      do i8=1,8
      xn(i8)=0.0d0
      if(i8.eq.5.or.i8.eq.6.or.i8.eq.7.or.i8.eq.8) then
      xn(i8)=-ez*nz
      end if
      end do

      do i=1,nx-1
      do j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0

      do m8=1,8
      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)
      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
      cpos = cpos + cterm
      else
      cneg = cneg + cterm
      end if

      end do
      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num
      end do
      end do; end do

      end if


c  x=nx y=ny edge


      i=nx
      y=ny

      do i8=1,8
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
      end do

      dn=d2
      if (dn.eq.nz) then
            dn = nz-1
      end if

      do k=d1,dn
      m=nxy*k
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)

      sum_num=0.0d0
      do m8=1,8
      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)

      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
      cpos = cpos + cterm
      else
      cneg = cneg + cterm
      end if

      end do
      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num

      end do
      end do


c  x=nx z=nz edge

      if (myrank.eq.nprocs-1) then

      i=nx
      k=nz

      do i8=1,8
      xn(i8)=0.0d0

      if(i8.eq.2.or.i8.eq.3) then
      xn(i8)=-ex*nx
      end if

      if(i8.eq.5.or.i8.eq.8) then
      xn(i8)=-ez*nz
      end if

      if(i8.eq.6.or.i8.eq.7) then
      xn(i8)=-ez*nz-ex*nx
      end if

      end do

      do j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+nx
      call m2ijk(m,ii,jj,kk,nx,ny,nx)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0

      do m8=1,8
      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)

      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
      cpos = cpos + cterm
      else
      cneg = cneg + cterm
      end if

      end do

      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num

      end do
      end do

c  y=ny z=nz edge

      j=ny
      k=nz

      do i8=1,8
      xn(i8)=0.0d0

      if(i8.eq.5.or.i8.eq.6) then
      xn(i8)=-ez*nz
      end if

      if(i8.eq.3.or.i8.eq.4) then
      xn(i8)=-ey*ny
      end if

      if(i8.eq.7.or.i8.eq.8) then
      xn(i8)=-ey*ny-ez*nz
      end if
      end do

      do i=1,nx-1
      m=nxy*(nz-1)+nx*(ny-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nx)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0

      do m8=1,8
      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)
      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)
      
	  if (cterm.ge.0.0d0) then
            cpos = cpos + cterm
      else
            cneg = cneg + cterm
      end if
 
      end do

      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num
      end do
      end do

c  x=nx y=ny z=nz corner

      i=nx
      j=ny
      k=nz

      do i8=1,8
      xn(i8)=0.0d0

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

      end do

      m=nx*ny*nz
      call m2ijk(m,ii,jj,kk,nx,ny,nx)

      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0

      do m8=1,8
      sum_num=sum_num+xn(m8)*dk(vox(ii,jj,kk),m8,mm)
      cterm=0.5d0*xn(m8)*dk(vox(ii,jj,kk),m8,mm)*xn(mm)

      if (cterm.ge.0.0d0) then
      cpos = cpos + cterm
      else
      cneg = cneg + cterm
      end if

      end do
      b(ipx,ipy,ipz) = b(ipx,ipy,ipz) + sum_num

      end do

c  End if for (myrank.eq.nprocs-1)

      end if

      c3 = cpos + cneg
      CALL MPI_ALLREDUCE(c3,C,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      if (myrank.eq.0) then
            write(*,*) "Final C = ", C
      end if

      f_time(4) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*)myrank,"Etime to calculate C & b= ",f_time(4)-f_time(3)
      end if

      if (nprocs.gt.1) then

c  RECV a new slice per node.

      ab = 0.0d0
      ba = b(:,:,d2+1)

      f_time(5) = MPI_Wtime(ierr)
      call t2b_dp(ab,ba,nx,ny)
      f_time(6) = MPI_Wtime(ierr)
      b(:,:,d1) = b(:,:,d1) + ab

      if (pflag.eq.1) then
      write(*,*) myrank, " B upddate: t2b time= ",f_time(6)-f_time(5)
      end if

c  botp = d1-1

      ab = 0.0
      ba = b(:,:,d1-1)

      f_time(7) = MPI_Wtime(ierr)
      call b2t_dp(ab,ba,nx,ny)
      f_time(8) = MPI_Wtime(ierr)
      b(:,:,d2) = b(:,:,d2) + ab

      if (pflag.eq.1) then
      write(*,*) myrank, " B upddate: b2t time= ",f_time(8)-f_time(7)
      end if

c  Update ghost layers

c  RECV a new slice per node.

      ab = b(:,:,d1)
      ba = b(:,:,d2)

      f_time(9) = MPI_Wtime(ierr)
      call t2b_dp(ab,ba,nx,ny)
      f_time(10) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank, "B ghost upddate:t2b time= ",
     $   f_time(10)-f_time(9)
      end if

      b(:,:,d1-1) = ab

      ab = b(:,:,d1)
      ba = b(:,:,d2)

      f_time(11) = MPI_Wtime(ierr)
      call b2t_dp(ab,ba,nx,ny)
      f_time(12) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank, "B ghost upddate:b2t time= ",
     $   f_time(12)-f_time(11)
      end if

      b(:,:,d2+1) = ba

      else

c  nprocs=1

      b(:,:,d1) = b(:,:,d1) + b(:,:,d2+1)
      b(:,:,d2) = b(:,:,d2) + b(:,:,d1-1)
      b(:,:,d2+1) = b(:,:,d1)
      b(:,:,d1-1) = b(:,:,d2)

      end if

      deallocate(ab)
      deallocate(ba)

      f_time(13) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank, "Femat_mpi elapsed time= ",
     $   f_time(13)-f_time(1)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine energy_mpi(u,dk,b,C,nx,ny,nz,d1,d2,gb,utot,vox)
      implicit none

      include 'mpif.h'

      integer nx,ny,nz,d1,d2,myrank,nprocs,ierr
      integer m3,ik,ij,ii
      integer pflag,nphase

      double precision u(nx,ny,d1-1:d2+1)
      double precision b(nx,ny,d1-1:d2+1)
      double precision gb(nx,ny,d1-1:d2+1)
      integer*2 vox(nx,ny,d1-1:d2+1)
      double precision e_time(24)

      double precision c,utot
      double precision dk(nphase,8,8)

      double precision dutot
      double precision ex,ey,ez

      common/list1/pflag,nphase
      common/list2/ex,ey,ez

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      e_time(1) = MPI_Wtime(ierr)

      gb = 0.0d0


c  After this call, gb is calculated and data slabs
c  are updated and passed.

      call gbah(gb,u,dk,vox,nx,ny,nz,d1,d2)


c  Now do the rest of the gb calculations that appear
c  in original "energy" subroutine.

c  utot will be a per processor value.
c  Do an MPI_ALLREDUCE on dutot
c  so each node will have the current updated version.

      dutot=0.0d0

      do ik=d1,d2
      do ij=1,ny
      do ii=1,nx

      dutot=dutot+0.5d0*u(ii,ij,ik)*gb(ii,ij,ik)+
     $   b(ii,ij,ik)*u(ii,ij,ik)

      end do; end do; end do

      call MPI_ALLREDUCE(dutot,utot,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)

      utot = utot + C

c  easier to add C here than before the above MPI call.

      gb = gb + b

      return
      end

c  **********************************************************

      subroutine dembx_mpi(nx,ny,nz,d1,d2,Lstep,gb,u,vox,h,
     $   gg,dk,gtest,ldemb,kkk)

      implicit none

      include 'mpif.h'

      integer nx,ny,nz,d1,d2,ldemb,kkk,ijk
      integer Lstep,myrank,nprocs,ierr
      integer pflag,nphase

      double precision dgg,gg,gglast,lambda,hAh2,hAh,gamma,gtest

      double precision u(nx,ny,d1-1:d2+1)
      double precision gb(nx,ny,d1-1:d2+1)
      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision dk(nphase,8,8)

      double precision Ah(nx,ny,d1-1:d2+1)
      double precision h(nx,ny,d1-1:d2+1)

      common/list1/pflag,nphase

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if(kkk.eq.1) then
            h=gb
      end if

c  Lstep counts the number of conjugate gradient steps taken in
c  each call to dembx

      Lstep=0

      do ijk=1,ldemb
      Lstep=Lstep+1

      Ah=0.0d0

      call gbah(Ah,h,dk,vox,nx,ny,nz,d1,d2)

      hAh = 0.0d0
      hAh2= 0.0d0

      hAh2 = SUM(h(:,:,d1:d2)*Ah(:,:,d1:d2))

      call MPI_ALLREDUCE(hAh2,hAh,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $   MPI_COMM_WORLD, ierr)

      lambda=gg/hAh
      u=u-lambda*h
      gb=gb-lambda*Ah
      gglast=gg
      gg=0.0d0

      dgg = SUM(gb(:,:,d1:d2)*gb(:,:,d1:d2))
      call MPI_ALLREDUCE(dgg,gg,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)

      if (gg.lt.gtest) goto 1000

      gamma = gg/gglast
      h = gb + gamma*h

      end do

1000  continue

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine current_mpi(nx,ny,nz,ns,sigma,vox,u,d1,d2)

      implicit none
      include 'mpif.h'

      integer nx,ny,nz,ns,d1,d2,nxy
      integer i,j,k,m,n,nn
      integer ifxa,ifya,pflag,nphase

      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision af(3,8)
      double precision u(nx,ny,d1-1:d2+1), uu(8)
      double precision sigma(nphase,3,3)

      double precision cur1,cur2,cur3,ex,ey,ez
      double precision currx,curry,currz
      double precision cuxxp,cuyyp,cuzzp

      integer myrank, ierr, nprocs
      integer status(MPI_STATUS_SIZE)

      common/list1/pflag,nphase
      common/list2/ex,ey,ez
      common/list3/currx,curry,currz
      common/list4/cuxxp,cuyyp,cuzzp

      nxy=nx*ny
	  
c  af is the average field matrix, average field in a pixel is af*u(pixel).
c  The matrix af relates the nodal voltages to the average field in the pixel.

c  Set up single element average field matrix
      af(1,1)=0.25d0
      af(1,2)=-0.25d0
      af(1,3)=-0.25d0
      af(1,4)=0.25d0
      af(1,5)=0.25d0
      af(1,6)=-0.25d0
      af(1,7)=-0.25d0
      af(1,8)=0.25d0
      af(2,1)=0.25d0
      af(2,2)=0.25d0
      af(2,3)=-0.25d0
      af(2,4)=-0.25d0
      af(2,5)=0.25d0
      af(2,6)=0.25d0
      af(2,7)=-0.25d0
      af(2,8)=-0.25d0
      af(3,1)=0.25d0
      af(3,2)=0.25d0
      af(3,3)=0.25d0
      af(3,4)=0.25d0
      af(3,5)=-0.25d0
      af(3,6)=-0.25d0
      af(3,7)=-0.25d0
      af(3,8)=-0.25d0

c  now compute current in each pixel
      currx=0.0d0
      curry=0.0d0
      currz=0.0d0

c  compute average field in each pixel

      do 470 k=d1,d2
      do 470 j=1,ny
      do 470 i=1,nx
      m=(k-1)*nxy+(j-1)*nx+i

      if ((i+1).GT.nx) then
            ifxa = 1
      else
            ifxa = i+1
      end if

      if ((j+1).GT.ny) then
            ifya = 1
      else
            ifya = j+1
      end if
	  
c  load in elements of 8-vector using pd. bd. conds.

      uu(1)= u(i,j,k)
      uu(2)= u(ifxa,j,k)
      uu(3)= u(ifxa,ifya,k)
      uu(4)= u(i,ifya,k)
      uu(5)= u(i,j,k+1)
      uu(6)= u(ifxa,j,k+1)
      uu(7)= u(ifxa,ifya,k+1)
      uu(8)= u(i,ifya,k+1)

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
      cur1=0.0d0
      cur2=0.0d0
      cur3=0.0d0

      do 465 n=1,8
      do 465 nn=1,3

      cur1=cur1+sigma(vox(i,j,k),1,nn)*af(nn,n)*uu(n)
      cur2=cur2+sigma(vox(i,j,k),2,nn)*af(nn,n)*uu(n)
      cur3=cur3+sigma(vox(i,j,k),3,nn)*af(nn,n)*uu(n)

465   continue

c  sum into the global average currents

      currx=currx+cur1
      curry=curry+cur2
      currz=currz+cur3

470   continue

      call MPI_ALLREDUCE(currx,cuxxp,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(curry,cuyyp,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(currz,cuzzp,1,MPI_DOUBLE_PRECISION,
     $   MPI_SUM,MPI_COMM_WORLD,ierr)


c  Volume average currents
      cuxxp=cuxxp/dfloat(ns)
      cuyyp=cuyyp/dfloat(ns)
      cuzzp=cuzzp/dfloat(ns)

      return
      end

c  **********************************************************

      subroutine gbah(om,uh,dk,vox,nx,ny,nz,d1,d2)

      implicit none
      include 'mpif.h'

      integer nx,ny,nz,d1,d2,mxy,pflag,nphase
      integer im,jm,km,ifxa,ifxb,ifya,ifyb
      integer myrank,nprocs,ierr

      double precision uh(nx,ny,d1-1:d2+1)
      double precision om(nx,ny,d1-1:d2+1)
      double precision gb_time(6)

      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision dk(nphase,8,8)

      common/list1/pflag,nphase

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      gb_time(1) = MPI_Wtime(ierr)

      om = 0.0d0

      do km=d1,d2
      do jm=1,ny
      do im=1,nx

      if ((im+1).GT.nx) then
            ifxa = 1
      else
            ifxa = im+1
      end if

      if ((im-1).LE.0) then
            ifxb = nx
      else
            ifxb = im-1
      end if

      if ((jm+1).GT.ny) then
            ifya = 1
      else
            ifya = jm+1
      end if

      if ((jm-1).LE.0) then
            ifyb = ny
      else
            ifyb = jm-1
      end if

c  SELF TERM

      om(im,jm,km) =

c  u(ib(m,1))
     & uh(im,ifya,km)*
     &(dk(vox(im,jm,km),1,4)
     &+dk(vox(ifxb,jm,km),2,3)
     &+dk(vox(im,jm,km-1),5,8)
     &+dk(vox(ifxb,jm,km-1),6,7) )+

c  u(ib(m,2))
     & uh(ifxa,ifya,km)*
     & (dk(vox(im,jm,km),1,3)+dk(vox(im,jm,km-1),5,7) )+
	 
c  u(ib(m,3))
     & uh(ifxa,jm,km)*(dk(vox(im,jm,km),1,2)
     &+ dk(vox(im,ifyb,km),4,3)
     &+ dk(vox(im,ifyb,km-1),8,7)
     &+ dk(vox(im,jm,km-1),5,6) ) +

c  u(ib(m,4))
     & uh(ifxa,ifyb,km)*(dk(vox(im,ifyb,km),4,2)
     &+ dk(vox(im,ifyb,km-1),8,6) ) +

c  u(ib(m,5))
     & uh(im,ifyb,km)*(dk(vox(ifxb,ifyb,km),3,2)
     & +dk(vox(im,ifyb,km),4,1)
     & +dk(vox(ifxb,ifyb,km-1),7,6)
     & +dk(vox(im,ifyb,km-1),8,5) ) +

c  u(ib(m,6))
     & uh(ifxb,ifyb,km)*(dk(vox(ifxb,ifyb,km),3,1)
     &+ dk(vox(ifxb,ifyb,km-1),7,5) ) +

c  u(ib(m,7))
     & uh(ifxb,jm,km)*(
     & dk(vox(ifxb,ifyb,km),3,4)
     &+dk(vox(ifxb,jm,km),2,1)
     &+dk(vox(ifxb,ifyb,km-1),7,8)
     &+dk(vox(ifxb,jm,km-1),6,5) ) +

c  u(ib(m,8))
     & uh(ifxb,ifya,km)*( dk(vox(ifxb,jm,km),2,4)
     &+dk(vox(ifxb,jm,km-1),6,8) ) +

c  u(ib(m,9))
     & uh(im,ifya,km-1)*(dk(vox(im,jm,km-1),5,4)
     &+ dk(vox(ifxb,jm,km-1),6,3) ) +

c  u(ib(m,10))
     & uh(ifxa,ifya,km-1)*(dk(vox(im,jm,km-1),5,3) )+

c  u(ib(m,11))
     & uh(ifxa,jm,km-1)*(dk(vox(im,ifyb,km-1),8,3)
     &+ dk(vox(im,jm,km-1),5,2) ) +

c  u(ib(m,12))
     & uh(ifxa,ifyb,km-1)*( dk(vox(im,ifyb,km-1),8,2) )+

c  u(ib(m,13))
     & uh(im,ifyb,km-1)*(dk(vox(im,ifyb,km-1),8,1)
     &+ dk(vox(ifxb,ifyb,km-1),7,2) ) +

c  u(ib(m,14))
     & uh(ifxb,ifyb,km-1)*( dk(vox(ifxb,ifyb,km-1),7,1) )+

c  u(ib(m,15))
     & uh(ifxb,jm,km-1)*(dk(vox(ifxb,ifyb,km-1),7,4)
     &+ dk(vox(ifxb,jm,km-1),6,1) )+

c  u(ib(m,16))
     &uh(ifxb,ifya,km-1)*( dk(vox(ifxb,jm,km-1),6,4) )+

c  u(ib(m,17))
     & uh(im,ifya,km+1)*(dk(vox(im,jm,km),1,8)
     &+ dk(vox(ifxb,jm,km),2,7) )+

c  u(ib(m,18))
     & uh(ifxa,ifya,km+1)*( dk(vox(im,jm,km),1,7) )+

c  u(ib(m,19))
     & uh(ifxa,jm,km+1)*(dk(vox(im,jm,km),1,6)
     &+ dk(vox(im,ifyb,km),4,7) ) +

c  u(ib(m,20))
     & uh(ifxa,ifyb,km+1)*( dk(vox(im,ifyb,km),4,6) )+

c  u(ib(m,21))
     & uh(im,ifyb,km+1)*(dk(vox(im,ifyb,km),4,5)
     &+ dk(vox(ifxb,ifyb,km),3,6) ) +

c  u(ib(m,22))
     & uh(ifxb,ifyb,km+1)*( dk(vox(ifxb,ifyb,km),3,5) )+

c  u(ib(m,23))
     & uh(ifxb,jm,km+1)*(dk(vox(ifxb,ifyb,km),3,8)
     &+ dk(vox(ifxb,jm,km),2,5) ) +

c  u(ib(m,24))
     & uh(ifxb,ifya,km+1)*( dk(vox(ifxb,jm,km),2,8) )+

c  u(ib(m,25))
     & uh(im,jm,km-1)*(dk(vox(ifxb,ifyb,km-1),7,3)
     &+ dk(vox(im,ifyb,km-1),8,4)
     &+ dk(vox(ifxb,jm,km-1),6,2)
     &+ dk(vox(im,jm,km-1),5,1) ) +

c  u(ib(m,26))
     & uh(im,jm,km+1)*(
     & dk(vox(ifxb,ifyb,km),3,7)
     &+dk(vox(im,ifyb,km),4,8)
     &+dk(vox(im,jm,km),1,5)
     &+dk(vox(ifxb,jm,km),2,6) ) +

c  u(ib(m,27))
     & uh(im,jm,km)* (dk(vox(im,jm,km),1,1)
     &+ dk(vox(ifxb,jm,km),2,2)
     &+ dk(vox(ifxb,ifyb,km),3,3)
     &+ dk(vox(im,ifyb,km),4,4)
     &+ dk(vox(im,jm,km-1),5,5)
     &+ dk(vox(ifxb,jm,km-1),6,6)
     &+ dk(vox(ifxb,ifyb,km-1),7,7)
     &+ dk(vox(im,ifyb,km-1),8,8) )

      end do; end do; end do

      gb_time(2) = MPI_Wtime(ierr)

c  Do top/bottom layer switch on matrix: om

      call z_ghost_dp(om,nx,ny,d1,d2)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine dpixel(nx,ny,nz,ns,pix)
      implicit none

      integer nx,ny,nz,ns,nphase,nxy
      integer i,j,k,m,pflag
      integer*2 pix(nx,ny,nz)
      integer*2 pix0

      common/list1/pflag,nphase

c  (USER) If you want to set up a test image inside the program, instead of
c  reading it in from a file, this should be done inside this subroutine.

      nxy=nx*ny
      do 200 k=1,nz
      do 200 j=1,ny
      do 200 i=1,nx
      m=nxy*(k-1)+nx*(j-1)+i
      read(9,*) pix(i,j,k)
200   continue

      do k=1,nz
      do j=1,ny
      do i=1,nx

      pix0 = pix(i,j,k)

      if(pix0.lt.1) then
            write(7,*) "Phase label in pix < 1--error at ",i,j,k
      end if
      if(pix0.gt.nphase) then
            write(7,*) "Phase label in pix > nphase--error at ",i,j,k
      end if

      end do; end do; end do


      return
      end

c  **********************************************************

      subroutine dassig(nx,ny,nz,prob,pix)
      implicit none

      integer nx,ny,nz,ns,nphase,ii,jj,kk,i,pflag

      integer*2 pix(nx,ny,nz)
      double precision prob(nphase)

      common/list1/pflag,nphase

      ns=nx*ny*nz
      prob=0.0d0

      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
      do i=1,nphase
            if(pix(ii,jj,kk).eq.i) then
                  prob(i)=prob(i)+1.0d0
            end if
      end do; end do
      end do; end do

      prob=prob/dfloat(ns)

      return
      end

c  **********************************************************

      subroutine ipxyz(mm,i,j,k,ipx,ipy,ipz,nx,ny,nz)

      implicit none
      integer mm,i,j,k,ipx,ipy,ipz,nx,ny,nz

      if (mm.le.4) then
            ipz=k
      else
            ipz=k+1
      end if

      if ((mm.eq.1).OR.(mm.eq.5)) then
            ipx=i
            ipy=j
      end if

      if ((mm.eq.2).OR.(mm.eq.6)) then
            ipx = i+1
            ipy=j

            if (i.ge.nx) then
                  ipx=1
            end if

      end if

      if ((mm.eq.3).OR.(mm.eq.7)) then
            ipx = i+1
            if (i.ge.nx) then
                  ipx=1
            end if
            ipy = j+1
            if (j.ge.ny) then
                  ipy=1
            end if
      end if

      if ((mm.eq.4).OR.(mm.eq.8)) then
            ipx = i
            ipy = j+1

            if (j.ge.ny) then
            ipy=1
         end if

      end if

      return
      end

c  **********************************************************

      subroutine m2ijk(inps,i,j,k,ni,nj,nk)

      implicit none
      integer inps,ns
      integer c
      integer kdiv,jdiv
      integer rj,rk
      integer i,j,k,ni,nj,nk

      ns=ni*nj
      kdiv=inps/ns
      c = ns*kdiv
      rk = inps-c

      if (rk.eq.0) then
            k=kdiv
            j=nj
            i=ni
      else
            k=kdiv+1
      end if

      if (k.ne.kdiv) then
	  
  	        jdiv=rk/ni
            c=jdiv*ni
            rj = rk-c

            if (rj.eq.0) then
                  j=jdiv
                  i=ni
            else
                  j=jdiv+1
                  i=rj
            end if

      end if

      return
      end

c  **********************************************************

      subroutine z_ghost_int(arr0,mx,my,d1,d2)

      implicit none

      include 'mpif.h'

      integer mx,my,mz,d1,d2

      integer*2 arr0(mx,my,d1-1:d2+1)
      integer*2, allocatable :: bot(:,:), top(:,:)

      integer myrank, ierr, nprocs
      integer status(MPI_STATUS_SIZE)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      allocate(bot(mx,my))
      allocate(top(mx,my))

c  Get new bottom ghost plane.

      bot = arr0(:,:,d1)
      top = arr0(:,:,d2)

      call t2b(bot,top,mx,my)

      arr0(:,:,d1-1) = bot

c  Get new top ghost plane

      bot = arr0(:,:,d1)
      top = arr0(:,:,d2)

      call b2t(bot,top,mx,my)

      arr0(:,:,d2+1) = top

      deallocate(bot)
      deallocate(top)

      return
      end

c  **********************************************************

      subroutine z_ghost_dp(arr0,mx,my,d1,d2)

      implicit none

      include 'mpif.h'

      integer mx,my,d1,d2

      double precision arr0(mx,my,d1-1:d2+1)

      double precision, allocatable :: bot(:,:), top(:,:)

      integer myrank, ierr, nprocs
      integer status(MPI_STATUS_SIZE)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      allocate(bot(mx,my))
      allocate(top(mx,my))

c  Get new bottom ghost plane.

      bot = arr0(:,:,d1)
      top = arr0(:,:,d2)
      call t2b_dp(bot,top,mx,my)
      arr0(:,:,d1-1) = bot

c  Get new top ghost plane

      bot = arr0(:,:,d1)
      top = arr0(:,:,d2)
      call b2t_dp(bot,top,mx,my)
      arr0(:,:,d2+1) = top
      deallocate(bot)
      deallocate(top)

      return
      end

c  **********************************************************

      subroutine t2b(b_layer,t_layer,nx,ny)

c  This is an INTEGER*2 subroutine.

c  Used for transferring: pix bottom2top layers

c  RECV a new t_layer (TOP layer) per node.

      implicit none

      include 'mpif.h'

      integer nx,ny,nxy
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)
      
      integer*2 b_layer(nx,ny), t_layer(nx,ny)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      nxy=nx*ny

      ides = mod(myrank+1,nprocs)
      isrc = mod(myrank+nprocs-1,nprocs)


      if (myrank.eq.nprocs-1) then
      call MPI_Irecv(b_layer,2*nxy, MPI_BYTE, isrc,
     &   9,MPI_COMM_WORLD, irequest, ierr)
      call mpi_send(t_layer,2*nxy,MPI_BYTE,ides,9,MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)

      else

      call mpi_recv(b_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,
     &   status,ierr)
      call mpi_send(t_layer,2*nxy,MPI_BYTE,ides,9,MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine b2t(b_layer,t_layer,nx,ny)


c  This is an INTEGER*2 subroutine.

c  Used for transferring: pix bottom2top layers

c  RECV a new t_layer (TOP layer) per node.

      implicit none

      include 'mpif.h'

      integer nx,ny,nxy
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)

      integer*2 b_layer(nx,ny), t_layer(nx,ny)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      nxy=nx*ny

      ides = mod(myrank+nprocs-1,nprocs)
      isrc = mod(myrank+1,nprocs)

      if (myrank.eq.nprocs-1) then
      call MPI_Irecv(t_layer,2*nxy, MPI_BYTE, isrc,
     &   9,MPI_COMM_WORLD, irequest, ierr)
      call mpi_send(b_layer,2*nxy,MPI_BYTE,ides,9,
     &   MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)

      else

      call mpi_recv(t_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,
     &   status,ierr)
      call mpi_send(b_layer,2*nxy,MPI_BYTE,ides,9,
     &   MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine t2b_dp(b_layer,t_layer,nx,ny)

c  This is a double precision subroutine.

c  Used for transferring: u,b,and om top2bottom layers

c  RECV a new b_layer (BOTTOM layer) per node.

      implicit none

      include 'mpif.h'

      integer nx,ny,mxy
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)
      double precision b_layer(nx,ny), t_layer(nx,ny)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      mxy=nx*ny

      ides = mod(myrank+1,nprocs)
      isrc = mod(myrank+nprocs-1,nprocs)

      if (myrank.eq.nprocs-1) then
      call mpi_irecv(b_layer,mxy,MPI_DOUBLE_PRECISION,isrc,9,
     &   MPI_COMM_WORLD,irequest,ierr)
      call mpi_send(t_layer,mxy,MPI_DOUBLE_PRECISION,ides,9,
     &   MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)

      else

      call mpi_recv(b_layer,mxy,MPI_DOUBLE_PRECISION,isrc,9,
     &   MPI_COMM_WORLD,status,ierr)
      call mpi_send(t_layer,mxy,MPI_DOUBLE_PRECISION,ides,9,
     &   MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine b2t_dp(b_layer,t_layer,nx,ny)


c  This is a double precision subroutine.

c  Used for transferring: u,b,and om bottom2top layers

c  RECV a new t_layer (TOP layer) per node.

      implicit none

      include 'mpif.h'

      integer nx,ny,mxy
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)

      double precision b_layer(nx,ny), t_layer(nx,ny)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      mxy=nx*ny

      ides = mod(myrank+nprocs-1,nprocs)
      isrc = mod(myrank+1,nprocs)

      if (myrank.eq.nprocs-1) then
      call mpi_Irecv(t_layer,mxy,MPI_DOUBLE_PRECISION,isrc,9,
     &   MPI_COMM_WORLD,irequest,ierr)
      call mpi_send(b_layer,mxy,MPI_DOUBLE_PRECISION,ides,9,
     &   MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)

      else

      call mpi_recv(t_layer,mxy,MPI_DOUBLE_PRECISION,isrc,9,
     &   MPI_COMM_WORLD,status,ierr)
      call mpi_send(b_layer,mxy,MPI_DOUBLE_PRECISION,ides,9,
     &   MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end