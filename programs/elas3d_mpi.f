c  *********************** elas3d_mpi.f **************************

c  This is the new MPI version of the elas3d.f code from
c  Section 9.3.2 of NISTIR 6269.

c  The main differences with this code compared to the serial
c  version are:

c  1. Removal of ib array.
c  2. Change of dimensionality on pix from pix(m) to pix(i,j,k)
c  Maximum value of m = nx*ny*nz (nx,ny,nz are the array dims).
c  3. All important arrays (pix,vox,gb,b,u,h,Ah) are dynamically allocated.

c  IN THIS VERSION:

c  The USER needs the following input:
c  (Search for occurences of "USER" in the code).

c  1. A 3-D pixel value data file with input & output names.
c  2. The values of the 3 dimensions: (nx,ny,nz)
c  3. The number of phases in the mixture: nphase
c  4. A convergence value: gtest
c  5. Initial values for shears and strains: exx,eyy,ezz,exy,exz,eyz
c  6. Values for DEMBX_MPI and how long it will run: kmax & ldemb

c  7. Flag for printing timing info for all data
c       passing MPI routines ( FEMAT_MPI, ENERGY_MPI, DEMBX)
c       from MAIN is called: pflag
c       pflag Values = 0,1 0=no timing i    nfo; 1=print timing info

c       pflag is a common value.

c       Timing info for the RELAXATION loop is not
c       influenced by the pflag and will always be printed.

c       User may edit the code to supress the printing.

c  8. Timing info stored in arrays namex X_time(i)
c       Where X=n,f,e ie.
c       n_time is in MAIN
c       f_time is in FEMAT_MPI
c       e_time is in ENERGY_MPI

c  NB: One also needs to insure that the values for
c        phasemod(i,j) are initialized correctly in
c        SUBROUTINE phasemod_init.

c  END of NEW comments.

c  BEGIN ORIGINAL comments.

c  BACKGROUND


c  This program solves the linear elastic equations in a
c  random linear elastic material, subject to an applied macroscopic strain,
c  using the finite element method. Each pixel in the 3-D digital
c  image is a cubic tri-linear finite element, having its own
c  elastic moduli tensor. Periodic boundary conditions are maintained.
c  In the comments below, (USER) means that this is a section of code that
c  the user might have to change for his particular problem. Therefore the
c  user is encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the energy
c  1/2 uAu + b u + C, where A is the Hessian matrix composed of the
c  stiffness matrices (dk) for each pixel/element, b is a constant vector
c  and C is a constant that are determined by the applied strain and
c  the periodic boundary conditions, and u is a vector of
c  all the displacements. The solution
c  method used is the conjugate gradient relaxation algorithm.
c  Other variables are: gb is the gradient = Au+b, h and Ah are
c  auxiliary variables used in the conjugate gradient algorithm (in dembx),
c  dk(n,i,j) is the stiffness matrix of the n’th phase, cmod(n,i,j) is
c  the elastic moduli tensor of the n’th phase, pix is a vector that gives
c  the phase label of each pixel, ib is a matrix that gives the labels of
c  the 27 (counting itself) neighbors of a given node, prob is the volume
c  fractions of the various phases,
c  strxx, stryy, strzz, strxz, stryz, and strxy are the six Voigt
c  volume averaged total stresses, and
c  sxx, syy, szz, sxz, syz, and sxy are the six Voigt
c  volume averaged total strains.

c  DIMENSIONS

c  The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
c  ns=nx*ny*nz, with three components, where the digital image of the
c  microstructure considered is a rectangular paralleliped, nx x ny x nz
c  in size. The arrays pix and ib are are also dimensioned to the system size.
c  The array ib has 27 components, for the 27 neighbors of a node.
c  Note that the program is set up at present to have at most 100
c  different phases. This can easily be changed, simply by changing
c  the dimensions of dk, prob, and cmod. The parameter nphase gives the
c  number of phases being considered in the problem.
c  All arrays are passed between subroutines using simple common statements.

c  STRONGLY SUGGESTED: READ THE MANUAL BEFORE USING PROGRAM!!
      implicit none
      include 'mpif.h'

c  (USER) Change the nx,ny,nz dimensions at the beginning.
c  All important arrays are dynamically allocated.

      integer*2,allocatable :: dat(:,:,:),datn(:,:,:)
      integer*2,allocatable :: pix(:,:,:),pixn(:,:,:)
      integer*2,allocatable :: vox(:,:,:)

      integer,allocatable :: d1s(:),d2s(:)

      double precision,allocatable :: b(:,:,:,:)
      double precision,allocatable :: gb(:,:,:,:)
      double precision,allocatable :: u(:,:,:,:)
      double precision,allocatable :: h(:,:,:,:)

      double precision,allocatable :: phasemod(:,:),prob(:)
      double precision,allocatable :: dk(:,:,:,:,:),cmod(:,:,:)
      
      double precision dgg,gg,utot,gtest,C
      double precision exx,eyy,ezz,exz,eyz,exy
      double precision x,y,z,saves
      double precision strxxp,stryyp,strzzp,strxyp,strxzp,stryzp
      double precision sxxp,syyp,szzp,sxyp,sxzp,syzp

      double precision bulk,shear,young,pois

      integer d1,d2,ns,sxip,kkk
      integer i,j,k,nx,ny,nz,nxy,nphase
      integer rem,sz,sized
      integer mxy

      integer npoints,micro
      integer kmax,ldemb,ltot,lstep
      integer pflag

      integer irank
      integer myrank,ierr,nprocs
      integer status(MPI_STATUS_SIZE)

      double precision starttime,endtime,start_npoint,end_npoint
      double precision kkk_start,kkk_end
      double precision elapsed_time,stress_loop
      double precision n_time(24)

      common/list1/pflag,nphase
      common/list2/exx,eyy,ezz,exz,eyz,exy
      common/list3/strxxp,stryyp,strzzp,strxyp,strxzp,stryzp
      common/list4/sxxp,syyp,szzp,sxyp,sxzp,syzp

      call MPI_INIT(ierr)

      starttime = MPI_Wtime(ierr)

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      if (myrank.eq.0) then
      write(*,*) "There are ",nprocs," processors running this job."
      end if

c  USER: Change nx,ny,nz values to match your data.

      nx=512
      ny=512
      nz=512
      nphase=2

      nxy=nx*ny
      ns=nx*ny*nz
      sz=nz/nprocs
      mxy=3*nx*ny 

      gtest=1.e-6*ns

      ! pflag=0 for no timing info printed.
      ! pflag=1 for timing info printed.
      pflag = 0

c  End this USER section.

      utot =0.0d0

c  USER: put phasemod definitions in
c        subroutine "phasemod_init".

      allocate(phasemod(nphase,2))
      call phasemod_init(phasemod)

      allocate( dk(nphase,8,3,8,3) )

      allocate( cmod(nphase,6,6) )
      allocate( prob(nphase) )
      
      if (myrank.eq.0) then
      
      allocate (d1s(0:nprocs-1))
      allocate (d2s(0:nprocs-1))
      
      do irank=0,nprocs-1
            d1s(irank)=irank*sz+1
            d2s(irank)=(irank+1)*sz
      end do
      
      rem = nz-nprocs*sz
      
      if (rem.ne.0) then
            do j=1,rem
                  irank=nprocs-rem+j-1
                  d1s(irank)=d1s(irank)+j-1
                  d2s(irank)=d2s(irank)+j
            end do
      end if
      
c  Send all d1s(i) and d2s(i) from ROOT
c  to NODE i & store into d1 & d2
      
      do i=0,nprocs-1
      call MPI_SEND(d1s(i),1,MPI_integer,i,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(d2s(i),1,MPI_integer,i,1,MPI_COMM_WORLD,ierr)
      end do
      
      end if
      
      call MPI_RECV(d1,1,MPI_integer,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(d2,1,MPI_integer,0,1,MPI_COMM_WORLD,status,ierr)
      write(*,*) "Rank#",myrank,"d1= ",d1," d2= ",d2
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
c  Allocate other arrays which need d1&d2 values.
      
      allocate (gb(nx,ny,d1-1:d2+1,3))
      gb=0.0d0
      allocate(b(nx,ny,d1-1:d2+1,3))
      b = 0.0d0
      
      allocate (u(nx,ny,d1-1:d2+1,3))
      allocate (h(nx,ny,d1-1:d2+1,3))

c  Want the ability to calculate on a series
c  of input files based on a value & some if statements.

c  Compute the average stress and strain in each microstructure.
c  (USER) npoints is the number of microstructures to use.

      npoints=1
      n_time(1) = MPI_Wtime(ierr)

      do micro=1,npoints

c  Allocate pix, so root can read it.

      if (myrank.eq.0) then
            allocate (pix(nx,ny,nz))
      end if

      start_npoint=MPI_Wtime(ierr)
      n_time(2) = MPI_Wtime(ierr)
      
      if (myrank.eq.0) then

c  (USER) Unit 9 is the microstructure input file,
c         Unit 7 is the results output file.
c         Get pix from the input file (unit=9).

            if (micro.eq.1) then
            open (unit=9,file='512_6.txt')
            open (unit=7,file='512_6_out.out')
            end if

            write(*,*) "MICRO = ",micro
            write(7,*) "MICRO = ",micro

c Finally... read in pix

      write(*,*) "call dpixel"
      call dpixel(nx,ny,nz,ns,pix)
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
                        call MPI_SEND(pixn,2*sxip,MPI_BYTE,i,7,
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
      write(*,*) myrank," time to get original data= ",
     $   n_time(4)-n_time(3)
      end if
      
      allocate(vox(nx,ny,d1-1:d2+1))
      vox = 0
      
c  Make the copy

      do k=d1,d2
            vox(:,:,k) = dat(:,:,k)
      end do
      deallocate(dat)

c  Call z_ghost_int to make Z ghost layers of INTEGER*2 values (aka vox).

      call z_ghost_int(vox,nx,ny,nz,d1,d2)

77    format(3(a5,i5,2x))
78    format(a,3(i5,2x))

c  Apply chosen strains as a homogeneous macroscopic strain
c  as the initial condition.

      if (myrank.eq.0) then
            call dassig(nx,ny,nz,prob,pix)
      
      do i=1,nphase
            write(7,9020) i,phasemod(i,1),phasemod(i,2)
9020  format("Phase ",i3," bulk = ",f12.6," shear = ",f12.6)
      end do

      do i=1,nphase
            write(7,9065) i,prob(i)
9065  format("Volume fraction of phase ",i3," is ",f8.5)
      end do
      
      call flush(7)
      deallocate(pix)
      
      end if
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c  (USER) Set applied strains.
c  Actual shear strain applied is exy, exz, and eyz as
c  given in the statements below. The engineering shear strain, by which
c  the shear modulus is usually defined, is twice these values.

      exx=0.1d0
      eyy=0.1d0
      ezz=0.1d0
      exz=0.1d0/2.d0
      eyz=0.1d0/2.d0
      exy=0.1d0/2.d0

      if (myrank.eq.0) then
      write(7,*) "Applied engineering strains"
      write(7,*) " exx eyy ezz exz eyz exy"
      write(7,*) exx,eyy,ezz,2.*exz,2.*eyz,2.*exy
      write(*,*) "Applied engineering strains"
      write(*,*) " exx eyy ezz exz eyz exy"
      write(*,*) exx," ",eyy," ",ezz," ",2.*exz," ",2.*eyz," ",2.*exy
      call flush(7)
      end if

c  Set up the elastic modulus variables, finite element stiffness matrices,
c  the constant, C, and vector, b, required for computing the energy.
c  (USER) If anisotropic elastic moduli tensors are used, these need to be
c  input in subroutine femat.

      n_time(9)=MPI_Wtime(ierr)
      call femat_mpi(nx,ny,nz,phasemod,d1,d2,vox,b,dk,C,cmod)
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
            u(i,j,k,1)=x*exx+y*exy+z*exz
            u(i,j,k,2)=x*exy+y*eyy+z*eyz
            u(i,j,k,3)=x*exz+y*eyz+z*ezz
      end do; end do; end do
      
      call z_ghost_dp(u,nx,ny,3,d1,d2)
           
c  RELAXATION LOOP
c  (USER) kmax is the maximum number of times dembx will be called, with
c  ldemb conjugate gradient steps performed during each call. The total
c  number of conjugate gradient steps allowed for a given elastic
c  computation is kmax*ldemb.
      kmax=40
      ldemb=50
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
      dgg = SUM(gb(:,:,d1:d2,:)*gb(:,:,d1:d2,:))
      call MPI_ALLREDUCE(dgg,gg,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      n_time(17)=MPI_Wtime(ierr)

      if (myrank.eq.0) then
            write(*,*) " Initial Energy = ",utot," gg = ",gg
            write(7,*) " Initial Energy = ",utot," gg = ",gg
            call flush(7)
      end if

      elapsed_time=0.0d0
      stress_loop=0.0d0

      n_time(18)=MPI_Wtime(ierr)
      do kkk=1,kmax
            kkk_start = MPI_Wtime(ierr)

c  call dembx_mpi to go into the conjugate gradient solver

      call dembx_mpi(nx,ny,nz,d1,d2,Lstep,gb,u,vox,h,gg,dk,gtest,ldemb,
     $   kkk)
      ltot=ltot+Lstep
      call energy_mpi(u,dk,b,C,nx,ny,nz,d1,d2,gb,utot,vox)

      kkk_end = MPI_Wtime(ierr)
      elapsed_time=elapsed_time+(kkk_end-kkk_start)

      if (myrank.eq.0) then
            write(7,*) "Energy = ",utot," gg = ",gg
            write(7,*) "Number of conjugate steps = ",ltot
            write(7,*) "Root took ",kkk_end-kkk_start," s for ",ltot,
     $         "conjugate steps."
            write(7,*) "Elapsed time=",elapsed_time," s for ",ltot,
     $         "conjugate steps."
     
            write(*,*) "Energy = ",utot," gg = ",gg
            write(*,*) "Number of conjugate steps = ",ltot
            write(*,*) "Root took ",kkk_end-kkk_start," s for ",ltot,
     $         "conjugate steps."
            write(*,*) "Elapsed time= ",elapsed_time," s for ",ltot,
     $         "conjugate steps."
            call flush(7)
      end if

c  Call energy_mpi to compute energy after dembx_mpi call. If gg < gtest,
c  this will be the final energy. If gg is still larger than gtest,
c  then this will give an intermediate energy with which to check how the
c  relaxation process is coming along.
     
c If relaxation process is finished, jump out of loop
      if(gg.le.gtest) goto 444
     
c  If relaxation process will continue, compute and output stresses
c  and strains as an additional aid to judge how the
c  relaxation procedure is progressing.
      
      n_time(19)=MPI_Wtime(ierr)
      call stress_mpi(nx,ny,nz,ns,u,vox,cmod,d1,d2)
      n_time(20)=MPI_Wtime(ierr)
     
      if (myrank.eq.0) then
            write(7,*) " stresses: xx,yy,zz,xz,yz,xy"
            write(7,*) strxxp,stryyp,strzzp,strxzp,stryzp,strxyp
            write(7,*) " strains: xx,yy,zz,xz,yz,xy"
            write(7,*) sxxp,syyp,szzp,sxzp,syzp,sxyp
            call flush(7)
      end if
      
      end do
      
      n_time(21)=MPI_Wtime(ierr)
444   call stress_mpi(nx,ny,nz,ns,u,vox,cmod,d1,d2)
      n_time(22)=MPI_Wtime(ierr)
     
      if (myrank.eq.0) then
            write(7,*) " stresses: xx,yy,zz,xz,yz,xy"
            write(7,*) strxxp,stryyp,strzzp,strxzp,stryzp,strxyp
            write(7,*) " strains: xx,yy,zz,xz,yz,xy"
            write(7,*) sxxp,syyp,szzp,sxzp,syzp,sxyp
            write(*,*) "Energy = ",utot," gg = ",gg
            write(*,*) "Number of conjugate steps = ",ltot
      end if

      bulk=(strxxp+stryyp+strzzp)/(sxxp+syyp+szzp)/3.0d0
      shear=(strxyp/sxyp+strxzp/sxzp+stryzp/syzp)/3.0d0
      young=9.d0*bulk*shear/(3.d0*bulk+shear)
      pois=(3.d0*bulk-2.d0*shear)/2.d0/(3.d0*bulk+shear)

      if (myrank.eq.0) then
            write(7,*) " bulk modulus = ",bulk
            write(7,*) " shear modulus = ",shear
            write(7,*) " Youngs modulus = ",young
            write(7,*) " Poissons ratio = ",pois

            write(*,*) " bulk modulus = ",bulk
            write(*,*) " shear modulus = ",shear
            write(*,*) " Youngs modulus = ",young
            write(*,*) " Poissons ratio = ",pois
      close(unit=9)
      close(unit=7)

      end if

c  Do another using loop var: npoints

      n_time(23) = MPI_Wtime(ierr)
      
      write(*,*) myrank," took ",n_time(23)-n_time(2),
     $   "s for npoints file ",micro
      
      deallocate(vox)
      
      end do

      n_time(24) = MPI_Wtime(ierr)
      write(*,*) myrank," took ",n_time(24)-n_time(1),"for all ",
     $   npoints," micro structures."
     
      endtime = MPI_Wtime(ierr)
      write(*,*) myrank," took ",endtime-starttime,"s in MAIN."

      CALL MPI_FINALIZE(ierr)
      
      end

c  **********************************************************

      subroutine femat_mpi(nx,ny,nz,phasemod,d1,d2,vox,b,dk,C,cmod)

      implicit none
      
      include 'mpif.h'
      
      integer i,ierr,nx,j,ny,nz
      integer d1,d2,myrank,nprocs
      integer ipx,ipy,ipz
      integer nxy,k,nm,ijk,mm,nn,ii,jj,kk,ll
      integer i3,i8,dn,m,m3,m8
      integer pflag,nphase

      integer status(MPI_STATUS_SIZE)
            
      integer*2 vox(nx,ny,d1-1:d2+1)
            
      double precision sum_num,cterm,cpos,cneg
      double precision c,c3,x,y,z
      double precision f_time(24)
            
      double precision dk(nphase,8,3,8,3)
      double precision dndx(8),dndy(8),dndz(8)
      double precision g(3,3,3),ck(6,6),cmu(6,6),cmod(nphase,6,6)
      double precision es(6,8,3),delta(8,3)
      double precision b(nx,ny,d1-1:d2+1,3)
      double precision,allocatable :: ab(:,:,:),ba(:,:,:)
      double precision exx,eyy,ezz,exz,eyz,exy

      double precision phasemod(nphase,2)
            
      common/list1/pflag,nphase
            
      common/list2/exx,eyy,ezz,exz,eyz,exy
            
      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )
            
      f_time(1) = MPI_Wtime(ierr)
      nxy=nx*ny
      
      allocate (ab(nx,ny,3))
      allocate (ba(nx,ny,3))
            
c  initialize stiffness matrices
            
      dk=0.0d0

c  set up elastic moduli matrices for each kind of element
c  ck and cmu are the bulk and shear modulus matrices, which need to be
c  weighted by the actual bulk and shear moduli

      ck(1,1)=1.0d0
      ck(1,2)=1.0d0
      ck(1,3)=1.0d0
      ck(1,4)=0.0d0
      ck(1,5)=0.0d0
      ck(1,6)=0.0d0
      ck(2,1)=1.0d0
      ck(2,2)=1.0d0
      ck(2,3)=1.0d0
      ck(2,4)=0.0d0
      ck(2,5)=0.0d0
      ck(2,6)=0.0d0
      ck(3,1)=1.0d0
      ck(3,2)=1.0d0
      ck(3,3)=1.0d0
      ck(3,4)=0.0d0
      ck(3,5)=0.0d0
      ck(3,6)=0.0d0
      ck(4,1)=0.0d0
      ck(4,2)=0.0d0
      ck(4,3)=0.0d0
      ck(4,4)=0.0d0
      ck(4,5)=0.0d0
      ck(4,6)=0.0d0
      ck(5,1)=0.0d0
      ck(5,2)=0.0d0
      ck(5,3)=0.0d0
      ck(5,4)=0.0d0
      ck(5,5)=0.0d0
      ck(5,6)=0.0d0
      ck(6,1)=0.0d0
      ck(6,2)=0.0d0
      ck(6,3)=0.0d0
      ck(6,4)=0.0d0
      ck(6,5)=0.0d0
      ck(6,6)=0.0d0

      cmu(1,1)=4.0d0/3.0d0
      cmu(1,2)=-2.0d0/3.0d0
      cmu(1,3)=-2.0d0/3.0d0
      cmu(1,4)=0.0d0
      cmu(1,5)=0.0d0
      cmu(1,6)=0.0d0
      cmu(2,1)=-2.0d0/3.0d0
      cmu(2,2)=4.0d0/3.0d0
      cmu(2,3)=-2.0d0/3.0d0
      cmu(2,4)=0.0d0
      cmu(2,5)=0.0d0
      cmu(2,6)=0.0d0
      cmu(3,1)=-2.0d0/3.0d0
      cmu(3,2)=-2.0d0/3.0d0
      cmu(3,3)=4.0d0/3.0d0
      cmu(3,4)=0.0d0
      cmu(3,5)=0.0d0
      cmu(3,6)=0.0d0
      cmu(4,1)=0.0d0
      cmu(4,2)=0.0d0
      cmu(4,3)=0.0d0
      cmu(4,4)=1.0d0
      cmu(4,5)=0.0d0
      cmu(4,6)=0.0d0
      cmu(5,1)=0.0d0
      cmu(5,2)=0.0d0
      cmu(5,3)=0.0d0
      cmu(5,4)=0.0d0
      cmu(5,5)=1.0d0
      cmu(5,6)=0.0d0
      cmu(6,1)=0.0d0
      cmu(6,2)=0.0d0
      cmu(6,3)=0.0d0
      cmu(6,4)=0.0d0
      cmu(6,5)=0.0d0
      cmu(6,6)=1.0d0
      
      do k=1,nphase
      do j=1,6
      do i=1,6
      
      cmod(k,i,j)=phasemod(k,1)*ck(i,j)+phasemod(k,2)*cmu(i,j)
      
      end do; end do; end do
            
c  set up Simpson’s integration rule weight vector
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
            
c  loop over the nphase kinds of pixels and Simpson’s rule quadrature
c  points in order to compute the stiffness matrices. Stiffness matrices
c  of trilinear finite elements are quadratic in x, y, and z, so that
c  Simpson’s rule quadrature gives exact results.

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
      
c  now build strain matrix
      
      es=0.0d0

      es(1,:,1)=dndx
      es(2,:,2)=dndy
      es(3,:,3)=dndz
      es(4,:,1)=dndz
      es(4,:,3)=dndx
      es(5,:,2)=dndz
      es(5,:,3)=dndy
      es(6,:,1)=dndy

      es(6,:,2)=dndx

c  Matrix multiply to determine value at (x,y,z), multiply by
c  proper weight, and sum_num into dk, the stiffness matrix
      
      f_time(2) = MPI_Wtime(ierr)
      
      do mm=1,3
      do nn=1,3
      do ii=1,8
      do jj=1,8

c  Define sum over strain matrices and elastic moduli matrix for
c  stiffness matrix
      sum_num=0.0d0
      do kk=1,6
      do ll=1,6
      sum_num=sum_num+es(kk,ii,mm)*cmod(ijk,kk,ll)*es(ll,jj,nn)
      end do; end do
      dk(ijk,ii,mm,jj,nn)=dk(ijk,ii,mm,jj,nn)+g(i,j,k)*sum_num/216.
      
      end do; end do; end do; end do
      end do; end do; end do; end do

      f_time(3) = MPI_Wtime(ierr)
      
      if (pflag.eq.1) then
      write(*,*) myrank,"time to calculate dk = ",f_time(3)-f_time(2)
      endif

c  Initialize b and C
      if (myrank.eq.0) then
            write(*,*) "Initializing b & C."
      end if

      b=0.0d0
      C=0.0d0
      c3=0.0d0

999   format(4(i4,1x),3(f9.6,1x))

c  x=nx face
      
      do i3=1,3
      do i8=1,8
      delta(i8,i3) = 0.0d0

      if(i8.eq.2.or.i8.eq.3.or.i8.eq.6.or.i8.eq.7) then
            delta(i8,1)=exx*nx
            delta(i8,2)=exy*nx
            delta(i8,3)=exz*nx
      end if
      
      end do; end do

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

      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m3=1,3
      do m8=1,8

      cterm =0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)
      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if
      
      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      end do; end do

c  Assign b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn) + sum_num

      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num

      end do; end do
      end do; end do

c  y=ny face

      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0d0
      if(i8.eq.3.or.i8.eq.4.or.i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      end do; end do

      do i=1,nx-1
      do k=d1,dn
      m=nxy*(k-1)+nx*(ny-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m3=1,3
      do m8=1,8

      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)

      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if
      
      end do; end do

      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num

      end do; end do
      end do; end do

c  Zface calcs

c  Only the last node does these series of calculations since
c  it contains all the necessary data therefore no data transfer
c  occurs.

      if (myrank.eq.nprocs-1) then
      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0d0
      if(i8.eq.5.or.i8.eq.6.or.i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      end do; end do

      do i=1,nx-1
      do j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m3=1,3
      do m8=1,8
      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)
      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if

      end do; end do
      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num
      end do; end do
      end do; end do

      end if

c  x=nx y=ny edge

      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.2.or.i8.eq.6) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if

      if(i8.eq.4.or.i8.eq.8) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.3.or.i8.eq.7) then
      delta(i8,1)=exy*ny+exx*nx
      delta(i8,2)=eyy*ny+exy*nx
      delta(i8,3)=eyz*ny+exz*nx
      end if
      end do; end do

      dn=d2
      if (dn.eq.nz) then
            dn = nz-1
      end if
      
      do k=d1,dn
      m=nxy*k
      call m2ijk(m,ii,jj,kk,nx,ny,nz)

      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)

      sum_num=0.0d0
      do m3=1,3
      do m8=1,8
      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)

      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if

      end do; end do
      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num

      end do; end do
      end do

c  x=nx z=nz edge

      if (myrank.eq.nprocs-1) then

      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0d0
      if(i8.eq.2.or.i8.eq.3) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
      if(i8.eq.5.or.i8.eq.8) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.6.or.i8.eq.7) then
      delta(i8,1)=exz*nz+exx*nx
      delta(i8,2)=eyz*nz+exy*nx
      delta(i8,3)=ezz*nz+exz*nx
      end if
      end do; end do

      do j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+nx
      call m2ijk(m,ii,jj,kk,nx,ny,nx)

      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m3=1,3
      do m8=1,8

      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)

      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if

      end do; end do

      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num

      end do; end do
      end do

c  y=ny z=nz edge

      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0d0
      if(i8.eq.5.or.i8.eq.6) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.3.or.i8.eq.4) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exy*ny+exz*nz
      delta(i8,2)=eyy*ny+eyz*nz
      delta(i8,3)=eyz*ny+ezz*nz
      end if
      end do; end do

      do i=1,nx-1
      m=nxy*(nz-1)+nx*(ny-1)+i
      call m2ijk(m,ii,jj,kk,nx,ny,nx)
      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)
      sum_num=0.0d0
      do m3=1,3
      do m8=1,8

      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)

      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if

      end do; end do
      
      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num
      end do; end do
      end do

c  x=nx y=ny z=nz corner

      do i3=1,3
      do i8=1,8
      delta(i8,i3)=0.0d0
      if(i8.eq.2) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
      if(i8.eq.4) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.5) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.8) then
      delta(i8,1)=exy*ny+exz*nz
      delta(i8,2)=eyy*ny+eyz*nz
      delta(i8,3)=eyz*ny+ezz*nz
      end if
      if(i8.eq.6) then
      delta(i8,1)=exx*nx+exz*nz
      delta(i8,2)=exy*nx+eyz*nz
      delta(i8,3)=exz*nx+ezz*nz
      end if
      if(i8.eq.3) then
      delta(i8,1)=exx*nx+exy*ny
      delta(i8,2)=exy*nx+eyy*ny
      delta(i8,3)=exz*nx+eyz*ny
      end if
      if(i8.eq.7) then
      delta(i8,1)=exx*nx+exy*ny+exz*nz
      delta(i8,2)=exy*nx+eyy*ny+eyz*nz
      delta(i8,3)=exz*nx+eyz*ny+ezz*nz
      end if
      end do; end do

      m=nx*ny*nz
      call m2ijk(m,ii,jj,kk,nx,ny,nx)
      do nn=1,3
      do mm=1,8
      call ipxyz(mm,ii,jj,kk,ipx,ipy,ipz,nx,ny,nz)

      sum_num=0.0d0
      do m3=1,3
      do m8=1,8
      sum_num=sum_num+delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)
      cterm=0.5d0*delta(m8,m3)*dk(vox(ii,jj,kk),m8,m3,mm,nn)*delta(mm,
     $   nn)

      if (cterm.ge.0.0d0) then
            cpos = cpos+cterm
      else
            cneg = cneg+cterm
      end if

      end do; end do
      b(ipx,ipy,ipz,nn) = b(ipx,ipy,ipz,nn)+sum_num

      end do; end do

c  End if for (myrank.eq.nprocs-1)

      end if

      c3 = cpos+cneg
      CALL MPI_ALLREDUCE(c3,C,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      if (myrank.eq.0) then
            write(*,*) "Final C = ",C
      end if

      f_time(4) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*)myrank,"Etime to calculate C & b= ",f_time(4)-f_time(3)
      end if

      if (nprocs.gt.1) then

c  RECV a new slice per node.

      ab = 0.0d0
      ba = b(:,:,d2+1,:)

      f_time(5) = MPI_Wtime(ierr)
      call t2b_dp(ab,ba,nx,ny,3)
      f_time(6) = MPI_Wtime(ierr)
      b(:,:,d1,:) = b(:,:,d1,:)+ab

      if (pflag.eq.1) then
      write(*,*) myrank," B upddate: t2b time= ",f_time(6)-f_time(5)
      end if

c  botp = d1-1

      ab = 0.0
      ba = b(:,:,d1-1,:)

      f_time(7) = MPI_Wtime(ierr)
      call b2t_dp(ab,ba,nx,ny,3)
      f_time(8) = MPI_Wtime(ierr)
      b(:,:,d2,:) = b(:,:,d2,:)+ab

      if (pflag.eq.1) then
      write(*,*) myrank," B upddate: b2t time= ",f_time(8)-f_time(7)
      end if

c  Update ghost layers

c  RECV a new slice per node.

      ab = b(:,:,d1,:)
      ba = b(:,:,d2,:)

      f_time(9) = MPI_Wtime(ierr)
      call t2b_dp(ab,ba,nx,ny,3)
      f_time(10) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank,"B ghost upddate:t2b time= ",
     $   f_time(10)-f_time(9)
      end if

      b(:,:,d1-1,:) = ab
      
      ab = b(:,:,d1,:)
      ba = b(:,:,d2,:)

      f_time(11) = MPI_Wtime(ierr)
      call b2t_dp(ab,ba,nx,ny,3)
      f_time(12) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank,"B ghost upddate:b2t time= ",
     $   f_time(12)-f_time(11)
      end if

      b(:,:,d2+1,:) = ba

      else

c  nprocs=1

      b(:,:,d1,:) = b(:,:,d1,:)+b(:,:,d2+1,:)
      b(:,:,d2,:) = b(:,:,d2,:)+b(:,:,d1-1,:)
      b(:,:,d2+1,:) = b(:,:,d1,:)
      b(:,:,d1-1,:) = b(:,:,d2,:)

      end if

      deallocate(ab)
      deallocate(ba)

      f_time(13) = MPI_Wtime(ierr)

      if (pflag.eq.1) then
      write(*,*) myrank,"Femat_mpi elapsed time= ",f_time(13)-f_time(1)
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

      double precision u(nx,ny,d1-1:d2+1,3)
      double precision b(nx,ny,d1-1:d2+1,3)
      double precision gb(nx,ny,d1-1:d2+1,3)
      integer*2 vox(nx,ny,d1-1:d2+1)
      double precision e_time(24)

      double precision c,c3,utot
      double precision dk(nphase,8,3,8,3)

      double precision dutot

      double precision exx,eyy,ezz,exz,eyz,exy

      common/list1/pflag,nphase
      common/list2/exx,eyy,ezz,exz,eyz,exy

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      e_time(1) = MPI_Wtime(ierr)

      dutot = 0.0d0
c  After this call, gb is calculated and data slabs
c  are updated and passed.

      call gbah(gb,u,dk,vox,nx,ny,nz,d1,d2)

c  Now do the rest of the gb calculations that appear
c  in original "energy" subroutine.

c  utot will be a per processor value.
c  Do an MPI_ALLREDUCE on dutot
c  so each node will have the current updated version.

      dutot=0.0d0
      do m3=1,3
      do ik=d1,d2
      do ij=1,ny
      do ii=1,nx

      dutot=dutot+0.5d0*u(ii,ij,ik,m3)*gb(ii,ij,ik,m3)+b(ii,ij,ik,
     $   m3)*u(ii,ij,ik,m3)
      end do; end do; end do; end do

      call MPI_ALLREDUCE(dutot,utot,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      utot = utot+C
c  easier to add C here than before the above MPI call.

      gb = gb+b

      return
      end

c  **********************************************************

      subroutine dembx_mpi(nx,ny,nz,d1,d2,Lstep,gb,u,vox,h,gg,dk,gtest,
     $   ldemb,kkk)

      implicit none

      include 'mpif.h'

      integer nx,ny,nz,d1,d2,ldemb,kkk,ijk
      integer Lstep,myrank,nprocs,ierr
      integer pflag,nphase

      double precision dgg,gg,gglast,lambda,hAh2,hAh,gamma,gtest

      double precision u(nx,ny,d1-1:d2+1,3)
      double precision gb(nx,ny,d1-1:d2+1,3)
      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision dk(nphase,8,3,8,3)

      double precision Ah(nx,ny,d1-1:d2+1,3)
      double precision h(nx,ny,d1-1:d2+1,3)

      common/list1/pflag,nphase

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

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

      hAh2 = SUM(h(:,:,d1:d2,:)*Ah(:,:,d1:d2,:))

      call MPI_ALLREDUCE(hAh2,hAh,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      lambda=gg/hAh
      u=u-lambda*h
      gb=gb-lambda*Ah

      gglast=gg

      gg=0.0d0

      dgg = SUM(gb(:,:,d1:d2,:)*gb(:,:,d1:d2,:))
      call MPI_ALLREDUCE(dgg,gg,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      if (gg.lt.gtest) goto 1000

      gamma = gg/gglast
      h = gb+gamma*h

      end do
1000  continue

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      return
      end

c  **********************************************************

      subroutine stress_mpi(nx,ny,nz,ns,u,vox,cmod,d1,d2)

      implicit none
      include 'mpif.h'
      
      integer nx,ny,ns,d1,d2
      integer ifxa,ifya
      integer pflag,nphase
      
      double precision u(nx,ny,d1-1:d2+1,3),uu(8,3)
      double precision dndx(8),dndy(8),dndz(8)
      double precision es(6,8,3),cmod(nphase,6,6)
      integer*2 vox(nx,ny,d1-1:d2+1)
      integer myrank,ierr,nprocs
      integer status(MPI_STATUS_SIZE)
      
      integer nz,nxy,i,j,k,m,mm,n3,n8,n
      
      double precision strxx,stryy,strzz,strxz,stryz,strxy
      double precision str11,str22,str33,str13,str23,str12
      double precision strxxp,stryyp,strzzp

      double precision strxzp,stryzp,strxyp

      double precision s11,s22,s33,s13,s23,s12
      double precision sxx,syy,szz,sxz,syz,sxy
      double precision sxxp,syyp,szzp,sxzp,syzp,sxyp
      double precision exx,eyy,ezz,exz,eyz,exy

      common/list1/pflag,nphase
      common/list2/exx,eyy,ezz,exz,eyz,exy
      common/list3/strxxp,stryyp,strzzp,strxyp,strxzp,stryzp
      common/list4/sxxp,syyp,szzp,sxyp,sxzp,syzp

            call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
            call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      nxy=nx*ny

c  set up single element strain matrix
c  dndx, dndy, and dndz are the components of the average strain
c  matrix in a pixel

      dndx(1)=-0.25d0
      dndx(2)=0.25d0
      dndx(3)=0.25d0
      dndx(4)=-0.25d0
      dndx(5)=-0.25d0
      dndx(6)=0.25d0
      dndx(7)=0.25d0
      dndx(8)=-0.25d0
      dndy(1)=-0.25d0
      dndy(2)=-0.25d0
      dndy(3)=0.25d0
      dndy(4)=0.25d0
      dndy(5)=-0.25d0
      dndy(6)=-0.25d0
      dndy(7)=0.25d0
      dndy(8)=0.25d0
      dndz(1)=-0.25d0
      dndz(2)=-0.25d0
      dndz(3)=-0.25d0
      dndz(4)=-0.25d0
      dndz(5)=0.25d0
      dndz(6)=0.25d0
      dndz(7)=0.25d0
      dndz(8)=0.25d0
c  Build averaged strain matrix, follows code in femat, but for average
c  strain over the pixel, not the strain at a point.
      
      es = 0.0d0
      es(1,:,1)=dndx
      es(2,:,2)=dndy
      es(3,:,3)=dndz
      es(4,:,1)=dndz
      es(4,:,3)=dndx
      es(5,:,2)=dndz
      es(5,:,3)=dndy
      es(6,:,1)=dndy
      es(6,:,2)=dndx

c  Compute components of the average stress and strain tensors in each pixel
      strxx=0.0d0
      stryy=0.0d0
      strzz=0.0d0
      strxz=0.0d0
      stryz=0.0d0
      strxy=0.0d0
      sxx=0.0d0
      syy=0.0d0
      szz=0.0d0
      sxz=0.0d0
      syz=0.0d0
      sxy=0.0d0

      strxxp=0.0d0
      stryyp=0.0d0
      strzzp=0.0d0
      strxzp=0.0d0
      stryzp=0.0d0
      strxyp=0.0d0
      sxxp=0.0d0
      syyp=0.0d0
      szzp=0.0d0
      sxzp=0.0d0
      syzp=0.0d0
      sxyp=0.0d0

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

      do mm=1,3
      uu(1,mm)= u(i,j,k,mm)
      uu(2,mm)= u(ifxa,j,k,mm)
      uu(3,mm)= u(ifxa,ifya,k,mm)
      uu(4,mm)= u(i,ifya,k,mm)
      uu(5,mm)= u(i,j,k+1,mm)
      uu(6,mm)= u(ifxa,j,k+1,mm)
      uu(7,mm)= u(ifxa,ifya,k+1,mm)
      uu(8,mm)= u(i,ifya,k+1,mm)
      end do

c  Correct for periodic boundary conditions, some displacements are wrong
c  for a pixel on a periodic boundary. Since they come from an opposite
c  face, need to put in applied strain to correct them.
      if(i.eq.nx) then
      uu(2,1)=uu(2,1)+exx*nx
      uu(2,2)=uu(2,2)+exy*nx
      uu(2,3)=uu(2,3)+exz*nx
      uu(3,1)=uu(3,1)+exx*nx
      uu(3,2)=uu(3,2)+exy*nx
      uu(3,3)=uu(3,3)+exz*nx
      uu(6,1)=uu(6,1)+exx*nx
      uu(6,2)=uu(6,2)+exy*nx
      uu(6,3)=uu(6,3)+exz*nx
      uu(7,1)=uu(7,1)+exx*nx
      uu(7,2)=uu(7,2)+exy*nx
      uu(7,3)=uu(7,3)+exz*nx
      end if
      if(j.eq.ny) then
      uu(3,1)=uu(3,1)+exy*ny
      uu(3,2)=uu(3,2)+eyy*ny
      uu(3,3)=uu(3,3)+eyz*ny
      uu(4,1)=uu(4,1)+exy*ny
      uu(4,2)=uu(4,2)+eyy*ny
      uu(4,3)=uu(4,3)+eyz*ny
      uu(7,1)=uu(7,1)+exy*ny
      uu(7,2)=uu(7,2)+eyy*ny
      uu(7,3)=uu(7,3)+eyz*ny
      uu(8,1)=uu(8,1)+exy*ny
      uu(8,2)=uu(8,2)+eyy*ny
      uu(8,3)=uu(8,3)+eyz*ny
      end if
      if(k.eq.nz) then
      uu(5,1)=uu(5,1)+exz*nz
      uu(5,2)=uu(5,2)+eyz*nz
      uu(5,3)=uu(5,3)+ezz*nz
      uu(6,1)=uu(6,1)+exz*nz
      uu(6,2)=uu(6,2)+eyz*nz
      uu(6,3)=uu(6,3)+ezz*nz
      uu(7,1)=uu(7,1)+exz*nz
      uu(7,2)=uu(7,2)+eyz*nz
      uu(7,3)=uu(7,3)+ezz*nz
      uu(8,1)=uu(8,1)+exz*nz
      uu(8,2)=uu(8,2)+eyz*nz
      uu(8,3)=uu(8,3)+ezz*nz
      end if

c  local stresses and strains in a pixel
      str11=0.0d0
      str22=0.0d0
      str33=0.0d0
      str13=0.0d0
      str23=0.0d0
      str12=0.0d0
      s11=0.0d0
      s22=0.0d0
      s33=0.0d0
      s13=0.0d0
      s23=0.0d0
      s12=0.0d0

      do n3=1,3
      do n8=1,8
      s11=s11+es(1,n8,n3)*uu(n8,n3)
      s22=s22+es(2,n8,n3)*uu(n8,n3)
      s33=s33+es(3,n8,n3)*uu(n8,n3)
      s13=s13+es(4,n8,n3)*uu(n8,n3)
      s23=s23+es(5,n8,n3)*uu(n8,n3)
      s12=s12+es(6,n8,n3)*uu(n8,n3)

      do n=1,6
            str11=str11+cmod(vox(i,j,k),1,n)*es(n,n8,n3)*uu(n8,n3)
            str22=str22+cmod(vox(i,j,k),2,n)*es(n,n8,n3)*uu(n8,n3)
            str33=str33+cmod(vox(i,j,k),3,n)*es(n,n8,n3)*uu(n8,n3)
            str13=str13+cmod(vox(i,j,k),4,n)*es(n,n8,n3)*uu(n8,n3)
            str23=str23+cmod(vox(i,j,k),5,n)*es(n,n8,n3)*uu(n8,n3)
            str12=str12+cmod(vox(i,j,k),6,n)*es(n,n8,n3)*uu(n8,n3)
      end do

      end do; end do

c  sum local strains and stresses into global values

      strxx=strxx+str11
      stryy=stryy+str22
      strzz=strzz+str33
      strxz=strxz+str13
      stryz=stryz+str23
      strxy=strxy+str12
      sxx=sxx+s11
      syy=syy+s22
      szz=szz+s33
      sxz=sxz+s13
      syz=syz+s23
      sxy=sxy+s12
470   continue

c  Now do MPI to gather all strNN and sNN terms,
c  add them at root, then do this final calculation
c  and write them to disk.

      call MPI_ALLREDUCE(strxx,strxxp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(stryy,stryyp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(strzz,strzzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(strxz,strxzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(strxy,strxyp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(stryz,stryzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(sxx,sxxp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(syy,syyp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(szz,szzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)
      
      call MPI_ALLREDUCE(sxz,sxzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(sxy,sxyp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)
      
      call MPI_ALLREDUCE(syz,syzp,1,MPI_double_precision,MPI_SUM,
     $   MPI_COMM_WORLD,ierr)

c  Now root has strxx,stryy, ... sxy
c  Let him write out data to disk after this Volume averaging

c  Volume average of global stresses and strains
      strxxp=strxxp/dfloat(ns)
      stryyp=stryyp/dfloat(ns)
      strzzp=strzzp/dfloat(ns)
      strxzp=strxzp/dfloat(ns)
      stryzp=stryzp/dfloat(ns)
      strxyp=strxyp/dfloat(ns)
      sxxp=sxxp/dfloat(ns)
      syyp=syyp/dfloat(ns)
      szzp=szzp/dfloat(ns)
      sxzp=sxzp/dfloat(ns)
      syzp=syzp/dfloat(ns)
      sxyp=sxyp/dfloat(ns)

      if (myrank.eq.0) then

      write(*,*) "strxxp = ",strxxp
      write(*,*) "stryyp = ",stryyp
      write(*,*) "strzzp = ",strzzp
      write(*,*) "strxyp = ",strxyp
      write(*,*) "strxzp = ",strxzp
      write(*,*) "stryzp = ",stryzp

      write(*,*) "sxxp = ",sxxp
      write(*,*) "syyp = ",syyp
      write(*,*) "szzp = ",szzp
      write(*,*) "sxyp = ",sxyp
      write(*,*) "sxzp = ",sxzp
      write(*,*) "syzp = ",syzp

      end if

      return
      end

c  **********************************************************

      subroutine gbah(om,uh,dk,vox,nx,ny,nz,d1,d2)
      implicit none

      include 'mpif.h'

      integer nx,ny,nz,d1,d2
      integer im,jm,km,j,ifxa,ifxb,ifya,ifyb
      integer myrank,nprocs,ierr
      integer pflag,nphase

      double precision uh(nx,ny,d1-1:d2+1,3)
      double precision om(nx,ny,d1-1:d2+1,3)
      double precision gb_time(6)
      
      integer*2 vox(nx,ny,d1-1:d2+1)

      double precision dk(nphase,8,3,8,3)

      common/list1/pflag,nphase

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

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

      do j=1,3

c  SELF TERM

      om(im,jm,km,j) =

c  u(ib(m,1),n)

     & SUM ( uh(im,ifya,km,:)*(dk(vox(im,jm,km),1,j,4,:)+dk(vox(ifxb,jm,
     &  km),2,j,3,:)+dk(vox(im,jm,km-1),5,j,8,:)+dk(vox(ifxb,jm,km-1),6,
     &  j,7,:) ))+

c  u(ib(m,2),n)
     & SUM ( uh(ifxa,ifya,km,:)*(dk(vox(im,jm,km),1,j,3,:)+dk(vox(im,jm,
     &  km-1),5,j,7,:) ))+

c  u(ib(m,3),n)

     & SUM ( uh(ifxa,jm,km,:)*(dk(vox(im,jm,km),1,j,2,:)+dk(vox(im,ifyb,
     &  km),4,j,3,:)+dk(vox(im,ifyb,km-1),8,j,7,:)+dk(vox(im,jm,km-1),5,
     &  j,6,:) ))+

c  u(ib(m,4),n)
     & SUM ( uh(ifxa,ifyb,km,:)*(dk(vox(im,ifyb,km),4,j,2,:)+dk(vox(im,
     &  ifyb,km-1),8,j,6,:) ))+

c  u(ib(m,5),n)
     & SUM ( uh(im,ifyb,km,:)*(dk(vox(ifxb,ifyb,km),3,j,2,:)+dk(vox(im,
     &  ifyb,km),4,j,1,:)+dk(vox(ifxb,ifyb,km-1),7,j,6,:)+dk(vox(im,
     &  ifyb,km-1),8,j,5,:) ))+

c  u(ib(m,6),n)
     & SUM ( uh(ifxb,ifyb,km,:)*(dk(vox(ifxb,ifyb,km),3,j,1,
     &  :)+dk(vox(ifxb,ifyb,km-1),7,j,5,:) ))+

c  u(ib(m,7),n)
     & SUM(uh(ifxb,jm,km,:)*(dk(vox(ifxb,ifyb,km),3,j,4,:)+dk(vox(ifxb,
     &  jm,km),2,j,1,:)+dk(vox(ifxb,ifyb,km-1),7,j,8,:)+dk(vox(ifxb,jm,
     &  km-1),6,j,5,:) ))+

c  u(ib(m,8),n)
     & SUM (uh(ifxb,ifya,km,:)*( dk(vox(ifxb,jm,km),2,j,4,
     &  :)+dk(vox(ifxb,jm,km-1),6,j,8,:) ))+
     
c  u(ib(m,9),n)
     & SUM ( uh(im,ifya,km-1,:)*(dk(vox(im,jm,km-1),5,j,4,
     &  :)+dk(vox(ifxb,jm,km-1),6,j,3,:) ))+

c  u(ib(m,10),n)
     & SUM ( uh(ifxa,ifya,km-1,:)*(dk(vox(im,jm,km-1),5,j,3,:) ))+

c  u(ib(m,11),n)
     & SUM ( uh(ifxa,jm,km-1,:)*(dk(vox(im,ifyb,km-1),8,j,3,
     &  :)+dk(vox(im,jm,km-1),5,j,2,:) ))+

c  u(ib(m,12),n)
     & SUM( uh(ifxa,ifyb,km-1,:)*( dk(vox(im,ifyb,km-1),8,j,2,:) ))+

c  u(ib(m,13),n)
     & SUM ( uh(im,ifyb,km-1,:)*(dk(vox(im,ifyb,km-1),8,j,1,
     &  :)+dk(vox(ifxb,ifyb,km-1),7,j,2,:) ))+

c  u(ib(m,14),n)
     & SUM( uh(ifxb,ifyb,km-1,:)*( dk(vox(ifxb,ifyb,km-1),7,j,1,:) ))+

c  u(ib(m,15),n)
     & SUM ( uh(ifxb,jm,km-1,:)*(dk(vox(ifxb,ifyb,km-1),7,j,4,
     &  :)+dk(vox(ifxb,jm,km-1),6,j,1,:) ))+

c  u(ib(m,16),n)
     &SUM(uh(ifxb,ifya,km-1,:)*( dk(vox(ifxb,jm,km-1),6,j,4,:) ))+

c  u(ib(m,17),n)
     & SUM ( uh(im,ifya,km+1,:)*(dk(vox(im,jm,km),1,j,8,:)+dk(vox(ifxb,
     &  jm,km),2,j,7,:) ))+

c  u(ib(m,18),n)
     & SUM (uh(ifxa,ifya,km+1,:)*( dk(vox(im,jm,km),1,j,7,:) ))+

c  u(ib(m,19),n)
     & SUM ( uh(ifxa,jm,km+1,:)*(dk(vox(im,jm,km),1,j,6,:)+dk(vox(im,
     &  ifyb,km),4,j,7,:) ))+

c  u(ib(m,20),n)
     & SUM (uh(ifxa,ifyb,km+1,:)*( dk(vox(im,ifyb,km),4,j,6,:) ))+

c  u(ib(m,21),n)
     & SUM ( uh(im,ifyb,km+1,:)*(dk(vox(im,ifyb,km),4,j,5,
     &  :)+dk(vox(ifxb,ifyb,km),3,j,6,:) ))+

c  u(ib(m,22),n)
     & SUM(uh(ifxb,ifyb,km+1,:)*( dk(vox(ifxb,ifyb,km),3,j,5,:) ))+

c  u(ib(m,23),n)
     & SUM ( uh(ifxb,jm,km+1,:)*(dk(vox(ifxb,ifyb,km),3,j,8,
     &  :)+dk(vox(ifxb,jm,km),2,j,5,:) ))+

c  u(ib(m,24),n)
     & SUM(uh(ifxb,ifya,km+1,:)*( dk(vox(ifxb,jm,km),2,j,8,:) ))+

c  u(ib(m,25),n)
     & SUM ( uh(im,jm,km-1,:)*(dk(vox(ifxb,ifyb,km-1),7,j,3,
     &  :)+dk(vox(im,ifyb,km-1),8,j,4,:)+dk(vox(ifxb,jm,km-1),6,j,2,
     &  :)+dk(vox(im,jm,km-1),5,j,1,:) ))+

c  u(ib(m,26),n)
     & SUM(uh(im,jm,km+1,:)*(dk(vox(ifxb,ifyb,km),3,j,7,:)+dk(vox(im,
     &  ifyb,km),4,j,8,:)+dk(vox(im,jm,km),1,j,5,:)+dk(vox(ifxb,jm,km),
     &  2,j,6,:) ))+
     
c  u(ib(m,27),n)
     & SUM( uh(im,jm,km,:)* (dk(vox(im,jm,km),1,j,1,:)+dk(vox(ifxb,jm,
     &  km),2,j,2,:)+dk(vox(ifxb,ifyb,km),3,j,3,:)+dk(vox(im,ifyb,km),4,
     &  j,4,:)+dk(vox(im,jm,km-1),5,j,5,:)+dk(vox(ifxb,jm,km-1),6,j,6,
     &  :)+dk(vox(ifxb,ifyb,km-1),7,j,7,:)+dk(vox(im,ifyb,km-1),8,j,8,
     &  :) ))

      end do
 
      end do; end do; end do

      gb_time(2) = MPI_Wtime(ierr)
      
      if (pflag.eq.1) then
      write(*,*)myrank,"Etime to calc gb/Ah=",gb_time(2)-gb_time(1)
      endif

c  Do top/bottom layer switch on matrix: om

      call z_ghost_dp(om,nx,ny,3,d1,d2)

      if (pflag.eq.1) then
      write(*,*)myrank,"Etime for t2b gb/Ah=",gb_time(4)-gb_time(3)
      write(*,*)myrank,"Etime for b2t gb/Ah=",gb_time(6)-gb_time(5)
      endif

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

      if(pix(i,j,k).eq.0) then
            pix(i,j,k)=46
      end if

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

      subroutine z_ghost_int(arr0,mx,my,mz,d1,d2)

      implicit none

      include 'mpif.h'
      
      integer mx,my,mz,d1,d2
      
      integer*2 arr0(mx,my,d1-1:d2+1)
      integer*2,allocatable :: bot(:,:),top(:,:)

      integer myrank,ierr,nprocs
      integer status(MPI_STATUS_SIZE)
      
      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

c  Make the Z Ghost

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

      subroutine z_ghost_dp(arr0,mx,my,mz,d1,d2)

      implicit none
      
      include 'mpif.h'
      
      integer mx,my,mz,d1,d2
      
      double precision arr0(mx,my,d1-1:d2+1,mz)
      
      double precision,allocatable :: bot(:,:,:),top(:,:,:)

      integer myrank,ierr,nprocs
      integer status(MPI_STATUS_SIZE)
      
      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )
      
c  Make the Z Ghost
      
      allocate(bot(mx,my,mz))
      allocate(top(mx,my,mz))
      
c  Get new bottom ghost plane.
      
      bot = arr0(:,:,d1,:)
      top = arr0(:,:,d2,:)
      call t2b_dp(bot,top,mx,my,3)
      arr0(:,:,d1-1,:) = bot
      
c  Get new top ghost plane
      
      bot = arr0(:,:,d1,:)
      top = arr0(:,:,d2,:)
      call b2t_dp(bot,top,mx,my,3)
      arr0(:,:,d2+1,:) = top
      deallocate(bot)
      deallocate(top)
      
      return
      end
      
c  **********************************************************

      subroutine t2b(b_layer,t_layer,nx,ny)
      
c  This is an INTEGER*2 subroutine.

c  Used for transferring: pix bottom2top layes

c  RECV a new t_layer (TOP layer) per node
      
      implicit none
      
      include 'mpif.h'

      integer nx,ny,nxy
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)
      
      integer*2 b_layer(nx,ny),t_layer(nx,ny)
      
      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )
      
      nxy=nx*ny
      
      ides = mod(myrank+1,nprocs)
      isrc = mod(myrank+nprocs-1,nprocs)
      
      if (myrank.eq.nprocs-1) then
      call MPI_Irecv(b_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,
     &   irequest,ierr)
      call mpi_send(t_layer,2*nxy,MPI_BYTE,ides,9,MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)
      else
            
      call mpi_recv(b_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,status,
     &   ierr)
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

      integer*2 b_layer(nx,ny),t_layer(nx,ny)

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      nxy=nx*ny

      ides = mod(myrank+nprocs-1,nprocs)
      isrc = mod(myrank+1,nprocs)

      if (myrank.eq.nprocs-1) then
      call MPI_Irecv(t_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,
     &   irequest,ierr)
      call mpi_send(b_layer,2*nxy,MPI_BYTE,ides,9,MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)
      
      else
      call mpi_recv(t_layer,2*nxy,MPI_BYTE,isrc,9,MPI_COMM_WORLD,status,
     &   ierr)
      call mpi_send(b_layer,2*nxy,MPI_BYTE,ides,9,MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine t2b_dp(b_layer,t_layer,nx,ny,i)
      
c  This is a double precision subroutine.

c  Used for transferring: u,b,and om top2bottom layers

c  RECV a new b_layer (BOTTOM layer) per node.
      
      implicit none
      
      include 'mpif.h'
      
      integer nx,ny,mxy,i
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)
      double precision b_layer(nx,ny,i),t_layer(nx,ny,i)

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      mxy=i*nx*ny

      ides = mod(myrank+1,nprocs)
      isrc = mod(myrank+nprocs-1,nprocs)

      if (myrank.eq.nprocs-1) then
      call mpi_irecv(b_layer,mxy,MPI_double_precision,isrc,9,
     &   MPI_COMM_WORLD,irequest,ierr)
      call mpi_send(t_layer,mxy,MPI_double_precision,ides,9,
     &   MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)
      else

      call mpi_recv(b_layer,mxy,MPI_double_precision,isrc,9,
     &   MPI_COMM_WORLD,status,ierr)
      call mpi_send(t_layer,mxy,MPI_double_precision,ides,9,
     &   MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      return
      end

c  **********************************************************

      subroutine b2t_dp(b_layer,t_layer,nx,ny,i)
      
c  This is a double precision subroutine.

c  Used for transferring: u,b,and om bottom2top layers

c  RECV a new t_layer (TOP layer) per node.

      implicit none

      include 'mpif.h'
      
      integer nx,ny,mxy,i
      integer ides,isrc,irequest
      integer myrank,nprocs,ierr
      integer status(MPI_STATUS_SIZE)

      double precision b_layer(nx,ny,i),t_layer(nx,ny,i)

      call MPI_COMM_RANK( MPI_COMM_WORLD,myrank,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,nprocs,ierr )

      mxy=i*nx*ny

      ides = mod(myrank+nprocs-1,nprocs)
      isrc = mod(myrank+1,nprocs)

      if (myrank.eq.nprocs-1) then
      call mpi_Irecv(t_layer,mxy,MPI_double_precision,isrc,9,
     &   MPI_COMM_WORLD,irequest,ierr)
      call mpi_send(b_layer,mxy,MPI_double_precision,ides,9,
     &   MPI_COMM_WORLD,ierr)
      call MPI_WAIT(irequest,status,ierr)
      else

      call mpi_recv(t_layer,mxy,MPI_double_precision,isrc,9,
     &   MPI_COMM_WORLD,status,ierr)
      call mpi_send(b_layer,mxy,MPI_double_precision,ides,9,
     &   MPI_COMM_WORLD,ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c  **********************************************************

      subroutine phasemod_init(phasemod)

c  USER: Put phasemod definitions here

      implicit none
      integer pflag,nphase,i
      double precision phasemod(nphase,2),saves
      
      common/list1/pflag,nphase
      
      phasemod = 0.0d0

	phasemod(1,1)=76.8
	phasemod(1,2)=32.0
	phasemod(2,1)=2.2
	phasemod(2,2)=0.0

! c  C3S
!       phasemod(1,1)=117.6d0
!       phasemod(1,2)=0.314d0
! c  C2S (same as C3S for now)
!       phasemod(2,1)=117.6d0
!       phasemod(2,2)=0.314d0
! c  C3A (same as C3S for now)
!       phasemod(3,1)=117.6d0
!       phasemod(3,2)=0.314d0
! c  C4AF (same as C3S for now)
!       phasemod(4,1)=117.6d0
!       phasemod(4,2)=0.314d0
! c  gypsum (use from paper with Sylvain)
!       phasemod(5,1)=45.7d0
!       phasemod(5,2)=0.33d0
! c  hemihydrate (same as gypsum for now)
!       phasemod(6,1)=0.5*(45.7d0+80.0d0)
!       phasemod(6,2)=0.5*(0.33d0+0.275d0)
! c  anhydrite (same as gypsum for now)
!       phasemod(7,1)=80.0d0
!       phasemod(7,2)=0.275d0
! c  pozzolan (no pozzolan)
!       phasemod(8,1)=0.0d0
!       phasemod(8,2)=0.0d0
! c  inert
!       phasemod(9,1)=0.0d0
!       phasemod(9,2)=0.0d0
! c  slag
!       phasemod(10,1)=0.0d0
!       phasemod(10,2)=0.0d0
! c  ASG flyash
!       phasemod(11,1)=0.0d0
!       phasemod(11,2)=0.0d0
! c  CAS2 fly ash
!       phasemod(12,1)=0.0d0
!       phasemod(12,2)=0.0d0
! c  CH
!       phasemod(13,1)=42. 3d0
!       phasemod(13,2)=0. 324d0
! c  C-S-H    
!       phasemod(14,1)=22.4d0
!       phasemod(14,2)=0.25d0
! c  C3AH6 (same as C-S-H for now)
!       phasemod(15,1)=phasemod(14,1)
!       phasemod(15,2)=phasemod(14,2)
! c  ettringite (from C3A) (1/3 gypsum for now)
!       phasemod(16,1)=phasemod(14,1)
!       phasemod(16,2)=phasemod(14,2)
! c  ettringite (from C4AF)
!       phasemod(17,1)=phasemod(16,1)
!       phasemod(17,2)=phasemod(16,2)
! c  Afm
!       phasemod(18,1)=phasemod(13,1)
!       phasemod(18,2)=phasemod(13,2)
! c  FH3 (same as C-S-H for now)
!       phasemod(19,1)=phasemod(14,1)
!       phasemod(19,2)=phasemod(14,2)
! c  pozzolanic C-S-H
!       phasemod(20,1)=phasemod(14,1)
!       phasemod(20,2)=phasemod(14,2)
! c  Slag C-S-H
!       phasemod(21,1)=phasemod(14,1)
!       phasemod(21,2)=phasemod(14,2)
! c  CaCl2 (in fly ash)
!       phasemod(22,1)=0.0d0
!       phasemod(22,2)=0.0d0
! c  Friedel Salt
!       phasemod(23,1)=0.0d0
!       phasemod(23,2)=0.0d0
! c  Stratlingite (from fly ash presence)
!       phasemod(24,1)=0.0d0
!       phasemod(24,2)=0.0d0
! c  Secondary gypsum (same modulus as regular gypsum)
!       phasemod(25,1)=phasemod(5,1)
!       phasemod(25,2)=phasemod(5,2)
! c  CaCO3
!       phasemod(26,1)=79.6d0
!       phasemod(26,2)=0.31d0
! c  Afmc
!       phasemod(27,1)=phasemod(13,1)
!       phasemod(27,2)=phasemod(13,2)
! c  Inert aggregate
!       phasemod(28,1)=0.0d0
!       phasemod(28,2)=0.0d0
! c  Absorbed gypsum (in C-S-H) treat as regular gypsum
!       phasemod(29,1)=phasemod(5,1)
!       phasemod(29,2)=phasemod(5,2)
! c  Fly ash
!       phasemod(30,1)=0.0d0
!       phasemod(30,2)=0.0d0
! c  C3A (fly ash)
!       phasemod(35,1)=0.0d0
!       phasemod(35,2)=0.0d0
! c  Empty porosity (no water)
!       phasemod(45,1)=0.0d0
!       phasemod(45,2)=0.0d0
! c  Water-filled porosity (change from label of zero in hydration program)
! c  input as bulk modulus (1) and shear modulus (2), preserve in do 1144 below
!       phasemod(46,1)=2.0d0
!       phasemod(46,2)=0.0d0
! c  Switched off phase for early age.
!       phasemod(88,1)=0.0d0
!       phasemod(88,2)=0.0d0

c  USER: end of phasemod defs

c  (USER) Program uses bulk modulus (1) and shear modulus (2), so transform
c  Young’s modulis (1) and Poisson’s ratio (2).
      
!       do 1144 i=1,nphase

!       if(i.eq.46) goto 1144

!       saves=phasemod(i,1)
!       phasemod(i,1)=phasemod(i,1)/3.d0/(1.d0-2.d0*phasemod(i,2))
!       phasemod(i,2)=saves/2.d0/(1.d0+phasemod(i,2))

! 1144  continue
      return
      end