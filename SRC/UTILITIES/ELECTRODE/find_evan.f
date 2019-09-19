       program find_evan
cccccc AL is the lattice of the large cell
cccccc ALw is the lattice of the electrode unit cell

      implicit double precision (a-h,o-z)

      parameter (nm=2000)
      parameter (nevan=2000)

      real*8, allocatable, dimension (:,:,:) :: phase

      complex*16, allocatable, dimension (:,:,:,:) :: ucw
      complex*16, allocatable, dimension (:,:,:) :: uc
      complex*16, allocatable, dimension (:,:,:) :: cphase
      complex*16, allocatable, dimension (:,:,:) :: ur
      complex*16, allocatable, dimension (:) :: ur_tmp


      real*8 E_evan(nevan),E_evanC(nevan)
      integer ist_evan(nevan),ikpt_evan(nevan),iGX_evan(nevan),
     &     iband_evan(nevan)
       real*8 dE_dk(400),ak_w(400)
      integer nline_w(400)

      real*8 E_linew(nm,nm)
      real*8 ALw(3,3),AL(3,3)

      integer ik1(400,10),ik2(400,10),ikw(400,10)
      real*8 Ekw(400,2)
      integer num1(400),num2(400),numw_new(nm)
      integer ikpt_st1w(400),ikpt_st2w(400),
     &   ikpt_st3w(400),
     &   i_st1w(400),
     &  i_st2w(400),i_st3w(400),
     &  numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)

      real*8 x_st1w(400),x_st2w(400),x_st3w(400)
      real*8 kk1,kk2,kk3,k1k2,k2k3,
     &     k2k1,k3k1,k3k2,k1k3


      complex*16 cc

      character*11 fileh


      ALw(:,1)=AL(:,1)*n1w*1.d0/n1
 
      pi=4*datan(1.d0)
      
      open(10,file="connect.input")
      rewind(10)
      read(10,*) 
      read(10,*) 
      read(10,*) nintep
      close(10)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Real band structure in GX
      open(14,file="E_line_W2K2")
      rewind(14)
      read(14,*) nstw,nkptw
      do j=1,nstw
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      close(14)

cccccc first, the evanescent state at the Gamma point
      num_evan=0
      do j=1,nstw
      if(ikpt_linew(1,j).eq.1.and.numw(j).ge.3*nintep) then    ! Gamma point, forget about short lines
      dEdk1=2*(E_linew(nintep+1,j)-E_linew(1,j))   ! curvature at Gamma (assuming a slop=0)
      dEdk2=E_linew(1,j)+E_linew(2*nintep+1,j)-2*E_linew(nintep+1,j) ! curvature calculated at 2

      if(dabs(dEdk1).lt.5*abs(dEdk2)) then       ! the Gamma curv similar to 2 curv, indeed, slop=0
      
      num_evan=num_evan+1
      E_evanC(num_evan)=dEdk1     
      E_evan(num_evan)=E_linew(1,j)
cccc if E_evanC>0, then evan state is at [-infinity,E_evan]
cccc if E_evane<0, then evan state is at [E_evan,infinity]
      ist_evan(num_evan)=ist_linew(1,j)
      ikpt_evan(num_evan)=ikpt_linew(1,j)
      iGX_evan(num_evan)=0
      iband_evan(num_evan)=j
      endif
      endif
      enddo
cccccccccccccccccccccccccccccccccccc
cccccc second, the evanescent state at the X point
      do j=1,nstw
       numwt=numw(j)
      if(ikpt_linew(numwt,j).eq.nkptw.and.
     & numwt.ge.3*nintep) then    ! X-point, forget about short lines

      dEdk1=2*(E_linew(numwt-nintep,j)-E_linew(numwt,j))   ! curvature at X (assuming a slop=0)
      dEdk2=E_linew(numwt,j)+E_linew(numwt-2*nintep,j)
     &     -2*E_linew(numwt-nintep,j) ! curvature calculated at 2

      if(dabs(dEdk1).lt.5*abs(dEdk2)) then       ! the Gamma curv similar to 2 curv, indeed, slop=0

      num_evan=num_evan+1
      E_evanC(num_evan)=dEdk1     
      E_evan(num_evan)=E_linew(numwt,j)
      ist_evan(num_evan)=ist_linew(numwt,j)
      ikpt_evan(num_evan)=ikpt_linew(numwt,j)
      iGX_evan(num_evan)=1
      iband_evan(num_evan)=j
       endif
       endif
       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc now, find max-min points in the middle 

      do j=1,nstw
      do i=2,numw(j)-1

      if(E_linew(i,j).gt.E_linew(i-1,j).and.
     &   E_linew(i,j).ge.E_linew(i+1,j)) then

      num_evan=num_evan+1
      E_evanC(num_evan)=E_linew(i-1,j)+E_linew(i+1,j)
     &              -2*E_linew(i,j)
      E_evan(num_evan)=E_linew(i,j)
      ist_evan(num_evan)=ist_linew(i,j)
      ikpt_evan(num_evan)=ikpt_linew(i,j)
      iGX_evan(num_evan)=2
      iband_evan(num_evan)=j
      endif

ccccccccccccccccccccccccccccccccccc
        if(E_linew(i,j).le.E_linew(i-1,j).and.
     &   E_linew(i,j).lt.E_linew(i+1,j)) then

      num_evan=num_evan+1
      E_evanC(num_evan)=E_linew(i-1,j)+E_linew(i+1,j)
     &              -2*E_linew(i,j)
      E_evan(num_evan)=E_linew(i,j)
      ist_evan(num_evan)=ist_linew(i,j)
      ikpt_evan(num_evan)=ikpt_linew(i,j)
      iGX_evan(num_evan)=2
      iband_evan(num_evan)=j
      endif
      enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(14,file="E_line_W2K2",position="append")
      write(14,*) "*************************"
      write(14,*) num_evan,"  evanescent states"
       do ii=1,num_evan
       write(14,606) ii,ikpt_evan(ii),ist_evan(ii),iGX_evan(ii),
     & iband_evan(ii), E_evan(ii),E_evanC(ii)
       enddo
       close(14)
606   format(5(i4,1x),2(E18.10,1x))
      stop
       end
