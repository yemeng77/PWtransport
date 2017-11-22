      subroutine fwdcpfft_comp_block(psi,nblock,nr1,nr2,nr3)
c------------------------------------------------------
c
c Written by A. Canning (CRAY-EPFL) 25th July 94  
c 
c
c output = psi(x,y,z)  wavefunction in fourier space load balanced
c                      as small sphere columns in x direction
c input = psi(z,y,x) real space grid with each PE having consecutive 
c                    sets of z columns.
c
c computes forward fft specifically for CP algo
c ie taking cube and going to two slices (one slice for gamma)
c then  to a sphere.  See fftprep.f for more details
c fftprep must be called once before this routine
c for setups.
       
      use fft_data
      use load_data
      use data

      implicit none

      include "mpif.h"

      real*8,allocatable,dimension(:) :: worknr1,worknr2,worknr3
      integer,allocatable,dimension(:):: idum_start1,idum_start2
      integer,allocatable,dimension(:):: idum_start3,idum_start4
      complex*16,allocatable,dimension(:) :: combuf1,combuf2
      complex*16,allocatable,dimension(:,:) :: psiy

      
      complex*16 psi(mr_n,nblock)

      integer nblock

      integer inode, nnodes, inode_tot, nnodes_tot
      integer nnodes_k, nnodes_b, inode_k, inode_b
      integer icolor_k, icolor_b, ikey_k, ikey_b, icolor, ikey
      integer num_group_k, num_group_b
      integer MPI_COMM_K1, MPI_COMM_K2, MPI_COMM_B1, MPI_COMM_B2
      integer MPI_COMM_K, MPI_COMM_N

      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &     nnodes_k, nnodes_b, inode_k, inode_b,
     &     icolor_k, icolor_b, ikey_k, ikey_b, icolor, ikey,
     &     num_group_k, num_group_b, 
     &     MPI_COMM_K1, MPI_COMM_K2, MPI_COMM_B1, MPI_COMM_B2,
     &     MPI_COMM_K, MPI_COMM_N


      integer mpistatus(mpi_status_size)
      integer ireq(nnodes)

c scalars used

      integer i,ib,ic,idum,ii,iloc,ilocadd,isc,isign,itar,itaradd,
     c        itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1,nr2,nr3,
     c        nr3u,iloc_dum,ierr,idum2,iblock,idum1,idum3,k
      integer nworknr1,nworknr2,nworknr3

      real*8 fac,s

c
      nworknr1 = 20000+2.28*nr1x
      nworknr2 = 20000+2.28*nr2x
      nworknr3 = 20000+2.28*nr3x

      allocate(worknr1(nworknr1))
      allocate(worknr2(nworknr2))
      allocate(worknr3(nworknr3))
      allocate(psiy(mr_n,nblock))
      allocate(combuf1(mr_n*nblock))
      allocate(combuf2(mr_n*nblock))
      allocate(idum_start1(nnodes))
      allocate(idum_start2(nnodes))
      allocate(idum_start3(nnodes))
      allocate(idum_start4(nnodes))


c
      isign = 1
c
c full FFT on the cube psi this is on the z dir
c layout is nz,ny,nx  colulmns of height nr3 complex
c reslt not put directly into psi (problem for scfft)
c due to different sizes of input and output
c

      do iblock=1, nblock

         do i = 1,ncolz
            ilocadd = 1+(i-1)*nr3
            call system_ccfft(0,isign,nr3,1.0d0,psi(ilocadd,iblock),
     &           psi(ilocadd,iblock),tabnr3fw,
     &           worknr3,0,ntabnr3,nworknr3)


         enddo

c
c now transpose nz,ny,nx to ny,nz,nx in the two slice mode
c into psiy. Each PE will have ncoly columns in psiy 
c


      idum1 = 1
      idum2 = 1
      do i = 1,nnodes
        idum_start1(i)=idum1-1
        idum_start2(i)=idum2-1
        idum2=idum2+ivpacn2(i)
        
         if(ivunpn2(i).ne.0) then
         if(i.ne.inode) then
         idum3=idum_start1(i)*nblock+ivunpn2(i)*(iblock-1)+1
            do j = 1,ivunpn2(i)
               combuf1(idum3) = psi(ivunp2(idum1),iblock)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
            enddo
         else
            idum1=idum1+ivunpn2(i)
         endif
         endif
       enddo
       enddo   ! iblock

cccccccccccccccccccccccccccccccccccccccccccccccccc
      do k = 1,nnodes-1
       i=mod(inode+k-1,nnodes)+1
       if(ivunpn2(i).ne.0) then
       call mpi_isend(combuf1(idum_start1(i)*nblock+1),
     &        ivunpn2(i)*nblock,
     &     mpi_double_complex,i-1,inode,mpi_comm_k,ireq(i),ierr)
       endif
      enddo

      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivpacn2(i).ne.0) then
       call mpi_recv(combuf2(idum_start2(i)*nblock+1),
     &        ivpacn2(i)*nblock,
     &        mpi_double_complex,i-1,i,mpi_comm_k,mpistatus,ierr)
        endif
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc
      do i = 1, nnodes
         if(ivunpn2(i).ne.0.and.i.ne.inode) then
         call mpi_wait(ireq(i), mpistatus, ierr)
         endif
      end do
ccccccccccccccccccccccccccccccccccccccc

      do 4000 iblock=1,nblock

      idum2=idum_start1(inode)+1
      idum1=idum_start2(inode)+1
      do j=1,ivpacn2(inode)
       psiy(ivpac2(idum1),iblock)=psi(ivunp2(idum2),iblock)
       idum2=idum2+1
       idum1=idum1+1
      enddo

ccccccccccccccccccccccccccccccccccccccccccccc
      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivpacn2(i).ne.0) then
       idum3=idum_start2(i)*nblock+ivpacn2(i)*(iblock-1)+1
       idum1=idum_start2(i)+1
            do j = 1,ivpacn2(i)
               psiy(ivpac2(idum1),iblock) = combuf2(idum3)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
           enddo
        endif
      enddo


c
c now do FFT's on the two slice
c
         do i = 1,2*ncoly
            ilocadd = 1+(i-1)*nr2
            call system_ccfft(0,isign,nr2,1.0d0,psiy(ilocadd,iblock),
     &           psiy(ilocadd,iblock),tabnr2fw,worknr2,
     &           0,ntabnr2,nworknr2)

         enddo

cccccccccccccccccccccccccccccccccccccccccccccccccc

c 
c now transpose back to format of program into psi 
c ie into x columns load balanced 
c 

      idum1 = 1
      idum2 = 1
      do i = 1,nnodes
      idum_start3(i)=idum1-1
      idum_start4(i)=idum2-1
      idum2=idum2+ivpacn1(i)
         
         if(ivunpn1(i).ne.0) then
         if(i.ne.inode) then
         idum3=idum_start3(i)*nblock+ivunpn1(i)*(iblock-1)+1
            do j = 1,ivunpn1(i)
               combuf1(idum3) = psiy(ivunp1(idum1),iblock)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
            enddo
          else
           idum1 = idum1 + ivunpn1(i)
          endif
          endif
       enddo

4000   continue
ccccccccccccccccccccccccccccccccccccc

      do k = 1,nnodes-1
        i=mod(inode+k-1,nnodes)+1
         if(ivunpn1(i).ne.0) then
         call mpi_isend(combuf1(idum_start3(i)*nblock+1),
     &        ivunpn1(i)*nblock,
     &        mpi_double_complex,i-1,inode,mpi_comm_k,ireq(i),ierr)
         endif
      enddo

      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivpacn1(i).ne.0) then
       call mpi_recv(combuf2(idum_start4(i)*nblock+1),
     &        ivpacn1(i)*nblock,
     &        mpi_double_complex,i-1,i,mpi_comm_k,mpistatus,ierr)
        endif
       enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      fac=1.d0/dfloat(nr1*nr2*nr3)

      do 5000 iblock=1,nblock    ! big do iblock

      idum1=idum_start3(inode)+1
      idum2=idum_start4(inode)+1
      do j=1,ivpacn1(inode)
       psi(ivpac1(idum2),iblock)=psiy(ivunp1(idum1),iblock)
       idum2=idum2+1
       idum1=idum1+1
      enddo
cccccccccccccccccccccccccccccc

      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivpacn1(i).ne.0) then
       idum3=idum_start4(i)*nblock+ivpacn1(i)*(iblock-1)+1
       idum1=idum_start4(i)+1
            do j = 1,ivpacn1(i)
               psi(ivpac1(idum1),iblock) = combuf2(idum3)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
            enddo
       endif
       enddo

c
c for psi in the x direction 
c
         do jcol = 1,ncol(inode)
            ilocadd = 1+(jcol-1)*nr1
            call system_ccfft(0,isign,nr1,1.0d0,psi(ilocadd,iblock),
     &           psi(ilocadd,iblock),tabnr1fw,worknr1,
     &           0,ntabnr1,nworknr1)

         enddo                  


         do i = 1,ncol(inode)*nr1
            psi(i,iblock) = psi(i,iblock)*fac
         enddo

5000    continue
ccccccccccccccccccccccccccccccc

      do i = 1, nnodes
         if(ivunpn1(i).ne.0.and.i.ne.inode) then
         call mpi_wait(ireq(i), mpistatus, ierr)
         endif
      end do
ccccccccccccccccccccccccccccccc


      deallocate(idum_start1)
      deallocate(idum_start2)
      deallocate(idum_start3)
      deallocate(idum_start4)
      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)
      deallocate(psiy)
      deallocate(combuf1)
      deallocate(combuf2)


      end  

