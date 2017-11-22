      subroutine invcpfft_comp_block(psi,nblock,nr1,nr2,nr3)

c     Written by A. Canning (CRAY-EPFL) 25th July 94  
c     
c     output = psi(z,y,x)  wavefunction in real space z columns 
c     input  = psi(x,y,z)  x columns of g vectors in load balalanced distr.
c     
c     computes inverse fft specifically for CP algo 
c     ie taking sphere and going to cylinder then two (one for gamma)
c     slices then a complete cube. 
c     see fftprep.f for more details
c     fftprep must be called once before this routine
c     for setups.
c     
c     last revised ADV 7/7/95 to parametrise the no. of 
c     chunks on each PE
c     AMC revised 21/6/97 hardwired for gamma point complex to real fft
c
c   invcfft_comp is slower than fwdcfft because the filling of zero 
cc  in fwdcfft, there is no need to fill the zero. This simple step of 
ccc filling the zero is expensive because the large number of point. 
ccc This filling to zero is more than 1/2 of the ccfft, this is a bit crazy. 

      use load_data
      use fft_data
      use data
      
      implicit none

      include "mpif.h"
c     
c     
      real*8,allocatable,dimension(:) :: worknr1,worknr2,worknr3
      integer,allocatable,dimension(:) :: idum_start1,idum_start2
      integer,allocatable,dimension(:) :: idum_start3,idum_start4
      complex*16,allocatable,dimension(:) :: combuf1,combuf2
      complex*16,allocatable,dimension(:,:) :: psiy

      real*8 time_000,time_001

      complex*16 psi(mr_n,nblock)

      integer ireq(nnodes)

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

c     scalars used

      integer i,ib,ic,idum,ii,ilocadd,isc,isign,itar,itaradd,
     c     itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1,nr2,nr3,
     c     nr3u,iloc_dum,ierr, idum2, iblock,idum1,idum3,k
      integer nworknr1,nworknr2,nworknr3


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
      isign = -1
c     
c     do the FFT's in place on the non-zero columns 
c     of psi
c     so this is the x dir FFT
c     

cccc do iblock the outer most loop, so the psi after ccfft for a given iblock might
cccc still in the cache, that will save a lot of time for the re-shaffling of the data

      do iblock=1, nblock

         do jcol = 1,ncol(inode)
            ilocadd = 1+(jcol-1)*nr1
            
            call system_ccfft(0,isign,nr1,1.0d0,psi(ilocadd,iblock),
     &           psi(ilocadd,iblock),tabnr1in,worknr1,0,ntabnr1,
     &           nworknr1)

         enddo

c     
c     transpose to the y slice  y col mode  in psiy
c     

      idum1 = 1
      idum2 = 1
      do i = 1,nnodes
      idum_start1(i)=idum1-1
      idum_start2(i)=idum2-1
      idum2=idum2+ivunpn1(i)
         
         if(ivpacn1(i).ne.0) then 
         if(i.ne.inode) then
         idum3=idum_start1(i)*nblock+ivpacn1(i)*(iblock-1)+1
            do j = 1,ivpacn1(i)
            combuf1(idum3) = psi(ivpac1(idum1),iblock)
            idum1=idum1+1
            idum3=idum3+1
            enddo
         else
             idum1=idum1+ivpacn1(i)
        endif
        endif
      enddo   
      enddo   ! iblock


      do k = 1,nnodes-1
       i=mod(inode+k-1,nnodes)+1
       if(ivpacn1(i).ne.0) then
       call mpi_isend(combuf1(idum_start1(i)*nblock+1),
     &        ivpacn1(i)*nblock,
     &        mpi_double_complex,i-1,inode,mpi_comm_k,ireq(i),ierr)
       endif
      enddo


ccccccccccccccccccccccccccccccccccccccccccccccc

cccccc We might want to combine preparing combuf1 with send, and receive with combuf2. 
cccccc basically, after combuf for each node is ready, then send. 
cccccc and after combuf from each node is received, then reassign. 
cccccc they both take time.
      
      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivunpn1(i).ne.0) then
       call mpi_recv(combuf2(idum_start2(i)*nblock+1),
     &        ivunpn1(i)*nblock,
     &        mpi_double_complex,i-1,i,mpi_comm_k,mpistatus,ierr)
       endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccc
cccc Now, we need to make sure combuf1 can be rewritten. This is wait for the mpi_isend
cccc This is needed, so we can write to combuf1 again
      do i = 1, nnodes
        if(ivpacn1(i).ne.0.and.i.ne.inode) then
        call mpi_wait(ireq(i), mpistatus, ierr)
        endif
      end do

cccccccccccccccccccccccccccccccccccccccc
      do 4000 iblock=1, nblock   ! big nblock 

         do ii = 1,ncoly*2*nr2
            psiy(ii,iblock) = 0.0   ! this simple step is expensive
         enddo

      idum2=idum_start1(inode)+1
      idum1=idum_start2(inode)+1
      do j=1,ivunpn1(inode)
       psiy(ivunp1(idum1),iblock)=psi(ivpac1(idum2),iblock)
       idum2=idum2+1
       idum1=idum1+1
      enddo

ccccccccccccccccc

      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivunpn1(i).ne.0) then
        idum3=idum_start2(i)*nblock+ivunpn1(i)*(iblock-1)+1
        idum1=idum_start2(i)+1
        do j=1,ivunpn1(i)
               psiy(ivunp1(idum1),iblock) = combuf2(idum3)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
        enddo
       endif
       enddo


c     
c     do FFT on the y direction on the two (one for gamma) slices 
c     this is the y dir FFT
c     each PE should have ncoly columns for the FFT
c     

         do i = 1,2*ncoly
            ilocadd = 1+(i-1)*nr2
            call system_ccfft(0,isign,nr2,1.0d0,psiy(ilocadd,iblock),
     &           psiy(ilocadd,iblock),tabnr2in,
     &           worknr2,0,ntabnr2,nworknr2)

         enddo
c     
c    transpose slice using mpi calls
c     zero values of psi which we do not put into
c

      
      idum1 = 1
      idum2 = 1
      do i = 1,nnodes
       idum_start3(i)=idum1-1
       idum_start4(i)=idum2-1
       idum2=idum2+ivunpn2(i)

       if(ivpacn2(i).ne.0) then
       if(i.ne.inode) then
       idum3=idum_start3(i)*nblock+ivpacn2(i)*(iblock-1)+1
            do j = 1,ivpacn2(i)
               combuf1(idum3) = psiy(ivpac2(idum1),iblock)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
            enddo
        else
             idum1=idum1 + ivpacn2(i)
        endif
        endif
      enddo

4000  continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k = 1,nnodes-1
         i=mod(inode+k-1,nnodes)+1
         if(ivpacn2(i).ne.0) then
         call mpi_isend(combuf1(idum_start3(i)*nblock+1),
     &      ivpacn2(i)*nblock,
     &        mpi_double_complex,i-1,inode,mpi_comm_k,ireq(i),ierr)
         endif
      enddo


      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivunpn2(i).ne.0) then
       idum1=idum_start4(i)*nblock+1
       call mpi_recv(combuf2(idum1),
     &        ivunpn2(i)*nblock,
     &        mpi_double_complex,i-1,i,mpi_comm_k,mpistatus,ierr)
       endif
       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccc
      do 5000 iblock=1, nblock     ! big do block

         do i = 1,ncolz
            idum = (i-1)*nr3
            do j = mgz+1,(nr3-mgz)
               psi(idum+j,iblock)=0.0     ! this simple step is expensive 
            enddo
         enddo

      idum1=idum_start3(inode)+1
      idum2=idum_start4(inode)+1
      do j=1,ivunpn2(inode)
      psi(ivunp2(idum2),iblock)=psiy(ivpac2(idum1),iblock)
      idum2=idum2+1
      idum1=idum1+1
      enddo
ccccccccccccccccccccccccccccccccccccccccccc
      do k = 1,nnodes-1
       i=mod(inode-k-1+nnodes,nnodes)+1
       if(ivunpn2(i).ne.0) then
       idum3=idum_start4(i)*nblock+ivunpn2(i)*(iblock-1)+1
       idum1=idum_start4(i)+1
            do j = 1,ivunpn2(i)
               psi(ivunp2(idum1),iblock) = combuf2(idum3)
               idum1 = idum1 + 1
               idum3 = idum3 + 1
            enddo
        endif
        enddo

c     
c     now do FFT's on z direction so each proc has nr1*nr2/nnodes 
c     columns of height nr3 complex numbers (real for gamma point)
c     
         do i = 1,ncolz
            ilocadd = 1+(i-1)*nr3
            call system_ccfft(0,isign,nr3,1.0d0,psi(ilocadd,iblock),
     &           psi(ilocadd,iblock),tabnr3in,worknr3,0,
     &           ntabnr3,nworknr3)

         enddo
   
5000   continue

ccccccccccccccccccccccccccccccccccccccccccccccc
      do i = 1, nnodes
         if(ivunpn2(i).ne.0.and.i.ne.inode) then
         call mpi_wait(ireq(i), mpistatus, ierr)
         endif
      end do
ccccccccccccccccccccccccccccccccccccccccccccccc



      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)
      deallocate(psiy)
      deallocate(combuf1)
      deallocate(combuf2)
      deallocate(idum_start1)
      deallocate(idum_start2)
      deallocate(idum_start3)
      deallocate(idum_start4)

      end
      


