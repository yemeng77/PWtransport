      subroutine  d3fft_comp_block(fg_n,fr_n,isign,kpt,nblock)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


***************************************************
****
****    Warning ! fr_n(mr_n) must be statically allocated
****    before the call of this subroutine using:
****    real*8 fr_n(1)
****    pointer(fr_n_p,fr_n)
****    call shpalloc(fr_n_p,mr_n,errcode,-1)
****
****
***************************************************
****    warning !!!
****  for  fr_n --> fg_n, the original fr_n will be destroyed
****  this FFT includes a fg_n(i), i.e the smooth Ecut mask
****  As a result, the wavefunction in real space fr_n(i) are
****  no longer normalized and orthogonal. Use other FFT to
****  get the "real" (without wg_n) real space wavefunction.
***************************************************

***************************************************
****  d3fft_real(fg,fr,isign)
****  isign=1: fr_n --> fg_n   (forward fft)
****  isign=-1: fg_n --> fr_n  (inverse fft)
****  normalization: \sum_i fr_n(i)^2 *vol/(n1*n2*n3)=1
****  \sum_i fg_n(i)^2  * 2 * vol = 1
***************************************************


      use fft_data
      use load_data
      use data

      implicit none

      include "param.escan_real"

      include "mpif.h"
      
      
      complex*16 fg_n(mg_nx,nblock)
      complex*16 fr_n(mr_n,nblock)   
      complex*16 cdum

      integer i,ig,indepg,n1_inv,n2_inv,indepg_d,isign,
     &     jjnode_dum,kpt,ierr,nblock,j
      


      if(isign.eq.-1) then

c     put into x col format for fft in workr_n

         do j=1, nblock
            do i = 1,ncol(inode)*n1
               fr_n(i,j) = dcmplx(0.d0,0.d0)
            enddo
         enddo

         do j=1, nblock
            do ig = 1, ngtotnod(inode,kpt)
               indepg=(jjcol(n2p_n(ig),n3p_n(ig))-1)*n1+n1p_n(ig)
               fr_n(indepg,j) = fg_n(ig,j)*wg_n(ig,kpt)
            enddo
         enddo

c         call mpi_barrier(MPI_COMM_K,ierr)
 
!         do j=1, nblock
!            call invcpfft_comp(fr_n(1,j),n1,n2,n3)
!         enddo
         call invcpfft_comp_block(fr_n,nblock,n1,n2,n3)
         
         return
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc end of inversed FFT 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      if(isign.eq.1) then

!         do j=1, nblock
!            call fwdcpfft_comp(fr_n(1,j),n1,n2,n3)
!         enddo
         call fwdcpfft_comp_block(fr_n,nblock,n1,n2,n3)

c     put back into load balanced g vector distribution
c     

         do j=1, nblock
            do ig = 1, ngtotnod(inode,kpt)
               indepg=(jjcol(n2p_n(ig),n3p_n(ig))-1)*n1+n1p_n(ig)
               fg_n(ig,j)=wg_n(ig,kpt)*fr_n(indepg,j)
            enddo
         enddo

c         call mpi_barrier(MPI_COMM_K,ierr)
         return
      endif


      return
      end


