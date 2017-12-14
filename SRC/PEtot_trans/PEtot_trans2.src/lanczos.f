      subroutine lanczos_eigen(ilocal,nline,
     &  vr,workr_n,kpt,Eref,AL,mxc
     &  eigen,ug_n)
****************************************
cc     Written by Meng Ye, December 13, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     Use Lanczos algorithm. T=V^dagger * (H-Eref)**2 * V,
cc     where V is a ng_n*mxc matrix V with orthonormal columns, 
cc     T is a mxc*mxc real symmetric tridiagonal matrix.
****************************************
*  Arguments
*  =========
*
*  D       the mxc diagonal elements of the tridiagonal matrix T
*          
*  E       the (mxc-1) subdiagonal elements of the tridiagonal matrix T
*
*  V       all elements of matrix V^dagger
*
*  =====================================================================

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

       real*8 AL(3,3)
c       complex*16 workr_n(mg_nx)
       complex*16 workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang

       real*8 D(mxc),E(mxc-1),Z(mxc,mxc)
       real*8, allocatable, dimension(:,:) :: work
       complex*16, allocatable, dimension (:,:) :: V
       complex*16, allocatable, dimension (:) :: v_temp,vh_temp,w_temp

       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek

       real*8 ran1

       ng_n=ngtotnod(inode,kpt)

       allocate(V(mxc,ng_n))
       allocate(v_temp(ng_n))
       allocate(vh_temp(ng_n))
       allocate(w_temp(ng_n))

************************************************
**** generate random inital v_temp
**** v_temp = v_temp / ||v_temp||
************************************************
       iranm=-2291-inode*3651
       do i=1,ng_n
       x = ran1(iranm)
       y = ran1(iranm)
       v_temp(i)= dcmplx(x-0.5d0,y-0.5d0)
       enddo

       s=0.d0
       do i=1,ng_n
       s=s+cdabs(v_temp(i))**2
       enddo
       
       call global_sumr(s)
       s=1/dsqrt(s*vol)

       do i=1,ng_n
       v_temp(i)=s*v_temp(i)
       enddo
      

      do 4000 nint=1,mxc

       if nint.ne.1 then
************************************************
**** beta = ||w_temp||
**** v_temp = w_temp / beta
**** beta is the nint-1 th element of E
**** v_temp is the nint th row of V
************************************************
       beta=0.d0

       do i=1,ng_n
       beta=beta+cdabs(w_temp(i))**2
       enddo
       call global_sumr(beta)
       beta=dsqrt(beta*vol)
       if (beta.eq.0.d0) exit
       do i=1,ng_n
       v_temp(i)=w_temp(i)/beta
       enddo
       E(nint-1)=beta
       endif

       do i=1,ng_n
       V(nint,i)=v_temp(i)
       enddo


************************************************
**** w_temp = (H-Eref)**2 * v_temp
************************************************
       call Hpsi_comp(v_temp,vh_temp,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       vh_temp(i)=vh_temp(i)-Eref*v_temp(i)
       enddo

       call Hpsi_comp(vh_temp,w_temp,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       w_temp(i)=w_temp(i)-Eref*vh_temp(i)
       enddo


************************************************
**** alpha = w_temp^dagger * v_temp
**** w_temp = w_temp - alpha * v_temp, when nint=1
**** w_temp = w_temp - alpha * v_temp - beta * v_temp(old), when nint>1
************************************************
       alpha=0.d0
       
       do i=1,ng_n
       alpha=alpha+dreal(dconjg(w_temp(i))*v_temp(i))
       enddo
       call global_sumr(alpha)
       alpha=dsqrt(alpha*vol)

       if (nint.eq.1) then
       do i=1,ng_n
       w_temp(i)=w_temp(i)-alpha*V(nint,i)
       enddo
       else
       do i=1,ng_n
       w_temp(i)=w_temp(i)-alpha*V(nint,i)-beta*V(nint-1,i)
       enddo
       endif

4000  continue


      deallocate(v_temp)
      deallocate(vh_temp)
      deallocate(w_temp)
***********************************************

************************************************
**** V=V*
************************************************
      do i=1,mxc
      do j=1,ng_n
      V(i,j)=dconjg(V(i,j))
      enddo
      enddo

      call detev('V',mxc,D,E,Z,mxc,work,info)

      deallocate(V)

      return
      end

