       subroutine occup(dE,itypeFermi,E_st,totNel,nkpt,
     & occ,E_NSC,TS,workr_n,Ef,islda,dV_bias)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'
      
      real*8 Ef,occ(mst,nkpt,islda),E_st(mst,nkpt,islda)
      real*8 occ_L(mst,nkpt,islda),occ_R(mst,nkpt,islda)
      complex*16 workr_n(mr_n)
      complex*16, allocatable, dimension(:,:) :: workr2_n
      real*8, allocatable, dimension (:) :: rho_n_tmp
      real*8 weight_L(mst,nkpt,islda),weight_L_tmp(mst,nkpt,islda)
      real*8 weight_R(mst,nkpt,islda),weight_R_tmp(mst,nkpt,islda)
      real*8 tmp_real2(nblock_band_mx)
      real*8 tmp_real(mst)
		

******************************************


      nblock=mx
      allocate(workr2_n(mr_n,nblock_band_mx))
	     
	     
      n=totNel/2+0.1d0
      if(n.gt.mx) then
      write(6,*) "n.gt.mx, stop", n,mx,mst
      stop
      endif


**********************************************************************
***** Calculate left and right weights ******************************************
**********************************************************************


*** get the charge density

*** calculate Left and Right contributions to the charge density, including mask function


	weight_L = 0.d0
	weight_R = 0.d0

c	tmp_real = 0.d0
c	tmp_real2 = 0.d0


	write(82,*) kpt_slda_dis(1)
	write(82,*) kpt_slda_dis(2)

	num_kpt_proc=kpt_slda_dis(2)-kpt_slda_dis(1)+1
	do 102 iislda=1,islda

	do kpt=1,nkpt

	 if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then      ! the big division of work

      call gen_G_comp(kpt,0)
      call fftprep_comp(n1,n2,n3)

      if(num_kpt_proc.gt.1) then
      call ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)
      endif

      call d3fft_comp_block(ug_n_bp,workr2_n,-1,kpt,nblock_band_mx)
      do m=1, nblock_band_mx
	 mm = band_dis(1) + m - 1
         mp = (nblock_band_mx*icolor_b) + m
	 do i=1, nr_n
	    weight_L(mp,kpt,iislda) = weight_L(mp,kpt,iislda)+
     &           + vmask_nL(i)*cdabs(workr2_n(i,m))**2
	    weight_R(mp,kpt,iislda) = weight_R(mp,kpt,iislda)+
     &           + (1.d0-vmask_nL(i))*cdabs(workr2_n(i,m))**2
	 enddo ! loop over nr_n
      enddo !loop over m states

ccc SUM OVER MPI_COMM_B1

	m2 = nblock_band_mx*icolor_b+1

       call mpi_allreduce(weight_L(m2,kpt,iislda),
     &      tmp_real2,nblock_band_mx,MPI_REAL8,
     &		MPI_SUM,MPI_COMM_B1,ierr)
      do m=1,nblock_band_mx
      mp = (nblock_band_mx*icolor_b) + m
      weight_L(mp,kpt,iislda) = tmp_real2(m)*vol/(n1*n2*n3)
      enddo


	m2 = nblock_band_mx*icolor_b+1

      call mpi_allreduce(weight_R(m2,kpt,iislda),
     &      tmp_real2,nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_B1,ierr)
      do m=1, nblock_band_mx
      mp = (nblock_band_mx*icolor_b) + m
      weight_R(mp,kpt,iislda) = tmp_real2(m)*vol/(n1*n2*n3)
      enddo

c	do m=1,nblock_band_mx
c	mp = (nblock_band_mx*icolor_b) + m
c	write(66,*) mp,weight_R(mp,kpt,iislda),
c     &	weight_L(mp,kpt,iislda)+weight_R(mp,kpt,iislda)
c	enddo


ccc SUM OVER MPI_COMM_B2

c      call mpi_allreduce(weight_L(1,kpt,iislda),tmp_real,
c     &      mx,MPI_REAL8,MPI_SUM,MPI_COMM_B2,ierr)
c	do m=1,mx
c	weight_L(m,kpt,iislda) = tmp_real(m)
c	enddo


c      call mpi_allreduce(weight_R(1,kpt,iislda),tmp_real,
c     &      mx,MPI_REAL8,MPI_SUM,MPI_COMM_B2,ierr)
c	do m=1,mx
c	weight_R(m,kpt,iislda) = tmp_real(m)
c	enddo

       endif     ! kpt_slda_dis(1),(2)
      enddo      ! kpt
102   continue   ! iislda


!      real*8 weight_L(mst,nkpt,islda),weight_L_tmp(mst,nkpt,islda)
!      real*8 weight_R(mst,nkpt,islda),weight_R_tmp(mst,nkpt,islda)

      do iislda=1,islda
         call mpi_allreduce(weight_L,weight_L_tmp,mst*nkpt*islda,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_B2,ierr)
         call mpi_allreduce(weight_L_tmp,weight_L,mst*nkpt*islda,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_K2,ierr)
         call mpi_allreduce(weight_R,weight_R_tmp,mst*nkpt*islda,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_B2,ierr)
         call mpi_allreduce(weight_R_tmp,weight_R,mst*nkpt*islda,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_K2,ierr)
      enddo


**********************************************************************
**********************************************************************
**********************************************************************



      Ef=E_st(n,1,1)
      if(mx.gt.n) Ef=(E_st(n,1,1)+E_st(n+1,1,1))/2

****** first, using three point Emax, Emin, then try (Emax+Emin)/2
****** only at the final steps, using the following derivatives.

      Emax=-10000.d0
      Emin=10000.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      if(E_st(m,kpt,iislda).gt.Emax) Emax=E_st(m,kpt,iislda)
      if(E_st(m,kpt,iislda).lt.Emin) Emin=E_st(m,kpt,iislda)
      enddo
      enddo
      enddo
      dV_bias_au = dV_bias/27.211396d0

      Emax=Emax+20*dE+dV_bias_au/2
      Emin=Emin-20*dE-dV_bias_au/2



      s1=2*totNel
      s2=0.d0

      do it=1,100

c      if(mod(it,2).eq.1) then 
      if(it.gt.3.and.s2.gt.0.d0.and.s1.lt.2*totNel) then 
      fac=(s1-totNel)/(s1-s2)
      Ef=Emin*fac+Emax*(1.d0-fac)
      else
      Ef=(Emax+Emin)/2
      endif

      if(dabs(s1-totNel).lt.1.D-12*totNel) then
      Ef=Emax
      goto 90
      endif

      if(dabs(s2-totNel).lt.1.D-12*totNel) then
      Ef=Emin
      goto 90
      endif



ccc Here now include the left and right occupations

	dV_bias_au = dV_bias/27.211396d0

      Ef_L = Ef - (dV_bias_au/2)	
      Ef_R = Ef + (dV_bias_au/2)	

c      write(6,*) 'Ef,Ef_L,Ef_R,dV_bias',Ef,Ef_L,Ef_R,dV_bias


      s_L=0.d0
      s_R=0.d0
      s_tot=0.d0
      s_test=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx

      y_L=(E_st(m,kpt,iislda)-Ef_L)/dE
      call Fermi(f_occ_L,S_occ_L,y_L,itypeFermi)
      s_L=s_L+f_occ_L*weighkpt(kpt)*weight_L(m,kpt,iislda)

      y_R=(E_st(m,kpt,iislda)-Ef_R)/dE
      call Fermi(f_occ_R,S_occ_R,y_R,itypeFermi)
      s_R=s_R+f_occ_R*weighkpt(kpt)*weight_R(m,kpt,iislda)

      y=(E_st(m,kpt,iislda)-Ef)/dE
      call Fermi(f_occ,S_occ,y,itypeFermi)
      s_test=s_test+f_occ*weighkpt(kpt)

      enddo
      enddo
      enddo

      s_tot = s_L + s_R
      s_tot=s_tot*2/islda
      s_test=s_test*2/islda
      if(s_tot.ge.totNel) then
      Emax=Ef
      s1=s_tot
      else
      Emin=Ef
      s2=s_tot
      endif

      if(inode.eq.1) then
      write(6,*) "test", it,s_tot,s_test,Ef
      endif

      enddo   ! do it=1,100


90    continue

500   format("s_L,s_R,s_tot,s_test", 
     &    3x,4(f7.2,1x))


      occ_charge=s_tot


ccccccc TS is the kT*entropy term for free energy, F=E-TS
ccccccc This term is needed in order for the total
ccccccc energy (free energy) be the local minimont at
ccccccc the Schrodinger equation.
cccccccccccc The current formula is only correct for Fermi-Dirac distribution

      dV_bias_au = dV_bias/27.211396d0

      Ef_L = Ef - (dV_bias_au/2)	
      Ef_R = Ef + (dV_bias_au/2)	

      TS=0.d0
      TS_L=0.d0
      TS_R=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx

      y_L=(E_st(m,kpt,iislda)-Ef_L)/dE
      call Fermi(yy_L,S_occ_L,y_L,itypeFermi)

      occ_L(m,kpt,iislda)=2.d0/islda*yy_L*weighkpt(kpt)

      TS_L=TS_L+2.d0/islda*weighkpt(kpt)*
     &	dE*S_occ_L*weight_L(m,kpt,iislda)

      y_R=(E_st(m,kpt,iislda)-Ef_R)/dE
      call Fermi(yy_R,S_occ_R,y_R,itypeFermi)

      occ_R(m,kpt,iislda)=2.d0/islda*yy_R*weighkpt(kpt)
      TS_R=TS_R+2.d0/islda*weighkpt(kpt)*
     &	dE*S_occ_R*weight_R(m,kpt,iislda)


      enddo
      enddo
      enddo
      TS = TS_L + TS_R
***********************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccc
***********************************************


*** get the charge density
      num_kpt_proc=kpt_slda_dis(2)-kpt_slda_dis(1)+1
      do 100 iislda=1,islda
      do i=1,nr_n
      rho_n(i,iislda)=0.d0
      enddo

      do kpt=1,nkpt

       if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then      ! the big division of work

      call gen_G_comp(kpt,0)
      call fftprep_comp(n1,n2,n3)

      if(num_kpt_proc.gt.1) then
      call ugIOBP(ug_n_bp,kpt,2,0,iislda,-1,nkpt,islda)
      endif

      call d3fft_comp_block(ug_n_bp,workr2_n,-1,kpt,nblock_band_mx)
      do m=1, nblock_band_mx
         mm = band_dis(1) + m - 1
  	 mp = (nblock_band_mx*icolor_b) + m
         do i=1, nr_n
            rho_n(i,iislda) = rho_n(i,iislda) 
     &   + (occ_L(mm,kpt,iislda)*vmask_nL(i)
     &                        *cdabs(workr2_n(i,m))**2)
     &  + (occ_R(mm,kpt,iislda)*(1.d0-vmask_nL(i))
     &                        *cdabs(workr2_n(i,m))**2)
         enddo
      enddo


       endif     ! kpt_slda_dis(1),(2)
      enddo
100   continue

      deallocate(workr2_n)
cccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(rho_n_tmp(nr_n))
      do iislda=1,islda
         call mpi_allreduce(rho_n(1,iislda),rho_n_tmp,nr_n,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_B2,ierr)
         call mpi_allreduce(rho_n_tmp,rho_n(1,iislda),nr_n,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_K2,ierr)         
      enddo
      deallocate(rho_n_tmp)
ccccccccccccccccccccccccccccccccccccccccccccccccc


      E_NSC=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
c      E_NSC=E_NSC+occ(m,kpt,iislda)*E_st(m,kpt,iislda)

      E_NSC=E_NSC+(occ_L(m,kpt,iislda)*weight_L(m,kpt,iislda)
     &  +occ_R(m,kpt,iislda)*weight_R(m,kpt,iislda))
     &		*E_st(m,kpt,iislda)
      occ(m,kpt,iislda)=occ_L(m,kpt,iislda)*weight_L(m,kpt,iislda)
     &   +occ_R(m,kpt,iislda)*weight_R(m,kpt,iislda)
      enddo
      enddo
      enddo

*********************************************

      return
      end


      
      subroutine Fermi(f_occ,S_occ,x,itypeFermi)
      implicit double precision(a-h,o-z)

cccccccccc   According to: 
cccc  (1) M. Methfessel, A.T. Paxton, Phys. Rev. B 40, 3616 (1989)
cccc  (2) G. Kresse, J. Furthmuller, Comp. Mat. Sci. 6, 15 (1996). 

      real*8 h(0:30)


      y=x

      if(itypeFermi.eq.1) then       !  Fermi-Diract smearing
      if(y.gt.100.d0) y=100.d0
      if(y.lt.-100.d0) y=-100.d0
      f_occ=1.d0/(dexp(y)+1.d0)
      S_occ=-((dabs(1.d0-f_occ)+1.D-30)*dlog(dabs(1.d0-f_occ)+
     & 1.D-30)+(dabs(f_occ)+1.D-30)*dlog(dabs(f_occ)+1.D-30))
      return
      endif

ccc Using M. Methfessel, A.T. Paxton, P.R.B, 40, 3616 (1989). 
      if(itypeFermi.ge.20.and.itypeFermi.lt.30) then 
      n=itypeFermi-20    ! the order of S_n function
      if(y.gt.9.d0) y=9.d0 
      if(y.lt.-9.d0) y=-9.d0
      f_occ=0.5d0*(1-derf(y))
      fg=exp(-y**2)

      h(0)=1.d0
      h(1)=2*y
      do i=1,2*n
      h(i+1)=2*y*h(i)-2*i*h(i-1)
      enddo

      a=1.d0/dsqrt(4*datan(1.d0))

      do i=1,n
      a=-a/(i*4)
      f_occ=f_occ+a*h(2*i-1)*fg
      enddo
      S_occ=0.5d0*a*h(2*n)*fg

      if(y.lt.-8.99) then
      f_occ=1.d0
      S_occ=0.d0
      endif

      return
      endif
ccccccccccccccccccccccccccccccccccccc

      write(6,*) "itypeFermi not defined, stop",itypeFermi
      stop
      end
      
      
      
