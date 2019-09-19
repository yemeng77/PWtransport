      subroutine wave_electrode_3kpts(ikpt_st1,ikpt_st2,
     & ikpt_st3,i_st1,i_st2,i_st3,x_st1,x_st2,x_st3,
     & uc,n1,n2,n3,fileh,iflag,cphase,phase,nintep)

ccc A. Garcia-Lekue, (August 2008) ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  this subroutine performs a 3-point interpolation
ccc  for running  electrode states
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


      implicit double precision (a-h,o-z)
      complex*16, allocatable, dimension (:,:,:) :: uc1,uc2,uc3
      complex*16 uc(n1,n2,n3)
      complex*16 cc1,cc2,cc
      complex*16 cphase(n1,n2,n3)
      real*8 phase(n1,n2,n3)
      character*7 fileh
***********************************************************

      allocate(uc1(n1,n2,n3))
      allocate(uc2(n1,n2,n3))
      allocate(uc3(n1,n2,n3))

      call wave_electrode_interp(ikpt_st1,i_st1,uc1,n1,n2,n3,
     &     fileh,nintep)
      call wave_electrode_interp(ikpt_st2,i_st2,uc2,n1,n2,n3,
     &     fileh,nintep)
      call wave_electrode_interp(ikpt_st3,i_st3,uc3,n1,n2,n3,
     &     fileh,nintep)

cccccccccccccccccccccccccccccccccccccccccccccccccc
      k12=ikpt_st2-ikpt_st1     ! k12 should be 1
      k13=ikpt_st3-ikpt_st1     ! k13 should be 2
      if(iflag.eq.1) then
      uc2=uc2*(cphase**k12)
      uc3=uc3*(cphase**k13)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1
      cc1=cc1+uc1(i,j,k)*dconjg(uc2(i,j,k))
      cc2=cc2+uc1(i,j,k)*dconjg(uc3(i,j,k))
      enddo
      enddo
      enddo
      cc1=cc1/abs(cc1)
      cc2=cc2/abs(cc2)
      x1=x_st1
      cc1=cc1*x_st2
      cc2=cc2*x_st3


      sum0=0.d0
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc(i,j,k)=uc1(i,j,k)*x1+cc1*uc2(i,j,k)+cc2*uc3(i,j,k)
      sum=sum+cdabs(uc(i,j,k))**2
      sum0=sum0+cdabs(uc1(i,j,k))**2
      enddo
      enddo
      enddo
      fact=dsqrt(sum0/sum)
      uc=uc*fact      ! normalize uc,this is needed, later, it is assumed uc is normalized
cccccccccccccccccccccccccccccccccccccccccccc

      if(iflag.eq.1) then
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc(i,j,k)=uc(i,j,k)*cdexp(dcmplx(0.d0,
     &  ((k12*x_st2)+(k13*x_st3))*phase(i,j,k)))    ! This is correct 1+ak=k12*x_st2+k13*x_st3
      enddo
      enddo
      enddo
      endif

      deallocate(uc1)
      deallocate(uc2)
      deallocate(uc3)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
