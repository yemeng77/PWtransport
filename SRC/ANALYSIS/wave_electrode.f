      subroutine wave_electrode(ikpt_st1,ikpt_st2,
     & i_st1,i_st2,x_st1,x_st2,uc,n1,n2,n3,fileh,iflag,
     & cphase,phase,nintep)

      implicit double precision (a-h,o-z)
      complex*16, allocatable, dimension (:,:,:) :: uc1,uc2
      complex*16 uc(n1,n2,n3)
      complex*16 cc1,cc2,cc3,cc
      complex*16 cphase(n1,n2,n3)
      real*8 phase(n1,n2,n3)
      character*7 fileh
***********************************************************

      allocate(uc1(n1,n2,n3))
      allocate(uc2(n1,n2,n3))

      call wave_electrode_interp(ikpt_st1,i_st1,uc1,n1,n2,n3,
     &     fileh,nintep)
      call wave_electrode_interp(ikpt_st2,i_st2,uc2,n1,n2,n3,
     &     fileh,nintep)

cccccccccccccccccccccccccccccccccccccccc
c***output  the density
       open(11,file="graph.int_elec")
       rewind(11)
       do i=1,n1
       sum=0.d0
       do k=1,n3
       do j=1,n2
       sum=sum+cdabs(uc2(i,j,k))**2
       enddo
       enddo
       sum=sum/(n2*n3)
c        write(11,*) i,dreal(uc2(i,2,3)/uc2(1,2,3)),
c     &          dreal(uc2(i,16,20)/uc2(1,2,3)),
c     &  abs(uc2(i,16,20)/uc2(1,2,3)),sum
        write(11,*) i,sum
          enddo
       close(11)
333    format(i6,2x,E18.5)
ccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccc
       if(iflag.eq.1) then
       uc2=uc2*cphase
       endif
cccccccccccccccccccccccccccccccccccccccccccccccccc

       cc1=dcmplx(0.d0,0.d0)
       do k=1,n3
       do j=1,n2
       do i=1,n1
       cc1=cc1+uc1(i,j,k)*dconjg(uc2(i,j,k))
       enddo
       enddo
       enddo
       cc1=cc1/cdabs(cc1)
       x1=x_st1
       cc1=cc1*x_st2


       do k=1,n3
       do j=1,n2
       do i=1,n1
       uc(i,j,k)=uc1(i,j,k)*x1+cc1*uc2(i,j,k)
       enddo
       enddo
       enddo

       if(iflag.eq.1) then
       do k=1,n3
       do j=1,n2
       do i=1,n1
       uc(i,j,k)=uc(i,j,k)*cdexp(dcmplx(0.d0,x_st2*phase(i,j,k)))
       enddo
       enddo
       enddo
       endif

       deallocate(uc1)
       deallocate(uc2)
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
