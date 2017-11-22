      subroutine wave_electrode(ikpt_st1,ikpt_st2,
     & i_st1,i_st2,x_st1,x_st2,uc,m1,m2,m3,fileh,iflag,
     & cphase,phase)

      implicit double precision (a-h,o-z)
      complex*16, allocatable, dimension (:,:,:) :: uc1,uc2
      complex*16, allocatable, dimension (:) :: uc_tmp1
      real*8, allocatable, dimension (:) :: ucR1,ucI1
      real*8, allocatable, dimension (:) :: ucR2,ucI2
      complex*16 uc(m1,m2,m3)
      complex*16 cc1,cc2,cc3,cc
      complex*16 cphase(m1,m2,m3)
      real*8 phase(m1,m2,m3)
      character*7 fileh
      character*20 filename
***********************************************************
      ikpt_100=ikpt_st1/100
      ikpt_10=(ikpt_st1-ikpt_100*100)/10
      ikpt_1=ikpt_st1-ikpt_100*100-ikpt_10*10
      filename=fileh//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)


c       write(6,*) '1st filename in wave_electrode4',
c     &          filename
c	write(6,*) ikpt_st1,ikpt_st2
c	write(6,*) x_st1,x_st2
c	write(6,*) iflag


      open(10,file=filename,form="unformatted")
      rewind(10)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(10) AL
      if(n1.ne.m1.or.n2.ne.m2.or.n3.ne.m3) then
      write(6,*) "n1,n2,n3.ne.m1,m2,m3,stop",n1,n2,n3,m1,m2,m3
      stop
      endif
      nr=n1*n2*n3
      nr_n=nr/nnodes

      allocate(uc_tmp1(nr_n))
      allocate(ucR1(nr_n))
      allocate(ucI1(nr_n))
      allocate(uc1(n1,n2,n3))

       do iislda=ispin_i,ispin_f
       do ist=1,i_st1
       do iread=1,nnodes

       if(ist.eq.i_st1) then
       read(10) (ucR1(i),i=1,nr_n),(ucI1(i),i=1,nr_n)
       do i=1,nr_n
       uc_tmp1(i)=dcmplx(ucR1(i),ucI1(i))
       enddo
       else
       read(10)       ! jump over this record, much faster. 
       endif

       if(ist.eq.i_st1) then
       do ii=1,nr_n
       jj=ii+(iread-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc1(i,j,k)=uc_tmp1(ii)	
       enddo
       endif
       enddo
       enddo
       enddo
       close(10)
       deallocate(uc_tmp1)
cccccccccccccccccccccccccccccccccccccccccccccccccc
       
***********************************************************
      ikpt_100=ikpt_st2/100
      ikpt_10=(ikpt_st2-ikpt_100*100)/10
      ikpt_1=ikpt_st2-ikpt_100*100-ikpt_10*10
      filename=fileh//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)


c	write(6,*) '2nd filename in wave_electrode4_run',
c     &          filename


      open(10,file=filename,form="unformatted")
      rewind(10)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(10) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes
      if(n1.ne.m1.or.n2.ne.m2.or.n3.ne.m3) then
      write(6,*) "n1,n2,n3.ne.m1,m2,m3,stop",n1,n2,n3,m1,m2,m3
      stop
      endif

      allocate(uc_tmp1(nr_n))
      allocate(ucR2(nr_n))
      allocate(ucI2(nr_n))
      allocate(uc2(n1,n2,n3))

       do iislda=ispin_i,ispin_f
       do ist=1,i_st2
       do iread=1,nnodes
       read(10) (ucR2(i),i=1,nr_n),(ucI2(i),i=1,nr_n)
       do i=1,nr_n
       uc_tmp1(i)=dcmplx(ucR2(i),ucI2(i))
       enddo

       if(ist.eq.i_st2) then
       do ii=1,nr_n
       jj=ii+(iread-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc2(i,j,k)=uc_tmp1(ii)
       enddo
       endif
       enddo
       enddo
       enddo
       close(10)

cccccccccccccccccccccccccccccccccccccccc
c***output  the density
       open(11,file="graph.int_elec")
       rewind(11)
       do i=1,n1
       sum=0.d0
       do k=1,n3
       do j=1,n2
       sum=sum+abs(uc2(i,j,k))**2
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
      deallocate(uc_tmp1)

       cc1=dcmplx(0.d0,0.d0)
       do k=1,n3
       do j=1,n2
       do i=1,n1
       cc1=cc1+uc1(i,j,k)*dconjg(uc2(i,j,k))
       enddo
       enddo
       enddo
       cc1=cc1/abs(cc1)
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
