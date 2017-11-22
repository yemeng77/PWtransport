      subroutine wave_electrode_3kpts(ikpt_st1,ikpt_st2,
     & ikpt_st3,i_st1,i_st2,i_st3,x_st1,x_st2,x_st3,
     & uc,m1,m2,m3,fileh,iflag,
     & cphase,phase)

ccc A. Garcia-Lekue, (August 2008) ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  this subroutine performs a 3-point interpolation
ccc  for running  electrode states
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


      implicit double precision (a-h,o-z)
      complex*16, allocatable, dimension (:,:,:) :: uc1,uc2
      complex*16, allocatable, dimension (:,:,:) :: uc3
      complex*16, allocatable, dimension (:) :: uc_tmp1
      real*8, allocatable, dimension (:) :: ucR1,ucI1
      real*8, allocatable, dimension (:) :: ucR2,ucI2
      real*8, allocatable, dimension (:) :: ucR3,ucI3
      complex*16 uc(m1,m2,m3)
      complex*16 cc1,cc2,cc
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

c	write(6,*) '1st filename in wave_electrode',
c     &          filename

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
       
      ikpt_100=ikpt_st2/100
      ikpt_10=(ikpt_st2-ikpt_100*100)/10
      ikpt_1=ikpt_st2-ikpt_100*100-ikpt_10*10
      filename=fileh//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)

c      write(6,*) '2nd filename in wave_electrode4_run',
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
      allocate(uc2(n1,n2,n3))
      allocate(ucR2(nr_n))
      allocate(ucI2(nr_n))

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

      deallocate(uc_tmp1)
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc hau nere kosetza

c      if (ikpt_st3.gt.101)goto 307

      ikpt_100=ikpt_st3/100
      ikpt_10=(ikpt_st3-ikpt_100*100)/10
      ikpt_1=ikpt_st3-ikpt_100*100-ikpt_10*10
      filename=fileh//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)


c      write(6,*) '3rd filename in wave_electrode4_run',
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
      allocate(uc3(n1,n2,n3))
      allocate(ucR3(nr_n))
      allocate(ucI3(nr_n))

       do iislda=ispin_i,ispin_f
      do ist=1,i_st3
      do iread=1,nnodes
       read(10) (ucR2(i),i=1,nr_n),(ucI2(i),i=1,nr_n)
       do i=1,nr_n
       uc_tmp1(i)=dcmplx(ucR2(i),ucI2(i))
       enddo


      if(ist.eq.i_st3) then
      do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      uc3(i,j,k)=uc_tmp1(ii)
      enddo
      endif
      enddo
      enddo
      enddo
      close(10)


      deallocate(uc_tmp1)

c  307 continue
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
