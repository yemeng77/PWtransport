      subroutine wave_electrode_interp(ikpt_st1,ikpt_st2,nst,
     & vec1,vec2,x_st1,x_st2,
     & uc,m1,m2,m3,fileh)

ccc A. Garcia-Lekue, (August 2008) ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  this subroutine performs a 3-point interpolation
ccc  for running  electrode states
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


      implicit double precision (a-h,o-z)

      complex*16, allocatable, dimension (:,:,:) :: uc1,uc2
      complex*16, allocatable, dimension (:) :: uc_tmp
      real*8, allocatable, dimension (:) :: ucR,ucI
      complex*16 uc(m1,m2,m3)
      complex*16 vec1(nst),vec2(nst)
      complex*16 cc
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

      allocate(uc_tmp(nr_n))
      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      allocate(uc1(n1,n2,n3))

      uc1=dcmplx(0.d0,0.d0)

      sum0=0.d0

      do iislda=ispin_i,ispin_f
      do ist=1,nst
      do iread=1,nnodes
       read(10) (ucR(i),i=1,nr_n),(ucI(i),i=1,nr_n)
       do i=1,nr_n
       uc_tmp(i)=dcmplx(ucR(i),ucI(i))
       enddo
       do ii=1,nr_n
        jj=ii+(iread-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        uc1(i,j,k)=uc1(i,j,k)+vec1(ist)*uc_tmp(ii)
        if(ist.eq.1) then
          sum0=sum0+cdabs(uc_tmp(ii))**2
        endif
       enddo
      enddo
      enddo
      enddo
      close(10)

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

      allocate(uc2(n1,n2,n3))

      uc2=dcmplx(0.d0,0.d0)

      do iislda=ispin_i,ispin_f
      do ist=1,nst
      do iread=1,nnodes
       read(10) (ucR(i),i=1,nr_n),(ucI(i),i=1,nr_n)
       do i=1,nr_n
       uc_tmp(i)=dcmplx(ucR(i),ucI(i))
       enddo
       do ii=1,nr_n
        jj=ii+(iread-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        uc2(i,j,k)=uc2(i,j,k)+vec2(ist)*uc_tmp(ii)
       enddo
      enddo
      enddo
      enddo
      close(10)

      deallocate(uc_tmp)
cccccccccccccccccccccccccccccccccccccccccccccccccc

      sum0=0.d0
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc(i,j,k)=x_st1*uc1(i,j,k)+x_st2*uc2(i,j,k)
      sum=sum+cdabs(uc(i,j,k))**2
      sum0=sum0+cdabs(uc1(i,j,k))**2
      enddo
      enddo
      enddo
      fact=dsqrt(sum0/sum)
      uc=uc*fact      ! normalize uc,this is needed, later, it is assumed uc is normalized
cccccccccccccccccccccccccccccccccccccccccccc
      deallocate(uc1)
      deallocate(uc2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
