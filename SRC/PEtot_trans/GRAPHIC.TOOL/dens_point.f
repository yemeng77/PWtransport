      program dens_point
cccccccccccccccccccccccccccccccccccccccccc
cccc This is a serial code to generate any point in the box from 
cccc the charge density on a regular grid. It uses the Fourier components, and FFT interpolation
cc     Written by Lin-Wang Wang
******************************************

      implicit double precision (a-h,o-z)


      parameter (lwork=50000)

      real*8 work(lwork)

      real*8, allocatable, dimension(:,:,:) :: vr
      real*8, allocatable, dimension(:,:,:) :: vr_r,vr_i
      real*8, allocatable, dimension(:) :: vr_tmp

      real*8 AL(3,3)

      complex*16 cai,cai2pi,cc,csum
      character*50 filename

       write(6,*) "input the name of dens (vr) file to be plot"
       read(5,*) filename
       open(11,file=filename,form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) AL
       write(6,*) AL(1,1),AL(2,1),AL(3,1)
       write(6,*) AL(1,2),AL(2,2),AL(3,2)
       write(6,*) AL(1,3),AL(2,3),AL(3,3)
       write(6,*) "n1,n2,n3=", n1,n2,n3

       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr(n1,n2,n3))
       allocate(vr_r(n1,n2,n3))
       allocate(vr_i(n1,n2,n3))

       do iread=1,nnodes
       read(11) vr_tmp

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       close(11)



       vr_r=vr
       vr_i=0.d0


c      naux=0                    ! the dcft3 will dynamically allocate memory
c      scale=1.d0/(n1*n2*n3)
c      call dcft3(vrc,n1,n1*n2,vrc,n1,n1*n2,n1,n2,n3,
c     &     -1,scale,aux,naux)

       call  cfft(n1,n2,n3,vr_r,vr_i,work,lwork,1)

ccccccccccccccccccccccccccccccccccccccccccccccc
         cai2pi=2*dcmplx(0.d0,4*atan(1.d0))
c         x1=
c         x2=
c         x3=
1000     continue
         write(6,*) "input x1,x2,x3"
         read(5,*) x1,x2,x3
         csum=dcmplx(0.d0,0.d0)
         do k=1,n3
            ak=0.d0        ! set the center point (n3/2) to zero
            fack=1.d0
            if(k.le.(n3+1)/2) ak=k-1
            if(k.gt.n3/2+1) ak=k-n3-1
            if(k.eq.n3/2.and.mod(n3,2).eq.0) fack=0.d0
         do j=1,n2
            aj=0.d0        ! set the center point (n3/2) to zero
            facj=1.d0
            if(j.le.(n2+1)/2) aj=j-1
            if(j.gt.n2/2+1) aj=j-n2-1
            if(j.eq.n2/2.and.mod(n2,2).eq.0) facj=0.d0
         do i=1,n1
            ai=0.d0        ! set the center point (n3/2) to zero
            faci=1.d0
            if(i.le.(n1+1)/2) ai=i-1
            if(i.gt.n1/2+1) ai=i-n1-1
            if(i.eq.n1/2.and.mod(n1,2).eq.0) faci=0.d0

          cc=faci*facj*fack*dcmplx(vr_r(i,j,k),vr_i(i,j,k))
          cc=cc/(n1*n2*n3)
          csum=csum+cc*exp(cai2pi*(x1*ai+x2*aj+x3*ak))
          enddo
          enddo
          enddo

         write(6,*) "csum=", csum
 
         i1=x1*n1+1.D-8+n1
         i1=mod(i1,n1)+1
         j1=x2*n2+1.D-8+n2
         j1=mod(j1,n2)+1
         k1=x3*n3+1.D-8+n3
         k1=mod(k1,n3)+1

         write(6,*) "vr=",vr(i1,j1,k1) 

         goto 1000


         stop
         end


