      subroutine write_wl(n1,n2,n3,n1w,nnodes,uc,AL,istate)

      implicit double precision (a-h,o-z)

      real*8 uc(n1w,n2,n3)
      real*8 AL(3,3)
      complex*16 cc
      complex*16, allocatable, dimension (:) :: ur_tmp

      nd=5
      pi=4*datan(1.d0)

      cc=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
        if(abs(uc(i,j,k)).gt.abs(cc)) cc=uc(i,j,k)
      enddo
      enddo
      enddo
      cc=dconjg(cc)/abs(cc)

      ii1=istate/10
      ii2=istate-ii1*10

      open(27,file="wr_test."//char(ii1+48)//char(ii2+48),
     &      form="unformatted")
      rewind(27)
      write(27) n1,n2,n3,nnodes
      write(27) AL

      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(ur_tmp(nr_n))

      do iread=1,nnodes
         do ii=1,nr_n
         jj=ii+(iread-1)*nr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         ur_tmp(ii)=dcmplx(0.d0,0.d0)
         if(i.ge.n1-n1w+1) then
           iw=i-(n1-n1w)
           fac=1.d0
           if(iw.gt.n1w+1-nd) fac=(1.d0-dsin((iw-n1w-1)*pi/2.d0/nd))/2
           if(iw.lt.nd+1) fac=(1.d0+dsin((iw-1)*pi/2.d0/nd))/2
           ur_tmp(ii)=(uc(iw,j,k)*cc)*fac
         endif
         if(i.le.nd) then
           iw=i+n1w
           fac=(1.d0-dsin((iw-n1w-1)*pi/2.d0/nd))/2
           ur_tmp(ii)=(uc(iw-n1w,j,k)*cc)*fac
         endif
         if(i.ge.n1-n1w-nd+2.and.i.le.n1-n1w) then
           iw=i-n1+n1w
           fac=(1.d0+dsin((iw-1)*pi/2.d0/nd))/2
           ur_tmp(ii)=(uc(n1w+iw,j,k)*cc)*fac
         endif
         enddo
         write(27) ur_tmp
      enddo
      close(27)

      deallocate(ur_tmp)

      end subroutine