      subroutine wave_system(uc,n1,n2,n3,nnodes,mstate,num_wr_dis)

      implicit double precision (a-h,o-z)

      complex*16, allocatable, dimension (:) :: uc_tmp
      complex*16 uc(n1,n2,n3,10)
      real*8 sum_st(10),sum_tmp(10)

      integer num_wr_dis(2)

      complex*16 cc
***********************************************************
      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(uc_tmp(nr_n))

      do iist=1,num_wr_dis(1)-1
      do iread=1,nnodes
      read(21)
      enddo
      enddo

      ist=0
      do iist=num_wr_dis(1),num_wr_dis(2)
      ist=ist+1
      do iread=1,nnodes
      read(21) (uc_tmp(i),i=1,nr_n)
       do ii=1,nr_n
       jj=ii+(iread-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc(i,j,k,ist)=uc_tmp(ii)
       enddo
      enddo
      enddo
      
      do iist=num_wr_dis(2)+1,mstate
      do iread=1,nnodes
      read(21)
      enddo
      enddo

      deallocate(uc_tmp)

cccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      nst=num_wr_dis(2)-num_wr_dis(1)+1
      do ist=1,nst
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
      sum=sum+cdabs(uc(i,j,k,ist))**2
      enddo
      enddo
      enddo
      sum_tmp(ist)=sum
      sum=1.d0/dsqrt(sum)
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc(i,j,k,ist)=sum*uc(i,j,k,ist)
      enddo
      enddo
      enddo
      enddo
c      write(6,*) "***the |uc|^2 of input Well(system) states,mst=",
c     &    mstate
c      write(6,888) (sum_tmp(ist),ist=1,mstate)
c888   format(12(E6.1,1x))
  
ccccccccccccccccccccccccccccccccc

      iplot=0
      if(iplot.eq.1) then
      open(10,file="graph.plot.Wellst")
      rewind(10)
      do i=1,n1
      do ll=1,mstate
      sum=0.d0
      do k=1,n3
      do j=1,n2
      sum=sum+cdabs(uc(i,j,k,ll))**2
      enddo
      enddo
      sum_st(ll)=sum/(n2*n3)
      enddo
      write(10,202) i,(sum_st(ll),ll=1,mstate)
      enddo
      close(10)
202   format(i5,2x,11(E10.4,1x))
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
