      subroutine wave_system_ym(ucL,ucR,nnposit1,nnposit2,
     &  n1w,n1,n2,n3,nnodes,num_wr_dis,nst)

      implicit double precision (a-h,o-z)

      integer num_wr_dis(2)
      complex*16, allocatable, dimension (:) :: uc_tmp
      complex*16 ucL(n1w,n2,n3,nst)   ! uc does not contain uc^*
      complex*16 ucR(n1w,n2,n3,nst)   ! uc does not contain uc^*
      complex*16 uc(n1,n2,n3)
***********************************************************
      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(uc_tmp(nr_n))
      
      do inum=1,num_wr_dis(1)-1
       do iread=1,nnodes
        read(21)
       enddo
      enddo

      ist=0
      do inum=num_wr_dis(1),num_wr_dis(2)
       ist=ist+1
       do iread=1,nnodes
        read(21) (uc_tmp(i),i=1,nr_n)
        do ii=1,nr_n
          jj=ii+(iread-1)*nr_n
          i=(jj-1)/(n2*n3)+1
          j=(jj-1-(i-1)*n2*n3)/n3+1
          k=jj-(i-1)*n2*n3-(j-1)*n3
          uc(i,j,k)=uc_tmp(ii)
        enddo
       enddo

       sum=0.d0
       do k=1,n3
       do j=1,n2
       do i=1,n1
         sum=sum+cdabs(uc(i,j,k))**2
       enddo
       enddo
       enddo
       sum=1.d0/dsqrt(sum)
       do k=1,n3
       do j=1,n2
       do i=1,n1
         uc(i,j,k)=sum*uc(i,j,k)
       enddo
       enddo
       enddo

       do k=1,n3
       do j=1,n2
       do i=1,n1w
         ucL(i,j,k,ist)=uc(i+nnposit1,j,k)
         ucR(i,j,k,ist)=uc(i+nnposit2,j,k)
       enddo
       enddo
       enddo
      enddo
      deallocate(uc_tmp)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
