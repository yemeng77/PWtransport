      program test
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      real*8 AL(3,3)
      complex*16, allocatable, dimension (:) :: uc_tmp
      complex*16, allocatable, dimension (:,:,:,:) :: uc1
      real*8, allocatable, dimension (:) :: ucR,ucI

      open(10,file="wr.K21.052",form="unformatted")
      rewind(10)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
c      read(10) AL(1,1),AL(2,1),AL(3,1) 	      
c      read(10) AL(1,2),AL(2,2),AL(3,2) 	      
c      read(10) AL(1,3),AL(2,3),AL(3,3) 	      
      read(10) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes

      nst=iw_f-iw_i+1
      write(6,*) 'nst',nst

      write(6,*) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      write(6,*) nr_n

      allocate(uc_tmp(nr_n))
      allocate(uc1(n1,n2,n3,nst))
      allocate(ucR(nr_n))
      allocate(ucI(nr_n))

      write(6,*) 'nr_n',nr_n

      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f
	
      do iproc=1,nnodes

      read(10) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
      do i=1,nr_n
      uc_tmp(i)=dcmplx(ucR(i),ucI(i))
      enddo
	if ((m.eq.1).and.(iproc.eq.1)) then
	do i=1,10
	write(65,*) i,ucR(i),ucI(i)
	enddo
	endif

      do ii=1,nr_n
      jj=ii+(iproc-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      uc1(i,j,k,m)=uc_tmp(ii)
      enddo

      enddo !end loop over iproc

      enddo !end loop over m
      enddo !end loop over iislda

c	deallocate(uc_tmp)
c	deallocate(ucR)
c	deallocate(ucI)

       close(10)


       write(6,*) 'Input state to be represented'
        read(5,*) jst

        open(11,file='graph.plot1D')
        rewind(11)
        do i=1,n1
        sum=0.d0
        do j=1,n2
        do k=1,n3
        sum=sum+abs(uc1(i,j,k,jst))**2
        enddo
        enddo
        sum=sum/(n2*n3)
        write(11,333) i,sum
        enddo
        close(11)
333     format(i6,2x,E18.5)


	deallocate(uc1)

	stop
	end

