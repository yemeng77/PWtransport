      program overlap_mpi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Lin-Wang Wang, Feb. 2011
cccccc   this program calculate the connectivity one after another k points for the band structure
cccccc   It also makes the necessary rotations for the degenerated states
cccccc modified by Meng Ye, Feb. 2018
cccccc by now, the number of process must be equal to the nkpt-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      include "mpif.h"

      integer status(MPI_STATUS_SIZE)

      real*8 AL(3,3),tmp_real
      complex*16 cc,cc_tmp
      complex*16, allocatable, dimension (:,:,:,:) :: uc
      complex*16, allocatable, dimension (:,:,:) :: uc_tmp
      real*8,allocatable,dimension(:) :: ucR,ucI
      complex*16, allocatable, dimension(:,:) :: cc_coeff

      character*20 fwr_all


      open(10,file="connect.input")
      read(10,*) nst,nkpt
      close(10)
      nst=nst-1
      nkpt0=nkpt-1

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)      
      inode=inode+1

      if(nnodes.ne.nkpt0) then
       if(inode.eq.1) write(6,*) "nnodes is not equal to nkpt-1",
     &   nnodes, nkpt
       call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      ikpt=inode
      kpt_dens=ikpt
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all="wr.new."//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)
      
      open(11,file=fwr_all,form="unformatted")
      rewind(11)
      read(11) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
      read(11) AL

      nr=n1*n2*n3
      nr_n=nr/nnodes 

      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      allocate(uc(n1,n2,n3,nst))
      allocate(uc_tmp(n1,n2,n3))
      allocate(cc_coeff(nst,nst))

      kpt_dens=ikpt+1
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all="wr.new."//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)

      if(inode.eq.1) then
        open(10,file="overlap.matrix",form="unformatted")
        rewind(10)
        write(10) ispin_i,ispin_f
        write(10) nst,nkpt
      endif

      open(12,file=fwr_all,form="unformatted")
      rewind(12)
      read(12) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(12) AL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f-1
      do iproc=1,nnodes
       read(12) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc(i,j,k,m)=dcmplx(ucR(ii),ucI(ii))
       enddo
      enddo !end loop over iproc
      enddo !end loop over m

      do m=iw_i,iw_f-1
      do iproc=1,nnodes
       read(11) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc_tmp(i,j,k)=dcmplx(ucR(ii),ucI(ii))
       enddo
      enddo !end loop over iproc

      do ist=1,nst
          cc=dcmplx(0.d0,0.d0)
          do k=1,n3
          do j=1,n2
          do i=1,n1
          cc=cc+uc(i,j,k,ist)*dconjg(uc_tmp(i,j,k))
          enddo
          enddo
          enddo
          cc_coeff(ist,m)=cc
      enddo
      enddo !end loop over m

      if(inode.ne.1) then
        call mpi_send(cc_coeff,nst*nst,MPI_DOUBLE_COMPLEX,0,100+iislda,
     &       MPI_COMM_WORLD,ierr)
      else
        write(10) iislda
        write(10) inode
        write(6,*) inode
        do ist2=1,nst
        write(10) (cc_coeff(ist1,ist2),ist1=1,nst)
        enddo
        do ikpt=2,nkpt0
          call mpi_recv(cc_coeff,nst*nst,MPI_DOUBLE_COMPLEX,ikpt-1,
     &       100+iislda,MPI_COMM_WORLD,status,ierr)
          write(10) ikpt
          write(6,*) ikpt
          do ist2=1,nst
          write(10) (cc_coeff(ist1,ist2),ist1=1,nst)
          enddo
        enddo
      endif
      enddo !end loop over iislda

      close(11)
      close(12)

      if(inode.eq.1) then
      close(10)
      endif

      deallocate(ucR)
      deallocate(ucI)
      deallocate(uc)
      deallocate(uc_tmp)
      deallocate(cc_coeff)

      call mpi_finalize(ierr)

      stop
      end