      program wr_interp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Meng Ye, Jnu. 2018
cccccc   This program use interpolation method to calculate the wave function between two near k-point. 
cccccc   It also calculate connect information between the insert k-point
cccccc   By now, the number of process must be equal to nkpt-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      include "mpif.h"

      integer status(MPI_STATUS_SIZE)

      real*8 AL(3,3)
      complex*16 cai,cc,c11,c12,c21,c22

      integer kpt_dis(2)

      real*8, allocatable, dimension (:,:) :: Eband,Eband_new,tmp
      real*8, allocatable, dimension (:) :: err_band,EE
      complex*16, allocatable, dimension (:,:,:,:) :: uc1,uc2
      complex*16, allocatable, dimension (:,:,:) :: cphase
      real*8, allocatable, dimension (:) :: ucR,ucI
      complex*16, allocatable, dimension (:,:) :: S_m,S,S11,S12,S21,S22
      complex*16, allocatable, dimension (:,:) :: hh
      complex*16, allocatable, dimension (:,:) :: coeff1,coeff2

      complex*16, allocatable, dimension (:) :: work
      real*8, allocatable, dimension (:) :: rwork

      character*20 fwr_all,fwr_head

      open(10,file="connect.input")
      read(10,*) nst,nkpt
      read(10,*) iwrnew
      read(10,*) nintp ! The number of points need to be inserted the two nearest k-point
      close(10)

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)      
      inode=inode+1

      nkpt0=nkpt-1
      if(nnodes.gt.nkpt0) then
       if(inode.eq.1) write(6,*) "nnodes is larger than nkpt-1, stop",
     &   nnodes, nkpt
       call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      nkpt_new=nkpt0*nintp+1

      allocate(Eband(nst,nkpt))
      allocate(err_band(nst))
      open(10,file="Ek.tmp")
      do ikpt=1,nkpt
       do j=1,4
       read(10,*)
       enddo
       read(10,*) (err_band(i),i=1,nst)
       read(10,*)
       read(10,*) (Eband(i,ikpt),i=1,nst)
      enddo
      close(10)
      deallocate(err_band)

      n_tmp1=nkpt0/nnodes
      n_tmp2=nkpt0-n_tmp1*nnodes
      if(inode.le.n_tmp2) then
        kpt_dis(1)=(inode-1)*(n_tmp1+1)+1
        kpt_dis(2)=inode*(n_tmp1+1)
      else
        kpt_dis(1)=(inode-1)*n_tmp1+n_tmp2+1
        kpt_dis(2)=inode*n_tmp1+n_tmp2
      endif

      if(iwrnew.eq.0) then
        fwr_head="wr.out"
      else
        fwr_head="wr.new."
      endif
      nst=nst-1

      kpt_dens=kpt_dis(1)
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all=trim(fwr_head)//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)      
      open(11,file=fwr_all,form="unformatted")
      rewind(11)
      read(11) n1,n2,n3,nnodes_t,ispin_i,ispin_f,iw_i,iw_f 
      read(11) AL

      nr=n1*n2*n3
      nr_n=nr/nnodes_t
      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      allocate(uc1(n1,n2,n3,nst))
      allocate(uc2(n1,n2,n3,nst))

      do m=1,nst
      do iproc=1,nnodes_t
       read(11) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc2(i,j,k,m)=dcmplx(ucR(ii),ucI(ii))
       enddo
      enddo !end loop over iproc
      enddo !end loop over m
      close(11)

      sum00=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
      sum00=sum00+cdabs(uc2(i,j,k,1))**2
      enddo
      enddo
      enddo

      allocate(cphase(n1,n2,n3))
      pi=4*datan(1.d0)
      cai=dcmplx(0.d0,1.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1
      cphase(i,j,k)=cdexp(cai*pi*(i-1.d0)/(nkpt0*nintp*n1))
      enddo
      enddo
      enddo

      allocate(S_m(nst,nst))
      allocate(S(nst,nst))
      allocate(S11(nst,nst))
      allocate(S12(nst,nst))
      allocate(S21(nst,nst))
      allocate(S22(nst,nst))
      allocate(coeff1(nst,nst))
      allocate(coeff2(nst,nst))
      allocate(hh(nst,nst))
      allocate(EE(nst))

      lwork=2*nst-1
      allocate(work(lwork))
      allocate(rwork(3*nst-2))
      allocate(Eband_new(nst,nkpt_new))
      Eband_new=0.d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 3000 ikpt=kpt_dis(1),kpt_dis(2)
      fwr_all="wr.interp"//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)      
      open(12,file=fwr_all,form="unformatted")
      rewind(12)
      fwr_all="wr.overlap"//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)      
      open(13,file=fwr_all,form="unformatted")
      rewind(13)

      uc1=uc2

      kpt_dens=ikpt+1
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all=trim(fwr_head)//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)
      open(11,file=fwr_all,form="unformatted")
      rewind(11)
      read(11) n1,n2,n3,nnodes_t,ispin_i,ispin_f,iw_i,iw_f
      read(11) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes_t
      do m=1,nst
      do iproc=1,nnodes_t
       read(11) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc2(i,j,k,m)=dcmplx(ucR(ii),ucI(ii))
       enddo
      enddo !end loop over iproc
      enddo
      close(11)

      do ist1=1,nst
      do ist2=1,nst
        cc=dcmplx(0.d0,0.d0)
        cc11=dcmplx(0.d0,0.d0)
        cc12=dcmplx(0.d0,0.d0)
        cc21=dcmplx(0.d0,0.d0)
        cc22=dcmplx(0.d0,0.d0)
        do k=1,n3
        do j=1,n2
        do i=1,n1
        cc=cc+uc1(i,j,k,ist1)*dconjg(uc2(i,j,k,ist2))
        cc11=cc11+uc1(i,j,k,ist1)*dconjg(uc1(i,j,k,ist2))*cphase(i,j,k)
        cc12=cc12+uc1(i,j,k,ist1)*dconjg(uc2(i,j,k,ist2))*cphase(i,j,k)
        cc21=cc21+uc2(i,j,k,ist1)*dconjg(uc1(i,j,k,ist2))*cphase(i,j,k)
        cc22=cc22+uc2(i,j,k,ist1)*dconjg(uc2(i,j,k,ist2))*cphase(i,j,k)
        enddo
        enddo
        enddo
        S_m(ist1,ist2)=cc/sum00
        S11(ist1,ist2)=cc11/sum00
        S12(ist1,ist2)=cc12/sum00
        S21(ist1,ist2)=cc21/sum00
        S22(ist1,ist2)=cc22/sum00
      enddo
      enddo

      write(13) S_m
      write(13) S11
      write(13) S12
      write(13) S21
      write(13) S22
      close(13)

      S=S_m

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Gram–Schmidt process, make S_m to be an unitary matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do m1=1,nst
        do m2=1,m1-1
          cc=dcmplx(0.d0,0.d0)
          do m3=1,nst
          cc=cc+S_m(m1,m3)*dconjg(S_m(m2,m3))
          enddo
          do m3=1,nst
          S_m(m1,m3)=S_m(m1,m3)-cc*S_m(m2,m3)
          enddo
        enddo
        sum=0.d0
        do m3=1,nst
        sum=sum+cdabs(S_m(m1,m3))**2
        enddo
        sum=dsqrt(1.d0/sum)
        do m3=1,nst
        S_m(m1,m3)=S_m(m1,m3)*sum
        enddo
      enddo

      do iintp=1,nintp-1
        w2=dble(iintp)/dble(nintp)
        w1=1.d0-w2
        do m1=1,nst
        do m2=m1,nst
          cc=dcmplx(0.d0,0.d0)
          do m3=1,nst
            cc=cc+Eband(m3,ikpt+1)*dconjg(S_m(m1,m3))*S_m(m2,m3)
          enddo
          hh(m1,m2)=w2*cc
        enddo
        hh(m1,m1)=hh(m1,m1)+w1*Eband(m1,ikpt)
        enddo

        call zheev('V','U',nst,hh,nst,EE,work,lwork,rwork,info)

        if(info.ne.0) then
        write(6,*) "Something wrong with zheev, stop. ikpt,iintp="
     &           ,ikpt,iintp
        call mpi_abort(MPI_COMM_WORLD,ierr)
        endif

        ikpt_new=(ikpt-1)*nintp+1+iintp
        Eband_new(:,ikpt_new)=EE(:)

        coeff1(:,:)=w1*hh(:,:)
        do ist=1,nst
        do m=1,nst
          cc=dcmplx(0.d0,0.d0)
          do m1=1,nst
          cc=cc+S_m(m1,m)*hh(m1,ist)
          enddo
          coeff2(m,ist)=w2*cc
        enddo
        enddo

C     Some test code!
C         do ist=1,nst
C           sum=0.d0
C           do k=1,n3
C           do j=1,n2
C           do i=1,n1
C           cc=dcmplx(0.d0,0.d0)
C           do m=1,nst
C           cc=cc+coeff1(m,ist)*uc1(i,j,k,m)-coeff2(m,ist)*uc2(i,j,k,m)
C           enddo
C           sum=sum+cdabs(cc)**2
C           enddo
C           enddo
C           enddo
C           write(6,*) "err=",dsqrt(sum/sum00)
C         enddo

        
        do ist=1,nst
        sum=0.d0
        do i=1,nst
        do j=1,nst
        sum=sum+2*dreal(coeff1(i,ist)*dconjg(coeff2(j,ist))*S(i,j))
        enddo
        sum=sum+cdabs(coeff1(i,ist))**2+cdabs(coeff2(i,ist))**2
        enddo
        sum=dsqrt(1.d0/sum)
        coeff1(:,ist)=coeff1(:,ist)*sum
        coeff2(:,ist)=coeff2(:,ist)*sum
        write(12) (coeff1(i,ist),i=1,nst)
        write(12) (coeff2(i,ist),i=1,nst)
        enddo

      enddo

      close(12)

3000  continue

      deallocate(ucR)
      deallocate(ucI)
      deallocate(uc1)
      deallocate(uc2)
      deallocate(cphase)
      deallocate(S_m)
      deallocate(S11)
      deallocate(S12)
      deallocate(S21)
      deallocate(S22)
      deallocate(hh)
      deallocate(EE)
      deallocate(work)
      deallocate(rwork)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      allocate(tmp(nst,nkpt_new))
      call mpi_allreduce(Eband_new,tmp,nst*nkpt_new,MPI_REAL8,
     &  MPI_SUM,MPI_COMM_WORLD,ierr)
      Eband_new=tmp
      deallocate(tmp)

      do i=1,nkpt
        Eband_new(:,(i-1)*nintp+1)=Eband(:,i)
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      deallocate(Eband)

      if(inode.eq.1) then
        open(22,file="Ek.interp")
        rewind(22)
        do ikpt=1,nkpt_new
        write(22,*) "*********************************"
        write(22,*) "*********************************"
        write(22,104) ikpt
        write(22,*) "eigen energies, in eV "
        write(22,103) (Eband_new(i,ikpt), i=1,nst)        
        enddo
        close(22)

        open(23,file="wr.interp",form="unformatted")
        rewind(23)
        write(23) nkpt_new,nst
        open(24,file="wr.overlap",form="unformatted")
        rewind(24)
        write(24) nkpt,nst

        do ikpt=1,nkpt0
        coeff1=dcmplx(0.d0,0.d0)
        coeff2=dcmplx(0.d0,0.d0)
        do ist=1,nst
         coeff1(ist,ist)=dcmplx(1.d0,0.d0)
         write(23) (coeff1(i,ist),i=1,nst)
         write(23) (coeff2(i,ist),i=1,nst)
        enddo
        kpt_dens=ikpt
        ikpt_100=kpt_dens/100
        ikpt_10=(kpt_dens-ikpt_100*100)/10
        ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
        fwr_all="wr.interp"//char(48+ikpt_100)//char(48+ikpt_10)//
     &     char(48+ikpt_1)      
        open(12,file=fwr_all,form="unformatted")
        rewind(12)
        fwr_all="wr.overlap"//char(48+ikpt_100)//char(48+ikpt_10)//
     &     char(48+ikpt_1)      
        open(13,file=fwr_all,form="unformatted")
        rewind(13)
        do i=1,5
          read(13) S
          write(24) S
        enddo
        close(13,status='delete')

        do iintp=1,nintp-1
          do ist=1,nst
          read(12) (coeff1(i,ist),i=1,nst)
          read(12) (coeff2(i,ist),i=1,nst)
          enddo
          do ist=1,nst
          write(23) (coeff1(i,ist),i=1,nst)
          write(23) (coeff2(i,ist),i=1,nst)
          enddo
        enddo
        close(12,status='delete')
        enddo
        coeff1=dcmplx(0.d0,0.d0)
        coeff2=dcmplx(0.d0,0.d0)
        do ist=1,nst
         coeff1(ist,ist)=dcmplx(1.d0,0.d0)
         write(23) (coeff1(i,ist),i=1,nst)
         write(23) (coeff2(i,ist),i=1,nst)
        enddo
        close(23)
        close(24)
      endif

103   format(5(f12.8,1x))
104   format("kpt= ", i5)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      deallocate(coeff1)
      deallocate(coeff2)
      deallocate(S)
      deallocate(Eband_new)
      
      call mpi_finalize(ierr)

      stop
      end
