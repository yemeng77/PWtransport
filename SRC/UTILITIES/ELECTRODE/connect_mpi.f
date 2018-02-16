      program connect
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
      complex*16 cc
      complex*16, allocatable, dimension (:,:,:,:) :: uc
      complex*16, allocatable, dimension (:,:,:) :: uc_tmp,cphase
      integer,allocatable, dimension(:,:) :: iconnect
      real*8,allocatable,dimension(:,:) :: connectX
      real*8,allocatable,dimension(:) :: ucR,ucI 
      real*8, allocatable,dimension(:,:) :: Eband,E_band,err_band
      integer, allocatable,dimension(:,:) :: ist_band,ikpt_band,
     &   iflag_band
      integer ind_g2s(100,1000),ind_s2g(1000),ind_s2g2(1000),nst_g(100)
      integer nkpt_band(1000),mgl(1000)
      complex*16 cai
      complex*16, allocatable, dimension(:,:) :: cc_coeff,cc_tmp
      complex*16, allocatable, dimension(:,:) :: rotate,rotate_tmp 

      character*20 fwr_all


      open(10,file="connect.input")
      read(10,*) nst,nkpt
      close(10)
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

      allocate(Eband(nst,nkpt))
      allocate(err_band(nst,nkpt))
      allocate(E_band(3*nst,nkpt))     ! nband can be 3*nst
      allocate(ist_band(3*nst,nkpt))
      allocate(ikpt_band(3*nst,nkpt))
      allocate(iflag_band(nst,nkpt))

      if(inode.eq.1) write(6,*) "copy report to Ek.tmp"

      open(10,file="Ek.tmp")
      do ikpt=1,nkpt
      do j=1,4
       read(10,*)
      enddo
       read(10,*) (err_band(i,ikpt),i=1,nst)
       read(10,*)
       read(10,*) (Eband(i,ikpt),i=1,nst)
      enddo
      close(10)

      allocate(iconnect(nst,nkpt))
      allocate(connectX(nst,nkpt))

      ikpt=inode
      kpt_dens=ikpt
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all="wr.out"//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)
      
      open(11,file=fwr_all,form="unformatted")
      rewind(11)
      read(11) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
      read(11) AL
      if(inode.eq.1) then
      open(10,file="wr.new.001",form="unformatted")
      rewind(10)        
      write(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
      write(10) AL
      write(6,*) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      endif

      nr=n1*n2*n3
      nr_n=nr/nnodes 

      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      allocate(uc(n1,n2,n3,nst))
      allocate(uc_tmp(n1,n2,n3))
      allocate(cphase(n1,n2,n3))

ccccccccccccccccccccccccccccc
      pi=4*datan(1.d0)
      cai=dcmplx(0.d0,1.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1
      cphase(i,j,k)=cdexp(cai*pi*(i-1.d0)/(nkpt0*n1))
      enddo
      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccc
      kpt_dens=ikpt+1
      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all="wr.out"//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)

      open(12,file=fwr_all,form="unformatted")
      rewind(12)
      read(12) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(12) AL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f
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
      enddo !end loop over iislda

      close(12)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
ccc In the future, it should use the energy to reduce this 
ccc calculation. 
ccccccccccccccccccccccccccccccccccccccccccc
cccccc first, figure out the number of  degenerated energy group exist in ikpt+1 
      num_g=0
      Eg=1.D+20
      do ist=1,nst
       if(abs(Eband(ist,ikpt+1)-Eg).gt.1.E-4) then  ! new E degenerated group
       num_g=num_g+1
       Eg=Eband(ist,ikpt+1)
       ind_in=1                     ! index inside the group
       else           ! inside the old group
       ind_in=ind_in+1
       endif
       ind_g2s(ind_in,num_g)=ist      ! index, from the group to the original state number
       ind_s2g(ist)=num_g             ! index, from the state to group number
       ind_s2g2(ist)=ind_in
       nst_g(num_g)=ind_in           ! the number of states inside a group
      enddo
cccccc number of group equals num_g


ccccc assuming all the states are normalized in the same way!
      sum00=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
       sum00=sum00+cdabs(uc(i,j,k,1))**2
      enddo
      enddo
      enddo
ccccccccccccccccccccccccccccccccc
      allocate(cc_coeff(nst,nst))
      allocate(cc_tmp(nst,nst))
      allocate(rotate(nst,nst))
      allocate(rotate_tmp(nst,nst))

      cc_tmp=dcmplx(0.d0,0.d0)
      
      do iislda=ispin_i,ispin_f
      do ist1=iw_i,iw_f
       do iproc=1,nnodes
        read(11) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
        if(inode.eq.1) write(10) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
        do ii=1,nr_n
        jj=ii+(iproc-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        uc_tmp(i,j,k)=dcmplx(ucR(ii),ucI(ii))
        enddo
       enddo !end loop over iproc
       do 301 ig=1,num_g
        if(abs(Eband(ist1,ikpt)-
     &     Eband(ind_g2s(1,ig),ikpt+1)).gt.0.2) goto 301
        do ind_gs=1,nst_g(ig)
          ist2=ind_g2s(ind_gs,ig)
          cc=dcmplx(0.d0,0.d0)
          do k=1,n3
          do j=1,n2
          do i=1,n1
          cc=cc+uc_tmp(i,j,k)*dconjg(uc(i,j,k,ist2))*cphase(i,j,k)
          enddo
          enddo
          enddo
          cc=cc/sum00
          cc_tmp(ist1,ist2)=cc
        enddo
301    continue
      enddo !end loop over ist1
      enddo !end loop over iislda
      
      if(inode.eq.1) then
        close(10)
        cc_coeff=cc_tmp
      else
        call mpi_recv(rotate,nst*nst,MPI_DOUBLE_COMPLEX,inode-2,100,
     &       MPI_COMM_WORLD,status,ierr)
        do ist2=1,nst
        do ist1=1,nst
          cc=dcmplx(0.d0,0.d0)
          do m=1,nst
          cc=cc+rotate(ist1,m)*cc_tmp(m,ist2)
          enddo
          cc_coeff(ist1,ist2)=cc
        enddo
        enddo
      endif

      rotate=dcmplx(0.d0,0.d0)
      do i=1,nst
      rotate(i,i)=dcmplx(1.d0,0.d0)
      enddo

      do 400 ist1=1,nst
cccccc note, for ikpt, it is one state (ist1) at a time. 
cccccc for ikpt+1, it is one group (ig) at a time.
       mg=0
       do 302 ig=1,num_g
        sum_g=0.d0
        if(abs(Eband(ist1,ikpt)-
     &    Eband(ind_g2s(1,ig),ikpt+1)).gt.0.2.or.
     &    nst_g(ig).eq.0) goto 302  ! only try connection within 0.2 eV,nst_g(ig).eq.0, the group has been used up

        do ind_gs=1,nst_g(ig)
        ist2=ind_g2s(ind_gs,ig)
        sum_g=sum_g+cdabs(cc_coeff(ist1,ist2))**2
        enddo    ! do ind_gs

       if(sum_g.gt.0.25) then    ! this is the group. the total projection must be larger than 0.25
       mg=ig
       goto 303
       endif

302    continue! do group
303    continue
      
       mgl(ist1)=mg

       if(mg.eq.0) then
       sum_g=0.d0
       connectX(ist1,ikpt)=0.d0    ! dead end
       iconnect(ist1,ikpt)=0       ! dead end
       goto 400
       endif

cccc construct the new state
       sum=0.d0
       E_tmp=0.d0
       cc_max=0.d0
       ind_gs_max=0
ccccccccccccccc
       do ind_gs=1,nst_g(mg)
       ist2=ind_g2s(ind_gs,mg)
       tmp_real=cdabs(cc_coeff(ist1,ist2))**2
       E_tmp=E_tmp+Eband(ist2,ikpt+1)*tmp_real
       sum=sum+tmp_real
       if(tmp_real.gt.cc_max) then
       cc_max=tmp_real
       ind_gs_max=ind_gs
       endif
       enddo
ccccccccccccccc
       E_tmp=E_tmp/sum     ! E_tmp is the energy, uc_tmp is the state         !
       cc_coeff(ist1,:)=cc_coeff(ist1,:)/dsqrt(sum)    ! renormalize the coeff, to be used below.
       connectX(ist1,ikpt)=sum

ccccccc  uc_tmp=cc_coeff(ist1,ist2)*uc2(ist2)

ccccccc  Now, we will store uc_tmp in ind_g2s(1,mg), E_tmp in Eband(ind_g2s(1,mg),ikpt+1)
ccccccc  We will reduce the number of states inside mg by 1. 
ccccccc  We will keep the remaining states inside mg (for the nst_g(mg)-1 states)


cccccccccccccccccccccccccccccccccccc
       if(nst_g(mg).eq.1) then   
       iconnect(ist1,ikpt)=ind_g2s(1,mg)         ! the connection information
       nst_g(mg)=nst_g(mg)-1          ! uc2 need not to be rewritten
       goto 400
       endif
ccccccc nst_g(mg).gt.1,  deal with the remaining states inside the group

cccc first, project out the uc_tmp from the remaining states
       do ind_gs=1,nst_g(mg)
        ist2=ind_g2s(ind_gs,mg)
        sum=0.d0
        cc_tmp(:,ist2)=cc_coeff(:,ist2)
        do i=1,nst_g(mg)
         m=ind_g2s(i,mg)
         rotate_tmp(ist2,m)=rotate(ist2,m)
         do j=1,nst_g(mg)
         n=ind_g2s(j,mg)
         rotate_tmp(ist2,m)=rotate_tmp(ist2,m)
     &      -dconjg(cc_coeff(ist1,ist2))*cc_coeff(ist1,n)*rotate(n,m)
         enddo
         sum=sum+cdabs(rotate_tmp(ist2,m))**2
         cc_tmp(:,ist2)=cc_tmp(:,ist2)
     &     -cc_coeff(ist1,ist2)*dconjg(cc_coeff(ist1,m))*cc_coeff(:,m)
        enddo
        sum=dsqrt(1.d0/sum)
        do i=1,nst_g(mg)
          m=ind_g2s(i,mg)
          rotate_tmp(ist2,m)=rotate_tmp(ist2,m)*sum
        enddo
        cc_tmp(:,ist2)=cc_tmp(:,ist2)*sum
       enddo

cccccc now, replace the ind_gs_max by the first one (i.e, get rid of the ind_gs_max one) 
       ist2_max=ind_g2s(ind_gs_max,mg)
       ist2_1=ind_g2s(1,mg)
       rotate_tmp(ist2_max,:)=rotate_tmp(ist2_1,:)
       cc_tmp(:,ist2_max)=cc_tmp(:,ist2_1)
       Eband(ist2_max,ikpt+1)=Eband(ist2_1,ikpt+1)

ccccc replace the first state in the group  by uc_tmp
       iconnect(ist1,ikpt)=ist2_1         ! the connection information
       do i=1,nst_g(mg)
        m=ind_g2s(i,mg)
        rotate_tmp(ist2_1,m)=dcmplx(0.d0,0.d0)
        do j=1,nst_g(mg)
          n=ind_g2s(j,mg)
          rotate_tmp(ist2_1,m)=rotate_tmp(ist2_1,m)
     &       +cc_coeff(ist1,n)*rotate(n,m)
        enddo
       enddo

       Eband(ist2_1,ikpt+1)=E_tmp

       do i=1,nst_g(mg)
        m=ind_g2s(i,mg)
        cc_coeff(:,m)=cc_tmp(:,m)
        do j=1,nst_g(mg)
        n=ind_g2s(j,mg)
        rotate(n,m)=rotate_tmp(n,m)
        enddo
       enddo
cccccccccccccccccccccccccccccccccccccc


cccccc now, reduce the nst_g(mg) by one, and re-order all the index
       nst_g(mg)=nst_g(mg)-1

       do ind_in=1,nst_g(mg) 
       ind_g2s(ind_in,mg)=ind_g2s(ind_in+1,mg)
       ist=ind_g2s(ind_in+1,mg)
       ind_s2g(ist)=mg
       ind_s2g2(ist)=ind_in
       enddo
ccccccccccccccc, now, do a Gram Schmidt orthogonalization among the nst_g(mg) states

       do i=1,nst_g(mg)
        m=ind_g2s(i,mg)
        do j=1,i-1
         n=ind_g2s(j,mg)
         cc=dcmplx(0.d0,0.d0)
         do l=1,nst
         cc=cc+rotate(m,l)*dconjg(rotate(n,l))
         enddo
         rotate(m,:)=rotate(m,:)-cc*rotate(n,:)
         cc_coeff(:,m)=cc_coeff(:,m)-dconjg(cc)*cc_coeff(:,n)
        enddo
        sum=0.d0
        do k=1,nst
        sum=sum+cdabs(rotate(m,k))**2
        enddo
        sum=dsqrt(1.d0/sum)
        rotate(m,:)=rotate(m,:)*sum
        cc_coeff(:,m)=cc_coeff(:,m)*sum
       enddo
 
cccccccccccccccccccccccccccccccccccccccccccccccc
400   continue
cccccccccccccccccccccccccccccccccccccccccccccc
cc       need to write out the wavefunction uc2, it has been changed, rotated
      if(inode.ne.nkpt0) then
      call mpi_send(rotate,nst*nst,MPI_DOUBLE_COMPLEX,inode,100,
     &     MPI_COMM_WORLD,ierr)
      endif

      fwr_all="wr.new."//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)
      open(10,file=fwr_all,form="unformatted")
      rewind(10)
      write(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      write(10) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes

      do iislda=ispin_i,ispin_f
       do ist=iw_i,iw_f
        uc_tmp=dcmplx(0.d0,0.d0)
        do m=1,nst
        uc_tmp(:,:,:)=uc_tmp(:,:,:)+rotate(ist,m)*uc(:,:,:,m)
        enddo
        do iproc=1,nnodes
         do ii=1,nr_n
         jj=ii+(iproc-1)*nr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         ucR(ii)=dreal(uc_tmp(i,j,k))
         ucI(ii)=dimag(uc_tmp(i,j,k))
         enddo
         write(10) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
        enddo !end loop over iproc
       enddo !end loop over m
      enddo !end loop over iislda

      close(10)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      deallocate(ucR)
      deallocate(ucI)
      deallocate(uc)
      deallocate(uc_tmp)
      deallocate(cphase)

      deallocate(cc_coeff)
      deallocate(cc_tmp)
      deallocate(rotate)
      deallocate(rotate_tmp)

ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      
      if(inode.ne.1) then
        call mpi_send(num_g,1,MPI_INTEGER,0,101,MPI_COMM_WORLD,ierr)
        call mpi_send(iconnect(1,ikpt),nst,MPI_INTEGER,0,102,
     &       MPI_COMM_WORLD,ierr)
        call mpi_send(mgl,nst,MPI_INTEGER,0,103,MPI_COMM_WORLD,ierr)
        call mpi_send(connectX(1,ikpt),nst,MPI_REAL8,0,104,
     &       MPI_COMM_WORLD,ierr)
      else
       open(10,file="connect.matrix.W2K2")
       rewind(10)
       write(10,*) nst,nkpt
       do ikpt=1,nkpt-1
        write(6,*) "ikpt=",ikpt
        if(ikpt.ne.1) then
        call mpi_recv(num_g,1,MPI_INTEGER,ikpt-1,101,MPI_COMM_WORLD,
     &       status,ierr)
        call mpi_recv(iconnect(1,ikpt),nst,MPI_INTEGER,ikpt-1,102,
     &       MPI_COMM_WORLD,status,ierr)
        call mpi_recv(mgl,nst,MPI_INTEGER,ikpt-1,103,MPI_COMM_WORLD,
     &       status,ierr)
        call mpi_recv(connectX(1,ikpt),nst,MPI_REAL8,ikpt-1,104,
     &       MPI_COMM_WORLD,status,ierr)
        endif
        write(6,*) "num_g=",num_g
        do ist11=1,nst
        write(10,700) ikpt,ist11,iconnect(ist11,ikpt),
     &   connectX(ist11,ikpt)
        write(6,444) ikpt,ist11,iconnect(ist11,ikpt),mgl(ist11),
     &   connectX(ist11,ikpt)
       enddo
       enddo
       close(10)
      endif
700   format(i4,2x,i4,4x,i4,4x,f9.6)
444   format("ikpt,ist11,iconn,mg,connX",4(i6,1x),f15.10)
      
      call mpi_barrier(MPI_COMM_WORLD,ierr)
ccccccccccccccccccccccccccccccccccccc
      if(inode.ne.1) then
       call mpi_send(Eband(1,ikpt+1),nst,MPI_REAL8,0,105,
     &       MPI_COMM_WORLD,ierr)
      else
       E_band=0.d0
ccccccccccc initial (begining) of the bands
ccccccccccc Note, this procedure is used, so it can start a band in the middle of the kpts

ccccc
       iflag_band=0
       iband=0
       do 800 ikpt11=1,nkpt-1
       if(ikpt11.ne.1) then
       call mpi_recv(Eband(1,ikpt11+1),nst,MPI_REAL8,ikpt11-1,105,
     &       MPI_COMM_WORLD,status,ierr)
       endif
       do 800 ist11=1,nst
        if(iflag_band(ist11,ikpt11).eq.1) goto 800
        iband=iband+1        ! the beginning of a new band

        if(iband.gt.3*nst) then
        write(6,*) "iband.gt.3*nst,stop"
        stop
        endif
        E_band(iband,1)=Eband(ist11,ikpt11)
        ist_band(iband,1)=ist11
        ikpt_band(iband,1)=ikpt11
        iflag_band(ist11,ikpt11)=1    ! used
        ist=ist11
        do ikpt=ikpt11,nkpt-1     ! do this just for a single band  
        if(connectX(ist,ikpt).lt.0.01) then
        nkpt_band(iband)=ikpt-ikpt11+1       ! end of the band
        goto 304
        endif
        
        ist=iconnect(ist,ikpt)
        E_band(iband,ikpt-ikpt11+2)=Eband(ist,ikpt+1)
        ist_band(iband,ikpt-ikpt11+2)=ist
        ikpt_band(iband,ikpt-ikpt11+2)=ikpt+1
        iflag_band(ist,ikpt+1)=1
        enddo
        
        nkpt_band(iband)=nkpt-ikpt11+1
304    continue
       if(nkpt_band(iband).lt.3) iband=iband-1      ! remove the short bands, at the top 
ccccccccccccccccccccccccccccccccccccc
800    continue
       nband=iband
ccccccccccccccccccccccccccccccccccccc
       open(11,file="band.W2K2")
       rewind(11)
ccccc Note, the k-point might not be correct here for plotting
       do i=1,nkpt
       write(11,110) (i-1.d0)/(nkpt-1),(E_band(j,i),j=1,nband)
       enddo
110    format(f5.3,2x,220(f11.5,1x))
       close(11)

       open(11,file="band.W2K2_plot")
       rewind(11)
       do i=1,nkpt
       write(11,110) (i-1.d0)/(nkpt-1),(Eband(j,i),j=1,nst)
       enddo
       close(11)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

       open(14,file="E_line_W2K2")
       write(14,*) nband,nkpt
       do j=1,nband
       write(14,*) j,nkpt_band(j)
       write(14,120) (E_band(j,i),i=1,nkpt_band(j))
       write(14,121) (ist_band(j,i),i=1,nkpt_band(j))
       write(14,121) (ikpt_band(j,i),i=1,nkpt_band(j))
       enddo
       close(14)
120    format(5(f10.6,1x))
121    format(5(i5,6x))
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      deallocate(Eband)
      deallocate(err_band)
      deallocate(E_band)
      deallocate(ist_band)
      deallocate(ikpt_band)
      deallocate(iflag_band)

      deallocate(iconnect)
      deallocate(connectX)
ccccc we will use the file as the data pass here, find_evan is 
ccccc  an independent code

      if(inode.eq.1) call find_evan()

      call mpi_finalize(ierr)

      stop
      end
