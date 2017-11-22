      program connect
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Lin-Wang Wang, Feb. 2011
cccccc   this program calculate the connectivity one after another k points for the band structure
cccccc   It also makes the necessary rotations for the degenerated states
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      real*8 AL(3,3)
      complex*16 cc
      complex*16, allocatable, dimension (:,:,:,:) :: uc1,uc2
      complex*16, allocatable, dimension (:,:,:) :: uc_tmp,cphase
      integer,allocatable, dimension(:,:) :: iconnect
      real*8,allocatable,dimension(:,:) :: connectX
      real*8,allocatable,dimension(:) :: ucR,ucR2 
      real*8,allocatable,dimension(:) :: ucI,ucI2 
      real*8, allocatable,dimension(:,:) :: Eband,E_band,err_band
      integer, allocatable,dimension(:,:) :: ist_band,ikpt_band,
     &   iflag_band
      integer ind_g2s(100,1000),ind_s2g(1000),ind_s2g2(1000),nst_g(100)
      integer nkpt_band(1000)
      complex*16 cc_coeff(100),cai

      character*20 fwr_all


      write(6,*) "input nst,nkpt (e.g., nkpt=101)"
      read(5,*) nst,nkpt
	nkpt0=nkpt-1
      allocate(Eband(nst,nkpt))
      allocate(err_band(nst,nkpt))
      allocate(E_band(3*nst,nkpt))     ! nband can be 3*nst
      allocate(ist_band(3*nst,nkpt))
      allocate(ikpt_band(3*nst,nkpt))
      allocate(iflag_band(nst,nkpt))

      write(6,*) "copy report to Ek.tmp"

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

      
      open(10,file="wr.out001",form="unformatted")
      rewind(10)
      open(11,file="wr.new.001",form="unformatted")
      rewind(11)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
      read(10) AL 	      
      write(11) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 
      write(11) AL 	      
      nr=n1*n2*n3
      nr_n=nr/nnodes

      write(6,*) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f 

      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      allocate(ucR2(nr_n))
      allocate(ucI2(nr_n))
      allocate(uc1(n1,n2,n3,nst))
      allocate(uc2(n1,n2,n3,nst))
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

      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f-1

      do iproc=1,nnodes

      read(10) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
      write(11) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)

       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc2(i,j,k,m)=dcmplx(ucR(ii),ucI(ii))
       enddo


      enddo

      enddo
      enddo

       close(10)
       close(11)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 3000 ikpt=1,nkpt-1

      
      write(6,*) "ikpt=",ikpt

       uc1=uc2

ccccccccccccccccccccccccccccccccccccccccccccccccccc
       kpt_dens=ikpt+1

      ikpt_100=kpt_dens/100
      ikpt_10=(kpt_dens-ikpt_100*100)/10
      ikpt_1=kpt_dens-ikpt_100*100-ikpt_10*10
      fwr_all="wr.out"//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)


      open(10,file=fwr_all,form="unformatted")
      rewind(10)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(10) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes

      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f-1

      do iproc=1,nnodes

      read(10) (ucR2(i),i=1,nr_n), (ucI2(i),i=1,nr_n)

       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc2(i,j,k,m)=dcmplx(ucR2(ii),ucI2(ii))
       enddo

      enddo !end loop over iproc

      enddo !end loop over m
      enddo !end loop over iislda


       close(10)
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

       write(6,*) "num_g=",num_g
cccccc number of group equals num_g


ccccc assuming all the states are normalized in the same way!
       sum00=0.d0
       do k=1,n3
       do j=1,n2
       do i=1,n1
       sum00=sum00+abs(uc1(i,j,k,1))**2
       enddo
       enddo
       enddo
ccccccccccccccccccccccccccccccccc


       do 400 ist11=1,nst

cccccc note, for ikpt, it is one state (ist11) at a time. 
cccccc for ikpt+1, it is one group (ig) at a time.
       
       mg=0
       do 301 ig=1,num_g
       
       sum_g=0.d0
       if(abs(Eband(ist11,ikpt)-
     &  Eband(ind_g2s(1,ig),ikpt+1)).gt.0.2.or.
     &  nst_g(ig).eq.0) goto 301  ! only try connection within 0.2 eV,nst_g(ig).eq.0, the group has been used up
   
       do ind_gs=1,nst_g(ig)
       ist2=ind_g2s(ind_gs,ig)
       cc=dcmplx(0.d0,0.d0)
       do k=1,n3
       do j=1,n2
       do i=1,n1
       cc=cc+uc1(i,j,k,ist11)*dconjg(uc2(i,j,k,ist2))*
     &  cphase(i,j,k)
       enddo
       enddo
       enddo
       cc_coeff(ind_gs)=cc/sum00     ! normalize

       sum_g=sum_g+cdabs(cc_coeff(ind_gs))**2

       enddo    ! do ind_gs

       if(sum_g.gt.0.25) then    ! this is the group. the total projection must be larger than 0.25
       mg=ig
       goto 302
       endif

301    continue   ! do group
302    continue

       if(mg.eq.0) sum_g=0.d0

       if(mg.eq.0) then
       connectX(ist11,ikpt)=0.d0    ! dead end
       iconnect(ist11,ikpt)=0       ! dead end
       goto 399
       endif
      

cccc construct the new state
       sum=0.d0
       E_tmp=0.d0
       cc_max=0.d0
       ind_gs_max=0
       uc_tmp=dcmplx(0.d0,0.d0)
ccccccccccccccc
       do ind_gs=1,nst_g(mg)

       ist2=ind_g2s(ind_gs,mg)
       E_tmp=E_tmp+Eband(ist2,ikpt+1)*cdabs(cc_coeff(ind_gs))**2
       sum=sum+cdabs(cc_coeff(ind_gs))**2

       if(cdabs(cc_coeff(ind_gs))**2.gt.cc_max) then
       cc_max=cdabs(cc_coeff(ind_gs))**2
       ind_gs_max=ind_gs
       endif

       uc_tmp(:,:,:)=uc_tmp(:,:,:)+cc_coeff(ind_gs)*uc2(:,:,:,ist2)

       enddo
ccccccccccccccc
       E_tmp=E_tmp/sum     ! E_tmp is the energy, uc_tmp is the state         !
       uc_tmp=uc_tmp/dsqrt(sum)     ! normalize it
       cc_coeff=cc_coeff/dsqrt(sum)    ! renormalize the coeff, to be used below. 
       connectX(ist11,ikpt)=sum

ccccccc  Now, we will store uc_tmp in ind_g2s(1,mg), E_tmp in Eband(ind_g2s(1,mg),ikpt+1)
ccccccc  We will reduce the number of states inside mg by 1. 
ccccccc  We will keep the remaining states inside mg (for the nst_g(mg)-1 states)


cccccccccccccccccccccccccccccccccccc
       if(nst_g(mg).eq.1) then   

       iconnect(ist11,ikpt)=ind_g2s(1,mg)         ! the connection information
       nst_g(mg)=nst_g(mg)-1          ! uc2 need not to be rewritten

       else    ! nst_g(mg).gt.1,  deal with the remaining states inside the group

cccc first, project out the uc_tmp from the remaining states
       do ind_gs=1,nst_g(mg)
       ist2=ind_g2s(ind_gs,mg)
       sum=0.d0
       cc=dcmplx(0.d0,0.d0)
       do k=1,n3
       do j=1,n2
       do i=1,n1
       uc2(i,j,k,ist2)=uc2(i,j,k,ist2)-dconjg(cc_coeff(ind_gs))*
     &  uc_tmp(i,j,k)

       sum=sum+cdabs(uc2(i,j,k,ist2))**2
       cc=cc+uc2(i,j,k,ist2)*dconjg(uc_tmp(i,j,k))
       enddo
       enddo
       enddo
       uc2(:,:,:,ist2)=uc2(:,:,:,ist2)*dsqrt(sum00/sum)
       enddo    ! ind_gs

cccccc now, replace the ind_gs_max by the first one (i.e, get rid of the ind_gs_max one) 
       ist2_max=ind_g2s(ind_gs_max,mg)
       ist2_1=ind_g2s(1,mg)
       uc2(:,:,:,ist2_max)=uc2(:,:,:,ist2_1)
       Eband(ist2_max,ikpt+1)=Eband(ist2_1,ikpt+1)

ccccc replace the first state in the group  by uc_tmp
       iconnect(ist11,ikpt)=ist2_1         ! the connection information
       uc2(:,:,:,ist2_1)=uc_tmp(:,:,:)
       Eband(ist2_1,ikpt+1)=E_tmp
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

       do ind_in1=1,nst_g(mg)
       ist1=ind_g2s(ind_in1,mg)
       
       do ind_in2=1,ind_in1-1
       ist2=ind_g2s(ind_in2,mg)

       cc=dcmplx(0.d0,0.d0)
       do k=1,n3
       do j=1,n2
       do i=1,n1
       cc=cc+uc2(i,j,k,ist1)*dconjg(uc2(i,j,k,ist2))
       enddo
       enddo
       enddo
       cc=cc/sum00
       uc2(:,:,:,ist1)=uc2(:,:,:,ist1)-cc*uc2(:,:,:,ist2)
       enddo       ! ind_in2

       sum=0.d0
       do k=1,n3
       do j=1,n2
       do i=1,n1
       sum=sum+cdabs(uc2(i,j,k,ist1))**2
       enddo
       enddo
       enddo
       uc2(:,:,:,ist1)=uc2(:,:,:,ist1)*dsqrt(sum00/sum)
       enddo   ! ind_in1

       endif  ! nst_g(mg).gt.1  
cccccccccccccccccccccccccccccccccccccccccccccccc

399    continue
       write(6,444) ikpt,ist11,iconnect(ist11,ikpt),mg,
     &   connectX(ist11,ikpt)
444    format("ikpt,ist11,iconn,mg,connX",4(i6,1x),f15.10)

400    continue     ! do ist11=1,nst
cccccccccccccccccccccccccccccccccccccccccccccc
cc       need to write out the wavefunction uc2, it has been changed, rotated
      fwr_all="wr.new."//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)

      open(10,file=fwr_all,form="unformatted")
      rewind(10)
      write(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      write(10) AL
      nr=n1*n2*n3
      nr_n=nr/nnodes

      do iislda=ispin_i,ispin_f
      do m=iw_i,iw_f-1

      do iproc=1,nnodes


       do ii=1,nr_n
       jj=ii+(iproc-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       ucR2(ii)=dreal(uc2(i,j,k,m))
       ucI2(ii)=dimag(uc2(i,j,k,m))
       enddo

      write(10) (ucR2(i),i=1,nr_n), (ucI2(i),i=1,nr_n)

      enddo !end loop over iproc

      enddo !end loop over m
      enddo !end loop over iislda

       close(10)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc

3000    continue     ! do ikpt=1,nst

       open(10,file="connect.matrix.W2K2")
       rewind(10)
       write(10,*) nst,nkpt
       do ikpt=1,nkpt-1 
       do ist=1,nst
       write(10,700) ikpt,ist,iconnect(ist,ikpt),
     &   connectX(ist,ikpt)
       enddo
       enddo
       close(10)
700    format(i4,2x,i4,4x,i4,4x,f9.6)


ccccccccccccccccccccccccccccccccccccc
       E_band=0.d0
ccccccccccc initial (begining) of the bands
ccccccccccc Note, this procedure is used, so it can start a band in the middle of the kpts

ccccc
       iflag_band=0
       iband=0
       do 800 ikpt11=1,nkpt-1
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
       goto 303
       endif

       ist=iconnect(ist,ikpt)
       E_band(iband,ikpt-ikpt11+2)=Eband(ist,ikpt+1)
       ist_band(iband,ikpt-ikpt11+2)=ist
       ikpt_band(iband,ikpt-ikpt11+2)=ikpt+1
       iflag_band(ist,ikpt+1)=1
       enddo
       nkpt_band(iband)=nkpt-ikpt11+1
303    continue
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
ccccc we will use the file as the data pass here, find_evan is 
ccccc  an independent code

       call find_evan()

       stop
       end
