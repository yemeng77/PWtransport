      subroutine find_wrl(Ew,n1w,n1,n2,n3,nnodes,
     &  dV,dE_evan,dk_same,AL,mstate)
cccccc AL is the lattice of the large cell
cccccc ALw is the lattice of the electrode unit cell

      implicit double precision (a-h,o-z)

      parameter (nm=400)

      real*8, allocatable, dimension (:,:,:) :: phase

      complex*16, allocatable, dimension (:,:,:) :: uc
      complex*16, allocatable, dimension (:,:,:) :: cphase
      complex*16, allocatable, dimension (:,:,:) :: ur
      complex*16, allocatable, dimension (:) :: ur_tmp


      complex*16 cphase_ucw(400)
      real*8 E_evan(500),E_evanC(500)
      integer ist_evan(500),ikpt_evan(500),iGX_evan(500),
     &     iband_evan(500)
      real*8 E_evan1(500),E_evan2(500)
      integer index(500)
      real*8 dE_min1(500),dE_min2(500),dE_min12(500)
       real*8 dE_dk(400),ak_w(400)

      real*8 E_linew(nm,nm)
      real*8 ALw(3,3),AL(3,3)

      integer ik1(400,10),ik2(400,10),ikw(400,10)
      real*8 Ekw(400,10)
      integer num1(400),num2(400),numw_new(nm)
      integer ikpt_st1w(400),ikpt_st2w(400),
     &   ikpt_st3w(400),
     &   i_st1w(400),
     &  i_st2w(400),i_st3w(400),
     &  numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)

      real*8 x_st1w(400),x_st2w(400),x_st3w(400)
      real*8 kk1,kk2,kk3,k1k2,k2k3,
     &     k2k1,k3k1,k3k2,k1k3


      complex*16 cc

      character*7 fileh


      ALw(:,1)=AL(:,1)*n1w*1.d0/n1
 
      dE_min=0.d0
      pi=4*datan(1.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	Ew1=Ew-(dV/2)
	Ew2=Ew+(dV/2)
	write(6,*) "Ew,Ew1,Ew2",Ew,Ew1,Ew2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Real band structure in GX
      open(14,file="E_line_W2K2")
      read(14,*) nstw,nkptw
      do j=1,nstw
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo

      read(14,*) 
      read(14,*) num_evanT
      do ii=1,num_evanT
      read(14,*) ii2,ikpt_evan(ii),ist_evan(ii),
     & iGX_evan(ii),iband_evan(ii),E_evan(ii),E_evanC(ii)
      enddo
      close(14)
ccccccccccccccccccccccccccccccccccccccccccccccccc

      allocate(cphase(n1w,n2,n3))
      allocate(phase(n1w,n2,n3))
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      phase(i,j,k)=pi*(i-1.d0)/(nkptw-1.d0)/n1w
      cphase(i,j,k)=cdexp(-dcmplx(0.d0,phase(i,j,k)))
      enddo
      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc End of preprocessing
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ik1=0
      ik2=0

ccccccccccccccccccccccccccccccccccccccc
ccccccccc  find the running waves

ccc find intersections with RBS for E-dV/2

      do j=1,nstw    ! go through different band-line

      ij=0
      do i=1,numw(j)-1    ! the number of k-points in each band-line

      if((E_linew(i,j).le.Ew1.and.E_linew(i+1,j).gt.Ew1).or.
     &  (E_linew(i+1,j).lt.Ew1.and.E_linew(i,j).ge.Ew1)) then
         ij=ij+1
         ik1(j,ij) = i
      endif
      enddo !end loop over i(kpts)

        num1(j)=ij
      enddo !end loop over j(bands)

ccc find intersections with RBS for E+dV/2

      do j=1,nstw    ! go through different band-line

      ij=0
      do i=1,numw(j)-1    ! the number of k-points in each band-line

      if((E_linew(i,j).le.Ew2.and.E_linew(i+1,j).gt.Ew2).or.
     &  (E_linew(i+1,j).lt.Ew2.and.E_linew(i,j).ge.Ew2)) then

         ij=ij+1
         ik2(j,ij) = i
      endif
      enddo 
         num2(j)=ij
      enddo !end loop over j(bands)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ikw=0
      Ekw=0.d0

      do j=1,nstw    ! go through different band-line
        n1k=num1(j)
        n2k=num2(j)
        nj=0

ccccccccc first, identify the same kpt wave from Ew1,Ew2
        if((n1k.ne.0).and.(n2k.ne.0)) then
        do i1=1,n1k
        do i2=1,n2k
         mk12=abs(ik1(j,i1)-ik2(j,i2))
	 dmk12=mk12/(nkptw-1.0)
         if(dmk12.lt.dk_same) then
         ik2(j,i2)=-100       ! will not be used later
         endif
	enddo 
	enddo 
        endif

      do i1=1,n1k
      nj=nj+1
      ikw(j,nj)=ik1(j,i1)
      Ekw(j,nj)=Ew1
      enddo
     
      do i2=1,n2k
      if(ik2(j,i2).gt.0) then
      nj=nj+1
      ikw(j,nj)=ik2(j,i2)
      Ekw(j,nj)=Ew2
      endif
      enddo
 
         numw_new(j)=nj
      enddo
      
201    format(3I5,e12.5)	
202    format(3I5,e12.5,1I5)	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      num_st=0

      do j=1,nstw    ! go through different band-line

      do ij=1,numw_new(j)    

      i = ikw(j,ij)

      Ewin=Ekw(j,ij)


      num_st=num_st+1


ccccccccc Now, find out the i1,i2,i3 for three point interpolation
ccccccccc Actually, for generating wrl, this might not be necessary, one point should
ccccccccc be good enough

      if(i.eq.1) then
       i33=i+2
      elseif(i.eq.(numw(j)-1)) then
       i33=i-1
      else
       if(abs(E_linew(i+2,j)-Ewin).gt.abs(E_linew(i-1,j)-Ewin)) then
          i33=i-1
       else
          i33=i+2
       endif
      endif

cccccccccccccccccccccccccccccccccccccccccc
cccc note, assuming (it is correct), ikpt_linew(i,j) and i has the same order. They can only be shifted by a number
      if(i33.eq.i+2) then
      i1=i
      i2=i+1
      i3=i+2
      else
      i1=i-1
      i2=i
      i3=i+1
      endif
ccccccccccccccccccccccccccccccccccccccccc
      E1=E_linew(i1,j)
      E2=E_linew(i2,j)
      E3=E_linew(i3,j)
ccccccccccccccccccccccccccccccccccccc
cccc assume E(k)=a1*k^2+a2*k+a3, and k=0 at i2, -1 at i1, 1 at i3 : E(-1)=E1,E(0)=E2,E(1)=E3
ccccc then
      a1=(E3+E1-2*E2)/2
      a2=(E3-E1)/2
      a3=E2
ccccc Then: E(k)=E1*(k^2-k)/2+E2*(1-k^2)+E3*(k^2+k)/2=E1*w1+E2*w2+E3*w3   ! the linear coeff for w1,w2,w3
ccccccccccccccccccccccccccc
cccc Then solve Ew=a1*ak^2+a2*ak+a3 to find ak
cccc ak is the kpoint distance from ikpt_linew(i2,j)
cccccccc
      ak1=(-a2+dsqrt(abs(a2**2-4*a1*(a3-Ewin))))/(2*a1)
      ak2=(-a2-dsqrt(abs(a2**2-4*a1*(a3-Ewin))))/(2*a1)
      if(abs(ak1).lt.abs(ak2)) then
      ak=ak1
      else
      ak=ak2
      endif

204    format(3I5,3(e12.5,1x))	


      if(abs(ak).gt.1) then
      write(6,*) "ak.gt.1, strange, stop",ak
      stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
      dE_dk3=(2*a1*ak+a2)/(pi/a11/(nkptw-1))
      dE_dk(num_st)=dE_dk3
      w1=(ak**2-ak)/2
      w2=1.d0-ak**2
      w3=(ak**2+ak)/2

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ikpt_st1w(num_st)=ikpt_linew(i1,j)
      ikpt_st2w(num_st)=ikpt_linew(i2,j) 
      ikpt_st3w(num_st)=ikpt_linew(i3,j) 

      i_st1w(num_st)=ist_linew(i1,j)
      i_st2w(num_st)=ist_linew(i2,j)
      i_st3w(num_st)=ist_linew(i3,j)

      x_st1w(num_st) = w1
      x_st2w(num_st) = w2
      x_st3w(num_st) = w3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ak_w(num_st)=ak+ikpt_linew(i2,j)

      cphase_ucw(num_st)=
     &  cdexp(dcmplx(0.d0,pi*(ak_w(num_st)-1.d0)/(nkptw-1)))

      enddo       ! different state within band
      enddo       !  different bands
ccccccccccccccccccccccccccccccccc


c      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
      write(6,*) "The number of running waves =", num_st
      write(6,301) (ikpt_st1w(i),i=1,num_st)
      write(6,302) (i_st1w(i),i=1,num_st)
c      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
301   format("ikpt= ",20(i4,1x)) 
302   format("i_st= ",20(i4,1x)) 

      num_run=num_st
ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   calculate the evanescence states

c      open(16,file="evanescent.out")
c      rewind(16) 
c      read(16,*) num_evanT
c      do ii=1,num_evanT
c      read(16,*) ii2,ikpt_evan(ii),ist_evan(ii),
c     & iGX_evan(ii),iband_evan(ii),E_evan(ii),E_evanC(ii)
c      enddo
c      close(16)



      ij=0
      do ii=1,num_evanT
      dE1=Ew1-E_evan(ii)
      dE2=Ew2-E_evan(ii)
      if((dabs(dE1).lt.dE_evan.and.
     &  dE1*E_evanC(ii).le.0.d0).or.
     & (dabs(dE2).lt.dE_evan.and.
     &  dE2*E_evanC(ii).le.0.d0)) then 
      ij=ij+1
      index(ij)=ii
      endif
      enddo

      num_evan=ij
       
 
ccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,num_evan
 
      ii=index(i)

      num_st = num_st+1

      ikpt_st1w(num_st)=ikpt_evan(ii)

      if(ikpt_st1w(num_st).lt.nkptw-1) then
      ikpt_st2w(num_st)=ikpt_st1w(num_st)+1
      ikpt_st3w(num_st)=ikpt_st1w(num_st)+2
      else
      ikpt_st2w(num_st)=ikpt_st1w(num_st)-1
      ikpt_st3w(num_st)=ikpt_st1w(num_st)-2
      endif
      i_st1w(num_st)=ist_evan(ii)
      i_st2w(num_st)=ist_evan(ii)
      i_st3w(num_st)=ist_evan(ii)
      x_st1w(num_st)=1.d0
      x_st2w(num_st)=0.d0
      x_st3w(num_st)=0.d0
      ak_w(num_st)=ikpt_evan(ii)
      dE_dk(num_st)=0.d0
      ak=ikpt_evan(ii)
      cphase_ucw(num_st)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))

      write(6,*) "YYYYYYYYYYYYYYYYYY"
      write(6,501)  i,ist_evan(ii),ikpt_evan(ii),
     &  iGX_evan(ii),E_evan(ii)

501    format("evanescent state:ind,ist,ikpt,E,iGX: ",
     & 	4(i3,1x),f9.6,1x)

      enddo !end loop over num_evan

cccccccccccccccccccccccccccccccccccccccccccccc
5000	continue
cccccccccccccccccccccccccccccccccccccccccccccc
      write(6,*) 'num_evan',num_evan
      num_tot=num_run+num_evan

      write(6,*) 'num_tot',num_tot      ! num of electrode states, counting both phi_i and phi_i^*

      mstate=num_tot

cccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccc

      fileh="wr.new."
      nr=n1*n2*n3
      nr_n=nr/nnodes

      allocate(uc(n1w,n2,n3))
      allocate(ur_tmp(nr_n))

      
      do 4000 istate=1,num_tot



      uc=dcmplx(0.d0,0.d0)

      call wave_electrode_3kpts_LW(ikpt_st1w(istate),
     & ikpt_st2w(istate),
     & ikpt_st3w(istate),
     & i_st1w(istate),i_st2w(istate),i_st3w(istate),
     & x_st1w(istate),x_st2w(istate),x_st3w(istate),
     & uc,n1w,n2,n3,fileh,1,cphase,phase)
ccccccccccccccccccccccccccccccccccc


ccc uc is the state that we will use to construct w_l (for each ist)
cccc This is a fixed parameter!!
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

       open(10,file="wr_test."//char(ii1+48)//char(ii2+48),
     &      form="unformatted")
       rewind(10)

       write(10) n1,n2,n3,nnodes
       write(10) AL

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
      if(iw.gt.n1w+1-nd)
     &       fac=(1.d0-dsin((iw-n1w-1)*pi/2.d0/nd))/2
      if(iw.lt.nd+1)
     &       fac=(1.d0+dsin((iw-1)*pi/2.d0/nd))/2
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

        write(10) ur_tmp
        enddo
       close(10)


4000   continue

      deallocate(uc)
      deallocate(ur_tmp)
      deallocate(phase)
      deallocate(cphase)

ccccccccccccccccccccccccccccccccccc
      return

      end subroutine find_wrl
