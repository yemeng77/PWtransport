      subroutine wave_decomp(ccy2_st,Ew,num_st,
     & num_run,idble,iposit,nnposit,ucw,
     & n1w,n1,n2,n3,cphase_ucw,dE_dk,ak_w,nline_w,
     & numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,
     & num_mx,imx,cphase,phase,num_stWell,uc,E_evan,E_evanC,
     & ist_evan,ikpt_evan,iGX_evan,num_evan,weight,a11,
     & num_iter_evan,dE_evan,cc_R)

ccccccccc  this subroutine generates the electrod state
ccccccccc  coefficients "ccy2_st" for "num_stWell" given input system wavefunction "uc". 
cccccccc  ccy2_st(i_electrode_st,j_Well_st) is the complex coeff. of the 
cccccccc  ith electrode state on the jth Well state. 
cccccccc  Note, i_electrode_st  ucw and ucw^* are in i*2-1, and i*2
cccccccc  the j_Well_st  well states uc and uc^* are in j, and j+num_stWell
cccccccc

      implicit double precision (a-h,o-z)

      parameter (nm=400)

      complex*16, allocatable, dimension (:,:,:) :: uc_test,
     &   uc_test2,S
      complex*16, allocatable, dimension (:,:,:,:) :: uc_R
      complex*16, allocatable, dimension (:,:) :: S_m,hh,vec1,vec2
      complex*16, allocatable, dimension (:) :: vec_tmp
      real*8, allocatable, dimension(:):: EE
      real*8, allocatable, dimension (:,:) :: E_st

      complex*16 ccy2_st(90,90),cc_R(90,90)
      complex*16 ucw(n1w,n2,n3,40)
      complex*16 cphase_ucw(400)
      complex*16 cphase(n1w,n2,n3)
      complex*16 uc(n1,n2,n3,50)
      real*8  phase(n1w,n2,n3)
      real*8 E_evan(2000),E_evanC(2000)
      integer ist_evan(2000),ikpt_evan(2000),iGX_evan(2000),i
     &     used_evan(2000)
       real*8 dE_dk(400),ak_w(400)
      real*8 weight(90),aI_tmp(400)
      integer nline_w(400)

      real*8 E_linew(nm,nm)

      integer ikpt_st1w(400),ikpt_st2w(400),
     &   i_st1w(400),i_st2w(400),
     &  numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm),
     &  num_mx(400),imx(10,400)

      real*8 x_st1w(400),x_st2w(400)
      real*8 kk1,kk2,kk3,k1k2,k2k3,
     &     k2k1,k3k1,k3k2,k1k3


      complex*16 ccm(200,200),ccA(200,200)
      complex*16 ccy(200,40),ccy2(200),ccy3(200)
      integer  ipiv(200),iposit(200),idble(200),ievan(200)

      complex*16 cc,cc1,cc2,cc3,cc_st

      character*7 fileh
 
      dE_min=0.d0
      write(6,*) "Ew=",Ew


      pi=4*datan(1.d0)

      iused_evan=0
      num_st=0
      iposit_next=1
      num_iter=0


      open(77,file="overlap.matrix",form="unformatted")
      rewind(77)
      read(77) ispin_i,ispin_f
      read(77) nst,nkpt
      read(77) iislda
      allocate(S(nst,nst,nkpt-1))
      do ikpt=1,nkpt-1
        read(77) ikptt
        do ist2=1,nst
        read(77) (S(ist1,ist2,ikpt),ist1=1,nst)
        enddo
      enddo
      close(77)
      allocate(S_m(nst,nst))


      allocate(E_st(nst,nkpt))
      E_st=0.d0
      do j=1,nstw    ! go through different band-line
      do i=1,numw(j)-1    ! the number of k-points in each band-line
       E_st(ist_linew(i,j),ikpt_linew(i,j))=E_linew(i,j)
      enddo
      enddo

      allocate(hh(nst,nst))
      allocate(EE(nst))
      allocate(vec_tmp(nst))
      allocate(vec1(nst,400))
      allocate(vec2(nst,400))


3000  continue
      num_iter=num_iter+1
ccccccccccccccccccccccccccccccccccccccc
      if(num_iter.gt.1) goto 1000
ccccccccc  find the running waves
      num_st_old=0

      do j=1,nstw    ! go through different band-line

      do i=1,numw(j)-1    ! the number of k-points in each band-line

      if((E_linew(i,j).le.Ew.and.E_linew(i+1,j).gt.Ew).or.
     &  (E_linew(i+1,j).lt.Ew.and.E_linew(i,j).ge.Ew)) then


      num_st=num_st+1
      iposit(num_st)=iposit_next
      idble(num_st)=2
      ievan(num_st)=0

      kpt=ikpt_linew(i,j)

      S_m(:,:)=S(:,:,kpt)
      do m1=1,nst
      do m2=1,m1-1
        cc=dcmplx(0.d0,0.d0)
        do m3=1,nst
           cc=cc+S_m(m3,m1)*dconjg(S_m(m3,m2))
        enddo
        do m3=1,nst
           S_m(m3,m1)=S_m(m3,m1)-cc*S_m(m3,m2)
        enddo
      enddo
        sum=0.d0
        do m3=1,nst
          sum=sum+cdabs(S_m(m3,m1))**2
        enddo
        sum=1/dsqrt(sum)
        do m3=1,nst
          S_m(m3,m1)=S_m(m3,m1)*sum
        enddo
      enddo  ! m1

      Ek1=E_linew(i,j)
      Ek2=E_linew(i+1,j)
      ak=(Ew-Ek1)/(Ek2-Ek1)
      ist_min=min(ist_linew(i,j),ist_linew(i+1,j))
      ist_max=max(ist_linew(i,j),ist_linew(i+1,j))

      tol_err=1.D-6

1007  continue
      call energy_interp(ak,nst,E_st(1,kpt),S_m,0,
     &  ist_min,ist_max,Ew,Ek,vec1(1,num_st))

      ak_mov=1.D-8
      call energy_interp(ak-ak_mov,nst,E_st(1,kpt),S_m,0,
     &  ist_min,ist_max,Ek,Ek1,vec_tmp)
      call energy_interp(ak+ak_mov,nst,E_st(1,kpt),S_m,0,
     &  ist_min,ist_max,Ek,Ek2,vec_tmp)
      dEdk=(Ek2-Ek1)/(2*ak_mov)

      if(dabs(Ek-Ew).gt.tol_err*1.D-1) then
        ak=ak+(Ew-Ek)/dEdk
        goto 1007
      endif

      do m1=1,nst
        do m2=1,nst
          S_m(m1,m2)=dconjg(S(m2,m1,kpt))
        enddo
      enddo
      do m1=1,nst
      do m2=1,m1-1
        cc=dcmplx(0.d0,0.d0)
        do m3=1,nst
           cc=cc+S_m(m3,m1)*dconjg(S_m(m3,m2))
        enddo
        do m3=1,nst
           S_m(m3,m1)=S_m(m3,m1)-cc*S_m(m3,m2)
        enddo
      enddo
        sum=0.d0
        do m3=1,nst
          sum=sum+cdabs(S_m(m3,m1))**2
        enddo
        sum=1/dsqrt(sum)
        do m3=1,nst
          S_m(m3,m1)=S_m(m3,m1)*sum
        enddo
      enddo  ! m1

      call energy_interp(ak,nst,E_st(1,kpt),S_m,1,
     &  ist_min,ist_max,Ek,Ek2,vec2(1,num_st))

      if(dabs(Ek2-Ew).gt.tol_err) then
        write(6,*) "Somthing is wrong with
     &   interpolation. akw,Ek1,Ek2,ak",dble(kpt)+ak,Ek,Ek2
        stop
      endif

      if(abs(ak).gt.1) then
      write(6,*) "ak.gt.1, strange, stop",ak
      stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
      ikpt_st1w(num_st)=ikpt_linew(i,j)
      ikpt_st2w(num_st)=ikpt_linew(i+1,j)  

      i_st1w(num_st)=ist_linew(i,j)
      i_st2w(num_st)=ist_linew(i+1,j)

      x_st1w(num_st) = 1.d0-ak
      x_st2w(num_st) = ak

      dE_dk(num_st)=dEdk/(pi/a11/(nkptw-1))
      ak_w(num_st)=ak+dble(kpt)

      nline_w(num_st)=j

      iposit_next=iposit_next+idble(num_st)

      cphase_ucw(num_st)=
     &  cdexp(dcmplx(0.d0,pi*(ak_w(num_st)-1.d0)/(nkptw-1)))

      endif
      enddo
      enddo

      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
      write(6,*) "The number of running waves =", num_st
      do i=1,num_st
        write(6,*) i_st1w(i),ak_w(i),dE_dk(i)
      enddo
      write(6,301) (ikpt_st1w(i),i=1,num_st)
      write(6,302) (i_st1w(i),i=1,num_st)
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
301   format("ikpt= ",20(i4,1x)) 
302   format("i_st= ",20(i4,1x)) 

      deallocate(S_m)
      deallocate(E_st)
      deallocate(hh)
      deallocate(EE)
      deallocate(vec_tmp)



      if(num_st.eq.0) then
      goto 1000      ! go straight to calculate the evanescence states
      else
      num_run=num_st
      goto 1001       ! jump over the following part for evanescence state calculation
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
1000  continue        ! calculate the evanescence states

      num_st_old=num_st

      dE_min=1000.d0
      do j=1,num_evan 

      dE=Ew-E_evan(j)
       if(iused_evan(j).eq.0.and.
     &  dE*E_evanC(j).le.0) then

       if(abs(dE).lt.dE_min) then
       dE_min=abs(dE)
       j_min=j
       endif
       endif
      enddo

       if(dE_min.gt.dE_evan) goto 5000      !  end
cccccccccccccccccccccccc
       iused_evan(j_min)=1
       num_st=num_st+1
      iposit(num_st)=iposit_next
       ievan(num_st)=1
      if(iGX_evan(j_min).eq.2) then
      idble(num_st)=2
      else
      idble(num_st)=1
      endif
      iposit_next=iposit_next+idble(num_st)
      
      ikpt_st1w(num_st)=ikpt_evan(j_min)
      if(ikpt_st1w(num_st).lt.nkptw) then
      ikpt_st2w(num_st)=ikpt_st1w(num_st)+1
      else
      ikpt_st2w(num_st)=ikpt_st1w(num_st)-1
      endif
      i_st1w(num_st)=ist_evan(j_min)
      i_st2w(num_st)=ist_evan(j_min)
      x_st1w(num_st)=1.d0
      x_st2w(num_st)=0.d0
      ak_w(num_st)=ikpt_evan(j_min)
      nline_w(num_st)=j_min
      dE_dk(num_st)=0.d0
      ak=ikpt_evan(j_min)
      cphase_ucw(num_st)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))

      write(6,*) "YYYYYYYYYYYYYYYYYY"
      write(*,*) "evanescent state,ind,ist,iGX,E_evan0,dE_min "
      write(6,501)  num_st-num_run,ist_evan(j_min),iGX_evan(j_min),
     &  E_evan(j_min),dE_min
      write(6,*) "YYYYYYYYYYYYYYYYYY"

501    format("evan_state,ind,ist,iGX,E_evan,dE_min ",
     & 	3(i3,1x),2(f9.6,1x))

cccccccccccccccccccccccccccccccccccccccccccccc
1001   continue       ! end calc. of running wave and evanescence states

ccccccccccccccccccccccccccccccc
      num_tot=iposit_next-1

      write(6,503) num_tot,num_run,num_tot-num_run,num_iter
503   format("num_tot,num_run,num_evan,num_iter",
     &  4(i3,1x))      

      if(num_iter.eq.1) then
      ccm=dcmplx(0.d0,0.d0)
      ccA=dcmplx(0.d0,0.d0)
      ccy=dcmplx(0.d0,0.d0)
      ccy2=dcmplx(0.d0,0.d0)
      ccy3=dcmplx(0.d0,0.d0)
      ipiv=0
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccc

      fileh="wr.new."
      do 1020 ii=num_st_old+1,num_st
      if(num_iter.eq.1) then
      call wave_electrode_interp(ikpt_st1w(ii),ikpt_st2w(ii),
     & nst,vec1(1,ii),vec2(1,ii),
     & x_st1w(ii),x_st2w(ii),
     & ucw(1,1,1,ii),n1w,n2,n3,fileh)
      else
      call wave_electrode(ikpt_st1w(ii),ikpt_st2w(ii),
     & i_st1w(ii),i_st2w(ii),x_st1w(ii),x_st2w(ii),
     & ucw(1,1,1,ii),n1w,n2,n3,fileh,1,cphase,phase)
      endif


        if(idble(ii).eq.1) then     ! make ucw real, Note, even for X-point evanescent states, we can make it real
        cc=dcmplx(0.d0,0.d0)
        sum=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        cc=cc+ucw(i,j,k,ii)**2
        sum=sum+abs(ucw(i,j,k,ii))**2
        enddo
        enddo
        enddo
        cc=1.d0/sqrt(cc/sum)
        sum1=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii)=dcmplx(dreal(cc*ucw(i,j,k,ii)),0.d0)
        sum1=sum1+abs(ucw(i,j,k,ii))**2
        enddo
        enddo
        enddo
        fact=dsqrt(sum/sum1)
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii)=fact*ucw(i,j,k,ii)
        enddo
        enddo
        enddo
       
        endif      ! idble(ii).eq.1


1020  continue

      deallocate(vec1)
      deallocate(vec2)
ccccccccccccccccccccccccccccccccccc
 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TEST ORTHOGONALITY
cccccccccc    test
        do ii1=num_st_old+1,num_st
        do ii2=1,ii1-1
        sum2=0.d0
        cc=dcmplx(0.d0,0.d0)
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        sum2=sum2+abs(ucw(i,j,k,ii2))**2
        cc=cc+ucw(i,j,k,ii1)*dconjg(ucw(i,j,k,ii2))
        enddo
        enddo
        enddo
        if(ii2.eq.1) sum000=sum2
        cc=cc/sum2
        orth1=abs(cc)

        if(abs(ak_w(ii1)-ak_w(ii2)).lt.1.E-3) then      ! ak_w, from 1 to 101, for the same ak, make sure they are orth
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii1)=ucw(i,j,k,ii1)-cc*ucw(i,j,k,ii2)
        enddo
        enddo
        enddo
        orth1=0.d0
        endif

        enddo    ! ii2

        sum=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        sum=sum+abs(ucw(i,j,k,ii1))**2
        enddo
        enddo
        enddo
        if(ii1.eq.1) sum000=sum
        fact=dsqrt(sum000/sum)
        sum=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii1)=ucw(i,j,k,ii1)*fact
        sum=sum+abs(ucw(i,j,k,ii1))**2
        enddo
        enddo
        enddo
        
        enddo     ! ii1

ccccccccccccccccccccccccccccccccccc
      do ii1=num_st_old+1,num_st
      do ii2=1,ii1
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc1=cc1+ucw(i,j,k,ii1)*dconjg(ucw(i,j,k,ii2))
      cc2=cc2+dconjg(ucw(i,j,k,ii1))*dconjg(ucw(i,j,k,ii2))
      enddo
      enddo
      enddo

      ccm(iposit(ii2),iposit(ii1))=cc1             ! increase the overlap matrix
      ccm(iposit(ii1),iposit(ii2))=dconjg(cc1)
      if(idble(ii1).eq.2) then
      ccm(iposit(ii2),iposit(ii1)+1)=cc2
      ccm(iposit(ii1)+1,iposit(ii2))=dconjg(cc2)
      endif
      if(idble(ii2).eq.2) then
      ccm(iposit(ii1),iposit(ii2)+1)=cc2
      ccm(iposit(ii2)+1,iposit(ii1))=dconjg(cc2)
      endif
      if(idble(ii2).eq.2.and.idble(ii1).eq.2) then
      ccm(iposit(ii2)+1,iposit(ii1)+1)=dconjg(cc1)
      ccm(iposit(ii1)+1,iposit(ii2)+1)=cc1
      endif 

      enddo
      enddo

       open(56,file='decomp_systemst.out')
ccccccccccccccccccccccccccccccccccccccccccccccc
      if(num_iter.eq.1) then
      allocate(uc_test(n1w,n2,n3))
      allocate(uc_test2(n1w,n2,n3))
      allocate(uc_R(n1w,n2,n3,num_stWell))
      endif

      num_bad=0
      do 600 num1=1,num_stWell
 
      do ii1=num_st_old+1,num_st
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc1=cc1+uc(i+nnposit,j,k,num1)*dconjg(ucw(i,j,k,ii1))
      cc2=cc2+uc(i+nnposit,j,k,num1)*ucw(i,j,k,ii1)
      enddo
      enddo
      enddo
      ccy(iposit(ii1),num1)=cc1
      if(idble(ii1).eq.2) then        ! not at the Gamma point, k.ne.0
      ccy(iposit(ii1)+1,num1)=cc2
      endif

      enddo


ccccccccccccccccccccccccccccccccccccccccc
ccc solve for: ccm(ii1,ii2)*ca(ii2)=ccy(ii1)
ccc then, uc= ucw(ii1)*ca(ii1)
cccccccccccccccccccccccccccccccccccccc
cccc lapack
      ccA=ccm
      do i=1,num_tot            ! num_tot, number of electrode states, counting both phi_i and phi_i^*
      ccy2(i)=ccy(i,num1)
      enddo

      call zgesv(num_tot,1,ccA,200,ipiv,ccy2,200,info)

      do i=1,num_tot
      cc=dcmplx(0.d0,0.d0)
      do j=1,num_tot
      cc=cc+ccm(i,j)*ccy2(j)
      enddo
      ccy3(i)=cc
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do i=1,num_tot
       ccy2_st(i,num1)=ccy2(i)
       enddo
       
       do i=1,num_st
       if(idble(i).eq.2) then
       ccy2_st(iposit(i),num1+num_stWell)=dconjg(ccy2(iposit(i)+1))         ! num1, from 1 to num_stWell, system states
       ccy2_st(iposit(i)+1,num1+num_stWell)=dconjg(ccy2(iposit(i)))
       endif
       if(idble(i).eq.1) then
       ccy2_st(iposit(i),num1+num_stWell)=dconjg(ccy2(iposit(i)))           ! this could be wrong for X point evanescent states, because it is not real
       endif

       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      uc_test=dcmplx(0.d0,0.d0)
      uc_test2=dcmplx(0.d0,0.d0)

      aIcurrent=0.d0
      do ii=1,num_st

      if(ievan(ii).eq.0) then
      fact2=1.d0
      else
      fact2=0.d0
      endif

      if(idble(ii).eq.2) then

      aIcurrent=aIcurrent+(cdabs(ccy2(iposit(ii)))**2-
     &   cdabs(ccy2(iposit(ii)+1))**2)*dE_dk(ii)


      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)+
     &  ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii))
      uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*(ccy2(iposit(ii))*
     &  ucw(i,j,k,ii)+ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii)))
      enddo
      enddo
      enddo
      else
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)
      uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*ccy2(iposit(ii))*
     &   ucw(i,j,k,ii)
      enddo
      enddo
      enddo
      endif

      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum1=0.d0
      sum2=0.d0
      diff=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      sum1=sum1+abs(uc_test(i,j,k))**2
      sum2=sum2+abs(uc(i+nnposit,j,k,num1))**2
      diff=diff+abs(uc(i+nnposit,j,k,num1)-uc_test(i,j,k))**2
ccccc THE CHANGE !
      uc_R(i,j,k,num1)=uc(i+nnposit,j,k,num1)-uc_test(i,j,k)      ! uc_test2 does not include the evanescent states
      enddo
      enddo
      enddo
      err=abs(1.d0-sum1/sum2)
      ibad=0
      if(err.gt.0.01.or.err*sum2.gt.0.001) then 
      num_bad=num_bad+1
      ibad=1
      endif
      if(num1.eq.1) then
      write(56,*) "***********************************************"
      write(56,*) "******** iteration", num_iter, "***************"
      endif
      write(56,306) num_iter,num1,sum1/sum2,
     &          aIcurrent,(abs(ccy2(i)),i=1,num_tot)

      
      weight(num1)=sum1/sum2
      aI_tmp(num1)=aIcurrent
      weight(num1+num_stWell)=sum1/sum2*0.99999999999d0   ! just to dist these two states
303   format("nnposit,ibad,sum,diff,abs_sum", 
     &  i4,2x,i2,f10.7,1x,2(E10.4,1x))
304   format("num_st,sum_fract,diff,sum_absolute,Current", 
     &  i4,2x,f10.7,1x,2(E10.4,1x),2x,E14.8)
305   format(i3,1x,f10.7,1x,E14.8,3x,20(E6.1,1x))
306   format(i3,i3,1x,f10.7,1x,E14.8,3x,20(E6.1,1x))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
600   continue


	if(num_bad.gt.0.and.num_iter.lt.num_iter_evan.and.
     &   dabs(dE_min).lt.dE_evan) goto 3000
5000  continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      sum=sum+abs(ucw(i,j,k,1))**2
      enddo
      enddo
      enddo

      fact=1.d0/sum

      do num1=1,num_stWell
      do num2=1,num_stWell
      cc=dcmplx(0.d0,0.d0)
      cc1=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc=cc+uc_R(i,j,k,num1)*dconjg(uc_R(i,j,k,num2))
      cc1=cc1+dconjg(uc_R(i,j,k,num1))*dconjg(uc_R(i,j,k,num2))
      enddo
      enddo
      enddo
      cc=cc*fact
      cc1=cc1*fact
      cc_R(num1,num2)=cc
      cc_R(num1+num_stWell,num2+num_stWell)=dconjg(cc)
      cc_R(num1+num_stWell,num2)=cc1
      cc_R(num1,num2+num_stWell)=dconjg(cc1)
      enddo
      enddo
      cc_R=dconjg(cc_R)     ! important
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      deallocate(uc_test)
ccccccc output is ccy2_st and ucw

      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
