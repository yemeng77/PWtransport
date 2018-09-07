      subroutine wave_decomp_ym(ccy2_st,Ew,num_st,num_tot,
     & idble,iposit,ievan,ucw,n1,n2,n3,dE_dk,ak_w,num_stWell,
     & ucLR,weight,cc_R,num_wr_dis,nst)

ccccccccc  this subroutine generates the electrod state
ccccccccc  coefficients "ccy2_st" for "num_stWell" given input system wavefunction "uc". 
cccccccc  ccy2_st(i_electrode_st,j_Well_st) is the complex coeff. of the 
cccccccc  ith electrode state on the jth Well state. 
cccccccc  Note, i_electrode_st  ucw and ucw^* are in i*2-1, and i*2
cccccccc  the j_Well_st  well states uc and uc^* are in j, and j+num_stWell
cccccccc modify by Meng Ye, Sep 2018
cccccccc  now read electrode states from file

      implicit double precision (a-h,o-z)

      parameter (nm=200)
      parameter (nsm=400)

      complex*16, allocatable, dimension (:,:,:) :: uc_test,uc_test2
      complex*16, allocatable, dimension (:,:,:,:) :: uc_R
      complex*16 ccy2_st(nsm,nsm),cc_R(nsm,nsm)
      complex*16 ucw(n1,n2,n3,nm)
      complex*16 ucLR(n1,n2,n3,nst)

      real*8 weight(nsm),aI_tmp(nm)
      real*8 dE_dk(nm)

      complex*16 ccm(nsm,nsm),ccA(nsm,nsm)
      complex*16 ccy(nsm,nsm),ccy2(nsm),ccy3(nsm)
      integer ipiv(nm),iposit(nm),idble(nm),ievan(nm)
      integer num_wr_dis(2)

      complex*16 cc,cc1,cc2,cc3,cc_st

      ccm=dcmplx(0.d0,0.d0)
      ccA=dcmplx(0.d0,0.d0)
      ccy=dcmplx(0.d0,0.d0)
      ccy2=dcmplx(0.d0,0.d0)
      ccy3=dcmplx(0.d0,0.d0)
      ipiv=0
      write(6,*) n1,n2,n3,num_st,nst
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ii1=1,num_st
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

ccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(uc_test(n1,n2,n3))
      allocate(uc_test2(n1,n2,n3))
      allocate(uc_R(n1,n2,n3,num_stWell))

      ist=0
      do 600 num1=num_wr_dis(1),num_wr_dis(2)
      ist=ist+1
      do ii=1,num_st
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1
      cc1=cc1+ucLR(i,j,k,ist)*dconjg(ucw(i,j,k,ii))
      cc2=cc2+ucLR(i,j,k,ist)*ucw(i,j,k,ii)
      enddo
      enddo
      enddo
      ccy(iposit(ii),num1)=cc1
      if(idble(ii).eq.2) ccy(iposit(ii)+1,num1)=cc2
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

      call zgesv(num_tot,1,ccA,nm,ipiv,ccy2,nm,info)

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
     &      cdabs(ccy2(iposit(ii)+1))**2)*dE_dk(ii)
          do k=1,n3
          do j=1,n2
          do i=1,n1
          uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)+
     &      ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii))
          uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*(ccy2(iposit(ii))*
     &      ucw(i,j,k,ii)+ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii)))
          enddo
          enddo
          enddo
        else
          do k=1,n3
          do j=1,n2
          do i=1,n1
          uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)
          uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*ccy2(iposit(ii))*
     &      ucw(i,j,k,ii)
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
      do i=1,n1
      sum1=sum1+abs(uc_test(i,j,k))**2
      sum2=sum2+abs(ucLR(i,j,k,ist))**2
      diff=diff+abs(ucLR(i,j,k,ist)-uc_test(i,j,k))**2
ccccc THE CHANGE !
      uc_R(i,j,k,num1)=ucLR(i,j,k,ist)-uc_test2(i,j,k)      ! uc_test2 does not include the evanescent states
      enddo
      enddo
      enddo
      err=abs(1.d0-sum1/sum2)
      ibad=0
      if(err.gt.0.01.or.err*sum2.gt.0.001) then 
      num_bad=num_bad+1
      ibad=1
      endif
      
      weight(num1)=sum1/sum2
      aI_tmp(num1)=aIcurrent
      weight(num1+num_stWell)=sum1/sum2*0.99999999999d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
600   continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1
        sum=sum+abs(ucw(i,j,k,1))**2
      enddo
      enddo
      enddo

      fact=1.d0/sum
      do num1=num_wr_dis(1),num_wr_dis(2)
      do num2=num_wr_dis(1),num_wr_dis(2)
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
      deallocate(uc_test2)
      deallocate(uc_R)
ccccccc output is ccy2_st and ucw

      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
