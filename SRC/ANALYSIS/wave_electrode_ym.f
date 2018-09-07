      subroutine wave_electrode_ym(Ew,num_st,num_run,num_tot,ftail,
     & idble,iposit,ievan,ucw,n1,n2,n3,cphase_ucw,dE_dk,ak_w,nline_w)

cccccccc modify by Meng Ye, Sep 2018
cccccccc read electrode states from file

      implicit double precision (a-h,o-z)

      parameter (nm=200)

      complex*16 ucw(n1,n2,n3,nm)
      real*8 dE_dk(nm),ak_w(nm)
      complex*16 cphase_ucw(nm)
      real*8, allocatable, dimension (:) :: ucR,ucI
      integer iposit(nm),idble(nm),ievan(nm),nline_w(nm)
      complex*16 cc

      character*7 ftail

      open(77,file="wr.run_"//trim(ftail),form="unformatted")
      rewind(77)
      read(77) n1t,n2t,n3t,nnodes,nst
      if(n1.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
        write(6,*) "n1,n2,n3 changed,stop",n1,n2,n3,n1t,n2t,n3t
        stop
      endif
      read(77)
      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      num_st=0
      iposit_next=1
      do ist=1,nst
        num_st=num_st+1
        iposit(num_st)=iposit_next
        idble(num_st)=2
        ievan(num_st)=0
        read(77) ikpt,i_st,nline_w(num_st),ak_w(num_st),dE_dk(num_st),
     &    cphase_ucw(num_st)
        do iproc=1,nnodes
          read(77) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
          do ii=1,nr_n
            jj=ii+(iproc-1)*nr_n
            i=(jj-1)/(n2*n3)+1
            j=(jj-1-(i-1)*n2*n3)/n3+1
            k=jj-(i-1)*n2*n3-(j-1)*n3
            ucw(i,j,k,num_st)=dcmplx(ucR(ii),ucI(ii))
          enddo 
        enddo
        iposit_next=iposit_next+idble(num_st)
      enddo
      deallocate(ucR)
      deallocate(ucI)
      close(77)

      num_run=num_st

      open(77,file="wr.evan_"//trim(ftail),form="unformatted")
      rewind(77)
      read(77) n1t,n2t,n3t,nnodes,nst
      if(n1.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
        write(6,*) "n1,n2,n3 changed,stop",n1,n2,n3,n1t,n2t,n3t
        stop
      endif
      read(77)
      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(ucR(nr_n))
      allocate(ucI(nr_n))
      do ist=1,nst
        num_st=num_st+1
        iposit(num_st)=iposit_next
        ievan(num_st)=1
        read(77) ikpt,i_st,iGX_evan,nline_w(num_st),ak_w(num_st),
     &    dE_dk(num_st),cphase_ucw(num_st)
        do iproc=1,nnodes
          read(77) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
          do ii=1,nr_n
            jj=ii+(iproc-1)*nr_n
            i=(jj-1)/(n2*n3)+1
            j=(jj-1-(i-1)*n2*n3)/n3+1
            k=jj-(i-1)*n2*n3-(j-1)*n3
            ucw(i,j,k,num_st)=dcmplx(ucR(ii),ucI(ii))
          enddo 
        enddo
        if(iGX_evan.eq.2) then
          idble(num_st)=2
        else
          idble(num_st)=1
        endif
        iposit_next=iposit_next+idble(num_st)
      enddo
      deallocate(ucR)
      deallocate(ucI)
      close(77)

      num_tot=iposit_next-1

      do ii=1,num_st
        if(idble(ii).eq.1) then     ! make ucw real, Note, even for X-point evanescent states, we can make it real
          cc=dcmplx(0.d0,0.d0)
          sum=0.d0
          do k=1,n3
          do j=1,n2
          do i=1,n1
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
      enddo

ccccccccccccccccccccccccccccccccccc
 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TEST ORTHOGONALITY
cccccccccc    test
      do ii1=1,num_st
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

      return
      end