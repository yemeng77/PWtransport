      program fix_band
cccccc Written by Meng Ye. Aug, 2018

      implicit double precision (a-h,o-z)

      parameter (nm=2000)

      real*8 E_linew(nm,nm)
      integer numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)
      integer numw_tmp,ist_tmp(nm),ikpt_tmp(nm)
      real*8 E_tmp(nm)

      integer, allocatable, dimension (:) :: iband_max,iband_min
      integer, allocatable, dimension (:) :: inum_max,inum_min
      real*8, allocatable, dimension (:) :: E_max,E_min

      complex*16, allocatable, dimension (:,:,:) :: S,U
      complex*16, allocatable, dimension (:,:,:) :: cc1,cc2
      complex*16, allocatable, dimension (:) :: cc1_tmp1,cc1_tmp2
      complex*16, allocatable, dimension (:) :: cc2_tmp1,cc2_tmp2
      complex*16, allocatable, dimension (:) :: tmp1,tmp2,tmp3,tmp4

      complex*16 phase
      integer dk

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Real band structure in GX
      open(14,file="E_line_W2K2")
      rewind(14)
      read(14,*) nstw,nkptw,nintep
      do j=1,nstw
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      read(14,*)
      read(14,*) nevan
      do ii=1,nevan
       read(14,*)
      enddo
      read(14,*)
      read(14,*) nevan
      allocate(iband_max(nevan))
      allocate(iband_min(nevan))
      allocate(inum_max(nevan))
      allocate(inum_min(nevan))
      allocate(E_max(nevan))
      allocate(E_min(nevan))
      do ii=1,nevan
        read(14,*) iit,inum_max(ii),iband_max(ii),
     &  inum_min(ii),iband_min(ii),E_max(ii),E_min(ii)
      enddo
      close(14)

      open(12,file="wr.overlap",form="unformatted")
      rewind(12)
      read(12) nkpt0,nst
      open(13,file="wr.rotation",form="unformatted")
      rewind(13)
      read(13) nkptt,nstt
      if(nkpt0.ne.nkptt.or.nst.ne.nstt) then
        write(6,*) "nkpt or nst in wr.overlap and wr.rotation are not 
     &        same,",nkpt,nst,nkptt,nstt
        stop
      endif
      if((nkpt0-1)*nintep+1.ne.nkptw) then
        write(6,*) "nkpt0, nkpt, nintep are not match,"
     &       ,nkpt0,nkptw,nintep
        stop
      endif

      open(11,file="wr.interp.new",form="unformatted")
      rewind(11)
      read(11)
      allocate(cc1(nst,nst,nkptw))
      allocate(cc2(nst,nst,nkptw))
      do ikpt=1,nkptw
        do ist=1,nst
          read(11) cc1(:,ist,ikpt)
          read(11) cc2(:,ist,ikpt)
        enddo
      enddo
      close(11)

      allocate(S(nst,nst,nkpt0-1))
      allocate(U(nst,nst,nkpt0-1))
      do ikpt=1,nkpt0-1
        read(12) S(:,:,ikpt)
        do i=1,4
          read(12)
        enddo
        read(13) U(:,:,ikpt)
      enddo
      close(12)
      close(13)

      allocate(tmp1(nst))
      allocate(tmp2(nst))
      allocate(tmp3(nst))
      allocate(tmp4(nst))
      allocate(cc1_tmp1(nst))
      allocate(cc2_tmp1(nst))
      allocate(cc1_tmp2(nst))
      allocate(cc2_tmp2(nst))

      E_tol=1.0d-2

      do ii=1,nevan
        if(E_min(ii)-E_max(ii).gt.E_tol) cycle

        i1=inum_max(ii)
        j1=iband_max(ii)
        E1=E_max(ii)
        i2=inum_min(ii)
        j2=iband_min(ii)
        E2=E_min(ii)
        kpt1=ikpt_linew(i1,j1)
        kpt2=ikpt_linew(i2,j2)    
        i11=i1-kpt1+kpt2
        E11=E_linew(i11,j1)
        dE=4.0d0*(E2-E11)
        do dk=1,nkptw
          dE1=E_linew(i2+dk,j2)-E_linew(i11-dk,j1)
          dE2=E_linew(i2-dk,j2)-E_linew(i11+dk,j1)
          if(dE1.ge.dE.and.dE2.ge.dE) goto 707
        enddo
707     continue
        
        if(i11+dk.gt.numw(j1)) dk=numw(j1)-i11
        if(i11-dk.lt.1) dk=i11-1
        if(i2+dk.gt.numw(j2)) dk=numw(j2)-i2
        if(i2-dk.lt.1) dk=i2-1

        ist1_1=ist_linew(i11-dk,j1)
        ist2_1=ist_linew(i2+dk,j2)
        ist1_2=ist_linew(i11+dk,j1)
        ist2_2=ist_linew(i2-dk,j2)
        ikpt1_1=ikpt_linew(i11-dk,j1)
        ikpt2_1=ikpt_linew(i2+dk,j2)
        ikpt1_2=ikpt_linew(i11+dk,j1)
        ikpt2_2=ikpt_linew(i2-dk,j2)

        cc1_tmp1=cc1(:,ist1_1,ikpt1_1)
        cc2_tmp1=cc2(:,ist1_1,ikpt1_1)
        cc1_tmp2=cc1(:,ist2_1,ikpt2_1)
        cc2_tmp2=cc2(:,ist2_1,ikpt2_1)

        ikpt1=(ikpt1_1-1)/nintep+1
        ikpt2=(ikpt2_1-1)/nintep+1
        do iikpt=ikpt1,ikpt2-1
          tmp1=dcmplx(0.d0,0.d0)
          tmp2=dcmplx(0.d0,0.d0)
          do m=1,nst
            do m1=1,nst
            tmp1(m)=tmp1(m)+U(m1,m,iikpt)*cc1_tmp1(m1)
            tmp2(m)=tmp2(m)+U(m1,m,iikpt+1)*cc2_tmp1(m1)
            enddo
          enddo
          cc1_tmp1=tmp1
          cc2_tmp1=tmp2
        enddo
        phase=dcmplx(0.d0,0.d0)
        do m1=1,nst
        do m2=1,nst
          phase=phase+cc1_tmp1(m1)*dconjg(cc2_tmp2(m2))*S(m1,m2,ikpt2)
     &      +cc2_tmp1(m1)*dconjg(cc1_tmp2(m2)*S(m2,m1,ikpt2))
        enddo
        phase=phase+cc1_tmp1(m1)*dconjg(cc1_tmp2(m1))
     &    +cc2_tmp1(m1)*dconjg(cc2_tmp2(m1))
        enddo
        phase=phase/cdabs(phase)

        do i=i11-dk+1,i11+dk-1
          w1=dble(i11+dk-i)/(2.0d0*dble(dk))
          w2=1.0d0-w1
          E_linew(i,j1)=w1*E_linew(i11-dk,j1)+w2*E_linew(i2+dk,j2)
          ikpt=ikpt_linew(i,j1)
          ist=ist_linew(i,j1)
          do m=1,nst
            tmp1(m)=w1*cc1_tmp1(m)+w2*cc1_tmp2(m)*phase
            tmp2(m)=w1*cc2_tmp1(m)+w2*cc2_tmp2(m)*phase
          enddo
          do iikpt=ikpt2-1,(ikpt-1)/nintep+1,-1
            tmp3=dcmplx(0.d0,0.d0)
            tmp4=dcmplx(0.d0,0.d0)
            do m=1,nst
              do m1=1,nst
              tmp3(m)=tmp3(m)+dconjg(U(m,m1,iikpt))*tmp1(m1)!U is an unitary matrix
              tmp4(m)=tmp4(m)+dconjg(U(m,m1,iikpt+1))*tmp2(m1)
              enddo
            enddo
            tmp1=tmp3
            tmp2=tmp4
          enddo
          iikpt=(ikpt-1)/nintep+1
          sum=0.d0
          do m1=1,nst
          do m2=1,nst
          sum=sum+2*dreal(tmp1(m1)*dconjg(tmp2(m2))*S(m1,m2,iikpt))
          enddo
          sum=sum+cdabs(tmp1(m1))**2+cdabs(tmp2(m1))**2
          enddo
          sum=dsqrt(1.d0/sum)
          do m=1,nst
          cc1(m,ist,ikpt)=tmp1(m)*sum
          cc2(m,ist,ikpt)=tmp2(m)*sum
          enddo
        enddo

        cc1_tmp1=cc1(:,ist1_2,ikpt1_2)
        cc2_tmp1=cc2(:,ist1_2,ikpt1_2)
        cc1_tmp2=cc1(:,ist2_2,ikpt2_2)
        cc2_tmp2=cc2(:,ist2_2,ikpt2_2)
        ikpt1=(ikpt1_2-1)/nintep+1
        ikpt2=(ikpt2_2-1)/nintep+1
        do iikpt=ikpt2,ikpt1-1
          tmp1=dcmplx(0.d0,0.d0)
          tmp2=dcmplx(0.d0,0.d0)
          do m=1,nst
            do m1=1,nst
            tmp1(m)=tmp1(m)+U(m1,m,iikpt)*cc1_tmp2(m1)
            tmp2(m)=tmp2(m)+U(m1,m,iikpt+1)*cc2_tmp2(m1)
            enddo
          enddo
          cc1_tmp2=tmp1
          cc2_tmp2=tmp2
        enddo
        phase=dcmplx(0.d0,0.d0)
        do m1=1,nst
        do m2=1,nst
          phase=phase+cc1_tmp2(m1)*dconjg(cc2_tmp1(m2))*S(m1,m2,ikpt1)
     &      +cc2_tmp2(m1)*dconjg(cc1_tmp1(m2)*S(m2,m1,ikpt1))
        enddo
        phase=phase+cc1_tmp2(m1)*dconjg(cc1_tmp1(m1))
     &    +cc2_tmp2(m1)*dconjg(cc2_tmp1(m1))
        enddo
        phase=phase/cdabs(phase)

        do i=i2-dk+1,i2+dk-1
          w1=dble(i2+dk-i)/(2.0d0*dble(dk))
          w2=1.0d0-w1
          E_linew(i,j2)=w1*E_linew(i2-dk,j2)+w2*E_linew(i11+dk,j1)
          ikpt=ikpt_linew(i,j2)
          ist=ist_linew(i,j2)
          do m=1,nst
            tmp1(m)=w1*cc1_tmp2(m)+w2*cc1_tmp1(m)*phase
            tmp2(m)=w1*cc2_tmp2(m)+w2*cc2_tmp1(m)*phase
          enddo
          do iikpt=ikpt1-1,(ikpt-1)/nintep+1,-1
            tmp3=dcmplx(0.d0,0.d0)
            tmp4=dcmplx(0.d0,0.d0)
            do m=1,nst
              do m1=1,nst
              tmp3(m)=tmp3(m)+dconjg(U(m,m1,iikpt))*tmp1(m1)!U is an unitary matrix
              tmp4(m)=tmp4(m)+dconjg(U(m,m1,iikpt+1))*tmp2(m1)
              enddo
            enddo
            tmp1=tmp3
            tmp2=tmp4
          enddo
          iikpt=(ikpt-1)/nintep+1
          sum=0.d0
          do m1=1,nst
          do m2=1,nst
          sum=sum+2*dreal(tmp1(m1)*dconjg(tmp2(m2))*S(m1,m2,iikpt))
          enddo
          sum=sum+cdabs(tmp1(m1))**2+cdabs(tmp2(m1))**2
          enddo
          sum=dsqrt(1.d0/sum)
          do m=1,nst
          cc1(m,ist,ikpt)=tmp1(m)*sum
          cc2(m,ist,ikpt)=tmp2(m)*sum
          enddo
        enddo
      
        E_tmp=E_linew(:,j1)
        ist_tmp=ist_linew(:,j1)
        ikpt_tmp=ikpt_linew(:,j1)
        numw_tmp=numw(j1)
        do i=dk,numw(j2)-i2
          E_linew(i11+i,j1)=E_linew(i2+i,j2)
          ikpt_linew(i11+i,j1)=ikpt_linew(i2+i,j2)
          ist_linew(i11+i,j1)=ist_linew(i2+i,j2)
        enddo
        numw(j1)=i11+numw(j2)-i2
        do i=dk,numw_tmp-i11
          E_linew(i2+i,j2)=E_tmp(i11+i)
          ikpt_linew(i2+i,j2)=ikpt_tmp(i11+i)
          ist_linew(i2+i,j2)=ist_tmp(i11+i)
        enddo
        numw(j2)=i2+numw_tmp-i11
        
        do iii=ii+1,nevan
          if(iband_max(iii).eq.j1.and.inum_max(iii).ge.i11+dk) then
            iband_max(iii)=j2
            inum_max(iii)=inum_max(iii)-i11+i2
          else if(iband_max(iii).eq.j2.and.inum_max(iii).ge.i2+dk) then
            iband_max(iii)=j1
            inum_max(iii)=inum_max(iii)-i2+i11
          endif
          if(iband_min(iii).eq.j1.and.inum_min(iii).ge.i11+dk) then
            iband_min(iii)=j2
            inum_min(iii)=inum_min(iii)-i11+i2
          else if(iband_min(iii).eq.j2.and.inum_min(iii).ge.i2+dk) then
            iband_min(iii)=j1
            inum_min(iii)=inum_min(iii)-i2+i11
          endif
        enddo

      enddo

      deallocate(tmp1)
      deallocate(tmp2)
      deallocate(tmp3)
      deallocate(tmp4)
      deallocate(cc1_tmp1)
      deallocate(cc2_tmp1)
      deallocate(cc1_tmp2)
      deallocate(cc2_tmp2)
      deallocate(S)
      deallocate(U)

      deallocate(iband_max)
      deallocate(iband_min)
      deallocate(inum_max)
      deallocate(inum_min)
      deallocate(E_max)
      deallocate(E_min)

      open(11,file="wr.interp.new.fix",form="unformatted")
      rewind(11)
      write(11) nkptw,nst,nkpt0,nintep
      do ikpt=1,nkptw
        do ist=1,nst
          write(11) cc1(:,ist,ikpt)
          write(11) cc2(:,ist,ikpt)
        enddo
      enddo
      close(11)

      deallocate(cc1)
      deallocate(cc2)

      open(14,file="E_line_W2K2.fix")
      rewind(14)
      write(14,*) nstw,nkptw,nintep
      do j=1,nstw
      write(14,*) j,numw(j)
      write(14,120) (E_linew(i,j),i=1,numw(j))
      write(14,121) (ist_linew(i,j),i=1,numw(j))
      write(14,121) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      close(14)

120   format(5(f10.6,1x))
121   format(5(i5,6x))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stop

      end