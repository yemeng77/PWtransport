      program find_evan_ym
cccccc Written by Meng Ye. Jul, 2018

      implicit double precision (a-h,o-z)

      parameter (nm=2000)
      parameter (nevan=2000)

      real*8 E_evan(nevan),E_evanC(nevan)
      integer ist_evan(nevan),ikpt_evan(nevan),iGX_evan(nevan),
     &     iband_evan(nevan)
      integer iband_max(nevan),iband_min(nevan)
      integer inum_max(nevan),inum_min(nevan)
      integer iflag_max(nevan),iflag_min(nevan)

      real*8 E_linew(nm,nm)
      integer numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)

      complex*16, allocatable, dimension (:,:,:) :: S,U
      complex*16, allocatable, dimension (:) :: cc1_kpt1,cc2_kpt1
      complex*16, allocatable, dimension (:) :: cc1_kpt2,cc2_kpt2
      complex*16, allocatable, dimension (:) :: tmp1,tmp2

      complex*16 cc

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
      close(14)

      open(11,file="wr.overlap",form="unformatted")
      rewind(11)
      read(11) nkpt,nst
      open(12,file="wr.rotation",form="unformatted")
      rewind(12)
      read(12) nkptt,nstt

      if(nkpt.ne.nkptt.or.nst.ne.nstt) then
        write(6,*) "nkpt or nst in wr.overlap and wr.rotation are not 
     &        same,",nkpt,nst,nkptt,nstt
        stop
      endif

      if((nkpt-1)*nintep+1.ne.nkptw) then
        write(6,*) "nkpt0, nkpt, nintep are not match,"
     &       ,nkpt,nkptw,nintep
        stop
      endif

      allocate(S(nst,nst,nkpt-1))
      allocate(U(nst,nst,nkpt-1))

      do ikpt=1,nkpt-1
        read(11) S(:,:,ikpt)
        do i=1,4
          read(11)
        enddo
        read(12) U(:,:,ikpt)
      enddo

      close(11)
      close(12)

cccccc first, the evanescent state at the Gamma point
      num_evan=0
      do j=1,nstw
      if(ikpt_linew(1,j).eq.1.and.numw(j).gt.2*nintep) then    ! Gamma point, forget about short lines
      dEdk1=2*(E_linew(nintep+1,j)-E_linew(1,j))   ! curvature at Gamma (assuming a slop=0)
      dEdk2=E_linew(1,j)+E_linew(2*nintep+1,j)-2*E_linew(nintep+1,j) ! curvature calculated at 2

      if(dabs(dEdk1).lt.5*abs(dEdk2)) then       ! the Gamma curv similar to 2 curv, indeed, slop=0
      
      num_evan=num_evan+1
      E_evanC(num_evan)=dEdk1     
      E_evan(num_evan)=E_linew(1,j)
cccc if E_evanC>0, then evan state is at [-infinity,E_evan]
cccc if E_evane<0, then evan state is at [E_evan,infinity]
      ist_evan(num_evan)=ist_linew(1,j)
      ikpt_evan(num_evan)=ikpt_linew(1,j)
      iGX_evan(num_evan)=0
      iband_evan(num_evan)=j
      endif
      endif
      enddo
cccccccccccccccccccccccccccccccccccc
cccccc second, the evanescent state at the X point
      do j=1,nstw
       numwt=numw(j)
      if(ikpt_linew(numwt,j).eq.nkptw.and.
     & numwt.gt.2*nintep) then    ! X-point, forget about short lines

      dEdk1=2*(E_linew(numwt-nintep,j)-E_linew(numwt,j))   ! curvature at X (assuming a slop=0)
      dEdk2=E_linew(numwt,j)+E_linew(numwt-2*nintep,j)
     &     -2*E_linew(numwt-nintep,j) ! curvature calculated at 2

      if(dabs(dEdk1).lt.5*abs(dEdk2)) then       ! the Gamma curv similar to 2 curv, indeed, slop=0

      num_evan=num_evan+1
      E_evanC(num_evan)=dEdk1     
      E_evan(num_evan)=E_linew(numwt,j)
      ist_evan(num_evan)=ist_linew(numwt,j)
      ikpt_evan(num_evan)=ikpt_linew(numwt,j)
      iGX_evan(num_evan)=1
      iband_evan(num_evan)=j
      endif
      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(14,file="E_line_W2K2",position="append")
      write(14,*) "*************************"
      write(14,*) num_evan," evanescent states at the Gamma and X point"
      do ii=1,num_evan
       write(14,606) ii,ikpt_evan(ii),ist_evan(ii),iGX_evan(ii),
     & iband_evan(ii), E_evan(ii),E_evanC(ii)
      enddo
606   format(5(i4,1x),2(E18.10,1x))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc now, find max-min points in the middle 
      num_evan=0

      num_max=0
      num_min=0
      iflag_max=0
      iflag_min=0

      do j=1,nstw
      do i=2,numw(j)-1

      if(E_linew(i,j).gt.E_linew(i-1,j).and.
     &   E_linew(i,j).ge.E_linew(i+1,j)) then
      num_max=num_max+1
      iband_max(num_max)=j
      inum_max(num_max)=i
      endif

ccccccccccccccccccccccccccccccccccc
      if(E_linew(i,j).le.E_linew(i-1,j).and.
     &   E_linew(i,j).lt.E_linew(i+1,j)) then
      num_min=num_min+1
      iband_min(num_min)=j
      inum_min(num_min)=i
      endif
      enddo
      enddo

      write(6,*) "The numbers of maximum and minimum points are",
     &           num_max,num_min
      
      kmax=floor(nkptw*0.05d0)
      tol=0.5d0

      allocate(cc1_kpt1(nst))
      allocate(cc2_kpt1(nst))
      allocate(cc1_kpt2(nst))
      allocate(cc2_kpt2(nst))
      allocate(tmp1(nst))
      allocate(tmp2(nst))

      do 2000 imax=1,num_max
        i1=inum_max(imax)
        j1=iband_max(imax)
        kpt1=ikpt_linew(i1,j1)
        E1=E_linew(i1,j1)
        do 1000 imin=1,num_min
          if(iflag_min(imin).eq.0) then
            i2=inum_min(imin)
            j2=iband_min(imin)
            kpt2=ikpt_linew(i2,j2)
            E2=E_linew(i2,j2)
            if(abs(kpt1-kpt2).le.kmax.and.E2.gt.E1) then
              i11=i1-kpt1+kpt2
              if(i11.gt.0.and.i11.le.numw(j1)) then
                E11=E_linew(i11,j1)
                if(E11.gt.E1) goto 1000
                dE=E2-E11
                if(dE.le.0.d0) goto 1000
                do k1=1,nkptw
                  if((i11-k1).le.0.or.(i2+k1).gt.numw(j2)) goto 1000
                  if(E_linew(i11-k1,j1).gt.E1.or.
     &               E_linew(i2+k1,j2).lt.E2) goto 1000
                  dE1=E_linew(i2+k1,j2)-E_linew(i11-k1,j1)
                  if(dE1.ge.4*dE) goto 301
                enddo
301             continue
                do k2=1,nkptw
                  if((i11+k2).gt.numw(j1).or.(i2-k2).le.0) goto 1000
                  if(E_linew(i11+k1,j1).gt.E1.or.
     &               E_linew(i2-k1,j2).lt.E2) goto 1000
                  dE2=E_linew(i2-k2,j2)-E_linew(i11+k1,j1)
                  if(dE2.ge.4*dE) goto 302
                enddo
302             continue

                ist1=ist_linew(i11-k1,j1)
                ist2=ist_linew(i2+k1,j2)
                ikpt1=ikpt_linew(i11-k1,j1)
                ikpt2=ikpt_linew(i2+k1,j2)
                call dot_product()
                if(cdabs(cc).lt.tol) goto 1000

                ist1=ist_linew(i11+k2,j1)
                ist2=ist_linew(i2-k2,j2)
                ikpt1=ikpt_linew(i11+k2,j1)
                ikpt2=ikpt_linew(i2-k2,j2)
                call dot_product()
                if(cdabs(cc).lt.tol) goto 1000
                
                iflag_max(imax)=imin
                iflag_min(imin)=imax
                num_evan=num_evan+1
                goto 2000
              endif
            endif
          endif
1000  continue
2000  continue

      deallocate(S)
      deallocate(U)

      write(14,*) "*************************"
      write(14,*) num_evan," evanescent states in the middle"
      ii=0
      do imax=1,num_max
      imin=iflag_max(imax)
      if(imin.ne.0) then
       ii=ii+1
       i1=inum_max(imax)
       j1=iband_max(imax)
       E1=E_linew(i1,j1)
       i2=inum_min(imin)
       j2=iband_min(imin)
       E2=E_linew(i2,j2)
       write(14,707) ii,i1,j1,i2,j2,E1,E2
      endif
      enddo
707   format(5(i4,1x),2(E18.10,1x))

      close(14)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stop
      contains

      subroutine dot_product()
      implicit double precision (a-h,o-z)

      open(10,file="wr.interp.new",form="unformatted")
      rewind(10)
      read(10)
      do ikpt=1,max(ikpt1,ikpt2)
        do ist=1,nst
          if(ikpt.eq.ikpt1.and.ist.eq.ist1) then
            read(10) cc1_kpt1
            read(10) cc2_kpt1
          else if(ikpt.eq.ikpt2.and.ist.eq.ist2) then
            read(10) cc1_kpt2
            read(10) cc2_kpt2
          else
            read(10)
            read(10)
          endif
        enddo
      enddo
      close(10)

      if(ikpt1.gt.ikpt2) then
        tmp1=cc1_kpt1
        tmp2=cc2_kpt1
        cc1_kpt1=cc1_kpt2
        cc2_kpt1=cc2_kpt2
        cc1_kpt2=tmp1
        cc2_kpt2=tmp2
        ikpt=ikpt1
        ikpt1=ikpt2
        ikpt2=ikpt
      endif

      ikpt1=(ikpt1-1)/nintep+1
      ikpt2=(ikpt2-1)/nintep+1

      do ikpt=ikpt1,ikpt2-1
        tmp1=dcmplx(0.d0,0.d0)
        tmp2=dcmplx(0.d0,0.d0)
        do m=1,nst
          do m1=1,nst
          tmp1(m)=tmp1(m)+U(m1,m,ikpt)*cc1_kpt1(m1)
          tmp2(m)=tmp2(m)+U(m1,m,ikpt+1)*cc2_kpt1(m1)
          enddo
        enddo
        cc1_kpt1=tmp1
        cc2_kpt1=tmp2
      enddo

      cc=dcmplx(0.d0,0.d0)
      do m1=1,nst
        do m2=1,nst
          cc=cc+cc1_kpt1(m1)*dconjg(cc2_kpt2(m2))*S(m1,m2,ikpt2)
     &      +cc2_kpt1(m1)*dconjg(cc1_kpt2(m2)*S(m2,m1,ikpt2))
        enddo
        cc=cc+cc1_kpt1(m1)*dconjg(cc1_kpt2(m1))
     &    +cc2_kpt1(m1)*dconjg(cc2_kpt2(m1))
      enddo

      return
      end subroutine dot_product


      end
