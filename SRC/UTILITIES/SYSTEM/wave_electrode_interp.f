      subroutine wave_electrode_interp(ikpt_st,i_st,uc,m1,m2,m3,
     & fileh,nintep)

      implicit double precision (a-h,o-z)
      complex*16, allocatable, dimension (:) :: cc1,cc2
      real*8, allocatable, dimension (:) :: ucR,ucI
      complex*16 uc(m1,m2,m3)
      integer nintep
      character*7 fileh
      character*20 filename
***********************************************************
      if(nintep.eq.1) then
        ikpt_st1=ikpt_st
      else
        open(12,file="wr.interp.new",form="unformatted")
        rewind(12)
        read(12) nkpt,nst,nkpt0,nintp
        if(nintp.ne.nintep) then
          write(6,*) "nintep .ne. nintp, stop",nintep,nintp
          stop
        endif
        allocate(cc1(nst))
        allocate(cc2(nst))
        do ikpt=1,ikpt_st-1
        do ist=1,nst
        read(12)
        read(12)
        enddo
        enddo
        do ist=1,i_st
        read(12) (cc1(i),i=1,nst)
        read(12) (cc2(i),i=1,nst)
        enddo
        close(12)
        ikpt_st1=(ikpt_st-1)/nintep+1
        ikpt_st2=ikpt_st-1-(ikpt_st1-1)*nintep
      endif

      ikpt_100=ikpt_st1/100
      ikpt_10=(ikpt_st1-ikpt_100*100)/10
      ikpt_1=ikpt_st1-ikpt_100*100-ikpt_10*10
      filename=trim(fileh)//char(48+ikpt_100)//char(48+ikpt_10)//
     &   char(48+ikpt_1)

      open(10,file=filename,form="unformatted")
      rewind(10)
      read(10) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
      read(10) AL
      if(n1.ne.m1.or.n2.ne.m2.or.n3.ne.m3) then
      write(6,*) "n1,n2,n3.ne.m1,m2,m3,stop",n1,n2,n3,m1,m2,m3
      stop
      endif
      nr=n1*n2*n3
      nr_n=nr/nnodes

      allocate(ucR(nr_n))
      allocate(ucI(nr_n))

      if(nintep.eq.1) then
        do ist=1,i_st-1
        do iread=1,nnodes
        read(10)       ! jump over this record, much faster. 
        enddo
        enddo
        do iread=1,nnodes
        read(10) (ucR(i),i=1,nr_n),(ucI(i),i=1,nr_n)
        do ii=1,nr_n
        jj=ii+(iread-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        uc(i,j,k)=dcmplx(ucR(ii),ucI(ii))
        enddo
        enddo
        deallocate(ucR)
        deallocate(ucI)
        close(10)
      else
        uc=dcmplx(0.d0,0.d0)
        do ist=1,nst
        do iread=1,nnodes
        read(10) (ucR(i),i=1,nr_n),(ucI(i),i=1,nr_n)
        do ii=1,nr_n
        jj=ii+(iread-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        uc(i,j,k)=uc(i,j,k)+cc1(ist)*dcmplx(ucR(ii),ucI(ii))
        enddo
        enddo
        enddo
        deallocate(ucR)
        deallocate(ucI)
        close(10)

        if(ikpt_st2.ne.0) then
          ikpt_st1=ikpt_st1+1
          ikpt_100=ikpt_st1/100
          ikpt_10=(ikpt_st1-ikpt_100*100)/10
          ikpt_1=ikpt_st1-ikpt_100*100-ikpt_10*10
          filename=trim(fileh)//char(48+ikpt_100)//char(48+ikpt_10)//
     &    char(48+ikpt_1)
          open(11,file=filename,form="unformatted")
          rewind(11)
          read(11) n1,n2,n3,nnodes,ispin_i,ispin_f,iw_i,iw_f
          read(11) AL
          if(n1.ne.m1.or.n2.ne.m2.or.n3.ne.m3) then
            write(6,*) "n1,n2,n3.ne.m1,m2,m3,stop",n1,n2,n3,m1,m2,m3
            stop
          endif
          nr=n1*n2*n3
          nr_n=nr/nnodes
          allocate(ucR(nr_n))
          allocate(ucI(nr_n))
          do ist=1,nst
          do iread=1,nnodes
          read(11) (ucR(i),i=1,nr_n),(ucI(i),i=1,nr_n)
          do ii=1,nr_n
          jj=ii+(iread-1)*nr_n
          i=(jj-1)/(n2*n3)+1
          j=(jj-1-(i-1)*n2*n3)/n3+1
          k=jj-(i-1)*n2*n3-(j-1)*n3
          uc(i,j,k)=uc(i,j,k)+cc2(ist)*dcmplx(ucR(ii),ucI(ii))
          enddo
          enddo
          enddo
          deallocate(ucR)
          deallocate(ucI)
          close(11)
        endif
        deallocate(cc1)
        deallocate(cc2)
      endif 

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
