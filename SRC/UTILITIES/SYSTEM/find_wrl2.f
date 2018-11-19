      program find_wrl2
cccccc AL is the lattice of the large cell

      implicit double precision (a-h,o-z)

      complex*16, allocatable, dimension (:,:,:) :: ur
      complex*16, allocatable, dimension (:) :: ur_tmp
      real*8, allocatable, dimension (:) :: urR,urI

      real*8, allocatable, dimension (:) :: err_st,E_st

      real*8 AL(3,3)

      complex*16 cc

      open(10,file="find_wrl.input")
      rewind(10)
      read(10,*) n1,n2,n3,n1w,nnodes
      read(10,*) Ew,dV,dE_evan,dk_same
      read(10,*) mst,num0
      close(10)

      open(21,file="wr.out001",form="unformatted")
      rewind(21)
      read(21) n1t,n2t,n3t,nnodesw,ispin_i,ispin_f,iw_i,iw_f
      read(21) AL

      if(n1.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
      write(6,*) "n1,n2,n3 changed in wr.new.***, stop",
     & n1,n2,n3,n1t,n2t,n3t
      stop
      endif

      nr=n1*n2*n3
      nr_n=nr/nnodes
      nr_nw=nr/nnodesw

      allocate(ur(n1,n2,n3))
      allocate(ur_tmp(nr_n))
      allocate(urR(nr_nw))
      allocate(urI(nr_nw))
      allocate(err_st(mst))
      allocate(E_st(mst))

      open(16,file="eigen_wg0")
      rewind(16)
      do i=1,4
        read(16,*)
      enddo
      read(16,*) (err_st(i),i=1,mst)
      read(16,*)
      read(16,*) (E_st(i),i=1,mst)
      close(16)

      dEmax=1.D-6*27.211396d0

      istate=num0

      do 4000 iw=iw_i,iw_f
      
      dE=a=dabs(E_st(iw)-Ew)
      if(dE.gt.dEmax) then
        do iread=1,nnodesw
          read(21)
        enddo
        goto 4000
      endif

      do iread=1,nnodesw
        read(21) (urR(i),i=1,nr_nw), (urI(i),i=1,nr_nw)
        do ii=1,nr_nw
          jj=ii+(iread-1)*nr_nw
          i=(jj-1)/(n2*n3)+1
          j=(jj-1-(i-1)*n2*n3)/n3+1
          k=jj-(i-1)*n2*n3-(j-1)*n3
          ur(i,j,k)=dcmplx(urR(ii),urI(ii))
        enddo
      enddo

      istate=istate+1
      nd=5
      pi=4*datan(1.d0)

      ii1=istate/10
      ii2=istate-ii1*10

      open(10,file="wr_test."//char(ii1+48)//char(ii2+48),
     &    form="unformatted")
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
      ur_tmp(ii)=ur(i,j,k)*fac
      endif

      if(i.le.nd) then
      iw=i+n1w
      fac=(1.d0-dsin((iw-n1w-1)*pi/2.d0/nd))/2
      ur_tmp(ii)=ur(i,j,k)*fac
      endif

      if(i.ge.n1-n1w-nd+2.and.i.le.n1-n1w) then
      iw=i-n1+n1w
      fac=(1.d0+dsin((iw-1)*pi/2.d0/nd))/2
      ur_tmp(ii)=ur(i,j,k)*fac
      endif

      enddo

        write(10) ur_tmp
        enddo

       write(10) iw
       close(10)


4000   continue

      deallocate(ur)
      deallocate(ur_tmp)
      deallocate(urR)
      deallocate(urI)
      deallocate(err_st)
      deallocate(E_st)

ccccccccccccccccccccccccccccccccccc
      stop

      end
