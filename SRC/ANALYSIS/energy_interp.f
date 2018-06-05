      subroutine energy_interp(ak,nst,E_st,S_m,idirect,
     &  ist1,ist2,Ew,eigen,vec)

      implicit double precision (a-h,o-z)

      real*8 ak,Ew,eigen
      integer nst,is1,ist2,idirect,ist
      complex*16 S_m(nst,nst),hh(nst,nst)
      complex*16 vec(nst)
      real*8 EE(nst)
      real*8 E_st(nst,2)

      complex*16, allocatable, dimension(:) :: work
      real*8, allocatable, dimension(:):: rwork

      hh=dcmplx(0.d0,0.d0)

      if(idirect.eq.0) then
        do m1=1,nst
        do m2=1,nst
          do m3=1,nst
          hh(m1,m2)=hh(m1,m2)+ak*E_st(m3,2)*S_m(m1,m3)*
     &      dconjg(S_m(m2,m3))
          enddo
        enddo
          hh(m1,m1)=hh(m1,m1)+(1-ak)*E_st(m1,1)
        enddo
      else
        do m1=1,nst
        do m2=m1,nst
          do m3=1,nst
          hh(m1,m2)=hh(m1,m2)+(1-ak)*E_st(m3,1)*S_m(m1,m3)*
     &      dconjg(S_m(m2,m3))
          enddo
        enddo
          hh(m1,m1)=hh(m1,m1)+ak*E_st(m1,2)
        enddo
      endif

      lwork=2*nst-1
      allocate(work(lwork))
      allocate(rwork(3*nst-2))
      call zheev('V','U',nst,hh,nst,EE,work,lwork,rwork,info)
      deallocate(work)
      deallocate(rwork)

      if(info.ne.0) then
        write(6,*) "Somthing is wrong with
     &   diagonalization of hh when ak=",ak
        stop
      endif

      err_min=1.D6
      do m=ist1,ist2
        if(dabs(EE(m)-Ew).lt.dabs(err_min)) then
          err_min=dabs(EE(m)-Ew)
          ist=m
        endif
      enddo

      eigen=EE(ist)
      vec(:)=hh(:,ist)

      return
      end