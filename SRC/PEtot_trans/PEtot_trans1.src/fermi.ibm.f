      subroutine Fermi(f_occ,S_occ,x,itypeFermi)
      implicit double precision(a-h,o-z)

cccccccccc   According to: 
cccc  (1) M. Methfessel, A.T. Paxton, Phys. Rev. B 40, 3616 (1989)
cccc  (2) G. Kresse, J. Furthmuller, Comp. Mat. Sci. 6, 15 (1996). 

      real*8 h(0:30)

      y=x

      if(itypeFermi.eq.1) then       !  Fermi-Diract smearing
      if(y.gt.100.d0) y=100.d0
      if(y.lt.-100.d0) y=-100.d0
      f_occ=1.d0/(dexp(y)+1.d0)
      S_occ=-((dabs(1.d0-f_occ)+1.D-30)*dlog(dabs(1.d0-f_occ)+
     & 1.D-30)+(dabs(f_occ)+1.D-30)*dlog(dabs(f_occ)+1.D-30))
      return
      endif

ccc Using M. Methfessel, A.T. Paxton, P.R.B, 40, 3616 (1989). 
      if(itypeFermi.ge.20.and.itypeFermi.lt.30) then 
      n=itypeFermi-20    ! the order of S_n function
      if(y.gt.9.d0) y=9.d0 
      if(y.lt.-9.d0) y=-9.d0
      f_occ=0.5d0*(1-erf(y))
      fg=exp(-y**2)

      h(0)=1.d0
      h(1)=2*y
      do i=1,2*n
      h(i+1)=2*y*h(i)-2*i*h(i-1)
      enddo

      a=1.d0/dsqrt(4*datan(1.d0))

      do i=1,n
      a=-a/(i*4)
      f_occ=f_occ+a*h(2*i-1)*fg
      enddo
      S_occ=0.5d0*a*h(2*n)*fg

      if(y.lt.-8.99) then
      f_occ=1.d0
      S_occ=0.d0
      endif

      return
      endif
ccccccccccccccccccccccccccccccccccccc

      write(6,*) "itypeFermi not defined, stop",itypeFermi
      stop
      end

