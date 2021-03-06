      subroutine system_scfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccc input: x  (real)
cccc output: y  (complex)
      implicit none
      integer init,isign,n,isys,ntable,nwork,i,info
      real*8 scale,table(*),work(*)
      real*8 x(*),y(*),fact1,fact2
ccccccc real to complex 1D fft

cccccc for XT3
c      if(init.eq.1) then
c      call dzfft(0,n,scale,0,0,table,0,0)
c      else
c      call dzfft(isign,n,scale,x,y,table,work,isys)
c      endif
cccccc for T3E
c      if(init.eq.1) then
c      call scfft(0,n,scale,0,0,table,0,0)
c      else
c      call scfft(isign,n,scale,x,y,table,work,isys)
c      endif
cccccccccccccccccccccc
cccccc for IBM SP2, essl  
c      if(init.eq.1) then
c      call drcft(1,0,0,0,0,
c     &    n,1,-isign,scale,table,ntable,0,0)
c      else
c      call drcft(0,x,0,y,0,n,1,-isign,scale,
c     &    table,ntable,work,nwork) 
c      endif
cccccccccccccccccccccccccccccccccccc
cccccc for Intel KML
       if(init.eq.1) then
       call dzfft1d(y,n,0,table)
       else
       do i=1,n
       y(i)=x(i)*scale
       enddo
       call dzfft1d(y,n,-1,table)    ! or -isign
cccc dzfft1d(y,n,1,table) and dzfft1d(y,n,-1,table) are the same !
       if(isign.eq.1) then
       do i=2,n+2,2
       y(i)=-y(i)
       enddo
       endif
       endif
cccccccccccccccccccccccccccccccccccccccccccc
ccccc    for AMD ACML lib
c       if(init.eq.1) then
c       call dzfft(0,n,y,table,info)
c       else
c       call dzfft(1,n,x,table,info)      ! Note, input x will be destroyed, but it is okay in PEtot
c          fact1=dsqrt(n*1.d0)*scale
c          fact2=-isign*fact1
c       do i=0,n/2
c         y(2*i+1)=fact1*x(i+1)
c        enddo
c        do i=1,(n-1)/2
c         y(2*i+2)=fact2*x(n-i+1)
c        enddo
c         y(2)=0.d0
c         if(mod(n,2).eq.0) y(n+2)=0.d0
c       endif
ccccccccccccccccccccccccccccccccccccccccc



      return
      end

 
      
     

