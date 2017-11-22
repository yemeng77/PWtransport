      subroutine system_scfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccc input: x  (real)   [1,n] ! Note at the return, x has been destroyed
cccc output: y  (complex) [1,n/2+1]
      implicit none
      integer init,isign,n,isys,ntable,nwork
      real*8 scale,table(*),work(*)
      real*8 x(*),y(*)
ccccccc real to complex 1D fft

cccccc for T3E
c      if(init.eq.1) then
c      call scfft(0,n,scale,0,0,table,0,0)
c      else
c      call scfft(isign,n,scale,x,y,table,work,isys)
c      endif
cccccccccccccccccccccc
cccccc for IBM SP2, essl  
      if(init.eq.1) then
      call drcft(1,0,0,0,0,
     &    n,1,-isign,scale,table,ntable,0,0)
      else
      call drcft(0,x,0,y,0,n,1,-isign,scale,
     &    table,ntable,work,nwork) 
      endif
cccccccccccccccccccccccccccccccccccc
cccccc for Intel KML
c       if(init.eq.1) then
c       call dzfft1d(y,n,0,table)
c       else
c       do i=1,n
c       y(i)=x(i)*scale
c       enddo
c       call dzfft1d(y,n,-1,table)    ! or -isign
ccccc dzfft1d(y,n,1,table) and dzfft1d(y,n,-1,table) are the same !
c       if(isign.eq.1) then
c       do i=2,n+2,2
c       y(i)=-y(i)
c       enddo
c       endif
c       endif
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

 
      
     

