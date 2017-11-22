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
       if(init.eq.1) then
       call dzfft(0,n,y,table,info)
       else
       call dzfft(1,n,x,table,info)      ! Note, input x will be destroyed, but it is okay in PEtot
          fact1=dsqrt(n*1.d0)*scale
          fact2=-isign*fact1
       do i=0,n/2
         y(2*i+1)=fact1*x(i+1)
        enddo
        do i=1,(n-1)/2
         y(2*i+2)=fact2*x(n-i+1)
        enddo
         y(2)=0.d0
         if(mod(n,2).eq.0) y(n+2)=0.d0
       endif
ccccccccccccccccccccccccccccccccccccccccc



      return
      end

 
      
     

