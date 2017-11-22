      subroutine ext(i)
c     
c  Stops program in case of errors or completion.
c
c  i is a stop parameter
c   000-099 main (0 is normal exit)
c   100-199 input
c   200-299 charge
c   300-399 vionic
c   400-499 velect
c   500-599 dsolv1
c   600-699 dsolv2 (including difnrl and difrel)
c   700-799 etotal
c   800-899 pseudo, pseudk, pseudt and pseudv
c
      if (i .ne. 0) write(6,10) i
 10   format('stop parameter =',i3)
      close (unit=1)
      close (unit=3)
      close (unit=5)
      close (unit=6)
      stop
      end
