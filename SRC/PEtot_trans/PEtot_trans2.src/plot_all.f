       program plot_all

ccccc This program plot all types of charge densities. 

       implicit double precision (a-h,o-z)

       real*8 AL(3,3),ALI(3,3)
       character*20 file1

       character*20 filename

       real*8, allocatable, dimension (:,:,:) :: vr
       real*8, allocatable, dimension (:,:,:) :: vr_xyz
       real*8, allocatable, dimension(:) :: vr_tmp

       write(6,*)  "input the type of input file"
       write(6,*) "1: graph.R (from Escan)" 
       write(6,*) "2: dens(or vr) for both PEtot and Escan"
       write(6,*) "3: charge_out (wavefunction from PEtot)"
       read(5,*) iflag

       call input_vr( )

       write(6,*) "n1,n2,n3=",n1,n2,n3

*******************************  
**** output  the density 

300    continue
c       write(6,*) "input kth"
c       read(5,*) kth
c       if(kth.lt.1.or.kth.gt.n3) stop
       write(6,*) "input ith"
       read(5,*) ith
       if(ith.lt.1.or.ith.gt.n1) stop


       open(12,file='graph.plot')
       rewind(12)
       do k=1,n3
       do j=1,n2
       write(12,100) vr(ith,j,k)
       enddo
       write(12,100)
       enddo
       close(12)

       goto 300
100    format(E14.5)

       stop

       contains 

       subroutine input_vr()
       implicit double precision (a-h,o-z)
*************************************************************
*************************************************************
       if(iflag.eq.1) then
       write(6,*) "input the graph.R (Escan) style file name"
       read(5,*) filename

       write(6,*)  "input the n1,n2,n3 in the graph.R file"
       read(5,*) n1,n2,n3

       allocate(vr(n1,n2,n3))

       write(6,*) "input the numth of the state to be plotted" 
       read(5,*) istate

       open(17,file=filename,form="unformatted")
       rewind(17) 
       do i=1,istate
       read(17) vr
       enddo
       close(17)

c       write(6,*) "input the name of atom.config file for the AL inform"
c       read(5,*) filename

c       open(10,file=filename)
c       rewind(10)
c       read(10,*) num_atom
c       read(10,*) AL(1,1),AL(2,1),AL(3,1)
c       read(10,*) AL(1,2),AL(2,2),AL(3,2)
c       read(10,*) AL(1,3),AL(2,3),AL(3,3)
c       close(10)

       endif    ! end for graph.R style input
*************************************************************
**********************************************************

       if(iflag.eq.2) then
       write(6,*) "input the name of dens (vr) file to be plot"
       read(5,*) filename
       open(11,file=filename,form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) AL

       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr(n1,n2,n3))

       do iread=1,nnodes
       read(11) vr_tmp

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       close(11)

       endif     ! end of dens(vr) input
********************************************************
********************************************************
       if(iflag.eq.3) then
       write(6,*) "input the name of charge_out (from PEtot)"
       read(5,*) file1
       open(11,file=file1)
       rewind(11)
       read(11,*) n1,n2,n3
       read(11,*) AL(1,1),AL(2,1),AL(3,1)
       read(11,*) AL(1,2),AL(2,2),AL(3,2)
       read(11,*) AL(1,3),AL(2,3),AL(3,3)
       read(11,*) nnodes

       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr(n1,n2,n3))

       do iread=1,nnodes
       read(11,*) (vr_tmp(i),i=1,nr_n)

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       close(11)
       endif       ! end of charge_out input
************************************************************
       return
       end subroutine input_vr

       end

      

