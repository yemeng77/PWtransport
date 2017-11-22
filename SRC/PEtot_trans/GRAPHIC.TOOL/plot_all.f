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
       write(6,*) "4: motif (from gen_motall)"
       read(5,*) iflag

       call input_vr( )

       write(6,*) "n1,n2,n3=",n1,n2,n3

*******************************  
**** output  the density 
cccccccccccccccc
cc       test
       write(6,*) "point(1); line(2); splot(3)"
       read(5,*) iflag

       if(iflag.eq.1) then
700    continue
       write(6,*) "input i,j,k"
       read(5,*) i,j,k
       write(6,*) i,j,k,vr(i,j,k)
       goto 700
       endif

       if(iflag.eq.2) then
600    continue
       write(6,*) "input i_orient(1,2,3) and ijk1,iik2"
       read(5,*) i_orient,ijk1,ijk2
       open(12,file="graph.plot")
       rewind(12)
       if(i_orient.eq.1) then
       do i=1,n1
       write(12,101) i, vr(i,ijk1,ijk2)
       enddo
       endif
       if(i_orient.eq.2) then
       do i=1,n2
       write(12,101) i, vr(ijk1,i,ijk2)
       enddo
       endif
       if(i_orient.eq.3) then
       do i=1,n3
       write(12,101) i, vr(ijk1,ijk2,i)
       enddo
       endif
       close(12)
       goto 600
       endif
101    format(i5,E12.5)
       
       
300    continue
       write(6,*) "input i_orient(1,2,3) and ijkth"
       read(5,*) i_orient,ijkth

       if(i_orient.eq.1) then
       open(12,file='graph.plot')
       rewind(12)
       do k=1,n3
       do j=1,n2
       write(12,100) vr(ijkth,j,k)
       enddo
       write(12,100)
       enddo
       close(12)
       endif

       if(i_orient.eq.2) then
       open(12,file='graph.plot')
       rewind(12)
       do k=1,n3
       do i=1,n1
       write(12,100) vr(i,ijkth,k)
       enddo
       write(12,100)
       enddo
       close(12)
       endif

       if(i_orient.eq.3) then
       open(12,file='graph.plot')
       rewind(12)
       do j=1,n2
       do i=1,n1
       write(12,100) vr(i,j,ijkth)
       enddo
       write(12,100)
       enddo
       close(12)
       endif


       goto 300
100    format(E14.5)

       stop

       contains 

       subroutine input_vr()
       implicit double precision (a-h,o-z)
       parameter (m=80)
       real*4 dens_m(-m:m,-m:m,-m:m)
       integer iat_neigh(10)
       real*8 x_neigh(3,10)
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
       write(6,*) AL(1,1),AL(2,1),AL(3,1)
       write(6,*) AL(1,2),AL(2,2),AL(3,2)
       write(6,*) AL(1,3),AL(2,3),AL(3,3)

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
**********************************************************

       if(iflag.eq.4) then
       write(6,*) "input the name of motif file to be plot"
       read(5,*) filename
       open(11,file=filename,form="unformatted")
       rewind(11)
       read(11) num_neigh, rad_box
       read(11) AL
       do ia=1,num_neigh
          read(11) iat_neigh(ia),x_neigh(1,ia),x_neigh(2,ia),
     &        x_neigh(3,ia)
       enddo
       read(11) mb,sum2
       write(6,*) "input the nth motif (e.g, for emotif)"
       read(5,*) nth
       do i=1,nth-1
       read(11)
       enddo
       read(11) dens_m
       close(11)

       n1=2*m+1
       n2=2*m+1
       n3=2*m+1
       allocate(vr(n1,n2,n3))

       do k=-m,m
          do j=-m,m
             do i=-m,m
                vr(i+m+1,j+m+1,k+m+1)=dble(dens_m(i,j,k))
             enddo
          enddo
       enddo

       endif     ! end of motif input
************************************************************
       return
       end subroutine input_vr

       end

      

