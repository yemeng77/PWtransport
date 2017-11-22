       program plot_all2plt

ccccc This program convert a unformatted dens (or vr), graph.R, or charge_out file 
ccccc into a formatted xxx.plt_f file for 3D isosurface plot 
ccccc using gOpenMol. This formatted .plt_f file needs to be
ccccc transformed into a unformatted .plt file using pltfile
ccccc within gopenmol (RUN), before be read in by the gopenmol-plot-contour
ccccc option. 

       implicit double precision (a-h,o-z)

       real*8 AL(3,3),ALI(3,3),AL0(3,3)
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

       call input_vr()

       write(6,*) "n1,n2,n3=",n1,n2,n3
       write(6,*) "*** current AL is ***"
       write(6,101) AL(1,1),AL(2,1),AL(3,1)
       write(6,101) AL(1,2),AL(2,2),AL(3,2)
       write(6,101) AL(1,3),AL(2,3),AL(3,3)
101    format(3(f12.5,1x))

************************************************
         write(6,*) "input 1: use AL as the box; 2: new box"
         read(5,*) iflag

         if(iflag.eq.1) then
         x10=0.d0
         x20=0.d0
         x30=0.d0
         AL0=AL
         else
         write(6,*)
     & "input the corner position of new box as x1,x2,x3 of AL"
         read(5,*) x10,x20,x30
         write(6,*) "input AL0(:,1) of the new box in a.u"
         read(5,*) AL0(1,1),AL0(2,1),AL0(3,1)
         write(6,*) "input AL0(:,2) of the new box in a.u"
         read(5,*) AL0(1,2),AL0(2,2),AL0(3,2)
         write(6,*) "input AL0(:,3) of the new box in a.u"
         read(5,*) AL0(1,3),AL0(2,3),AL0(3,3)
         endif

        write(6,*) 
     & "It assumes the box is orth. If AL is not, plot is distorted" 

        write(6,*) 
     & "to overlay the atom plot correctly, same box must be used as in
     & conf2xyz.f" 
        write(6,*) 
     & "and output in conf2xyz.f must use x1,x2,x3 option, not xyz"


       xmax=dsqrt(AL0(1,1)**2+AL0(2,1)**2+AL0(3,1)**2)
       ymax=dsqrt(AL0(1,2)**2+AL0(2,2)**2+AL0(3,2)**2)
       zmax=dsqrt(AL0(1,3)**2+AL0(2,3)**2+AL0(3,3)**2)

       xtmp=dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
       ytmp=dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
       ztmp=dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)
       ratio=max(n1/xtmp,n2/ytmp,n3/ztmp)

       write(6,*) "input m1,m2,m3 (the grid for the box)"
       if(iflag.eq.1) then
       write(6,*) "suggested m1,m2,m3: ", n1,n2,n3
       else
       m1=xmax*ratio
       m2=ymax*ratio
       m3=zmax*ratio
       write(6,*) "suggested m1,m2,m3: ", m1,m2,m3
       endif

       read(5,*) m1,m2,m3

       allocate(vr_xyz(m1,m2,m3))

       call get_ALI(AL,ALI)

       do k=1,m3
       x3=(k-1.d0)/m3
       do j=1,m2
       x2=(j-1.d0)/m2
       do i=1,m1
       x1=(i-1.d0)/m1

       x=AL0(1,1)*x1+AL0(1,2)*x2+AL0(1,3)*x3
       y=AL0(2,1)*x1+AL0(2,2)*x2+AL0(2,3)*x3
       z=AL0(3,1)*x1+AL0(3,2)*x2+AL0(3,3)*x3

       ai1=(ALI(1,1)*x+ALI(2,1)*y+ALI(3,1)*z+x10+10)*n1+1
       aj1=(ALI(1,2)*x+ALI(2,2)*y+ALI(3,2)*z+x20+10)*n2+1
       ak1=(ALI(1,3)*x+ALI(2,3)*y+ALI(3,3)*z+x30+10)*n3+1
       i1=ai1
       j1=aj1
       k1=ak1


       fi=i1+1-ai1
       fi2=1.d0-fi
       fj=j1+1-aj1
       fj2=1.d0-fj
       fk=k1+1-ak1
       fk2=1.d0-fk

       i1m=mod(i1-1+20*n1,n1)+1
       i1p=mod(i1+20*n1,n1)+1
       j1m=mod(j1-1+20*n2,n2)+1
       j1p=mod(j1+20*n2,n2)+1
       k1m=mod(k1-1+20*n3,n3)+1
       k1p=mod(k1+20*n3,n3)+1

       vr_xyz(i,j,k)=fi*fj*fk*vr(i1m,j1m,k1m)+
     &   fi2*fj*fk*vr(i1p,j1m,k1m)
     &  +fi*fj2*fk*vr(i1m,j1p,k1m)+fi*fj*fk2*vr(i1m,j1m,k1p)
     &  +fi2*fj2*fk*vr(i1p,j1p,k1m)+fi2*fj*fk2*vr(i1p,j1m,k1p)
     &  +fi*fj2*fk2*vr(i1m,j1p,k1p)+fi2*fj2*fk2*vr(i1p,j1p,k1p)


       enddo
       enddo
       enddo


       xmax=xmax*0.529177
       ymax=ymax*0.529177
       zmax=zmax*0.529177

       write(6,*) 
     & "output written in graph.plt_f, to be converted into" 
       write(6,*)
     & "graph.plt by running Pltfile(conversion) within gopenmol" 

       open(12,file='graph.plt_f')
       rewind(12)
       write(12,99) 3,3
       write(12,98) m3,m2,m1

99     format(2(i2,2x))
98     format(3(i5,2x))
       write(12,200) 0.,zmax,0.,ymax,0.,xmax
200    format(6(E13.5,1x))
       write(12,100) (((vr_xyz(i,j,k),i=1,m1),j=1,m2),k=1,m3)

100    format(6(E13.5,1x))
       close(12)


       deallocate(vr)

       stop

       contains 
************************************************

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

       write(6,*) "input the name of atom.config file for the AL inform"
       read(5,*) filename

       open(10,file=filename)
       rewind(10)
       read(10,*) num_atom
       read(10,*) AL(1,1),AL(2,1),AL(3,1)
       read(10,*) AL(1,2),AL(2,2),AL(3,2)
       read(10,*) AL(1,3),AL(2,3),AL(3,3)
       close(10)

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

      

