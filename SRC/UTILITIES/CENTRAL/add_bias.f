       program add_bias
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc This program adds a bias voltage potential to an existing
ccc system potential. In particular, we have added a potential 
ccc V/2*sin(pi*z/L) to get the  potential for a given bias V. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       implicit double precision (a-h,o-z)

       real*8 AL_cent(3,3),AL(3,3),x_cent(3,200)
       integer iat_cent(200)
       real*8, allocatable, dimension (:,:,:) :: vr_cent
       real*8, allocatable, dimension(:) :: vr_tmp

cccccccccccccccccccccccccccccccccccccccccccc
       open(11,file="xatom.system_vac")
       rewind(11)
       read(11,*) natom_cent
       read(11,*) AL_cent(1,1),AL_cent(2,1),AL_cent(3,1)
       read(11,*) AL_cent(1,2),AL_cent(2,2),AL_cent(3,2)
       read(11,*) AL_cent(1,3),AL_cent(2,3),AL_cent(3,3)
       do i=1,natom_cent
       read(11,*) iat_cent(i),x_cent(1,i),x_cent(2,i),x_cent(3,i)
       enddo
       close(11)

       open(11,file="vr.system_vac",form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) AL

       diff=0.d0
       do j=1,3
       do i=1,3
       diff=diff+abs(AL(i,j)-AL_cent(i,j))
       enddo
       enddo
       if(diff.gt.1.D-5) then
       write(6,*) "AL.ne.AL_cent, stop", diff
       stop
       endif


       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr_cent(n1,n2,n3))

       do iread=1,nnodes
       read(11) vr_tmp

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr_cent(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       close(11)
       write(6,*) "n1,n2,n3",n1,n2,n3
100    format(3(f10.6,1x))
       n1c=n1
       deallocate(vr_tmp)

c Represent planar averaged of vr.Chain
        open(53,file="graphCenter.plot")
        rewind(53)
        do i=1,n1
        sum=0.d0
        do j=1,n2
        do k=1,n3
        sum = sum +(vr_cent(i,j,k))
        enddo
        enddo
	xfrac=(i-1.0)/(n1-1.0)
 	xcent=(xfrac-0.5)*AL_cent(1,1)
        write(53,*) i,xcent,sum/(n2*n3)
        enddo
	close(53)
cccccccccccccccccccccccccccccccccccccccccccccccccc
       pi=4*datan(1.d0)
       write(6,*) "input bias voltage dV (eV)"
       read(5,*) dV
       dV=dV/27.211396d0
       do i=1,n1c
       x=(i-n1c/2.d0-1.d0)*pi/n1c
       ddV=dV*dsin(x)/2
       do k=1,n3
       do j=1,n2
       vr_cent(i,j,k)=vr_cent(i,j,k)+ddV
       enddo
       enddo
       enddo
cccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccc
       nnodes=64
       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))

	write(6,*)  "New potential output in vr.system_new"

       open(11,file="vr.system_new",form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL

       do iread=1,nnodes

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr_tmp(ii)=vr_cent(i,j,k)

       enddo
       write(11) vr_tmp
       enddo
       close(11)
cccccccccccccccccccccccccccccccccccccccccccc
200    format(3(f15.8,1x))
300    format(i4,2x,3(f14.10,1x),2x,3(i2,1x))

c Represent planar averaged of vr_system
        open(51,file="graphVr.plot")
        rewind(51)
        do i=1,n1
        sum=0.d0
        do j=1,n2
        do k=1,n3
        sum = sum +(vr_cent(i,j,k))
        enddo
        enddo
	xfrac=(i-1.0)/(n1-1.0)
        xcent=(xfrac-0.5)*AL(1,1)
        write(51,*) i,xcent,sum/(n2*n3)
        enddo
        close(51)



       stop
       end

