       program gen_Vext

ccccc This program plot all types of charge densities. 

       implicit double precision (a-h,o-z)

       real*8 AL(3,3),ALI(3,3)

       character*20 filename

       real*8, allocatable, dimension (:,:,:) :: vr
       real*8, allocatable, dimension (:,:,:) :: vr_xyz
       real*8, allocatable, dimension(:) :: vr_tmp

       write(6,*) "input the name of atom.config file"
       read(5,*) filename
       open(10,file=filename)
       rewind(10)
       read(10,*) nat
       read(10,*) AL(1,1),AL(2,1),AL(3,1)
       read(10,*) AL(1,2),AL(2,2),AL(3,2)
       read(10,*) AL(1,3),AL(2,3),AL(3,3)
       close(10)
       write(6,*) "--------- AL ----------"
       write(6,100) AL(1,1),AL(2,1),AL(3,1)
       write(6,100) AL(1,2),AL(2,2),AL(3,2)
       write(6,100) AL(1,3),AL(2,3),AL(3,3)
       write(6,*) "--------- AL ----------"
100    format(3(f12.6,1x))
       write(6,*) "input n1,n2,n3,nnodes"
       read(6,*) n1,n2,n3,nnodes
       write(6,*) "vext= v_jump*[(x-x1)+0.5] for x < x1"
       write(6,*) "vext= v_jump*[(x-x1)-0.5] for x < x1"
       write(6,*) "input v_jump(a.u.), x1(0:1)"
       read(6,*) v_jump,x1


       open(11,file="Vext_gen",form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL

       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr(n1,n2,n3))

       do iread=1,nnodes

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       x=(i-1.d0)/n1
       if(x.lt.x1) then
       vr_tmp(ii)=v_jump*(x-x1+0.5d0)
       else
       vr_tmp(ii)=v_jump*(x-x1-0.5d0)
       endif
       enddo

       write(11) vr_tmp

       enddo

       close(11)
       write(6,*) "potential written in Vext_gen"

       stop
       end
