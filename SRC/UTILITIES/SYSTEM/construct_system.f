       program construct_system
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc This program constructs the system potential vr.system and xatom.system file
ccc from the central potential and electrode potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       implicit double precision (a-h,o-z)

       real*8 AL_cent(3,3),AL(3,3),AL_el(3,3),AL_tot(3,3)
       real*8 x_tot(3,2000),x_cent(3,1000),x_el(3,200)
       integer iat_tot(2000),iat_cent(1000),iat_el(200)

       real*8, allocatable, dimension (:,:,:) :: vr_cent,vr_el
       real*8, allocatable, dimension(:) :: vr_tmp

       character*30 fileV1, fileV2, fileXat1, fileXat2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       open(12,file="construct.input")
       rewind(12)
       read(12,*) dV
       read(12,*) ipL, ipR
       read(12,*) nnL, nnR, iadd_L,dcut
       read(12,*) fileV1,fileV2      ! central, electrode
       read(12,*) fileXat1,fileXat2  ! central, electrode
       read(12,*) nnodes_sys         ! new nnodes for system (to allow more nodes)
       close(12)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       write(6,*) "Reading central region atomic configuration"
       open(11,file=fileXat1)
       rewind(11)
       read(11,*) natom_cent
       read(11,*) AL_cent(1,1),AL_cent(2,1),AL_cent(3,1)
       read(11,*) AL_cent(1,2),AL_cent(2,2),AL_cent(3,2)
       read(11,*) AL_cent(1,3),AL_cent(2,3),AL_cent(3,3)
       do i=1,natom_cent
       read(11,*) iat_cent(i),x_cent(1,i),x_cent(2,i),x_cent(3,i)
       enddo
       close(11)

       write(6,*) "Reading central region system potential"
       open(11,file=fileV1,form="unformatted")
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
       write(6,*) "AL.ne.AL_cent, stop"
       write(6,*) AL_cent
       write(6,*) AL
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
100    format(3(f10.6,1x))
       n1c=n1
       n2c=n2
       n3c=n3
       nnodes_c=nnodes
       deallocate(vr_tmp)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(6,*) "Reading the electrode atomic configuration"
       open(11,file=fileXat2)
       rewind(11)
       read(11,*) natom_el
       read(11,*) AL_el(1,1),AL_el(2,1),AL_el(3,1)
       read(11,*) AL_el(1,2),AL_el(2,2),AL_el(3,2)
       read(11,*) AL_el(1,3),AL_el(2,3),AL_el(3,3)
       do i=1,natom_el
       read(11,*) iat_el(i),x_el(1,i),x_el(2,i),x_el(3,i)
       enddo
       close(11)

       write(6,*) "Reading the electrode system potential"
       open(11,file=fileV2,form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) AL

       diff=0.d0
       do j=1,3
       do i=1,3
       diff=diff+abs(AL(i,j)-AL_el(i,j))
       enddo
       enddo
       if(diff.gt.1.D-5) then
       write(6,*) "AL.ne.AL_el, stop"
       write(6,*) AL_el
       write(6,*) AL
       stop
       endif

       diff=0.d0
       do i=1,3
       diff=diff+abs(AL_el(i,1)-AL_cent(i,1)*n1/n1c)
       diff=diff+abs(AL_el(i,2)-AL_cent(i,2))
       diff=diff+abs(AL_el(i,3)-AL_cent(i,3))
       enddo
       if(diff.gt.1.D-3) then
       write(6,*) "AL_el.ne.AL_cent,y,z, stop"
       write(6,*) AL_cent
       write(6,*) AL_el
       stop
       endif

       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr_el(n1,n2,n3))

       do iread=1,nnodes
       read(11) vr_tmp

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr_el(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       close(11)
       n1e=n1
       n2e=n2
       n3e=n3
ccccccccccccccccccc
       if(n2e.ne.n2c.or.n3e.ne.n3c) then
       write(6,*) "n2e,n2c,n3e,n3c not the same",
     &  n2e,n2c, n3e,n3c
       stop
       endif

       deallocate(vr_tmp)


ccc Bias voltage in atomic units
        dV = dV/27.211396d0
cccccc in the construction, the left electrode will be
cccccc lowered by dV/2, while the right electrode will be raised by dV/2
cccccc and the central part will be shifted to match the left and right
cccccc ipL-1 of central is the n1e of electrode
cccccc ipR+1 of central is the 1 of electrode
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        ave_el_n1e=0.d0
        ave_el_1=0.d0
        ave_ipLm1=0.d0
        ave_ipRp1=0.d0
        do j=1,n2
        do k=1,n3
        ave_el_n1e=ave_el_n1e+vr_el(n1e,j,k)
        ave_el_1=ave_el_1+vr_el(1,j,k)
        ave_ipLm1=ave_ipLm1+vr_cent(ipL-1,j,k)
        ave_ipRp1=ave_ipRp1+vr_cent(ipR+1,j,k)
	enddo
	enddo
        ave_el_n1e=ave_el_n1e/(n2*n3)
        ave_el_1=ave_el_1/(n2*n3)
        ave_ipLm1=ave_ipLm1/(n2*n3)
        ave_ipRp1=ave_ipRp1/(n2*n3)

        shift=(ave_el_n1e-dV/2+ave_el_1+dV/2-
     &   ave_ipLm1-ave_ipRp1)/2

        err=abs((ave_el_n1e-dV/2)-(ave_ipLm1+shift))
     &   +abs((ave_el_1+dV/2)-(ave_ipRp1+shift))

        write(6,*) "connection err (eV)", err*27.211396

	
cccc iadd_L is an extension of the connection at the left end. 
cccc The purpose is to get a proper n1 in order to do FFT
cccc it can be positive or negative, and should be smaller than n1e 

        n1=(nnL+nnR)*n1e+ipR-ipL+1+iadd_L
        write(6,*) "n1,n2,n3",n1,n2,n3

        AL_tot(:,1)=AL_el(:,1)*n1/n1e
        AL_tot(:,2)=AL_el(:,2)
        AL_tot(:,3)=AL_el(:,3)

cccccccccccccccccccccccccccccccccccccccccc
!       nnodes=nnodes_c
       nnodes=nnodes_sys      ! new nnodes for system (to allow more nodes)
       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))

       write(6,*) "output potential in vr.system"

       open(11,file="vr.system",form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL_tot

       do iread=1,nnodes

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       if(i.le.n1e*nnL+iadd_L) then
       i1=mod(i-1-iadd_L+2*n1e,n1e)+1
       vr_tmp(ii)=vr_el(i1,j,k)-dV/2
       endif

       if(i.gt.n1e*nnL+iadd_L.and.
     &  i.le.n1e*nnL+iadd_L+ipR-ipL+1) then
       i1=i-n1e*nnL-iadd_L+ipL-1
       vr_tmp(ii)=vr_cent(i1,j,k)+shift
       endif

       if(i.gt.n1e*nnL+iadd_L+ipR-ipL+1) then
       i1=i-(n1e*nnL+iadd_L+ipR-ipL+1)
       i1=mod(i1-1,n1e)+1
       vr_tmp(ii)=vr_el(i1,j,k)+dV/2
       endif

       enddo
       write(11) vr_tmp
       enddo
       close(11)
cccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccc

       do i=1,natom_el 
       if(x_el(1,i).gt.1.d0-1.D-8) x_el(1,i)=x_el(1,i)-1.d0
       if(x_el(1,i).lt.-1.D-8) x_el(1,i)=x_el(1,i)+1.d0
       enddo

       do i=1,natom_cent
       if(x_cent(1,i).gt.1.d0-1.D-8) 
     &          x_cent(1,i)=x_cent(1,i)-1.d0
       if(x_cent(1,i).lt.-1.D-8) x_cent(1,i)=x_cent(1,i)+1.d0
       enddo


       xcL=(ipL-1.d0)/n1c-1.D-5
       xcR=(ipR-1.d0)/n1c-1.D-5
       xcLL=(ipL-1.d0)/n1c
       xcRR=(ipR-1.d0)/n1c


       xmax=0.d0
       do i=1,natom_el
       if(x_el(1,i).gt.xmax) xmax=x_el(1,i)
       enddo

       dis=(1-xmax)*sqrt(AL_el(1,1)**2+AL_el(2,1)**2+AL_el(3,1)**2)
       dcut=-dis+dcut
       xcut=dcut/sqrt(AL_tot(1,1)**2+AL_tot(2,1)**2+AL_tot(3,1)**2)


       num=0
       if(iadd_L.eq.0) then
       do ii=1,nnL    
       do i=1,natom_el
       xtmp=(x_el(1,i)*n1e+n1e*(ii-1))*1.d0/n1
       num=num+1
       x_tot(1,num)=xtmp
       x_tot(2,num)=x_el(2,i)
       x_tot(3,num)=x_el(3,i)
       iat_tot(num)=iat_el(i)
       enddo
       enddo
       else
       do ii=0,nnL     ! start from 0, instead of 1 for the possible iadd_L
       do i=1,natom_el
       xtmp=(x_el(1,i)*n1e+n1e*(ii-1)+iadd_L)*1.d0/n1
       if(xtmp.ge.xcut) then
       num=num+1
       x_tot(1,num)=xtmp
       x_tot(2,num)=x_el(2,i)
       x_tot(3,num)=x_el(3,i)
       iat_tot(num)=iat_el(i)
       endif
       enddo
       enddo
       endif



       do i=1,natom_cent
       if(x_cent(1,i).gt.xcL.and.x_cent(1,i).le.xcR) then
       num=num+1
       xtmp=((x_cent(1,i)-xcLL)*n1c+n1e*nnL+iadd_L)*1.d0/n1
       x_tot(1,num)=xtmp
       x_tot(2,num)=x_cent(2,i)
       x_tot(3,num)=x_cent(3,i)
       iat_tot(num)=iat_cent(i)
       endif
       enddo
       

       do ii=1,nnR
       do i=1,natom_el
       num=num+1
       xtmp=(n1e*nnL+ipR-ipL+1+x_el(1,i)*n1e+
     &                    n1e*(ii-1)+iadd_L)*1.d0/n1
       x_tot(1,num)=xtmp
       x_tot(2,num)=x_el(2,i)
       x_tot(3,num)=x_el(3,i)
       iat_tot(num)=iat_el(i)
       enddo
       enddo
       
     
       write(6,*) "Out put atomic config xatom.system ",num

       open(11,file="xatom.system")
       rewind(11)

       write(11,*) num
       write(11,200) AL_tot(1,1),AL_tot(2,1),AL_tot(3,1)
       write(11,200) AL_tot(1,2),AL_tot(2,2),AL_tot(3,2)
       write(11,200) AL_tot(1,3),AL_tot(2,3),AL_tot(3,3)
       do i=1,num
       write(11,300) iat_tot(i),x_tot(1,i),x_tot(2,i),x_tot(3,i),0,0,0
       enddo
       close(11)

200    format(3(f20.13,1x))
300    format(i4,2x,3(f14.10,1x),2x,3(i2,1x))

       stop
       end

