      program gen_mask
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This programs generates a mask function later used in the
c calculation of the self-consistent potential of the system
c Caution: this utility is valid only for cubic unit cells!
c          For other cells, modification should be made.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      real*8 AL(3,3)
      real*8, allocatable, dimension (:) ::  xp, yp, zp 
      real*8, allocatable, dimension (:,:) ::  xatom
      real*8, allocatable, dimension (:) :: fw_L, vr_tmp
      real*8, allocatable, dimension (:,:,:) :: w_L 
      integer  iatom(1000) 
      integer  imov_at(3,1000)
      character*30 filename

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Read input data file	
        open(18,file="mask.input")
        rewind(18)
	read(18,*) n1,n2,n3,nnodes
 	read(18,*) xpL,xpR
        read(18,*) filename
	close(18)

        open(10,file=filename)
        read(10,*) natom
        read(10,*) AL(1,1),AL(2,1),AL(3,1)
        read(10,*) AL(1,2),AL(2,2),AL(3,2)
        read(10,*) AL(1,3),AL(2,3),AL(3,3)
        close(10)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       allocate(fw_L(n1))

	iL = int(xpL*n1)+1
	iR = int(xpR*n1)+1

        do i=1,iL
        fw_L(i)=1.d0
        enddo

        do i=iR,n1
        fw_L(i)=0.d0
        enddo
        
        pi=4*datan(1.d0)

        do i=iL+1,iR-1
        x=(i-iL)*1.d0/(iR-iL)
        y=dcos(x*pi/2)**2
        if(x.le.0.5) then
        y=1.d0-y
        y=y**2
        y=1.d0-y/(2*0.5**2)
        else
        y=y**2/(2*0.5**2)
        endif
        fw_L(i)=y
        enddo

        open(12,file="mask.dat")
	do i = 1,n1
	write(12,*) i,fw_L(i)
	enddo
	close(12)

cccccccccccccccccccccccccccccccccccccccccc
       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))

	write(6,*) "Mask function is output in mask.central"

       open(11,file="mask.central",form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL

       do iread=1,nnodes

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr_tmp(ii)=fw_L(i)

       enddo
       write(11) vr_tmp
       enddo
       close(11)
cccccccccccccccccccccccccccccccccccccccccccc


	stop
	end









