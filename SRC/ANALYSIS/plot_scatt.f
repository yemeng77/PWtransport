       program plot_vr_rho_3D 

ccc A. Garcia-Lekue, (August 2008) ccc
ccccccccccccccccccccccccccccccccccccccccccccccc
ccc this short routine is used to plot
ccc density and potential; in 1D (graph.plot1D)
ccc or 2D  (graph.plot,graph.plot_2)
ccccccccccccccccccccccccccccccccccccccccccccccc

       implicit double precision (a-h,o-z)


       character*20 filename

       complex*16, allocatable, dimension (:,:,:) :: vr
       complex*16, allocatable, dimension(:) :: vr_tmp
       integer num_st(200)


       call input_vr( )

       write(6,*) "n1,n2,n3=",n1,n2,n3

*******************************  
**** output  the density 


       open(12,file='graph.plot1D')
       rewind(12)
	do i=1,n1
	sum=0.d0
	do j=1,n2
	do k=1,n3
	sum=sum+abs(vr(i,j,k))**2
	if ((j.eq.1).and.(k.eq.1)) then
 	write(46,*) i,j,k,vr(i,j,k)
	endif 
	enddo
	enddo
	sum=sum/(n2*n3)	
	write(12,*) i,sum
	enddo
	close(12)

	write(6,*) 'Input ith (from 1 to n1)',n1
	read(5,*) ith


       open(13,file='graph.plot3D_C')
       rewind(13)
       open(12,file='graph.plot3D')
       rewind(12)
       open(14,file='graph.plot3D_sq')
       rewind(14)
       do j=1,n2
       do k=1,n3
       j1=mod(j+n2/2-1,n2)+1
       k1=mod(k+n3/2-1,n3)+1
       R1=dreal(vr(ith,j,k))
       R2=dimag(vr(ith,j,k))
       R3=abs(vr(ith,j,k))**2
       if(abs(R1).lt.1.D-40) R1=0.d0
       if(abs(R2).lt.1.D-40) R2=0.d0
       write(13,*) j1,k1,R1,R2
       write(12,*) j,k,R1,R2
       write(14,*) j,k,R3
       enddo
       write(13,*)
       write(12,*)
       write(14,*)
       enddo
       close(12)
       close(13)
       close(14)

100    format(E14.5)

       stop

       contains 

       subroutine input_vr()
       implicit double precision (a-h,o-z)
*************************************************************
*************************************************************

       open(11,file="scatte_st_3D.out",form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) numE
  
 
	write(6,*) 'Total number of energy states',numE
	write(6,*) 'Input index of energy to be plotted'
        read(5,*) iE_in 


       nr=n1*n2*n3
       nr_n=nr/nnodes
       allocate(vr_tmp(nr_n))
       allocate(vr(n1,n2,n3))

ccccccccc  roll over to the starting position in the file 11
	do  iE=1,iE_in-1
	read(11) jE,num_st(iE)
        do ist=1,num_st(iE)
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

        enddo !end loop over num_st(iE)
	enddo !end loop over (1,iE_in-1)


	read(11) jE,num_st_in
	write(6,*) 'Total number of scattering states',num_st_in
	write(6,*) 'Input index of scattering state to be plotted'
        read(5,*) ist_in 

	do ist=1,ist_in
	write(6,*) 'read in for ist',ist
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

	enddo !last vr(i,j,k) corresponds to ist_in

       close(11)

	deallocate(vr_tmp)
********************************************************
********************************************************
       return
       end subroutine input_vr

       end

      

