       program plot_all

ccccc This program plot all types of charge densities. 

       implicit double precision (a-h,o-z)

       real*8 AL(3,3),ALI(3,3)
       character*50 file1

       character*50 filename

       real*8, allocatable, dimension (:,:,:) :: vr
       real*8, allocatable, dimension (:,:,:) :: vr_xyz
       real*8, allocatable, dimension(:) :: vr_tmp

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
       write(6,*) "input the full name of group control file"
       read(5,*) filename
       write(6,*) "input the ifrag and num_proc_group"
       read(5,*) ifrag,num_proc_group

       open(10,file=filename)
       read(10,*) num_frag,num_tot_rec,LREC 
       do i=1,ifrag
       read(10,*) natom,n1,n2,n3,totNel,ic1,ic2,ic3,ns1,
     &  ns2,ns3,ind_frag,ipos_rec,num_rec
       enddo
       close(10)
       write(6,*) "natom,n1,n2,n3,totNel,ic1,ic2,ic3,ns1,ns2,ns3,"//
     &   "ind_frag,ipos_rec,num_rec"


       write(6,*) natom,n1,n2,n3,totNel,ic1,ic2,ic3,ns1,ns2,ns3,
     &  ind_frag,ipos_rec,num_rec
        
       
       LREC8=LREC*8
       write(6,*) "input the full name of group dens file"
       read(5,*) filename


       open(11,file=filename,form="unformatted",access="direct",
     &   recl=LREC8,iostat=ierr)
       rewind(11)

       nnodes=num_proc_group
       nr=n1*n2*n3
       nr_n=nr/nnodes
       num_rec_proc=num_rec/num_proc_group
       LREC_end=nr_n-(num_rec_proc-1)*LREC
       write(6,*) "ipos_rec,num_rec_proc,num_proc_group,LREC,LREC_end",
     & ipos_rec,num_rec_proc,num_proc_group,LREC,LREC_end

       allocate(vr_tmp(LREC))
       allocate(vr(n1,n2,n3))

       do iread=1,num_proc_group
       do irec=1,num_rec_proc          ! this do loop to replace the original read(*) vr_tmp2
       length=LREC
       if(irec.eq.num_rec_proc) length=LREC_end
       ipos=ipos_rec+(iread-1)*num_rec_proc+irec-1
       read(11,rec=ipos) vr_tmp(1:length)

       do ii=1,length
       jj=(iread-1)*nr_n+(irec-1)*LREC+ii
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       vr(i,j,k)=vr_tmp(ii)
       enddo
       enddo
       enddo

       close(11)

********************************************************
************************************************************
       return
       end subroutine input_vr

       end

      

