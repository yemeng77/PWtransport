       program plot_DOS
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

****************************************
****  It stores the wavefunction in G space, only in half
****  of the E_cut sphere. 
******************************************

      implicit double precision (a-h,o-z)

       parameter (matom=4000,mtype=6,nE=500)
******************************************
       real*8 AL(3,3),AL_t(3,3)
***********************************************
*************************************************
       real*8 xatom(3,matom),w(matom),weight(3,mtype)
       integer iatom(matom),imov_at(3,matom),iiatom(mtype),
     &    ityatom(matom)
       integer numref(matom),iref_start(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)

ccccccccc not really used  
       integer smatr(3,3,48),nrot
       integer iCGmth0(100),iCGmth1(100),iscfmth0(100),iscfmth1(100)
       real*8 FermidE0(100),FermidE1(100),totNel
       integer itypeFermi0(100),itypeFermi1(100)
       integer kpt_dens(2),ispin_dens(2),iw_dens(2)
       integer niter0,nline0,niter1,nline1
       integer icoul
       real*8 xcoul(3)
       integer ipsp_type(mtype),nref_type(mtype)
       character*20 vwr_atom(mtype),fforce_out,fdens_out
       character*20 fwg_in(2),fwg_out(2),frho_in(2),frho_out(2),
     & fvr_in(2),fvr_out(2),f_tmp,fxatom_out,fvext_in
       character*20 f_xatom,sym_file,kpt_file
cccccccccccccccccccccccccccccc
       real*8, allocatable, dimension(:):: weighkpt,akx,aky,akz 
       real*8, allocatable, dimension(:,:,:):: E_st
       complex*16, allocatable, dimension(:,:) :: beta_psi_tmp
       complex*16, allocatable, dimension(:,:,:,:) :: beta_psi
       real*8 DOS(nE),DOS_ref(nE,3,mtype)
       character*9 fname
 
**************************************************
**************************************************

       open(9,file='etot.input',status='old',action='read',iostat=ierr) 
       if(ierr.ne.0) then
       write(6,*) "file ***etot.input*** does not exist, stop"
       stop
       endif
       read(9,*)
       read(9,*,iostat=ierr) i1, f_xatom
       if(ierr.ne.0) call error_stop(i1)
       call readf_xatom()
       read(9,*,iostat=ierr) i1, n1,n2,n3,n1L,n2L,n3L
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, islda,igga
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, Ecut,Ecut2,Ecut2L,Smth
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, icoul,xcoul(1),xcoul(2),xcoul(3)
       if(ierr.ne.0) call error_stop(i1)
       if(islda.ne.1.and.islda.ne.2) then
       write(6,*) "islda must be 1 (lda) or 2 (slda), stop", islda
       stop
       endif
       if(igga.ne.0.and.igga.ne.1) then
       write(6,*) "igga must be 0 (no gga) or 1 (gga), stop", igga
       stop
       endif
       read(9,*,iostat=ierr) i1, iwg_in,(fwg_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, iwg_out,(fwg_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, irho_in,(frho_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, irho_out,(frho_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, ivr_in,(fvr_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, ivr_out,(fvr_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, ivext_in,fvext_in
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1,idens_out,kpt_dens(1),kpt_dens(2),
     &  ispin_dens(1),ispin_dens(2),iw_dens(1),iw_dens(2),
     &  fdens_out
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, iforce,fforce_out
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, isym,sym_file
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, ikpt_yno,kpt_file
       if(ierr.ne.0) call error_stop(i1)
         
       read(9,*,iostat=ierr) i1, totNel,mx,tolug,tolE
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, niter0,nline0
       if(ierr.ne.0) call error_stop(i1)
         do i=1,niter0
         read(9,*,iostat=ierr) iCGmth0(i),iscfmth0(i),FermidE0(i),
     &   itypeFermi0(i)
       if(ierr.ne.0) call error_stop(i1)
         FermidE0(i)=FermidE0(i)/27.211396d0
	 enddo
       read(9,*,iostat=ierr) i1, num_mov,tolforce,dtstart,
     &    dd_limit,fxatom_out,imv_cont
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, niter1,nline1
       if(ierr.ne.0) call error_stop(i1)
         do i=1,niter1
         read(9,*,iostat=ierr) iCGmth1(i),iscfmth1(i),FermidE1(i),
     &    itypeFermi1(i)
       if(ierr.ne.0) call error_stop(i1)
         FermidE1(i)=FermidE1(i)/27.211396d0
	 enddo
       read(9,*,iostat=ierr) i1, ilocal
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, rcut
       if(ierr.ne.0) call error_stop(i1)
       read(9,*,iostat=ierr) i1, ntype
       if(ierr.ne.0) call error_stop(i1)
       do ia=1,ntype
       read(9,*,iostat=ierr)    vwr_atom(ia),ipsp_type(ia)
       if(ierr.ne.0) call error_stop(i1)
       enddo
       close(9)

       ipsp_all=1
       do ia=1,ntype
       if(ipsp_type(ia).eq.2) ipsp_all=2
       enddo
       if(ilocal.eq.1) ipsp_all=1
       
*************************************************
**** change Ecut,Eref to A.U
*************************************************
       Ecut=Ecut/2
       Ecut2=Ecut2/2
       Ecut2L=Ecut2L/2
*************************************************
       nrot=1
*************************************************

       nr=n1*n2*n3

*************************************************
       nh1=n1/2+1
*************************************************
        do ia=1,ntype
        if(ipsp_type(ia).eq.1) then
        call readvwr_head()
        else
        call readusp_head(vwr_atom(ia),iiatom(ia),nref_type(ia))
        endif
        enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            nref_tot=0
            mref=0
           do ia=1,natom
            iref_start(ia)=nref_tot
            iitype=0
            do itype=1,ntype
            if(iatom(ia).eq.iiatom(itype)) iitype=itype
            enddo
            if(iitype.eq.0) then
            write(6,*) "itype not found, stop", iatom(ia),ia
            stop
            endif
            ityatom(ia)=iitype
            numref(ia)=nref_type(iitype)
            nref_tot=nref_tot+nref_type(iitype)
            if(nref_type(iitype).gt.mref) mref=nref_type(iitype)
        enddo

       nref_tot_tmp=nref_tot
       natom_tmp=natom
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       open(23,file="eigen_all.store",form="unformatted")
       rewind(23)
       read(23) islda,nkpt,mx,nref_tot,natom,nnodes

cccc Note, usually, eigen_all.store is from previous step calculation, its nref_tot is smaller than in the bpsiifilexxxx

       nref_tot=nref_tot_tmp       

       if(nref_tot_tmp.ne.nref_tot.or.
     &  natom_tmp.ne.natom) then
       write(6,*) "nref_tot,natom,changed, stop", 
     &  nref_tot_tmp,nref_tot,natom_tmp,natom
       stop
       endif
       write(6,*) "islda,nkpt,mx=",islda,nkpt,mx
       allocate(weighkpt(nkpt))
       allocate(akx(nkpt))
       allocate(aky(nkpt))
       allocate(akz(nkpt))
       allocate(E_st(mx,nkpt,islda))
       do iislda=1,islda
       do kpt=1,nkpt
       read(23) iislda_tmp,kpt_tmp,weighkpt(kpt),
     &       akx(kpt),aky(kpt),akz(kpt)
       read(23) (E_st(i,kpt,iislda),i=1,mx)
       enddo
       enddo
       close(23)

       
       mx_n=mx/nnodes+1
       allocate(beta_psi_tmp(nref_tot,mx_n))
       allocate(beta_psi(nref_tot,mx,nkpt,islda))

       do iislda=1,islda
       do kpt=1,nkpt

       fname="bpsiiofil"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

        do inode=1,nnodes
        read(10) beta_psi_tmp
        do im=1,mx_n
        m=im+(inode-1)*mx_n
        if(m.le.mx) then
        beta_psi(:,m,kpt,iislda)=beta_psi_tmp(:,im)


        endif
        enddo
        enddo
        close(10)
        enddo
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Emin=100000
       Emax=-100000
       do iislda=1,islda
       do kpt=1,nkpt
       do m=1,mx
       if(E_st(m,kpt,iislda).lt.Emin) Emin=E_st(m,kpt,iislda)
       if(E_st(m,kpt,iislda).gt.Emax) Emax=E_st(m,kpt,iislda)
       enddo
       enddo
       enddo
       write(6,*) "Emin,Emax=",Emin,Emax

       write(6,*) "input E_b (eV), Gaussian broadening"
       read(5,*) E_b

       Emin=Emin-4*E_b
       Emax=Emax+4*E_b


       dE=(Emax-Emin)/nE
       nb=2*E_b/dE

       DOS=0.d0
       DOS_ref=0.d0

       do iislda=1,islda
       do kpt=1,nkpt
       do m=1,mx
       
cccccccccccccccccccccccccccccccc
       weight=0.d0

       iref=0
       do ia=1,natom
        itype=ityatom(ia) 

       do ll=1,is_ref(itype)
       iref=iref+1
       weight(1,itype)=weight(1,itype)+
     &     abs(beta_psi(iref,m,kpt,iislda))**2*w(ia)

c        write(6,*) "beta",m,kpt,iref,iref,beta_psi(iref,m,kpt,iislda)
       enddo

       do ll=1,ip_ref(itype)*3
       iref=iref+1
       weight(2,itype)=weight(2,itype)+
     &     abs(beta_psi(iref,m,kpt,iislda))**2*w(ia)

c        write(6,*) "beta",m,kpt,iref,beta_psi(iref,m,kpt,iislda)
       enddo

       do ll=1,id_ref(itype)*5
       iref=iref+1
       weight(3,itype)=weight(3,itype)+
     &     abs(beta_psi(iref,m,kpt,iislda))**2*w(ia)

c        write(6,*) "beta",m,kpt,iref,iref,beta_psi(iref,m,kpt,iislda)
       enddo

       enddo    ! ia

ccccccccccccccccccccccccccccccccccccccc
       Ep=E_st(m,kpt,iislda)
       iEp=(Ep-Emin)/dE+1

       do i=-nb,nb
       fact=exp(-(i*2.d0/nb)**2)
       ip=iEp+i
       DOS(ip)=DOS(ip)+fact*weighkpt(kpt)

       do itype=1,ntype
       DOS_ref(ip,1,itype)=DOS_ref(ip,1,itype)+
     &  fact*weighkpt(kpt)*weight(1,itype)
       DOS_ref(ip,2,itype)=DOS_ref(ip,2,itype)+
     &  fact*weighkpt(kpt)*weight(2,itype)
       DOS_ref(ip,3,itype)=DOS_ref(ip,3,itype)+
     &  fact*weighkpt(kpt)*weight(3,itype)
       enddo

       enddo   ! i
cccccccccccccccccccccccccccccccccccccccccccc

       enddo
       enddo
       enddo
cccccccccccccccccccccccccccccccccccccccccccc
ccccc s,p,d, need to rescale by the (psi*dV(r))^2 factor
       

       open(10,file="graph.DOS")
       rewind(10)
       do i=1,nE
       write(10,200) Emin+(i-1)*dE, DOS(i),
     &   (DOS_ref(i,1,itype),DOS_ref(i,2,itype),
     &   DOS_ref(i,3,itype),itype=1,ntype)
       enddo
       close(10)

200    format(E12.5,2x,E10.4,2x,15(E10.4,1x))



*************************************************
      stop
      contains

*******************************************
      subroutine readvwr_head()

      implicit double precision (a-h,o-z)


      open(10,file=vwr_atom(ia),status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
       write(6,*) "vwr_file ***",filename,"*** does not exist, stop"
       stop
      endif
      read(10,*) nrr_t,ic_t,iiatom(ia),zatom_t,iloc_t,occ_s_t,
     &  occ_p_t,occ_d_t
      read(10,*) is_ref_t,ip_ref_t,id_ref_t,
     &  is_TB_t,ip_TB_t,id_TB_t
      close(10)


      if(iloc_t.eq.1) is_ref_t=0
      if(iloc_t.eq.2) ip_ref_t=0
      if(iloc_t.eq.3) id_ref_t=0
ccccccccccccccccccccccccccccccccccccccccccccccc
      is_ref_t=1
      ip_ref_t=1
      id_ref_t=1
cccccccccccccccccccccccccccccccccccccccccccccccccc



      nref_type(ia)=is_ref_t+ip_ref_t*3+id_ref_t*5
      is_ref(ia)=is_ref_t
      ip_ref(ia)=ip_ref_t
      id_ref(ia)=id_ref_t

     

      return
      end subroutine readvwr_head



      subroutine readf_xatom()
      implicit double precision (a-h,o-z)
      real*8, allocatable, dimension(:,:) :: xatom_tmp
      integer, allocatable, dimension(:,:) :: imov_at_tmp
      integer, allocatable, dimension(:)   ::  iatom_tmp
      real*8, allocatable, dimension(:)   ::  w_tmp
   

       write(6,*) "input: 0 (total DOS) or 1 (partial DOS)"
       read(5,*) ipart_DOS
       if(ipart_DOS.eq.1) then
       write(6,*)  "For atom based partial DOS, there must be an"//
     &  " extra last column (atomic weight) in ", f_xatom 
       endif
   
       open(10,file=f_xatom,status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
        write(6,*) "f_xatom file **",f_xatom, " ***does not exist, stop"
        stop
       endif
       
       rewind(10)
       read(10,*) natom
      
       if(natom.gt.matom) then
       write(6,*) "natom.gt.matom, increase matom in data.f, stop"
       stop
       endif

       allocate(xatom_tmp(3,natom))
       allocate(imov_at_tmp(3,natom))
       allocate(iatom_tmp(natom))
       allocate(w_tmp(natom))

       read(10,*) (AL(i,1),i=1,3)
       read(10,*) (AL(i,2),i=1,3)
       read(10,*) (AL(i,3),i=1,3)
       if(ipart_DOS.eq.0) then
       do i=1,natom
       read(10,*) iatom_tmp(i),xatom_tmp(1,i),
     & xatom_tmp(2,i),xatom_tmp(3,i),
     & imov_at_tmp(1,i),imov_at_tmp(2,i),imov_at_tmp(3,i)
       w_tmp(i)=1.d0
       enddo
       else
       do i=1,natom
       read(10,*) iatom_tmp(i),xatom_tmp(1,i),
     & xatom_tmp(2,i),xatom_tmp(3,i),
     & imov_at_tmp(1,i),imov_at_tmp(2,i),imov_at_tmp(3,i),w_tmp(i)
       enddo
       endif

       close(10)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Now, re-arrange xatom, so the same atoms are consequentive together. 
cccc This is useful to speed up the getwmask.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       ii=0
100    continue
       ncount=0
       itype_tmp=-2
       do i=1,natom
       if(itype_tmp.eq.-2.and.iatom_tmp(i).ne.-1)
     &  itype_tmp=iatom_tmp(i)

       if(iatom_tmp(i).eq.itype_tmp) then
       ii=ii+1
       ncount=ncount+1
       iatom(ii)=iatom_tmp(i)
       iatom_tmp(i)=-1
       xatom(:,ii)=xatom_tmp(:,i)
       imov_at(:,ii)=imov_at_tmp(:,i)
       w(ii)=w_tmp(i)
       endif
       enddo
       if(ncount.gt.0) goto 100
       if(ii.ne.natom) then
       write(6,*) "something wrong to rearrange xatom, stop"
       stop
       endif

       deallocate(xatom_tmp)
       deallocate(imov_at_tmp)
       deallocate(iatom_tmp)
       deallocate(w_tmp)
       
      return
      end subroutine readf_xatom
**************************************************

      subroutine error_stop(i1)
      implicit double precision(a-h,o-z)
      write(6,*) "error in etot.input, line=",i1
      stop

      return
      end subroutine error_stop

      end
      
      

