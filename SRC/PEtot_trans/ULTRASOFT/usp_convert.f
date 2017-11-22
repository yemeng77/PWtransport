c       subroutine w_line_usp(vwr_atom,zatom,iiatom,
c     &  Ealpha,occ_t,isNL,Dij0,Qij,qijrad,
c     &  icore,is_ref,
c     &  ip_ref,id_ref,qi,wq,qi2,vq,
c     &  rhoq,rhocq,vqT,ri,amr,lll,nbeta,rcut_q1t,
c     &  rcut_q2t,qfuncLM0,r_at,a_r,b_r)

         program usp_convert
c----------------------------------------------------------------------------
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

      implicit double precision(a-h,o-z)
c      include "param.escan_real"

      parameter(mnq=2000)


      real*8 rhoc(2000),vlocT(2000),vloc(2000)
      real*8 ch_at(2000),vw(2000),qfuncLM0(12,2000)
      real*8 r_at(2000)

 
      real*8 qi2(mnq),vq(mnq),rhoq(mnq)
      real*8 rhocq(mnq)
      real*8 vqT(mnq)

cccccccc wq has 1/rmask for real space implementation
cccccccc    and has no 1/rmask for q space implementation

      real*8 qi(mnq),wq(mnq,8)
      real*8 ri(201),amr(201)
      real*8 Dij0(32,32),Qij(32,32),qijrad(0:6,32,32)

      integer  iiatom

      real*8  zatom,Ealpha
      real*8  occ_s,occ_p,occ_d,occ_t
      integer isNL(3),icore
      integer is_ref,ip_ref,id_ref
      integer is_TB,ip_TB,id_TB


      character*50  vwr_atom,vwr_atom2
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
c

cccc idim3 must be 8, otherwise, many of the other dimensions in other subroutines
cccc need to be changed !
      parameter( idim1 = 1000      , idim2 = 26         )
      parameter( idim3 = 8        , idim4 = 5          )
      parameter( idim5 = 4         , idim6 = 2*idim5-1  )
      parameter( idim7 = 4         , idim8 = 20         )

c     idim1  .ge.  no. of points in the radial mesh
c     idim2  .ge.  no. of shells in the atom
c     idim3  .ge.  no. of beta functions in vanderbilt potential
c     idim4  .ge.  no. of valence states for pseudization
c     idim5  .ge.  lmax + 1, lmax is maximum l value of potential
c     idim7  .ge.  no. of reference states for each angular momentum value
c     idim8  .ge.  no. of terms in taylor expansion of pseudized q_i,j



c
c.....radial mesh information
      dimension r(idim1),rab(idim1),rinner(idim6)
c.....charge densities
      dimension rsatom(idim1),rspsco(idim1)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vlocSR(idim1),vloc0(idim1)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3)
c.....version number and date
      dimension iver(3),idmy(3)
c.....logic for logarithmic mesh types
      logical tlog
c.....wave functions
      dimension snl(idim1,idim2)
c
c.....title
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
c
      integer stderr

      real*8 dum(idim1),qrad(0:6,idim3,idim3)
c      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
c--------------------------------------------------------------------------
c
c       r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c
        write(6,*) "input the usp name"
        read(5,*) vwr_atom
        write(6,*) "input the txt output usp name"
        read(5,*) vwr_atom2
       

        iops=10
        iops2=11
        iout=6
        stderr=6
      

c       save data so compatibility checks can be run later
c        mold  = mesh
c        r2    = r(2)
c        rmesh = r(mesh)
c
        open( unit = iops , file = vwr_atom , status = 'old' ,
     +    form = 'unformatted' )
        open( unit = iops2 , file = vwr_atom2 )
c
        read (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
        write(iops2,*) (iver(i),i=1,3),(idmy(i),i=1,3)
c
        if ( iver(1) .gt. 99 .or. iver(1) .lt. 1 .or. 
     +       iver(2) .gt. 9  .or. iver(2) .lt. 0 .or.
     +       iver(3) .gt. 9  .or. iver(3) .lt. 0       ) then
c
c         handle errors caused by pre-version number files
c
          write(stderr,*) '***error reading pseudopotential file'
          write(stderr,*) 'file does not contain version number'
          write(stderr,*) 'taking corrective action - will assume:'
c
          iver(1) = 2
          iver(2) = 1
          iver(3) = 0
c
          idmy(1) = 0
          idmy(2) = 0
          idmy(3) = 0
c
          write(stderr,90) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
   90     format('pseudopotential program version',i3,'.',
     +      i1,'.',i1,'   date:',i3,' - ',i2,' - ',i4)
c
          rewind (iops)

          stop     ! I hope this will not occure
c
        endif
c
c       set the mesh type
        if ( iver(1) .eq. 1 ) then
          tlog = .false.
        else
          tlog = .true.
        endif
c
c
        read (iops) title,zps,zvps,exftps,nvalps,mesh,etot
        write (iops2,*) title,zps,zvps,exftps,nvalps,mesh,etot
        read (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        write (iops2,*) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        read (iops) keyps,ifpcor,rinner(1)
        write (iops2,*) keyps,ifpcor,rinner(1)
c
        if ( keyps .ne. 3 ) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) 'keyps =',keyps,' out of programed range'
          call exit(1)
        endif
c
        if ( iver(1) .ge. 3 ) then
          read(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
          write(iops2,*) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
c
        if (10*iver(1)+iver(2).ge.51) then
          read (iops) (rinner(i),i=1,nang*2-1)
          write (iops2,*) (rinner(i),i=1,nang*2-1)
        else
          if (nang.gt.1) then
            do i=2,2*nang-1
              rinner(i)=rinner(1)
            end do
          endif
        endif
c
        if ( iver(1) .ge. 4 ) then
          read(iops) irelps
          write(iops2,*) irelps
        endif
c
c
c       set the number of angular momentum terms in q_ij to read in
        if ( iver(1) .eq. 1 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 5
        elseif ( iver(1) .eq. 2 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 2 * nvalps - 1
        else
          nlc = 2 * nang - 1
        endif
c
c
        read (iops) (rc(i),i=1,nang)
        write (iops2,*) (rc(i),i=1,nang)
        read (iops) nbeta,kkbeta
        write (iops2,*) nbeta,kkbeta
        do 100 j=1,nbeta
          read (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          write (iops2,*) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 100 k=j,nbeta
            read (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
            write (iops2,*) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
            ddd0(k,j)=ddd0(j,k)
            ddd (k,j)=ddd (j,k)
            qqq (k,j)=qqq (j,k)
            do 100 i=1,kkbeta
              qfunc(i,k,j)=qfunc(i,j,k)
  100   continue
        do 110 i=1,nqf
        do 110 lp = 1,nlc
          qfcoef(i,lp,k,j)=qfcoef(i,lp,j,k)
  110   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          read (iops) (iptype(j),j=1,nbeta),npf,ptryc
          write (iops2,*) (iptype(j),j=1,nbeta),npf,ptryc
        endif
c
        read (iops) rcloc,(vloc0(i),i=1,mesh)
        write (iops2,*) rcloc,(vloc0(i),i=1,mesh)
c
c       set index arrays nbl and nbl0
        do 150 lp = 1,nang
          nbl(lp) = 0
  150   continue
        do 160 ib = 1,nbeta
          lp = lll(ib) + 1
          nbl(lp) = nbl(lp) + 1
  160   continue
        lmaxps = lll(nbeta)
        nbl0(1) = 0
        do 170 i = 1,lmaxps
          nbl0(i+1) = nbl0(i) + nbl(i)
  170   continue
c
c       possible readin of the frozen core correction
        if (ifpcor.gt.0) then
           read(iops) rpcor
           write(iops2,*) rpcor
           read (iops) (rspsco(i),i=1,mesh)
           write (iops2,*) (rspsco(i),i=1,mesh)
        endif
c
        read (iops) (vlocSR(i),i=1,mesh)
        write (iops2,*) (vlocSR(i),i=1,mesh)
        read (iops) (rsatom(i),i=1,mesh)
        write (iops2,*) (rsatom(i),i=1,mesh)
c
c       with the logarithmic mesh potentials grid information included
        if ( tlog ) then
          read (iops) (r(i),i=1,mesh)
          write (iops2,*) (r(i),i=1,mesh)
          read (iops) (rab(i),i=1,mesh)
          write (iops2,*) (rab(i),i=1,mesh)
        endif
        if (iver(1) .ge. 6) then
c           nchi = nvales        ! warning, used before defined ! use nvalps ?
           nchi = nvalps        ! warning, used before defined ! use nvalps ?
           if (iver(1) .ge. 7) then
           read (iops) nchi
           write (iops2,*) nchi
           endif
           ncores=0
           read (iops) ((snl(i,j), i=1,mesh),j=ncores+1,ncores+nchi)
           write (iops2,*) ((snl(i,j), i=1,mesh),j=ncores+1,ncores+nchi)
        endif
c
        close (iops)
        close (iops2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   END OF THE READ IN
cccccccccccccccccccccccccccccccccccccccccccccc
          stop

          end

          
