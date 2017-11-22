       program gen_vext_den

ccccc This program plot all types of charge densities. 

       implicit double precision (a-h,o-z)

       real*8 AL(3,3),ALt(3,3),ALI(3,3)
       real*8 x(3,1000),zatom(10),zzatom(1000)
       integer iat(1000),itype(10)

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
       do i=1,nat
       read(10,*) iat(i),x(1,i),x(2,i),x(3,i)
       enddo
       close(10)

       call get_ALI(AL,ALI)

       write(6,*) "--------- AL ----------"
       write(6,100) AL(1,1),AL(2,1),AL(3,1)
       write(6,100) AL(1,2),AL(2,2),AL(3,2)
       write(6,100) AL(1,3),AL(2,3),AL(3,3)
       write(6,*) "--------- AL ----------"
100    format(3(f12.6,1x))

       ntype=0
       do i=1,nat
       do j=1,ntype
       if(itype(j).eq.iat(i))  goto 20
       enddo
       ntype=ntype+1
       itype(ntype)=iat(i)
20     continue
       enddo
ccccccccccccccc
       write(6,*) "ntype=",ntype
       do j=1,ntype
       write(6,*) "input val. el charge (positive) for atom type=", 
     &  itype(j)
       read(5,*) zatom(j)
       enddo
cccccccccccccccccccccccccccc
       totNel=0.d0
       do i=1,nat
       do j=1,ntype
       if(itype(j).eq.iat(i)) zzatom(i)=zatom(j)
       enddo
       totNel=totNel+zzatom(i)
       enddo
       write(6,*) "totNel=", totNel
cccccccccccccccccccccccccccc
      
       write(6,*) "input the name of dens file"
       read(5,*) filename
       open(11,file=filename,form="unformatted")
       rewind(11)
       read(11) n1,n2,n3,nnodes
       read(11) ALt
 
       dd=0.d0
       do j=1,3
       do i=1,3
       dd=dd+dabs(AL(i,j)-ALt(i,j))
       enddo
       enddo
       if(dd.gt.1.D-3) then
       write(6,*) "AL from two inputs are not the same, stop"
       write(6,*) ALt
       stop
       endif

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
cccccccccccccccccccccccccccccccccccccccccccccccc
       write(6,*) "input x1"
       read(5,*) x1
       sum=0.d0
       sum2=0.d0
       do k=1,n3
       do j=1,n2
       do i=1,n1
       xt=(i-1.d0)/n1
       if(xt.gt.x1) xt=xt-1.d0
       sum=sum+vr(i,j,k)*xt
       sum2=sum2+vr(i,j,k)
       enddo
       enddo
       enddo
       sum=sum/sum2*totNel
       do i=1,nat
       xt=x(1,i)
       if(xt.gt.x1) xt=xt-1.d0
       sum=sum-zzatom(i)*xt
       enddo

       vol=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))
     &     +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))
     &     +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))
       vol=dabs(vol)
       h=AL(1,1)*ALI(1,1)+AL(2,1)*ALI(2,1)+AL(3,1)*ALI(3,1)
       h=h/dsqrt(ALI(1,1)**2+ALI(2,1)**2+ALI(3,1)**2)
       sum=sum*h**2/vol      ! general formula
       pi=4*datan(1.d0)

       v_jump=4*pi*sum

       write(6,*) "v_jump(a.u)=",v_jump
       
       open(11,file="vext_dens",form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL

       sum=0.d0
       do iread=1,nnodes

       do ii=1,nr_n

       jj=ii+(iread-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       xt=(i-1.d0)/n1
       if(xt.lt.x1) then
       vr_tmp(ii)=v_jump*(xt-x1+0.5d0)
       else
       vr_tmp(ii)=v_jump*(xt-x1-0.5d0)
       endif
       sum=sum+vr_tmp(ii)*vr(i,j,k)
       enddo

       write(11) vr_tmp

       enddo

       close(11)
       write(6,*) "potential written in vext_dens"
ccccccccccccccccccccccccccccccccccccccc
       sum=sum*vol/(n1*n2*n3)
       sum1=0.d0
       do i=1,nat
       xt=x(1,i)
       if(xt.lt.x1) then
       vt=v_jump*(xt-x1+0.5d0)
       else
       vt=v_jump*(xt-x1-0.5d0)
       endif
       sum1=sum1-zzatom(i)*vt
       enddo

       write(6,200) sum,sum1,sum1+sum
200    format("E_rhoV,E_IV,E_tot=",3(E14.7,1x))

       stop
       end



      subroutine get_ALI(AL,ALI)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


      real*8 AL(3,3),ALI(3,3)
      real*8 tmp(3)

      do i=1,3
         do j=1,3
            ALI(j,i)=AL(i,j)
         enddo
         tmp(i)=1
      enddo

      call gaussj(ALI,3,3,tmp,1,1)

*****************************************
*     *  \sum_i AL(i,j1)*ALI(i,j2)= \delta_j1,j2
*     *  2*pi*ALI(i,j) is the jth G vector unit in (i)x,y,z components
*****************************************
      return
      end

       SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
cccccccccc this is double precision gauss program
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=200)
      real*8 A(NP,NP),B(NP,MP)
      integer IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG= 0.d0
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix 1'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW

        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.d0) then
        PAUSE 'Singular matrix. 2'
        endif
        PIVINV=1.d0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.d0
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.d0
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

      
      
      
      
      
  
     

