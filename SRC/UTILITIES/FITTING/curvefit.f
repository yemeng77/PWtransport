       program curvefit
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc This program reads in the transmission coefficients from tmp [T(k(E_i))]
ccc and interpolates them to the whole k=[0,1] region (\Gamma-X region).
ccc It then calculates the conductance in the energy region of interest. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       parameter (nm=101)
       implicit double precision (a-h,o-z)
       real*8 ak(200),w(200),a(200),Ew(200)
       real*8 aa(0:200,0:200),y(0:200),yy(0:200,200)
       integer ipiv(0:200)
       real*8 akall(1500),Tall(1500),Eall(1500),
     &  Ref(1500),Asum(1500),afL(1500),afR(1500) 
       real*8 afL2(1500),afR2(1500)
       integer iband(1500),ind_b(80),ibflag(200)
       real*8 ak_max(200),ak_min(200)
       character*20 filename

       integer ist_linew(nm,nm),ikpt_linew(nm,nm),numw(nm)
       real*8 E_linew(nm,nm)
       real*8 Et(1001),Tsum(1001)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(10,file="curvefit.input")
	rewind(10)
	read(10,*) E_Fermi
	read(10,*) dV
	read(10,*) alpha
	read(10,*) E1 
	read(10,*) E2
	close(10)	
	write(6,*) "Fermi energy(eV): ",E_Fermi
	write(6,*) "Bias voltage(eV): ",dV
	write(6,*) "Energy min and max:", E2,E1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(14,file="E_line_W2K2")
      read(14,*) nstw,nkptw
      do j=1,nstw
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      close(14)
cccccccccccccccccccccccccccccccccc

       open(10,file="tmp")
       rewind(10)
       read(10,*) mm
       do ll=1,mm
       read(10,406) ia,iband(ll),akall(ll),Eall(ll),
     & Tall(ll),Ref(ll),Asum(ll),afL(ll),afR(ll),
     & afL2(ll),afR2(ll),
     &	qq1,qq2,qq3,qq4,n1,n2

       Tall(ll)=Tall(ll)/Asum(ll)
       enddo
       close(10)
405   format(2(i4,1x),f10.6,2x,f12.6,2x,3(E13.6,1x),3x,
     &  2(E13.6,1x),3x,4(E8.2,1x))
406   format(2(i4,1x),f10.6,1x,f12.6,2x,3(E17.10,1x),2x,
     &  2(E9.2,1x),2x,2(E9.2,1x),3x,4(E8.2,1x),1x,2(i3,1x))

ccccccccccccccccccccccccccccccccccccccccccccccc

       ind_b=0 
       yy=0.d0
       ibflag=0

1000   continue

       ib=0
       do ll=1,mm
       if(ind_b(iband(ll)).eq.0) then
       ib=iband(ll)
       ind_b(ib)=1
       goto 300
       endif
       enddo

300    continue

       if(ib.eq.0) goto 2000

      open(25,file="tmp_filtered")

       nn=0
       ak_max(ib)=0.d0
       ak_min(ib)=1.d0
       do ll=1,mm
       if(iband(ll).eq.ib) then


c       if((abs(Asum(ll)-1.d0).lt.0.05.and.
c     &   afR(ll).gt.0.95.and.
c     &  (afL(ll).gt.0.85.or.(Tall(ll).lt.1.E-5.and.afL(ll).gt.0.1))
c     &  .and.Tall(ll).gt.1.E-10).or.ib.eq.11) then    ! special  

	if(abs(Asum(ll)-1.d0).lt.0.005) then

       write(25,1006) iband(ll),akall(ll),Eall(ll),
     & Tall(ll),Ref(ll),Asum(ll)

1006   format(1I5,5(E13.6,1x))

        nn=nn+1
        ak(nn)=akall(ll)
        w(nn)=dlog(abs(Tall(ll)))
        Ew(nn)=Eall(ll)
        if(ak(nn).gt.ak_max(ib)) ak_max(ib)=ak(nn)
        if(ak(nn).lt.ak_min(ib)) ak_min(ib)=ak(nn)
        endif
       endif
       enddo
       sum=sum/nn


       if(nn.eq.0) goto 1000

c       alpha=1.d0


       aa=0.d0
       y=0.d0

       do 100 ll=1,nn

       ww=1.d0*dexp(w(ll)/5)

c       akx=ak(ll)*99.999999
       akx=ak(ll)*199.999999
       i1=akx
       x1=i1+1-akx
       x2=1-x1
       i2=i1+1

       aa(i1,i1)=aa(i1,i1)+2*x1**2*ww
       aa(i1,i2)=aa(i1,i2)+2*x1*x2*ww
       y(i1)=y(i1)+2*w(ll)*x1*ww

       aa(i2,i2)=aa(i2,i2)+2*x2**2*ww
       aa(i2,i1)=aa(i2,i1)+2*x2*x1*ww
       y(i2)=y(i2)+2*w(ll)*x2*ww

100    continue

cccccccccccccccccccccccccccccccccccc
c       do 101 i=1,99
       do 101 i=1,199

       aa(i-1,i-1)=aa(i-1,i-1)+2*alpha
       aa(i-1,i+1)=aa(i-1,i+1)+2*alpha
       aa(i-1,i)=aa(i-1,i)-4*alpha
      
       aa(i+1,i-1)=aa(i+1,i-1)+2*alpha
       aa(i+1,i+1)=aa(i+1,i+1)+2*alpha
       aa(i+1,i)=aa(i+1,i)-4*alpha

       aa(i,i-1)=aa(i,i-1)-4*alpha
       aa(i,i+1)=aa(i,i+1)-4*alpha
       aa(i,i)=aa(i,i)+8*alpha
101    continue
cccccccccccccccccccccccccccccccccccccc

c       call dgesv(101,1,aa,101,ipiv,y,101,info)
       call dgesv(201,1,aa,201,ipiv,y,201,info)

       if(ib.lt.10) then
       filename="graph."//char(48+ib)
       else
       ibt=ib/10
       filename="graph."//char(48+ibt)//char(48+ib-ibt*10)
       endif
   
       open(10,file=filename)
       rewind(10)
       do i=0,200
       write(10,*) i/200.d0,y(i),dexp(y(i))
       enddo
106   format(3(E13.6,1x))
       close(10)

       if(ib.lt.10) then
       filename="graphT."//char(48+ib)
       else
       ibt=ib/10
       filename="graphT."//char(48+ibt)//char(48+ib-ibt*10)
       endif

       open(10,file=filename)
       rewind(10)
       do ll=1,nn
       write(10,105) ak(ll),w(ll),dexp(w(ll)),Ew(ll)
       enddo
105   format(4(E13.6,1x))
       close(10)


        do i=0,200
        yy(i,ib)=y(i)
        enddo

        ibflag(ib)=1

       goto 1000
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccc
2000   continue


       do 3000 i1=1,1001
       E=E2+(E1-E2)*(i1-1)/1000.d0

       sum=0.d0
       do j=1,nstw

       do i=1,numw(j)-1
       if((E_linew(i,j).le.E.and.E_linew(i+1,j).gt.E).
     & or.(E_linew(i+1,j).lt.E.and.E_linew(i,j).ge.E)) then


        if(ibflag(j).eq.0) then
        write(6,*) "point not found (kpt,band), stop", i,j
        stop
        endif

         ik1=(ikpt_linew(i,j)-1)*200/(nkptw-1)    ! covert 50 grid to 100 grid
         ik2=ik1+2


       x1=(E_linew(i+1,j)-E)/(E_linew(i+1,j)-E_linew(i,j))
       x2=1.d0-x1

       ya=yy(ik1,j)*x1+yy(ik2,j)*x2


       if(ik1/200.d0.lt.ak_min(j).or.ik2/200.d0.gt.ak_max(j)) then
       if(ik1/200.d0.lt.ak_min(j)) ik_m=ak_min(j)*200.d0
       if(ik2/200.d0.gt.ak_max(j)) ik_m=ak_max(j)*200.d0
       ya_m=yy(ik_m,j)
       if(ya.lt.ya_m-2.d0) then
       ya=ya_m-2.d0       ! special treatment, don't trust the extrapolation
       write(6,*) "warning, outside,ya.lt.ya_m-2, ib,E,ak,ya",
     &      j,E,ik1/200.d0,ya
       endif

       if(ya.gt.ya_m+1.d0) then
       ya=ya_m+1.d0       
       write(6,*) "warning, outside,ya.gt.ya_m+1, ib,E,ak,ya",
     &      j,E,ik1/200.d0,ya
       write(37,*) "warning, outside,ya.gt.ya_m+1, ib,E,ak,ya",
     &      j,E,ik1/200.d0,ya
       endif  
       endif
       
       if(ya.ge.0.d0) then
       write(6,*) "warning, ya>0, set to 0", j,i,ya,E
       ya=0.d0
       endif


       sum=sum+dexp(ya)
       endif
       enddo
       enddo
       Et(i1)=E
       Tsum(i1)=sum
3000   continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccc total current J= G_0 *\int_Emin^Emax dE *T(E)
cccc G_0=2*e^2/h = (12.9 K Omu)^-1
cccc here, Emin=E_Fermi-dV, Emax=E_Fermi
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc here, all the energies are calliborated with right hand electrode level
ccccccccccccccccccccccccccccccccccccc

       Emax=E_Fermi
       Emin=Emax-dV

       sum=0.d0
       do i=1,1001
       if(Et(i).ge.Emin.and.Et(i).le.Emax) then
        if(Et(i-1).lt.Emin) then
        sum=sum+Tsum(i)*(Et(i)-Emin)
        endif
        if(Et(i+1).gt.Emax) then
        sum=sum+Tsum(i)*(Emax-Et(i))
        endif
        if(Et(i+1).le.Emax) then
        sum=sum+(Et(i+1)-Et(i))*(Tsum(i+1)+Tsum(i))/2
        endif
       endif
       enddo
   
       write(6,*) "J(in unit of G0*eV)=",sum
       

       open(10,file="graph.T")
       rewind(10)
       do i1=1,1001
       write(10,*) Et(i1),Et(i1)-E_Fermi,Tsum(i1)
       enddo
107   format(3(E13.6,1x))
       close(10)
       

       stop
       end

       

      

       

        
       
       

