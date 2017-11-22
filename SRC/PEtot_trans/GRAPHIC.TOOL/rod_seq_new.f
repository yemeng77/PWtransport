      program rod_seq
      implicit double precision (a-h,o-z)

      parameter (natom=200000)

      real*8 xatom(3,natom),watom(natom)
      real*8 ALL(3,3),ALI(3,3),AL2(3,3),AL2I(3,3)
      real*8 xh(3,3)

      integer iatom(natom),iseq(1000),index(natom)
      integer neigh(natom)

      real*8 xyz_nb(3,4,natom)
      real*8 xyz_mb1(3),xyz_mb2(3),xyz_mb(3,2)
      real*8, allocatable,  dimension(:,:) ::  xyz_a,x123_a
      integer,  allocatable,  dimension(:) ::  iat_a,nneigh_a
      integer,  allocatable,  dimension(:,:) ::  ineigh_a


ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
cccc In this version, we only deal with ZB and WZ structures, please check sequence
cccc Now, the surface Cd-H, Se-H bond directions in ZB and WZ are the same. 
ccccccccccccccccccccccccccccccccccccccccccccc

      open(10,file="sequence")
      rewind(10)
      read(10,*) nseq,nc,dc,hh,iZB     ! nc: center, dc: rod diam, hh: rod high
      do i=1,nseq
      read(10,*) iseq(i)
      enddo
      close(10)
      if(iZB.eq.1) write(6,*) "calc. ZB structure"
      if(iZB.eq.2) write(6,*) "calc. WZ structure"

      dc_st=dc
      hh_st=hh

  

      if(iZB.eq.1) then
      alat=11.49165277

      ALL(1,1)=alat/dsqrt(2.d0)
      ALL(2,1)=0.d0
      ALL(3,1)=0.d0

      ALL(1,2)=-alat/dsqrt(2.d0)/2
      ALL(2,2)=dsqrt(3.d0)*alat/dsqrt(2.d0)/2
      ALL(3,2)=0.d0

      ALL(1,3)=0.d0
      ALL(2,3)=0.d0
      ALL(3,3)=alat/dsqrt(3.d0)

      xh(1,1)=0.d0
      xh(2,1)=0.d0
      xh(3,1)=3.d0/8       ! the position at the center of the vertical bond

      xh(1,2)=1.d0/3.d0
      xh(2,2)=2.d0/3.d0
      xh(3,2)=3.d0/8

      xh(1,3)=2.d0/3.d0
      xh(2,3)=1.d0/3.d0
      xh(3,3)=3.d0/8

      dh0=3.d0/8.d0    ! ideal, when i1-i2-i3, i1.ne.i3, ZB structure
      endif

      if(iZB.eq.2) then      ! experimental WZ structure
      alat=11.49165277

      ALL(1,1)=alat/dsqrt(2.d0)
      ALL(2,1)=0.d0
      ALL(3,1)=0.d0

      ALL(1,2)=-alat/dsqrt(2.d0)/2
      ALL(2,2)=dsqrt(3.d0)*alat/dsqrt(2.d0)/2
      ALL(3,2)=0.d0

      ALL(1,3)=0.d0
      ALL(2,3)=0.d0
c      ALL(3,3)=alat/dsqrt(3.d0)
      ALL(3,3)=13.24887d0/2             ! experimental length

      xh(1,1)=0.d0
      xh(2,1)=0.d0
      xh(3,1)=3.d0/8       ! the position at the center of the vertical bond

      xh(1,2)=1.d0/3.d0
      xh(2,2)=2.d0/3.d0
      xh(3,2)=3.d0/8

      xh(1,3)=2.d0/3.d0
      xh(2,3)=1.d0/3.d0
      xh(3,3)=3.d0/8

c      dh1=0.3685       ! Wurtzite, fitted when i1-i2-i3, i1.eq.i3 
      dh1=0.36790       ! Wurtzite,  experimental value

      facH_v=alat/dsqrt(3.d0)/ALL(3,3)*0.375d0/dh1
      facH_nv=alat/dsqrt(3.d0)/ALL(3,3)*
     &         (0.5d0-0.375d0)/(0.5d0-dh1)
      endif



ccccccccccccccccccccccccccccccccccccccccccccccccc

      dc=dc*ALL(3,3)
      hh=hh*ALL(3,3)

      nn1=2*dc/(alat/dsqrt(2.d0))*2/dsqrt(3.d0)+0.5+2
      nn2=nn1
      nn3=(2*dc+hh)/ALL(3,3)+0.5+2


       nn1=nn1+2
       nn2=nn2+1
       nn3=nn3+1

c      write(6,*) "nn3,nn1,nn2=",nn3,nn1,nn2
c      write(6,*) "unit: one double layer;xy unit cell length"  
c      write(6,*) "input the new nn3,nn1,nn2"
c      read(5,*) nn3,nn1,nn2
      

      xc=ALL(1,1)*nn1/2+ALL(1,2)*nn2/2+ALL(1,3)*nc
      yc=ALL(2,1)*nn1/2+ALL(2,2)*nn2/2+ALL(2,3)*nc
      zc=ALL(3,1)*nn1/2+ALL(3,2)*nn2/2+ALL(3,3)*nc
      zc1=zc-hh/2
      zc2=zc+hh/2

      num_a=(1+nn1)*(1+nn2)*nseq*2  
      allocate(xyz_a(3,num_a))
      allocate(x123_a(3,num_a))
      allocate(iat_a(num_a))
      allocate(nneigh_a(num_a))
      allocate(ineigh_a(4,num_a))
      

      num=0
      do 100 k=2,nseq-1     !  one, double layer
      do 100 ik=-1,1,2      ! one double layer has two atoms, one for Cd, one for Se
      if(iseq(k-1).eq.iseq(k+1)) dh=dh1     ! Wurtrize
      if(iseq(k-1).ne.iseq(k+1)) dh=dh0     ! ZB situation
      do 100 i=0,nn1
      do 100 j=0,nn2
      x1=xh(1,iseq(k))+i
      x2=xh(2,iseq(k))+j
      x3=xh(3,iseq(k))+k+ik*dh
      x=ALL(1,1)*x1+ALL(1,2)*x2+ALL(1,3)*x3
      y=ALL(2,1)*x1+ALL(2,2)*x2+ALL(2,3)*x3
      z=ALL(3,1)*x1+ALL(3,2)*x2+ALL(3,3)*x3
      num=num+1
      if(ik.eq.-1) iat_a(num)=34
      if(ik.eq.1) iat_a(num)=48
      x123_a(1,num)=x1/nn1
      x123_a(2,num)=x2/nn2
      x123_a(3,num)=(x3-nc+nn3/2.d0)/nn3
      xyz_a(1,num)=x
      xyz_a(2,num)=y
      xyz_a(3,num)=z
100   continue

ccccccccccccccccccccccccccccccccccccccccc
      num_tot=num
      ddc=(1.25*alat*dsqrt(3.d0)/4)**2

      do i=1,num_tot
      num=0
      do j=1,num_tot
      if(i.ne.j.and.iat_a(i).ne.iat_a(j)) then
      dd=(xyz_a(1,i)-xyz_a(1,j))**2+(xyz_a(2,i)-xyz_a(2,j))**2
     &  +(xyz_a(3,i)-xyz_a(3,j))**2
      if(dd.lt.ddc) then
      num=num+1
      ineigh_a(num,i)=j
      endif
      endif
      enddo
      nneigh_a(i)=num
      if(num.gt.4) then
      write(6,*) "num of neigh.gt.4, stop", i, num
      stop
      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccc


      num=0
      do 120 i=1,num_tot
      x=xyz_a(1,i)
      y=xyz_a(2,i)
      z=xyz_a(3,i)

      dd0=(x-xc)**2+(y-yc)**2
      dd1=(x-xc)**2+(y-yc)**2+(z-zc1)**2
      dd2=(x-xc)**2+(y-yc)**2+(z-zc2)**2

      if((z.ge.zc1.and.z.le.zc2.and.dd0.lt.dc**2).or.
     &  (z.gt.zc2.and.dd2.lt.dc**2).or.
     &  (z.lt.zc1.and.dd1.lt.dc**2)) then 

      num=num+1
      if(num.gt.natom) then
      write(6,*) "num.gt.natom,stop",num,natom
      stop
      endif
      xatom(1,num)=x123_a(1,i)
      xatom(2,num)=x123_a(2,i)
      xatom(3,num)=x123_a(3,i)
      iatom(num)=iat_a(i)
      index(num)=i
      watom(num)=1.d0
      endif
120   continue

      write(6,*) "num inside the QD geometry=",num


      do i=1,3
      ALL(i,1)=ALL(i,1)*nn1
      ALL(i,2)=ALL(i,2)*nn2
      ALL(i,3)=ALL(i,3)*nn3
      enddo
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc remove the site which has only one bond
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ddc=(1.25*alat*dsqrt(3.d0)/4)**2
1000  continue
      do 60 i=1,num
      neigh(i)=0
      do 59 j=1,num
      if(j.eq.i) goto 59
      dx1=xatom(1,j)-xatom(1,i)
      dx2=xatom(2,j)-xatom(2,i)
      dx3=xatom(3,j)-xatom(3,i)
      dx=ALL(1,1)*dx1+ALL(1,2)*dx2+ALL(1,3)*dx3
      dy=ALL(2,1)*dx1+ALL(2,2)*dx2+ALL(2,3)*dx3
      dz=ALL(3,1)*dx1+ALL(3,2)*dx2+ALL(3,3)*dx3
      dd=dx**2+dy**2+dz**2
      if(dd.lt.ddc) then
      neigh(i)=neigh(i)+1
      if(neigh(i).gt.4) then
      write(6,*) "neigh.gt.4, stop", i, neigh(i)
      stop
      endif
      nb=neigh(i)
      xyz_nb(1,nb,i)=dx
      xyz_nb(2,nb,i)=dy
      xyz_nb(3,nb,i)=dz
      endif

59    continue
60    continue

      nchangetot=0
      num2=0
      do i=1,num
      if(neigh(i).gt.1) then
      num2=num2+1
      xatom(1,num2)=xatom(1,i)
      xatom(2,num2)=xatom(2,i)
      xatom(3,num2)=xatom(3,i)
      iatom(num2)=iatom(i)
      watom(num2)=watom(i)
      neigh(num2)=neigh(i)
      index(num2)=index(i)
      xyz_nb(:,:,num2)=xyz_nb(:,:,i)
      else
      nchangetot=nchangetot+1
      endif
      enddo
cccccccccccccccccccccccc
      num=num2
      if(nchangetot.gt.0) goto 1000
      write(6,*) "num after remove one bond Cd/Se: ",num


      call get_ALI(ALL,ALI)


      write(6,*) "add surface Cd (1) or not (0)" 
      read(5,*) iadd_Cd

      if(iadd_Cd.eq.0) goto 900
ccccccccccccccccccccccccccccccccccccccc
cccc  add the Cd atom on the surface Se atom
ccccccccccccccccccccccccccccccccccccccc
      num2=num
      do i=1,num
      if(iatom(i).eq.34.and.neigh(i).lt.4) then 
      i2=index(i)
      if(nneigh_a(i2).ne.4) then
      write(6,*) "nneigh_a.ne.4, increase nn1,nn2,nn3",nneigh_a(i2),i2,i
      stop
      endif
      do ii=1,nneigh_a(i2)
       j2=ineigh_a(ii,i2)
       iflag=0
       do j=1,num2    ! add
       if(index(j).eq.j2) iflag=1
       enddo
       if(iflag.eq.0)  then      ! this Cd atom should be included, but not included yet
       num2=num2+1
      xatom(1,num2)=x123_a(1,j2)
      xatom(2,num2)=x123_a(2,j2)
      xatom(3,num2)=x123_a(3,j2)
      iatom(num2)=iat_a(j2)
      watom(num2)=1.d0
      index(num2)=j2
cc      neigh(num2)=neigh(i)     ! there is no information for the neigh for the new one!
       endif
      enddo
      endif
      enddo
ccccccccccccccccccccccccccccccccccccc
      num=num2
      write(6,*) "num Cd+Se atom after add surface Cd num=",num
900   continue
ccccccccccccccccccccccccccccccccccccccc
cccc  add the surface H atom 
ccccccccccccccccccccccccccccccccccccccc
      num2=num
      dx_Cd=0.d0
      dy_Cd=0.d0
      dz_Cd=0.d0
      dx_Se=0.d0
      dy_Se=0.d0
      dz_Se=0.d0
      do 90 i=1,num
ccccccccc  first, recalculate the neigh
      i2=index(i)

      if(nneigh_a(i2).ne.4) then
      write(6,*) "nneigh_a.ne.4, increase nn1,nn2,nn3",nneigh_a(i2),i2,i
      stop
      endif

      do 90 ii=1,nneigh_a(i2)     ! check for each of its neighbore
       j2=ineigh_a(ii,i2)
       iflag=0
       do j=1,num    
       if(index(j).eq.j2) iflag=1
       enddo
       if(iflag.eq.0)  then      ! This neighbore is not in the QD, attach a H on it
       num2=num2+1
         
       if(iZB.eq.1) then
       xatom(1,num2)=x123_a(1,i2)+(x123_a(1,j2)-x123_a(1,i2))*0.5d0
       xatom(2,num2)=x123_a(2,i2)+(x123_a(2,j2)-x123_a(2,i2))*0.5d0
       xatom(3,num2)=x123_a(3,i2)+(x123_a(3,j2)-x123_a(3,i2))*0.5d0
       endif
       if(iZB.eq.2) then
ccccccccc  make the WZ surface has the same as the ZB surface Cd-H, Se-H bond orientations and length
       xatom(1,num2)=x123_a(1,i2)+(x123_a(1,j2)-x123_a(1,i2))*0.5d0
       xatom(2,num2)=x123_a(2,i2)+(x123_a(2,j2)-x123_a(2,i2))*0.5d0
         if(abs(x123_a(1,j2)-x123_a(1,i2))*nn1+
     &         abs(x123_a(2,j2)-x123_a(2,i2))*nn2.gt.0.001) then
ccccccccc  this is not a vertical bond
       xatom(3,num2)=x123_a(3,i2)+(x123_a(3,j2)-x123_a(3,i2))*
     &               0.5d0*facH_nv
         else
ccccccccc  this is a vertical bond
       xatom(3,num2)=x123_a(3,i2)+(x123_a(3,j2)-x123_a(3,i2))*
     &               0.5d0*facH_v
         endif
       endif

       if(iat_a(i2).eq.34) iatom(num2)=105
       if(iat_a(i2).eq.48) iatom(num2)=115
       watom(num2)=1.d0
       index(num2)=j2
cc      neigh(num2)=neigh(i)     ! there is no information for the neigh for the new one!
cccccccccccccccccccccc TEST TEST
        x1=xatom(1,num2)-x123_a(1,i2)
        x2=xatom(2,num2)-x123_a(2,i2)
        x3=xatom(3,num2)-x123_a(3,i2)
        x=ALL(1,1)*x1+ALL(1,2)*x2+ALL(1,3)*x3
        y=ALL(2,1)*x1+ALL(2,2)*x2+ALL(2,3)*x3
        z=ALL(3,1)*x1+ALL(3,2)*x2+ALL(3,3)*x3
        if(iat_a(i2).eq.34) then
        dx_Se=dx_Se+x
        dy_Se=dy_Se+y
        dz_Se=dz_Se+z
        endif
        if(iat_a(i2).eq.48) then
        dx_Cd=dx_Cd+x
        dy_Cd=dy_Cd+y
        dz_Cd=dz_Cd+z
        endif
        
       endif 
90     continue
ccccccccccccccccccccccccccccccccccccc

      write(6,*) "total number of atoms: ",num2
      write(6,*) "surface H dxyz summation" 
      write(6,*) "Cd: dx,dy,dz ", dx_Cd,dy_Cd,dz_Cd
      write(6,*) "Se: dx,dy,dz ", dx_Se,dy_Se,dz_Se

cccccccc  change a box, to orthogonal box for 3DF calculation

      AL2=0.d0
      if(iZB.eq.1) then
      AL2(3,1)=alat/dsqrt(3.d0)
      else
      AL2(3,1)=13.24887d0/2             ! experimental length
      endif

      AL2(1,2)=alat*3/(2*dsqrt(2.d0))
      AL2(2,2)=-alat*dsqrt(3.d0)/(2*dsqrt(2.d0))

      AL2(2,3)=alat*dsqrt(3.d0)/dsqrt(2.d0)

      write(6,*) "nn1,nn2,nn3=",nn3-2,(nn1-1)/dsqrt(3.d0),
     &      (nn1-1)/dsqrt(3.d0)
      write(6,*) "unit: one double layer;xy unit cell length"  
      write(6,*) "input the new nn1,nn2,nn3"
      read(5,*) nn1,nn2,nn3

      AL2(:,1)=AL2(:,1)*nn1
      AL2(:,2)=AL2(:,2)*nn2
      AL2(:,3)=AL2(:,3)*nn3

      call get_ALI(AL2,AL2I)

      do i=1,num2
      x1=xatom(1,i)-0.5
      x2=xatom(2,i)-0.5
      x3=xatom(3,i)-0.5
      x=ALL(1,1)*x1+ALL(1,2)*x2+ALL(1,3)*x3
      y=ALL(2,1)*x1+ALL(2,2)*x2+ALL(2,3)*x3
      z=ALL(3,1)*x1+ALL(3,2)*x2+ALL(3,3)*x3
      xatom(1,i)=AL2I(1,1)*x+AL2I(2,1)*y+AL2I(3,1)*z+0.5
      xatom(2,i)=AL2I(1,2)*x+AL2I(2,2)*y+AL2I(3,2)*z+0.5
      xatom(3,i)=AL2I(1,3)*x+AL2I(2,3)*y+AL2I(3,3)*z+0.5
      enddo


      num_Cd=0
      num_Se=0
      num_115=0
      num_105=0
      do i=1,num2
      if(iatom(i).eq.34) num_Se=num_Se+1
      if(iatom(i).eq.48) num_Cd=num_Cd+1
      if(iatom(i).eq.115) num_115=num_115+1
      if(iatom(i).eq.105) num_105=num_105+1
      enddo

      write(6,*) "num_Cd,Se,115,105: ", 
     & num_Cd,num_Se,num_115,num_105




cc      call reorder(iatom,xatom,num2) 

      open(10,file="xatom.config")
      rewind(10)
      write(10,*) num2,num_Cd,num_Se,num_115,num_105,dc_st,hh_st,iZB,
     &  nn1,nn2,nn3
      write(10,200) AL2(1,1),AL2(2,1),AL2(3,1)
      write(10,200) AL2(1,2),AL2(2,2),AL2(3,2)
      write(10,200) AL2(1,3),AL2(2,3),AL2(3,3)
      do i=1,num2
      write(10,300) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i),
     &   0,0,0
      enddo

      close(10)

200   format(3(f14.7,1x))
300   format(i4,2x,3(f13.9,1x),2x,3(i3,1x))
      stop
      end

      
