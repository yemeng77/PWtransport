      program calc_IV
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc This program reads in the transmission T(E) for a given number
ccc of bias voltages [one graph.T file for each bias voltage]
ccc and calculates the I.vs.V, and dI/dV.vs.V curves
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      real*8 E0(1001,0:4),T(1001,0:4)
      real*8 V(501),aI(501),aR(501)

      mV=4

      do iV=0,mV
      open(11,file="graph.T."//char(48+iV)//"eV")
      rewind(11)
      do i=1,1001
      read(11,*) E1,E0(i,iV),T(i,iV)
      enddo
      close(11)
      enddo
      

      do 100 i=1,501
      dV=(i-1)*4/500.d0
      iV=dV*mV/4*0.99999

      x=dV-iV

      if(iV.gt.0.and.iV.lt.mV-1) then
      ym1=-0.5*x+x**2-0.5*x**3
      y0=1.d0-2.5*x**2+1.5*x**3
      y1=0.5*x+2*x**2-1.5*x**3
      y2=-0.5*x**2+0.5*x**3
      endif
    
      if(iV.eq.0) then
      ym1=0.d0
      y0=1.d0-1.5*x+0.5*x**2
      y1=2*x-x**2
      y2=-0.5*x+0.5*x**2
      endif

      if(iV.eq.mV-1) then
      ym1=-0.5*x+0.5*x**2
      y0=1.d0-x**2
      y1=0.5*x+0.5*x**2
      y2=0.d0
      endif

      sum=0.d0
      do j=1,1000

      E=-dV*(j-0.5d0)/1000.d0

      if(iV.gt.0) then
      ak=(E-E0(1,iV-1))/(E0(1001,iV-1)-E0(1,iV-1))*1000+1
      k=ak*0.9999999999999d0
      z1=k+1-ak
      z2=1.d0-z1
      sum=sum+(z1*T(k,iV-1)+z2*T(k+1,iV-1))*ym1
      endif

      ak=(E-E0(1,iV))/(E0(1001,iV)-E0(1,iV))*1000+1
      k=ak*0.9999999999999d0
      z1=k+1-ak
      z2=1.d0-z1
      sum=sum+(z1*T(k,iV)+z2*T(k+1,iV))*y0

      ak=(E-E0(1,iV+1))/(E0(1001,iV+1)-E0(1,iV+1))*1000+1
      k=ak*0.9999999999999d0
      z1=k+1-ak
      z2=1.d0-z1
      sum=sum+(z1*T(k,iV+1)+z2*T(k+1,iV+1))*y1

      if(iV.lt.mV-1) then
      ak=(E-E0(1,iV+2))/(E0(1001,iV+2)-E0(1,iV+2))*1000+1
      k=ak*0.9999999999999d0
      z1=k+1-ak
      z2=1.d0-z1
      sum=sum+(z1*T(k,iV+2)+z2*T(k+1,iV+2))*y2
      endif
      

      enddo
      sum=sum*dV/1000.d0
      V(i)=dV
      aI(i)=sum
100   continue

      do i=2,500
      aR(i)=(aI(i+1)-aI(i-1))/(V(i+1)-V(i-1))
      enddo
      aR(1)=aR(2)
      aR(501)=aR(500)

      open(10,file="IV.graph")
      rewind(10)
      do i=1,501
      write(10,300) V(i),aI(i),aR(i)
      enddo
      close(10)
300   format(3(E13.6,1x))
      stop
      end


      
      
      
