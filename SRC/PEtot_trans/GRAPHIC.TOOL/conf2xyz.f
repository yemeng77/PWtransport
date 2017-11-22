         program conf2xyz
         implicit double precision (a-h,o-z)
         real*8 AL(3,3),AL0(3,3),ALI(3,3)
         character*50 fname_in,fname_out
         character*2 atomtype(100000),atom
         real*8 xyz(3,100000)

         write(6,*) "this program converts .config to .xyz file"

         write(6,*) "input the name of input xxx.config file" 
         read(5,*) fname_in
         write(6,*) "input the name of output xxx.xyz file"
         read(5,*) fname_out


         open(10,file=fname_in)
         rewind(10)
         read(10,*) num
         read(10,*) AL(1,1),AL(2,1),AL(3,1)
         read(10,*) AL(1,2),AL(2,2),AL(3,2)
         read(10,*) AL(1,3),AL(2,3),AL(3,3)


         write(6,*) "***** current AL is *******"
         write(6,201) AL(1,1),AL(2,1),AL(3,1)
         write(6,201) AL(1,2),AL(2,2),AL(3,2)
         write(6,201) AL(1,3),AL(2,3),AL(3,3)
201      format(3(f12.5,1x))

         write(6,*) "input 1: use AL as the box; 2: new box" 
         read(5,*) iflag

         if(iflag.eq.1) then
         x10=0.d0
         x20=0.d0
         x30=0.d0
         AL0=AL
         else
         write(6,*) 
     & "input the corner position of new box as x1,x2,x3 of AL"
         read(5,*) x10,x20,x30
         write(6,*) "input AL0(:,1) of the new box in a.u"
         read(5,*) AL0(1,1),AL0(2,1),AL0(3,1)
         write(6,*) "input AL0(:,2) of the new box in a.u"
         read(5,*) AL0(1,2),AL0(2,2),AL0(3,2)
         write(6,*) "input AL0(:,3) of the new box in a.u"
         read(5,*) AL0(1,3),AL0(2,3),AL0(3,3)
         endif

         write(6,*) 
     &"input 1: output x,y,z coord; 2 output x1,x2,x3 of the box"
         write(6,*) 
     & "Note: for x1,x2,x3, the plot will assume the box is orth"
         write(6,*) 
     & "To overlay with charge in gopenmol, x1,x2,x3 need to be used"
         read(5,*) iflag2

         call get_ALI(AL,ALI)
         

         x1max=-10
         x2max=-10
         x3max=-10
         x1min=10
         x2min=10
         x3min=10
         do i=0,1
         do j=0,1
         do k=0,1
         x=AL0(1,1)*i+AL0(1,2)*j+AL0(1,3)*k
         y=AL0(2,1)*i+AL0(2,2)*j+AL0(2,3)*k
         z=AL0(3,1)*i+AL0(3,2)*j+AL0(3,3)*k
        
         x1=x*ALI(1,1)+y*ALI(2,1)+z*ALI(3,1)+x10
         x2=x*ALI(1,2)+y*ALI(2,2)+z*ALI(3,2)+x20
         x3=x*ALI(1,3)+y*ALI(2,3)+z*ALI(3,3)+x30
         if(x1.gt.x1max) x1max=x1
         if(x2.gt.x2max) x2max=x2
         if(x3.gt.x3max) x3max=x3
         if(x1.lt.x1min) x1min=x1
         if(x2.lt.x2min) x2min=x2
         if(x3.lt.x3min) x3min=x3
         enddo
         enddo
         enddo
         i1max=x1max+0.99999d0
         i2max=x2max+0.99999d0
         i3max=x3max+0.99999d0
         i1min=x1min-0.99999d0
         i2min=x2min-0.99999d0
         i3min=x3min-0.99999d0
         write(6,*) "i1,2,3max: ", i1max,i2max,i3max
         write(6,*) "i1,2,3min: ", i1min,i2min,i3min


         call get_ALI(AL0,ALI)

      xmax=dsqrt(AL0(1,1)**2+AL0(2,1)**2+AL0(3,1)**2)
      ymax=dsqrt(AL0(1,2)**2+AL0(2,2)**2+AL0(3,2)**2)
      zmax=dsqrt(AL0(1,3)**2+AL0(2,3)**2+AL0(3,3)**2)
         
         num_count=0
         do 300 i=1,num
         read(10,*) iatom, x1t,x2t,x3t


         do 200 i1=i1min,i1max
         do 200 i2=i2min,i2max
         do 200 i3=i3min,i3max

         x1=x1t-x10+i1
         x2=x2t-x20+i2
         x3=x3t-x30+i3
         
         x=AL(1,1)*x1+AL(1,2)*x2+AL(1,3)*x3
         y=AL(2,1)*x1+AL(2,2)*x2+AL(2,3)*x3
         z=AL(3,1)*x1+AL(3,2)*x2+AL(3,3)*x3

         x1n=x*ALI(1,1)+y*ALI(2,1)+z*ALI(3,1)
         x2n=x*ALI(1,2)+y*ALI(2,2)+z*ALI(3,2)
         x3n=x*ALI(1,3)+y*ALI(2,3)+z*ALI(3,3)


         if(x1n.ge.0.d0.and.x1n.le.1.d0.and.x2n.ge.0.d0.and.
     &  x2n.le.1.d0.and.x3n.ge.0.d0.and.x3n.le.1.d0) then


         num_count=num_count+1
         call atom_name(atomtype(num_count),iatom)
         if(iflag2.eq.2) then
         xyz(1,num_count)=x1n*xmax
         xyz(2,num_count)=x2n*ymax
         xyz(3,num_count)=x3n*zmax
         else
         x=AL0(1,1)*x1n+AL0(1,2)*x2n+AL0(1,3)*x3n
         y=AL0(2,1)*x1n+AL0(2,2)*x2n+AL0(2,3)*x3n
         z=AL0(3,1)*x1n+AL0(3,2)*x2n+AL0(3,3)*x3n
         xyz(1,num_count)=x
         xyz(2,num_count)=y
         xyz(3,num_count)=z
         endif

         endif

200     continue
300     continue

         write(6,*) "num of atom in the plot", num_count

         open(11,file=fname_out)
         rewind(11)
         write(11,*) num_count
         write(11,*) 
         af=0.529177
         do i=1,num_count
         write(11,100) atomtype(i),xyz(1,i)*af,xyz(2,i)*af,xyz(3,i)*af
         enddo
         close(11)
100     format(a2,3x,3(f13.7,1x))
         stop
         end



         subroutine atom_name(atomtype,iatom)
         implicit double precision (a-h,o-z)
         character*2 atomtype
         integer iatom

         atomtype="XX"

         if(iatom.eq.1) atomtype="H"
         if(iatom.eq.2) atomtype="He"
         if(iatom.eq.5) atomtype="B"
         if(iatom.eq.6) atomtype="C"
         if(iatom.eq.7) atomtype="N"
         if(iatom.eq.8) atomtype="O"
         if(iatom.eq.9) atomtype="F"
         if(iatom.eq.12) atomtype="Mg"
         if(iatom.eq.13) atomtype="Al"
         if(iatom.eq.14) atomtype="Si"
         if(iatom.eq.15) atomtype="P"
         if(iatom.eq.16) atomtype="S"
         if(iatom.eq.17) atomtype="Cl"
         if(iatom.eq.20) atomtype="Ca"
         if(iatom.eq.25) atomtype="Mn"
         if(iatom.eq.26) atomtype="Fe"
         if(iatom.eq.27) atomtype="Co"
         if(iatom.eq.28) atomtype="Ni"
         if(iatom.eq.29) atomtype="Cu"
         if(iatom.eq.30) atomtype="Zn"
         if(iatom.eq.30) atomtype="Zn"
         if(iatom.eq.31) atomtype="Ga"
         if(iatom.eq.32) atomtype="Ge"
         if(iatom.eq.33) atomtype="As"
         if(iatom.eq.34) atomtype="Se"
         if(iatom.eq.35) atomtype="Br"
         if(iatom.eq.47) atomtype="Ag"
         if(iatom.eq.48) atomtype="Cd"
         if(iatom.eq.49) atomtype="In"
         if(iatom.eq.50) atomtype="Sn"
         if(iatom.eq.51) atomtype="Sb"
         if(iatom.eq.52) atomtype="Te"
         if(iatom.eq.53) atomtype="I"
         if(iatom.eq.79) atomtype="Au"
         if(iatom.eq.105) atomtype="H"
         if(iatom.eq.115) atomtype="H"
         if(iatom.eq.101) atomtype="H"
         if(iatom.eq.102) atomtype="H"
         if(iatom.eq.100) atomtype="H"

         if(atomtype.eq."XX") then
         write(6,*) "atom type ", iatom, "not found in the program"
         write(6,*) "add this type in the program"
         stop
         endif

         return
         end




        

        
      
        
