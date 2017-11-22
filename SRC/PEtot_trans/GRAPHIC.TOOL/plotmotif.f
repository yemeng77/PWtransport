       program plotmotif
       implicit double precision (a-h,o-z)
       real*4,allocatable,dimension(:,:,:) :: dens_m
       real*8 AL_mbox(3,3)
       integer iat_neigh(10)
       real*8 x_neigh(3,10)
       character*2 atomtype

       character*30 filename
       write(6,*) "input the motif file name"
       read(5,*) filename
       open(10,file=filename,form="unformatted")
       rewind(10)
       read(10) num_neigh,rad_box
       read(10) AL_mbox
       do ia=1,num_neigh
       read(10) iat_neigh(ia),x_neigh(1,ia),x_neigh(2,ia),
     &  x_neigh(3,ia)
       enddo
       read(10) m,sum2
       allocate(dens_m(-m:m,-m:m,-m:m))
       read(10) dens_m
       close(10)

       open(12,file='motif.plt_f')
       rewind(12)
       write(12,99) 3,3
       write(12,98) m+1,m+1,m+1
        
99     format(2(i2,2x))
98     format(3(i5,2x))
       zmax=AL_mbox(3,3)*0.529177/2
       ymax=AL_mbox(2,2)*0.529177/2
       xmax=AL_mbox(1,1)*0.529177/2
       write(12,200)-zmax,zmax,-ymax,ymax,-xmax,xmax
       write(6,200)-zmax,zmax,-ymax,ymax,-xmax,xmax
200    format(6(E13.5,1x))
       write(12,100) (((dens_m(i,j,k),i=-m,m,2),j=-m,m,2),
     &      k=-m,m,2)
100    format(6(E13.5,1x))
       close(12)

       open(11,file="motif.xyz")
       rewind(11)
       write(11,*) num_neigh
c       write(11,*)
       af=0.529177
       do i=1,num_neigh
       x=AL_mbox(1,1)*x_neigh(1,i)+AL_mbox(1,2)*x_neigh(2,i)+
     &   AL_mbox(1,3)*x_neigh(3,i)
       y=AL_mbox(2,1)*x_neigh(1,i)+AL_mbox(2,2)*x_neigh(2,i)+
     &   AL_mbox(2,3)*x_neigh(3,i)
       z=AL_mbox(3,1)*x_neigh(1,i)+AL_mbox(3,2)*x_neigh(2,i)+
     &   AL_mbox(3,3)*x_neigh(3,i)
       if(i.eq.1) then
       x0=x
       y0=y
       z0=z
       endif
       call atom_name(atomtype,iat_neigh(i))
       write(11,500) atomtype ,(x-x0)*af,(y-y0)*af,(z-z0)*af
       enddo
         close(11)
500     format(a2,3x,3(f13.7,1x))
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


      

       
   
