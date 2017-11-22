program test_zgemm

 double complex, dimension(5,2) :: a 
 double complex, dimension(5,3) :: b
 double complex, dimension(2,3) :: c
 integer :: i, j

 do i=1, 2
    do j=1, 5
       a(j,i) = (1.d0,0.d0)*j + (i-1)*(5.d0,0.d0)
    enddo
 enddo

 do i=1, 3
    do j=1, 5
       b(j,i) = (0.1d0,0.d0)*j + (i-1)*(0.5d0,0.d0)    
    enddo
 enddo

 print*, a(:,1)
 print*, a(:,2)
 print*, ''

 print*, b(:,1)
 print*, b(:,2)
 print*, b(:,3)
 print*, ''

 call zgemm('c', 'n', 2, 3, 4, (1.d0,0.d0), a, 5, b, 5, (0.d0,0.d0), c, 2)

 print*, c(1,:)
 print*, c(2,:)
  
end program test_zgemm

