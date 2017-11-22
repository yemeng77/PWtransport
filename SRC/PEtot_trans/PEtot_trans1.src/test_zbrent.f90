program test_zbrent
  implicit none
  double precision :: f, zbrent
  external f, zbrent

  print*,  zbrent(f,-0.5d0,1.5d0,1.d-5)
 
  print*, f(0.5d0) 

end program test_zbrent




