      subroutine map_distort(wv_n,wvd_n,mr_n)
      implicit double precision (a-h,o-z)
      include "param.escan_real"

      real*8 wv_n(mr_n),wvd_n(mr_n)

      real*8 xgrid(3,mr_n) 
