         subroutine global_sumc(sum_local_pe)
c--------------------------------------------------------

      implicit none
c
      include "mpif.h"
      complex*16    sum_local_pe,res
      integer ierr

      integer inode, nnodes, inode_tot, nnodes_tot
      integer nnodes_k, nnodes_b, inode_k, inode_b
      integer icolor_k, icolor_b, ikey_k, ikey_b, icolor, ikey
      integer num_group_k, num_group_b
      integer MPI_COMM_K1, MPI_COMM_K2, MPI_COMM_B1, MPI_COMM_B2
      integer MPI_COMM_K, MPI_COMM_N
      
      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &     nnodes_k, nnodes_b, inode_k, inode_b,
     &     icolor_k, icolor_b, ikey_k, ikey_b, icolor, ikey,
     &     num_group_k, num_group_b, 
     &     MPI_COMM_K1, MPI_COMM_K2, MPI_COMM_B1, MPI_COMM_B2,
     &     MPI_COMM_K, MPI_COMM_N
c     
      call mpi_allreduce(sum_local_pe,res,1,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)

      sum_local_pe = res
      return
      end


