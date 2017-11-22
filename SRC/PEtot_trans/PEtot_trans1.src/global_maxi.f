c------------------------------------------------------
      subroutine global_maxi(inum)
c------------------------------------------------------
c output=inum 
c this routine calculates the global max of inum 
c answer put on each PE
c
c
      implicit none
c
      include "mpif.h"
      
      integer inum,inum_tmp,ierr

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

      
      call mpi_allreduce(inum,inum_tmp,1,
     &  MPI_INTEGER,MPI_MAX,MPI_COMM_K,ierr)

      inum=inum_tmp

      return
      end

