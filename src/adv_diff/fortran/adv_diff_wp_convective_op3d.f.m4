define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes r = u.grad(q)
c
c     where u is vector valued face centered velocity
c     q is cell centered with depth d
c     returns r_data at cell centeres
c     computes grad(q) using weno + wave propagation
c     interpolation coefficients and weights must be provided
c     currently only works for interp orders 3 (k=2) and 5 (k=3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adv_diff_wp_convective_op3d(
     &            q_data, q_gcw,
     &            u_data_0, u_data_1, u_data_2, u_gcw,
     &            r_data, r_gcw, depth,
     &            ilower0, ilower1, ilower2,
     &            iupper0, iupper1, iupper2,
     &            dx,
     &            interp_coefs, smooth_weights, k)
      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      INTEGER depth

      INTEGER q_gcw
      REAL q_data(CELL3d(ilower,iupper,q_gcw),0:(depth-1))

      REAL s_data_0(FACE3d0(ilower,iupper,0),0:1)
      REAL s_data_1(FACE3d1(ilower,iupper,0),0:1)
      REAL s_data_2(FACE3d2(ilower,iupper,0),0:1)

      INTEGER u_gcw
      REAL u_data_0(FACE3d0(ilower,iupper,u_gcw))
      REAL u_data_1(FACE3d1(ilower,iupper,u_gcw))
      REAL u_data_2(FACE3d2(ilower,iupper,u_gcw))

      INTEGER r_gcw
      REAL r_data(CELL3d(ilower,iupper,r_gcw),0:(depth-1))

      REAL dx(0:2)

      INTEGER k, j
      REAL interp_coefs(0:k,0:(k-1))
      REAL smooth_weights(0:(k-1))

      INTEGER i0, i1, i2
      REAL inv_dx, inv_dy, inv_dz

      do j=0,(depth-1)
      call reconstruct_data_on_patch_3d(q_data(:,:,:,j), q_gcw,
     &             s_data_0, s_data_1, s_data_2, 0,
     &             ilower0, ilower1, ilower2,
     &             iupper0, iupper1, iupper2,
     &             interp_coefs, smooth_weights, k)
      inv_dx = 1.d0/dx(0)
      inv_dy = 1.d0/dx(1)
      inv_dz = 1.d0/dx(2)
      do i2 = ilower2, iupper2
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
         r_data(i0,i1,i2,j) =
     &     inv_dx*(max(u_data_0(i0,i1,i2),0.d0)*
     &     (s_data_0(i0,i1,i2,1)-s_data_0(i0,i1,i2,0))
     &     + min(u_data_0(i0+1,i1,i2),0.d0)*
     &     (s_data_0(i0+1,i1,i2,1)-s_data_0(i0+1,i1,i2,0))
     &     + 0.5d0*(u_data_0(i0+1,i1,i2)+u_data_0(i0,i1,i2))*
     &     (s_data_0(i0+1,i1,i2,0)-s_data_0(i0,i1,i2,1)))

         r_data(i0,i1,i2,j) = r_data(i0,i1,i2,j) +
     &     inv_dy*(max(u_data_1(i1,i2,i0),0.d0)*
     &     (s_data_1(i1,i2,i0,1)-s_data_1(i1,i2,i0,0))
     &     + min(u_data_1(i1+1,i2,i0),0.d0)*
     &     (s_data_1(i1+1,i2,i0,1)-s_data_1(i1+1,i2,i0,0))
     &     + 0.5d0*(u_data_1(i1+1,i2,i0)+u_data_1(i1,i2,i0))*
     &     (s_data_1(i1+1,i2,i0,0)-s_data_1(i1,i2,i0,1)))

         r_data(i0,i1,i2,j) = r_data(i0,i1,i2,j) +
     &     inv_dz*(max(u_data_2(i2,i0,i1),0.d0)*
     &     (s_data_2(i2,i0,i1,1)-s_data_2(i2,i0,i1,0))
     &     + min(u_data_2(i2+1,i0,i1),0.d0)*
     &     (s_data_2(i2+1,i0,i1,1)-s_data_2(i2+1,i0,i1,0))
     &     + 0.5d0*(u_data_2(i2+1,i0,i1)+u_data_2(i2,i0,i1))*
     &     (s_data_2(i2+1,i0,i1,0)-s_data_2(i2,i0,i1,1)))
      enddo; enddo; enddo; enddo
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Reconstructs data on patches using a weno scheme
c       the convex and interpolation weights must be supplied
c
c       q_data is cell centered with depth 1
c       r_data_* are face centered with depth 2
c         and return the values reconstructed from each side
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reconstruct_data_on_patch_3d(q_data, q_gcw,
     &            r_data_0, r_data_1, r_data_2, r_gcw,
     &            ilower0, ilower1, ilower2,
     &            iupper0, iupper1, iupper2,
     &            interp_coefs, smooth_weights, k)

      implicit none
      INTEGER k
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      INTEGER q_gcw
      REAL q_data(CELL3d(ilower,iupper,q_gcw))

      INTEGER r_gcw
      REAL r_data_0(FACE3d0(ilower,iupper,r_gcw),0:1)
      REAL r_data_1(FACE3d1(ilower,iupper,r_gcw),0:1)
      REAL r_data_2(FACE3d2(ilower,iupper,r_gcw),0:1)

      REAL interp_coefs(0:k,0:(k-1))
      REAL smooth_weights(0:(k-1))

      REAL interp_values_p(0:(k-1))
      REAL interp_values_n(0:(k-1))
      REAL smooth_id_p(0:(k-1))
      REAL smooth_id_n(0:(k-1))
      REAL weights_p(0:(k-1))
      REAL weights_n(0:(k-1))
      REAL interp_values(0:(k-1))
      REAL smooth_id(0:(k-1))
      REAL weights(0:(k-1))

      INTEGER i0, i1, i2
      INTEGER j, r

      REAL eps, total, alpha
      eps = 1.0d-7

c     X-Direction
      do i2 = ilower2, iupper2; do i1 = ilower1, iupper1
        do i0=ilower0,iupper0+1
          do r=0,k-1
            interp_values_p(r) = 0.d0
            interp_values_n(r) = 0.d0
            do j=0,k-1
              interp_values_p(r) = interp_values_p(r)
     &          + interp_coefs(r,j)*q_data(i0-r+j,i1,i2)
              interp_values_n(r) = interp_values_n(r)
     &          + interp_coefs(r+1,j)*q_data(i0-1-r+j,i1,i2)
            enddo
          enddo
          smooth_id_p(0) = 13.d0/12.d0*(q_data(i0,i1,i2)
     &         -2.d0*q_data(i0+1,i1,i2)+q_data(i0+2,i1,i2))**2
     &      + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0+1,i1,i2)
     &         +q_data(i0+2,i1,i2))**2
          smooth_id_p(1) = 13.d0/12.d0*(q_data(i0-1,i1,i2)
     &        -2.d0*q_data(i0,i1,i2)+q_data(i0+1,i1,i2))**2
     &      + 0.25d0*(q_data(i0-1,i1,i2)-q_data(i0+1,i1,i2))**2
          smooth_id_p(2) = 13.d0/12.d0*(q_data(i0-2,i1,i2)
     &        -2.d0*q_data(i0-1,i1,i2)+q_data(i0,i1,i2))**2
     &      + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0-1,i1,i2)
     &        +q_data(i0-2,i1,i2))**2

          smooth_id_n(0) = 13.d0/12.d0*(q_data(i0-1,i1,i2)
     &        -2.d0*q_data(i0,i1,i2)+q_data(i0+1,i1,i2))**2
     &      + 0.25d0*(3.d0*q_data(i0-1,i1,i2)-4.d0*q_data(i0,i1,i2)
     &        +q_data(i0+1,i1,i2))**2
          smooth_id_n(1) = 13.d0/12.d0*(q_data(i0-2,i1,i2)
     &        -2.d0*q_data(i0-1,i1,i2)+q_data(i0,i1,i2))**2
     &      + 0.25d0*(q_data(i0-2,i1,i2)-q_data(i0,i1,i2))**2
          smooth_id_n(2) = 13.d0/12.d0*(q_data(i0-3,i1,i2)
     &        -2.d0*q_data(i0-2,i1,i2)+q_data(i0-1,i1,i2))**2
     &      + 0.25d0*(3.d0*q_data(i0-1,i1,i2)
     &        -4.d0*q_data(i0-2,i1,i2)+q_data(i0-3,i1,i2))**2

          total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(k-1-j)/(eps+smooth_id_p(j))**2
             total = total + alpha
             weights_p(j) = alpha
           enddo
           do j=0,k-1
             weights_p(j) = weights_p(j)/total
           enddo
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(j)/(eps+smooth_id_n(j))**2
             total = total + alpha
             weights_n(j) = alpha
           enddo
           do j=0,k-1
             weights_n(j) = weights_n(j)/total
           enddo
           r_data_0(i0,i1,i2,0) = 0.d0
           r_data_0(i0,i1,i2,1) = 0.d0
           do r=0,k-1
             r_data_0(i0,i1,i2,0) = r_data_0(i0,i1,i2,0)
     &               + weights_n(r)*interp_values_n(r)
             r_data_0(i0,i1,i2,1) = r_data_0(i0,i1,i2,1)
     &               + weights_p(r)*interp_values_p(r)
           enddo
         enddo
       enddo; enddo

c       SECOND INTERPOLANT
c       Y DIRECTION
c     Interpolate in other direction
      do i2 = ilower2, iupper2; do i1 = ilower1, iupper1+1
        do i0=ilower0,iupper0
          do r=0,k-1
            interp_values_p(r) = 0.d0
            interp_values_n(r) = 0.d0
            do j=0,k-1
              interp_values_p(r) = interp_values_p(r)
     &          + interp_coefs(r,j)*q_data(i0,i1-r+j,i2)
              interp_values_n(r) = interp_values_n(r)
     &          + interp_coefs(r+1,j)*q_data(i0,i1-1-r+j,i2)
            enddo
          enddo
          smooth_id_p(0) = 13.d0/12.d0*(q_data(i0,i1,i2)
     &          -2.d0*q_data(i0,i1+1,i2)+q_data(i0,i1+2,i2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0,i1+1,i2)
     &          +q_data(i0,i1+2,i2))**2
          smooth_id_p(1) = 13.d0/12.d0*(q_data(i0,i1-1,i2)
     &          -2.d0*q_data(i0,i1,i2)+q_data(i0,i1+1,i2))**2
     &       + 0.25d0*(q_data(i0,i1-1,i2)-q_data(i0,i1+1,i2))**2
          smooth_id_p(2) = 13.d0/12.d0*(q_data(i0,i1-2,i2)
     &          -2.d0*q_data(i0,i1-1,i2)+q_data(i0,i1,i2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0,i1-1,i2)
     &          +q_data(i0,i1-2,i2))**2

          smooth_id_n(0) = 13.d0/12.d0*(q_data(i0,i1-1,i2)
     &          -2.d0*q_data(i0,i1,i2)+q_data(i0,i1+1,i2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1-1,i2)-4.d0*q_data(i0,i1,i2)
     &          +q_data(i0,i1+1,i2))**2
          smooth_id_n(1) = 13.d0/12.d0*(q_data(i0,i1-2,i2)
     &          -2.d0*q_data(i0,i1-1,i2)+q_data(i0,i1,i2))**2
     &       + 0.25d0*(q_data(i0,i1-2,i2)-q_data(i0,i1,i2))**2
          smooth_id_n(2) = 13.d0/12.d0*(q_data(i0,i1-3,i2)
     &          -2.d0*q_data(i0,i1-2,i2)+q_data(i0,i1-1,i2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1-1,i2)-4.d0*q_data(i0,i1-2,i2)
     &          +q_data(i0,i1-3,i2))**2

          total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(k-1-j)/(eps+smooth_id_p(j))**2
             total = total + alpha
             weights_p(j) = alpha
           enddo
           do j=0,k-1
             weights_p(j) = weights_p(j)/total
           enddo
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(j)/(eps+smooth_id_n(j))**2
             total = total + alpha
             weights_n(j) = alpha
           enddo
           do j=0,k-1
             weights_n(j) = weights_n(j)/total
           enddo
           r_data_1(i1,i2,i0,0) = 0.d0
           r_data_1(i1,i2,i0,1) = 0.d0
           do r=0,k-1
             r_data_1(i1,i2,i0,0) = r_data_1(i1,i2,i0,0)
     &               + weights_n(r)*interp_values_n(r)
             r_data_1(i1,i2,i0,1) = r_data_1(i1,i2,i0,1)
     &               + weights_p(r)*interp_values_p(r)
           enddo
         enddo
       enddo; enddo

c       THIRD INTERPOLANT
c       Z DIRECTION
c     Interpolate in other direction
      do i2 = ilower2, iupper2+1; do i1 = ilower1, iupper1
        do i0=ilower0,iupper0
          do r=0,k-1
            interp_values_p(r) = 0.d0
            interp_values_n(r) = 0.d0
            do j=0,k-1
              interp_values_p(r) = interp_values_p(r)
     &          + interp_coefs(r,j)*q_data(i0,i1,i2-r+j)
              interp_values_n(r) = interp_values_n(r)
     &          + interp_coefs(r+1,j)*q_data(i0,i1,i2-1-r+j)
            enddo
          enddo
          smooth_id_p(0) = 13.d0/12.d0*(q_data(i0,i1,i2)
     &          -2.d0*q_data(i0,i1,i2+1)+q_data(i0,i1,i2+2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0,i1,i2+1)
     &          +q_data(i0,i1,i2+2))**2
          smooth_id_p(1) = 13.d0/12.d0*(q_data(i0,i1,i2-1)
     &          -2.d0*q_data(i0,i1,i2)+q_data(i0,i1,i2+1))**2
     &       + 0.25d0*(q_data(i0,i1,i2-1)-q_data(i0,i1,i2+1))**2
          smooth_id_p(2) = 13.d0/12.d0*(q_data(i0,i1,i2-2)
     &          -2.d0*q_data(i0,i1,i2-1)+q_data(i0,i1,i2))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2)-4.d0*q_data(i0,i1,i2-1)
     &          +q_data(i0,i1,i2-2))**2

          smooth_id_n(0) = 13.d0/12.d0*(q_data(i0,i1,i2-1)
     &          -2.d0*q_data(i0,i1,i2)+q_data(i0,i1,i2+1))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2-1)-4.d0*q_data(i0,i1,i2)
     &          +q_data(i0,i1,i2+1))**2
          smooth_id_n(1) = 13.d0/12.d0*(q_data(i0,i1,i2-2)
     &          -2.d0*q_data(i0,i1,i2-1)+q_data(i0,i1,i2))**2
     &       + 0.25d0*(q_data(i0,i1,i2-2)-q_data(i0,i1,i2))**2
          smooth_id_n(2) = 13.d0/12.d0*(q_data(i0,i1,i2-3)
     &          -2.d0*q_data(i0,i1,i2-2)+q_data(i0,i1,i2-1))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1,i2-1)-4.d0*q_data(i0,i1,i2-2)
     &          +q_data(i0,i1,i2-3))**2

          total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(k-1-j)/(eps+smooth_id_p(j))**2
             total = total + alpha
             weights_p(j) = alpha
           enddo
           do j=0,k-1
             weights_p(j) = weights_p(j)/total
           enddo
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(j)/(eps+smooth_id_n(j))**2
             total = total + alpha
             weights_n(j) = alpha
           enddo
           do j=0,k-1
             weights_n(j) = weights_n(j)/total
           enddo
           r_data_2(i2,i0,i1,0) = 0.d0
           r_data_2(i2,i0,i1,1) = 0.d0
           do r=0,k-1
             r_data_2(i2,i0,i1,0) = r_data_2(i2,i0,i1,0)
     &               + weights_n(r)*interp_values_n(r)
             r_data_2(i2,i0,i1,1) = r_data_2(i2,i0,i1,1)
     &               + weights_p(r)*interp_values_p(r)
           enddo
         enddo
       enddo; enddo
      endsubroutine

