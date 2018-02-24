define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
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
      subroutine adv_diff_wp_convective_op2d(
     &            q_data, q_gcw,
     &            u_data_0, u_data_1, u_gcw,
     &            r_data, r_gcw, d,
     &            ilower0, ilower1,
     &            iupper0, iupper1,
     &            dx,
     &            interp_coefs, smooth_weights, k)

      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER d

      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw),0:(d-1))

      REAL s_data_0(FACE2d0(ilower,iupper,0),0:1)
      REAL s_data_1(FACE2d1(ilower,iupper,0),0:1)

      INTEGER u_gcw
      REAL u_data_0(FACE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(FACE2d1(ilower,iupper,u_gcw))

      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:(d-1))

      REAL dx(0:1)

      INTEGER k, j
      INTEGER i0, i1
      REAL interp_coefs(0:k,0:(k-1))
      REAL smooth_weights(0:(k-1))


      do j=0,(d-1)
      call reconstruct_data_on_patch_2d(q_data(:,:,j), q_gcw,
     &             s_data_0, s_data_1, 0,
     &             ilower0, ilower1, iupper0, iupper1,
     &             interp_coefs, smooth_weights, k)

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
         r_data(i0,i1,j) =
     &     1.d0/dx(0)*(max(u_data_0(i0,i1),0.d0)*
     &     (s_data_0(i0,i1,1)-s_data_0(i0,i1,0))
     &     + min(u_data_0(i0+1,i1),0.d0)*
     &     (s_data_0(i0+1,i1,1)-s_data_0(i0+1,i1,0))
     &     + 0.5d0*(u_data_0(i0+1,i1)+u_data_0(i0,i1))*
     &     (s_data_0(i0+1,i1,0)-s_data_0(i0,i1,1)))

         r_data(i0,i1,j) = r_data(i0,i1,j) +
     &     1.d0/dx(1)*(max(u_data_1(i1,i0),0.d0)*
     &     (s_data_1(i1,i0,1)-s_data_1(i1,i0,0))
     &     + min(u_data_1(i1+1,i0),0.d0)*
     &     (s_data_1(i1+1,i0,1)-s_data_1(i1+1,i0,0))
     &     + 0.5d0*(u_data_1(i1+1,i0)+u_data_1(i1,i0))*
     &     (s_data_1(i1+1,i0,0)-s_data_1(i1,i0,1)))
        enddo
      enddo
      enddo
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
      subroutine reconstruct_data_on_patch_2d(q_data, q_gcw,
     &            r_data_0, r_data_1, r_gcw,
     &            ilower0, ilower1, iupper0, iupper1,
     &            interp_weights, smooth_weights, k)

       INTEGER k
       INTEGER ilower0, iupper0
       INTEGER ilower1, iupper1

       INTEGER q_gcw
       REAL q_data(CELL2d(ilower,iupper,q_gcw))

       INTEGER r_gcw
       REAL r_data_0(FACE2d0(ilower,iupper,r_gcw),0:1)
       REAL r_data_1(FACE2d1(ilower,iupper,r_gcw),0:1)

       REAL interp_weights(0:k,0:(k-1))
       REAL smooth_weights(0:k-1)

       REAL interp_values_p(0:k-1)
       REAL interp_values_n(0:k-1)
       REAL smooth_id_p(0:k-1)
       REAL smooth_id_n(0:k-1)
       REAL weights_p(0:k-1)
       REAL weights_n(0:k-1)
       REAL s_vals_x((ilower0-k):(iupper0+k),
     &               (ilower1):(iupper1))
       REAL s_vals_y((ilower1-k):(iupper1+k),
     &               (ilower0):(iupper0))

       INTEGER i0, i1
       INTEGER j, r

       REAL eps, total, alpha

      eps = 1.0d-7
c     X DIRECTION
      do i1 = ilower1, iupper1
        do i0=ilower0,iupper0+1
          do r=0,k-1
            interp_values_p(r) = 0.d0
            interp_values_n(r) = 0.d0
            do j=0,k-1
              interp_values_p(r) = interp_values_p(r)
     &          + interp_weights(r,j)*q_data(i0-r+j,i1)
              interp_values_n(r) = interp_values_n(r)
     &          + interp_weights(r+1,j)*q_data(i0-1-r+j,i1)
            enddo
          enddo
          smooth_id_p(0) = 13.d0/12.d0*(q_data(i0,i1)
     &         -2.d0*q_data(i0+1,i1)+q_data(i0+2,i1))**2
     &      + 0.25d0*(3.d0*q_data(i0,i1)-4.d0*q_data(i0+1,i1)
     &         +q_data(i0+2,i1))**2
          smooth_id_p(1) = 13.d0/12.d0*(q_data(i0-1,i1)
     &        -2.d0*q_data(i0,i1)+q_data(i0+1,i1))**2
     &      + 0.25d0*(q_data(i0-1,i1)-q_data(i0+1,i1))**2
          smooth_id_p(2) = 13.d0/12.d0*(q_data(i0-2,i1)
     &        -2.d0*q_data(i0-1,i1)+q_data(i0,i1))**2
     &      + 0.25d0*(3.d0*q_data(i0,i1)-4.d0*q_data(i0-1,i1)
     &        +q_data(i0-2,i1))**2

          smooth_id_n(0) = 13.d0/12.d0*(q_data(i0-1,i1)
     &        -2.d0*q_data(i0,i1)+q_data(i0+1,i1))**2
     &      + 0.25d0*(3.d0*q_data(i0-1,i1)-4.d0*q_data(i0,i1)
     &        +q_data(i0+1,i1))**2
          smooth_id_n(1) = 13.d0/12.d0*(q_data(i0-2,i1)
     &        -2.d0*q_data(i0-1,i1)+q_data(i0,i1))**2
     &      + 0.25d0*(q_data(i0-2,i1)-q_data(i0,i1))**2
          smooth_id_n(2) = 13.d0/12.d0*(q_data(i0-3,i1)
     &        -2.d0*q_data(i0-2,i1)+q_data(i0-1,i1))**2
     &      + 0.25d0*(3.d0*q_data(i0-1,i1)
     &        -4.d0*q_data(i0-2,i1)+q_data(i0-3,i1))**2

          total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(k-1-j)
     &           /(eps+smooth_id_p(j))**2
             total = total + alpha
             weights_p(j) = alpha
           enddo
           do j=0,k-1
             weights_p(j) = weights_p(j)/total
           enddo
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(j)
     &           /(eps+smooth_id_n(j))**2
             total = total + alpha
             weights_n(j) = alpha
           enddo
           do j=0,k-1
             weights_n(j) = weights_n(j)/total
           enddo
           r_data_0(i0,i1,0) = 0.d0
           r_data_0(i0,i1,1) = 0.d0
           do r=0,k-1
             r_data_0(i0,i1,0) = r_data_0(i0,i1,0)
     &               + weights_n(r)*interp_values_n(r)
             r_data_0(i0,i1,1) = r_data_0(i0,i1,1)
     &               + weights_p(r)*interp_values_p(r)
           enddo
         enddo
       enddo

c      Y DIRECTION
       do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
           do r=0,k-1
             interp_values_p(r) = 0.d0
             interp_values_n(r) = 0.d0
             do j=0,k-1
               interp_values_p(r) = interp_values_p(r)
     &          + interp_weights(r,j)*q_data(i0,i1+j-r)
               interp_values_n(r) = interp_values_n(r)
     &          + interp_weights(r+1,j)*q_data(i0,i1-1+j-r)
             enddo
           enddo
           smooth_id_p(0) = 13.d0/12.d0*(q_data(i0,i1)
     &        -2.d0*q_data(i0,i1+1)+q_data(i0,i1+2))**2
     &        +0.25d0*(3.d0*q_data(i0,i1)
     &        -4.d0*q_data(i0,i1+1)+q_data(i0,i1+2))**2
           smooth_id_p(1) = 13.d0/12.d0*(q_data(i0,i1-1)
     &        -2.d0*q_data(i0,i1)+q_data(i0,i1+1))**2
     &        +0.25d0*(q_data(i0,i1-1)-q_data(i0,i1+1))**2
           smooth_id_p(2) = 13.d0/12.d0*(q_data(i0,i1-2)
     &        -2.d0*q_data(i0,i1-1)+q_data(i0,i1))**2
     &        +0.25d0*(q_data(i0,i1-2)
     &        -4.d0*q_data(i0,i1-1)+3.d0*q_data(i0,i1))**2

           smooth_id_n(0) = 13.d0/12.d0*(q_data(i0,i1-1)
     &        -2.d0*q_data(i0,i1)+q_data(i0,i1+1))**2
     &        +0.25d0*(3.d0*q_data(i0,i1-1)
     &        -4.d0*q_data(i0,i1)+q_data(i0,i1+1))**2
           smooth_id_n(1) = 13.d0/12.d0*(q_data(i0,i1-2)
     &        -2.d0*q_data(i0,i1-1)+q_data(i0,i1))**2
     &        +0.25d0*(q_data(i0,i1-2)-q_data(i0,i1))**2
           smooth_id_n(2) = 13.d0/12.d0*(q_data(i0,i1-3)
     &        -2.d0*q_data(i0,i1-2)+q_data(i0,i1-1))**2
     &        +0.25d0*(q_data(i0,i1-3)
     &        -4.d0*q_data(i0,i1-2)+3.d0*q_data(i0,i1-1))**2
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(k-1-j)
     &           /((eps+smooth_id_p(j))**2)
             total = total + alpha
             weights_p(j) = alpha
           enddo
           do j=0,k-1
             weights_p(j) = weights_p(j)/total
           enddo
           total = 0.d0
           do j=0,k-1
             alpha = smooth_weights(j)
     &           /((eps+smooth_id_n(j))**2)
             total = total + alpha
             weights_n(j) = alpha
           enddo
           do j=0,k-1
             weights_n(j) = weights_n(j)/total
           enddo
           r_data_1(i1,i0,0) = 0.d0
           r_data_1(i1,i0,1) = 0.d0
           do r=0,k-1
             r_data_1(i1,i0,0) = r_data_1(i1,i0,0)
     &         + weights_n(r)*interp_values_n(r)
             r_data_1(i1,i0,1) = r_data_1(i1,i0,1)
     &         + weights_p(r)*interp_values_p(r)
           enddo
         enddo
       enddo
       end subroutine
