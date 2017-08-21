c
c     Routines to compute quantities related to variable coefficient
c     generalized Laplace operators.
c
c     Created on 12 Jul 2016 by Nishant Nangia
c
c     Copyright (c) 2002-2014, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes (f0,f1,f2) = alpha div mu (2*tau)
c     where tau = 1/2(grad(u) + grad(u)^T) is passed in
c
c     Computes the side-centered variable coefficient generalized
c     Laplacian, with cell and edge-centered coefficient mu, cell-centered
c     tauxx, tauyy, tauzz, and edge-centered tauxy, tauxz, tauyz
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vclaplace3d(
     &     f0,f1,f2,f_gcw,
     &     mu_cc,mu_cc_gcw,
     &     mu_ec0, mu_ec1, mu_ec2 ,mu_ec_gcw,
     &     tauxx_cc,tauyy_cc,tauzz_cc,tau_cc_gcw,
     &     tauxy_ec,tauxz_ec,tauyz_ec,tau_ec_gcw,
     &     alpha,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER f_gcw,mu_cc_gcw,mu_ec_gcw
      INTEGER tau_cc_gcw,tau_ec_gcw

      REAL alpha

      REAL mu_cc(CELL3d(ilower,iupper,mu_cc_gcw))
      REAL tauxx_cc(CELL3d(ilower,iupper,tau_cc_gcw))
      REAL tauyy_cc(CELL3d(ilower,iupper,tau_cc_gcw))
      REAL tauzz_cc(CELL3d(ilower,iupper,tau_cc_gcw))
      
      REAL mu_ec0(EDGE3d0(ilower,iupper,mu_ec_gcw))
      REAL mu_ec1(EDGE3d1(ilower,iupper,mu_ec_gcw))
      REAL mu_ec2(EDGE3d2(ilower,iupper,mu_ec_gcw))
      REAL tauyz_ec(EDGE3d0(ilower,iupper,tau_ec_gcw))
      REAL tauxz_ec(EDGE3d1(ilower,iupper,tau_ec_gcw))
      REAL tauxy_ec(EDGE3d2(ilower,iupper,tau_ec_gcw))
      
      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL f0(SIDE3d0(ilower,iupper,f_gcw))
      REAL f1(SIDE3d1(ilower,iupper,f_gcw))
      REAL f2(SIDE3d2(ilower,iupper,f_gcw))
      
      
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL fac0,fac1,fac2
      REAL d_dx,d_dy,d_dz
c
c     Compute the discrete divergence of mu (tau).
c
      fac0 = alpha*2.d0/(dx(0))
      fac1 = alpha*2.d0/(dx(1))
      fac2 = alpha*2.d0/(dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
c              Compute each component of divergence
               d_dx = fac0*(mu_cc(i0,i1,i2)*tauxx_cc(i0,i1,i2)
     &                     -mu_cc(i0-1,i1,i2)*tauxx_cc(i0-1,i1,i2))
               
               d_dy = fac1*(mu_ec2(i0,i1+1,i2)*tauxy_ec(i0,i1+1,i2)
     &                     -mu_ec2(i0,i1,i2)*tauxy_ec(i0,i1,i2))
     
               d_dz = fac2*(mu_ec1(i0,i1,i2+1)*tauxz_ec(i0,i1,i2+1)
     &                     -mu_ec1(i0,i1,i2)*tauxz_ec(i0,i1,i2))
     
c              Set force to be summation of divergence components
               f0(i0,i1,i2) = d_dx + d_dy + d_dz
            enddo
         enddo
      enddo

      
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
c              Compute each component of divergence
               d_dx = fac0*(mu_ec2(i0+1,i1,i2)*tauxy_ec(i0+1,i1,i2)
     &                     -mu_ec2(i0,i1,i2)*tauxy_ec(i0,i1,i2))
               
               d_dy = fac1*(mu_cc(i0,i1,i2)*tauyy_cc(i0,i1,i2)
     &                     -mu_cc(i0,i1-1,i2)*tauyy_cc(i0,i1-1,i2))
     
               d_dz = fac2*(mu_ec0(i0,i1,i2+1)*tauyz_ec(i0,i1,i2+1)
     &                     -mu_ec0(i0,i1,i2)*tauyz_ec(i0,i1,i2))
     
c              Set force to be summation of divergence components
               f1(i0,i1,i2) = d_dx + d_dy + d_dz
            enddo
         enddo
      enddo
      
      
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
c              Compute each component of divergence
               d_dx = fac0*(mu_ec1(i0+1,i1,i2)*tauxz_ec(i0+1,i1,i2)
     &                     -mu_ec1(i0,i1,i2)*tauxz_ec(i0,i1,i2))
               
               d_dy = fac1*(mu_ec0(i0,i1+1,i2)*tauyz_ec(i0,i1+1,i2)
     &                     -mu_ec0(i0,i1,i2)*tauyz_ec(i0,i1,i2))
     
               d_dz = fac2*(mu_cc(i0,i1,i2)*tauzz_cc(i0,i1,i2)
     &                     -mu_cc(i0,i1,i2-1)*tauzz_cc(i0,i1,i2-1))
     
c              Set force to be summation of divergence components
               f2(i0,i1,i2) = d_dx + d_dy + d_dz
            enddo
         enddo
      enddo
c
      return
      end
c
