c
c     Routines to compute misc math operations on patches.
c
c     Created on 05 Jan 2004 by Boyce Griffith
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
c     Computes U = alpha V.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multiply13d(
     &     U,U_gcw,
     &     alpha,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw

      REAL alpha

      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the linear sum.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = alpha*V(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = alpha V + beta W.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multiplyadd13d(
     &     U,U_gcw,
     &     alpha,
     &     V,V_gcw,
     &     beta,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw,W_gcw

      REAL alpha,beta

      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL W(CELL3d(ilower,iupper,W_gcw))
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the linear sum.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = alpha*V(i0,i1,i2) + beta*W(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = A V.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multiply23d(
     &     U,U_gcw,
     &     A,A_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,A_gcw,V_gcw

      REAL A(CELL3d(ilower,iupper,A_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the linear sum.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = A(i0,i1,i2)*V(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = A V + beta W.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multiplyadd23d(
     &     U,U_gcw,
     &     A,A_gcw,
     &     V,V_gcw,
     &     beta,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,A_gcw,V_gcw,W_gcw

      REAL beta

      REAL A(CELL3d(ilower,iupper,A_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL W(CELL3d(ilower,iupper,W_gcw))
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the linear sum.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = A(i0,i1,i2)*V(i0,i1,i2) + beta*W(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = A V + B W.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multiplyadd33d(
     &     U,U_gcw,
     &     A,A_gcw,
     &     V,V_gcw,
     &     B,B_gcw,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,A_gcw,V_gcw,B_gcw,W_gcw

      REAL A(CELL3d(ilower,iupper,A_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL B(CELL3d(ilower,iupper,B_gcw))
      REAL W(CELL3d(ilower,iupper,W_gcw))
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the linear sum.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = A(i0,i1,i2)*V(i0,i1,i2)
     &              + B(i0,i1,i2)*W(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = |V|_1.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine pwl1norm3d(
     &     U,U_gcw,
     &     V,V_gcw,V_depth,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw,V_depth

      REAL V(CELL3d(ilower,iupper,V_gcw),0:V_depth-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      REAL n
      INTEGER i0,i1,i2,d
c
c     Compute the pointwise norm.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               n = 0.d0

               do d = 0,V_depth-1
                  n = n + dabs(V(i0,i1,i2,d))
               enddo

               U(i0,i1,i2) = n
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = |V|_2.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine pwl2norm3d(
     &     U,U_gcw,
     &     V,V_gcw,V_depth,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw,V_depth

      REAL V(CELL3d(ilower,iupper,V_gcw),0:V_depth-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      REAL n
      INTEGER i0,i1,i2,d
c
c     Compute the pointwise norm.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               n = 0.d0

               do d = 0,V_depth-1
                  n = n + V(i0,i1,i2,d)*V(i0,i1,i2,d)
               enddo

               U(i0,i1,i2) = dsqrt(n)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes U = |V|_oo.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine pwmaxnorm3d(
     &     U,U_gcw,
     &     V,V_gcw,V_depth,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw,V_depth

      REAL V(CELL3d(ilower,iupper,V_gcw),0:V_depth-1)
c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
c
c     Local variables.
c
      REAL n
      INTEGER i0,i1,i2,d
c
c     Compute the pointwise norm.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               n = 0.d0

               do d = 0,V_depth-1
                  n = dmax1(n,dabs(V(i0,i1,i2,d)))
               enddo

               U(i0,i1,i2) = n
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Helper subroutines for fast sweeping algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     Get neighboring points
      subroutine getneighbors(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     a,b,c,
     &     hx,hy,hz,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER i0,i1,i2
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw
      REAL dx(0:NDIM-1)

c
c     Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL a,b,c
      REAL hx,hy,hz
      
c     Get values for grid spacings and neighboring points.
      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
      b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
      c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
            
c     Take care of Dirichlet boundaries.
      if (a < 0.d0) then
         a  = 0.d0
         hx = hx/2.d0
      endif
            
      if (b < 0.d0) then
         b  = 0.d0
         hy = hy/2.d0
      endif
               
      if (c < 0.d0) then
         c  = 0.d0
         hz = hz/2.d0
      endif
      
      
      
      return
      end
c
c     Sort the neighboring points
      subroutine sortneighbors(
     &     a,b,c,
     &     hx,hy,hz,
     &     a1,a2,a3,
     &     h1,h2,h3)
c
      implicit none
c
c     Input.
c
      REAL a,b,c
      REAL hx,hy,hz
      
c
c     Output.
c
      REAL a1,a2,a3
      REAL h1,h2,h3
      
c     Sort a,b,c
      if (a .le. b) then
         if (a .le. c) then
            if(b .le. c) then
               a1 = a; a2 = b; a3 = c
               h1 = hx; h2 = hy; h3 = hz
            else
               a1 = a; a2 = c; a3 = b
               h1 = hx; h2 = hz; h3 = hy
            endif
         else
            a1 = c; a2 = a; a3 = b
            h1 = hz; h2 = hx; h3 = hy
         endif
      else 
         if (b .le. c) then
            if (a .le. c) then
               a1 = b; a2 = a; a3 = c
               h1 = hy; h2 = hx; h3 = hz
            else
               a1 = b; a2 = c; a3 = a
               h1 = hy; h2 = hz; h3 = hx
            endif
         else
            a1 = c; a2 = b; a3 = a
            h1 = hz; h2 = hy; h3 = hx
         endif
      endif
      
      return
      end
c      
c     Find dbar
      subroutine computedbar(
     &     dbar,
     &     a1,a2,a3,
     &     h1,h2,h3)
c
      implicit none

c
c     Input.
c
      REAL a1,a2,a3
      REAL h1,h2,h3
      
c
c     Output.
c
      REAL dbar

c
c     Local variables.
c
      REAL Q,R,S
      REAL dtil
      
      
c     Algorithm to find dbar
      dtil = a1 + h1
      if (dtil .le. a2) then
         dbar = dtil
      else
         Q = h1*h1 + h2*h2
         R = -2.d0*(h2*h2*a1 + h1*h1*a2)
         S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
         dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
      endif
              
      if (dtil .lt. a3) then
         dbar = dtil
              
      else
         Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
         R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
         S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
         dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
         dbar = dtil
      endif
                    
      return
      end
               

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out fast sweeping algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fastsweep3d(
     &     U,U_gcw,
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
      INTEGER U_gcw

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    a,b,c
      REAL    a1,a2,a3
      REAL    hx,hy,hz
      REAL    h1,h2,h3
      REAL    dbar
c
c     Fast-sweeping algorithm
c

      
c     Do the eight sweeping directions
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = ilower2,iupper2
         do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = iupper2,ilower2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = iupper2,ilower2,-1
         do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = iupper2,ilower2,-1
         do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = iupper2,ilower2,-1
         do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
              
c              Get values for grid spacings and neighboring points.
               call getneighbors (U,U_gcw,
     &                            ilower0,iupper0,
     &                            ilower1,iupper1,
     &                            ilower2,iupper2,
     &                            i0,i1,i2,a,b,c,hx,hy,hz,dx)
               
c             Sort a,b,c.
              call sortneighbors (a,b,c,
     &                            hx,hy,hz,
     &                            a1,a2,a3,
     &                            h1,h2,h3)
              
c             Compute dbar
              call computedbar(dbar,a1,a2,a3,h1,h2,h3)
                    
                        
               U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      
      return
      end
c
c