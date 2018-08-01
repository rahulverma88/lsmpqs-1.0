c***********************************************************************
c
c  File:        qss_spatial_derivatives2d.f
c***********************************************************************
c  Algorithms to compute the 2D spatial derivative terms in the level 
c  set equation.
c 
c***********************************************************************

c***********************************************************************
c  Computes the 2D variational normal velocity term, but with constant
c  theta.
c***********************************************************************
      subroutine qss2dSetVarNorm(
     &  mask,
     &  mask_x,
     &  mask_y,
     &  var_a,
     &  a0,
     &  ca,
     &  ilo_gb, 
     &  ihi_gb,
     &  jlo_gb, 
     &  jhi_gb,
     &  ilo_fb,
     &  ihi_fb,
     &  jlo_fb,
     &  jhi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real mask(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)    
      real mask_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real mask_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real var_a(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)      
      real dx
      real ca
      real a0
      integer i, j
      real h1, h2, s, C2, C1, eps_var, absgrad, x
      eps_var = 1.5 * dx      
      C1 = 0.04
      C2 = 1
      
c       loop over included cells {
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
              call qssHeaviside(h1, -mask(i,j), eps_var)
              call qssHeaviside(h2, mask(i,j), eps_var)
              absgrad = mask_x(i,j) * mask_x(i,j) 
     &                     + mask_y(i,j) * mask_y(i,j)
                absgrad = SQRT(absgrad)
                x = mask(i,j)*mask(i,j) + absgrad*dx*absgrad*dx
                x = SQRT(x)
                s = mask(i,j)/x
                var_a(i,j) = h1 * a0 - s * h2 * C2 * COS(ca)*absgrad
          enddo
        enddo
c       } end loop over grid
        return
      end
c***********************************************************************

c***********************************************************************
c  Computes the 2D variational normal velocity term, but with spatially
c  varying theta.
c***********************************************************************
      subroutine qss2dSetVarNormTheta(
     &  mask,
     &  mask_x,
     &  mask_y,
     &  var_a,
     &  a0,
     &  ca,
     &  ilo_gb, 
     &  ihi_gb,
     &  jlo_gb, 
     &  jhi_gb,
     &  ilo_fb,
     &  ihi_fb,
     &  jlo_fb,
     &  jhi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real mask(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)    
      real mask_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real mask_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real var_a(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)      
      real dx
      real ca(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
     
      real a0
      integer i, j
      real h1, h2, s, C2, C1, eps_var, absgrad, x
      eps_var = 1.5 * dx      
      C1 = 0.04
      C2 = 1
      
c       loop over included cells {
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb
              call qssHeaviside(h1, -mask(i,j), eps_var)
              call qssHeaviside(h2, mask(i,j), eps_var)
              absgrad = mask_x(i,j) * mask_x(i,j) 
     &                     + mask_y(i,j) * mask_y(i,j)
                absgrad = SQRT(absgrad)
                x = mask(i,j)*mask(i,j) + absgrad*dx*absgrad*dx
                x = SQRT(x)
                s = mask(i,j)/x
                var_a(i,j) = h1 * a0 - s * h2 * C2 *
     &            COS(ca(i,j))*absgrad
          enddo
        enddo
c       } end loop over grid
        return
      end
c***********************************************************************

c***********************************************************************
c  Computes the 2D variational curvature and convective velocity terms
c***********************************************************************
      subroutine qss2dSetVarCurvAdv(
     &  mask,
     &  mask_x,
     &  mask_y,
     &  var_b,
     &  vel_x,
     &  vel_y,
     &  b_max_over_dx,
     &  max_U_over_dx,
     &  const_b,
     &  ilo_gb, 
     &  ihi_gb,
     &  jlo_gb, 
     &  jhi_gb,
     &  ilo_fb,
     &  ihi_fb,
     &  jlo_fb,
     &  jhi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real mask(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)      
      real mask_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real mask_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real var_b(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)      
      real vel_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real vel_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb) 
      real dx
      
      real const_b
      
      integer i, j
      real h1, h2, s, C2, C1, eps_var, x, absgrad
      real b_max_over_dx, max_U_over_dx, b_max, U_over_dX_cur
      real inv_dx, inv_dy

      eps_var = 1.5 * dx      
      C1 = const_b
      C2 = 1
      
c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0
      
c     initialize b_max to -1
      b_max = -1.0d0  

c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dx
      
c       loop over included cells {
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
                call qssHeaviside(h1, -mask(i,j), eps_var)
                call qssHeaviside(h2, mask(i,j), eps_var)
                absgrad = mask_x(i,j) * mask_x(i,j) +
     &                     mask_y(i,j) * mask_y(i,j) 
                absgrad = SQRT(absgrad)
                x = mask(i,j)*mask(i,j) + absgrad*dx*absgrad*dx
                x = SQRT(x)
                s = mask(i,j)/x
                var_b(i,j) = h1*C1
c                h2*C1
                
                if (var_b(i,j) > b_max) then
                    b_max = var_b(i,j)
                endif
                
                vel_x(i,j) = s * C2 * mask_x(i,j) * h2
                vel_y(i,j) = s * C2 * mask_y(i,j) * h2
                
                U_over_dX_cur = abs(vel_x(i,j))*inv_dx
     &                        + abs(vel_y(i,j))*inv_dy

                if (U_over_dX_cur .gt. max_U_over_dX) then
                  max_U_over_dX = U_over_dX_cur
                endif  
            enddo
          enddo
c       } end loop over grid

          b_max_over_dx = 2 * b_max * 
     &               (inv_dx*inv_dx + inv_dy*inv_dy )
          return
      end
c***********************************************************************

c***********************************************************************
c  Computes the right hand side of the 2D level set equation, 
c  with second order accurate stencil for  variational a, b and 
c  convective velocity terms
c***********************************************************************
      subroutine qss2dGetRHS2nd(
     &  phi,
     &  vel_n,
     &  b,
     &  vel_x, vel_y,
     &  max_H_over_dX,
     &  lse_rhs,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      
c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus
      real phi_x_minus, phi_y_minus
      
c     upwinded derivatives
      real phi_x_adv, phi_y_adv
    
c     central difference derivatives
      real phi_x, phi_y
      
c     second order central difference derivatives
      real phi_xx, phi_yy, phi_xy
      
      real curv, grad_mag2
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_n(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real b(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real vel_y(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb)
      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real max_H_over_dX
      real H_over_dX_cur
      real max_dx_sq
      real norm_grad_phi, phi_x_cur, phi_y_cur
      real temp_norm, temp_curv, temp_vel

      real dx_factor, dy_factor

c     compute denominator values, for curvature term
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      
c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c           phi_x_plus
            D1_x = phi(i,j) - phi(i-1,j)
            D1_x_plus = phi(i+1,j) - phi(i,j)
            D2_x = phi(i+1,j) + phi(i-1,j) - 2*phi(i,j)
            D2_x_plus = phi(i+2,j) + phi(i,j) - 2*phi(i+1,j)
            D2_x_minus = phi(i-2,j) + phi(i,j) - 2*phi(i-1,j)
            
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x)*inv_dx
            else
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x_plus)*inv_dx
            endif

c           phi_x_minus
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = (D1_x 
     &                           + half*D2_x_minus)*inv_dx
            else
              phi_x_minus = (D1_x
     &                           + half*D2_x)*inv_dx
            endif

c           phi_y_plus
            D1_y = phi(i,j) - phi(i,j-1)
            D1_y_plus = phi(i,j+1) - phi(i,j)
            D2_y = phi(i,j+1) + phi(i,j-1) - 2*phi(i,j)
            D2_y_plus = phi(i,j+2) + phi(i,j) - 2*phi(i,j+1)
            D2_y_minus = phi(i,j-2) + phi(i,j) - 2*phi(i,j-1)
            
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y)*inv_dy
            else
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y_plus)*inv_dy
            endif

c           phi_y_minus
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = (D1_y 
     &                           + half*D2_y_minus)*inv_dy
            else
              phi_y_minus = (D1_y
     &                           + half*D2_y)*inv_dy
            endif
            
            vel_n_cur = vel_n(i,j)
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi


c             compute contribution to lse_rhs(i,j) 
              temp_norm = - vel_n_cur*sqrt(norm_grad_phi_sq)

            endif
        
c           get max_H_over_dX for calculating stable dt
c           initialize max_H_over_dX to -1
            max_H_over_dX = -1.0d0

c           compute max_dx_sq
            max_dx_sq = max(dx,dy)
            max_dx_sq = max(dx,dy) * max(dx,dy)
      
            phi_x_cur = max(abs(phi_x_plus),
     &                          abs(phi_x_minus))
            phi_y_cur = max(abs(phi_y_plus),
     &                          abs(phi_y_minus))
            norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                            + phi_y_cur*phi_y_cur 
     &                            + max_dx_sq )

            H_over_dX_cur = abs(vel_n(i,j)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

            if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
            endif
        
            
    
c           compute advective term contributions
c              phi_x
              if (vel_x(i,j).ge.0) then
                phi_x_adv = phi_x_minus
              else
                phi_x_adv = phi_x_plus
              endif

c              phi_y
              if (vel_y(i,j).ge.0) then
                phi_y_adv = phi_y_minus
              else
                phi_y_adv = phi_y_plus
              endif    
              
              
            temp_vel = - ( vel_x(i,j)*phi_x_adv
     &                       + vel_y(i,j)*phi_y_adv )
     
     
c         compute curvature term
            phi_x = (phi(i+1,j) - phi(i-1,j))*dx_factor
            phi_y = (phi(i,j+1) - phi(i,j-1))*dy_factor
               
            phi_xx = (phi(i+2,j) + phi(i-2,j)
     &               - 2*phi(i,j))*dx_factor*dx_factor
            phi_yy = (phi(i,j+2) + phi(i,j-2)
     &               - 2*phi(i,j))*dy_factor*dy_factor
        
            phi_xy = (phi(i+1,j+1) + phi(i-1,j-1)
     &               - phi(i-1,j+1) - phi(i+1,j-1))
     &               *dx_factor*dy_factor
    
c           compute squared magnitude of gradient
            grad_mag2 = phi_x * phi_x
     &                + phi_y * phi_y
            if (grad_mag2 .lt. zero_tol) then
              curv = 0.d0
            else
              curv = phi_xx*phi_y*phi_y 
     &             +   phi_yy*phi_x*phi_x  
     &             - 2*phi_xy*phi_x*phi_y
              curv = curv / grad_mag2 
            endif

            temp_curv = b(i,j)*curv

	    lse_rhs(i,j) = temp_norm + temp_curv + temp_vel    
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c  Compute the right hand side of the 2D level set equation, with 
c  second order accurate stencils, and constant zero contact angle 
c  everywhere - so a and b are not varying
c  and the convective velocity term does not exist.
c***********************************************************************
      subroutine qss2dGetRHS2ndconst(
     &  phi,
     &  vel_n,
     &  b,
     &  max_H_over_dX,
     &  lse_rhs,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      
c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus
      real phi_x_minus, phi_y_minus
      
c     upwinded derivatives
      real phi_x_adv, phi_y_adv
    
c     central difference derivatives
      real phi_x, phi_y
      
c     second order central difference derivatives
      real phi_xx, phi_yy, phi_xy
      
      real curv, grad_mag2
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_n
      real b
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb)
      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real max_H_over_dX
      real H_over_dX_cur
      real max_dx_sq
      real norm_grad_phi, phi_x_cur, phi_y_cur
      real temp_norm, temp_curv, temp_vel

      real dx_factor, dy_factor

c     compute denominator values, for curvature term
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      
c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c           phi_x_plus
            D1_x = phi(i,j) - phi(i-1,j)
            D1_x_plus = phi(i+1,j) - phi(i,j)
            D2_x = phi(i+1,j) + phi(i-1,j) - 2*phi(i,j)
            D2_x_plus = phi(i+2,j) + phi(i,j) - 2*phi(i+1,j)
            D2_x_minus = phi(i-2,j) + phi(i,j) - 2*phi(i-1,j)
            
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x)*inv_dx
            else
              phi_x_plus = (D1_x_plus 
     &                          - half*D2_x_plus)*inv_dx
            endif

c           phi_x_minus
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = (D1_x 
     &                           + half*D2_x_minus)*inv_dx
            else
              phi_x_minus = (D1_x
     &                           + half*D2_x)*inv_dx
            endif

c           phi_y_plus
            D1_y = phi(i,j) - phi(i,j-1)
            D1_y_plus = phi(i,j+1) - phi(i,j)
            D2_y = phi(i,j+1) + phi(i,j-1) - 2*phi(i,j)
            D2_y_plus = phi(i,j+2) + phi(i,j) - 2*phi(i,j+1)
            D2_y_minus = phi(i,j-2) + phi(i,j) - 2*phi(i,j-1)
            
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y)*inv_dy
            else
              phi_y_plus = (D1_y_plus 
     &                          - half*D2_y_plus)*inv_dy
            endif

c           phi_y_minus
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = (D1_y 
     &                           + half*D2_y_minus)*inv_dy
            else
              phi_y_minus = (D1_y
     &                           + half*D2_y)*inv_dy
            endif
            
            vel_n_cur = vel_n
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi


c             compute contribution to lse_rhs(i,j) 
              temp_norm = - vel_n_cur*sqrt(norm_grad_phi_sq)

            endif
        
c           get max_H_over_dX for calculating stable dt
c           initialize max_H_over_dX to -1
            max_H_over_dX = -1.0d0

c           compute max_dx_sq
            max_dx_sq = max(dx,dy)
            max_dx_sq = max(dx,dy) * max(dx,dy)
      
            phi_x_cur = max(abs(phi_x_plus),
     &                          abs(phi_x_minus))
            phi_y_cur = max(abs(phi_y_plus),
     &                          abs(phi_y_minus))
            norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                            + phi_y_cur*phi_y_cur 
     &                            + max_dx_sq )

            H_over_dX_cur = abs(vel_n) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

            if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
            endif
        
     
c         compute curvature term
            phi_x = (phi(i+1,j) - phi(i-1,j))*dx_factor
            phi_y = (phi(i,j+1) - phi(i,j-1))*dy_factor
               
            phi_xx = (phi(i+2,j) + phi(i-2,j)
     &               - 2*phi(i,j))*dx_factor*dx_factor
            phi_yy = (phi(i,j+2) + phi(i,j-2)
     &               - 2*phi(i,j))*dy_factor*dy_factor
        
            phi_xy = (phi(i+1,j+1) + phi(i-1,j-1)
     &               - phi(i-1,j+1) - phi(i+1,j-1))
     &               *dx_factor*dy_factor
    
c           compute squared magnitude of gradient
            grad_mag2 = phi_x * phi_x
     &                + phi_y * phi_y
            if (grad_mag2 .lt. zero_tol) then
              curv = 0.d0
            else
              curv = phi_xx*phi_y*phi_y 
     &             +   phi_yy*phi_x*phi_x  
     &             - 2*phi_xy*phi_x*phi_y
              curv = curv / grad_mag2 
            endif

            temp_curv = b*curv

	    lse_rhs(i,j) = temp_norm + temp_curv + temp_vel    
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c  Computes the right hand side of the 2D level set equation, 
c  with third order accurate stencil for  variational a, b and 
c  convective velocity terms
c***********************************************************************
      subroutine qss2dgetRHS3rd(
     &  phi,
     &  vel_n,
     &  b,
     &  vel_x, vel_y,
     &  max_H_over_dX,
     &  lse_rhs,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      
c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus
      real phi_x_minus, phi_y_minus
      
c     upwinded derivatives
      real phi_x_adv, phi_y_adv
    
c     central difference derivatives
      real phi_x, phi_y
      
c     second order central difference derivatives
      real phi_xx, phi_yy, phi_xy
      
      real curv, grad_mag2
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_n(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real b(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real vel_y(ilo_phi_gb:ihi_phi_gb,
     &           jlo_phi_gb:jhi_phi_gb)
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb)
      real D1_x, D1_x_plus, D1_x_plus_plus, D1_x_plus_plus_plus
      real D1_x_minus, D1_x_minus_minus
      real D1_y, D1_y_plus, D1_y_plus_plus, D1_y_plus_plus_plus
      real D1_y_minus, D1_y_minus_minus
       
      real D2_x, D2_x_plus, D2_x_minus, D2_x_minus_minus
      real D2_x_plus_plus
      real D2_y, D2_y_plus, D2_y_minus, D2_y_minus_minus
      real D2_y_plus_plus
      
      real D3_x, D3_x_plus, D3_x_plus_plus, D3_x_minus
      real D3_y, D3_y_plus, D3_y_plus_plus, D3_y_minus
      
      
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      parameter (sixth=1.d0/6.d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real max_H_over_dX
      real H_over_dX_cur
      real max_dx_sq
      real norm_grad_phi, phi_x_cur, phi_y_cur
      real temp_norm, temp_curv, temp_vel

      real dx_factor, dy_factor

c     compute denominator values, for curvature term
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      
c     { begin loop over grid 
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            D1_x = phi(i,j) - phi(i-1,j)
            D1_x_plus = phi(i+1,j) - phi(i,j)
            D1_x_minus = phi(i-1,j) - phi(i-2,j)
            D1_x_minus_minus = phi(i-2,j) - phi(i-3,j)
            D1_x_plus_plus = phi(i+2,j) - phi(i+1,j)
            D1_x_plus_plus_plus = phi(i+3,j) - phi(i+2,j)
            
            D2_x = D1_x_plus - D1_x
            D2_x_plus = D1_x_plus_plus - D1_x_plus
            D2_x_minus = D1_x - D1_x_minus
            D2_x_plus_plus = D1_x_plus_plus_plus - D1_x_plus_plus
            D2_x_minus_minus = D1_x_minus - D1_x_minus_minus
            
            D3_x = D2_x - D2_x_minus
            D3_x_plus = D2_x_plus - D2_x
            D3_x_plus_plus = D2_x_plus_plus - D2_x_plus
            D3_x_minus = D2_x_minus - D2_x_minus_minus       
            
c           { begin calculation of phi_x_plus
            phi_x_plus = D1_x_plus
  
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = phi_x_plus
     &                          - half*D2_x 
              if (abs(D3_x).lt.abs(D3_x_plus)) then
                phi_x_plus = phi_x_plus 
     &                            - sixth*D3_x
              else
                phi_x_plus = phi_x_plus 
     &                            - sixth*D3_x_plus
              endif
            else
              phi_x_plus = phi_x_plus
     &                          - half*D2_x_plus
              if (abs(D3_x_plus).lt.abs(D3_x_plus_plus)) then
                phi_x_plus = phi_x_plus
     &                            + third*D3_x_plus
              else
                phi_x_plus = phi_x_plus 
     &                            + third*D3_x_plus_plus
              endif
            endif

c           divide phi_x_plus by dx
            phi_x_plus = phi_x_plus*inv_dx

c           } end calculation of phi_x_plus

c           { begin calculation of phi_x_minus
            phi_x_minus = D1_x
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = phi_x_minus 
     &                           + half*D2_x_minus 
              if (abs(D3_x_minus).lt.abs(D3_x)) then
                phi_x_minus = phi_x_minus 
     &                             + third*D3_x_minus
              else
                phi_x_minus = phi_x_minus 
     &                             + third*D3_x
              endif
            else
              phi_x_minus = phi_x_minus 
     &                           + half*D2_x 
              if (abs(D3_x).lt.abs(D3_x_plus)) then
                phi_x_minus = phi_x_minus 
     &                            - sixth*D3_x
              else
                phi_x_minus = phi_x_minus 
     &                            - sixth*D3_x_plus
              endif
            endif

c           divide phi_x_minus by dx
            phi_x_minus = phi_x_minus*inv_dx

c           } end calculation of phi_x_minus


            D1_y = phi(i,j) - phi(i,j-1)
            D1_y_plus = phi(i,j+1) - phi(i,j)
            D1_y_minus = phi(i,j-1) - phi(i,j-2)
            D1_y_minus_minus = phi(i,j-2) - phi(i,j-3)
            D1_y_plus_plus = phi(i,j+2) - phi(i,j+1)
            D1_y_plus_plus_plus = phi(i,j+3) - phi(i,j+2)
            
            D2_y = D1_y_plus - D1_y
            D2_y_plus = D1_y_plus_plus - D1_y_plus
            D2_y_minus = D1_y - D1_y_minus
            D2_y_plus_plus = D1_y_plus_plus_plus - D1_y_plus_plus
            D2_y_minus_minus = D1_y_minus - D1_y_minus_minus
            
            D3_y = D2_y - D2_y_minus
            D3_y_plus = D2_y_plus - D2_y
            D3_y_plus_plus = D2_y_plus_plus - D2_y_plus
            D3_y_minus = D2_y_minus - D2_y_minus_minus       
            
c           { begin calculation of phi_y_plus
            phi_y_plus = D1_y_plus
  
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = phi_y_plus
     &                          - half*D2_y 
              if (abs(D3_y).lt.abs(D3_y_plus)) then
                phi_y_plus = phi_y_plus 
     &                            - sixth*D3_y
              else
                phi_y_plus = phi_y_plus 
     &                            - sixth*D3_y_plus
              endif
            else
              phi_y_plus = phi_y_plus
     &                          - half*D2_y_plus
              if (abs(D3_y_plus).lt.abs(D3_y_plus_plus)) then
                phi_y_plus = phi_y_plus
     &                            + third*D3_y_plus
              else
                phi_y_plus = phi_y_plus 
     &                            + third*D3_y_plus_plus
              endif
            endif

c           divide phi_y_plus by dx
            phi_y_plus = phi_y_plus*inv_dy

c           } end calculation of phi_y_plus

c           { begin calculation of phi_y_minus
            phi_y_minus = D1_y
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = phi_y_minus 
     &                           + half*D2_y_minus 
              if (abs(D3_y_minus).lt.abs(D3_y)) then
                phi_y_minus = phi_y_minus 
     &                             + third*D3_y_minus
              else
                phi_y_minus = phi_y_minus 
     &                             + third*D3_y
              endif
            else
              phi_y_minus = phi_y_minus 
     &                           + half*D2_y 
              if (abs(D3_y).lt.abs(D3_y_plus)) then
                phi_y_minus = phi_y_minus 
     &                            - sixth*D3_y
              else
                phi_y_minus = phi_y_minus 
     &                            - sixth*D3_y_plus
              endif
            endif

c           divide phi_y_minus by dx
            phi_y_minus = phi_y_minus*inv_dx

c           } end calculation of phi_y_minus

            vel_n_cur = vel_n(i,j)
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi


c             compute contribution to lse_rhs(i,j,k) 
              temp_norm = - vel_n_cur*sqrt(norm_grad_phi_sq)

            endif
        
c           get max_H_over_dX for calculating stable dt
c           initialize max_H_over_dX to -1
            max_H_over_dX = -1.0d0

c           compute max_dx_sq
            max_dx_sq = max(dx,dy)
            max_dx_sq = max(dx,dy) * max(dx,dy)
      
            phi_x_cur = max(abs(phi_x_plus),
     &                          abs(phi_x_minus))
            phi_y_cur = max(abs(phi_y_plus),
     &                          abs(phi_y_minus))
            norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                            + phi_y_cur*phi_y_cur 
     &                            + max_dx_sq )

            H_over_dX_cur = abs(vel_n(i,j)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

            if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
            endif
        
            
    
c           compute advective term contributions
c              phi_x
              if (vel_x(i,j).ge.0) then
                phi_x_adv = phi_x_minus
              else
                phi_x_adv = phi_x_plus
              endif

c              phi_y
              if (vel_y(i,j).ge.0) then
                phi_y_adv = phi_y_minus
              else
                phi_y_adv = phi_y_plus
              endif    
              
            temp_vel = - ( vel_x(i,j)*phi_x_adv
     &                       + vel_y(i,j)*phi_y_adv )
     
     
c         compute curvature term
            phi_x = (phi(i+1,j) - phi(i-1,j))*dx_factor
            phi_y = (phi(i,j+1) - phi(i,j-1))*dy_factor
   
            phi_xx = (phi(i+2,j) + phi(i-2,j)
     &               - 2*phi(i,j))*dx_factor*dx_factor

            phi_yy = (phi(i,j+2) + phi(i,j-2)
     &               - 2*phi(i,j))*dy_factor*dy_factor
        
            phi_xy = (phi(i+1,j+1) + phi(i-1,j-1)
     &               - phi(i-1,j+1) - phi(i+1,j-1))
     &               *dx_factor*dy_factor
    
c           compute squared magnitude of gradient
            grad_mag2 = phi_x * phi_x
     &                + phi_y * phi_y
            if (grad_mag2 .lt. zero_tol) then
              curv = 0.d0
            else
              curv = phi_xx*phi_y*phi_y 
     &             +   phi_yy*phi_x*phi_x  
     &             - 2*phi_xy*phi_x*phi_y
              curv = curv / grad_mag2 
            endif

            temp_curv = b(i,j)*curv

	    lse_rhs(i,j) = temp_norm + temp_curv + temp_vel    
          enddo
        enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c  Compute the right hand side of the 2D level set equation, with 
c  third order accurate stencils, and constant zero contact angle 
c  everywhere - so a and b are not varying
c  and the convective velocity term does not exist.
c***********************************************************************
      subroutine qss2dgetRHS3rdConst(
     &  phi,
     &  vel_n,
     &  b,
     &  max_H_over_dX,
     &  lse_rhs,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none
      
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      
c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus
      real phi_x_minus, phi_y_minus
      
c     upwinded derivatives
      real phi_x_adv, phi_y_adv
    
c     central difference derivatives
      real phi_x, phi_y
      
c     second order central difference derivatives
      real phi_xx, phi_yy, phi_xy
      
      real curv, grad_mag2
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_n
      real b
      real lse_rhs(ilo_phi_gb:ihi_phi_gb,
     &             jlo_phi_gb:jhi_phi_gb)
      real D1_x, D1_x_plus, D1_x_plus_plus, D1_x_plus_plus_plus
      real D1_x_minus, D1_x_minus_minus
      real D1_y, D1_y_plus, D1_y_plus_plus, D1_y_plus_plus_plus
      real D1_y_minus, D1_y_minus_minus
       
      real D2_x, D2_x_plus, D2_x_minus, D2_x_minus_minus
      real D2_x_plus_plus
      real D2_y, D2_y_plus, D2_y_minus, D2_y_minus_minus
      real D2_y_plus_plus
      
      real D3_x, D3_x_plus, D3_x_plus_plus, D3_x_minus
      real D3_y, D3_y_plus, D3_y_plus_plus, D3_y_minus
      
      
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      parameter (sixth=1.d0/6.d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)
      real vel_n_cur
      real norm_grad_phi_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real max_H_over_dX
      real H_over_dX_cur
      real max_dx_sq
      real norm_grad_phi, phi_x_cur, phi_y_cur
      real temp_norm, temp_curv, temp_vel

      real dx_factor, dy_factor

c     compute denominator values, for curvature term
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      
c     { begin loop over grid 
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            D1_x = phi(i,j) - phi(i-1,j)
            D1_x_plus = phi(i+1,j) - phi(i,j)
            D1_x_minus = phi(i-1,j) - phi(i-2,j)
            D1_x_minus_minus = phi(i-2,j) - phi(i-3,j)
            D1_x_plus_plus = phi(i+2,j) - phi(i+1,j)
            D1_x_plus_plus_plus = phi(i+3,j) - phi(i+2,j)
            
            D2_x = D1_x_plus - D1_x
            D2_x_plus = D1_x_plus_plus - D1_x_plus
            D2_x_minus = D1_x - D1_x_minus
            D2_x_plus_plus = D1_x_plus_plus_plus - D1_x_plus_plus
            D2_x_minus_minus = D1_x_minus - D1_x_minus_minus
            
            D3_x = D2_x - D2_x_minus
            D3_x_plus = D2_x_plus - D2_x
            D3_x_plus_plus = D2_x_plus_plus - D2_x_plus
            D3_x_minus = D2_x_minus - D2_x_minus_minus       
            
c           { begin calculation of phi_x_plus
            phi_x_plus = D1_x_plus
  
            if (abs(D2_x).lt.abs(D2_x_plus)) then
              phi_x_plus = phi_x_plus
     &                          - half*D2_x 
              if (abs(D3_x).lt.abs(D3_x_plus)) then
                phi_x_plus = phi_x_plus 
     &                            - sixth*D3_x
              else
                phi_x_plus = phi_x_plus 
     &                            - sixth*D3_x_plus
              endif
            else
              phi_x_plus = phi_x_plus
     &                          - half*D2_x_plus
              if (abs(D3_x_plus).lt.abs(D3_x_plus_plus)) then
                phi_x_plus = phi_x_plus
     &                            + third*D3_x_plus
              else
                phi_x_plus = phi_x_plus 
     &                            + third*D3_x_plus_plus
              endif
            endif

c           divide phi_x_plus by dx
            phi_x_plus = phi_x_plus*inv_dx

c           } end calculation of phi_x_plus

c           { begin calculation of phi_x_minus
            phi_x_minus = D1_x
            if (abs(D2_x_minus).lt.abs(D2_x)) then
              phi_x_minus = phi_x_minus 
     &                           + half*D2_x_minus 
              if (abs(D3_x_minus).lt.abs(D3_x)) then
                phi_x_minus = phi_x_minus 
     &                             + third*D3_x_minus
              else
                phi_x_minus = phi_x_minus 
     &                             + third*D3_x
              endif
            else
              phi_x_minus = phi_x_minus 
     &                           + half*D2_x 
              if (abs(D3_x).lt.abs(D3_x_plus)) then
                phi_x_minus = phi_x_minus 
     &                            - sixth*D3_x
              else
                phi_x_minus = phi_x_minus 
     &                            - sixth*D3_x_plus
              endif
            endif

c           divide phi_x_minus by dx
            phi_x_minus = phi_x_minus*inv_dx

c           } end calculation of phi_x_minus


            D1_y = phi(i,j) - phi(i,j-1)
            D1_y_plus = phi(i,j+1) - phi(i,j)
            D1_y_minus = phi(i,j-1) - phi(i,j-2)
            D1_y_minus_minus = phi(i,j-2) - phi(i,j-3)
            D1_y_plus_plus = phi(i,j+2) - phi(i,j+1)
            D1_y_plus_plus_plus = phi(i,j+3) - phi(i,j+2)
            
            D2_y = D1_y_plus - D1_y
            D2_y_plus = D1_y_plus_plus - D1_y_plus
            D2_y_minus = D1_y - D1_y_minus
            D2_y_plus_plus = D1_y_plus_plus_plus - D1_y_plus_plus
            D2_y_minus_minus = D1_y_minus - D1_y_minus_minus
            
            D3_y = D2_y - D2_y_minus
            D3_y_plus = D2_y_plus - D2_y
            D3_y_plus_plus = D2_y_plus_plus - D2_y_plus
            D3_y_minus = D2_y_minus - D2_y_minus_minus       
            
c           { begin calculation of phi_y_plus
            phi_y_plus = D1_y_plus
  
            if (abs(D2_y).lt.abs(D2_y_plus)) then
              phi_y_plus = phi_y_plus
     &                          - half*D2_y 
              if (abs(D3_y).lt.abs(D3_y_plus)) then
                phi_y_plus = phi_y_plus 
     &                            - sixth*D3_y
              else
                phi_y_plus = phi_y_plus 
     &                            - sixth*D3_y_plus
              endif
            else
              phi_y_plus = phi_y_plus
     &                          - half*D2_y_plus
              if (abs(D3_y_plus).lt.abs(D3_y_plus_plus)) then
                phi_y_plus = phi_y_plus
     &                            + third*D3_y_plus
              else
                phi_y_plus = phi_y_plus 
     &                            + third*D3_y_plus_plus
              endif
            endif

c           divide phi_y_plus by dx
            phi_y_plus = phi_y_plus*inv_dy

c           } end calculation of phi_y_plus

c           { begin calculation of phi_y_minus
            phi_y_minus = D1_y
            if (abs(D2_y_minus).lt.abs(D2_y)) then
              phi_y_minus = phi_y_minus 
     &                           + half*D2_y_minus 
              if (abs(D3_y_minus).lt.abs(D3_y)) then
                phi_y_minus = phi_y_minus 
     &                             + third*D3_y_minus
              else
                phi_y_minus = phi_y_minus 
     &                             + third*D3_y
              endif
            else
              phi_y_minus = phi_y_minus 
     &                           + half*D2_y 
              if (abs(D3_y).lt.abs(D3_y_plus)) then
                phi_y_minus = phi_y_minus 
     &                            - sixth*D3_y
              else
                phi_y_minus = phi_y_minus 
     &                            - sixth*D3_y_plus
              endif
            endif

c           divide phi_y_minus by dx
            phi_y_minus = phi_y_minus*inv_dx

c           } end calculation of phi_y_minus

            vel_n_cur = vel_n
            if (abs(vel_n_cur) .ge. zero_tol) then

c             { begin Godunov selection of grad_phi

              if (vel_n_cur .gt. 0.d0) then
                norm_grad_phi_sq = max(max(phi_x_minus,0.d0)**2,
     &                                 min(phi_x_plus,0.d0)**2 )
     &                           + max(max(phi_y_minus,0.d0)**2,
     &                                 min(phi_y_plus,0.d0)**2 )
              else
                norm_grad_phi_sq = max(min(phi_x_minus,0.d0)**2,
     &                                 max(phi_x_plus,0.d0)**2 )
     &                           + max(min(phi_y_minus,0.d0)**2,
     &                                 max(phi_y_plus,0.d0)**2 )
              endif

c             } end Godunov selection of grad_phi


c             compute contribution to lse_rhs(i,j,k) 
              temp_norm = - vel_n_cur*sqrt(norm_grad_phi_sq)

            endif
        
c           get max_H_over_dX for calculating stable dt
c           initialize max_H_over_dX to -1
            max_H_over_dX = -1.0d0

c           compute max_dx_sq
            max_dx_sq = max(dx,dy)
            max_dx_sq = max(dx,dy) * max(dx,dy)
      
            phi_x_cur = max(abs(phi_x_plus),
     &                          abs(phi_x_minus))
            phi_y_cur = max(abs(phi_y_plus),
     &                          abs(phi_y_minus))
            norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                            + phi_y_cur*phi_y_cur 
     &                            + max_dx_sq )

            H_over_dX_cur = abs(vel_n) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

            if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
            endif
        
     
c         compute curvature term
            phi_x = (phi(i+1,j) - phi(i-1,j))*dx_factor
            phi_y = (phi(i,j+1) - phi(i,j-1))*dy_factor
   
            phi_xx = (phi(i+2,j) + phi(i-2,j)
     &               - 2*phi(i,j))*dx_factor*dx_factor

            phi_yy = (phi(i,j+2) + phi(i,j-2)
     &               - 2*phi(i,j))*dy_factor*dy_factor
        
            phi_xy = (phi(i+1,j+1) + phi(i-1,j-1)
     &               - phi(i-1,j+1) - phi(i+1,j-1))
     &               *dx_factor*dy_factor
    
c           compute squared magnitude of gradient
            grad_mag2 = phi_x * phi_x
     &                + phi_y * phi_y
            if (grad_mag2 .lt. zero_tol) then
              curv = 0.d0
            else
              curv = phi_xx*phi_y*phi_y 
     &             +   phi_yy*phi_x*phi_x  
     &             - 2*phi_xy*phi_x*phi_y
              curv = curv / grad_mag2 
            endif

            temp_curv = b*curv

	    lse_rhs(i,j) = temp_norm + temp_curv 
          enddo
        enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c  Computes the 2D central differencing gradient. Applied to curvature
c  term calculation
c***********************************************************************
      subroutine qss2dCentralGradOrder2(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer i,j
      real dx_factor, dy_factor
c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     { begin loop over grid 
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            phi_x(i,j) = (phi(i+1,j) - phi(i-1,j))*dx_factor
            phi_y(i,j) = (phi(i,j+1) - phi(i,j-1))*dy_factor
   
          enddo
        enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c  Computes the 2D signed linear extrapolation boundary conditions
c***********************************************************************
      subroutine qss2dSignedLinearExtrapolation(
     &  phi,
     &  ilo_gb, ihi_gb, jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  bdry_location_idx)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb, jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      integer bdry_location_idx
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables       
      integer i,j, zero
      parameter (zero = 0)
      real s, abs_diff, dist, slope
      real one
      parameter (one = 1.0d0)

      if (bdry_location_idx .eq. 0) then
c     { extrapolate data in x-direction at lower end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          s = sign(one,phi(ilo_fb,j))
          abs_diff = abs(phi(ilo_fb,j) - phi(ilo_fb+1,j))
          slope = s*abs_diff
          do i = ilo_gb, ilo_fb-1
            dist = ilo_fb - i
            phi(i,j) = phi(ilo_fb,j) + slope*dist
          enddo
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at lower end
       
      elseif (bdry_location_idx .eq. 1) then
c     { extrapolate data in x-direction at upper end

c       { begin j loop
        do j = jlo_gb, jhi_gb
          s = sign(one,phi(ihi_fb,j))
          abs_diff = abs(phi(ihi_fb,j) - phi(ihi_fb-1,j))
          slope = s*abs_diff
          do i = ihi_fb+1, ihi_gb
            dist = i - ihi_fb
            phi(i,j) = phi(ihi_fb,j) + slope*dist
          enddo 
        enddo  
c       } end j loop

c     } end extrapolate data in x-direction at upper end

      elseif (bdry_location_idx .eq. 2) then
c     { extrapolate data in y-direction at lower end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          s = sign(one,phi(i,jlo_fb))
          abs_diff = abs(phi(i,jlo_fb) - phi(i,jlo_fb+1))
          slope = s*abs_diff
          do j = jlo_gb, jlo_fb-1
            dist = jlo_fb - j 
            phi(i,j) = phi(i,jlo_fb) + slope*dist
          enddo
        enddo 
c       } end i loop

c     } end extrapolate data in y-direction at lower end

      elseif (bdry_location_idx .eq. 3) then
c     { extrapolate data in y-direction at upper end

c       { begin i loop
        do i = ilo_gb, ihi_gb
          s = sign(one,phi(i,jhi_fb))
          abs_diff = abs(phi(i,jhi_fb) - phi(i,jhi_fb-1))
          slope = s*abs_diff
          do j = jhi_fb+1, jhi_gb
            dist = j - jhi_fb
            phi(i,j) = phi(i,jhi_fb) + slope*dist
          enddo 
        enddo 
c       } end i loop

c     } end extrapolate data in y-direction at lower end

      endif

      return
      end
c } end subroutine
c***********************************************************************



