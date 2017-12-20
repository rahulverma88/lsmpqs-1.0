c***********************************************************************
c
c  File:        qss_reinitialization2d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision: 149 $
c  Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
c  Description: F77 routines for reinitialization of 2d level set functions
c
c***********************************************************************

c***********************************************************************
c The algorithms and notation in these subroutines closely follows
c the discussion in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  qss2dComputeReinitializationEqnRHS() computes the right-hand side of 
c  the reinitialization equation using a Godunov scheme to select the 
c  numerical discretization of the sgn(phi) |grad(phi)| term.  
c  Forward (plus) and backward (minus) spatial derivatives used in 
c  the Godunov calculation must be supplied by the user.
c
c  Arguments:
c    reinit_rhs (out):       right-hand side of reinitialization 
c                            equation
c    phi (in):               level set function at current iteration
c                            of reinitialization process
c    phi0 (in):              level set function at initial iteration
c                            iteration of reinitialization process
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    use_phi0_for_sgn (in):  flag to specify whether phi0 should be
c                            used in the computation of sgn(phi).
c                              0 = use phi (do NOT use phi0)
c                              1 = use phi0
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c  NOTES:
c   (1) if use_phi0_for_sgn is not equal to 0 or 1, the default
c       behavior of qss2dComputeReinitializationEqnRHS() is to use 
c       phi (i.e. equivalent to setting use_phi0_for_sgn to 0)
c
c***********************************************************************
      subroutine ComputeReinitializationEqnRHS2d(
     &  reinit_rhs,
     &  phi,
     &  phi0,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy,
     &  use_phi0_for_sgn)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real reinit_rhs(ilo_gb:ihi_gb,
     &                jlo_gb:jhi_gb)
      real phi(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
      real phi0(ilo_gb:ihi_gb,
     &          jlo_gb:jhi_gb)
      real dx, dy
      integer use_phi0_for_sgn
      real phi_cur
      integer DIM
      parameter (DIM=2)
      real grad_phi_plus_cur(1:DIM)
      real grad_phi_minus_cur(1:DIM)
      real grad_phi_star(1:DIM)
      integer i,j
      integer dir
      real sgn_phi
      real norm_grad_phi_sq
      real dx_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real one
      parameter (one=1.d0)

c     plus and minus upwind derivatives
      real phi_x_plus, phi_y_plus
      real phi_x_minus, phi_y_minus
      
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
      
      real inv_dx, inv_dy
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      
c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy



c     set value of dx_sq to be square of max{dx,dy}
      dx_sq = max(dx,dy)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method
c----------------------------------------------------
c     { begin condition on use_phi0_for_sgn

c       -----------------------------------------------
c       use phi in computation of smoothed sgn function
c       -----------------------------------------------
c       { begin loop over grid
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
            
            if (use_phi0_for_sgn .ne. 1) then

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

c             cache phi and spatial derivative approximations
              phi_cur = phi(i,j)
              grad_phi_plus_cur(1) = phi_x_plus
              grad_phi_plus_cur(2) = phi_y_plus
              grad_phi_minus_cur(1) = phi_x_minus
              grad_phi_minus_cur(2) = phi_y_minus
  
c             { begin Godunov selection of grad_phi
              do dir=1,DIM

                if (phi_cur .gt. 0.d0) then
                  grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                  grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                          0.d0)
                else
                  grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                  grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
                endif

                grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                   grad_phi_minus_cur(dir)) 

              enddo
c             } end Godunov selection of grad_phi

c             compute reinit_rhs(i,j) using smoothed sgn(phi)
              if (abs(phi_cur) .ge. zero_tol) then
                norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                           + grad_phi_star(2)*grad_phi_star(2)
                sgn_phi = phi_cur
     &                  / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)
                reinit_rhs(i,j) = 
     &            sgn_phi*(one - sqrt(norm_grad_phi_sq))
              else
                reinit_rhs(i,j) = 0.d0
              endif
        else

c       ------------------------------------------------
c       use phi0 in computation of smoothed sgn function 
c       ------------------------------------------------
  
c             cache phi and spatial derivative approximations
              phi_cur = phi0(i,j)
              grad_phi_plus_cur(1) = phi_x_plus
              grad_phi_plus_cur(2) = phi_y_plus
              grad_phi_minus_cur(1) = phi_x_minus
              grad_phi_minus_cur(2) = phi_y_minus
  
c             { begin Godunov selection of grad_phi
              do dir=1,DIM

                if (phi_cur .gt. 0.d0) then
                  grad_phi_plus_cur(dir) = max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                  grad_phi_minus_cur(dir) = max(grad_phi_minus_cur(dir),
     &                                          0.d0)
                else
                  grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                  grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
                endif

                grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                   grad_phi_minus_cur(dir)) 

              enddo
c             } end Godunov selection of grad_phi

c             compute reinit_rhs(i,j) using smoothed sgn(phi)
              if (abs(phi_cur) .ge. zero_tol) then
                norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                           + grad_phi_star(2)*grad_phi_star(2)
                sgn_phi = phi_cur / sqrt(phi_cur*phi_cur + dx_sq)
                reinit_rhs(i,j) = 
     &            sgn_phi*(one - sqrt(norm_grad_phi_sq))
              else
                reinit_rhs(i,j) = 0.d0
              endif

             endif
c     } end condition on use_phi0_for_sgn

            enddo
          enddo
c       } end loop over grid


      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  qss2dComputeOrthogonalizationEqnRHS() computes the right-hand side 
c  of the orthogonalization equation:
c
c    phi_t + grad(phi) dot { sgn(psi)/|grad(psi)| grad(psi) } = 0
c
c  Upwinding is used to select whether the forward (plus) or backward
c  (minus) spatial derivative should be used for grad(phi).  grad(psi)
c  is computed by averaging the forward and backward spatial derivatives
c  for grad(psi).  Forward and backward spatial derivatives used in the 
c  calculation must be supplied by the user.
c
c  Arguments:
c    othro_rhs (out):        right-hand side of orthogonalization
c                            equation
c    psi (in):               data array for psi
c    phi_*_plus (in):        forward spatial derivatives for grad(phi)
c    phi_*_minus (in):       backward spatial derivatives for grad(phi)
c    psi_*_plus (in):        forward spatial derivatives for grad(psi)
c    psi_*_minus (in):       backward spatial derivatives for grad(psi)
c    *_gb (in):              index range for ghostbox
c    *_fb (in):              index range for fillbox
c
c***********************************************************************
      subroutine qss2dComputeOrthogonalizationEqnRHS(
     &  ortho_rhs,
     &  ilo_rhs_gb, ihi_rhs_gb, 
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  psi,
     &  ilo_psi_gb, ihi_psi_gb, 
     &  jlo_psi_gb, jhi_psi_gb,
     &  psi_x_plus, psi_y_plus,
     &  ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb, 
     &  jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb,
     &  psi_x_minus, psi_y_minus,
     &  ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb, 
     &  jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, 
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_psi_gb, ihi_psi_gb
      integer jlo_psi_gb, jhi_psi_gb
      integer ilo_grad_psi_plus_gb, ihi_grad_psi_plus_gb
      integer jlo_grad_psi_plus_gb, jhi_grad_psi_plus_gb
      integer ilo_grad_psi_minus_gb, ihi_grad_psi_minus_gb
      integer jlo_grad_psi_minus_gb, jhi_grad_psi_minus_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real ortho_rhs(ilo_rhs_gb:ihi_rhs_gb,
     &               jlo_rhs_gb:jhi_rhs_gb)
      real psi(ilo_psi_gb:ihi_psi_gb,
     &         jlo_psi_gb:jhi_psi_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real psi_x_plus(ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb)
      real psi_y_plus(ilo_grad_psi_plus_gb:ihi_grad_psi_plus_gb,
     &                jlo_grad_psi_plus_gb:jhi_grad_psi_plus_gb)
      real psi_x_minus(ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                 jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb)
      real psi_y_minus(ilo_grad_psi_minus_gb:ihi_grad_psi_minus_gb,
     &                 jlo_grad_psi_minus_gb:jhi_grad_psi_minus_gb)
      real dx, dy
      real dx_sq
      integer DIM
      parameter (DIM=2)
      real grad_psi_star(1:DIM)
      real norm_grad_psi
      real sgn_psi
      integer i,j
      real psi_zero_tol, grad_psi_zero_tol
      parameter (psi_zero_tol=1.e-5)
      parameter (grad_psi_zero_tol=1.e-5)
      real one, half
      parameter (one=1.d0,half=0.5d0)

c     set value of dx_sq to be square of max{dx,dy}
      dx_sq = max(dx,dy)
      dx_sq = dx_sq*dx_sq

c----------------------------------------------------
c      compute RHS of orthogonalization equation
c      using upwinding
c----------------------------------------------------

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         average forward and backward derivatives of 
c         grad(psi) to get grad_psi_star
          grad_psi_star(1) = half 
     &                     * (psi_x_plus(i,j)+ psi_x_minus(i,j))
          grad_psi_star(2) = half 
     &                     * (psi_y_plus(i,j)+ psi_y_minus(i,j))
 
c         compute norm of grad(psi)
          norm_grad_psi = grad_psi_star(1)*grad_psi_star(1)
     &                  + grad_psi_star(2)*grad_psi_star(2)
          norm_grad_psi = sqrt(norm_grad_psi)

c         compute ortho_rhs(i,j) using upwinding on sgn_psi*grad(psi) 
          if ( (abs(psi(i,j)) .gt. psi_zero_tol) .and.
     &         (norm_grad_psi .gt. grad_psi_zero_tol) ) then

c           CASE: nontrivial psi and grad(psi) 

            sgn_psi = psi(i,j)/sqrt(psi(i,j)*psi(i,j) + dx_sq)

            ortho_rhs(i,j) = -1.d0/norm_grad_psi 
     &        * ( max(sgn_psi*grad_psi_star(1),0.d0)
     &                       *phi_x_minus(i,j)
     &          + min(sgn_psi*grad_psi_star(1),0.d0)
     &                        *phi_x_plus(i,j)
     &          + max(sgn_psi*grad_psi_star(2),0.d0)
     &                        *phi_y_minus(i,j)
     &          + min(sgn_psi*grad_psi_star(2),0.d0)
     &                        *phi_y_plus(i,j) )
          else

c           CASE: grad(psi) = 0 CASE: psi = 0 

            ortho_rhs(i,j) = 0.d0

          endif

        enddo
      enddo
c     } end loop over grid


      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************

      subroutine qss2dComputeDistanceForSubcellFix(
     &  distance0,
     &  phi0,
     &  ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_phi0_gb, ihi_phi0_gb, jlo_phi0_gb, jhi_phi0_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real distance0(ilo_phi0_gb:ihi_phi0_gb,
     &               jlo_phi0_gb:jhi_phi0_gb)
      real phi0(ilo_phi0_gb:ihi_phi0_gb,
     &          jlo_phi0_gb:jhi_phi0_gb)
      real dx, dy, max_dx
      real zero_tol
      parameter (zero_tol=1.e-5)
      real large_distance_flag,  zero
      parameter (zero=0.d0)
      real near_plus_x, near_minus_x
      real near_plus_y, near_minus_y
      logical near
      real d1, d2, d3, delta
      integer i,j
      
      max_dx = max(dx,dy);
      large_distance_flag = -1000.d0*max_dx;
      
c       { begin loop over grid
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           determine whether we're near the interface
c           dumb/quick way to do it - this should be
c           computed once and provided to the function
            near_plus_x = phi0(i+1,j)*phi0(i,j);
            near_minus_x = phi0(i-1,j)*phi0(i,j);
	    near_plus_y = phi0(i,j+1)*phi0(i,j);
            near_minus_y = phi0(i,j-1)*phi0(i,j);

            near = (near_plus_x  .lt. zero) .or. 
     &             (near_minus_x .lt. zero) .or.
     &             (near_plus_y  .lt. zero) .or. 
     &             (near_minus_y .lt. zero);
             
	   if( near .eqv. .false. ) then
	     distance0(i,j) = large_distance_flag;
	   else        
              d1 = abs(phi0(i+1,j) - phi0(i-1,j));
	      d2 = abs(phi0(i+1,j) - phi0(i,j));
	      d3 = abs(phi0(i,j) - phi0(i-1,j));
	      delta = max(d1,d2,d3,zero_tol);
	      
	      d1 = abs(phi0(i,j+1) - phi0(i,j-1));
	      d2 = abs(phi0(i,j+1) - phi0(i,j));
	      d3 = abs(phi0(i,j) - phi0(i,j-1));
	      delta = max(d1,d2,d3,delta);
	                
              distance0(i,j) = 2*phi0(i,j)/delta;     
	   endif
	   
	 enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine qss2dComputeReinitializationEqnRHSSubcellFixOrder1(
     &  reinit_rhs,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  phi0, distance0,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real reinit_rhs(ilo_phi_gb:ihi_phi_gb,
     &                jlo_phi_gb:jhi_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi0(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real distance0(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_x_plus
      real phi_y_plus
      real phi_x_minus
      real phi_y_minus
      real dx, dy, max_dx
      real phi_cur
      integer DIM
      parameter (DIM=2)
      real grad_phi_plus_cur(1:DIM)
      real grad_phi_minus_cur(1:DIM)
      real grad_phi_star(1:DIM)
      integer i,j
      integer dir
      real sgn_phi0
      real norm_grad_phi_sq
      real dx_sq
      real zero_tol
      parameter (zero_tol=1.e-5)
      real one, zero, large_distance_flag
      parameter (one=1.d0,  zero=0.d0)

      real D1_x, D1_x_plus, D2_x, D2_x_plus, D2_x_minus
      real D1_y, D1_y_plus, D2_y, D2_y_plus, D2_y_minus
      
      real inv_dx, inv_dy
      
c     set value of dx_sq to be square of max{dx,dy}
      max_dx = max(dx,dy)
      large_distance_flag = -1000.d0*max_dx;
      
c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      
      dx_sq = max_dx*max_dx

c----------------------------------------------------
c      compute RHS of reinitialization equation
c      using Godunov's method with subcell fix near interface
c----------------------------------------------------

c       { begin loop over grid
       do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c           phi_x_plus, minus
            D1_x = phi(i,j) - phi(i-1,j)
            D1_x_plus = phi(i+1,j) - phi(i,j)
            phi_x_plus = (D1_x_plus)*inv_dx
            phi_x_minus = (D1_x)*inv_dx            
            
c           phi_y_plus, minus
            D1_y = phi(i,j) - phi(i,j-1)
            D1_y_plus = phi(i,j+1) - phi(i,j)
            phi_y_plus = (D1_y_plus)*inv_dy
            phi_y_minus = (D1_y)*inv_dy
            
                       
	      if (abs(phi0(i,j)) .gt. zero_tol) then
	        sgn_phi0 = sign(one,phi0(i,j));
	      else
	        sgn_phi0 = zero;
	      endif
     
          if( distance0(i,j) .eq. large_distance_flag ) then   
c           Away from interface we do standard Godunov gradient selection
c             cache phi and spatial derivative approximations
             phi_cur = phi(i,j)
             grad_phi_plus_cur(1) = phi_x_plus
             grad_phi_plus_cur(2) = phi_y_plus
             grad_phi_minus_cur(1) = phi_x_minus
             grad_phi_minus_cur(2) = phi_y_minus

c             { begin Godunov selection of grad_phi
             do dir=1,DIM
               if (phi_cur .gt. 0.d0) then
                grad_phi_plus_cur(dir)=max(-grad_phi_plus_cur(dir),
     &                                         0.d0)
                grad_phi_minus_cur(dir)=max(grad_phi_minus_cur(dir),
     &                                          0.d0)
               else
                 grad_phi_plus_cur(dir) = 
     &              max(grad_phi_plus_cur(dir), 0.d0)
                 grad_phi_minus_cur(dir) = 
     &              max(-grad_phi_minus_cur(dir), 0.d0)
               endif

               grad_phi_star(dir) = max(grad_phi_plus_cur(dir),
     &                                  grad_phi_minus_cur(dir)) 

             enddo
c             } end Godunov selection of grad_phi    

c           compute reinit_rhs(i,j,k) - no need for smoothing since
c           we're away from the interface
            if (abs(phi_cur) .gt. zero_tol) then 
               norm_grad_phi_sq = grad_phi_star(1)*grad_phi_star(1)
     &                          + grad_phi_star(2)*grad_phi_star(2)
               
	           reinit_rhs(i,j) = sgn_phi0 - sgn_phi0*
     &                  	            sqrt(norm_grad_phi_sq)
            else
               reinit_rhs(i,j) = 0.d0
            endif
          else

c        Close to the interface make sure you don't use information from
c        the other side of the interface (Russo/Smereka 2000. subcell fix)
	      reinit_rhs(i,j) = distance0(i,j) - 
     &                       sgn_phi0*abs(phi(i,j))/max_dx;
	      endif
	      
        enddo
       enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
