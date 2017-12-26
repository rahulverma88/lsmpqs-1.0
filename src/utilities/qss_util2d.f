c***********************************************************************
      subroutine qss2dMaxNormDiff(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb)
      real max_norm_diff
      real next_diff
      integer i,j


c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib) 
     &                   - field2(ilo_ib,jlo_ib))

c       loop over included cells { 
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                next_diff = abs(field1(i,j) - field2(i,j))
                if (next_diff .gt. max_norm_diff) then
                  max_norm_diff = next_diff
                endif
 
            enddo
          enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c***********************************************************************

      subroutine qss2dMaxNormDiffLocal(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  local_zone)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb)
      real max_norm_diff
      real next_diff
      real local_zone
      integer i,j


c     initialize max_norm_diff
      max_norm_diff = -1d0
c     abs( field1(ilo_ib,jlo_ib) 
c     &                   - field2(ilo_ib,jlo_ib))

c       loop over included cells { 
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
                if (abs(field1(i,j)) < local_zone) then
                  next_diff = abs(field1(i,j) - field2(i,j))
                  if (next_diff .gt. max_norm_diff) then
                    max_norm_diff = next_diff
                  endif
                endif
            enddo
          enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine imposeMaskParallel2d(
     &  phi_masked,
     &  mask,
     &  phi,
     &  overlap,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  cur_jlo_gb, cur_jhi_gb)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer cur_jlo_gb, cur_jhi_gb
      
      real overlap
      real phi_masked(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
      real mask(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
      real phi(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
     
      integer i, j
      
        do j=cur_jlo_gb,cur_jhi_gb
          do i=ilo_gb,ihi_gb
            
            if ((mask(i,j) - overlap) .gt. phi(i,j)) then
               phi_masked(i,j) = mask(i,j) - overlap
            else
               phi_masked(i,j) = phi(i,j)
            endif
          enddo
        enddo       
      return
      end
      
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine copyDataParallel2d(
     &  phi_copy,
     &  phi_orig,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  cur_jlo_gb, cur_jhi_gb)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer cur_jlo_gb, cur_jhi_gb
      
      real overlap
      real phi_copy(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
      real phi_orig(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb)
     
      integer i, j, k
      
        do j=cur_jlo_gb,cur_jhi_gb
          do i=ilo_gb,ihi_gb
                phi_copy(i,j) = phi_orig(i,j)
          enddo
        enddo
       
      return
      end
      
c } end subroutine
c***********************************************************************
c***********************************************************************
c
c  lsm3dVolumeRegionPhiLessThanZero() computes the volume of the 
c  region where the level set function is less than 0.  
c
c  Arguments:
c    volume (out):          volume of the region where phi < 0
c    phi (in):              level set function
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine qss2dVolumeRegionPhiLessThanZero(
     &  volume,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy

c     initialize volume to zero
      volume = 0.0d0
           
c       loop over included cells {
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
                phi_cur = phi(i,j)
                phi_cur_over_epsilon = phi_cur/epsilon
    
                if (phi_cur .lt. -epsilon) then
                  volume = volume + dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 0.5d0*(1 - phi_cur_over_epsilon
     &                                   - one_over_pi
     &                                   * sin(pi*phi_cur_over_epsilon))
                  volume = volume + one_minus_H*dV
                endif

           enddo
          enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dVolumeRegionPhiGreaterThanZero() computes the volume of the 
c  region where the level set function is greater than 0.  
c
c  Arguments:
c    volume (out):          volume of the region where phi > 0
c    phi (in):              level set function
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine qss2dVolumeRegionPhiGreaterThanZero(
     &  volume,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy
      dV = dx * dy 

c     initialize volume to zero
      volume = 0.0d0

c       loop over included cells {
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  volume = volume + dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5*(1 + phi_cur_over_epsilon 
     &                       + one_over_pi*sin(pi*phi_cur_over_epsilon))
                  volume = volume + H*dV
                endif

            enddo
          enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************

      subroutine imposeTrapVelocity2d( 
     &  phi_w,
     &  phi_nw, 
     &  normal_velocity, 
     &  curvature_coeff,
     &  vel_x, 
     &  vel_y,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine      
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real phi_w(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb) 
      real phi_nw(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb) 
      real normal_velocity(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)    
      real curvature_coeff(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real vel_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
      real vel_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb)
     
      real dx
      integer i, j, k
      real h1, h2, eps_var, x1, x2, x
      eps_var = 1.5 * dx    
      
c       loop over included cells {
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
                x1 = -phi_nw(i,j)
                x2 = -phi_w(i,j) 
                call qssHeaviside(h1, x1, eps_var)
                call qssHeaviside(h2, x2, eps_var)
                x = normal_velocity(i,j)
                normal_velocity(i,j) = h1 * h2 * x
c                if (phi_w(i,j) > eps_var) then
c                    normal_velocity(i,j) = 0
c                    curvature_coeff(i,j) = 0
c                    vel_x(i,j) = 0
c                    vel_y(i,j) = 0
c                endif
                
                x = curvature_coeff(i,j)
                curvature_coeff(i,j) = h1 * h2 * x
                
                x = vel_x(i,j)
                vel_x(i,j) = h1 * h2 * x
                
                x = vel_y(i,j)
                vel_y(i,j) = h1 * h2 * x
                
            enddo
          enddo
c       } end loop over grid
        return
      end
c***********************************************************************

        
     
     
