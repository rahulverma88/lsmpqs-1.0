c***********************************************************************
      subroutine qss3dMaxNormDiff(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  klo_field1_gb, khi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  klo_field2_gb, khi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer klo_field1_gb, khi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer klo_field2_gb, khi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb,
     &            klo_field1_gb:khi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb,
     &            klo_field2_gb:khi_field2_gb)
      real max_norm_diff
      real next_diff
      integer i,j,k


c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib,klo_ib) 
     &                   - field2(ilo_ib,jlo_ib,klo_ib))

c       loop over included cells { 
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                next_diff = abs(field1(i,j,k) - field2(i,j,k))
                if (next_diff .gt. max_norm_diff) then
                  max_norm_diff = next_diff
                endif
 
            enddo
          enddo
        enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c***********************************************************************

      subroutine qss3dMaxNormDiffLocal(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  klo_field1_gb, khi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  klo_field2_gb, khi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  local_zone)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer klo_field1_gb, khi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer klo_field2_gb, khi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb,
     &            klo_field1_gb:khi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb,
     &            klo_field2_gb:khi_field2_gb)
      real max_norm_diff
      real next_diff
      integer i,j,k
      real local_zone

c     initialize max_norm_diff
      max_norm_diff = -1d0

c       loop over included cells { 
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
                if (abs(field1(i,j,k)) < local_zone) then
                    next_diff = abs(field1(i,j,k) - field2(i,j,k))
                    if (next_diff .gt. max_norm_diff) then
                        max_norm_diff = next_diff
                    endif
                endif
            enddo
          enddo
        enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c*************


c***********************************************************************
      subroutine imposeMaskParallel(
     &  phi_masked,
     &  mask,
     &  phi,
     &  overlap,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  cur_klo_gb, cur_khi_gb)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer cur_klo_gb, cur_khi_gb
      
      real overlap
      real phi_masked(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real mask(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real phi(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
     
      integer i, j, k
      
      do k=cur_klo_gb,cur_khi_gb
        do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb
            
            if ((mask(i,j,k) - overlap) .gt. phi(i,j,k)) then
               phi_masked(i,j,k) = mask(i,j,k)
            else
               phi_masked(i,j,k) = phi(i,j,k)
            endif
          enddo
        enddo
      enddo
       
      return
      end
      
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine imposeMaskParallelVar(
     &  phi_masked,
     &  mask,
     &  phi,
     &  overlap,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  cur_klo_gb, cur_khi_gb)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer cur_klo_gb, cur_khi_gb
      
      real overlap(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real phi_masked(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real mask(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real phi(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
     
      integer i, j, k
      
      do k=cur_klo_gb,cur_khi_gb
        do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb
            
            if ((mask(i,j,k) - overlap(i,j,k)) .gt. phi(i,j,k)) then
               phi_masked(i,j,k) = mask(i,j,k)
            else
               phi_masked(i,j,k) = phi(i,j,k)
            endif
          enddo
        enddo
      enddo
       
      return
      end
      
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine copyDataParallel3d(
     &  phi_copy,
     &  phi_orig,
     &  ilo_gb, ihi_gb, 
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  cur_klo_gb, cur_khi_gb)
c***********************************************************************
c { begin subroutine
      implicit none
      
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer cur_klo_gb, cur_khi_gb
      
      real overlap
      real phi_copy(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
      real phi_orig(ilo_gb:ihi_gb,
     &         jlo_gb:jhi_gb,
     &         klo_gb:khi_gb)
     
      integer i, j, k
      
      do k=cur_klo_gb,cur_khi_gb
        do j=jlo_gb,jhi_gb
          do i=ilo_gb,ihi_gb
                phi_copy(i,j,k) = phi_orig(i,j,k)
          enddo
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
      subroutine qss3dVolumeRegionPhiLessThanZero(
     &  volume,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize volume to zero
      volume = 0.0d0
           
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
                phi_cur = phi(i,j,k)
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
      subroutine qss3dVolumeRegionPhiGreaterThanZero(
     &  volume,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize volume to zero
      volume = 0.0d0

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j,k)
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
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
c***********************************************************************
      subroutine qssHeaviside( val, x, eps )
c***********************************************************************
c { begin subroutine      
      implicit none
      
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      
      real x, val, eps, val1, val2, val3
      
      call basicHeaviside(val1, x+eps)
      call basicHeaviside(val2, eps-x)
      call basicHeaviside(val3, x-eps)
      
      val = (0.5d0 + x/(2*eps) + 
     &      (0.5d0*one_over_pi)*sin((pi*x)/(eps)))*val1*val2 + val3
      
      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine basicHeaviside( val, x )
c***********************************************************************
c { begin subroutine      
      implicit none
      
      real x, val
      
      if ( x .lt. 0) then
        val = 0d0
      elseif (x .gt. 0) then
        val = 1d0
      elseif (x .eq. 0) then
        val = 0.5d0
      endif
      
      return
      end
c } end subroutine
c***********************************************************************
        
c***********************************************************************

      subroutine imposeTrapVelocity3d( 
     &  phi_w,
     &  phi_nw, 
     &  normal_velocity, 
     &  curvature_coeff,
     &  vel_x, 
     &  vel_y,
     &  vel_z,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine      
      implicit none

      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real phi_w(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb) 
      real phi_nw(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb) 
      real normal_velocity(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)    
      real curvature_coeff(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real vel_x(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real vel_y(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)
      real vel_z(ilo_gb:ihi_gb,
     &           jlo_gb:jhi_gb,
     &           klo_gb:khi_gb)      
      real dx
      integer i, j, k
      real h1, h2, eps_var, x
      eps_var = 1.5 * dx    
      
c       loop over included cells {
        do k=klo_fb,khi_fb
          do j=jlo_fb,jhi_fb
            do i=ilo_fb,ihi_fb
c                call qssHeaviside(h1, -phi_nw(i,j,k), eps_var)
c                call qssHeaviside(h2, -phi_w(i,j,k), eps_var)
c                x = normal_velocity(i,j,k)
c                normal_velocity(i,j,k) = h1 * h2 * x
                if (phi_w(i,j,k) > eps_var) then
                    normal_velocity(i,j,k) = 0
                    curvature_coeff(i,j,k) = 0
                    vel_x(i,j,k) = 0
                    vel_y(i,j,k) = 0
                    vel_z(i,j,k) = 0
                endif                   
c               x = curvature_coeff(i,j,k)
c                curvature_coeff(i,j,k) = h1 * h2 * x
c                
c                x = vel_x(i,j,k)
c                vel_x(i,j,k) = h1 * h2 * x
                
c                x = vel_y(i,j,k)
c                vel_y(i,j,k) = h1 * h2 * x
                
c                x = vel_z(i,j,k)
c                vel_z(i,j,k) = h1 * h2 * x
            enddo
          enddo
        enddo
c       } end loop over grid
        return
      end
c***********************************************************************


        
     
     