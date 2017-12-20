c***********************************************************************
      subroutine qss3dRK1Step(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  klo_u_next_gb, khi_u_next_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer klo_u_next_gb, khi_u_next_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb,
     &                        klo_u_next_gb:khi_u_next_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      integer i, j, k
      real dt

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            u_next(i,j,k) = u_cur(i,j,k) + dt*rhs(i,j,k)

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine qss3dTVDRK2Stage2(
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  klo_u_stage2_gb, khi_u_stage2_gb,
     &  u_stage1,
     &  ilo_u_stage1_gb, ihi_u_stage1_gb,
     &  jlo_u_stage1_gb, jhi_u_stage1_gb,
     &  klo_u_stage1_gb, khi_u_stage1_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer klo_u_stage2_gb, khi_u_stage2_gb
      integer ilo_u_stage1_gb, ihi_u_stage1_gb
      integer jlo_u_stage1_gb, jhi_u_stage1_gb
      integer klo_u_stage1_gb, khi_u_stage1_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb,
     &                          klo_u_stage2_gb:khi_u_stage2_gb)
      real u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb,
     &                          jlo_u_stage1_gb:jhi_u_stage1_gb,
     &                          klo_u_stage1_gb:khi_u_stage1_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      integer i, j, k
      real dt

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            u_stage2(i,j,k) = 0.75d0*u_cur(i,j,k) 
     &                      + 0.25d0*(u_stage1(i,j,k) + dt*rhs(i,j,k))

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine qss3dTVDRK3Stage3(
     &  u_next,
     &  ilo_u_next_gb, ihi_u_next_gb,
     &  jlo_u_next_gb, jhi_u_next_gb,
     &  klo_u_next_gb, khi_u_next_gb,
     &  u_stage2,
     &  ilo_u_stage2_gb, ihi_u_stage2_gb,
     &  jlo_u_stage2_gb, jhi_u_stage2_gb,
     &  klo_u_stage2_gb, khi_u_stage2_gb,
     &  u_cur,
     &  ilo_u_cur_gb, ihi_u_cur_gb,
     &  jlo_u_cur_gb, jhi_u_cur_gb,
     &  klo_u_cur_gb, khi_u_cur_gb,
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  klo_rhs_gb, khi_rhs_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb,
     &  dt)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_u_next_gb, ihi_u_next_gb
      integer jlo_u_next_gb, jhi_u_next_gb
      integer klo_u_next_gb, khi_u_next_gb
      integer ilo_u_stage2_gb, ihi_u_stage2_gb
      integer jlo_u_stage2_gb, jhi_u_stage2_gb
      integer klo_u_stage2_gb, khi_u_stage2_gb
      integer ilo_u_cur_gb, ihi_u_cur_gb
      integer jlo_u_cur_gb, jhi_u_cur_gb
      integer klo_u_cur_gb, khi_u_cur_gb
      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer klo_rhs_gb, khi_rhs_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb, klo_fb, khi_fb
      real u_next(ilo_u_next_gb:ihi_u_next_gb,
     &                        jlo_u_next_gb:jhi_u_next_gb,
     &                        klo_u_next_gb:khi_u_next_gb)
      real u_stage2(ilo_u_stage2_gb:ihi_u_stage2_gb,
     &                          jlo_u_stage2_gb:jhi_u_stage2_gb,
     &                          klo_u_stage2_gb:khi_u_stage2_gb)
      real u_cur(ilo_u_cur_gb:ihi_u_cur_gb,
     &                       jlo_u_cur_gb:jhi_u_cur_gb,
     &                       klo_u_cur_gb:khi_u_cur_gb)
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &                     jlo_rhs_gb:jhi_rhs_gb,
     &                     klo_rhs_gb:khi_rhs_gb)
      integer i, j, k
      real dt
      real one_third, two_thirds
      parameter (one_third = 1.d0/3.d0)
      parameter (two_thirds = 2.d0/3.d0)

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

            u_next(i,j,k) = one_third*u_cur(i,j,k)
     &                    + two_thirds*( u_stage2(i,j,k) 
     &                                 + dt*rhs(i,j,k) )

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
