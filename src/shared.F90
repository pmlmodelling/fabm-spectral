#include "fabm_driver.h"

module fabm_spectral_shared

   use fabm_types, only: rk

   implicit none

   public

   real(rk), parameter :: pi = 3.14159265358979323846_rk
   real(rk), parameter :: deg2rad = pi / 180._rk
   real(rk), parameter :: rad2deg = 180._rk / pi

   real(rk), parameter :: Planck = 6.62606957e-34_rk  ! Planck constant (m2 kg/s)
   real(rk), parameter :: lightspeed = 299792458_rk   ! Speed of light (m/s)
   real(rk), parameter :: Avogadro = 6.02214129e23_rk ! Avogadro constant (/mol)

   interface interp
      module procedure interp_0d
      module procedure interp_1d
      module procedure interp_0d_scalar
      module procedure interp_1d_scalar
   end interface

contains

   subroutine calculate_integral_weights(xl, xr, n, x, w)
      integer,  intent(in)  :: n
      real(rk), intent(in)  :: x(n), xl, xr
      real(rk), intent(out) :: w(n)

      integer  :: i
      real(rk) :: deltax, f

      w = 0._rk

      if (x(1) > xr) then
         ! All points to right of desired range.
         w(1) = xr - xl
         return
      elseif (x(n) < xl) then
         ! All points to left of desired range.
         w(n) =  xr - xl
         return
      end if

      do i = 2, n
         deltax = x(i) - x(i-1)
         if (x(i) >= xl .and. x(i) <= xr) then
            ! Right-hand point in desired range.
            if (x(i-1) >= xl) then
               ! Whole segment in range
               ! Integral: 0.5*(y1+y2)*(x2-x1)
               w(i-1:i) = w(i-1:i) + 0.5_rk * deltax
            else
               ! Segment crosses left boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xl-x1) = y1*(1-(xl-x1)/(x2-x1)) + y2*(xl-x1)/(x2-x1)
               ! integral = 0.5*(yb+y2)*(x2-xl)
               ! yb+y2 = y1*(1-(xl-x1)/(x2-x1)) + y2*(1+(xl-x1)/(x2-x1))
               f = (xl - x(i-1)) / deltax
               w(i-1) = w(i-1) + (1.0_rk - f) * (x(i) - xl) * 0.5_rk
               w(i  ) = w(i  ) + (1.0_rk + f) * (x(i) - xl) * 0.5_rk
            end if
         elseif (x(i) > xr .and. x(i-1) < xr) then
            ! Right-hand point beyond desired range, left-hand point before right range boundary.
            if (x(i-1) >= xl) then
               ! Segment crosses right boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xr-x1) = y1*(1-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               ! integral = 0.5*(y1+yb)*(xr-x1)
               ! y1+yb = y1*(2-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               f = (xr - x(i-1)) / deltax
               w(i-1) = w(i-1) + (2.0_rk - f) * (xr - x(i-1)) * 0.5_rk
               w(i  ) = w(i  ) + (         f) * (xr - x(i-1)) * 0.5_rk
            else
               ! Segment crosses both boundaries
               ! y at centre: yc = y1 + (y2-y1)/(x2-x1)*(0.5*(xl+xr)-x1) = y1*(1-(0.5*(xl+xr)-x1)/(x2-x1)) + y2*(0.5*(xl+xr)-x1)/(x2-x1)
               ! integral: (xr-xl)*yc
               f = (0.5_rk*(xl+xr)-x(i-1))/deltax
               w(i-1) = w(i-1) + (1._rk-f)*(xr-xl)
               w(i  ) = w(i  ) + (      f)*(xr-xl)
            end if
         end if
      end do
      if (x(1)>xl) w(1) = w(1) + (x(1)-xl)
      if (x(n)<xr) w(n) = w(n) + (xr-x(n))
   end subroutine

   subroutine interp_0d(nsource,x,y,ntarget,targetx,targety)
      ! 1D interpolation, extrapolates beyond boundaries
      integer,                    intent(in)  :: nsource,ntarget
      real(rk),dimension(nsource),intent(in)  :: x,y
      real(rk),dimension(ntarget),intent(in)  :: targetx
      real(rk),dimension(ntarget),intent(out) :: targety
      integer :: i,j
      real(rk) :: frac

      i = 1
      do j = 1,ntarget
         do while (i+1<nsource)
            if (x(i+1)>=targetx(j)) exit
            i = i+1
         end do
         frac = (targetx(j)-x(i))/(x(i+1)-x(i))
         targety(j) = y(i) + frac*(y(i+1)-y(i))
      end do
   end subroutine

   subroutine interp_0d_scalar(nsource,x,y,targetx,targety)
      ! 1D interpolation, extrapolates beyond boundaries
      integer,                    intent(in)  :: nsource
      real(rk),dimension(nsource),intent(in)  :: x,y
      real(rk),                   intent(in)  :: targetx
      real(rk),                   intent(out) :: targety
      integer :: i
      real(rk) :: frac

      i = 1
      do while (i+1<nsource)
         if (x(i+1)>=targetx) exit
         i = i+1
      end do
      frac = (targetx-x(i))/(x(i+1)-x(i))
      targety = y(i) + frac*(y(i+1)-y(i))
   end subroutine

   subroutine interp_1d(m,nsource,x,y,ntarget,targetx,targety)
      ! 1D interpolation, extrapolates beyond boundaries
      integer,                      intent(in)  :: nsource,ntarget,m
      real(rk),dimension(nsource),  intent(in)  :: x
      real(rk),dimension(m,nsource),intent(in)  :: y
      real(rk),dimension(ntarget),  intent(in)  :: targetx
      real(rk),dimension(m,ntarget),intent(out) :: targety
      integer :: i,j
      real(rk) :: frac

      i = 1
      do j = 1,ntarget
         do while (i+1<nsource)
            if (x(i+1)>=targetx(j)) exit
            i = i+1
         end do
         frac = (targetx(j)-x(i))/(x(i+1)-x(i))
         targety(:,j) = y(:,i) + frac*(y(:,i+1)-y(:,i))
      end do
   end subroutine

   subroutine interp_1d_scalar(m,nsource,x,y,targetx,targety)
      ! 1D interpolation, extrapolates beyond boundaries
      integer,                      intent(in)  :: nsource,m
      real(rk),dimension(nsource),  intent(in)  :: x
      real(rk),dimension(m,nsource),intent(in)  :: y
      real(rk),                     intent(in)  :: targetx
      real(rk),dimension(m),        intent(out) :: targety
      integer :: i
      real(rk) :: frac

      i = 1
      do while (i+1<nsource)
         if (x(i+1)>=targetx) exit
         i = i+1
      end do
      frac = (targetx-x(i))/(x(i+1)-x(i))
      targety(:) = y(:,i) + frac*(y(:,i+1)-y(:,i))
   end subroutine

end module
