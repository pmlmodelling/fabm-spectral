#include "fabm_driver.h"

module fabm_ozone

   ! Spectral light model based on Bird 1984, Bird & Riordan 1986, Gregg &  Carder 1990, Gregg & Casey 2009
   ! Copyright PML 2018

   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   type,extends(type_base_model), public :: type_ozone
      type (type_horizontal_dependency_id)          :: id_lon, id_lat
      type (type_global_dependency_id)              :: id_yearday
      type (type_horizontal_diagnostic_variable_id) :: id_O3
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ozone), intent(inout), target :: self
      integer,            intent(in)            :: configunit

      call self%register_dependency(self%id_lon, standard_variables%longitude)
      call self%register_dependency(self%id_lat, standard_variables%latitude)
      call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)

      call self%register_diagnostic_variable(self%id_O3, 'O3', 'kg m-2', 'atmospheric ozone', &
        standard_variable=type_horizontal_standard_variable('atmosphere_mass_content_of_ozone', 'kg m-2'), &
        source=source_do_surface)
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ozone), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: longitude, latitude, yearday, O3

      _HORIZONTAL_LOOP_BEGIN_
          _GET_HORIZONTAL_(self%id_lon, longitude)
          _GET_HORIZONTAL_(self%id_lat, latitude)
          _GET_GLOBAL_(self%id_yearday, yearday)
          
          ! Use Van Heuklon (1979) to account for seasonal and geographical variation.
          ! Cf Bird 1984 p 466 reports O3=0.344 atm-cm
          O3 = estimate_ozone(longitude, latitude, yearday)

          ! From matm-cm to mol m-2 to kg m-2
          O3 = O3 * 0.4462e-3_rk * 48._rk / 1000

          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O3, O3)
      _HORIZONTAL_LOOP_END_
   end subroutine

   function estimate_ozone(longitude, latitude, yearday) result(O3)
      ! Estimate ozone concentration based on Van Heuklon (1979)
      real(rk), intent(in) :: longitude, latitude, yearday
      real(rk) :: O3
     
      real(rk), parameter :: pi = 3.14159265358979323846_rk
      real(rk), parameter :: deg2rad = pi / 180._rk

      ! Model parameters given in Van Heuklon 1979, Table 1
      real(rk), parameter :: D = 0.9865_rk, G = 20._rk, J = 235._rk
      real(rk), parameter :: A_N = 150._rk, beta_N = 1.28_rk, C_N = 40._rk, F_N = -30._rk,    H_N = 3._rk, I_NE = 20._rk, I_NW = 0._rk
      real(rk), parameter :: A_S = 100._rk, beta_S = 1.5_rk,  C_S = 30._rk, F_S = 152.625_rk, H_S = 2._rk, I_S = -75._rk

      real(rk) :: A, beta, C, E, F, H, I, phi, lambda

      E = yearday
      phi = latitude
      lambda = -180._rk + mod(longitude + 180._rk, 360._rk)

      if (latitude > 0._rk) then
         ! Northern hemisphere
         A = A_N
         beta = beta_N
         C = C_N
         F = F_N
         H = H_N
         if (lambda > 0._rk) then
            ! Eastern hemisphere
            I = I_NE
         else
            ! Western hemisphere
            I = I_NW
         end if
      else
         ! Southern hemisphere
         A = A_S
         beta = beta_S
         C = C_S
         F = F_S
         H = H_S
         I = I_S
      end if

      ! Van Heuklon (1979), Eq 4 - note conversion from degrees to radians
      O3 = J + (A + C * sin(D * (E + F) * deg2rad) + G * sin(H * (lambda + I) * deg2rad)) * sin(beta * phi * deg2rad)**2
   end function

end module