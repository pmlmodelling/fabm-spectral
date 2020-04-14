#include "fabm_driver.h"

module fabm_relhum

   ! Copyright PML 2018

   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   type,extends(type_base_model), public :: type_relhum
      type (type_horizontal_dependency_id)          :: id_airtemp, id_sphum, id_airpres
      type (type_horizontal_diagnostic_variable_id) :: id_relhum
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_relhum), intent(inout), target :: self
      integer,             intent(in)            :: configunit

      call self%register_dependency(self%id_sphum, standard_variables%surface_specific_humidity)
      call self%register_dependency(self%id_airtemp, standard_variables%surface_temperature)
      call self%register_dependency(self%id_airpres, standard_variables%surface_air_pressure)

      call self%register_diagnostic_variable(self%id_relhum, 'relhum', '-', 'relative_humidity', &
        standard_variable=type_horizontal_standard_variable('relative_humidity', '-'), source=source_do_surface)
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_relhum), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      ! Lowe 1977 coefficients applicable for air temperature in degrees Celsius (taken from GOTM's humidity.F90)
      real(rk), parameter :: a1=6.107799961
      real(rk), parameter :: a2=4.436518521e-1
      real(rk), parameter :: a3=1.428945805e-2
      real(rk), parameter :: a4=2.650648471e-4
      real(rk), parameter :: a5=3.031240396e-6
      real(rk), parameter :: a6=2.034080948e-8
      real(rk), parameter :: a7=6.136820929e-11

      ! From GOTM's airsea_variables.F90
      real(rk), parameter :: const06=0.62198
      
      real(rk) :: qa, ta, airp, relhum, ea_sat, ea

      _HORIZONTAL_LOOP_BEGIN_
          _GET_HORIZONTAL_(self%id_sphum, qa)      ! Specific humidity (kg kg-1)
          _GET_HORIZONTAL_(self%id_airtemp, ta)    ! 2m air temperature (degrees Celsius)
          _GET_HORIZONTAL_(self%id_airpres, airp)  ! Surface air pressure (Pa)

          ! Get saturation vapor pressure from air temperature (Lowe 1977; taken from GOTM's humidity.F90)
          ea_sat = a1 +ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*(a6+ta*a7)))))
          ea_sat = ea_sat * 100.0 ! Conversion millibar --> Pascal

          ! from specific humidity qa to vapour pressure ea
          ! GOTM gives inverse as qa = const06*ea/(airp-0.377*ea)
          ea = qa * airp / (const06 + 0.377*qa)
          relhum = ea / ea_sat

          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_relhum, relhum)
      _HORIZONTAL_LOOP_END_
   end subroutine

end module