#include "fabm_driver.h"

module fabm_spectral

   ! Spectral light model based on Bird 1984, Bird & Riordan 1986, Gregg &  Carder 1990, Gregg & Casey 2009
   ! Copyright PML 2018

   use fabm_types
   use fabm_standard_variables
   use fabm_expressions, only: temporal_mean

   implicit none

   private

   public :: slingo, nlambda_slingo, lambda_slingo

   type,extends(type_base_model), public :: type_spectral
      type (type_diagnostic_variable_id) :: id_swr, id_uv, id_par, id_par_E
      type (type_horizontal_diagnostic_variable_id) :: id_swr_sf, id_par_sf, id_uv_sf, id_par_E_sf, id_swr_dif_sf
      type (type_horizontal_diagnostic_variable_id) :: id_swr_sf_w, id_par_sf_w, id_uv_sf_w, id_par_E_sf_w
      type (type_horizontal_diagnostic_variable_id) :: id_alpha_a, id_beta_a, id_omega_a, id_F_a
      type (type_horizontal_dependency_id) :: id_lon, id_lat, id_cloud, id_wind_speed, id_airpres, id_relhum, id_lwp, id_O3, id_WV, id_mean_wind_speed, id_visibility, id_air_mass_type
      type (type_global_dependency_id) :: id_yearday
      type (type_dependency_id) :: id_h
      integer :: nlambda
      type (type_horizontal_diagnostic_variable_id), dimension(:), allocatable :: id_surface_band_dir, id_surface_band_dif
      type (type_diagnostic_variable_id), dimension(:), allocatable :: id_band_dir, id_band_dif
      real(rk), dimension(:), allocatable :: lambda, lambda_bounds, par_weights, par_E_weights, swr_weights, uv_weights, F
      real(rk), dimension(:), allocatable :: exter
      real(rk), dimension(:), allocatable :: a_o, a_u, a_v, tau_r
      real(rk), dimension(:), allocatable :: a_w, b_w
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

   real(rk), parameter :: pi = 3.14159265358979323846_rk
   real(rk), parameter :: deg2rad = pi / 180._rk
   real(rk), parameter :: rad2deg = 180._rk / pi

#include "slingo_const.inc"
#include "water_const.inc"
#include "oasim_const.inc"
#include "astm_const.inc"
#include "birdrior1986_const.inc"

   interface interp
      module procedure interp_0d
      module procedure interp_1d
      module procedure interp_0d_scalar
      module procedure interp_1d_scalar
   end interface

contains

   subroutine initialize(self, configunit)
      class (type_spectral), intent(inout), target :: self
      integer,               intent(in)            :: configunit

      integer :: l
      integer :: lambda_method
      real(rk) :: lambda_min, lambda_max
      logical :: save_spectra
      character(len=8) :: strwavelength, strindex

      integer, parameter :: exter_source = 2

      ! Coefficients for wavelength dependence of foam reflectance (Eqs A11, A12 in Gregg & Casey 2009)
      real(rk), parameter :: a0 = 0.9976_rk
      real(rk), parameter :: a1 = 0.2194_rk
      real(rk), parameter :: a2 = 0.0554_rk
      real(rk), parameter :: a3 = 0.0067_rk
      real(rk), parameter :: b0 = 5.026_rk
      real(rk), parameter :: b1 = -0.0114_rk
      real(rk), parameter :: b2 = 9.552e-6_rk
      real(rk), parameter :: b3 = -2.698e-9_rk

      real(rk), parameter :: Planck = 6.62606957e-34     ! Planck constant (m2 kg/s)
      real(rk), parameter :: lightspeed = 299792458_rk   ! Speed of light (m/s)
      real(rk), parameter :: Avogadro = 6.02214129e23_rk ! Avogadro constant (/mol)

      real(rk) :: log_T_w

      call self%get_parameter(lambda_method, 'lambda_method', '', 'choice of wavebands (0: custom range, 1: OASIM)', default=1)
      select case (lambda_method)
      case (0)
         call self%get_parameter(self%nlambda, 'nlambda', '', 'number of wavebands')
         call self%get_parameter(lambda_min, 'lambda_min', 'nm', 'minimum wavelength', minimum=0._rk)
         call self%get_parameter(lambda_max, 'lambda_max', 'nm', 'maximum wavelength', minimum=0._rk)
         allocate(self%lambda(self%nlambda), self%lambda_bounds(self%nlambda + 1))
         do l = 1, self%nlambda + 1
            self%lambda_bounds(l) = lambda_min + (l - 1) * (lambda_max - lambda_min) / self%nlambda
         end do
         self%lambda = (self%lambda_bounds(1:self%nlambda) + self%lambda_bounds(2:self%nlambda + 1)) / 2
      case (1)
         self%nlambda = nlambda_oasim
         allocate(self%lambda(self%nlambda), self%lambda_bounds(self%nlambda + 1))
         self%lambda(:) = lambda_oasim
      end select
      call self%get_parameter(save_spectra, 'save_spectra', '', 'save full spectral irradiance', default=.false.)

      ! Find wavelength bounds of photosynthetically active radiation
      allocate(self%par_weights(self%nlambda))
      allocate(self%swr_weights(self%nlambda))
      allocate(self%uv_weights(self%nlambda))
      allocate(self%par_E_weights(self%nlambda))
      call calculate_integral_weights(400._rk, 700._rk, self%nlambda, self%lambda, self%par_weights)
      call calculate_integral_weights(300._rk, 4000._rk, self%nlambda, self%lambda, self%swr_weights)
      call calculate_integral_weights(300._rk, 400._rk, self%nlambda, self%lambda, self%uv_weights)
      self%par_E_weights(:) = self%par_weights * self%lambda /(Planck*lightspeed)/Avogadro*1e-3_rk ! divide by 1e9 to go from nm to m, multiply by 1e6 to go from mol to umol

      call self%register_dependency(self%id_lon, standard_variables%longitude)
      call self%register_dependency(self%id_lat, standard_variables%latitude)
      call self%register_dependency(self%id_wind_speed, standard_variables%wind_speed)
      call self%register_dependency(self%id_cloud, standard_variables%cloud_area_fraction)
      call self%register_dependency(self%id_airpres, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
      call self%register_dependency(self%id_relhum, type_horizontal_standard_variable('relative_humidity', '-'))
      call self%register_dependency(self%id_lwp, type_horizontal_standard_variable('atmosphere_mass_content_of_cloud_liquid_water', 'kg m-2'))
      call self%register_dependency(self%id_O3, type_horizontal_standard_variable('atmosphere_mass_content_of_ozone', 'kg m-2'))
      call self%register_dependency(self%id_WV, type_horizontal_standard_variable('atmosphere_mass_content_of_water_vapor', 'kg m-2'))
      call self%register_dependency(self%id_visibility, type_horizontal_standard_variable('visibility_in_air', 'm'))
      call self%register_dependency(self%id_air_mass_type, type_horizontal_standard_variable('aerosol_air_mass_type', '-'))
      call self%register_dependency(self%id_mean_wind_speed, temporal_mean(self%id_wind_speed, period=86400._rk, resolution=3600._rk))

      !call self%register_diagnostic_variable(self%id_swr, 'swr',     'W/m^2',     'shortwave radiation')
      !call self%register_diagnostic_variable(self%id_uv,  'uv',      'W/m^2',     'ultraviolet radiative flux')
      !call self%register_diagnostic_variable(self%id_par, 'par',     'W/m^2',     'photosynthetically active radiation')
      !call self%register_diagnostic_variable(self%id_par_E, 'par_E', 'umol/m^2/s','photosynthetic photon flux density')
      call self%register_diagnostic_variable(self%id_swr_sf, 'swr_sf',     'W/m^2',     'downward shortwave radiation in air')
      call self%register_diagnostic_variable(self%id_swr_dif_sf, 'swr_dif_sf',     'W/m^2', 'downward diffuse shortwave radiation in air')
      call self%register_diagnostic_variable(self%id_uv_sf,  'uv_sf',      'W/m^2',     'downward ultraviolet radiative flux in air')
      call self%register_diagnostic_variable(self%id_par_sf, 'par_sf',     'W/m^2',     'downward photosynthetically active radiation in air')
      call self%register_diagnostic_variable(self%id_par_E_sf, 'par_E_sf', 'umol/m^2/s','downward photosynthetic photon flux density in air')
      call self%register_diagnostic_variable(self%id_swr_sf_w, 'swr_sf_w',     'W/m^2',     'downward shortwave radiation in water')
      call self%register_diagnostic_variable(self%id_uv_sf_w,  'uv_sf_w',      'W/m^2',     'downward ultraviolet radiative flux in water')
      call self%register_diagnostic_variable(self%id_par_sf_w, 'par_sf_w',     'W/m^2',     'downward photosynthetically active radiation in water')
      call self%register_diagnostic_variable(self%id_par_E_sf_w, 'par_E_sf_w', 'umol/m^2/s','downward photosynthetic photon flux density in water')

      call self%register_diagnostic_variable(self%id_alpha_a, 'alpha_a', '-', 'aerosol Angstrom exponent')
      call self%register_diagnostic_variable(self%id_beta_a, 'beta_a', '-', 'aerosol scale factor for optical thickness')
      call self%register_diagnostic_variable(self%id_omega_a, 'omega_a', '-', 'aerosol single scattering albedo')
      call self%register_diagnostic_variable(self%id_F_a, 'F_a', '-', 'aerosol forward scattering probability')

      ! Interpolate absorption and scattering spectra to user wavelength grid
      allocate(self%a_w(self%nlambda), self%b_w(self%nlambda))
      call interp(nlambda_w, lambda_w, a_w, self%nlambda, self%lambda, self%a_w)
      call interp(nlambda_w, lambda_w, b_w, self%nlambda, self%lambda, self%b_w)

      allocate(self%exter(self%nlambda), self%a_o(self%nlambda), self%a_v(self%nlambda), self%a_u(self%nlambda), self%tau_r(self%nlambda))
      if (exter_source == 1) then
         call interp(nlambda_oasim, lambda_oasim, ET_oasim, self%nlambda, self%lambda, self%exter)
      else
         call interp(nlambda_astm, lambda_astm, ET_astm, self%nlambda, self%lambda, self%exter)
      end if
      call interp(nlambda_oasim, lambda_oasim, a_o_oasim, self%nlambda, self%lambda, self%a_o)
      call interp(nlambda_oasim, lambda_oasim, a_v_oasim, self%nlambda, self%lambda, self%a_v)
      call interp(nlambda_oasim, lambda_oasim, a_u_oasim, self%nlambda, self%lambda, self%a_u)
      !call interp(nlambda_oasim, lambda_birdrior1986, a_u_birdrior1986, self%nlambda, self%lambda, self%a_u)

      ! Rayleigh optical thickness (Eq 2 Bird 1984, Eq 4 Bird & Riordan 1986, Eq 15 Gregg & Carder 1990)
      !call interp(nlambda_oasim, lambda_oasim, tau_r_oasim, self%nlambda, self%lambda, self%tau_r)
      self%tau_r = 1.0_rk / (115.6406_rk * (self%lambda/1000)**4 - 1.335_rk * (self%lambda/1000)**2)

      ! Protect against negative absorption coefficients produced by linear extrapolation
      self%a_o = max(0._rk, self%a_o)
      self%a_v = max(0._rk, self%a_v)
      self%a_u = max(0._rk, self%a_u)

      ! Wavelength dependence of foam reflectance (Eqs A11, A10 Gregg & Casey 2009)
      allocate(self%F(self%nlambda))
      do l = 1, self%nlambda
         if (self%lambda(l) < 900) then
            ! Note: close to 900 nm the expression below returns negative values.
            ! We clip to zero in line with Fig A1 Gregg & Casey 2009
            log_T_w = -(self%a_w(l) + 0.5_rk * self%b_w(l))
            self%F(l) = max(0._rk, a0 + a1 * log_T_w + a2 * log_T_w**2 + a3 * log_T_w**3)
         else
            ! Note: above 1700 nm the expression below returns negative values.
            ! We clip to zero in line with Fig A1 Gregg & Casey 2009
            self%F(l) = max(0._rk, b0 + b1 * self%lambda(l) + b2 * self%lambda(l)**2 + b3 * self%lambda(l)**3)
         end if
      end do

      if (save_spectra) then
         allocate(self%id_surface_band_dir(self%nlambda))
         allocate(self%id_surface_band_dif(self%nlambda))
         allocate(self%id_band_dir(self%nlambda))
         allocate(self%id_band_dif(self%nlambda))
         do l = 1, self%nlambda
            if (self%lambda(l) < 1000._rk) then
               write(strwavelength, '(f5.1)') self%lambda(l)
            else
               write(strwavelength, '(f6.1)') self%lambda(l)
            end if
            write(strindex, '(i0)') l
            call self%register_diagnostic_variable(self%id_surface_band_dir(l), 'dir_sf_band' // trim(strindex), 'W/m2/nm', 'downward direct irradiance in air @ ' // trim(strwavelength) // ' nm')
            call self%register_diagnostic_variable(self%id_surface_band_dif(l), 'dif_sf_band' // trim(strindex), 'W/m2/nm', 'downward diffuse irradiance in air @ ' // trim(strwavelength) // ' nm')
            !call self%register_diagnostic_variable(self%id_band_dir(l), 'dir_band' // trim(strindex), 'W/m2/nm', 'direct irradiance @ ' // trim(strwavelength) // ' nm')
            !call self%register_diagnostic_variable(self%id_band_dif(l), 'dif_band' // trim(strindex), 'W/m2/nm', 'diffuse irradiance @ ' // trim(strwavelength) // ' nm')
         end do
      end if
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_spectral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      !real(rk), parameter :: W0 = 0.928_rk  ! Single scattering albedo of the aerosol (Bird 1984 Eq 15)
      !real(rk), parameter :: fa = 0.82_rk   ! Forward to total scattering ratio of the aerosol (Bird 1984 Eq 15)
      real(rk), parameter :: pres0 = 101300._rk ! reference air pressure in Pa (Bird 1984 p461, Bird & Riordan p89)
      !real(rk), parameter :: ga = 0.05_rk ! ground albedo
      real(rk), parameter :: H_oz = 22._rk ! Height of maximum ozone concentration (km)
      real(rk), parameter :: r_e = (10._rk + 11.8_rk) / 2 ! equivalent radius of cloud drop size distribution (um) based on mean of Kiehl et al. & Han et al. (cf OASIM)

      real(rk) :: longitude, latitude, yearday, cloud_cover, wind_speed, airpres, relhum, LWP, water_vapour, WV, WM, visibility, AM
      real(rk) :: days, hour, theta, costheta, alpha_a, beta_a
      integer :: l
      real(rk) :: M, M_prime, M_oz
      real(rk) :: O3
      real(rk), dimension(self%nlambda) :: direct, diffuse, spectrum  ! Spectra at top of the water column (with refraction and reflection accounted for)
      real(rk) :: par_J, swr_J, uv_J, par_E, F_a, omega_a
      real(rk), dimension(self%nlambda) :: tau_a, T_a, T_oz, T_w, T_u, T_r, T_aa, T_as
      real(rk), dimension(self%nlambda) :: T_g, T_dclr, T_sclr, T_dcld, T_scld
      real(rk), dimension(self%nlambda) :: rho_d, rho_s

      _HORIZONTAL_LOOP_BEGIN_
          _GET_HORIZONTAL_(self%id_lon, longitude)
          _GET_HORIZONTAL_(self%id_lat, latitude)
          _GET_GLOBAL_(self%id_yearday, yearday)
          _GET_HORIZONTAL_(self%id_cloud, cloud_cover)      ! Cloud cover (fraction, 0-1)
          _GET_HORIZONTAL_(self%id_wind_speed, wind_speed)  ! Wind speed @ 10 m above surface (m/s)
          _GET_HORIZONTAL_(self%id_airpres, airpres)        ! Surface air pressure (Pa)
          _GET_HORIZONTAL_(self%id_relhum, relhum)          ! Relative humidity (-)
          _GET_HORIZONTAL_(self%id_lwp, LWP)                ! Cloud liquid water content (kg m-2)
          _GET_HORIZONTAL_(self%id_O3, O3)                  ! Ozone content (kg m-2)
          _GET_HORIZONTAL_(self%id_wv, water_vapour)        ! Total precipitable water vapour (kg m-2) - equivalent to mm
          _GET_HORIZONTAL_(self%id_mean_wind_speed, WM)     ! Daily mean wind speed @ 10 m above surface (m/s)
          _GET_HORIZONTAL_(self%id_visibility, visibility)  ! Visibility (m)
          _GET_HORIZONTAL_(self%id_air_mass_type, AM)       ! Aerosol air mass type (1: open ocean, 10: continental)

          WV = water_vapour / 10                ! from kg m-2 to cm
          O3 = O3 * (1000 / 48._rk) / 0.4462_rk ! from kg m-2 to mol m-2, then from mol m-2 to atm cm (Basher 1982)
          days = floor(yearday) + 1.0_rk
          hour = mod(yearday, 1.0_rk) * 24.0_rk

          ! Calculate zenith angle (in radians)
          theta = zenith_angle(days, hour, longitude, latitude)
          theta = min(theta, 0.5_rk * pi)  ! Restrict the input zenith angle between 0 and pi/2
          costheta = cos(theta)

          ! Atmospheric path length, a.k.a. relative air mass (Eq 3 Bird 1984, Eq 5 Bird & Riordan 1986, Eq 13 Gregg & Carder 1990, Eq A5 in Casey & Gregg 2009)
          ! Note this should always exceed 1, but as it is an approximation it does not near theta -> 0 (Tomasi et al. 1998 p14). Hence the max operator.
          M = max(1._rk, 1._rk / (costheta + 0.15_rk * (93.885_rk - theta * rad2deg)**(-1.253_rk)))

          ! Pressure-corrected atmospheric path length (Eq A6 Casey & Gregg 2009)
          M_prime = M * airpres / pres0

         ! Atmospheric path length for ozone
         ! Eq 10, Bird & Riordan 1986, Eq 14 Gregg & Carder 1990; NB 6370 is the earth's radius in km
         ! See also Tomasi et al. 1998 Eq 5
         M_oz = (1._rk + H_oz / 6370._rk) / sqrt(costheta**2 + 2 * H_oz / 6370._rk)

         ! Transmittance due to ozone absorption (Eq 8 Bird 1984, Eq 9 Bird & Riordan 1986, Eq 17 Gregg & Carder 1990)
         T_oz = exp(-self%a_o * O3 * M_oz)

         ! Transmittance due to water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
         ! Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
         T_w = exp((-0.2385_rk * self%a_v * WV * M) / (1._rk + 20.07_rk * self%a_v * WV * M)**0.45_rk)

         ! Transmittance due to uniformly mixed gas absorption - SHOULD use pressure corrected airmass
         ! Eq 10 Bird 1984, Eq 11 Bird and Riordan 1986, Eq 18 Gregg & Carder 1990
         ! Bird & Riordan use 118.93 rather than 118.3, but state 118.3 should be used in the future.
         T_u = exp(-1.41_rk * self%a_u * M_prime / (1._rk + 118.3_rk * self%a_u * M_prime)**0.45_rk)

         ! Transmittance terms that apply for both cloudy and clear skies.
         T_g = T_oz * T_w * T_u

         ! -------------------
         ! clear skies part
         ! -------------------

         ! Transmittance due to Rayleigh scattering (use precomputed optical thickness)
         T_r = exp(-M_prime * self%tau_r)

         ! Transmittance due to aerosol absorption (Eq 26 Gregg & Carder 1990)
         call navy_aerol_model(AM, WM, wind_speed, relhum, visibility, costheta, alpha_a, beta_a, F_a, omega_a)
         tau_a = beta_a * self%lambda**(-alpha_a)
         T_a = exp(-tau_a * M)

         ! Direct transmittance
         T_dclr = T_r * T_a

         ! Separate absorption and scattering components of aerosol transmittance
         T_aa = exp(-(1._rk - omega_a) * tau_a * M)
         T_as = exp(-omega_a * tau_a * M)
         T_sclr = T_aa * 0.5_rk * (1._rk - T_r**0.95_rk) + T_r**1.5_rk * T_aa * F_a * (1 - T_as)

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_alpha_a, alpha_a)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_beta_a, beta_a)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_omega_a, omega_a)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_F_a, F_a)

         ! -------------------
         ! cloudy skies part
         ! -------------------

         ! Transmittance due to absorption and scattering by clouds
         call slingo(costheta, max(0._rk, LWP) * 1000, r_e, self%nlambda, self%lambda, T_dcld, T_scld)

         ! Diffuse and direct irradiance streams (Eqs 1, 2 Gregg & Casey 2009)
         ! These combine terms for clear and cloudy skies, weighted by cloud cover fraction
         direct = self%exter * costheta * T_g * ((1._rk - cloud_cover) * T_dclr + cloud_cover * T_dcld)
         diffuse = self%exter * costheta * T_g * ((1._rk - cloud_cover) * T_sclr + cloud_cover * T_scld)

         spectrum = direct + diffuse
         par_J = sum(self%par_weights * spectrum)
         swr_J = sum(self%swr_weights * spectrum)
         par_E = sum(self%par_E_weights * spectrum)
         uv_J  = sum(self%uv_weights * spectrum)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_sf, par_J)  ! Photosynthetically Active Radiation (W/m2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_E_sf,par_E) ! Photosynthetically Active Radiation (umol/m2/s)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_swr_sf,  swr_J) ! Total shortwave radiation (W/m2) [up to 4000 nm]
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uv_sf, uv_J)    ! UV (W/m2)

         swr_J = sum(self%swr_weights * diffuse)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_swr_dif_sf, swr_J)

         if (allocated(self%id_surface_band_dir)) then
            do l = 1, self%nlambda
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_band_dir(l), direct(l))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_band_dif(l), diffuse(l))
            end do
         end if

         ! Sea surface reflectance
         call reflectance(self%nlambda, self%F, theta, wind_speed, rho_d, rho_s)

         ! Incorporate the loss due to reflectance
         direct = direct * (1._rk - rho_d)
         diffuse = diffuse * (1._rk - rho_s)
         spectrum = direct + diffuse

         par_J = sum(self%par_weights * spectrum)
         swr_J = sum(self%swr_weights * spectrum)
         par_E = sum(self%par_E_weights * spectrum)
         uv_J  = sum(self%uv_weights * spectrum)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_sf_w, par_J)  ! Photosynthetically Active Radiation (W/m2)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_E_sf_w,par_E) ! Photosynthetically Active Radiation (umol/m2/s)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_swr_sf_w,  swr_J) ! Total shortwave radiation (W/m2) [up to 4000 nm]
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uv_sf_w, uv_J)    ! UV (W/m2)

      _HORIZONTAL_LOOP_END_

   end subroutine do_surface

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

   real(rk) function zenith_angle(days, hour, dlon, dlat)
       real(rk), intent(in)  :: days, hour
       real(rk), intent(in)  :: dlon, dlat

       real(rk), parameter       :: yrdays = 365.24_rk
       real(rk)                  :: th0, th02, th03, sundec
       real(rk)                  :: thsun, coszen
       real(rk)                  :: rlon, rlat

       ! from now on everything in radians
       rlon = deg2rad * dlon
       rlat = deg2rad * dlat

       ! Sun declination from Fourier expansion of Spencer (1971, Search 2:172).
       th0 = 2._rk*pi*days/yrdays
       th02 = 2._rk*th0
       th03 = 3._rk*th0
       sundec =  0.006918_rk - 0.399912_rk*cos(th0) + 0.070257_rk*sin(th0) &
               - 0.006758_rk*cos(th02) + 0.000907*sin(th02)                &
               - 0.002697_rk*cos(th03) + 0.001480*sin(th03)

       ! sun hour angle :
       thsun = (hour-12._rk)*15._rk*deg2rad + rlon

       ! cosine of the solar zenith angle :
       coszen = sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)

       zenith_angle = acos(coszen)
   end function zenith_angle

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

   subroutine navy_aerol_model(AM, WM, W, RH, V, costheta, alpha, beta, F_a, omega_a)
      ! AM: air-mass type (1 = marine aerosol-dominated, 10 continental aerosol-dominated)
      ! WM: wind speed averaged over past 24 h (m s-1)
      ! W: instantaneous wind speed (m s-1)
      ! RH: relative humidity (-)
      ! V: visibility (m)
      ! costheta: cosine of zenith angle
      ! alpha: Angstrom exponent. NB optical thickness tau = beta * lambda**(-alpha)
      ! beta: scale factor for optical thickness
      ! F_a: forward scattering probability
      ! omega_a: single scattering albedo
      real(rk), intent(in) :: AM, WM, W, RH, V, costheta
      real(rk), intent(out) :: alpha, beta, F_a, omega_a

      real(rk), parameter :: R = 0.05_rk
      integer, parameter :: nsize = 3, gridsize = 3
      real(rk), parameter :: r_o(nsize) = (/0.03_rk, 0.24_rk, 2.0_rk/)
      real(rk), parameter :: H_a = 1000._rk ! Aerosol scale height (m) Gregg & Carder 1990 p1665
      real(rk), parameter :: r_grid(gridsize) = (/0.1_rk, 1._rk, 10._rk/)

      real(rk) :: relhum, A(nsize), f, gamma, tau_a550
      real(rk) :: dNdr(gridsize), y(gridsize), x(gridsize)
      integer :: i
      real(rk) :: c_a550
      real(rk) :: B1, B2, B3, cos_theta_bar

      ! Impose upper limit on relative humidity (RH = 1 causes division by 0 in expression for f below)
      relhum = min(0.999_rk, RH)

      ! Amplitude functions for aerosol components (Eqs 21-23 Gregg & Carder 1990)
      ! Units are number of particles per cubic centimeter per micrometer (Gathman 1983)
      A(1) = 2000 * AM * AM
      A(2) = max(0.5_rk, 5.866_rk * (WM - 2.2_rk))
      A(3) = max(1.4e-5_rk, 0.01527_rk * (W - 2.2_rk) * R)

      ! function relating particle size to relative humidity (Eq 24 Gregg & Carder 1990)
      f = ((2._rk - relhum) / (6 * (1._rk - relhum)))**(1._rk / 3._rk)

      ! Particle density at different radii (Gregg & Carder 1990 p 1665)
      do i = 1, gridsize
          dNdr(i) = sum(A * exp(-log(r_grid(i) / f / r_o)**2) / f)
      end do

      ! Assume Junge distribution (i.e., particle density is power law of radius): dN/dr = C r^gamma
      ! Estimate gamma with least squares.
      ! Note: least-squares estimate of slope is x-y covariance / variance of x
      ! Due to the choice of r_grid (0.1, 1, 10), the sum of x_i = log10 r_i is 0.
      ! As a result, we do not need to subtract x_mean*y_mean and x_mean^2 to compute (co)variances.
      x = log10(r_grid)
      y = log10(dNdr)
      gamma = sum(x * y) / sum(x**2)

      ! Calculate Angstrom exponent from exponent of Junge distribution (Eq 26 Gregg & Carder 1990).
      ! Junge distribution: dN/d(ln r) = r dN/dr = C r^(-v)
      ! Thus, dN/dr = C r^(-v-1)
      ! The relation to gamma defined above: gamma = -v - 1. Thus, v = -gamma - 1
      ! The relationship between Junge distribution and Angstrom exponent is alpha = v - 2 (e.g. Tomasi et al. 1983)
      ! Thus, alpha = -gamma - 3
      alpha = -(gamma + 3)

      ! Estimate concentration parameter (Eqs 28, 29 Gregg & Carder 1990)
      c_a550 = 3.91_rk / V
      tau_a550 = c_a550 * H_a
      beta = tau_a550 * 550._rk**alpha

      ! Asymmetry parameter (Eq 35 Gregg & Carder 1990) - called alpha in Gregg & Casey 2009
      cos_theta_bar = -0.1417_rk * min(max(0._rk, alpha), 1.2_rk) + 0.82_rk

      ! Forward scattering probability (Eqs 31-34 Gregg & Carder 1990)
      ! NB B1-B3 are A, B, C in Gregg & Casey 2009 Eqs 3-6
      ! NB B1-B3 are AFS, BFS, ALG in Bird & Riordan 1986 Eqs 22-26
      B3 = log(1._rk - cos_theta_bar)
      B1 = B3 * (1.4590_rk + B3 * ( 0.1595_rk + 0.4129_rk * B3))
      B2 = B3 * (0.0783_rk + B3 * (-0.3824_rk - 0.5874_rk * B3))
      F_a = 1._rk - 0.5_rk * exp((B1 + B2 * costheta) * costheta)

      ! Single scattering albedo (Eq 36 Gregg & Carder 1990)
      omega_a = (-0.0032_rk * AM + 0.972_rk) * exp(3.06e-2_rk * relhum)
   end subroutine

   subroutine reflectance(nlambda, F, theta, W, rho_d, rho_s)
      integer, intent(in) :: nlambda
      real(rk), intent(in) :: F(nlambda), theta, W
      real(rk), intent(out) :: rho_d(nlambda), rho_s(nlambda)

      ! Coefficients for foam reflectance-wind speed relationship (p1666 Gregg and Carder 1990)
      real(rk), parameter :: D1 = 2.2e-5_rk
      real(rk), parameter :: D2 = 4.0e-4_rk
      real(rk), parameter :: D3 = 4.5e-5_rk
      real(rk), parameter :: D4 = 4.0e-5_rk

      ! Air density (g/m3)
      real(rk), parameter :: rho_a = 1.2e3_rk

      ! Refractive index of seawater, Gregg and Carder 1990 p1667
      real(rk), parameter :: n_w = 1.341_rk

      real(rk) :: C_D, tau, rho_f_W, rho_f(nlambda)
      real(rk) :: theta_r, b, rho_dsp, rho_ssp

      ! Drag coefficient (Eqs 42, 43 Gregg and Carder 1990)
      ! Constants match Trenberth et al. 1989 p1508 J Clim. However, they use different wind speed thresholds
      ! and a constant value between 3 and 10 m s-1. Theirs is also continuous, whereas the expression below is discontinuous.
      if (W <= 0._rk) then
         C_D = 0._rk
      else if (W <= 7._rk) then
         C_D = (0.62_rk + 1.56_rk / W) * 1.e-3_rk
      else
         C_D = (0.49_rk + 0.065_rk * W) * 1.e-3_rk
      end if

      ! Wavelength-independent foam reflectance [affects direct and diffuse light] (Eqs 39-41 Gregg and Carder 1990)
      ! Reformulated to make dependence on surface stress (tau, units seem to be 10-3 m2/s2) explicit.
      tau = rho_a * C_D * W**2
      if (W <= 4._rk) then
         rho_f_W = 0._rk
      elseif (W <= 7._rk) then
         rho_f_W = D1 * tau - D2
      else
         rho_f_W = D3 * tau - D4 * W**2
      end if

      ! Final wavelength and wind speed-dependent foam reflectance
      rho_f = rho_f_W * F

      ! Calculate zenith angle (radians) inside the water, taking refraction into account: Snell's law
      theta_r = asin(sin(theta) / n_w)

      ! Direct light specular component
      if (theta * rad2deg >= 40._rk .and. W > 2._rk) then
         ! Wind speed dependant (Eqs 46, 47 Gregg and Carder 1990)
         b = -7.14e-4_rk * W + 0.0618_rk
         rho_dsp =  0.0253_rk * exp(b * (theta * rad2deg - 40._rk))
      else
         ! Fresnel's Law, Eq 44 Gregg and Carder 1990 - note: contains typo (internal 1/2), see Kirk 3rd ed 2011, p 46
         rho_dsp = 0.5_rk * (sin(theta - theta_r)**2 / sin(theta + theta_r)**2 + tan(theta - theta_r)**2 / tan(theta + theta_r)**2)
      end if 

      ! Diffuse light specular component (p 1667 Gregg and Carder 1990)
      if (W > 4._rk) then
         rho_ssp = 0.057_rk
      else
         rho_ssp = 0.066_rk
      end if

      rho_d = rho_dsp + rho_f
      rho_s = rho_ssp + rho_f
   end subroutine

   subroutine slingo(mu0, LWP, r_e, nlambda, lambda, T_dcld, T_scld)
      ! mu0: cosine of zenith angle
      ! LWP: liquid water path (g m-2)
      ! r_e: equivalent radius of drop size distribution (um)
      real(rk), intent(in) :: mu0, LWP, r_e, lambda(nlambda)
      integer, intent(in) :: nlambda
      real(rk), intent(out), dimension(nlambda) :: T_dcld, T_scld

      real(rk), dimension(nlambda_slingo) :: tau, omega, g
      real(rk), dimension(nlambda_slingo) :: beta0, beta_mu0, f, U2
      real(rk), dimension(nlambda_slingo) :: alpha1, alpha2, alpha3, alpha4, epsilon, M, E, gamma1, gamma2
      real(rk), dimension(nlambda_slingo) :: one_minus_omega_f, gamma_denom, DIF_denom
      real(rk), dimension(nlambda_slingo) :: T_DB, R_DIF, T_DIF, R_DIR, T_DIR
      real(rk), parameter :: U1 = 7._rk / 4._rk

      ! Slingo 1989 Eqs 1-3
      ! cloud optical depth, single scatter albedo, asymmetry parameter
      ! derived for r_e = [4.2, 16.6]
      ! NB OASIM uses a default r_e that equals the mean of 10 um (Kiehl et al., 1998 J. Clim.) and 11.8 um (Han et al., 1994 J. Clim.)
      ! NB Stephens 1978 J Atmos Sciences Eq 10 links tau directly to LWP, which allows estimation of r_e as e.g. in Slingo 1989 section #4
      tau = LWP * (a_slingo + b_slingo / r_e)
      omega = 1._rk - (c_slingo + d_slingo * r_e)
      g = e_slingo + f_slingo * r_e

      ! fraction of scattered diffuse radation which is scattered into the backward hemisphere (Slingo 1989 Eq 6)
      beta0 = 3._rk / 7._rk * (1._rk - g)

      ! fraction of scattered direct radation which is scattered into the backward hemisphere (Slingo 1989 Eq 7)
      beta_mu0 = 0.5_rk - 0.75_rk * mu0 * g / (1._rk + g)

      ! fraction of scattered direct flux which emerges at zenith angles close to that of the incident beam (Slingo 1989 Eq 8)
      f = g * g

      ! reciprocals of the effective cosine for diffuse upward and downward fluxes (Slingo 1989 Eqs 9 and 10)
      ! JB 2018-09-12: Setting lower limit of U2 to 1, given that cosines cannot exceed 1 (thus its reciprocal cannot be lower than 1)      
      U2 = max(1._rk, U1 * (1._rk - (1._rk - omega) / 7._rk / omega / beta0))

      alpha1 = U1 * (1._rk - omega * (1._rk - beta0))
      alpha2 = U2 * omega * beta0
      alpha3 = (1 - f) * omega * beta_mu0
      alpha4 = (1 - f) * omega * (1 - beta_mu0)
      epsilon = sqrt(max(0._rk, alpha1 * alpha1 - alpha2 * alpha2))
      M = alpha2 / (alpha1 + epsilon)
      E = exp(-epsilon * tau)
      one_minus_omega_f = 1._rk - omega * f
      gamma_denom = one_minus_omega_f**2 - epsilon**2 * mu0**2
      gamma1 = ( one_minus_omega_f * alpha3 - mu0 * (alpha1 * alpha3 + alpha2 * alpha4)) / gamma_denom
      gamma2 = (-one_minus_omega_f * alpha4 - mu0 * (alpha1 * alpha4 + alpha2 * alpha3)) / gamma_denom

      ! Transmissivity for the direct solar beam (Slingo 1989 Eq 20)
      T_DB = exp(-one_minus_omega_f * tau / mu0)

      ! Diffuse reflectivity for diffuse incident radiation (Slingo 1989 Eq 21)
      DIF_denom = 1._rk - E**2 * M**2
      R_DIF = M * (1._rk - E**2) / DIF_denom

      ! Diffuse transmissivity for diffuse incident radiation (Slingo 1989 Eq 22)
      T_DIF = E * (1._rk - M**2) / DIF_denom

      ! Diffuse reflectivity for direct incident radiation (Slingo 1989 Eq 21)
      R_DIR = max(0._rk, -gamma2 * R_DIF - gamma1 * T_DB * T_DIF + gamma1)

      ! Diffuse transmissivity for direct incident radiation (Slingo 1989 Eq 21)
      T_DIR = min(1._rk, -gamma2 * T_DIF - gamma1 * T_DB * R_DIF + gamma2 * T_DB)

      call interp(nlambda_slingo, lambda_slingo, T_DB, nlambda, lambda, T_dcld)
      call interp(nlambda_slingo, lambda_slingo, T_DIR, nlambda, lambda, T_scld)
   end subroutine

end module
