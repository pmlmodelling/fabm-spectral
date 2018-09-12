#include "fabm_driver.h"

module fabm_spectral

   ! A simple tracer model with support for sinking/floating and temperature-dependent decay.
   ! The temperature dependence of decay can be described with a Q10 formulation or an Arrhenius formulation.
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types

   implicit none

   private

   type,extends(type_base_model), public :: type_spectral
      type (type_horizontal_diagnostic_variable_id) :: id_swr, id_uv, id_par, id_par_E, id_O3, id_swr_sf, id_par_sf, id_uv_sf, id_par_E_sf
      type (type_horizontal_dependency_id) :: id_lon, id_lat, id_cloud, id_wind_speed, id_airpres
      type (type_global_dependency_id) :: id_yearday
      type (type_dependency_id) :: id_h
      integer :: nlambda
      type (type_horizontal_diagnostic_variable_id), dimension(:), allocatable :: id_surface_band_dir, id_surface_band_dif
      type (type_diagnostic_variable_id), dimension(:), allocatable :: id_band_dir, id_band_dif
      real(rk), dimension(:), allocatable :: lambda, par_weights, swr_weights, uv_weights
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

   integer, parameter  :: nlambda_bird = 39   ! Number of wavelengths looped over in atmosphere (Bird) code
   real(rk), parameter :: lambda_bird(nlambda_bird) = &
    (/300._rk, 305._rk, 310._rk, 315._rk, 320._rk, 325._rk, 330._rk, 335._rk, &   ! Table 1,  Bird Riordan 1986
      340._rk, 345._rk, 350._rk, 360._rk, 370._rk, 380._rk, 390._rk,  &                                       
      400._rk, 410._rk, 420._rk, 430._rk, 440._rk, 450._rk, 460._rk, 470._rk,  &                  
      480._rk, 490._rk, 500._rk, 510._rk, 520._rk, 530._rk, 540._rk, 550._rk,  &
      570._rk, 593._rk, 610._rk, 630._rk, 656._rk, 667.6_rk, 690._rk, 710._rk/)

   real(rk), parameter :: pi = 3.14159265358979323846_rk
   integer :: i500   ! Index of first wavelength > 500 nm
   real(rk) :: c(nlambda_bird, 7)

   interface interp
      module procedure interp_0d
      module procedure interp_1d
      module procedure interp_0d_scalar
      module procedure interp_1d_scalar
   end interface

contains

   subroutine initialize_module()
      ! Diffuse correction factor values Table 5,  Bird 1984
      real(rk), parameter :: c_ori(7, 7) = reshape( &
       (/0.88_rk, 1.08_rk, 1.11_rk, 1.04_rk, 1.15_rk, 1.12_rk, 1.32_rk,  &
         0.93_rk, 1.07_rk, 1.13_rk, 1.05_rk, 1.00_rk, 0.96_rk, 1.12_rk,  &
         1.02_rk, 1.11_rk, 1.18_rk, 1.09_rk, 1.00_rk, 0.96_rk, 1.07_rk,  &
         1.23_rk, 1.19_rk, 1.24_rk, 1.11_rk, 0.99_rk, 0.94_rk, 1.02_rk,  &
         2._rk,   1.51_rk, 1.46_rk, 1.24_rk, 1.06_rk, 0.99_rk, 1.10_rk,  &
         4._rk,   1.97_rk, 1.70_rk, 1.34_rk, 1.07_rk, 0.96_rk, 0.90_rk,  &
         6.3_rk,  3.76_rk, 2.61_rk, 1.72_rk, 1.22_rk, 1.04_rk, 0.80_rk/), (/7, 7/), order=(/2, 1/))

      ! Wavelengths at which diffuse correction factors are given (Table 5,  Bird 1984)
      real(rk), parameter :: clambda(7) = (/300._rk, 350._rk, 400._rk, 450._rk, 500._rk, 550._rk, 710._rk/)   ! Wavelengths at which diffuse correction factor is known; Table 5, Bird 1984

      ! Find index for first wavelength > 500 nm
      do i500 = 1, nlambda_bird
         if (lambda_bird(i500) > 500._rk) exit
      end do

      ! Interpolate diffuse correction factors to the wavelengths that we will resolve.
      call interp(7, 7, clambda, c_ori, nlambda_bird, lambda_bird, c)

      ! Transpose diffuse correction factor matrix to have wavelength dimension first.
      ! This enables runtime interpolation in second dimension (zenith angle).
      c = reshape(c, (/nlambda_bird, 7/), order=(/2, 1/))
   end subroutine initialize_module

   subroutine initialize(self,configunit)
      class (type_spectral), intent(inout), target :: self
      integer,               intent(in)            :: configunit

      integer :: l
      real(rk) :: lambda_min, lambda_max
      logical :: save_spectra
      character(len=8) :: strwavelength

      call initialize_module()

      call self%get_parameter(self%nlambda, 'nlambda', '', 'number of wavelengths')
      call self%get_parameter(lambda_min, 'lambda_min', 'nm', 'minimum wavelength')
      call self%get_parameter(lambda_max, 'lambda_max', 'nm', 'maximum wavelength')
      allocate(self%lambda(self%nlambda))
      do l = 1, self%nlambda
         self%lambda(l) = lambda_min + (l - 1) * (lambda_max - lambda_min) / (self%nlambda - 1)
      end do
      call self%get_parameter(save_spectra, 'save_spectra', '', 'save full spectral irradiance', default=.false.)

      ! Find wavelength bounds of photosynthetically active radiation
      allocate(self%par_weights(self%nlambda))
      allocate(self%swr_weights(self%nlambda))
      allocate(self%uv_weights(self%nlambda))
      call calculate_integral_weights(400._rk, 700._rk, self%nlambda, self%lambda, self%par_weights)
      call calculate_integral_weights(0._rk, 4000._rk, self%nlambda, self%lambda, self%swr_weights)
      call calculate_integral_weights(0._rk, 400._rk, self%nlambda, self%lambda, self%uv_weights)

      call self%register_dependency(self%id_lon, standard_variables%longitude)
      call self%register_dependency(self%id_lat, standard_variables%latitude)
      call self%register_dependency(self%id_wind_speed, standard_variables%wind_speed)
      call self%register_dependency(self%id_cloud, standard_variables%cloud_area_fraction)
      call self%register_dependency(self%id_airpres, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)

      call self%register_diagnostic_variable(self%id_swr, 'swr',     'W/m^2',     'shortwave radiation')
      call self%register_diagnostic_variable(self%id_uv,  'uv',      'W/m^2',     'ultraviolet radiative flux')
      call self%register_diagnostic_variable(self%id_par, 'par',     'W/m^2',     'photosynthetically active radiation')
      call self%register_diagnostic_variable(self%id_par_E, 'par_E', 'umol/m^2/s','photosynthetic photon flux density')
      call self%register_diagnostic_variable(self%id_swr_sf, 'swr_sf',     'W/m^2',     'net shortwave radiation in air')
      call self%register_diagnostic_variable(self%id_uv_sf,  'uv_sf',      'W/m^2',     'net ultraviolet radiative flux in air')
      call self%register_diagnostic_variable(self%id_par_sf, 'par_sf',     'W/m^2',     'net photosynthetically active radiation in air')
      call self%register_diagnostic_variable(self%id_par_E_sf, 'par_E_sf', 'umol/m^2/s','net photosynthetic photon flux density in air')
      call self%register_diagnostic_variable(self%id_O3, 'O3', 'atm cm', 'atmospheric ozone concentration')

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
            call self%register_diagnostic_variable(self%id_surface_band_dir(l), 'irradiance_dir_sf_' // trim(strwavelength), 'W/m2/nm', 'net direct irradiance in air @ ' // trim(strwavelength) // ' nm')
            call self%register_diagnostic_variable(self%id_surface_band_dif(l), 'irradiance_dif_sf_' // trim(strwavelength), 'W/m2/nm', 'net diffuse irradiance in air @ ' // trim(strwavelength) // ' nm')
            call self%register_diagnostic_variable(self%id_band_dir(l), 'irradiance_dir_' // trim(strwavelength), 'W/m2/nm', 'direct irradiance @ ' // trim(strwavelength) // ' nm')
            call self%register_diagnostic_variable(self%id_band_dif(l), 'irradiance_dif_' // trim(strwavelength), 'W/m2/nm', 'diffuse irradiance @ ' // trim(strwavelength) // ' nm')
         end do
      end if
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_spectral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk), parameter :: airden = 1.2e3_rk   ! Air density (g/m3), Text p1666 GregCard90
      real(rk), parameter :: d1 = 2.2e-5_rk, d2 = 4.0e-4_rk, d3 = 4.5e-5_rk, d4 = 4.0e-5_rk   ! Coefficients given by Text p1666 GregCard90

      real(rk), parameter :: W0 = 0.928_rk  ! Single scattering albedo of the aerosol (Bird 1984 Eq 15)
      real(rk), parameter :: fa = 0.82_rk   ! Forward to total scattering ratio of the aerosol (Bird 1984 Eq 15)
      real(rk), parameter :: pres0 = 1013._rk ! reference air pressure in mbar (Bird 1984 p461, Bird & Riordan p89)
      real(rk), parameter :: ga = 0.05_rk ! ground albedo

      ! Exter_rk,  mean irradiance at the top of the atmosphere (i.e. includes correction for Earth-sun distance and eccentricity; Table 1_rk,  Bird 1984)
      real(rk), parameter :: exter(nlambda_bird) = (/535.9_rk, 558.3_rk, 622._rk, 692.7_rk, 715.1_rk, 832.9_rk, 961.9_rk, 931.9_rk, 900.6_rk, 911.3_rk, 975.5_rk, 975.9_rk, 1119.9_rk, 1103.8_rk, 1033.8_rk,  &
         1479.1_rk, 1701.3_rk, 1740.4_rk, 1587.2_rk, 1837.0_rk, 2005.0_rk, 2043.0_rk,  &            
         1987.0_rk, 2027.0_rk, 1896.0_rk, 1909.0_rk, 1927.0_rk, 1831.0_rk, 1891.0_rk,  &
         1898.0_rk, 1892.0_rk, 1840.0_rk, 1768.0_rk, 1728.0_rk, 1658.0_rk, 1524.0_rk,  &
         1531.0_rk, 1420.0_rk, 1399.0_rk/)

      ! Zenith angles at which diffuse correction factor is known; Table 5, Bird 1984
      real(rk), parameter :: cd(7) = (/0._rk, 37._rk, 48.19_rk, 60._rk, 70._rk, 75._rk, 80._rk/)

      real(rk), parameter :: Planck = 6.62606957e-34     ! Planck constant (m2 kg/s)
      real(rk), parameter :: lightspeed = 299792458_rk   ! Speed of light (m/s)
      real(rk), parameter :: Avogadro = 6.02214129e23_rk ! Avogadro constant (/mol)
      
      real(rk) :: longitude, latitude, yearday, cloud_cover, wind_speed, pres
      real(rk) :: days, hour, zen, sunbet
      real(rk) :: toz(nlambda_bird), tw(nlambda_bird), ta(nlambda_bird), tr(nlambda_bird), tu(nlambda_bird)
      integer :: l
      real(rk) :: zend, airmass, ros(nlambda_bird), dir(nlambda_bird), &
                  c2(nlambda_bird), xx, Ir, Ia, Ig, O3, &
                  dif(nlambda_bird), foam, cdrag, b, dirspec, difspec, zenw
      real(rk), dimension(self%nlambda) :: direct, diffuse, spectrum  ! Spectra at top of the water column (with refraction and reflection accounted for)
      real(rk) :: par_J, swr_J, uv_J, par_E

      _HORIZONTAL_LOOP_BEGIN_
          _GET_HORIZONTAL_(self%id_lon, longitude)
          _GET_HORIZONTAL_(self%id_lat, latitude)
          _GET_GLOBAL_(self%id_yearday, yearday)
          _GET_HORIZONTAL_(self%id_cloud, cloud_cover)      ! Cloud cover (fraction, 0-1)
          _GET_HORIZONTAL_(self%id_wind_speed, wind_speed)  ! Wind speed @ 10 m above surface (m/s)
          _GET_HORIZONTAL_(self%id_airpres, pres)           ! Surface air pressure (Pa)

          days = floor(yearday) + 1.0_rk
          hour = mod(yearday, 1.0_rk) * 24.0_rk

          ! Calculate zenith angle and solar noon altitude (both in radians)
          call solar_angles(days, hour, longitude, latitude, zen, sunbet)

          zen = min(zen, 0.5_rk * pi)  ! Restrict the input zenith angle between 0 and pi/2
          zend = zen * 180._rk / pi 

          ! Bird 1984 p 466 reports O3=0.344 atm-cm. We use Van Heuklon (1979) to account for seasonal and geographical variation.
          O3 = estimate_ozone(longitude, latitude, yearday) * 0.001_rk
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O3, O3)

          ! Compute transmittance at relative air mass of 1.9
          call calculate_atmospheric_transmittance(1.9_rk, zen, O3, toz, tw, ta, tr, tu)         ! Following transmissivities sufficient for near-UV and visible wavelengths (Gregg and Carder, 1990)
          
          ! Calculate rho_s (air albedo)
          do l = 1, nlambda_bird
             ros(l) = toz(l) * tu(l) * tw(l) * (ta(l) * (1._rk - tr(l)) * 0.5_rk + tr(l) * (1._rk - ta(l)) * 0.22_rk * W0)      ! Eq (15) Bird 1984
          end do

          ! Calculate relative air mass (atmospheric path length)
          airmass = 1._rk / (cos(zen) + 0.15_rk * (93.885_rk - zend)**(-1.253_rk))          ! Eq 3 Bird 1984, Eq 5 Bird & Riordan 1986, Eq 13 Gregg & Carder 1990 note this is not pressure corrected.  
          if (airmass < 1._rk) airmass = 1._rk
          airmass = airmass * pres/100/pres0 ! NB converting surface pressure "pres" from Pa to mbar

          call calculate_atmospheric_transmittance(airmass, zen, O3, toz, tw, ta, tr, tu)

          ! Total Direct Normal Irradiance
          do L = 1, nlambda_bird
             dir(l) = exter(l) * tr(l) * ta(l) * tw(l) * toz(l) * tu(l)         ! Eq (1) Bird 1984, Calculating total direct normal irradiance from extraterrestrial irradiance and transmissivities at 24 wavelengths.
          end do

          ! Diffuse Irradiance

          ! Interpolating diffuse correction factor for the zenith angle of interest
          call interp(nlambda_bird, 7, cd, c, zend, c2)

          ! Calculating total diffuse irradiance (following Bird 1984: assumes independent scattering which results in underestimation of scattered irradiance at Z>60°, Bird and Riordan 1986)             
          do l = 1, nlambda_bird
             xx = exter(l) * cos(zen) * toz(l) * tw(l) * tu(l)                    
             Ir = xx * ta(l) * (1._rk - tr(l)) * 0.5_rk                   ! Eq (12) Bird 1984
             Ia = xx * tr(l) * (1._rk - ta(l)) * W0 * fa                  ! Eq (13) Bird 1984

             ! Removed Eq (14) Bird 1984: set to zero in Greg and Carder 90 text p1661 (makes no difference to rms)
             Ig = (dir(l) * cos(zen) + (Ir + Ia) * c2(l)) * ros(l) * ga / (1.D0 - ga * ros(l))   ! Eq (14) Bird 1984

             dif(l) = (Ir + Ia) * c2(l) + Ig                                 ! Eq (11) Bird 1984
          end do   

          ! Interpolate irradiance at Bird wavelengths to desired wavelenths
          call interp(nlambda_bird, lambda_bird, dir, self%nlambda, self%lambda, direct)
          call interp(nlambda_bird, lambda_bird, dif, self%nlambda, self%lambda, diffuse)

          ! exter in W/m2/um - divide by 1000 to convert to W/m2/nm  
          direct  = direct * 1e-3_rk
          diffuse = diffuse * 1e-3_rk

          ! Atmosphere-Ocean Refraction 

          ! Calculate zenith angle (radians) inside the water, taking refraction into account
          zenw = asin(sin(zen) / 1.341_rk)              ! Refractive index of seawater = 1.341, Gregg and Carder 1990

          ! Sea-Surface Reflection

          ! Setting the drag coefficient
          if (wind_speed <= 0._rk) then
             cdrag = 0._rk
          else if (wind_speed <= 7._rk)   then            
             cdrag = (0.62_rk + 1.56_rk / wind_speed) * 1.e-3_rk    ! Eq 42 Gregg and Carder 1990
          else                                     
             cdrag = (0.49_rk + 0.065_rk * wind_speed) * 1.e-3_rk   ! Eq 43 Gregg and Carder 1990
          end if

          ! Foam component: direct and diffuse light
          if (wind_speed <= 4._rk) then
             foam = 0._rk                                            ! Eq 39 Gregg and Carder 1990
          elseif (wind_speed <= 7._rk) then
             foam = d1 * airden * cdrag * wind_speed**2 - d2         ! Eq 40 Gregg and Carder 1990
          else
             foam = (d3 * airden * cdrag - d4) * wind_speed**2       ! Eq 41 Gregg and Carder 1990
          end if

          ! Direct light specular component
          if (zend >= 40._rk .and. wind_speed > 2._rk) then
             b = -7.14e-4_rk * wind_speed + 0.0618_rk              ! Eq 47 Gregg and Carder 1990
             dirspec =  0.0253_rk * exp(b * (zend - 40._rk))       ! Eq 46 Gregg and Carder 1990
          else
             ! Reflectance according to Fresnel's Law
             dirspec = 0.5_rk * (sin(zen - zenw)**2 / sin(zen + zenw)**2 + tan(zen - zenw)**2 / tan(zen + zenw)**2) ! Eq 44 Gregg and Carder 1990 - note: contains typo (internal 1/2), see Kirk 3rd ed 2011, p 46
          end if 

          ! Diffuse light specular component
          if (wind_speed > 4._rk) then               ! Text p 1667 GregCard90   
             difspec = 0.057_rk                     
          else
             difspec = 0.066_rk                  
          end if

          ! Apply reflection
          do l = 1, self%nlambda
             direct(l)  = direct(l) * (1._rk - (dirspec + foam))
             diffuse(l) = diffuse(l) * (1._rk - (difspec + foam))
          end do

          if (allocated(self%id_surface_band_dir)) then
             do l = 1, self%nlambda
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_band_dir(l), direct(l))
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_band_dif(l), diffuse(l))
             end do
          end if

          spectrum = direct + diffuse
          par_J = sum(self%par_weights * spectrum)
          swr_J = sum(self%swr_weights * spectrum)
          par_E = sum(self%par_weights * spectrum * self%lambda)
          uv_J  = sum(self%uv_weights * spectrum)
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_sf, par_J) ! Photosynthetically Active Radiation (W/m2)
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par_E_sf,par_E/(Planck*lightspeed)/Avogadro*1e-3_rk) ! Photosynthetically Active Radiation (umol/m2/s) - divide by 1e9 to go from nm to m, multiply by 1e6 to go from mol to umol
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_swr_sf,  swr_J) ! Total shortwave radiation (W/m2) [up to 4000 nm]
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uv, uv_J)       ! UV (W/m2)

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

   subroutine solar_angles(days,hour,dlon,dlat,zenith_angle,solar_noon_altitude)
       real(rk), intent(in)  :: days, hour
       real(rk), intent(in)  :: dlon, dlat
       real(rk), intent(out) :: zenith_angle, solar_noon_altitude

       real(rk), parameter       :: deg2rad = pi/180._rk
       real(rk), parameter       :: rad2deg = 180._rk/pi

       real(rk), parameter       :: yrdays = 365.24_rk
       real(rk)                  :: th0, th02, th03, sundec
       real(rk)                  :: thsun, coszen
       real(rk)                  :: rlon, rlat

       ! from now on everything in radians
       rlon = deg2rad*dlon
       rlat = deg2rad*dlat

       ! Sun declination from Fourier expansion of Spencer (1971, Search 2:172).
       th0 = 2._rk*pi*days/yrdays
       th02 = 2._rk*th0
       th03 = 3._rk*th0
       sundec =  0.006918_rk - 0.399912_rk*cos(th0) + 0.070257_rk*sin(th0) &
               - 0.006758_rk*cos(th02) + 0.000907*sin(th02)                &
               - 0.002697_rk*cos(th03) + 0.001480*sin(th03)

       ! sun hour angle :
       thsun = (hour-12._rk)*15._rk*deg2rad + rlon

       solar_noon_altitude = 0.5_rk*pi-rlat+sundec

       ! cosine of the solar zenith angle :
       coszen = sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)

       zenith_angle = acos(coszen)
   end subroutine solar_angles

   subroutine calculate_atmospheric_transmittance(airmass, zen, O3, toz, tw, ta, tr, tu)
      real(rk), intent(in) :: airmass, zen, O3
      real(rk), intent(out) :: toz(nlambda_bird), tr(nlambda_bird), ta(nlambda_bird), tw(nlambda_bird), tu(nlambda_bird)

      real(rk), parameter :: beta1  = 0.13238_rk  ! Alpha values from text p463, Bird 1984; beta values from Table 3, Bird 1984. Alpha/beta1 for 300–500 nm wavelengths, alpha/beta2 for 500–4000 nm wavelengths.
      real(rk), parameter :: alpha1 = 1.0274_rk
      real(rk), parameter :: beta2  = 0.116981_rk
      real(rk), parameter :: alpha2 = 1.206_rk
      real(rk), parameter :: w = 2._rk   ! Total precipitable water vapor (cm); from Table 5, Leckner 1978
      real(rk), parameter :: ho = 22._rk ! Height of maximum ozone concentration (km)
      !real(rk), parameter :: O3 = 0.344_rk ! Ozone amount (atm-cm), value taken from Bird 1984 p 466

      real(rk) :: mo
      integer :: l

      real(rk), parameter :: lampth(nlambda_bird) = (/0.3_rk, 0.305_rk, 0.31_rk, 0.315_rk, 0.32_rk, 0.325_rk, 0.33_rk, 0.335_rk, 0.34_rk, 0.345_rk, 0.35_rk, 0.36_rk, 0.37_rk, 0.38_rk, 0.39_rk,  &      ! Wavelengths in um; Table 1_rk,  Bird 1984
         0.4_rk,  0.41_rk,  0.42_rk, 0.43_rk, 0.44_rk,  0.45_rk,   0.46_rk, 0.47_rk,    &         
         0.48_rk, 0.49_rk,  0.5_rk,  0.51_rk, 0.52_rk,  0.53_rk,   0.54_rk, 0.55_rk,  &
         0.57_rk, 0.593_rk, 0.61_rk, 0.63_rk, 0.656_rk, 0.6676_rk, 0.69_rk, 0.71_rk/)

      real(rk), parameter :: AV(nlambda_bird) = (/0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,  &                     ! Water vapour absorption coefficients; Table 1_rk,  Bird Riordan 1986
         0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,  &      
         0.075_rk, 0._rk, 0._rk, 0._rk, 0._rk, 0.016_rk, 0.0125_rk/)

      real(rk), parameter :: AO(nlambda_bird) = (/10.0_rk, 4.80_rk, 2.70_rk, 1.35_rk, 0.8_rk, 0.38_rk, 0.16_rk, 0.075_rk, 0.04_rk, 0.019_rk, 0.007_rk, 0._rk, 0._rk, 0._rk, 0._rk,  &               ! Ozone absorption coefficients; Table 1_rk,  Bird Riordan 1984
         0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0.003_rk, 0.006_rk, 0.009_rk, 0.014_rk, 0.021_rk, 0.030_rk, 0.040_rk,  &   
         0.048_rk, 0.063_rk, 0.075_rk, 0.085_rk, 0.120_rk, 0.119_rk, 0.120_rk, 0.090_rk, 0.065_rk,  &
         0.051_rk, 0.028_rk, 0.018_rk/)

      real(rk), parameter :: AU(nlambda_bird) = (/0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,  &                     ! Water vapour absorption coefficients; Table 1_rk,  Bird Riordan 1986
         0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,  &      
         0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0.15_rk, 0._rk/)

      ! A) Raleigh scattering
      do l = 1, nlambda_bird
         tr(l) = exp(-airmass / (115.6406_rk * lampth(l)**4 - 1.335_rk * lampth(l)**2))   ! Eq 2 Bird 1984, Eq 4 Bird & Riordan 1986, Eq 15 Gregg & Carder 1990
      end do

      ! B) Aerosol scattering and absorption
      !  Turbidity is assumed to be 0.27 (a function of lambda, eq (4) of Bird 1984)
      do l = 1, i500 - 1                           
         ta(l) = exp(-beta1 * lampth(l)**(-alpha1) * airmass) ! Eq (6), Bird 1984. For 400–<500 nm, using alpha1 and beta1 turbidiy coefficient values
      end do
      do l = i500, nlambda_bird                               
         ta(l) = exp(-beta2 * lampth(l)**(-alpha2) * airmass) ! Eq (6), Bird 1984. For 500–710 nm, using alpha2 and beta2 turbidiy coefficient values
      end do

      ! C) Water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
      do l = 1, nlambda_bird 
         tw(l) = exp((-0.2385_rk * av(l) * w * airmass) / (1._rk + 20.07_rk * av(l) * w * airmass)**0.45_rk) ! Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
      end do

      ! D) Ozone and uniformly mixed gas absorption

      ! Ozone
      !mo = 35._rk/sqrt(1224._rk*cos(zen)**2 + 1._rk)               ! Eq 9, Bird 1984 
      mo = (1._rk + ho / 6370._rk) / sqrt(cos(zen)**2 + 2 * ho / 6370._rk)  ! Eq 10, Bird & Riordan 1986, Eq 14 Gregg & Carder 1990; NB 6370 is the earth's radius in km
      do l = 1, nlambda_bird
         toz(l) = exp(-ao(l) * O3 * mo)                             ! Eq 8 Bird 1984, Eq 9 Bird & Riordan 1986, Eq 17 Gregg & Carder 1990
      end do

      ! Uniformly mixed gas - SHOULD use pressure corrected airmass
      do L = 1, nlambda_bird
         tu(L) = exp(-1.41_rk * au(L) * airmass / (1._rk + 118.3_rk * au(L) * airmass)**0.45_rk)      ! Eq 10 Bird 1984, Eq 11, Bird and Riordan 1986, Eq 18 Gregg & Carder 1990        
      end do

   end subroutine calculate_atmospheric_transmittance

   function estimate_ozone(longitude, latitude, yearday) result(O3)
      ! Estimate ozone concentration based on Van Heuklon (1979)
      real(rk), intent(in) :: longitude, latitude, yearday
      real(rk) :: O3
      
      real(rk), parameter       :: deg2rad = pi/180._rk

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

   subroutine navy_aerols_model(AM, WM, W, RH, V, M, theta, nlambda, lambda, T_a)
      ! AM: air-mass type (1 = marine aerosol-dominated, 10 continental aerosol-dominated)
      ! WM: wind speed averaged over past 24 h (m s-1)
      ! W: instantaneous wind speed (m s-1)
      ! RH: relative humidity (-)
      ! V: visibility (m)
      ! M: air mass
      ! theta: zenith angle (radians)
      ! lambda: wave lengths (nm)
      real(rk), intent(in) :: AM, WM, W, RH, V, M, theta
      integer, intent(in) :: nlambda
      real(rk), intent(in) :: lambda(nlambda)
      real(rk), intent(out) :: T_a(nlambda)

      real(rk), parameter :: R = 0.05_rk
      integer, parameter :: nr = 3, nr_eval = 3
      real(rk), parameter :: r_o(nr) = (/0.03_rk, 0.24_rk, 2.0_rk/)
      real(rk), parameter :: H_a = 1000._rk ! Aerosol scale height (m) Gregg & Carder 1990 p1665
      real(rk), parameter :: r_eval(nr_eval) = (/0.1_rk, 1._rk, 10._rk/)

      real(rk) :: A(nr), f, gamma, alpha, beta
      real(rk) :: y(nr_eval), x(nr_eval), ymean, xmean
      real(rk) :: tau_alpha(nlambda)
      integer :: i
      real(rk) :: c_a550, tau_a550
      real(rk) :: B1, B2, B3, cos_theta_bar, F_a, omega_a

      ! Amplitude functions for aerosol components (Eqs 21-23 Gregg & Carder 1990)
      A(1) = 2000 * AM * AM
      A(2) = max(0.5_rk, 5.866_rk * (WM - 2.2_rk))
      A(3) = max(1.4e-5_rk, 0.01527_rk * (W - 2.2_rk) * R)

      ! function relating particle size to relative humidity (Eq 24 Gregg & Carder 1990)
      f = ((2._rk - RH) / (6 * (1._rk - RH)))**(1._rk / 3._rk)

      ! Estimate gamma with least squares
      do i = 1, nr_eval
          y(i) = log(sum(A*exp(-log(r_eval(i)/f/r_o)**2)/f))
      end do
      x = log(r_eval)
      xmean = sum(x) / nr_eval
      ymean = sum(y) / nr_eval

      ! Slope is x-y covariance / variance of x
      gamma = sum((y - ymean) * (x - xmean)) / sum((x - xmean)**2)

      ! Angstrom exponent (Eq 26 Gregg & Carder 1990)
      alpha = -(gamma + 3)

      ! Estimate concentration parameter (Eqs 28, 29 Gregg & Carder 1990)
      c_a550 = 3.91_rk / V
      tau_a550 = c_a550 * H_a
      beta = tau_a550 / 550._rk**(-alpha)

      ! Extinction coefficient
      tau_alpha = beta * lambda**(-alpha)

      ! Transmittance (Eq 26 Gregg & Carder 1990)
      T_a = exp(-tau_alpha * M)

      ! Asymmetry parameter (Eq 35 Gregg & Carder 1990) - called alpha in Gregg & Casey 2009
      cos_theta_bar = -0.1417_rk * min(max(0._rk, alpha), 1.2_rk) + 0.82_rk

      ! Forward scattering probability (NB B1-B3 are A-C in Gregg & Casey 2009 Eqs 3-6)
      B3 = log(1._rk - cos_theta_bar)
      B1 = B3 * (1.459_rk  + B3 * ( 0.1595_rk + 0.4129_rk * B3))
      B2 = B3 * (0.0783_rk + B3 * (-0.3824_rk - 0.5874_rk * B3))
      F_a = 1._rk - 0.5_rk * exp(B1 + B2 * cos(theta) * cos(theta))

      ! Single scattering albedo (Eq 36 Gregg & Carder 1990)
      omega_a = (-0.0032_rk * AM + 0.972_rk) * exp(3.06e-2_rk * RH)
   end subroutine
end module
