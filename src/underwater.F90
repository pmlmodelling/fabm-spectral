#include "fabm_driver.h"

module fabm_spectral_underwater

   ! Spectral light model based on Bird 1984, Bird & Riordan 1986, Gregg &  Carder 1990, Gregg & Casey 2009
   ! Copyright PML 2018

   use fabm_types
   use fabm_standard_variables
   use fabm_expressions, only: temporal_mean
   use fabm_particle

   use fabm_spectral_shared

   implicit none

   private

   type type_iop
      type (type_dependency_id) :: id_c         ! concentration (could be chl, carbon, or something else, but its product with a or b below should return units m-1)
      real(rk), dimension(:), allocatable :: a  ! specific absorption (m-1 concentration-1)
      real(rk), dimension(:), allocatable :: b  ! specific total scattering (m-1 concentration-1)
      real(rk) :: b_b                           ! ratio of backscattering to total scattering (dimensionless)
   end type

   type,extends(type_particle_model), public :: type_spectral_underwater
      type (type_diagnostic_variable_id) :: id_swr, id_uv, id_par, id_par_E, id_par_E_scalar, id_par_J_scalar, id_par_E_dif, id_swr_abs, id_secchi
      type (type_dependency_id) :: id_h
      type (type_surface_dependency_id), allocatable :: id_direct_sf(:), id_diffuse_sf(:)
      type (type_surface_dependency_id) :: id_costheta_r
      integer :: nlambda
      type (type_diagnostic_variable_id), dimension(:), allocatable :: id_band_dir, id_band_dif, id_a_band, id_b_band, id_Kd
      type (type_diagnostic_variable_id), dimension(:,:), allocatable :: id_a_iop
      real(rk), dimension(:), allocatable :: lambda, lambda_bounds, par_weights, par_E_weights, swr_weights, uv_weights, F, lambda_out

      real(rk), dimension(:), allocatable :: a_w, b_w, a_w_out, b_w_out
      type (type_iop), allocatable :: iops(:)
      integer :: l490_l
      integer :: spectral_output
      logical :: save_Kd
   contains
      procedure :: initialize
      procedure :: do_column
   end type

#include "oasim_const.inc"
#include "water_const.inc"
#include "phyto_const.inc"

contains

   subroutine initialize(self, configunit)
      class (type_spectral_underwater), intent(inout), target :: self
      integer,               intent(in)            :: configunit

      integer :: l
      integer :: lambda_method
      integer :: n_iop, i_iop, iop_type
      real(rk) :: lambda_min, lambda_max
      integer :: nlambda_out
      character(len=8) :: strwavelength, strindex, strindex2
      real(rk) :: lambda_ref_iop, a_star_iop, S_iop, b_star_iop, eta_iop

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

      call self%get_parameter(n_iop, 'n_iop', '', 'number of inherent optical properties (IOPs)', default=0)
      allocate(self%iops(n_iop))
      do i_iop = 1, n_iop
         allocate(self%iops(i_iop)%a(self%nlambda))
         allocate(self%iops(i_iop)%b(self%nlambda))
         write(strindex, '(i0)') i_iop
         call self%get_parameter(iop_type, 'iop'//trim(strindex)//'_type', '', 'type of IOP '//trim(strindex)//' (1: diatoms, 2: chlorophytes, 3: cyanobacteria, 4: coccolithophorids, 5: dinoflagellates, 6: detritus, 8: CDOC, 9: OM with custom absorption/scattering)', minimum=1, maximum=9)
         select case (iop_type)
         case (1) ! diatoms
            call interp(size(lambda_diatoms), lambda_diatoms, a_diatoms, self%nlambda, self%lambda, self%iops(i_iop)%a)
            call interp(size(lambda_diatoms), lambda_diatoms, b_diatoms, self%nlambda, self%lambda, self%iops(i_iop)%b)
            self%iops(i_iop)%b_b = 0.002 ! Gregg & Rousseau 2016 but originally Morel 1988
         case (2) ! chlorophytes
            call interp(size(lambda_chlorophytes), lambda_chlorophytes, a_chlorophytes, self%nlambda, self%lambda, self%iops(i_iop)%a)
            call interp(size(lambda_chlorophytes), lambda_chlorophytes, b_chlorophytes, self%nlambda, self%lambda, self%iops(i_iop)%b)
            self%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Ahn et al. 1992 as reported in Gregg & Rousseau 2016
         case (3) ! cyanobacteria
            call interp(size(lambda_cyanobacteria), lambda_cyanobacteria, a_cyanobacteria, self%nlambda, self%lambda, self%iops(i_iop)%a)
            call interp(size(lambda_cyanobacteria), lambda_cyanobacteria, b_cyanobacteria, self%nlambda, self%lambda, self%iops(i_iop)%b)
            self%iops(i_iop)%b_b = 0.0032 ! Gregg & Rousseau 2016 but originally Ahn et al. 1992
         case (4) ! coccolithophorids
            call interp(size(lambda_coccolithophores), lambda_coccolithophores, a_coccolithophores, self%nlambda, self%lambda, self%iops(i_iop)%a)
            call interp(size(lambda_coccolithophores), lambda_coccolithophores, b_coccolithophores, self%nlambda, self%lambda, self%iops(i_iop)%b)
            self%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Morel 1988 as reported in Gregg & Rousseau 2016
         case (5) ! dinoflagellates
            call interp(size(lambda_dinoflagellates), lambda_dinoflagellates, a_dinoflagellates, self%nlambda, self%lambda, self%iops(i_iop)%a)
            call interp(size(lambda_dinoflagellates), lambda_dinoflagellates, b_dinoflagellates, self%nlambda, self%lambda, self%iops(i_iop)%b)
            self%iops(i_iop)%b_b = 0.0029 ! Gregg & Rousseau 2016 but originally Morel 1988
         case (6) ! detritus (small, as described in Gregg & Rousseau 2016)
            ! Parameters below match the small organic detritus parametrization of Gallegos et al. 2011 (table 2) - the latter also offers a parametrizaton for large detritus.
            ! Note: Gallegos et al. 2011 constants are specific to dry weight! Did Gregg & Rousseau 2016 misinterpret them as specific to carbon weight?
            ! If so: Babin et al 2003 state 2.6 g DW per g C is representative for suspended OM
            ! NB 12.0107 converts from mg-1 to mmol-1
            ! JB 14/1/2019: adding factor 2.6 (DW/C) discussed above as that seems like to make L4/WCO match much better.
            self%iops(i_iop)%a(:) = 8e-5_rk * exp(-0.013_rk * (self%lambda - 440_rk)) * 12.0107_rk * 2.6_rk
            self%iops(i_iop)%b(:) = 0.00115_rk * (550._rk / self%lambda)**0.5_rk * 12.0107_rk * 2.6_rk
            self%iops(i_iop)%b_b = 0.005_rk
         !case (7) ! PIC
         !   self%iops(i_iop)%a(:) = 0
         !   call interp(size(), lambda_w, a_w, self%nlambda, self%lambda, self%iops(i_iop)%b)
         !   self%iops(i_iop)%b_b = 0.01_rk
         case (8) ! CDOC
            ! NB 12.0107 converts from mg-1 to mmol-1
            self%iops(i_iop)%a(:) = 2.98e-4_rk * exp(-0.014_rk * (self%lambda - 443_rk)) * 12.0107_rk
            self%iops(i_iop)%b(:) = 0
            self%iops(i_iop)%b_b = 0
         case (9) ! Custom carbon-specific absorption and total scattering spectra
            ! NB 12.0107 converts from mg-1 to mmol-1
            call self%get_parameter(a_star_iop, 'a_star_iop'//trim(strindex), 'm2/mg C', 'carbon-mass-specific absorption coefficient for IOP '//trim(strindex)//' at reference wavelength', minimum=0._rk, default=0._rk)
            if (a_star_iop /= 0._rk) then
               call self%get_parameter(lambda_ref_iop, 'lambda_a_iop'//trim(strindex), 'nm', 'reference wavelength for absorption by IOP '//trim(strindex))
               call self%get_parameter(S_iop, 'S_iop'//trim(strindex), '-', 'exponent of absorption spectrum for IOP '//trim(strindex), minimum=0._rk)
               self%iops(i_iop)%a(:) = a_star_iop * exp(-S_iop * (self%lambda - lambda_ref_iop)) * 12.0107_rk
            else
               self%iops(i_iop)%a(:) = 0
            end if

            call self%get_parameter(b_star_iop, 'b_star_iop'//trim(strindex), 'm2/mg C', 'carbon-mass-specific scattering coefficient for IOP '//trim(strindex)//' at reference wavelength', minimum=0._rk, default=0._rk)
            if (b_star_iop /= 0._rk) then
               call self%get_parameter(lambda_ref_iop, 'lambda_b_iop'//trim(strindex), 'nm', 'reference wavelength for scattering by IOP '//trim(strindex))
               call self%get_parameter(eta_iop, 'eta_iop'//trim(strindex), '-', 'exponent of scattering spectrum for IOP '//trim(strindex), minimum=0._rk)
               call self%get_parameter(self%iops(i_iop)%b_b, 'b_b_iop'//trim(strindex), '-', 'backscattering-to-total-scattering ratio for IOP '//trim(strindex), minimum=0._rk)
               self%iops(i_iop)%b(:) = b_star_iop * (lambda_ref_iop / self%lambda)**eta_iop * 12.0107_rk
            else
               self%iops(i_iop)%b(:) = 0
               self%iops(i_iop)%b_b = 0
            end if
         end select

         ! Protect against negative coefficients caused by extrapolation beyond source spectrum boundaries.
         self%iops(i_iop)%a(:) = max(self%iops(i_iop)%a, 0._rk)
         self%iops(i_iop)%b(:) = max(self%iops(i_iop)%b, 0._rk)

         ! Link to concentration metric to allow us to convert *specific* absorption/scattering into actual absorption and scattering (in m-1)
         if (iop_type >=1 .and. iop_type <= 5) then
            ! Phytoplankton: chlorophyll-specific absorption and scattering
            call self%register_dependency(self%iops(i_iop)%id_c, 'iop' // trim(strindex) // '_chl', 'mg Chl m-3', 'chlorophyll in IOP ' // trim(strindex))
            call self%request_coupling_to_model(self%iops(i_iop)%id_c, 'iop' // trim(strindex), type_bulk_standard_variable(name='total_chlorophyll'))
         else
            ! POM/DOM/PIC: carbon-specific absorption and scattering
            call self%register_dependency(self%iops(i_iop)%id_c, 'iop' // trim(strindex) // '_c', 'mmol C m-3', 'carbon in IOP ' // trim(strindex))
            call self%request_coupling_to_model(self%iops(i_iop)%id_c, 'iop' // trim(strindex), standard_variables%total_carbon)
         end if
      end do

      ! Find wavelength bounds of photosynthetically active radiation
      allocate(self%par_weights(self%nlambda))
      allocate(self%swr_weights(self%nlambda))
      allocate(self%uv_weights(self%nlambda))
      allocate(self%par_E_weights(self%nlambda))
      call calculate_integral_weights(400._rk, 700._rk, self%nlambda, self%lambda, self%par_weights)
      call calculate_integral_weights(300._rk, 4000._rk, self%nlambda, self%lambda, self%swr_weights)
      call calculate_integral_weights(300._rk, 400._rk, self%nlambda, self%lambda, self%uv_weights)
      self%par_E_weights(:) = self%par_weights * self%lambda /(Planck*lightspeed)/Avogadro*1e-3_rk ! divide by 1e9 to go from nm to m, multiply by 1e6 to go from mol to umol

      call self%register_dependency(self%id_h, standard_variables%cell_thickness)

      ! Scalar downwelling irradiance within the water column
      call self%register_diagnostic_variable(self%id_swr,     'swr',     'W/m^2',      'downwelling shortwave flux', standard_variable=standard_variables%downwelling_shortwave_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_uv,      'uv',      'W/m^2',      'downwelling ultraviolet radiative flux', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par,     'par',     'W/m^2',      'downwelling photosynthetic radiative flux', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_E,   'par_E',   'umol/m^2/s', 'downwelling photosynthetic photon flux', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_J_scalar,'par_J_scalar','W/m^2', 'scalar downwelling photosynthetic radiative flux', standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_E_scalar,'par_E_scalar','umol/m^2/s', 'scalar downwelling photosynthetic photon flux', source=source_do_column)
      call self%register_diagnostic_variable(self%id_swr_abs, 'swr_abs', 'W/m^2',      'absorption of shortwave energy in layer', standard_variable=standard_variables%net_rate_of_absorption_of_shortwave_energy_in_layer, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_E_dif, 'par_E_dif', 'W/m^2',      'diffusive downwelling photosynthetic photon flux', source=source_do_column)
      !call self%register_diagnostic_variable(self%id_secchi,  'secchi',  'm',          'Secchi depth (1.7/Kd 490)', standard_variable=standard_variables%secchi_depth, source=source_do_column)

      ! Interpolate absorption and scattering spectra to user wavelength grid
      allocate(self%a_w(self%nlambda), self%b_w(self%nlambda))
      call interp(nlambda_w, lambda_w, a_w, self%nlambda, self%lambda, self%a_w)
      call interp(nlambda_w, lambda_w, b_w, self%nlambda, self%lambda, self%b_w)

      call self%get_parameter(self%spectral_output, 'spectral_output', '', 'spectral output (0: none, 1: full, 2: selected wavelengths)', default=0, minimum=0, maximum=2)
      select case (self%spectral_output)
      case (1)
         allocate(self%lambda_out(self%nlambda))
         self%lambda_out(:) = self%lambda
      case (2)
         call self%get_parameter(nlambda_out, 'nlambda_out', '', 'number of wavebands for spectral output', minimum=1)
         allocate(self%lambda_out(nlambda_out))
         do l = 1, nlambda_out
            write(strindex, '(i0)') l
            call self%get_parameter(self%lambda_out(l), 'lambda' // trim(strindex) // '_out', '', 'output wavelength ' // trim(strindex))
         end do
         allocate(self%a_w_out(nlambda_out), self%b_w_out(nlambda_out))
         call interp(nlambda_w, lambda_w, a_w, nlambda_out, self%lambda_out, self%a_w_out)
         call interp(nlambda_w, lambda_w, b_w, nlambda_out, self%lambda_out, self%b_w_out)
      end select
      call self%get_parameter(self%save_Kd, 'save_Kd', '', 'compute attenuation', default=.false.)

      if (allocated(self%lambda_out)) then
         allocate(self%id_band_dir(size(self%lambda_out)))
         allocate(self%id_band_dif(size(self%lambda_out)))
         allocate(self%id_a_iop(size(self%lambda_out), size(self%iops)))
         allocate(self%id_a_band(size(self%lambda_out)))
         allocate(self%id_b_band(size(self%lambda_out)))
         if (self%save_Kd) allocate(self%id_Kd(size(self%lambda_out)))
         do l = 1, size(self%lambda_out)
            if (self%lambda_out(l) < 1000._rk) then
               write(strwavelength, '(f5.1)') self%lambda_out(l)
            else
               write(strwavelength, '(f6.1)') self%lambda_out(l)
            end if
            write(strindex, '(i0)') l
            !call self%register_diagnostic_variable(self%id_band_dir(l), 'dir_band' // trim(strindex), 'W/m2/nm', 'direct irradiance @ ' // trim(strwavelength) // ' nm', source=source_do_column)
            !call self%register_diagnostic_variable(self%id_band_dif(l), 'dif_band' // trim(strindex), 'W/m2/nm', 'diffuse irradiance @ ' // trim(strwavelength) // ' nm', source=source_do_column)
            if (self%save_Kd) call self%register_diagnostic_variable(self%id_Kd(l), 'Kd_band' // trim(strindex), 'm-1', 'attenuation @ ' // trim(strwavelength) // ' nm', source=source_do_column)
            do i_iop = 1, size(self%iops)
               write(strindex2, '(i0)') i_iop
               call self%register_diagnostic_variable(self%id_a_iop(l, i_iop), 'a_iop' // trim(strindex2) // '_band' // trim(strindex), '1/m', 'absorption by IOP ' // trim(strindex2) // ' @ ' // trim(strwavelength) // ' nm', source=source_do_column)
            end do
            call self%register_diagnostic_variable(self%id_a_band(l), 'a_band' // trim(strindex), 'm-1', 'total absorption excluding water @ ' // trim(strwavelength) // ' nm', source=source_do_column)
            call self%register_diagnostic_variable(self%id_b_band(l), 'b_band' // trim(strindex), 'm-1', 'total scattering excluding water @ ' // trim(strwavelength) // ' nm', source=source_do_column)
         end do
      end if
      
      allocate(self%id_diffuse_sf(self%nlambda))
      allocate(self%id_direct_sf(self%nlambda))
      do l = 1, self%nlambda
         if (self%lambda(l) < 1000._rk) then
            write(strwavelength, '(f5.1)') self%lambda(l)
         else
            write(strwavelength, '(f6.1)') self%lambda(l)
         end if
         write(strindex, '(i0)') l
         call self%register_dependency(self%id_direct_sf(l), 'direct_band' // trim(strindex), 'W/m2/nm', 'downward direct irradiance in water @ ' // trim(strwavelength) // ' nm')
         call self%register_dependency(self%id_diffuse_sf(l), 'diffuse_band' // trim(strindex), 'W/m2/nm', 'downward diffuse irradiance in water @ ' // trim(strwavelength) // ' nm')
         call self%request_coupling(self%id_direct_sf(l), type_surface_standard_variable('downward_direct_irradiance_in_water_at_' // trim(strwavelength) // '_nm', 'W m-2 nm-1'))
         call self%request_coupling(self%id_diffuse_sf(l), type_surface_standard_variable('downward_diffuse_irradiance_in_water_at_' // trim(strwavelength) // '_nm', 'W m-2 nm-1'))
      end do
      call self%register_dependency(self%id_costheta_r, 'costheta_r', '-', 'cosine of zenith angle of underwater direct irradiance')
   end subroutine initialize

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_spectral_underwater), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk), parameter :: mcosthetas = 0.831_rk        ! Mean of cosine of angle of diffuse radiation in water, assuming all angular contributions equal in air (Sathyendranath and Platt 1989, p 191) NB Ackleson et al 1994 use 0.9
      real(rk), parameter :: r_s = 1.5_rk                 ! Shape factor representing mean backscatter coefficient of diffuse irradiance (r_d in Ackleson et al. 1994 p 7487)

      integer :: l
      real(rk), dimension(self%nlambda) :: direct, diffuse, spectrum, Kd  ! Spectra at top of the water column (with refraction and reflection accounted for)
      real(rk) :: par_J, swr_J, uv_J, par_E

      real(rk), dimension(self%nlambda) :: a, b, b_b, a_iop
      real(rk), dimension(self%nlambda) :: f_att_d, f_att_s, f_prod_s
      real(rk), allocatable :: spectrum_out(:)
      integer :: i_iop
      real(rk) :: c_iop, h, swr_top, costheta_r, dir_frac

      if (self%spectral_output == 2) allocate(spectrum_out(size(self%lambda_out)))

      do l = 1, self%nlambda
         _GET_SURFACE_(self%id_direct_sf(l), direct(l))
         _GET_SURFACE_(self%id_diffuse_sf(l), diffuse(l))
      end do
      _GET_SURFACE_(self%id_costheta_r, costheta_r)
      spectrum = direct + diffuse
      swr_J = sum(self%swr_weights * spectrum)

      _DOWNWARD_LOOP_BEGIN_
         ! Save downwelling shortwave flux at top of the layer
         swr_top = swr_J

         ! Compute absorption, total scattering and backscattering in current layer from IOPs
         a = self%a_w
         b = self%b_w
         b_b = 0.5_rk * self%b_w
         do i_iop = 1, size(self%iops)
            _GET_(self%iops(i_iop)%id_c, c_iop)
            a_iop = c_iop * self%iops(i_iop)%a
            select case (self%spectral_output)
            case (1)
               do l = 1, self%nlambda
                  _SET_DIAGNOSTIC_(self%id_a_iop(l, i_iop), a_iop(l))
               end do
            case (2)
               call interp(self%nlambda, self%lambda, a_iop, size(self%lambda_out), self%lambda_out, spectrum_out)
               do l = 1, size(self%lambda_out)
                  _SET_DIAGNOSTIC_(self%id_a_iop(l, i_iop), spectrum_out(l))
               end do
            end select
            a = a + a_iop
            b = b + c_iop * self%iops(i_iop)%b
            b_b = b_b + c_iop * self%iops(i_iop)%b_b * self%iops(i_iop)%b
         end do

         ! Transmissivity of direct/diffuse attentuation and conversion from direct to diffuse - for one half of the layer
         _GET_(self%id_h, h)
         f_att_d = exp(-0.5_rk * (a + b) * h / costheta_r)         ! Gregg & Rousseau 2016 Eq 8
         f_att_s = exp(-0.5_rk * (a + r_s * b_b) * h / mcosthetas) ! Gregg & Rousseau 2016 Eq 9
         f_prod_s = exp(-0.5_rk * a * h / costheta_r) - f_att_d    ! Gregg & Rousseau 2016 Eq 14 but not accounting for backscattered fraction

         if (self%save_Kd .and. self%spectral_output /= 0) then
            do l = 1, self%nlambda
               !Kd(l) = 2 * (log(direct(l) + diffuse(l)) - log(direct(l) * (f_att_d(l) + f_prod_s(l)) + diffuse(l) * f_att_s(l))) / h
               if (direct(l) + diffuse(l) > 0) then
                  ! Direct and diffuse stream
                  dir_frac = direct(l) / (direct(l) + diffuse(l))
                  Kd(l) = - 2 *log(dir_frac * (f_att_d(l) + f_prod_s(l)) + (1.0_rk - dir_frac) * f_att_s(l)) / h
               else
                  ! Only diffuse stream
                  Kd(l) = (a(l) + r_s * b_b(l)) / mcosthetas
               end if
            end do
         end if

         ! From top to centre of layer
         direct = direct * f_att_d
         diffuse = diffuse * f_att_s + direct * f_prod_s
         spectrum = direct + diffuse

         par_J = sum(self%par_weights * spectrum)
         swr_J = sum(self%swr_weights * spectrum)
         par_E = sum(self%par_E_weights * spectrum)
         uv_J  = sum(self%uv_weights * spectrum)
         _SET_DIAGNOSTIC_(self%id_par, par_J)   ! Photosynthetically Active Radiation (W/m2)
         _SET_DIAGNOSTIC_(self%id_par_E, par_E) ! Photosynthetically Active Radiation (umol/m2/s)
         _SET_DIAGNOSTIC_(self%id_swr,  swr_J)  ! Total shortwave radiation (W/m2) [up to 4000 nm]
         _SET_DIAGNOSTIC_(self%id_uv, uv_J)     ! UV (W/m2)
         _SET_DIAGNOSTIC_(self%id_par_E_dif, sum(self%par_E_weights * diffuse)) ! Diffuse Photosynthetically Active photon flux (umol/m2/s)

         ! Compute scalar PAR as experienced by phytoplankton
         spectrum = direct / costheta_r + diffuse / mcosthetas
         par_J = sum(self%par_weights * spectrum)
         par_E = sum(self%par_E_weights * spectrum)
         _SET_DIAGNOSTIC_(self%id_par_J_scalar, par_J) ! Scalar Photosynthetically Active Radiation (W/m2)
         _SET_DIAGNOSTIC_(self%id_par_E_scalar, par_E) ! Scalar Photosynthetically Active photon flux (umol/m2/s)

         ! From centre to bottom of layer
         direct = direct * f_att_d
         diffuse = diffuse * f_att_s + direct * f_prod_s
         spectrum = direct + diffuse

         ! Save spectrally resolved outputs
         select case (self%spectral_output)
         case (1)
            do l = 1, self%nlambda
               _SET_DIAGNOSTIC_(self%id_a_band(l), a(l) - self%a_w(l))
               _SET_DIAGNOSTIC_(self%id_b_band(l), b(l) - self%b_w(l))
                if (self%save_Kd) _SET_DIAGNOSTIC_(self%id_Kd(l), Kd(l))
            end do
         case (2)
            call interp(self%nlambda, self%lambda, a, size(self%lambda_out), self%lambda_out, spectrum_out)
            do l = 1, size(self%lambda_out)
               _SET_DIAGNOSTIC_(self%id_a_band(l), spectrum_out(l) - self%a_w_out(l))
            end do
            call interp(self%nlambda, self%lambda, b, size(self%lambda_out), self%lambda_out, spectrum_out)
            do l = 1, size(self%lambda_out)
               _SET_DIAGNOSTIC_(self%id_b_band(l), spectrum_out(l) - self%b_w_out(l))
            end do
            if (self%save_Kd) then
               call interp(self%nlambda, self%lambda, Kd, size(self%lambda_out), self%lambda_out, spectrum_out)
               do l = 1, size(self%lambda_out)
                  _SET_DIAGNOSTIC_(self%id_Kd(l), spectrum_out(l))
               end do
            end if
         end select

         ! Compute remaining downwelling shortwave flux and from that, absorption [heating]
         swr_J = sum(self%swr_weights * spectrum)
         _SET_DIAGNOSTIC_(self%id_swr_abs, swr_top - swr_J)
         !_SET_DIAGNOSTIC_(self%id_secchi, 0._rk)
      _VERTICAL_LOOP_END_

      ! Put remaining shortwave in bottom layer
      ! (assumes all light is absorbed by sediment and injected into water column)
      _MOVE_TO_BOTTOM_
      _SET_DIAGNOSTIC_(self%id_swr_abs, swr_J)
   end subroutine do_column

end module
