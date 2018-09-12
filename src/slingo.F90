module slingo_clouds

   use fabm_types, only: rk

   implicit none

   private

   public :: slingo, nslingo

#include "slingo_const.inc"

contains
   
   subroutine slingo(mu0, LWP, r_e, T_DB, T_DIF, T_DIR, R_DIF, R_DIR)
      ! mu0: consine of zenith angle
      ! LWP: liquid water path (g m-2)
      ! r_e: equivalent radius of drop size distribution (um)
      real(rk), intent(in) :: mu0, LWP, r_e
      real(rk), intent(out), dimension(nslingo) :: T_DB, T_DIF, T_DIR, R_DIF, R_DIR

      real(rk), dimension(nslingo) :: tau, omega, g
      real(rk), dimension(nslingo) :: beta0, beta_mu0, f, U2
      real(rk), dimension(nslingo) :: alpha1, alpha2, alpha3, alpha4, epsilon, M, E, gamma1, gamma2
      real(rk), dimension(nslingo) :: one_minus_omega_f, gamma_denom, DIF_denom
      real(rk), parameter :: U1 = 7._rk / 4._rk

      ! Slingo 1989 Eqs 1-3
      ! cloud optical depth, single scatter albedo, asymmetry parameter
      ! derived for r_e = [4.2, 16.6]
      ! NB OASIM uses a default r_e that equals the mean of 10 um (Kiehl et al., 1998 J. Clim.) and 11.8 um (Han et al., 1994 J. Clim.)
      ! NB Stephens 1978 J Atmos Sciences Eq 10 links tau directly to LWP, which allows estimation of r_e as e.g. in Slingo 1989 section #4
      tau = LWP * (slingo_a + slingo_b / r_e)
      omega = 1._rk - (slingo_c + slingo_d * r_e)
      g = slingo_e + slingo_f * r_e

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
   end subroutine

   end module slingo_clouds