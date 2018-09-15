module spectral_test
   use iso_c_binding, only: c_double
   use fabm_spectral
   use fabm_types, only: rk

   implicit none

   contains

   subroutine run_spectral_test(mu0, LWP, r_e, T_DB, T_DIF, T_DIR, R_DIF, R_DIR) bind(c)
      !DEC$ ATTRIBUTES DLLEXPORT :: run_spectral_test
      real(c_double), value :: mu0, LWP, r_e
      real(c_double), dimension(nlambda_slingo) :: T_DB, T_DIF, T_DIR, R_DIF, R_DIR

      integer :: l

      call slingo(real(mu0, rk), real(LWP, rk), real(r_e, rk), nlambda_slingo, lambda_slingo, T_DB, T_DIR)
   end subroutine

end module