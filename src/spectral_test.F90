module spectral_test
   use iso_c_binding, only: c_float
   use slingo_clouds
   use fabm_types, only: rk

   implicit none

   contains

   subroutine run_spectral_test(mu0, LWP, r_e) bind(c)
      !DEC$ ATTRIBUTES DLLEXPORT :: run_spectral_test
      real(c_float), value :: mu0, LWP, r_e

      real(rk), dimension(nslingo) :: T_DB, T_DIF, T_DIR
      integer :: l

      call slingo(real(mu0, rk), real(LWP, rk), real(r_e, rk), T_DB, T_DIF, T_DIR)
      do l = 1, nslingo
         write (*,*) T_DB(l), T_DIF(l), T_DIR(l)
      end do
   end subroutine

end module