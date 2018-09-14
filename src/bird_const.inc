   integer, parameter  :: nlambda_bird = 39   ! Number of wavelengths looped over in atmosphere (Bird) code
   real(rk), parameter :: lambda_bird(nlambda_bird) = &
    (/300._rk, 305._rk, 310._rk, 315._rk, 320._rk, 325._rk, 330._rk, 335._rk, &   ! Table 1,  Bird Riordan 1986
      340._rk, 345._rk, 350._rk, 360._rk, 370._rk, 380._rk, 390._rk,  &                                       
      400._rk, 410._rk, 420._rk, 430._rk, 440._rk, 450._rk, 460._rk, 470._rk,  &                  
      480._rk, 490._rk, 500._rk, 510._rk, 520._rk, 530._rk, 540._rk, 550._rk,  &
      570._rk, 593._rk, 610._rk, 630._rk, 656._rk, 667.6_rk, 690._rk, 710._rk/)

      ! Exter_rk,  mean irradiance at the top of the atmosphere (i.e. includes correction for Earth-sun distance and eccentricity; Table 1_rk,  Bird 1984)
      real(rk), parameter :: exter(nlambda_bird) = (/535.9_rk, 558.3_rk, 622._rk, 692.7_rk, 715.1_rk, 832.9_rk, 961.9_rk, 931.9_rk, 900.6_rk, 911.3_rk, 975.5_rk, 975.9_rk, 1119.9_rk, 1103.8_rk, 1033.8_rk,  &
         1479.1_rk, 1701.3_rk, 1740.4_rk, 1587.2_rk, 1837.0_rk, 2005.0_rk, 2043.0_rk,  &            
         1987.0_rk, 2027.0_rk, 1896.0_rk, 1909.0_rk, 1927.0_rk, 1831.0_rk, 1891.0_rk,  &
         1898.0_rk, 1892.0_rk, 1840.0_rk, 1768.0_rk, 1728.0_rk, 1658.0_rk, 1524.0_rk,  &
         1531.0_rk, 1420.0_rk, 1399.0_rk/)

      ! Zenith angles at which diffuse correction factor is known; Table 5, Bird 1984
      real(rk), parameter :: cd(7) = (/0._rk, 37._rk, 48.19_rk, 60._rk, 70._rk, 75._rk, 80._rk/)
