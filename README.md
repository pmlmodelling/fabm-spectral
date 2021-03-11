This is a [FABM](https://fabm.net) implementation of the Ocean-Atmosphere Spectral Irradiance Model (OASIM), which simulates the propagation of spectrally resolved direct and diffuse irradiance in atmosphere and oceans. The spectral range of the model is configurable and can include any wavelength between 200 and 4000 nm. Thus it can include total shortwave, ultraviolet, and photosynthetically active radiation, among others.

This implementation is designed for *online coupling* to hydrodynamic-biogeochemical models through FABM. Both atmosphere and ocean components run online as part of the coupled simulation and are directly forced by bulk meteorological properties. The modelled underwater irradiance can drive biogeochemical processes and, in turn, biogeochemical variables can affect (spectrally resolved) light absorption and scattering. The model also calculates light absorption per water layer, which can be fed back to the physical model as a heating term.

## How to use

Through FABM, this code can be used in combination with many hydrodynamic and biogeochemical models. For more information about use with specific hydrodynamic models, please visit [the FABM wiki](https://fabm.net/wiki).

When building FABM, you will need explicitly configure it to include the spectral model. To do this, provide the following arguments to cmake:

`-DFABM_INSTITUTES="<OTHER_MODELS>;spectral" -DFABM_SPECTRAL_BASE=<SPECTRALDIR>`.

Here, `<OTHER_MODELS>` is a semi-colon-separated list of additional FABM packages you want to build, e.g., ersem.
NB the double quotes are _required_ to prevent your shell (e.g., bash) from interpreting the semi-colon as "end of command".
`<SPECTRALDIR>` is the root of the spectral code repository (i.e., the directory where this README file is found).

This code comes with several test cases in the `testcases` directory. Most of these combine the spectral module with the [ERSEM](https://ersem.com) ecosystem model.

## Implementation

An early version of this codebase was developed in 2013 based upon the work of [Bird (1984)](https://doi.org/10.1016/0038-092X(84)90260-3), [Bird & Riordan (1986)](https://doi.org/10.1175/1520-0450(1986)025<0087:SSSMFD>2.0.CO;2) and [Gregg & Carder (1990)](https:/doi.org/10.4319/lo.1990.35.8.1657). The implementation was completely revised in 2019-2021 and is now based on:

* Gregg, W. W., & Casey, N. W. (2009). Skill assessment of a spectral ocean-atmosphere radiative model. Journal of Marine Systems, 76(1–2), 49–63. doi: [10.1016/j.jmarsys.2008.05.007](https://doi.org/10.1016/j.jmarsys.2008.05.007).
* Gregg, W. W., & Rousseaux, C. S. (2016). Directional and Spectral Irradiance in Ocean Models: Effects on Simulated Global Phytoplankton, Nutrients, and Primary Production. Frontiers in Marine Science, 3. doi: [10.3389/fmars.2016.00240](https://doi.org/10.3389/fmars.2016.00240).
* Gregg, W. W., & Rousseaux, C. S. (2017). Simulating PACE Global Ocean Radiances. Frontiers in Marine Science, 4. doi: [10.3389/fmars.2017.00060](https://doi.org/10.3389/fmars.2017.00060).

The ocean component does not currently resolve upwelling irradiance and therefore is simpler than the version described by Gregg & Rousseaux.

## Publications

This work was described and used in:

* Skákala, J., Bruggeman, J., Brewin, R. J. W., Ford, D. A., & Ciavatta, S. (2020). Improved Representation of Underwater Light Field and Its Impact on Ecosystem Dynamics: A Study in the North Sea. Journal of Geophysical Research: Oceans, 125(7). doi: [10.1029/2020JC016122](https://doi.org/10.1029/2020JC016122)