add_library(fabm_models_spectral OBJECT
            spectral_model_library.F90
            spectral.F90
           )

add_dependencies(fabm_models_spectral fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015 - PML
#-----------------------------------------------------------------------
