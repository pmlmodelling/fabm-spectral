add_library(fabm_models_spectral OBJECT
            spectral_model_library.F90
            spectral.F90
            slingo.F90
           )

target_include_directories(fabm_models_spectral PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(fabm_models_spectral fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015 - PML
#-----------------------------------------------------------------------
