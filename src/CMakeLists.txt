add_library(fabm_models_spectral OBJECT
            spectral_model_library.F90
            spectral.F90
            ozone.F90
            relhum.F90
           )

target_include_directories(fabm_models_spectral PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(fabm_models_spectral fabm_base)

add_library(spectral_test SHARED EXCLUDE_FROM_ALL
            spectral_test.F90
            $<TARGET_OBJECTS:fabm_models_spectral>
           )
target_link_libraries(spectral_test fabm)
target_include_directories(spectral_test PRIVATE $<TARGET_PROPERTY:fabm_models_spectral,INTERFACE_INCLUDE_DIRECTORIES>)

#-----------------------------------------------------------------------
# Copyright (C) 2018 - PML
#-----------------------------------------------------------------------
