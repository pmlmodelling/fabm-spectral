module spectral_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use fabm_spectral
   use fabm_ozone
   use fabm_relhum

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: spectral_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('spectral'); allocate(type_spectral::model)
         case ('ozone'); allocate(type_ozone::model)
         case ('relhum'); allocate(type_relhum::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select
   end subroutine create

end module
