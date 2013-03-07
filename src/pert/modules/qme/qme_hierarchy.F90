!
!
!
!
#include "util_allocation.h"

#define SUPERINDEX_FROM_K_L(k, l, lmax) ((l) + (lmax)*((k)-1))
#define K_FROM_SUPERINDEX(superindex, lmax) (((superindex) - 1) / (lmax) + 1)
#define L_FROM_SUPERINDEX(superindex, lmax) (mod(((superindex)-1) , (lmax)) + 1)

#define PALLOCATABLE pointer
!allocatable
!pointer
#define PALLOCATED   associated
!allocated
!associated
#define PNULLIFY(x) nullify(x)
! leave empty if not using pointers

!
!
!
!
module qme_hierarchy

	use prepare
	use twod

	use resources_qme
	use std_types
	use nakajima_zwanzig_shared

	use numer_ode
	use numer_fft
	use numer_matrix
	use sci_misc

	use util_allocation

	implicit none

 	! declarations

    private::write_phi_config_file

 	contains

    subroutine fill_evolution_superoperator_hierarchy(type)
        character, intent(in) :: type
        integer(i4b) :: i, j


    end subroutine fill_evolution_superoperator_hierarchy

    subroutine write_phi_config_file
        integer(i4b) :: i, j


    end subroutine write_phi_config_file

end module qme_hierarchy

