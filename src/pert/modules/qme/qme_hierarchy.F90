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

 	character(len=64), parameter :: external_dir = "external", config_filename = "config.prm", &
 	                    out_filename = 'hierarchy_out.dat'

    private::write_phi_config_file
    private::run_phi

 	contains

    subroutine fill_evolution_superoperator_hierarchy(type)
        character, intent(in) :: type
        integer(i4b) :: i, j

        if(type == 'O') then
            do i=1, N1_from_type('O')
                call write_phi_config_file(i,1,type)
                call run_phi()
            end do
        else if(type == 'E') then
            do i=1,N1_from_type('E')
            do j=1,N1_from_type('E')
                call write_phi_config_file(i,j,type)
                call run_phi()
            end do
            end do
        end if
    end subroutine fill_evolution_superoperator_hierarchy

    subroutine write_phi_config_file(i0, j0, type)
        integer(i4b), intent(in)       :: i0, j0
        character, intent(in)          :: type

        integer(i4b)                   :: i, j
        character(len=256)             :: buff, buff2

        if(type /= 'O' .and. type /= 'E') then
            call print_error_message(-1, 'unknown type in write_phi_config_file')
        end if

        open(unit=32,file=trim(trim(out_dir)//trim('/../')//trim(external_dir)//'/'//trim(config_filename) ) , err=32)

        if(.false.)then
32          write(buff,*) trim('Creating directory '//trim(external_dir))
            call print_log_message(adjustl(trim(buff)), 5)

            write(buff,*) trim(trim('mkdir')//' '//trim(out_dir)//trim('/../')//trim(external_dir) )
            call system(trim(buff))
            open(unit=32,file=trim(trim(out_dir)//trim('/../')//trim(external_dir)//'/'//trim(config_filename) ) , err=42)
        end if

        write(buff,'(A,I1,A)') '(A,I', int(ceiling(0.01 + log(real( N1_from_type('E') + 1 ))/log(10.0))) ,')'
        !write(buff2,'(A,I1,A)') '(A,I', int(ceiling(0.01 + log(real( N1_from_type('E') + 1 ))/log(10.0))) ,')'

        write(32,'(A)') '# Number of sites'
        write(32,trim(buff)) 'NumStates=',          N1_from_type('E')+1
        write(32,'(A)') '# Number of baths'
        write(32,trim(buff)) 'NumCouplingTerms=',   N1_from_type('E')+1
        write(32,'(A)') 'HierarchyTruncation=4'
        write(32,'(A)') 'MatsubaraTerms=1'
        write(32,'(A)') 'OutputFile='//adjustl(trim(out_filename))
        write(32,'(A)') '#T in Kelvin'
        write(buff2,'(F10.3)') temp
        write(32,'(A)') trim('Temperature='//adjustl(trim(buff2)))
        write(32,'(A)') '# Truncation scheme (default:1) Time local truncation'
        write(32,'(A)') '#TimeLocal=0'
        write(32,'(A,I5)') '# H in [cm^-1], using RWA ', int(rwa*Energy_internal_to_cm)

        write(32,'(A)') 'Hamiltonian:'

        ! row with optical coherences and ground state
        write(32,'(A)', advance='no') '0'
        do i=1,N1_from_type('E')
            write(32,'(A)', advance='no') ', 0'
        end do
        write(32,'(A)') ' '

        do i=1, N1_from_type('E')
            write(32,'(A)', advance='no') '0' ! for optical coherence
        do j=1, N1_from_type('E')
            if(i /= j) then
                write(buff2,'(F10.3)') iblocks(1,1)%sblock%J(i,j)*Energy_internal_to_cm
                write(32,'(A)', advance='no') trim(', ' // adjustl(trim(buff2)))
            else
                write(buff2,'(F10.3)') (iblocks(1,1)%sblock%en(i) - rwa)*Energy_internal_to_cm
                write(32,'(A)', advance='no') trim(', ' // adjustl(trim(buff2)))
            end if
        end do
            write(32,'(A)') ' '
        end do


        write(32,'(A)') 'InitialDensityMatrix:'
        do j=1,N1_from_type('E')+1

        write(32,'(A)', advance='no') '0'
        do i=1,N1_from_type('E')
            if((type == 'O' .and. i0 == i .and. j == 1) .or. &
                type == 'E' .and. i0 == i .and. j0 == j-1) then
                write(32,'(A)', advance='no') ', 1'
            else
                write(32,'(A)', advance='no') ', 0'
            end if
        end do
        write(32,'(A)') ' '

        end do


        write(32,'(A)') '# cutoff frequencies in ps^-1]'
        write(32,'(A)') 'gamma:'
        write(32,'(A)', advance='no') '10'
        do j=1, N1_from_type('E')
            write(buff2,'(F10.3)') 10.0
            write(32,'(A)', advance='no') trim(', ' // adjustl(trim(buff2)))
        end do
        write(32,'(A)') ' '



        write(32,'(A)') '# reorganization energies in cm^-1'
        write(32,'(A)') 'lambda:'

        write(32,'(A)', advance='no') '1' ! for ground state lambda
        do j=1, N1_from_type('E')
            ! for a lack of fantasy, I use reorganization energies defined by
            ! imaginary part of first integral of correlation functions, this is,
            ! of course wrong and should be replaced by a better solution
            write(buff2,'(F10.3)')  -aimag(all_hoft(iblocks(1,1)%sblock%gindex(j))%gg( &
                          size(all_hoft(iblocks(1,1)%sblock%gindex(j))%gg))           &
                         ) * Energy_internal_to_cm
            write(32,'(A)', advance='no') trim(', ' // adjustl(trim(buff2)))
        end do
        write(32,'(A)') ' '


        write(32,'(A)') '# Timesteps in ps'
        write(buff2,'(F10.3)') dt
        write(32,'(A)') trim('Timestep='//adjustl(trim(buff2)))
        write(32,'(A)') '# Run time in ps'
        write(buff2,'(F10.3)') dt*grid_Nt
        write(32,'(A)') trim('Time='//adjustl(trim(buff2)))

    close(32)


    return

42  write(*,*) trim(trim(out_dir)//trim('/../')//trim(external_dir)//'/'//trim(config_filename) ), ' could not be created, exiting...'
    stop

    end subroutine write_phi_config_file

    subroutine run_phi
        !call system('pwd')
        call system('phi '//trim(trim(out_dir)//trim('/../')//trim(external_dir)//'/'//trim(config_filename) )//' rk4 2')
    end subroutine run_phi

end module qme_hierarchy

