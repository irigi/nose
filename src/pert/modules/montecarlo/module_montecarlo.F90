!
!
!
!
#include "util_allocation.h"
!#include "omp_lib.h"

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
module module_montecarlo

    use prepare

    use std_types

!    use omp_lib
    use qch_lib
    use sci_misc

    use numer_ode
    use numer_fft
    use numer_matrix
    use sci_misc
    use resources

    use nakajima_zwanzig_shared

    use util_allocation
    use resources_montecarlo

    implicit none

     ! declarations
     integer(i4b), private    ::    TRAJECTORIES, RUNS, MICROGROUP

     integer(i4b), private    ::    STEPS = -1
     integer(i4b), private    ::    Nl                                 ! number of one-excitons
     real(dp), dimension(:,:), pointer, private   ::  J_coupl          ! pointer to couplings
     real(dp), private        ::  timeStep                             ! timestep of the method
     real(dp), private        ::  jumps_in_one_run                     ! jumps in average in one run
     real(dp), private        ::  p_of_no_jump, p_of_no_jump_measured  ! probability of no jump (important for restarting)
     complex(dpc), dimension(:,:,:), allocatable, private     :: rho, rho_coherent
     real(dp), parameter, private :: CoherentFactorMultiplierToLowerJumpProb = 1.0_dp

!     integer(i1b), dimension(:,:,:), allocatable :: trajectory_depository
!     integer(i4b), dimension(:,:), allocatable   :: factor_depository

    complex, dimension(:,:,:,:,:), allocatable   :: gg_MC, hh_MC, cc_MC

    complex, dimension(:), allocatable     :: MC_polar_1
    real, dimension(:), allocatable        :: MC_spect_abs

    logical, parameter :: only_coherences = .true.

    real(dp), private :: dom, oma
    integer(i4b), private :: NFFT, padfac
    integer, parameter, private :: NFFT_basic = 1024

    logical, private :: MC_polar_1_collected


    public::do_montecarlo_work
    public::collect_montecarlo_data
    public::init_monte_carlo
    public::clean_montecarlo

!    private::collect_rho
    private::init_random_seed

    private::jump_probability_144
    private::jump_probability_444
    private::jump_phase_change
    private::jump_probability_total_IJ
    private::jump_probability_total_

    private::generate_trajectory
    private::generate_trajectory_re
    private::generate_trajectory_cmplx

    private::goft_site
    private::hoft_site
    private::init_goft_general
    private::deinit_goft_general
    private::goft_general
    private::dgoft_general
    private::goft_general_2
    private::dgoft_general_2
    private::goft_general_4
    private::dgoft_general_4
    private::hoft_general
    private::Vh_general
    private::hV_general
    private::VV_general

    private::calculate_Gfactor_from_trajectory_history
    private::calculate_Gfactor_from_trajectory_history_general_basis
    private::calculate_Ifactor_from_trajectory_history
    private::next_step_of_trajectory
    private::number_of_trajectory_steps
    private::perform_montecarlo
    private::create_spect_abs

    private::get_coherent_dynamics
    private::seam_superoperators_U

    private::write_gofts

    interface goft_general
        module procedure goft_general_2
        module procedure goft_general_4
    end interface
    interface dgoft_general
        module procedure dgoft_general_2
        module procedure dgoft_general_4
    end interface

    interface generate_trajectory
        module procedure generate_trajectory_re
        module procedure generate_trajectory_cmplx
    end interface

    interface jump_probability
        module procedure jump_probability_144
        module procedure jump_probability_444
    end interface
    interface jump_probability_total
        module procedure jump_probability_total_
        module procedure jump_probability_total_IJ
    end interface

    contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_montecarlo_work(err)
        integer, intent(out)    :: err

        integer(i4b)            :: r,s,k,l,i, closest_i, closest_j
        character(len=128)      :: cbuff
        character :: type

        complex(dpc), dimension(:,:), allocatable     :: rho_in
        real(dp), dimension(:), allocatable           :: vector

        Nl = iblocks(1,1)%eblock%N1

        if(only_coherences) then
            type = 'O'
        else
            type = 'E'
        end if

        call init_monte_carlo()
        call init_goft_general()
        call write_gofts()

        if(debug_gamma > 1e-12) then
            write(cbuff, '(A F6.4 A)') "Using debugging Gamma = ",debug_gamma," - result is unphysical"
            call print_warning_message(cbuff, -1)
        end if
        if(vynechani_G_ifu) then
            write(cbuff, '(A F6.4 A)') "WARNING! FAST_G METHOD CONTAINS SERIOUS ERROR, USE ONLY FOR TESTING!"
            call print_warning_message(cbuff, -1)
        end if

        write(cbuff, '(A L)') "Optical coherences ", only_coherences
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "G-functions ", g_functions
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Exciton basis ", use_exciton_basis
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Fixed Seed ", fixed_seed
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Modified Unraveling 2 ", modified_unraveling2
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Depository ", depository
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "FastG ", vynechani_G_ifu
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Exciton basis unraveling ", exciton_basis_unraveling
        call print_log_message(cbuff, 5)

        if(load_evops) then

        !*************************************************************
        !  Reading evops if possible
        !*************************************************************
        if (only_coherences) then
            if (parallel_id == 0) then
                ALLOCATE(rho_in, (N1_from_type(type), N2_from_type(type)))
                ALLOCATE(vector, (Nl*Nl*2+1))

                if(.not. use_exciton_basis) then
                    open(unit=11,file=trim(file_join(out_dir,"Ueg_local.dat")), err=130)
                else
                    open(unit=11,file=trim(file_join(out_dir,"Ueg_exciton.dat")), err=130)
                end if

                write(cbuff,'(i3)') Nl*Nl*2+1

                cbuff = '('//trim(cbuff)//'f15.8)' ! forming a descriptor for outputting row by row

                do i = 1, RUNS*STEPS

                  read(11,cbuff, err=110, end=120) vector(:)
                  write(*,*) vector(:)
                  write(*,*)

                  do l=1,Nl
                  do k=1,Nl
                    rho_in(k,1) = vector(1+2*k-1 + 2*Nl*(l-1)) + vector(1+2*k + 2*Nl*(l-1))*cmplx(0,1)
                    write(*,*) rho_in(k,1)
                    write(*,*)
                  end do

!                    if(use_exciton_basis) then
!                        call operator_to_exc(rho_in,'O')
!                    end if

                    evops(1,1)%Ueg(l, 1, :, 1, i) = rho_in(:,1)

                  end do

                end do

                goto 100

                130 call print_error_message(-1, "error reading Ueg_local(_exciton).dat")
                110 close(unit=11) !call print_error_message(-1, "error reading Ueg_local.dat")
                120 close(unit=11)
                100 call print_log_message("evops read successfuly",5)

                DEALLOCATE(rho_in)
                DEALLOCATE(vector)
            end if
        end if

        else                            ! not load evops

        !*************************************************************
        ! Calculation of evolution superops
        !*************************************************************
        do r=N1_from_type(type),1,-1
        do s=1,N2_from_type(type)

            call perform_montecarlo(r,s)

            do k=1,N1_from_type(type)
            do l=1,N2_from_type(type)

            rho(k,l,1) = 0.0_dp
            rho(r,s,1) = 1.0_dp

            rho_coherent(k,l,1) = 0.0_dp
            rho_coherent(r,s,1) = 1.0_dp

            if(maxval(abs(rho(:,:,min(STEPS, size(rho,3)) ))) > 0.05) then
                call print_warning_message("Optical coherences do not die sufficiently, there will be artefacts in spectra, choose more steps!",-1)
            end if

            ! saving result into evolution superoperator
            do i=1,Nt(1)
                if(i > Nt(1) .or. i < 1) then
                    call print_error_message(-1, 'index error in fill_evolution_superoperator_nakajima_zwangig')
                end if

                closest_i = i*(dt*gt(1))/timeStep + 0.5_dp

                if(closest_i < 1) then
                    closest_i = 1
                end if

                if(closest_i <= size(rho,3)) then

                    if(only_coherences) then
                        evops(1,1)%Ueg(k,l, r,s, INT((i-1)) + 1) = rho(k,l,closest_i)
                    else
                        evops(1,1)%Uee(k,l, r,s, INT((i-1)) + 1) = rho(k,l,closest_i)
                    endif

                else
                    if(only_coherences) then
                        evops(1,1)%Ueg(k,l, r,s, INT((i-1)) + 1) = 0.0_dp
                    else
                        evops(1,1)%Uee(k,l, r,s, INT((i-1)) + 1) = 0.0_dp
                    endif
                end if

                !write(*,*) evops(1,1)%Ueg(k,l, r,s, INT((i-1)) + 1), closest_i

            end do

            if(only_coherences) then
                write(cbuff,'(a,i1,i1,i1,i1,a)') 'Ueg(',k,l,r,s,') filled'
            else
                write(cbuff,'(a,i1,i1,i1,i1,a)') 'Uee(',k,l,r,s,') filled'
            end if

            call print_log_message(trim(cbuff),5)

            end do
            end do

        end do
        end do

        ! seamlessly extend the superoperators
        do r=2, RUNS
            closest_i = (r-1)*STEPS*timeStep/(dt*gt(1)) + 1 + 0.5_dp
            closest_j = r*STEPS/(dt*gt(1))*timeStep + 1 + 0.5_dp
            closest_j = closest_j - 1

            write(*,*) 'seam indices', closest_i, closest_j

            if(only_coherences) then
                call seam_superoperators_U(closest_i, closest_j, 'O')
            end if
        end do

        ! transformation into excitonic picture
        if(use_exciton_basis .and. (.not. exciton_basis_unraveling) ) then
            do i=1,Nt(1)
                if(only_coherences) then
                    call superops_to_exc(evops(1,1)%Ueg(:,:, :,:, INT((i-1)) + 1), 'O')
                else
                    call superops_to_exc(evops(1,1)%Uee(:,:, :,:, INT((i-1)) + 1), 'E')
                endif

            end do
        end if

        end if                            !***************************************************


        !*************************************************************
        ! Output of linear absorption spectrum
        !*************************************************************
        if (resources_output_contains("spect_abs")) then

            ! output linear absorption spectrum

            call create_polar_1()
            call create_spect_abs()

            if (parallel_id == 0) then
                open(unit=11,file=trim(file_join(out_dir,"spect_abs.dat")))

                do i = 1, NFFT
                     oma = ((-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

                    write(11,*) oma, MC_spect_abs(i)

                end do

                close(unit=11)
            end if

        end if

        call deinit_goft_general()

    end subroutine do_montecarlo_work

    subroutine collect_montecarlo_data(err)
        integer, intent(out) :: err
        complex(dpc), dimension(:,:), allocatable :: rho_out, rho_out2
        character :: type

        integer   :: k,i,l, fs
        real(dp)  :: oma, sum

        character(len = 300)     :: buffer
        character(len = 300)     :: name
        character(len = 800)     :: caption_line ! with intention of proceding output data with Origin

        err = 0

        if(only_coherences) then
            type = 'O'
        else
            type = 'E'
        end if

        ALLOCATE(rho_out, (N1_from_type(type), N2_from_type(type)))
        ALLOCATE(rho_out2, (Nl, Nl))

!        call collect_rho()

        !*************************************************************
        !  Outputting populations
        !*************************************************************
        if (resources_output_contains("populations") .or. resources_output_contains("local_populations")) then

            ! output of the populations

            if (parallel_id == 0) then
                if(use_exciton_basis) then
                    open(unit=11,file=trim(file_join(out_dir,"populations.dat")))
                else
                    open(unit=11,file=trim(file_join(out_dir,"local_populations.dat")))
                end if

                !-------------------- Organizes the caption line for Origin
                if(origin) then
                    caption_line = 'time'
                    do i = 1, Nl
                        write(buffer,'(i1)') i
                        buffer = 'Population'//trim(buffer)//'_Re'//' Population'//trim(buffer)//'_Im '
                        caption_line = trim(caption_line)//' '//trim(buffer)
                      end do
                      write(11,'(a, /)') trim(caption_line)
                  end if

                !--------------------


                write(buffer,'(i3)') Nl+1
                buffer = '('//trim(buffer)//'f15.8)' ! forming a descriptor for outputting row by row

                do i = 1, RUNS*STEPS

                  if(only_coherences) then
                      cycle
                  end if

                  write(11,buffer,advance='no') (i-1)*timeStep

                  sum = 0
                  if(only_coherences) then
                      rho_out(:,1) = rho(:,1,i)
                  else
                      rho_out = rho(:,:,i)
                  end if
                  if(use_exciton_basis .and. (.not. exciton_basis_unraveling) ) then
                      call operator_to_exc(rho_out,type)
                  end if

                  do k=1,Nl
                      write(11,buffer,advance='no') rho_out(k,k)
                      sum = sum + rho(k,k,i)
                  end do
                  write(11,buffer,advance='no') sum

                  write(11,buffer)

                end do

                close(unit=11)
            end if

        end if

        !*************************************************************
        !  Outputting coherences
        !*************************************************************
        if (resources_output_contains("coherences") .or. resources_output_contains("local_coherences")) then

            ! output of the coherences

            if (parallel_id == 0) then
                if(use_exciton_basis) then
                    open(unit=11,file=trim(file_join(out_dir,"coherences.dat")))
                else
                    open(unit=11,file=trim(file_join(out_dir,"local_coherences.dat")))
                end if

                !-------------------- Organizes the caption line for Origin
                if(origin) then
                    caption_line = 'time'
                    do k = 1, Nl
                    do l = 1, Nl
                        if(k == l .and. .not. only_coherences .or. only_coherences .and. l > 1) then
                            cycle
                        end if

                        write(buffer,'(i1i1)') k,l
                        buffer = 'Coherence'//trim(buffer)//'_Re'//' Coherence'//trim(buffer)//'_Im '
                        caption_line = trim(caption_line)//' '//trim(buffer)
                      end do
                      end do
                      write(11,'(a, /)') trim(caption_line)
                  end if

                !--------------------


                if(debug_gamma > 1e-12 .and. only_coherences .and. Nl == 2) then
                    write(buffer,'(i3)') 2*Nl+1
                else
                    write(buffer,'(i3)') Nl+1
                end if
                buffer = '('//trim(buffer)//'f15.8)' ! forming a descriptor for outputting row by row

                do i = 1, RUNS*STEPS

                  write(11,buffer,advance='no') (i-1)*timeStep

                  if(only_coherences) then
                      rho_out(:,1) = rho(:,1,i)
                  else
                      rho_out = rho(:,:,i)
                  end if
                  if(use_exciton_basis .and. (.not. exciton_basis_unraveling) ) then
                      call operator_to_exc(rho_out,type)
                  end if

                  do k=1,Nl
                  do l=1,Nl
                      if(k == l .and. .not. only_coherences .or. only_coherences .and. l > 1) then
                            cycle
                    end if

                      if(debug_gamma > 1e-12 .and. only_coherences .and. Nl == 2) then
                          !write(11,buffer,advance='no') rho_out(k,l), get_dimer_gamma_evol(k, (i-1)*timeStep)
                          rho_out2 = 0.0
                          rho_out2(1,1) = 1.0
                          call get_coherent_dynamics(rho_out2,(i-1)*timeStep)
                          write(11,buffer,advance='no') rho_out(k,l), rho_out2(k,l)
                      else
                          write(11,buffer,advance='no') rho_out(k,l)
                      end if
                  end do
                  end do

                  write(11,buffer)

                end do

                close(unit=11)
            end if

        end if

        !*************************************************************
        !  Outputting dens_opt_coh
        !*************************************************************
        if(only_coherences) then
            rho_out(:,1) = rho(:,1,1)
            if (resources_output_contains(NOSE_RDM_B01)) then
                call write_time_evolutions('O',.false.,.false.,rho_out)
            endif
            if (resources_output_contains(NOSE_RDM_B01_ABS)) then
                call write_time_evolutions('O',.true.,.false.,rho_out)
            endif
            if (resources_output_contains(NOSE_RDM_B01_CONJG)) then
                call write_time_evolutions('O',.false.,.true.,rho_out)
            endif

            call write_evolution_operators('O')
            !call write_redfield_tensor('O')
        else
            call write_evolution_operators('E')
        end if

        !*************************************************************
        !  Outputting evolution operators (mandatory)
        !*************************************************************
        if (only_coherences .and. .not. load_evops) then
            if (parallel_id == 0) then

                if(.not. use_exciton_basis) then
                    open(unit=11,file=trim(file_join(out_dir,"Ueg_local.dat")))
                else
                    open(unit=11,file=trim(file_join(out_dir,"Ueg_exciton.dat")))
                end if

                write(buffer,'(i3)') Nl+1

                buffer = '('//trim(buffer)//'f15.8)' ! forming a descriptor for outputting row by row

                do i = 1, RUNS*STEPS

                  write(11,buffer,advance='no') (i-1)*timeStep

                  do l=1,Nl
                  do k=1,Nl

                    rho_out(:,1) = evops(1,1)%Ueg(l, 1, :, 1, i)
!                    if(use_exciton_basis) then
!                        call operator_from_exc(rho_out,'O')
!                    end if

                      write(11,buffer,advance='no') rho_out(k,1)
                  end do
                  end do

                  write(11,buffer)

                end do

                close(unit=11)
            end if
        end if

        DEALLOCATE(rho_out)
        DEALLOCATE(rho_out2)
    end subroutine collect_montecarlo_data

    subroutine init_monte_carlo()
!        if(abs(gt(1)-1) > 1e-3) then
!            call print_error_message(-1, "gt(1) is unequal to 1 -- this is not properly handled by Montecarlo module")
!        end if

        if ((resources_output_contains("coherences") .and. resources_output_contains("local_coherences")) .or. &
                (resources_output_contains("populations") .and. resources_output_contains("local_populations")) .or. &
                (resources_output_contains("coherences") .and. resources_output_contains("local_populations")) .or. &
                (resources_output_contains("populations") .and. resources_output_contains("local_coherences"))  ) then

            call print_error_message(-1, "I can't print both local and exciton output (choose coherences or local_coherences, not both).")
        else if(resources_output_contains("local_coherences") .or. resources_output_contains("local_populations")) then
            use_exciton_basis = .false.
        else
            use_exciton_basis = .true.
        end if

        if(fixed_seed) then
            call init_random_seed(0)
        else
            call init_random_seed()
        end if

        TRAJECTORIES = 10000*Nt(2)
!        TRAJECTORIES_STORED = min(TRAJECTORIES, 1000)
        RUNS = 1 !gt(2)
        jumps_in_one_run = gt(3)

        MICROGROUP = Nt(3)

        Nl = iblocks(1,1)%eblock%N1
        STEPS = Nt(1)
        J_coupl => iblocks(1,1)%sblock%J
        !Nl = size(J_coupl,1)

        if(jump_probability_total() < 1e-5) then
            timeStep = (dt*gt(1))/RUNS
        else
            timeStep = jumps_in_one_run/jump_probability_total()/STEPS

            if(timeStep > (dt*gt(1))/RUNS) then
                timeStep = (dt*gt(1))/RUNS
            else

            end if
        end if

        jumps_in_one_run = jump_probability_total()*timeStep*STEPS

        p_of_no_jump = (1.0_dp - jump_probability_total()*timeStep)**STEPS


        ALLOCATE(rho,(Nl,Nl,STEPS*RUNS))
        ALLOCATE(rho_coherent,(Nl,Nl,STEPS*RUNS))
        rho = 0.0_dp
        rho_coherent = 0.0_dp

!        ALLOCATE(trajectory_depository,(TRAJECTORIES_STORED,2,STEPS*RUNS))
!        trajectory_depository = 0


            ALLOCATE(evops(1,1)%Ueg,(Nl, 1,Nl, 1,Nt(1)) )

            ALLOCATE(evops(1,1)%Uee,(Nl,Nl,Nl,Nl,Nt(1)) )

            ALLOCATE(evops(1,1)%Ugg,(1,1,1,1,Nt(1)) )

        evops(1,1)%Ueg = 0.0_dp
        evops(1,1)%Uee = 0.0_dp
        evops(1,1)%Ugg = 1.0_dp

        if((dt*gt(1))/timeStep > 10) then
            call print_warning_message("TimeStep in configuration file too big, g(t) can have too rough step", -1)
        end if
    end subroutine init_monte_carlo

    subroutine clean_montecarlo()
        DEALLOCATE(rho)
        DEALLOCATE(rho_coherent)
!        DEALLOCATE(trajectory_depository)
    end subroutine clean_montecarlo

    subroutine perform_montecarlo(i0, j0)
        integer(i4b), intent(in) :: i0, j0

        integer(i4b) :: run, i, j, a, b, jj, m, i_m, aa, bb, debug_int
        logical      :: debug_bool
        real(dp)     :: r, norm

        integer(i1b), dimension(2,STEPS*RUNS, MICROGROUP) :: draha
!@!        integer(i1b), dimension(TRAJECTORIES_STORED,2,STEPS*RUNS) :: depository_tmp
        complex(dpc), dimension(STEPS*RUNS, MICROGROUP)   :: factor, Gfactor, Ifactor
        complex(dpc), dimension(Nl,Nl,STEPS*RUNS, MICROGROUP) :: rho_micro, rho_coherent_micro

        complex(dpc) :: factor_in
!@!        complex(dpc), dimension(Nl, Nl) :: rho_init
        character(len=256) :: buff
        integer(i4b), dimension(3) :: time

        rho = 0.0_dp
        rho_coherent = 0.0_dp

        factor_in = 1.0_dp
        Gfactor = cmplx(1,0,dpc)

        write(buff,'(f12.3)') jump_probability_total()*timeStep
        buff = 'Averagely ' // trim(buff) // ' jumps in one run'
        call print_log_message(trim(buff),5)

        do run=1, RUNS

            call itime(time)
            write(buff, '(A I2 A I2.2 A I2.2 A I6 A F6.4)') '  time ',time(1), ':', time(2), ':', time(3) , &
                            ', STEPS', STEPS, ', timeStep ' , timeStep
            call print_log_message(trim(buff), 5)

            rho_micro = 0.0_dp
            rho_coherent_micro = 0.0_dp

            do i_m=1, TRAJECTORIES/MICROGROUP
            !$OMP PARALLEL DO PRIVATE(i,r,a,b)
            do m=1, MICROGROUP
                i = (i_m-1)*MICROGROUP+m ! denotes number of trajectories

                !write(*,*) m, i_m, MICROGROUP


                call random_number(r)
                a = INT(r * Nl + 1.0_dp)
                call random_number(r)
                b = INT(r * Nl + 1.0_dp)

                if(only_coherences) then
                    b = 1
                end if

!@!                if(run > 1) then
!@!                    ! preparation of initial condition
!@!                    rho_init(:,:) = 0.0_dp
!@!                    rho_init(i0,j0) = 1.0_dp
!@!                    !!!call get_coherent_dynamics(rho_init(:,:), (STEPS*(run-1)+1)*timeStep)
!@!
!@!                    if(rho_micro_N(run-1,i_m) == 0) then
!@!                        cycle
!@!                    end if
!@!
!@!                    rho_init(:,:) = rho_micro(:,:,run-1,i_m)/p_of_no_jump/rho_micro_N(run-1,i_m)
!@!                    rho_init(:,:) = rho_micro(:,:,run-1,i_m)/p_of_no_jump_measured/MICROGROUP
!@!
!@!                    if(abs(rho_init(a,b)) < 1e-6) then
!@!                        cycle
!@!                    end if
!@!
!@!
!@!                    if(m == 1) then
!@!                        write(*,*) 'p_of_no_jump,MICROGROUP',p_of_no_jump, p_of_no_jump_measured, MICROGROUP
!@!                        write(*,*) 'rho_init', rho_init
!@!                        write(*,*) 'rho_micro - source', rho_micro(:,:,run-1,i_m)
!@!                        write(*,*) 'rho_micro_N', rho_micro_N(run-1,i_m)
!@!                    end if
!@!
!@!
!@!                    call generate_trajectory(draha(:,:,m), factor, rho_init(a,b), a, b, i0, j0, run)
!@!
!@!                else
!@!                    rho_init(:,:) = 0.0_dp
!@!                    rho_init(i0,j0) = 1.0_dp

                    call generate_trajectory(draha(:,:,m), factor(:,m), 1.0_dp, i0, j0, i0, j0, run)

                    if(i_m == TRAJECTORIES/MICROGROUP .and. m == MICROGROUP) then
                        p_of_no_jump_measured = abs(rho(i0,j0,1))/i
                        write(*,*) 'p_of_no_jump_measured', p_of_no_jump_measured
                    end if
!@!                end if

                if(mod(i,1000) == 0) then
                            call itime(time)
                            write(buff, '(A I2 A I2.2 A I2.2 A I6 A I6 A I6)') '  time ',time(1), ':', time(2), ':', time(3) ,  &
                             ', STEPS/1000: ',i/1000+(run-1)*TRAJECTORIES/1000, ' ', i_m, ' ', m
                            call print_log_message(trim(buff), 5)
                end if

                if(g_functions) then
                if(exciton_basis_unraveling) then
                    call calculate_Gfactor_from_trajectory_history_general_basis(draha(:,:,m),cmplx(1,0,dp),Gfactor(:,m))
                else
                    call calculate_Gfactor_from_trajectory_history(draha(:,:,m),cmplx(1,0,dp),Gfactor(:,m))
                end if
                end if

                call calculate_Ifactor_from_trajectory_history(draha(:,:,m),cmplx(1,0,dp),Ifactor(:,m))
                Ifactor(:,m) = Ifactor(:,m)/Ifactor(STEPS*(run-1)+1,m)    ! coherent at the restart point
                Ifactor(1:STEPS*(run-1),m) = 0.0_dp                ! zero before it

                do j=1, STEPS
                    if(maxval(draha(:,j+STEPS*(run-1),m)) > 0 .and. minval(draha(:,j+STEPS*(run-1),m)) <= Nl) then
                        rho_micro(draha(1,j+STEPS*(run-1),m),draha(2,j+STEPS*(run-1),m),j+STEPS*(run-1),m) =             &
                        rho_micro(draha(1,j+STEPS*(run-1),m),draha(2,j+STEPS*(run-1),m),j+STEPS*(run-1),m)               &
                        + factor(j+STEPS*(run-1),m)*conjg(Gfactor(j+STEPS*(run-1),m))*Ifactor(j+STEPS*(run-1),m)

                        rho_coherent_micro(draha(1,j+STEPS*(run-1),m),draha(2,j+STEPS*(run-1),m),j+STEPS*(run-1),m) =    &
                        rho_coherent_micro(draha(1,j+STEPS*(run-1),m),draha(2,j+STEPS*(run-1),m),j+STEPS*(run-1),m)      &
                        + factor(j+STEPS*(run-1),m)*Ifactor(j+STEPS*(run-1),m)
                    end if
                end do

!@!                if(i <= TRAJECTORIES_STORED) then
!@!                    depository_tmp(i,:,:) = draha(:,:,m)
!@!                end if

            end do
            !$OMP END PARALLEL DO
            end do

            write(*,*) 'A'

            ! collect parallel results
            do m=1, MICROGROUP
              rho(:,:,:) = rho(:,:,:) +  rho_micro(:,:,:,m)
              rho_coherent(:,:,:) = rho_coherent(:,:,:) +  rho_coherent_micro(:,:,:,m)
            end do

            if(max(1+STEPS*(run-1),1+STEPS*(run-1) + 5) > size(rho_coherent,3)) then
                call print_error_message(-1,"error 42")
            end if

            !norm = abs(maxval(abs(rho_coherent(:,1,1+STEPS*(run-1):STEPS*run))))
            norm = abs(maxval(abs(rho_coherent(:,:,1+STEPS*(run-1):1+STEPS*(run-1) + 5)))) ! do not include the end
            write(*,*) 'norm:',norm
            if(norm == 0.0) then
              call print_error_message(-1,'Montecarlo zero norm!')
            end if
            !do j=STEPS, 1, -1
                !if(.not. only_coherences) then
                !    rho(:,:,j+STEPS*(run-1)) = rho(:,:,j+STEPS*(run-1))/trace(rho(:,:,j+STEPS*(run-1)))
                !    rho_coherent(:,:,j+STEPS*(run-1)) = rho_coherent(:,:,j+STEPS*(run-1))/trace(rho_coherent(:,:,j+STEPS*(run-1)))
                !else
                !    do jj=1, Nl
                !       rho(jj,1,j+STEPS*(run-1)) = rho(jj,1,j+STEPS*(run-1))/norm
                !       rho_coherent(jj,1,j+STEPS*(run-1)) = rho_coherent(jj,1,j+STEPS*(run-1))/norm
                !    end do
                !end if
            !end do

             rho = rho/norm
             rho_coherent = rho_coherent/norm


!            trajectory_depository = depository_tmp
!            trajectory_depository = 0
!            call collect_rho(run)
!            write(*,*) 'run',run

        end do

        !write(*,*) 'fIfGf_cor', fIfGf_cor/TRAJECTORIES
        !write(*,*)
        !write(*,*) 'fIfGf_cor - GI', fIfGf_cor/TRAJECTORIES - (Gf_cor*fIf_cor)/TRAJECTORIES/TRAJECTORIES
        !write(*,*)
        !write(*,*) 'Gf_cor', Gf_cor/TRAJECTORIES
        !write(*,*)
        !write(*,*) 'f_cor', f_cor/TRAJECTORIES
        !write(*,*)
        !write(*,*) 'If_cor', If_cor/TRAJECTORIES
        !write(*,*)
        !write(*,*) 'fIf_cor', fIf_cor/TRAJECTORIES
        !write(*,*)

    end subroutine perform_montecarlo

    recursive function goft_site(m,n,t) result(dgoft)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: dgoft
        integer    (i4b)            :: a, t_index
        character(len=300)          :: buff

        if(t < 0) then
            dgoft = conjg(goft_site(n,m,-t))
            return
        end if

        t_index = INT(t/dt)+1

        dgoft = 0.0_dp

        if(m == n) then

            if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex)) then
                write(buff,'(i2)') m
                buff = "m="//trim(buff)//" exceeds Nl1 size in goft_site()"
                call print_error_message(-1,buff)
                stop
            end if

            if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
                write(*,*) t, t_index, m,  size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)
                call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
                stop
            end if

            dgoft = all_goft(iblocks(1,1)%sblock%gindex(m))%gg(t_index)

        end if

    end function goft_site

    recursive function hoft_site(m,n,t1,t2) result(dgoft)
    ! hoft(a,a,t1,t2) = hoft(a,a,-t2,-t1)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: dgoft

        dgoft = 0.0_dp

        if(m == n) then
            dgoft = goft_site(m,m,t1) - goft_site(m,m,t1-t2) + goft_site(m,m,-t2)
        end if

    end function hoft_site

    subroutine deinit_goft_general()

        DEALLOCATE(gg_MC)
        DEALLOCATE(hh_MC)
        DEALLOCATE(cc_MC)

    end subroutine deinit_goft_general

    subroutine init_goft_general()
        integer(i4b) :: i,j,k,l,r,s,t_index, NNN
        real(dp), dimension(:,:), pointer :: SSS, SS1

        NNN = size(iblocks(1,1)%eblock%S1,1)

        ALLOCATE(gg_MC, (NNN, NNN, NNN, NNN, size(all_goft(1)%gg,1)) )
        ALLOCATE(hh_MC, (NNN, NNN, NNN, NNN, size(all_hoft(1)%gg,1)) )
        ALLOCATE(cc_MC, (NNN, NNN, NNN, NNN, size(all_coft(1)%gg,1)) )

        gg_MC = 0.0_dp
        hh_MC = 0.0_dp
        cc_MC = 0.0_dp

        SSS => iblocks(1,1)%eblock%SS
        SS1 => iblocks(1,1)%eblock%S1

        if(exciton_basis_unraveling) then

            ! init of exciton goft and hoft
            do t_index=1,size(all_goft(iblocks(1,1)%sblock%gindex(i))%gg)

            do i=1,NNN
            do j=1,NNN
            do k=1,NNN
            do l=1,NNN

            do r=1,NNN
            do s=1,NNN
                gg_MC(i,j,k,l,t_index) =     SS1(i,r) * all_goft(iblocks(1,1)%sblock%gindex(r))%gg(t_index) * SSS(r,j) &
                                                * SS1(k,s) * all_goft(iblocks(1,1)%sblock%gindex(s))%gg(t_index) * SSS(s,l)

                hh_MC(i,j,k,l,t_index) =     SS1(i,r) * all_hoft(iblocks(1,1)%sblock%gindex(r))%gg(t_index) * SSS(r,j) &
                                                * SS1(k,s) * all_hoft(iblocks(1,1)%sblock%gindex(s))%gg(t_index) * SSS(s,l)

                cc_MC(i,j,k,l,t_index) =     SS1(i,r) * all_coft(iblocks(1,1)%sblock%gindex(r))%gg(t_index) * SSS(r,j) &
                                                * SS1(k,s) * all_coft(iblocks(1,1)%sblock%gindex(s))%gg(t_index) * SSS(s,l)
            end do
            end do

            end do
            end do
            end do
            end do

            end do

        else

            do i=1,NNN
            do t_index=1,size(all_hoft(iblocks(1,1)%sblock%gindex(i))%gg)
                gg_MC(i,i,i,i,t_index) = all_goft(iblocks(1,1)%sblock%gindex(i))%gg(t_index)
                hh_MC(i,i,i,i,t_index) = all_hoft(iblocks(1,1)%sblock%gindex(i))%gg(t_index)
                cc_MC(i,i,i,i,t_index) = all_coft(iblocks(1,1)%sblock%gindex(i))%gg(t_index)
            end do
            end do
        end if

    end subroutine init_goft_general

    ! this is basicaly ordinary correlation function in general basis, so
    ! it has four idices, since it is quantity of second order in SBC
    recursive function VV_general(m,n,   i,j,t1,t2) result(res)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n,i,j
        complex(dpc)                :: res

        integer(i4b)                :: t_index
        real(dp)                    :: t

        t = t1 - t2

        if(t < 0) then
            res = conjg(VV_general(i,j,m,n,t2,t1))
            return
        end if

        t_index = INT(t/dt)+1

        res = 0.0_dp

        !!! what about terms ii-jj, are they in this term too, or are they in g-fctions of excitons?
        if(exciton_basis_unraveling) then
            res = cc_MC(m,n,i,j,t_index)
        end if

    end function VV_general

    function goft_general_2(m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = goft_general_4(m,m,n,n,t)
    end function goft_general_2

    recursive function goft_general_4(i,j,m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: i,j,m,n
        complex(dpc)                :: res
        integer (i4b)               :: a, t_index
        character(len=300)          :: buff

        if(t < 0) then
            res = conjg(goft_general_4(m,n,i,j,-t))
            return
        end if

        t_index = INT(t/dt)+1

        res = 0.0_dp

!        if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex) .or. n < 1 .or. n > size(iblocks(1,1)%sblock%gindex)) then
!            write(buff,'(i2)') m
!            buff = "m="//trim(buff)//" exceeds Nl1 size in goft_site()"
!            call print_error_message(-1,buff)
!            stop
!        end if
!
!        if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
!            write(*,*) t, t_index, m,  size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)
!            call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
!            stop
!        end if

        res = gg_MC(i,j,m,n,t_index)

    end function goft_general_4

    function dgoft_general_2(m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = dgoft_general_4(m,m,n,n,t)
    end function dgoft_general_2

    recursive function dgoft_general_4(i,j,m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n,i,j
        complex(dpc)                :: res
        integer (i4b)               :: a, t_index
        character(len=300)          :: buff

        if(t < 0) then
            res = -conjg(dgoft_general_4(m,n,i,j,-t))
            return
        end if

        t_index = INT(t/dt)+1

        res = 0.0_dp

!        if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex) .or. n < 1 .or. n > size(iblocks(1,1)%sblock%gindex)) then
!            write(buff,'(i2)') m
!            buff = "m="//trim(buff)//" exceeds Nl1 size in goft_site()"
!            call print_error_message(-1,buff)
!            stop
!        end if
!
!        if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
!            write(*,*) t, t_index, m,  size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)
!            call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
!            stop
!        end if

        res = hh_MC(i,j,m,n,t_index)

    end function dgoft_general_4

    function hoft_general(m,n,t1,t2) result(res)
    ! hoft(a,a,t1,t2) = hoft(a,a,-t2,-t1)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = 0.0_dp

        res = goft_general(m,n,t1) - goft_general(m,n,t1-t2) + goft_general(m,n,-t2)

    end function hoft_general

    recursive function hV_general(i,   m,n,th,tV) result(res)
        real(dp), intent(in)        :: th,tV
        integer(i4b), intent(in)    :: m,n,i
        complex(dpc)                :: res

        res = 0.0_dp

        res = goft_general(i,i,m,n,th) - dgoft_general(i,i,m,n,-tV) + dgoft_general(i,i,m,n,th-tV)

    end function hV_general

    recursive function Vh_general(m,n,   i,tV,th) result(res)
        real(dp), intent(in)        :: th,tV
        integer(i4b), intent(in)    :: m,n,i
        complex(dpc)                :: res

        res = 0.0_dp

        res = dgoft_general(m,n,i,i,tV) + goft_general(m,n,i,i,-th) - dgoft_general(m,n,i,i,tV-th)

    end function Vh_general

    !
    ! Defines the jump probability at given circumstances
    !
    function jump_probability_444(state_where_I_am,state_where_to_go,side) result(res)
        integer(i4b), intent(in)   :: state_where_I_am
        integer(i4b), intent(in)   :: side, state_where_to_go
        real(dp)                   :: res

        res = 0.0_dp

        if(exciton_basis_unraveling) then
            if(g_functions) then
                res = 1.0_dp
            end if
            ! it would be probably optimal to jump proportionaly to reorganization energy,
            ! but it will require some thinking since
        else
            ! technically, it should (probably) be state_where_I_am,state_where_to_go for side = 1
            ! and the other way for side = 2, but since the J is symmetric, it doesn't matter
            res = J_coupl(state_where_I_am,state_where_to_go)*CoherentFactorMultiplierToLowerJumpProb
        end if

    end function jump_probability_444

    ! just for interface-purposes
    function jump_probability_144(state_where_I_am,state_where_to_go,side) result(res)
        integer(i1b), intent(in)                 :: state_where_I_am
        integer(i4b), intent(in)                 :: side, state_where_to_go
        real(dp)                                 :: res

        integer(i4b)                             :: a

        a = state_where_I_am

        res = jump_probability_444(a,state_where_to_go,side)
    end function jump_probability_144

    !
    ! Total jump probability at given circumstances
    !
    function jump_probability_total_IJ(i0, j0) result(cumulative_probability)
        integer(i4b), intent(in)   :: i0, j0
        real(dp)                   :: cumulative_probability

        integer(i4b)               :: a

        cumulative_probability = 0.0_dp
        do a=1,Nl
            cumulative_probability = cumulative_probability + jump_probability(i0, a, 1)
        end do
        if(.not. only_coherences) then
            do a=1,Nl
                cumulative_probability = cumulative_probability + jump_probability(j0, a, 2)
            end do
        end if
    end function

    function jump_probability_total_() result(cumulative_probability)
        real(dp)                   :: cumulative_probability

        integer(i4b)               :: a

        cumulative_probability = 0.0_dp
        do a=1,Nl
            cumulative_probability = cumulative_probability + jump_probability_total_IJ(a,a)
        end do

    end function

    !
    ! Defines the jump phase change at given circumstances
    !
    function jump_phase_change(trajectory,side,i,a) result(res)
        integer(i1b), intent(in), dimension(:,:) :: trajectory
        integer(i4b), intent(in)                 :: side, i,a
        real(dp)                                 :: res

        res = 0.0_dp

        if(exciton_basis_unraveling) then
            ! phase change is the same as in site basis, because of the i/hbar in Schrodinger Eq.
            if(side == 1) then
                res = cmplx(0,-1,dpc)/CoherentFactorMultiplierToLowerJumpProb
            else if(side == 2) then
                res = cmplx(0,+1,dpc)/CoherentFactorMultiplierToLowerJumpProb
            end if
        else
            if(side == 1) then
                res = cmplx(0,-1,dpc)/CoherentFactorMultiplierToLowerJumpProb
            else if(side == 2) then
                res = cmplx(0,+1,dpc)/CoherentFactorMultiplierToLowerJumpProb
            end if
        end if
    end function jump_phase_change    !

    !
    ! Generates trajectory with STEPS steps and related pure-coupling-factors
    !
    subroutine generate_trajectory_re(trajectory, factor_out, factor_in, i0, j0, i00, j00, run)
        integer(i1b), intent(out), dimension(:,:) :: trajectory
        integer(i4b), intent(in)    :: i0, j0, run, i00, j00
        real(dp), intent(in)        :: factor_in
        complex(dpc), intent(out), dimension(:)    :: factor_out

        complex(dpc) :: ffactor

        ffactor = factor_in

        call generate_trajectory_cmplx(trajectory, factor_out, ffactor, i0, j0, i00, j00, run)
    end subroutine generate_trajectory_re

    recursive subroutine generate_trajectory_cmplx(trajectory, factor_out, factor_in, i0, j0, i00, j00, run)
        integer(i1b), intent(out), dimension(:,:) :: trajectory
        integer(i4b), intent(in)    :: i0, j0, run, i00, j00 ! beginning of run-th run and beginning of 1-st run
        complex(dpc), intent(in)    :: factor_in
        complex(dpc), intent(out), dimension(:)    :: factor_out

        integer(i4b) :: i, a, b, side, inside_run, u, v
        real(dp) :: cumulative_random, cumulative_probability

        if(size(trajectory,2) /= STEPS*RUNS .or. size(trajectory,1) /= 2 .or. &
            size(factor_out,1) /= STEPS*RUNS .or. run < 1 .or. run > RUNS) then
            call print_error_message(-1, "dimension error in generate_trajectory()")
        end if

        if((.not.(i0 == -1 .and. j0 == -1)) .or. depository) then
            ! if depository is on or we are in first recursion, we zero the trajectory and factor
            factor_out = 0.0_dp
            trajectory = 0

            trajectory(1, STEPS*(run-1) + 1) = i0
            trajectory(2, STEPS*(run-1) + 1) = j0
        end if

!@!        if(run > 1) then
!@!            ! we add random histories from depository to trajectory, backwards
!@!            if(depository) then
!@!                ! in this case, i00 and j00 are ignored - they are contained in depository
!@!                do inside_run=run-1,1,-1
!@!                    call random_number(cumulative_random)
!@!                    a = int(cumulative_random * TRAJECTORIES_STORED + 1.0_dp)
!@!                    b = 0 ! too high b signalizes there is not enough trajectories in that element
!@!
!@!                        do while(.not.(                                       &
!@!                     trajectory_depository(a, 1, STEPS*inside_run) ==         &
!@!                         trajectory(1, STEPS*inside_run + 1)                  &
!@!                         .and.                                                &
!@!                     trajectory_depository(a, 2, STEPS*inside_run) ==         &
!@!                         trajectory(2, STEPS*inside_run + 1)) )
!@!
!@!                        if(a < 1 .or. a > TRAJECTORIES_STORED) then
!@!                            call print_error_message(-1, "wrong 'a' generated in generate_trajectory()")
!@!                        end if
!@!                        if(b > TRAJECTORIES_STORED) then
!@!                            call print_log_message("no valid history to add in generate_trajectory()", 5)
!@!                            exit
!@!                        end if
!@!
!@!                        call random_number(cumulative_random)
!@!                        a = int(cumulative_random * TRAJECTORIES_STORED + 1.0_dp)
!@!                        b = b + 1
!@!                    end do
!@!
!@!                    if(b > TRAJECTORIES_STORED) then
!@!                        trajectory = 0
!@!                    else
!@!                        trajectory(:,STEPS*(inside_run-1)+1:STEPS*inside_run) = &
!@!                            trajectory_depository(a,:,STEPS*(inside_run-1)+1:STEPS*inside_run)
!@!                    end if
!@!                end do
!@!            else
!@!                if(run - 1 == 1) then
!@!                    call generate_trajectory_cmplx(trajectory, factor_out, factor_in, i00, j00, i00, j00, run-1)
!@!                end if
!@!
!@!                if(.not.(i0 == -1 .and. j0 == -1) .and. run - 1 > 1) then
!@!                    do while(.not.(                                                &
!@!                     i0 ==    trajectory(1, STEPS*run)                             &
!@!                         .and.                                                     &
!@!                     j0 ==    trajectory(2, STEPS*run)    ) )
!@!
!@!                         call generate_trajectory_cmplx(trajectory, factor_out, factor_in, -1, -1, i00, j00, run-1)
!@!                end do
!@!                end if
!@!
!@!            end if
!@!        end if

        factor_out(STEPS*(run-1)+1) = factor_in
        trajectory(1,STEPS*(run-1)+1) = i0
        trajectory(2,STEPS*(run-1)+1) = j0

        do i=1+STEPS*(run-1), STEPS-1+STEPS*(run-1)
            factor_out(i+1) = factor_out(i)
            trajectory(1,i+1) = trajectory(1,i)
            trajectory(2,i+1) = trajectory(2,i)

            call random_number(cumulative_random)
            cumulative_probability = 0.0_dp

            do a=1,Nl
            do side=1,2
                if(only_coherences .and. side == 2) then
                    exit
                end if

                cumulative_probability = cumulative_probability + jump_probability(trajectory(side,i),a,side)*timeStep

                if(cumulative_random < cumulative_probability) then
                    trajectory(side,i+1) = a
                    factor_out(i+1) = factor_out(i)*jump_phase_change(trajectory,side,i,a)
                    if(.not. modified_unraveling2) then
                        factor_out(1:i) = 0.0_dp
                    end if
                    goto 42
                end if

            end do
            end do

42          continue
        end do

        ! normalization according to the jump probability
        if(modified_unraveling2) then
            cumulative_probability = jump_probability_total(i0,j0)*timeStep

            do i=1+STEPS*(run-1), STEPS-1+STEPS*(run-1)
                factor_out(i) = factor_out(i) / (1.0 - cumulative_probability)**(i-STEPS*(run-1))
            end do
        end if

    end subroutine generate_trajectory_cmplx

    subroutine calculate_Ifactor_from_trajectory_history(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in)     :: trajectory
        complex(dpc), intent(in)                     :: initial_factor
        complex(dpc), dimension(:), intent(out)      :: factor_out

        real(dp)      :: time
        integer(i4b)  :: i, j, k

        if(size(trajectory,1) /= 2 .or. size(factor_out) /= size(trajectory,2)) then
            call print_error_message(-1, "dimension error in calculate_Ifactor_from_trajectory_history()")
        end if

        if(maxval(trajectory) > Nl) then
            call print_error_message(-1, "value error in calculate_Ifactor_from_trajectory_history()")
        end if

        factor_out = 0.0_dp
        factor_out(1) = initial_factor
        if(exciton_basis_unraveling) then

        ! exciton basis
        do i=2,size(factor_out)
            if(minval(trajectory(:,i)) < 1) then
                exit
            else
                if(.not. only_coherences) then
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (iblocks(1,1)%eblock%en(trajectory(1,i)) - iblocks(1,1)%eblock%en(trajectory(2,i))) )
                else
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (iblocks(1,1)%eblock%en(trajectory(1,i)) - rwa))
                end if
            end if
        end do

        else

        ! local basis
        do i=2,size(factor_out)
            if(minval(trajectory(:,i)) < 1) then
                exit
            else
                if(.not. only_coherences) then
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (iblocks(1,1)%sblock%en(trajectory(1,i)) - iblocks(1,1)%sblock%en(trajectory(2,i))) )
                else
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (iblocks(1,1)%sblock%en(trajectory(1,i)) - rwa) - debug_gamma*timeStep)
                end if
            end if
        end do

        end if ! exciton_basis unraveling and others

    end subroutine calculate_Ifactor_from_trajectory_history

     ! Let we have U^+_A(t1)U^+_B(t2)U^+_C(t3)...U_c(t3)U_b(t2)U_a(t1)W, then the cummulant expansion is derived
     ! from it starts as (1-i h_A(t1)-g_A(t1))^*(1-i h_B(t1)-g_B(t1))(1-i h_B(t1+t2)-g_B(t1+t2))^*
     ! (1-i h_C(t1+t2)-g_C(t1+t2)) ... (1-i h_Z(t1+..+tn)-g_Z(t1+..+tn))^* (1-i h_z(t1+..+tn)-g_z(t1+..+tn))
     ! (1-i h_z(t1+..+t(n-1))-g_z(t1+..+t(n-1)))^* (1-i h_y(t1+..+t(n-1))-g_y(t1+..+t(n-1)))
    subroutine calculate_Gfactor_from_trajectory_history(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in)  :: trajectory
        complex(dpc), intent(in)                  :: initial_factor
        complex(dpc), dimension(:), intent(out)   :: factor_out

        ! we overallocate jump-fields not to have to always calculate their dimensions
        integer(i4b), parameter :: overallocation = 10000
        integer(i4b), dimension(0:overallocation)    ::    &
                            jumps_index, jumps_tovalue_ket, jumps_tovalue_bra
        real(dp), dimension(0:overallocation)        ::    &
                            jumps_time
        logical, dimension(0:overallocation)         ::    &
                            jumps_ket

        integer(i4b) ::     i, j, k, number_of_jumps, last_index, tovalue,        &
                            fromvalue, ket_i, bra_i, ket_j, bra_j

        real(dp)     :: time_i, time_j, time
        logical      :: last_round_i
        complex(dpc) :: additive_factor

        complex(dpc) :: additive_factor_not_last_turn_z_minuleho
        integer(i4b) :: max_i_minuleho
        logical      :: lze_pouzit_AFNLTZM

        logical, parameter :: debug_G = .false.

        if(size(trajectory,1) /= 2 .or. size(factor_out) /= size(trajectory,2)) then
            call print_error_message(-1, "dimension error in calculate_factor_from_trajectory_history()")
        end if

        ! we prepare fields representing the jumps
        jumps_index(0)          = 1
        jumps_time(0)           = 0.0_dp
        jumps_tovalue_ket(0)    = trajectory(1,1)
        jumps_tovalue_bra(0)    = trajectory(2,1)

        tovalue = 1
        fromvalue = 1
        i = 1
        do while(.true.)
            if(tovalue == -1 .or. fromvalue == -1) then
                exit
            end if
            if(i > overallocation) then
                call print_error_message(-1,'overallocation insufficient in calculate_Gfactor_from_trajectory_history()')
            end if

            call next_step_of_trajectory(trajectory, jumps_index(i-1), jumps_index(i), fromvalue, tovalue, jumps_ket(i))

            if(jumps_ket(i)) then
                jumps_tovalue_ket(i) = tovalue
                jumps_tovalue_bra(i) = jumps_tovalue_bra(i-1)
            else
                jumps_tovalue_ket(i) = jumps_tovalue_ket(i-1)
                jumps_tovalue_bra(i) = tovalue
            end if

            jumps_time(i) = (jumps_index(i)-1)*timeStep
            i = i + 1
        end do

        number_of_jumps = i-2
        max_i_minuleho = -1


        ! now we perform evaluation of brackets of g and h-functions

        ! g-functions come first
        factor_out = 0.0_dp

        do k=1, size(factor_out)

            time = (k-1)*timeStep

            if(trajectory(1,k) < 1 .or. trajectory(1,k) > Nl .or. trajectory(2,k) < 1 .or. trajectory(2,k) > Nl) then
                return
            end if

            factor_out(k) = initial_factor
            additive_factor = 0.0_dp

            if(debug_G) then
                write(*,*)
                write(*,*) 'jumps_tovalue_ket', jumps_tovalue_ket(:5)
                write(*,*) 'jumps_tovalue_bra', jumps_tovalue_bra(:5)
                write(*,*) 'jumps_time', jumps_time(:5)
            end if

            ! last operator has time t instead of time of last jump
            do i=1, number_of_jumps + 1
                if(jumps_time(i) >= time .or. i > number_of_jumps) then
                    last_round_i = .true.
                    time_i = time
                    ket_i = trajectory(1,k)
                    bra_i = trajectory(2,k)

                    if(max_i_minuleho == i) then
                        lze_pouzit_AFNLTZM = vynechani_G_ifu
                    else
                        lze_pouzit_AFNLTZM = .false.
                        additive_factor_not_last_turn_z_minuleho = -9999
                        ! if not correctly set later, it will screw up everything and make the mistake obvious
                    end if
                    max_i_minuleho = i

                else
                    last_round_i = .false.
                    time_i = jumps_time(i)
                    ket_i = jumps_tovalue_ket(i)
                    bra_i = jumps_tovalue_bra(i)

                    if(lze_pouzit_AFNLTZM) then
                        additive_factor = additive_factor_not_last_turn_z_minuleho
                        cycle ! we wait until we are in the last round
                    end if
                end if

                ! there are always two 1hg-brackets for both bra and ket for each jump-time, with site i and (i-1).
                ! exception is the last round, where the bra and ket meet into two brackets instead of four

                if(.not. last_round_i) then
                    ! two ket-brackets
                    additive_factor = additive_factor+(-goft_site(ket_i,ket_i,time_i))
                        if(debug_G) then
                            write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if
                    additive_factor = additive_factor+&
                            conjg((-goft_site(jumps_tovalue_ket(i-1),jumps_tovalue_ket(i-1),time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1)
                        end if

                    ! two bra-brackets
                    if(.not. only_coherences) then
                        additive_factor = additive_factor+conjg((-goft_site(bra_i,bra_i,time_i)))
                            if(debug_G) then
                                write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                        additive_factor = additive_factor+&
                                (-goft_site(jumps_tovalue_bra(i-1),jumps_tovalue_bra(i-1),time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i,jumps_tovalue_bra(i-1)
                            end if
                    end if


                    ! hh between bra and ket brackets of the same i (4 terms => 6 pairs)
                    if(ket_i == bra_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H1', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H2', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, jumps_tovalue_bra(i-1)
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == bra_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H3', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1), bra_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(i-1),jumps_tovalue_bra(i-1),&
                                                                    time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H4', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1), jumps_tovalue_bra(i-1)
                            end if
                    endif
                    if(ket_i == jumps_tovalue_ket(i-1)) then
                        additive_factor = additive_factor+(hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H6', real(additive_factor), k, i, time, time_i, last_round_i, bra_i, bra_i
                            end if
                    endif

                ! one additional factor where bra and kets meet
                else
                    additive_factor = additive_factor+conjg((-goft_site(ket_i,ket_i,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_g--', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if

                    if(.not. only_coherences) then
                        additive_factor = additive_factor+(-goft_site(bra_i,bra_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg--', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                    end if

                    if(bra_i == ket_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_h (last)', real(additive_factor), k, i, time, time_i, last_round_i
                            end if
                    end if
                end if


                do j=1, i-1
!                    if(last_round_i .and. j == i-1) then
!                        exit
!                    end if

                    time_j = jumps_time(j)
                    ket_j = jumps_tovalue_ket(j)
                    bra_j = jumps_tovalue_bra(j)

                    ! hh between all combinations of i-j pairs
                    if(ket_i == ket_j .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor1', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == ket_j) then
                        additive_factor = additive_factor+(hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor2', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == ket_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor3', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == ket_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor4', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_ket(j-1) .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,+time_i))
                            if(debug_G) then
                                write(*,*) 'factor5', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_ket(j-1)) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor6', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_ket(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor7', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == jumps_tovalue_ket(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,+time_i))
                            if(debug_G) then
                                write(*,*) 'factor8', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_bra(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor9', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_bra(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor10', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_bra(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_j,time_i))            !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor11', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == jumps_tovalue_bra(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_j,+time_i))        !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor12', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == bra_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(bra_j,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor13', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == bra_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(bra_j,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor14', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == bra_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(bra_j,bra_j,time_j,time_i))                          !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor15', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == bra_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_j,bra_j,time_j,time_i))                          !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor16', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                end do

                ! At the time of jump, code is not correct. It is safer to exclude these cases than fix it.
                if(abs(time-jumps_time(i)) < timeStep/2.0_dp) then
                    additive_factor = -9999999.0

                    if(debug_G) then
                        write(*,*) 'case excluded', real(additive_factor)
                    end if
                end if

                if(real(additive_factor) > 0 .and. debug_G) then
                    write(*,*) 'zbyly-faktor', real(additive_factor)
                end if

                if(last_round_i) then
                    exit
                else
                    if(.not. lze_pouzit_AFNLTZM) then
                        additive_factor_not_last_turn_z_minuleho = additive_factor
                    end if
                end if
            end do

            factor_out(k) = factor_out(k)*exp(additive_factor)

        end do

        !write(*,*) 'factor', factor_out

    end subroutine calculate_Gfactor_from_trajectory_history

     ! Let we have U^+_A(t1) V_AB U^+_B(t2) V_BC U^+_C(t3)...U_c(t3)  V_cb  U_b(t2) V_ba U_a(t1)W,
     ! then the cummulant expansion is derived from it starts as
     ! (1-i h_A(t1)-g_A(t1))^* V_AB(t1) (1-i h_B(t1)-g_B(t1))(1-i h_B(t1+t2)-g_B(t1+t2))^*  V_BC(t1+t2)
     ! (1-i h_C(t1+t2)-g_C(t1+t2)) ... (1-i h_Z(t1+..+tn)-g_Z(t1+..+tn))^* (1-i h_z(t1+..+tn)-g_z(t1+..+tn))
     ! (1-i h_z(t1+..+t(n-1))-g_z(t1+..+t(n-1)))^* V_zy(t1+..+t(n-1) ) (1-i h_y(t1+..+t(n-1))-g_y(t1+..+t(n-1)))
    subroutine calculate_Gfactor_from_trajectory_history_general_basis(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in) :: trajectory
        complex(dpc), intent(in)                    :: initial_factor
        complex(dpc), dimension(:), intent(out) :: factor_out

        ! we overallocate jump-fields not to have to always calculate their dimensions
        integer(i4b), parameter :: overallocation = 10000
        integer(i4b), dimension(0:overallocation)   ::  &
                            jumps_index, jumps_tovalue_ket, jumps_tovalue_bra
        real(dp), dimension(0:overallocation)       ::  &
                            jumps_time
        logical, dimension(0:overallocation)        ::  &
                            jumps_ket

        integer(i4b) ::     i, j, k, number_of_jumps, last_index, tovalue,      &
                            fromvalue, ket_i, bra_i, ket_j, bra_j,  ket_i_, bra_i_, ket_j_, bra_j_

        real(dp)     :: time_i, time_j, time
        logical      :: last_round_i
        complex(dpc) :: additive_factor

        complex(dpc) :: additive_factor_not_last_turn_z_minuleho
        integer(i4b) :: max_i_minuleho
        logical      :: lze_pouzit_AFNLTZM

        logical, parameter :: debug_G = .false.

        if(size(trajectory,1) /= 2 .or. size(factor_out) /= size(trajectory,2)) then
            call print_error_message(-1, "dimension error in calculate_factor_from_trajectory_history()")
        end if

        ! we prepare fields representing the jumps
        jumps_index(0)              = 1
        jumps_time(0)            = 0.0_dp
        jumps_tovalue_ket(0)    = trajectory(1,1)
        jumps_tovalue_bra(0)    = trajectory(2,1)

        tovalue = 1
        fromvalue = 1
        i = 1
        do while(.true.)
            if(tovalue == -1 .or. fromvalue == -1) then
                exit
            end if
            if(i > overallocation) then
                call print_error_message(-1,'overallocation insufficient in calculate_Gfactor_from_trajectory_history()')
            end if

            call next_step_of_trajectory(trajectory, jumps_index(i-1), jumps_index(i), fromvalue, tovalue, jumps_ket(i))

            if(jumps_ket(i)) then
                jumps_tovalue_ket(i) = tovalue
                jumps_tovalue_bra(i) = jumps_tovalue_bra(i-1)
            else
                jumps_tovalue_ket(i) = jumps_tovalue_ket(i-1)
                jumps_tovalue_bra(i) = tovalue
            end if

            jumps_time(i) = (jumps_index(i)-1)*timeStep
            i = i + 1
        end do

        number_of_jumps = i-2
        max_i_minuleho = -1


        ! now we perform evaluation of brackets of g and h-functions

        ! g-functions come first
        factor_out = 0.0_dp

        do k=1, size(factor_out)

            time = (k-1)*timeStep

            if(trajectory(1,k) < 1 .or. trajectory(1,k) > Nl .or. trajectory(2,k) < 1 .or. trajectory(2,k) > Nl) then
                return
            end if

            factor_out(k) = initial_factor
            additive_factor = 0.0_dp

            if(debug_G) then
                write(*,*)
                write(*,*) 'jumps_tovalue_ket', jumps_tovalue_ket(:5)
                write(*,*) 'jumps_tovalue_bra', jumps_tovalue_bra(:5)
                write(*,*) 'jumps_time', jumps_time(:5)
            end if

            ! last operator has time t instead of time of last jump
            do i=1, number_of_jumps + 1
                if(jumps_time(i) >= time .or. i > number_of_jumps) then
                    last_round_i = .true.
                    time_i = time
                    ket_i = trajectory(1,k)
                    bra_i = trajectory(2,k)

                    if(max_i_minuleho == i) then
                        lze_pouzit_AFNLTZM = vynechani_G_ifu
                    else
                        lze_pouzit_AFNLTZM = .false.
                        additive_factor_not_last_turn_z_minuleho = -9999
                        ! if not correctly set later, it will screw up everything and make the mistake obvious
                    end if
                    max_i_minuleho = i

                else
                    last_round_i = .false.
                    time_i = jumps_time(i)
                    ket_i = jumps_tovalue_ket(i)
                    bra_i = jumps_tovalue_bra(i)
                    ket_i_ = jumps_tovalue_ket(i-1)
                    bra_i_ = jumps_tovalue_bra(i-1)

                    if(lze_pouzit_AFNLTZM) then
                        additive_factor = additive_factor_not_last_turn_z_minuleho
                        cycle ! we wait until we are in the last round
                    end if
                end if

                ! there are always two 1hg-brackets for both bra and ket for each jump-time, with site i and (i-1).
                ! exception is the last round, where the bra and ket meet into two brackets instead of four
                !
                ! The general ordering is ket_, ket, bra, bra_
                if(.not. last_round_i) then
                    ! two ket-brackets
                    additive_factor = additive_factor+(-goft_general(ket_i,ket_i,time_i))
                        if(debug_G) then
                            write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if
                    additive_factor = additive_factor+&
                                conjg((-goft_general(ket_i_,ket_i_,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_
                        end if

                    ! two bra-brackets
                    if(.not. only_coherences) then
                        additive_factor = additive_factor+conjg((-goft_general(bra_i,bra_i,time_i)))
                            if(debug_G) then
                                write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                        additive_factor = additive_factor+&
                                (-goft_general(jumps_tovalue_bra(i-1),jumps_tovalue_bra(i-1),time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i,jumps_tovalue_bra(i-1)
                            end if
                    end if

                    ! hh between bra and ket brackets of the same i (4 terms: C(4,2) => 6 pairs)
                        ! ket_i == ket_i_
                        additive_factor = additive_factor + (hoft_general(ket_i_,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            end if

                    if(.not. only_coherences) then
                        ! ket_i == bra_i
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H1', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i
                            end if

                        ! ket_i == bra_i_
                        additive_factor = additive_factor + (-hoft_general(ket_i,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H2', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i_
                            end if

                        ! ket_i_ == bra_i
                        additive_factor = additive_factor + (-hoft_general(ket_i_,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H3', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_, bra_i
                            end if

                        ! ket_i_ == bra_i_
                        additive_factor = additive_factor + (hoft_general(ket_i_,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H4', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_, bra_i_
                            end if

                        ! bra_i == bra_i_
                        additive_factor = additive_factor + (hoft_general(bra_i,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H6', real(additive_factor), k, i, time, time_i, last_round_i, bra_i, bra_i
                            end if
                    endif

                    ! hV, Vh between bra and ket brackets of the same i (4 h-terms, 2 V-terms => 8 pairs)
                        ! ket_i == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   ket_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_i_,   ket_i_,ket_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    if(.not. only_coherences) then
                        ! bra_i == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(bra_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_i_,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if



                    ! VV term
                        additive_factor = additive_factor + (VV_general(ket_i_,ket_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    end if


                ! one additional factor where bra and kets meet
                else ! not last round
                    additive_factor = additive_factor+conjg((-goft_general(ket_i,ket_i,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_g--', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if

                    if(.not. only_coherences) then
                        additive_factor = additive_factor+(-goft_general(bra_i,bra_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg--', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                    end if

                    ! bra_i == ket_i
                    if(.not. only_coherences) then
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_i,time_i,time_i))

                            if(debug_G) then
                                write(*,*) 'factor_h (last)', real(additive_factor), k, i, time, time_i, last_round_i
                            end if
                    end if

                    ! no V-terms here, I guess
                end if


                do j=1, i-1
                    time_j = jumps_time(j)
                    ket_j = jumps_tovalue_ket(j)
                    bra_j = jumps_tovalue_bra(j)
                    ket_j_ = jumps_tovalue_ket(j-1)
                    bra_j_ = jumps_tovalue_bra(j-1)

                    ! hh between all combinations of i-j pairs (8 terms: C(8,2) - 2*C(4,2)  => 16 pairs)
                    ! the correct ordering is:  ket_j, ket_j_, ket_i, ket_i_, bra_i_, bra_i, bra_j_, bra_j
                        ! ket_i_ == ket_j
                        additive_factor = additive_factor + (hoft_general(ket_j,ket_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor2', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                        ! ket_i_ == ket_j_
                        additive_factor = additive_factor + (-hoft_general(ket_j_,ket_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor6', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                    if(.not. only_coherences) then
                        ! bra_i_ == ket_j_
                        additive_factor = additive_factor + (hoft_general(ket_j_,bra_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor8', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == ket_j
                        additive_factor = additive_factor + (-hoft_general(ket_j,bra_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor4', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i_ == bra_j_
                        additive_factor = additive_factor + (hoft_general(ket_i_,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor10', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == bra_j_
                        additive_factor = additive_factor + (-hoft_general(bra_i_,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor12', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i_ == bra_j
                        additive_factor = additive_factor + (-hoft_general(ket_i_,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor14', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == bra_j
                        additive_factor = additive_factor + (hoft_general(bra_i_,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor16', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                    endif
                    if(.not. last_round_i) then
                        ! ket_i == ket_j
                        additive_factor = additive_factor + (-hoft_general(ket_j,ket_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor1', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == ket_j_
                        additive_factor = additive_factor + (hoft_general(ket_j_,ket_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor5', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(.not. only_coherences .and. .not. last_round_i) then
                        ! bra_i == ket_j
                        additive_factor = additive_factor + (hoft_general(ket_j,bra_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor3', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == ket_j_
                        additive_factor = additive_factor + (-hoft_general(ket_j_,bra_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor7', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == bra_j_
                        additive_factor = additive_factor + (-hoft_general(ket_i,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor9', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == bra_j_
                        additive_factor = additive_factor + (hoft_general(bra_i,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor11', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == bra_j
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor13', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == bra_j
                        additive_factor = additive_factor + (-hoft_general(bra_i,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor15', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif


                    ! hV, Vh between bra and ket brackets of the same i (8 h-terms, 4 V-terms, 8 x 4 - 2*8 => 16 pairs)
                        ! ket_i == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   ket_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_j,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_j_,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    if(.not. only_coherences) then
                        ! bra_j == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_j,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_j,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_j,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_j_,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   bra_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(bra_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(bra_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(ket_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(ket_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if



                    ! VV terms C(4,2) - 2 = 4
                        ! V(ket_j_,ket_j) == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (VV_general(ket_j_,ket_j,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(bra_j,bra_j_) == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (VV_general(ket_i_,ket_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(ket_j_,ket_j) == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (VV_general(ket_j_,ket_j,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(bra_j,bra_j_) == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (VV_general(bra_i,bra_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    end if

                end do

                ! At the time of jump, code is not correct. It is safer to exclude these cases than fix it.
                if(abs(time-jumps_time(i)) < timeStep/2.0_dp) then
                    additive_factor = -9999999.0

                    if(debug_G) then
                        write(*,*) 'case excluded', real(additive_factor)
                    end if
                end if

                if(real(additive_factor) > 0 .and. debug_G) then
                    write(*,*) 'zbyly-faktor', real(additive_factor)
                end if

                if(last_round_i) then
                    exit
                else
                    if(.not. lze_pouzit_AFNLTZM) then
                        additive_factor_not_last_turn_z_minuleho = additive_factor
                    end if
                end if
            end do

            factor_out(k) = factor_out(k)*exp(additive_factor)

        end do

    end subroutine calculate_Gfactor_from_trajectory_history_general_basis

    pure subroutine next_step_of_trajectory(trajectory, from_index, tindex, fromvalue, tovalue, ket)
        integer(i1b), dimension(:,:), intent(in) :: trajectory
        integer(i4b), intent(in)                 :: from_index
        integer(i4b), intent(out)                :: tindex, fromvalue, tovalue
        logical, intent(out)                     :: ket

        integer(i4b)    :: i

        if(from_index < 0) then
            return
        end if

        do i=from_index+1, size(trajectory,2)
            if(trajectory(1,i) /= trajectory(1,i-1)) then
                ket = .true.
                fromvalue = trajectory(1,i-1)
                tovalue = trajectory(1,i)
                tindex = i
                return
            elseif(trajectory(2,i) /= trajectory(2,i-1)) then
                ket = .false.
                fromvalue = trajectory(2,i-1)
                tovalue = trajectory(2,i)
                tindex = i
                return
            end if
        end do

        fromvalue = -1
        tovalue   = -1
        tindex = size(trajectory,2)
    end subroutine next_step_of_trajectory

    pure function number_of_trajectory_steps(trajectory) result(nav)
        integer(i4b) :: nav

        integer(i1b), dimension(:,:), intent(in)    :: trajectory
        integer(i4b)                                :: tindex, fromvalue, tovalue, from_index
        logical                                     :: ket

        tovalue = 1
        fromvalue = 1
        from_index = 1
        nav = -1

        do while(tovalue /= -1)
            nav = nav + 1
            call next_step_of_trajectory(trajectory, from_index, tindex, fromvalue, tovalue, ket)
            from_index = tindex
        end do
    end function number_of_trajectory_steps





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! Calculate linear absorption spectrum from polarization
    !
    subroutine create_spect_abs()

        complex(dpc), dimension(:),   allocatable    :: sig
        complex(dpc), dimension(:,:), allocatable    :: dat
        real(dp), dimension(:), allocatable          :: rs

        integer :: on, ch, i
        real(dp) :: oma, rr

        padfac = 4
        NFFT = (2**padfac)*NFFT_basic

        allocate(sig(NFFT))
        allocate(dat(1,NFFT))
        allocate(rs(NFFT))
        allocate(MC_spect_abs(NFFT))
        rs            = 0.0_dp
        MC_spect_abs = 0.0_dp

        dom = 2.0_dp*PI_D/(NFFT*(dt*gt(1)))

        if (.not.MC_polar_1_collected) then
            call collect_polar_1()
        end if

        sig = 0.0_dp
        sig(1:Nt(1)) = MC_polar_1(1:Nt(1))

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)


        do i = 1, NFFT
            oma = (-NFFT/2 + i)*dom + rwa

!                if (i == 0) then
!                    rs(i) = abs(oma)*real(sig(NFFT))
!                else
                    rs(i) = abs(oma)*real(sig(i))
!                end if


        end do

        rr = 1 !maxval(rs(:))

        MC_spect_abs = rs * (dt*gt(1))

    end subroutine create_spect_abs

    !
    ! Calculation of the first order polarization
    !
    subroutine create_polar_1()
        integer(i4b)     :: i,j,k,b, mi,ni
        integer          :: kk, kb
        ! this will be a local pointer to the global transition dipole moments
        real(dp), dimension(:,:), pointer         :: dd, dx,dy,dz
        ! orientation averaging factors
        real(dp), dimension(:,:), allocatable     :: as_orfact
        real(dp) :: vdni, vdmi

        if (.not.allocated(MC_polar_1)) then
            allocate(MC_polar_1(1:Nt(1)))
            MC_polar_1 = 0.0_dp
        end if

        call resources_rewind_blocks()

        ! loop over blocks
        kk = 0
        kb = 1
        do

            if (use_exciton_basis) then

                dx => iblocks(1,1)%eblock%dx
                   dy => iblocks(1,1)%eblock%dy
                dz => iblocks(1,1)%eblock%dz

            else

                dx => iblocks(1,1)%sblock%dx
                dy => iblocks(1,1)%sblock%dy
                dz => iblocks(1,1)%sblock%dz

            end if


            ! this just temporary, it will be put somewhere else
            allocate(as_orfact(N1,N1))

            do ni = 1,N1
                do mi = 1,N1
                    vdmi = sqrt(dx(mi,1)*dx(mi,1) + dy(mi,1)*dy(mi,1) + dz(mi,1)*dz(mi,1))
                    vdni = sqrt(dx(ni,1)*dx(ni,1) + dy(ni,1)*dy(ni,1) + dz(ni,1)*dz(ni,1))
                    if ((vdmi==0.0_dp).or.(vdni==0.0_dp)) then
                      as_orfact(ni,mi) = 0.0_dp
                    else
                      as_orfact(ni,mi) = (1.0_dp/3.0_dp)* &
                        (dx(ni,1)*dx(mi,1) + dy(ni,1)*dy(mi,1) + dz(ni,1)*dz(mi,1))/ &
                        (vdni*vdmi)
                    end if
                end do
            end do


            if (use_exciton_basis) then
                dd => iblocks(1,1)%eblock%dd
            else
                dd => iblocks(1,1)%sblock%dd
            end if

            ! loop over time
            do i = 1, Nt(1)
                ! loop over excited states
                do k = 1, N1
                    do j = 1, N1
                        MC_polar_1(i) = MC_polar_1(i) + (dd(k,1)*dd(j,1))*as_orfact(k,j)*evops(kb,kb)%Ueg(k,1,j,1,i) !gcohs%C(i,k+kk)
                    end do
                end do
            end do

! no deallocation
!            deallocate(dd)
            deallocate(as_orfact)

            if (.not.resources_have_next_block()) exit

            kk = kk + N1
            kb = kb + 1
            call resources_next_block()


        end do

    end subroutine create_polar_1

    !
    ! Collect polarizations from all processes
    !
    subroutine collect_polar_1()
        integer :: N
        character(len=256) :: buff

        N = N_realizations_local
        call parallel_average_complex_array(MC_polar_1,N)
        MC_polar_1_collected = .true.

        write(buff,'(a,i5,a)') "Polar_1  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_polar_1

    subroutine get_coherent_dynamics(rho,t)
        complex(dpc), dimension(:,:), intent(inout)    :: rho
        real(dp), intent(in)                           :: t

        integer(i4b) :: i,j
        complex(dpc), dimension(size(rho,1),1)         :: rho_O

        if(only_coherences) then
            rho_O(:,1) = rho(:,1)
            call operator_to_exc(rho_O,'O')
            rho(:,1) = rho_O(:,1)
        else
            call operator_to_exc(rho,'E')
        end if

        do i=1,size(rho,1)
        do j=1,size(rho,2)

            if(.not. only_coherences) then
                    rho(i,j) = rho(i,j)*exp(cmplx(0,-1,dpc)*t*&
                      (iblocks(1,1)%eblock%en(i) - iblocks(1,1)%eblock%en(j)) )
                else
                    rho(i,j) = rho(i,j)*exp(cmplx(0,-1,dpc)*t*&
                      (iblocks(1,1)%eblock%en(i) - rwa) - debug_gamma*t)
            end if

        end do
        end do

        if(only_coherences) then
            rho_O(:,1) = rho(:,1)
            call operator_from_exc(rho_O,'O')
            rho(:,1) = rho_O(:,1)
        else
            call operator_from_exc(rho,'E')
        end if

    end subroutine get_coherent_dynamics

    subroutine seam_superoperators_U(index_from, index_to, type)
        integer(i4b), intent(in)    :: index_from, index_to
        character, intent(in)       :: type

        complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: AA, UU, U
        integer(i4b) :: i

        if(    type == 'O' .and. (index_from > index_to .or. index_from <= 1 .or. index_to > size(evops(1,1)%Ueg, 5) ) .or. &
            type == 'E' .and. (index_from > index_to .or. index_from <= 1 .or. index_to > size(evops(1,1)%Uee, 5) ) .or. &
            .not.(type == 'O' .or. type == 'E') ) then

            write(*,*) size(evops(1,1)%Ueg, 5), index_from, index_to
            call print_error_message(-1, "error in seam_superoperators_U")
            stop
        end if

        !debug
        do i=1, size(evops(1,1)%Ueg,5)
            write(*,*) '  ', i,evops(1,1)%Ueg(1,1,1,1,i)
        end do

        if(type == 'O') then
            call superops_4indexed_to_2indexed(evops(1,1)%Ueg(:,:,:,:,index_from), U, type)
            call inv(U, UU)
            call superops_4indexed_to_2indexed(evops(1,1)%Ueg(:,:,:,:,index_from-1), U, type)

            AA = matmul(U,UU)

            write(*,*) 'A', AA

            do i=index_from,index_to
                call superops_4indexed_to_2indexed(evops(1,1)%Ueg(:,:,:,:,i), UU, type)
                UU = matmul(AA,UU)
                call superops_2indexed_to_4indexed(UU, evops(1,1)%Ueg(:,:,:,:,i), type)
            end do
        elseif(type == 'E') then

        endif


    end subroutine seam_superoperators_U

    subroutine write_time_evolutions(type,absolute_value,conjugate_coherences,rho0)
        ! writes time evolution of the density matrix elements according to the type,
        ! type = 'E', 'O', '2'. If type == '2', arrays are allocated as 0:N1+N1*(N-1)/2,
        ! otherwise 0:N1

          character, intent(in)                         :: type
          logical    , intent(in)                       :: absolute_value, conjugate_coherences
          complex(dpc), dimension(:,:),intent(in)       :: rho0

        real(dp), dimension(:,:), allocatable           :: dd
          real(dp)                                      :: s, x, magnitude
          complex(dpc)                                  :: s_c
          integer                                       :: a,b,i,j,k,l,Ublock

          complex(dpc), dimension(:,:,:), allocatable   :: rr
         complex(dpc), dimension(:,:,:,:,:), pointer    :: actual_U
          character(len=10)    :: cha,chb
          character(len=64)    :: name

          if(.not. (type == 'E' .or. type == 'O' .or. &
              size(rho0,1) == N1_from_type(type) .or. size(rho0,2) == N2_from_type(type))) then
             call print_error_message(-1,  'type error in write_time_evolutions')
             stop
          end if

        Ublock = 1

        allocate(rr(N1_from_type(type),N2_from_type(type),Nt(1)))
          rr = 0.0d0

        ! We set indices range according to block we evaluate. Because rho0 is
        ! whole density matrix, while evolution operators are only from particular
        ! block, offset is set between these indices.
        if (type == '2') then
            actual_U => evops(Ublock,Ublock)%Ufe
        else if (type == 'E') then
            actual_U => evops(Ublock,Ublock)%Uee
        else if (type == 'O') then
            actual_U => evops(Ublock,Ublock)%Ueg
        end if

          do a = 1, N1_from_type(type)
          do b = 1, N2_from_type(type)

          write(cha,'(i1)') a
          write(chb,'(i1)') b

        if (type == 'E') then
            name = 'dens_E_'//trim(cha)//'_'//trim(chb)//'.dat'
          else if (type == 'O') then
              name = 'dens_O_'//trim(cha)//'_'//trim(chb)//'.dat'
          else if (type == '2') then
              name = 'dens_2_'//trim(cha)//'_'//trim(chb)//'.dat'
          end if

          if(absolute_value .eqv. .true.) then
              name = 'abs_' // name
          endif

          if(conjugate_coherences .eqv. .true.) then
              name = 'conjg_' // name
          endif

        open(unit=11,file=trim(file_join(out_dir,trim(name))))

          do i = 1, Nt(1)

              do k = 1, N1_from_type(type)
            do l = 1, N2_from_type(type)

                if(conjugate_coherences .eqv. .true.) then
                    !!! averaging coherences
                     rr(a,b,i) = rr(a,b,i) + (actual_U(a,b,k,l,i)*rho0(k,l)            &
                                           + conjg(actual_U(b,a,k,l,i)*rho0(k,l) )     &
                                            )/2

                else
                     rr(a,b,i) = rr(a,b,i) + actual_U(a,b,k,l,i)*rho0(k,l)
                 end if

               end do
             end do

            if (absolute_value .eqv. .true.) then
                write(11,*) (i-1)*(dt*gt(1)), abs(rr(a,b,i))
            else
                write(11,*) (i-1)*(dt*gt(1)), real(rr(a,b,i)), aimag(rr(a,b,i))
              endif


          end do

          close(unit=11)

          end do
          end do

          if (type == 'E') then
              name = 'dens_trace.dat'
              open(11,file=trim(file_join(out_dir,trim(name))))

              do i = 1, Nt(1)
                 s_c = 0.0d0
                 do a = 1, N1
                    s_c = s_c + rr(a,a,i)
                 end do

                 write(11,*) (i-1)*(dt*gt(1)), real(s_c), aimag(s_c)
              end do
        end if

          close(11)


          deallocate(rr)

    end subroutine write_time_evolutions


    !*************************************************************
    !  Writing out correlation functions
    !*************************************************************

    subroutine write_gofts()
        integer (i4b)        :: i,j
        character(len=4)     :: number_goft
        character(len=100)   :: name

        do i=1,size(all_goft)

        if(i < 10) then
            write(number_goft,'(i1)')   i
        else if (i < 100) then
            write(number_goft,'(i2)')   i
        else
            write(number_goft,'(i3)')   i
        endif

        ! goft
        name = trim('site_goft') // trim(number_goft) // '.dat'

        open(UNIT=22, FILE = trim(file_join(out_dir,trim(name))))

        j = 1
        do while (j <= size(all_goft(i)%gg))
            write(22,*) dt*(j-1),' ',real(all_goft(i)%gg(j)),' ',aimag(all_goft(i)%gg(j))
            j = j + 1
        end do

        close(UNIT=22)

        ! hoft
        name = trim('site_hoft') // trim(number_goft) // '.dat'

        open(UNIT=22, FILE = trim(file_join(out_dir,trim(name))))

        j = 1
        do while (j <= size(all_hoft(i)%gg))
            write(22,*) dt*(j-1),' ',real(all_hoft(i)%gg(j)),' ',aimag(all_hoft(i)%gg(j))
            j = j + 1
        end do

        close(UNIT=22)

        ! coft
        name = trim('site_coft') // trim(number_goft) // '.dat'

        open(UNIT=22, FILE = trim(file_join(out_dir,trim(name))))

        j = 1
        do while (j <= size(all_goft(i)%gg))
            write(22,*) dt*(j-1),' ',real(all_coft(i)%gg(j)),' ',aimag(all_coft(i)%gg(j))
            j = j + 1
        end do

        close(UNIT=22)

        end do
    end subroutine write_gofts


end module module_montecarlo

