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
	use resources

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

 	character(len=64), parameter, private :: external_dir = "external", config_filename = "config.prm", &
 	                    out_filename = 'hierarchy_out.dat', config_filename_arend = "configA.prm"




    integer(i4b), private:: Nsys ! system size
    integer(i4b), private:: tmax = 6 ! depth of the hierarchy
    integer(i4b), private:: Ntimestept1   = 0 ! number of time steps during t1 in outer loop
    integer(i4b), private:: Ntimestept1in = 0 ! number of time steps during t1 in inner loop
    integer(i4b), private:: Ntimestept2   = 0 ! number of time steps during t2 in outer loop
    integer(i4b), private:: Ntimestept2in = 0 ! number of time steps during t2 in inner loop
    integer(i4b), private:: Ntimestept3   = 0 ! number of time steps during t3 in outer loop
    integer(i4b), private:: Ntimestept3in = 0 ! number of time steps during t3 in inner loop

    real(dp), private :: central_frequency = 0
    real(dp), private :: radiation_temperature = 5500
    real(dp), private :: max_light_time    = 1e9
    real(dp), private :: global_relax_rate = 0
    real(dp), private :: bloch_strength    = 0
    logical, private  :: gaussian_pulse    = .false.
    logical, private  :: normalize_trace   = .false.
    logical, private  :: complex_CF        = .false.
    logical, private  :: bloch_term        = .false.
    logical, private  :: light_hierarchy   = .false.
    logical, private  :: noise_term        = .false.
    integer(i4b)      :: realizations      = 0

    logical, parameter, private:: calculatereph = .true.
    logical, parameter, private:: calculatenonreph = .true.
    logical, parameter, private:: exciton_basis = .true.

    integer:: Nind ! number of indices in hierarchy
    complex(dpc), allocatable, private:: rho1(:,:), prhox1(:,:), prhodx1(:,:)
    complex(dpc), allocatable, private:: rho2(:,:,:), prhox2(:,:,:), prhodx2(:,:,:)
    complex(dpc), allocatable, private:: rho3s(:,:), prhox3s(:,:), prhodx3s(:,:)
    complex(dpc), allocatable, private:: rho3(:,:,:), prhox3(:,:,:), prhodx3(:,:,:)
    complex(dpc), allocatable, private:: rhoC(:,:,:), prhoxC(:,:,:), prhodxC(:,:,:)
    complex(dpc), allocatable, private:: HS(:,:), HS2(:,:)
    complex(dpc), allocatable, private:: V(:, :,:), V2(:,:,:), mu(:,:)
    complex(dpc), allocatable, private:: opLeft2(:,:,:), opRight2(:,:,:), opLRLeft2(:,:,:,:), opLRRight2(:,:,:,:)
    complex(dpc), allocatable, private:: opPlusLeft2(:,:,:,:), opPlusRight2(:,:,:,:), opMinLeft2(:,:,:,:), opMinRight2(:,:,:,:)
    !complex(dpc), allocatable, private:: opLeftC(:,:,:), opRightC(:,:,:), opLRLeftC(:,:,:,:), opLRRightC(:,:,:,:)
    !complex(dpc), allocatable, private:: opPlusLeftC(:,:,:,:), opPlusRightC(:,:,:,:), opMinLeftC(:,:,:,:), opMinRightC(:,:,:,:)
    complex(dpc), allocatable, private:: opLeft3(:,:,:), opRight3(:,:,:), opLRLeft3(:,:,:,:), opLRRight3(:,:,:,:)
    complex(dpc), allocatable, private:: opPlusLeft3(:,:,:,:), opPlusRight3(:,:,:,:), opMinLeft3(:,:,:,:), opMinRight3(:,:,:,:)
    complex(dpc), allocatable, private:: opLeft1(:,:,:), opPlusLeft1(:,:,:,:), opMinLeft1(:,:,:,:)
    complex(dpc), allocatable, private:: signal(:,:), signalpar(:,:), signalper(:,:)
    real(dp), allocatable, private:: lambda(:), beta(:), LLambda(:), Dtrans(:), Dlong(:)
    complex(dpc), private:: pc, ps, rhosum, mem, orcoeffpar
    complex(dpc), parameter, private:: izero = dcmplx(0.0, 0.0)
    complex(dpc), parameter, private:: iconst = dcmplx(0.0, 1.0)

    integer(i4b), private :: Nhier = 1! number of density matrices in hierarchy
    integer(i4b), private :: tt, nin, s, s2, nnt1, nt1in, nnt3, nnt2, w1, w2

    integer(i4b), allocatable, private :: perm(:,:)
    integer(i4b), allocatable, private :: permplus(:,:), permmin(:,:)

    real(dp), allocatable, private :: CholeskyCF(:,:)
    real(dp), allocatable, private :: oonoise(:)

    integer(i4b), private:: tierstart, tier, kk1, kk2, currentindex, nnn, n, m, nnp, mp, w, wp
    integer(i4b), private:: dir1, dir2, dir3, dir4, pol
    logical, private:: permexists
    integer(i4b), allocatable, private:: currentperm(:)

    ! polarization components
    ! during t1, we need delta = x (1), y (2) and z (3)
    ! during t2, we have delta LLambda = xx (1), xy (2), xz (3), yx (4), yy (5),
    ! yz (6), zx (7), zy (8), zz (9)
    ! during t2, we then need to multiply these to form
    ! delta LLambda beta = xxx (1), xxy (2), xxz(3), xyx (4), xyy(5), xzx(6),
    ! xzz (7), yxx (8), yxy (9), yyx (10), yyy (11), yyz (12), yzy (13),
    ! yzz (14), zxx (15), zxz (16), zyy (17), zyz (18), zzx (19), zzy (20),
    ! zzz (21)

    integer(i4b), parameter, private:: dirinda(21) = (/ 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3 /)
    integer(i4b), parameter, private:: dirindb(21) = (/ 1, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 3 /)
    integer(i4b), parameter, private:: dirindc(21) = (/ 1, 2, 3, 1, 2, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 2, 3, 1, 2, 3 /)
    integer(i4b), parameter, private:: dirindd(21) = (/ 1, 2, 3, 2, 1, 3, 1, 2, 1, 1, 2, 3, 3, 2, 3, 1, 3, 2, 1, 2, 3 /)

    ! xxxx, xxyy, xxzz, xyxy, xyyx, xzxz, xzzx, yxxy, yxyx, yyxx,
    ! yyyy, yyzz, yzyz, yzzy, zxxz, zxzx, zyyz, zyzy, zzxx, zzyy, zzzz

    real(dp), parameter, private:: orcoeffZZYY(21) = (/ 2.0/30, 4.0/30, 4.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, -1.0/30, 4.0/30, 2.0/30, 4.0/30, -1.0/30, -1.0/30,-1.0/30, -1.0/30, -1.0/30, -1.0/30, 4.0/30, 4.0/30, 2.0/30 /)

    real(dp), parameter, private:: orcoeffZZZZ(21) = (/ 6.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 6.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 2.0/30, 6.0/30 /)

    real(dp), parameter, private:: orcoeffZYZY(21) = (/ 2.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 2.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 4.0/30, -1.0/30, 4.0/30, -1.0/30, -1.0/30, 2.0/30 /)


    private::write_phi_config_file
    private::run_phi
    private::fill_evolution_superoperator_hierarchy_phi

    private::arend_init
    private::arend_allocate
    private::arend_deallocate
    private::arend_fill_parameters
    private::arend_init2
    private::arend_benchmark_E
    private::arend_light_hierarchy

    private::arend_initmult1
    private::arend_initmult2
    private::arend_initmult3

    private::arend_propagate1reph
    private::arend_propagate1nonreph
    private::arend_propagate2
    private::arend_propagate3
    private::arend_propagateC

    private::light_CF
    private::close_files
    private::open_files
    !private::arend_mancal_valkunas_quantum_light
    private::arend_mancal_valkunas_quantum_light2
    private::arend_bloch_equations_CW
    private::arend_bloch_equations_noise
    private::read_config_file

    private::generate_noise
    private::init_generate_noise

 	contains

 	complex(dpc) function light_CF(t1,t2) result(res)
 	    real(dp), intent(in)    :: t1, t2
 	    real(dp) :: gamma, beta, TTT

 	    gamma = 1/tau_of_projector
 	    beta = 1/kB_intK/radiation_temperature
 	    TTT = t1 - t2

        if(t1 + 1e-6 > t2 .and. t1 >= -1e-6) then
          if(gaussian_pulse) then
 	        res = exp(-(t1 - t2)*(t1 - t2)*gamma*gamma)
 	      elseif(complex_CF) then
 	        central_frequency = 0
 	        res = exp(-abs(TTT)*gamma)*(1 + sign(1.0_dp,TTT)*cmplx(0,1) * tan(gamma * beta / 2) )
 	      else
 	        res = exp(-abs(t1 - t2)*gamma)
 	      end if
 	        res = res * cos(abs(t1 - t2)*central_frequency)
 	    else
 	        res = 0.0_dp
 	    end if

 	    if(max(t1,t2) > max_light_time) then
 	        res = 0.0_dp
 	    end if
 	end function light_CF

 	subroutine read_config_file
 	    character(len=256)  :: buff = ""
 	    real(dp)            :: value = 0.0_dp
 	    integer(i4b)        :: i = 0

        open(unit=32,file=trim(trim(out_dir)//trim('/../')//trim(external_dir)//'/'//trim(config_filename_arend) ) , err=32, status='old')

        tmax = 6

        do while(i == 0)
          read(32, *, iostat=i) buff, value

          if(trim(adjustl(buff)) == 'tier') then
            write(*,*) buff, int(value)
            tmax = int(value)

          elseif(trim(adjustl(buff)) == 'maxLightTime') then
            write(*,*) buff, value
            max_light_time = value

          elseif(trim(adjustl(buff)) == 'centralFrequency') then
            write(*,*) buff, value/Energy_internal_to_cm, value, 'in cm'
            central_frequency = value/Energy_internal_to_cm

          elseif(trim(adjustl(buff)) == 'globalRelaxRate') then
            write(*,*) buff, value
            global_relax_rate = value

!            ! for closed system, always set to zero (often mistake)
!            if(submethod1 == 'C') then
!              global_relax_rate = 0.0_dp
!              write(*,*) buff, 'changed to 0.00000'
!            end if

          elseif(trim(adjustl(buff)) == 'blochElectricFieldStrength') then
            write(*,*) buff, value
            bloch_strength = value

          elseif(trim(adjustl(buff)) == 'radiationTemperature') then
            write(*,*) buff, value
            radiation_temperature = value

          elseif(trim(adjustl(buff)) == 'gaussianPulse') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              gaussian_pulse = .true.
            else
              gaussian_pulse = .false.
            end if

          elseif(trim(adjustl(buff)) == 'complexCF') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              complex_CF = .true.
            else
              complex_CF = .false.
            end if

          elseif(trim(adjustl(buff)) == 'normalizeTrace') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              normalize_trace = .true.
            else
              normalize_trace = .false.
            end if

          elseif(trim(adjustl(buff)) == 'blochTerm') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              bloch_term = .true.
            else
              bloch_term = .false.
            end if

          elseif(trim(adjustl(buff)) == 'lightHierarchy') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              light_hierarchy = .true.
            else
              light_hierarchy = .false.
            end if

          elseif(trim(adjustl(buff)) == 'noiseTerm') then
            write(*,*) buff, int(value)
            if(int(value) == 1) then
              noise_term = .true.
            else
              noise_term = .false.
            end if

          elseif(trim(adjustl(buff)) == 'noiseRealizations') then
            write(*,*) buff, int(value)
            realizations = int(value)

          else
            write(*,*) 'unrecognised:', buff, value
          end if
        end do

        close(32)

        return

32      call print_warning_message("couldn't read the supplementary config file, using default values" ,5)
 	end subroutine read_config_file

    subroutine fill_evolution_superoperator_hierarchy(type)
        character, intent(in) :: type
        integer(i4b) :: i, j

        call print_log_message("fill_evolution_superoperator_hierarchy called",5)
        if(submethod1 == 'E') then
        else
            call arend_main()
        end if
        !stop
    end subroutine fill_evolution_superoperator_hierarchy

    subroutine arend_main()
      character(len=128) :: buff
      integer(i4b)       :: b

      call read_config_file()
      call arend_init()
      call arend_allocate()
      call arend_fill_parameters()
      call arend_init2()

      !call arend_benchmark_E(1,2)

      !do b=1, Nsys
      !  call arend_benchmark_O(s)
      !end do

      write(buff,*) 'Submethod used:',submethod1
      call print_log_message(adjustl(trim(buff)), 5)

      if(bloch_term) then
        write(*,*) central_frequency, rwa,central_frequency-rwa
        call arend_bloch_equations_CW(submethod1)
      elseif(noise_term) then
        call arend_bloch_equations_noise(submethod1)
      elseif(light_hierarchy) then
        call arend_light_hierarchy(submethod1)
      else
        call arend_mancal_valkunas_quantum_light2(submethod1)
      end if



      call arend_deallocate()

    end subroutine arend_main

    subroutine arend_mancal_valkunas_quantum_light2(type)
      character, intent(in)  :: type
      character(len=128)     :: buff, buff2
      complex(dpc), dimension(0:Nsys, 0:Nsys)              :: rhotmp
      complex(dpc), dimension(0:Nsys, 0:Nsys, Ntimestept2) :: rho_physical
      complex(dpc), dimension(Nhier+1, 0:Nsys, 0:Nsys)     :: rhotmp2
      real(dp)     :: time1, time2, time
      complex(dpc) :: intFactor
      integer(i4b) :: nnt

      call arend_initmult1()
      call arend_initmult2()
      rho_physical = 0.0_dp

      ! XXXXXXXXXXXXXXXx
        dir1 = 1
        dir2 = 1
        dir3 = 1
        dir4 = 1
      ! XXXXXXXXXXXXXXXx

      write(*,*) 'mu_x =',mu(:,1)
      write(*,*) 'mu_y =',mu(:,2)
      write(*,*) 'mu_z =',mu(:,3)

      rhoC = 0.0_dp

      ! Initial condition
      do nnn = 1, 1!Nhier
        do s=1, Nsys
          rhoC(nnn, s, 0) = mu(s, dir1)
        end do
      end do

      do nnt1 = 1, Ntimestept1
        if(submethod1 == 'I' .and. nnt1 > 1) then
          exit
        end if

        time1 = (nnt1-1)*Ntimestept1in*dt
        intFactor = ( exp( -dt*Ntimestept1in*rwa*cmplx(0,1) ) -  1.0_dp )/rwa*cmplx(0.0,1.0) / (dt*Ntimestept1in)

        write(*,*) '   nnt1 = ', nnt1
        write(*,*) rhoC(1, 1:, 0)
        call flush()

        ! store the time evolution
        rhotmp2 = rhoC

        ! initialize t2
        do nnn = 1, 1!Nhier
          do s = 1, Nsys
            do s2 = 1, Nsys
              rhoC(nnn, s, s2) = rhoC(nnn, s, s2)  + mu(s, dir2) * rhoC(nnn, s2, 0)* exp( -time1*rwa*cmplx(0,1) ) * intFactor          &
                                                   + conjg(mu(s2, dir2) * rhoC(nnn, s, 0) * exp( -time1*rwa*cmplx(0,1) ) ) * intFactor

              if(nnt == 1) then
                rhoC(nnn, s, s2) = rhoC(nnn, s, s2) / 2
              end if
            end do
          end do
        end do

        ! propagate t2
        do nnt2 = 1, Ntimestept2-nnt1
          time2 = (nnt2-1)*Ntimestept2in*dt

          do s = 1, Nsys
          do s2 = 1, Nsys

          if(submethod2 == 'D') then
            nnt = nnt1+nnt2-1
            time = (nnt-1)*Ntimestept1in*dt
            rho_physical(s,s2,nnt) = rho_physical(s,s2,nnt) + rhoC(1,s,s2)*light_CF(time - time2, time - time2 - time1)*dt*Ntimestept1in
          else
            do nnt = max(nnt2+nnt1-1, 1), Ntimestept2
              time = (nnt-1)*Ntimestept1in*dt
              rho_physical(s,s2,nnt) = rho_physical(s,s2,nnt) + rhoC(1,s,s2)*light_CF(time - time2, time - time2 - time1)*dt*Ntimestept1in*dt*Ntimestept1in
            end do
          end if

          end do
          end do

          do nin = 1, Ntimestept2in
            call arend_propagateC(dt, time1+time2)
          end do
        end do

                  if(mod(nnt1,10) == 1 .or. Ntimestept1 == nnt1) then
                  call open_files()

                  ! print outcome
                  do nnt=1,Ntimestept2
                      rhotmp(:,:) = rho_physical(:,:,nnt) + transpose(conjg(rho_physical(:,:,nnt)))

                      if(exciton_basis) then
                        call operator_to_exc(rhotmp(1:,1:),'E')
                        !call operator_to_exc(rhotmp(0,1:),'O')
                        !call operator_to_exc(rhotmp(1:,0),'O')
                      end if

                      if(normalize_trace) then
                        rhotmp = rhotmp / (trace(rhotmp) + 1e-100_dp)
                      end if

                      do s=1, Nsys
                      do s2=1, Nsys
                        write(s+Nsys*s2+10,*) dt*Ntimestept2in*(nnt), real(rhotmp(s,s2)), aimag(rhotmp(s,s2))
                      end do
                      end do

                      if(maxval(abs(rhotmp)) > 1e-6 ) then
                        write(10,*) dt*Ntimestept2in*(nnt), entropy(rhotmp/trace(rhotmp) )
                      end if
                  end do

                  call close_files()
                  end if

        ! restore the time evolution
        rhoC = rhotmp2

        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagateC(dt, time1)
        end do
      end do ! over nnt1

    end subroutine arend_mancal_valkunas_quantum_light2


!    subroutine arend_mancal_valkunas_quantum_light(type)
!      character, intent(in)  :: type
!      character(len=128)     :: buff, buff2
!      complex(dpc), dimension(Nsys, Nsys)              :: rhotmp
!      complex(dpc), dimension(Nsys, Nsys, Ntimestept2) :: rho_physical
!      integer(i4b) :: nnt
!
!      call arend_initmult1()
!      call arend_initmult2()
!      rho_physical = 0.0_dp
!
!      ! XXXXXXXXXXXXXXXx
!        dir1 = 1
!        dir2 = 1
!        dir3 = 1
!        dir4 = 1
!      ! XXXXXXXXXXXXXXXx
!
!      write(*,*) 'mu_x =',mu(:,1)
!      write(*,*) 'mu_y =',mu(:,2)
!      write(*,*) 'mu_z =',mu(:,3)
!
!      ! Initial condition
!      do nnn = 1, 1!Nhier
!        do s=1, Nsys
!          rho1(nnn, s) = mu(s, dir1)
!        end do
!      end do
!
!      do nnt1 = 1, Ntimestept1
!        if(submethod1 == 'I' .and. nnt1 > 1) then
!          exit
!        end if
!
!        ! initialize t2
!        rho2 = 0.0_dp
!        do nnn = 1, 1!Nhier
!          do s = 1, Nsys
!            do s2 = 1, Nsys
!              ! check the sign of rwa here
!              rho2(nnn, s, s2) = rho2(nnn, s, s2) + mu(s, dir2) * rho1(nnn, s2)* exp( -nnt1*Ntimestept1in*dt*rwa*cmplx(0,1) ) &
!                                                  + conjg(mu(s2, dir2) * rho1(nnn, s) * exp( -nnt1*Ntimestept1in*dt*rwa*cmplx(0,1) ) )
!            end do
!          end do
!        end do
!
!        !write(*,*) '   nnt1 = ', nnt1
!        call flush()
!
!        write(*,*) rho1(1, :)
!
!        ! propagate t2
!        do nnt2 = 1, Ntimestept2-nnt1
!          do s = 1, Nsys
!          do s2 = 1, Nsys
!
!          if(submethod2 == 'D') then
!            nnt = nnt1+nnt2-1
!            rho_physical(s,s2,nnt) = rho_physical(s,s2,nnt) + rho2(1,s,s2)*light_CF((nnt-nnt2)*dt*Ntimestept1in, (nnt-nnt1-nnt2+1)*dt*Ntimestept1in)*dt*Ntimestept1in
!          else
!            do nnt = max(nnt2+nnt1-1, 1), Ntimestept2
!              rho_physical(s,s2,nnt) = rho_physical(s,s2,nnt) + rho2(1,s,s2)*light_CF((nnt-nnt2)*dt*Ntimestept1in, (nnt-nnt1-nnt2+1)*dt*Ntimestept1in)*dt*Ntimestept1in
!            end do
!          end if
!
!          end do
!          end do
!
!          !! XXXXXXXXXX
!          !rho_physical(:,:,nnt2) = rho2(1,:,:)
!
!          !if(real(rho2(1,1,1)) < 0) then
!          !  write(*,*) nnt1, nnt2, rho2(1,1,1), abs(rho1(1, :))
!          !end if
!
!          do nin = 1, Ntimestept2in
!            call arend_propagate2(dt)
!          end do
!        end do
!
!                  if(mod(nnt1,10) == 1 .or. Ntimestept1 == nnt1) then
!                  call open_files()
!
!                  ! print outcome
!                  do nnt=1,Ntimestept2
!                      rhotmp(:,:) = rho_physical(:,:,nnt)! + transpose(conjg(rho_physical(:,:,nnt)))
!                      if(exciton_basis) then
!                        call operator_to_exc(rhotmp(:,:),'E')
!                      end if
!                      do s=1, Nsys
!                      do s2=1, Nsys
!                        write(s+Nsys*s2+10,*) dt*Ntimestept2in*(nnt-1), real(rhotmp(s,s2)), aimag(rhotmp(s,s2))
!                      end do
!                      end do
!
!                      if(maxval(abs(rhotmp)) > 1e-6 ) then
!                        write(10,*) dt*Ntimestept2in*(nnt-1), entropy(rhotmp/trace(rhotmp))
!                      end if
!                  end do
!
!                  call close_files()
!                  end if
!
!        ! propagate t1
!        do nin = 1, Ntimestept1in
!          call arend_propagate1reph(dt)
!        end do
!      end do ! over nnt1
!
!    end subroutine arend_mancal_valkunas_quantum_light

    subroutine arend_bloch_equations_CW(type)
      character, intent(in)  :: type
      character(len=128)     :: buff, buff2
      complex(dpc), dimension(0:Nsys, 0:Nsys)              :: rhotmp
      complex(dpc), dimension(Nhier+1, 0:Nsys, 0:Nsys)     :: rhotmp2
      real(dp)     :: time1, time2, time

      call arend_initmult1()
      call arend_initmult2()
      bloch_term = .true.

      ! XXXXXXXXXXXXXXXx
        dir1 = 1
        dir2 = 1
        dir3 = 1
        dir4 = 1
      ! XXXXXXXXXXXXXXXx

      write(*,*) 'mu_x =',mu(:,1)
      write(*,*) 'mu_y =',mu(:,2)
      write(*,*) 'mu_z =',mu(:,3)

      rhoC = 0.0_dp

      ! Initial condition -- ground state
      do nnn = 1, 1!Nhier
        do s=1, Nsys
          rhoC(nnn, 0, 0) = 1
        end do
      end do

      call open_files()

      do nnt1 = 1, Ntimestept1
        time1 = (nnt1-1)*Ntimestept1in*dt

        write(*,*) '   nnt1 = ', nnt1
        write(*,*) rhoC(1, 1:, 0)
        call flush()

                  ! print outcome
                      rhotmp(:,:) = rhoC(1,:,:) !+ transpose(conjg(rho_physical(:,:,nnt)))

                      if(exciton_basis) then
                        call operator_to_exc(rhotmp(1:,1:),'E')
                        !call operator_to_exc(rhotmp(0,1:),'O')
                        !call operator_to_exc(rhotmp(1:,0),'O')
                      end if

                      if(normalize_trace) then
                        rhotmp = rhotmp / (trace(rhotmp(1:,1:)) + 1e-100_dp)
                      end if

                      do s=1, Nsys
                      do s2=1, Nsys
                        write(s+Nsys*s2+10,*) time1, real(rhotmp(s,s2)), aimag(rhotmp(s,s2))
                      end do
                      end do

                      if(maxval(abs(rhotmp)) > 1e-6 ) then
                        write(10,*) time1, entropy(rhotmp/trace(rhotmp) )
                      end if



        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagateC(dt, time1 + (nin-1)*dt)
        end do
      end do ! over nnt1

      call close_files()

    end subroutine arend_bloch_equations_CW

    subroutine arend_light_hierarchy(type)
      character, intent(in)  :: type
      character(len=128)     :: buff, buff2
      complex(dpc), dimension(Nsys, Nsys)              :: rhotmp
      complex(dpc), dimension(Nhier+1, Nsys, Nsys)     :: rhotmp2
      real(dp)     :: time1, time2, time
      integer(i4b) :: aa,bb

      call arend_initmult1()
      call arend_initmult2()
      light_hierarchy = .true.

      ! XXXXXXXXXXXXXXXx
        dir1 = 1
        dir2 = 1
        dir3 = 1
        dir4 = 1
      ! XXXXXXXXXXXXXXXx

      write(*,*) 'mu_x =',mu(:,1)
      write(*,*) 'mu_y =',mu(:,2)
      write(*,*) 'mu_z =',mu(:,3)

      do aa = 1, Nsys
      do bb = 1, Nsys
        write(*,*) 'VV',V(aa,:,bb)
      end do
      write(*,*)
      end do

      rho2 = 0.0_dp

      ! Initial condition -- ground state
      do nnn = 1, 1!Nhier
        rho2(nnn, Nsys, Nsys) = 1
      end do

      call open_files()

      do nnt1 = 1, Ntimestept1
        time1 = (nnt1-1)*Ntimestept1in*dt

        write(*,*) '   nnt1 = ', nnt1
        write(*,*) rho2(1, 1, :)
        write(*,*) rho2(1, 2, :)
        write(*,*) rho2(1, Nsys, :)
        call flush()

                  ! print outcome
                      rhotmp(:,:) = rho2(1,:,:) !+ transpose(conjg(rho_physical(:,:,nnt)))

                      if(exciton_basis) then
                        call operator_to_exc(rhotmp(1:(Nsys-1),1:(Nsys-1)),'E')
                        !call operator_to_exc(rhotmp(0,1:),'O')
                        !call operator_to_exc(rhotmp(1:,0),'O')
                      end if

                      if(normalize_trace) then
                        rhotmp = rhotmp / (trace(rhotmp(1:(Nsys-1),1:(Nsys-1))) + 1e-100_dp)
                      end if

                      do s=1, Nsys-1
                      do s2=1, Nsys-1
                        write(s+Nsys*s2+10,*) time1, real(rhotmp(s,s2)), aimag(rhotmp(s,s2))
                      end do
                      end do

                      if(maxval(abs(rhotmp(1:(Nsys-1),1:(Nsys-1)))) > 1e-6 ) then
                        write(10,*) time1, entropy(rhotmp(1:(Nsys-1),1:(Nsys-1))/trace(rhotmp(1:(Nsys-1),1:(Nsys-1))) )
                      end if



        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagate2(dt)
        end do
      end do ! over nnt1

      call close_files()

    end subroutine arend_light_hierarchy


    subroutine arend_bloch_equations_noise(type)
      character, intent(in)  :: type
      character(len=128)     :: buff, buff2
      complex(dpc), dimension(0:Nsys, 0:Nsys)              :: rhotmp
      complex(dpc), dimension(0:Nsys, 0:Nsys, Ntimestept1) :: rho_physical
      complex(dpc), dimension(Nhier+1, 0:Nsys, 0:Nsys)     :: rhotmp2
      real(dp)     :: time1, time2, time
      integer(i4b) :: nnoise
      logical      :: write_this_loop


      ALLOCATE(oonoise, (NtimestepT1*NtimestepT1In))

      call arend_initmult1()
      call arend_initmult2()
      rho_physical = 0.0_dp
      bloch_term = .false.
      noise_term = .true.

      call init_generate_noise()

      ! XXXXXXXXXXXXXXXx
        dir1 = 1
        dir2 = 1
        dir3 = 1
        dir4 = 1
      ! XXXXXXXXXXXXXXXx

      write(*,*) 'mu_x =',mu(:,1)
      write(*,*) 'mu_y =',mu(:,2)
      write(*,*) 'mu_z =',mu(:,3)

    do nnoise=1, realizations
      write(*,*) 'realization', nnoise

      if(mod(nnoise, 250) == 1) then
        write_this_loop = .true.
      else
        write_this_loop = .false.
      end if

      call generate_noise(oonoise)

      rhoC = 0.0_dp

      ! Initial condition -- ground state
      do nnn = 1, 1!Nhier
        do s=1, Nsys
          rhoC(nnn, 0, 0) = 1
        end do
      end do

      if(write_this_loop) &
        call open_files()

      do nnt1 = 1, Ntimestept1
        time1 = (nnt1-1)*Ntimestept1in*dt

        !write(*,*) '   nnt1 = ', nnt1
        !write(*,*) rhoC(1, 1:, 0)
        call flush()

                  ! print outcome
                      rhotmp(:,:) = rhoC(1,:,:) !+ transpose(conjg(rho_physical(:,:,nnt)))

                      if(exciton_basis) then
                        call operator_to_exc(rhotmp(1:,1:),'E')
                        !call operator_to_exc(rhotmp(0,1:),'O')
                        !call operator_to_exc(rhotmp(1:,0),'O')
                      end if

                      if(normalize_trace) then
                        rhotmp = rhotmp / (trace(rhotmp(1:,1:)) + 1e-100_dp)
                      end if

                      rho_physical(:,:,nnt1) = rho_physical(:,:,nnt1) + rhotmp(:,:)

                      if(write_this_loop) then
                      do s=1, Nsys
                      do s2=1, Nsys
                        write(s+Nsys*s2+10,*) time1, real(rho_physical(s,s2,nnt1)), aimag(rho_physical(s,s2,nnt1))
                      end do
                      end do

                      if(maxval(abs(rho_physical(:,:,nnt1))) > 1e-6 ) then
                        write(10,*) time1, entropy(rho_physical(:,:,nnt1)/trace(rho_physical(:,:,nnt1)) )
                      end if
                      end if



        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagateC(dt, time1)
        end do
      end do ! over nnt1

      if(write_this_loop) &
        call close_files()

    end do ! over nnoise

      DEALLOCATE(oonoise)

    end subroutine arend_bloch_equations_noise

    subroutine arend_benchmark_E(i0,j0)
      integer(i4b), intent(in) :: i0, j0

      character(len=128) :: buff, buff2
      complex(dpc), dimension(Nsys, Nsys) :: rhotmp

      do s=1, Nsys
      do s2=1, Nsys
        write(buff, '(I2)') s
        buff = repeat( '0', max(2-len_trim(adjustl(buff)), 0)  ) // adjustl(buff)

        write(buff2, '(I2)') s2
        buff2 = repeat( '0', max(2-len_trim(adjustl(buff2)), 0)  ) // adjustl(buff2)

        buff = 'rhoE'//trim(buff)//'-'//trim(buff2)//'.dat'
        call flush()
        open(unit=s+Nsys*s2+10,file=trim(file_join(out_dir,adjustl(trim(buff)))))
      end do
      end do

      call arend_initmult1()
      call arend_initmult2()
      !call arend_initmult3()

      ! Initial condition
      do nnn = 1, 1!Nhier
        do s=1, Nsys
                              if(s == i0) then
                                rho1(nnn, s) = 1.0
                              else
                                rho1(nnn, s) = 0.0
                              end if
          !rho1(nnn, s) = mu(s, dir1)
        end do
                              if(exciton_basis) then
                                call operator_from_exc(rho1(nnn,:),'O')
                              end if
      end do


      do nnt1 = 2, Ntimestept1
        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagate1reph(dt)
        end do
      end do ! over nnt1

      ! initialize t2
      rho2 = 0.0_dp
      do nnn = 1, 1!Nhier
                              if(exciton_basis) then
                                call operator_to_exc(rho1(nnn,:),'O')
                              end if
        do s = 1, Nsys
          do s2 = 1, Nsys
                              if(s == i0 .and. s2 == j0) then
                                rho2(nnn,s,s2) = rho1(nnn, s)
                              else
                                rho2(nnn,s,s2) = 0.0
                              end if
            !rho2(nnn, s, s2) = mu(s, dir2) * rho1(nnn, s2)
          end do
        end do
                              if(exciton_basis) then
                                call operator_from_exc(rho1(nnn,:),'O')
                                call operator_from_exc(rho2(nnn,:,:),'E')
                              end if
      end do

      ! propagate t2
      do nnt2 = 1, Ntimestept2
        rhotmp(:,:) = rho2(1,:,:)
        if(exciton_basis) then
          !call operator_to_exc(rho2(1,:,:),'E')
          call operator_to_exc(rhotmp(:,:),'E')
        end if
        do s=1, Nsys
        do s2=1, Nsys
          write(s+Nsys*s2+10,*) dt*Ntimestept2in*(nnt2-1)/1e15, real(rhotmp(s,s2)), aimag(rhotmp(s,s2))
        end do
        end do
        !if(exciton_basis) then
        !  call operator_from_exc(rho2(1,:,:),'E')
        !end if


        do nin = 1, Ntimestept2in
          call arend_propagate2(dt)
        end do
      end do

      do s=1, Nsys
      do s2=1, Nsys
        close(s+Nsys*s2+10)
      end do
      end do
    end subroutine arend_benchmark_E

    subroutine arend_benchmark_O(i0)
      integer(i4b), intent(in) :: i0

      character(len=128) :: buff

      do s=1, Nsys
        write(buff, '(I2)') s
        buff = repeat( '0', max(2-len_trim(adjustl(buff)), 0)  ) // adjustl(buff)

        buff = 'rhoO'//trim(buff)//'-G.dat'
        open(unit=s+10,file=trim(file_join(out_dir,adjustl(trim(buff)))))
      end do

      call arend_initmult1()
      !call arend_initmult2()
      !call arend_initmult3()

      ! Initial condition
      do nnn = 1, 1!Nhier
        do s=1, Nsys
                              if(s == i0) then
                                rho1(nnn, s) = 1.0
                              else
                                rho1(nnn, s) = 0.0
                              end if
          !rho1(nnn, s) = mu(s, dir1)
        end do
      end do


      do nnt1 = 1, Ntimestept1
        write(*,'(A6,I2,A12,I5,A3,I5)') "pol = ", pol, " / 21; t1 = ", nnt1, " / " , Ntimestept1
        call flush()

        do s=1, Nsys
            write(s+10,*) dt*Ntimestept2in*(nnt1-1)/1e15, real(rho1(1,s)) , aimag(rho1(1,s))

            if(mod(nnt1, Ntimestept1in) == 1) then
              evops(1,1)%Ueg(i0,1, s,1, INT((nnt1-1)/Ntimestept1in) + 1) = rho1(1,s)
            end if
        end do

        ! propagate t1
        do nin = 1, Ntimestept1in
          call arend_propagate1reph(dt)
        end do

      end do ! over nnt1

      do s=1, Nsys
        close(s+10)
      end do
    end subroutine arend_benchmark_O


    subroutine open_files()
      character(len=128)     :: buff, buff2

      do s=1, Nsys
      do s2=1, Nsys
        write(buff, '(I2)') s
        buff = repeat( '0', max(2-len_trim(adjustl(buff)), 0)  ) // adjustl(buff)

        write(buff2, '(I2)') s2
        buff2 = repeat( '0', max(2-len_trim(adjustl(buff2)), 0)  ) // adjustl(buff2)

        buff = 'rhoE'//trim(buff)//'-'//trim(buff2)//'.dat'
        call flush()
        open(unit=s+Nsys*s2+10,file=trim(file_join(out_dir,adjustl(trim(buff)))))
      end do
      end do

      ! entropy file
      buff = 'rhoEntropy.dat'
      call flush()
      open(unit=10,file=trim(file_join(out_dir,adjustl(trim(buff)))))
    end subroutine open_files

    subroutine close_files()
      do s=1, Nsys
      do s2=1, Nsys
        close(s+Nsys*s2+10)
      end do
      end do

      !entropy file
      close(10)
    end subroutine close_files

    subroutine arend_init()
      character(len=256) :: buff

      if(light_hierarchy) then ! include ground state into system
        Nsys = N1_from_type('E') + 1
      else
        Nsys = N1_from_type('E')
      end if
      Nind = Nsys

      do tt = 1, tmax
        Nhier = Nhier + numpermt(tt)
      end do

      Ntimestept1 = Nt(1)
      Ntimestept2 = Nt(2)
      Ntimestept3 = Nt(3)
      Ntimestept1in = gt(1)
      Ntimestept2in = gt(2)
      Ntimestept3in = gt(3)

      if(Ntimestept1in /= Ntimestept2in) then
        call print_error_message(-1,'Ntimestept1in /= Ntimestept2in, internal error')
      end if


      write(buff,*) 'number of elements in hierarchy = ', Nhier
      call print_log_message(buff,5)
    end subroutine arend_init

    subroutine arend_allocate()
      ALLOCATE(rho1,(Nhier+1, Nsys)) ! +1: store zeros in the last density matrix
      ALLOCATE(prhodx1,(Nhier+1, Nsys))
      ALLOCATE(prhox1,(Nhier+1, Nsys))
      rho1 = 0.0
      prhodx1 = 0.0
      prhox1 = 0.0

      ALLOCATE(rho2,(Nhier+1, Nsys, Nsys)) ! +1: store zeros in the last density matrix
      ALLOCATE(prhodx2,(Nhier+1, Nsys, Nsys))
      ALLOCATE(prhox2,(Nhier+1, Nsys, Nsys))
      rho2 = 0.0
      prhodx2 = 0.0
      prhox2 = 0.0

      ALLOCATE(rhoC,(Nhier+1, 0:Nsys, 0:Nsys))     ! +1: store zeros in the last density matrix
      ALLOCATE(prhodxC,(Nhier+1, 0:Nsys, 0:Nsys))  ! 0: -- include ground state at zero index
      ALLOCATE(prhoxC,(Nhier+1, 0:Nsys, 0:Nsys))
      rhoC = 0.0
      prhodxC = 0.0
      prhoxC = 0.0

      ALLOCATE(rho3,(Nhier+1, (Nsys*(Nsys-1))/2, Nsys)) ! +1: store zeros in the last density matrix
      ALLOCATE(prhodx3,(Nhier+1, (Nsys*(Nsys-1))/2, Nsys))
      ALLOCATE(prhox3,(Nhier+1, (Nsys*(Nsys-1))/2, Nsys))
      rho3 = 0.0
      prhodx3 = 0.0
      prhox3 = 0.0

      ALLOCATE(rho3s,(Nhier+1, Nsys)) ! +1: store zeros in the last density matrix
      ALLOCATE(prhodx3s,(Nhier+1, Nsys))
      ALLOCATE(prhox3s,(Nhier+1, Nsys))
      rho3s = 0.0
      prhodx3s = 0.0
      prhox3s = 0.0

      ALLOCATE(HS,(Nsys, Nsys))
      Hs = 0.0

      ALLOCATE(V,(Nsys, Nsys, Nsys))
      V = 0.0

      ALLOCATE(HS2,((Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      HS2 = 0.0

      ALLOCATE(V2,(Nsys, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      V2 = 0.0


      ALLOCATE(beta,(Nsys))
      ALLOCATE(LLambda,(Nsys))
      ALLOCATE(lambda,(Nsys))
      ALLOCATE(Dlong,(Nsys))
      ALLOCATE(Dtrans,(Nsys))


      ALLOCATE(perm,(Nhier, Nind))
      ALLOCATE(currentperm,(Nind))

      ALLOCATE(opLeft1,(Nhier, Nsys, Nsys))
      ALLOCATE(opPlusLeft1,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opMinLeft1,(Nhier, Nsys, Nsys, Nsys))


      ALLOCATE(opLeft2,(Nhier, Nsys, Nsys))
      ALLOCATE(opRight2,(Nhier, Nsys, Nsys))

      ALLOCATE(opLRLeft2,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opLRRight2,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opPlusLeft2,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opPlusRight2,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opMinLeft2,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opMinRight2,(Nhier, Nsys, Nsys, Nsys))

      ALLOCATE(opLeft3,(Nhier, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      ALLOCATE(opRight3,(Nhier, Nsys, Nsys))
      ALLOCATE(opLRLeft3,(Nhier, Nsys, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      ALLOCATE(opLRRight3,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opPlusLeft3,(Nhier, Nsys, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      ALLOCATE(opPlusRight3,(Nhier, Nsys, Nsys, Nsys))
      ALLOCATE(opMinLeft3,(Nhier, Nsys, (Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))
      ALLOCATE(opMinRight3,(Nhier, Nsys, Nsys, Nsys))

      ALLOCATE(permplus,(Nhier, Nsys))

      ALLOCATE(permmin,(Nhier, Nsys))

      ALLOCATE(signal,(Ntimestept1, Ntimestept3))
      ALLOCATE(signalpar,(Ntimestept1, Ntimestept3))
      ALLOCATE(signalper,(Ntimestept1, Ntimestept3))

      ALLOCATE(mu,(Nsys, 3))
    end subroutine arend_allocate

    subroutine arend_deallocate()
      DEALLOCATE(rho1)
      DEALLOCATE(prhodx1)
      DEALLOCATE(prhox1)

      DEALLOCATE(rho2)
      DEALLOCATE(prhodx2)
      DEALLOCATE(prhox2)

      DEALLOCATE(rhoC)
      DEALLOCATE(prhodxC)
      DEALLOCATE(prhoxC)

      DEALLOCATE(rho3)
      DEALLOCATE(prhodx3)
      DEALLOCATE(prhox3)

      DEALLOCATE(rho3s)
      DEALLOCATE(prhodx3s)
      DEALLOCATE(prhox3s)

      DEALLOCATE(HS)

      DEALLOCATE(V)

      DEALLOCATE(HS2)

      DEALLOCATE(V2)

      DEALLOCATE(beta)
      DEALLOCATE(LLambda)
      DEALLOCATE(lambda)
      DEALLOCATE(Dlong)
      DEALLOCATE(Dtrans)

      DEALLOCATE(perm)
      DEALLOCATE(currentperm)

      DEALLOCATE(opLeft1)
      DEALLOCATE(opPlusLeft1)
      DEALLOCATE(opMinLeft1)

      DEALLOCATE(opLeft2)
      DEALLOCATE(opRight2)

      DEALLOCATE(opLRLeft2)
      DEALLOCATE(opLRRight2)
      DEALLOCATE(opPlusLeft2)
      DEALLOCATE(opPlusRight2)
      DEALLOCATE(opMinLeft2)
      DEALLOCATE(opMinRight2)

      DEALLOCATE(opLeft3)
      DEALLOCATE(opRight3)
      DEALLOCATE(opLRLeft3)
      DEALLOCATE(opLRRight3)
      DEALLOCATE(opPlusLeft3)
      DEALLOCATE(opPlusRight3)
      DEALLOCATE(opMinLeft3)
      DEALLOCATE(opMinRight3)

      DEALLOCATE(permplus)

      DEALLOCATE(permmin)

      DEALLOCATE(signal)
      DEALLOCATE(signalpar)
      DEALLOCATE(signalper)

      DEALLOCATE(mu)
    end subroutine arend_deallocate


    subroutine arend_fill_parameters()
      character(len=256) :: buff

      ! parameters for system-bath coupling
      do s = 1, Nsys
        if(s == Nsys .and. light_hierarchy) then
          lambda(s) = bloch_strength
          LLambda(s) = 1.0/tau_of_projector
          beta(s) = 1/kB_intK/radiation_temperature
          Dlong(s) = 0.0_dp
          Dtrans(s) = 0.0_dp
        else
          lambda(s) = igofts(iblocks(1,1)%sblock%gindex(s))%goft%params(1,1)
          LLambda(s) = 1/igofts(iblocks(1,1)%sblock%gindex(s))%goft%params(2,1)
          beta(s) = 1/kB_intK/temp
          Dlong(s) = 1.0_dp
          Dtrans(s) = 0.0_dp

          if(submethod1 == 'C') then
            Dlong(s) = 0.0_dp
          end if
        end if

        write(buff,*) ';lambda  ', lambda(s)*Energy_internal_to_cm
        call print_log_message(adjustl(trim(buff)),5)
        write(buff,*) ';LLambda ', LLambda(s)
        call print_log_message(adjustl(trim(buff)),5)
        write(buff,*) ';beta    ', 1/beta(s)*Energy_internal_to_cm/0.69503568_dp
        call print_log_message(adjustl(trim(buff)),5)
        write(buff,*) ';Dlong   ', Dlong(s)
        call print_log_message(adjustl(trim(buff)),5)
        write(buff,*) ';Dtrans  ', Dtrans(s)
        call print_log_message(adjustl(trim(buff)),5)
      end do


      ! system Hamiltonian
      HS = 0.0_dp
      do s=1, Nsys
      do s2=1, Nsys
        if(max(s,s2) == Nsys .and. light_hierarchy) then
          if(s == s2) then
            HS(s,s) = - rwa
          end if
        else
          if(s == s2) then
            HS(s,s) = iblocks(1,1)%sblock%en(s) - rwa
          else
            HS(s,s2) = iblocks(1,1)%sblock%J(s,s2)
          end if
        end if
      end do
      end do

      ! transition dipoles, (Hilbert space index, xyz-index)
      mu = 0.0_dp
      do s=1, Nsys
        if(light_hierarchy .and. s == Nsys) then
          cycle
        end if

          mu(s,1) = current_s_block%dx(s,1)
          mu(s,2) = current_s_block%dy(s,1)
          mu(s,3) = current_s_block%dz(s,1)
      end do

      ! system-bath coupling, no correlation
      do s=1, Nsys
        V(s,s,s) = dcmplx(Dlong(s))
      end do

      if(light_hierarchy) then
        do s=1, Nsys-1
        do s2=1, Nsys-1
          V(Nsys,s,Nsys) = mu(s,1)
          V(Nsys,Nsys,s) = mu(s,1)
        end do
        end do
      end if

!      ! 2-quantum Hamiltonian and system-bath coupling operators
!      do n = 1, Nsys
!        do m = (n+1), Nsys
!          w = (n-1)*Nsys + m - (n*(n+1))/2
!
!          do nnp = 1, Nsys
!            do mp = (nnp+1), Nsys
!              wp = (nnp-1)*Nsys + mp - (nnp*(nnp+1))/2
!
!          if (n == nnp) then
!                 HS2(w, wp) = HS2(w, wp) + HS(m, mp)
!                 do s = 1, Nsys
!               V2(s, w, wp) = V2(s, w, wp) + V(s, m, mp)
!                 end do
!              end if
!          if (n == mp) then
!                 HS2(w, wp) = HS2(w, wp) + HS(m, nnp)
!                 do s = 1, Nsys
!               V2(s, w, wp) = V2(s, w, wp) + V(s, m, nnp)
!                 end do
!              end if
!          if (m == nnp) then
!                 HS2(w, wp) = HS2(w, wp) + HS(n, mp)
!                 do s = 1, Nsys
!               V2(s, w, wp) = V2(s, w, wp) + V(s, n, mp)
!                 end do
!              end if
!          if (m == mp) then
!                 HS2(w, wp) = HS2(w, wp) + HS(n, nnp)
!                 do s = 1, Nsys
!               V2(s, w, wp) = V2(s, w, wp) + V(s, n, nnp)
!                 end do
!              end if
!
!            end do
!           end do
!        end do
!      end do

    end subroutine arend_fill_parameters

    subroutine arend_init2()
      character(len=128) :: buff

      ! build index
      tierstart = 1 ! first element of the current tier
      currentindex = 2

      currentperm = 0.0
      perm = 0.0

      do tier = 1, tmax
        tierstart = tierstart + numpermt(tier-1)
        do kk1 = currentindex-numpermt(tier-1), currentindex-1 ! loop over all elements in the previous tier
          do kk2 = 1, Nind ! try to add a ball at position kk2
              currentperm(:) = perm(kk1, :)
              currentperm(kk2) = currentperm(kk2) + 1
              permexists = .False.
              do nnn = tierstart, currentindex-1 ! check if it is already in the list
                 if ( ALL(perm(nnn, :) == currentperm) ) then
                    permexists = .True.
                    !write(*,*) "permutation ", currentperm, " exists"
                 end if
              end do
              if (permexists .eqv. .False.) then ! if not, add it to the list
                perm(currentindex, :) = currentperm(:)
                currentindex = currentindex + 1
                !write(*,*) "new permutation ", currentperm
              end if

           end do
         end do
      end do

      !do s = 1, size(perm,1)
      !  do nnn = 1, size(perm,2)
      !    write(*,'(I2A)', advance='no') perm(s,nnn),' '
      !  end do
      !  write(*,*)
      !end do
      !write(*,*)

      ! cache plus and minus
      do nnn = 1, Nhier
        do s = 1, Nsys
          permplus(nnn, s) = nplus(nnn, s)
          permmin(nnn, s) = nmin(nnn, s)
        end do
      end do

      !! debug text to understand the tiers and permutations
      !do s = 1, size(perm,1)
      !  do nnn = 1, size(perm,2)
      !    write(*,'(I2A)', advance='no') perm(s,nnn),' '
      !  end do
      !
      !  write(*,'(A)',advance='no') '   '
      !
      !  do s2 = 1, size(perm,2)
      !  if(permplus(s,s2) < 1 .or. permplus(s,s2) > size(perm,1)) then
      !    cycle
      !  end if
      !
      !  do nnn = 1, size(perm,2)
      !    write(*,'(I2A)', advance='no') perm(permplus(s,s2),nnn),' '
      !  end do
      !  end do !s2
      !
      !  write(*,'(A)',advance='no') '   '
      !
      !  do s2 = 1, size(perm,2)
      !  if(permmin(s,s2) < 1 .or. permmin(s,s2) > size(perm,1)) then
      !    cycle
      !  end if
      !
      !  do nnn = 1, size(perm,2)
      !    write(*,'(I2A)', advance='no') perm(permmin(s,s2),nnn),' '
      !  end do
      !  end do !s2
      !
      !  write(*,*)
      !end do
      !write(*,*)

      call print_log_message("index complete",5)
    end subroutine arend_init2

    subroutine arend_2D()
      ! REPHASING

      if (calculatereph) then

      signalpar(:,:) = 0.0
      signalper(:,:) = 0.0

      ! loop over polarizations
      do pol = 1, 1 ! 21

      dir1 = dirinda(pol)
      dir2 = dirindb(pol)
      dir3 = dirindc(pol)
      dir4 = dirindd(pol)
      orcoeffpar = orcoeffZZZZ(pol)


      ! Initial condition
      do nnn = 1, Nhier
        do s=1, Nsys
          rho1(nnn, s) = mu(s, dir1)
        end do
      end do


      do nnt1 = 1, Ntimestept1
        write(*,'(A6,I2,A12,I5,A3,I5)') "pol = ", pol, " / 21; t1 = ", nnt1, " / " , Ntimestept1
        call flush()

      ! GB

      ! initialize t3
        do nnn = 1, Nhier
          mem = 0.0
          do s = 1, Nsys
            mem = mem + mu(s, dir2) * rho1(nnn, s)
          end do
          do s = 1, Nsys
            rho3s(nnn, s) = mu(s, dir3) * mem
          end do
        end do

       do nnt3 = 1, Ntimestept3

           do s = 1, Nsys
              signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
           end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3s(dt)
           end do


       end do



      ! SE and IA

        ! initialize t2
        do nnn = 1, Nhier
          do s = 1, Nsys
            do s2 = 1, Nsys
              rho2(nnn, s, s2) = mu(s, dir2) * rho1(nnn, s2)
            end do
          end do
        end do

       ! propagate t2
       do nnt2 = 1, Ntimestept2
         do nin = 1, Ntimestept2in
           call arend_propagate2(dt)
         end do
       end do

      ! SE

       ! initialize t3
        rho3s = 0.0
        do nnn = 1, Nhier
          do s = 1, Nsys
            do s2 = 1, Nsys
              rho3s(nnn, s) = rho3s(nnn, s) + mu(s2, dir3) * rho2(nnn, s, s2)
            end do
          end do
        end do

       do nnt3 = 1, Ntimestept3

           ! calculate signal
           do s = 1, Nsys
              signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
            end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3s(dt)
           end do


       end do




      ! IA
        ! initialize t3
        rho3 = 0.0
        do nnn = 1, Nhier
          do s = 1, Nsys
           do w1 = 1, Nsys
             do w2 = (w1+1), Nsys
                w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
                rho3(nnn, w, s) = rho3(nnn, w, s) + mu(w2, dir3) * rho2(nnn, w1, s)
                rho3(nnn, w, s) = rho3(nnn, w, s) + mu(w1, dir3) * rho2(nnn, w2, s)
              end do
            end do
          end do
        end do


       do nnt3 = 1, Ntimestept3

           ! calculate signal
           do w1 = 1, Nsys
             do w2 = (w1+1), Nsys
                w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
                signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) - orcoeffpar * mu(w2, dir4) * rho3(1, w, w1)
                signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) - orcoeffpar * mu(w1, dir4) * rho3(1, w, w2)
              end do
            end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3(dt)
           end do


       end do

       ! propagate t1
       do nin = 1, Ntimestept1in
         call arend_propagate1reph(dt)
       end do


      end do

      end do ! pol

      call print_log_message("rephasing calculation complete",5)

!      !output
!      open (55, File='reph.par.re' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), real(signalpar(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='reph.per.re' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), real(signalper(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='reph.par.im' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), imag(signalpar(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='reph.per.im' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), imag(signalper(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)



      end if ! calculatereph



      ! NON-REPHASING

      if (calculatenonreph) then

      signalpar(:,:) = 0.0
      signalper(:,:) = 0.0

      ! loop over polarizations
      do pol = 1,1
      dir1 = dirinda(pol)
      dir2 = dirindb(pol)
      dir3 = dirindc(pol)
      dir4 = dirindd(pol)
      orcoeffpar = orcoeffZZZZ(pol)



      ! Initial condition
      do nnn = 1, Nhier
        do s=1, Nsys
          rho1(nnn, s) = mu(s, dir1)
        end do
      end do


      do nnt1 = 1, Ntimestept1
        write(*,'(A6,I2,A12,I5,A3,I5)') "pol = ", pol, " / 21; t1 = ", nnt1, " / " , Ntimestept1

      ! GB

      ! initialize t3
        do nnn = 1, Nhier
          mem = 0.0
          do s = 1, Nsys
            mem = mem + mu(s, dir2) * rho1(nnn, s)
          end do
          do s = 1, Nsys
            rho3s(nnn, s) = mu(s, dir3) * mem
          end do
        end do

       do nnt3 = 1, Ntimestept3

           do s = 1, Nsys
              signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
           end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3s(dt)
           end do


       end do



      ! SE and IA

       ! initialize t2
        do nnn = 1, Nhier
          do s = 1, Nsys
            do s2 = 1, Nsys
              rho2(nnn, s, s2) = mu(s2, dir2) * rho1(nnn, s)
            end do
          end do
        end do

       ! propagate t2
       do nnt2 = 1, Ntimestept2
         do nin = 1, Ntimestept2in
           call arend_propagate2(dt)
         end do
       end do

      ! SE

       ! initialize t3
        rho3s = 0.0
        do nnn = 1, Nhier
          do s = 1, Nsys
            do s2 = 1, Nsys
              rho3s(nnn, s) = rho3s(nnn, s) + mu(s2, dir3) * rho2(nnn, s, s2)
            end do
          end do
        end do

       do nnt3 = 1, Ntimestept3

           ! calculate signal
           do s = 1, Nsys
              signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) + orcoeffpar * mu(s, dir4) * rho3s(1, s)
           end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3s(dt)
           end do


       end do




      ! IA
        ! initialize t3
        rho3 = 0.0
        do nnn = 1, Nhier
          do s = 1, Nsys
           do w1 = 1, Nsys
             do w2 = (w1+1), Nsys
                w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
                rho3(nnn, w, s) = rho3(nnn, w, s) + mu(w2, dir3) * rho2(nnn, w1, s)
                rho3(nnn, w, s) = rho3(nnn, w, s) + mu(w1, dir3) * rho2(nnn, w2, s)
              end do
            end do
          end do
        end do


       do nnt3 = 1, Ntimestept3

           ! calculate signal
           do w1 = 1, Nsys
             do w2 = (w1+1), Nsys
                w = (w1-1)*Nsys - (w1*(w1+1))/2 + w2
                signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) - orcoeffpar * mu(w2, dir4) * rho3(1, w, w1)
                signalpar(nnt1, nnt3) = signalpar(nnt1, nnt3) - orcoeffpar * mu(w1, dir4) * rho3(1, w, w2)
              end do
            end do

           ! propagate t3
           do nin = 1, Ntimestept3in
             call arend_propagate3(dt)
           end do


       end do

       ! propagate t1
       do nin = 1, Ntimestept1in
         call arend_propagate1nonreph(dt)
       end do


      end do

      end do ! pol

      call print_log_message("non-rephasing calculation complete",5)

!      !output
!      open (55, File='nonreph.par.re' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), real(signalpar(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='nonreph.per.re' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), real(signalper(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='nonreph.par.im' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), imag(signalpar(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)

!      open (55, File='nonreph.per.im' )
!      do nnt1 = 1, Ntimestept1
!        do nnn = 1, Ntimestept3
!        write(55,'(F12.6)', advance='no'), imag(signalper(nnt1, nnn))
!        end do
!        write(55,*)
!      end do
!      close(55)



      end if ! calculatenonreph
    end subroutine arend_2D

    function numpermt(tier) ! number of permutations in tier
      integer(i4b), intent(in):: tier
      integer(i4b) :: numpermt, tmp
          numpermt = int(factorial(Nind+tier-1)/(factorial(Nind-1)*factorial(tier)))
    end function numpermt

    ! return the index of given permutation
    function permindex (permutation)
       integer(i4b):: permindex
       integer(i4b), intent(in):: permutation(Nind)
       logical:: isfound
       integer(i4b):: findex, tt, fend

            isfound = .False.
            findex = 1
            do tt = 1, sum(permutation)-1 ! find the tier where this permutation should be
                    findex = findex + numpermt (tt)
            end do
            fend = findex + numpermt(sum(permutation))+1
            if (findex >= Nhier) then ! the permutation does not exist
                    isfound = .True.
                    findex = 0
            end if
            do while ((isfound .eqv. .false.) .and. (findex < fend))  ! actually, findex >= fend should never happen if all permutations are present in perm
                                                     ! but we include it for use when perm can be pruned.
              if (ALL(perm(findex,:) == permutation(:))) then
                 isfound = .True.
              end if
              findex = findex + 1
            end do

            if (isfound .eqv. .false.) then
                    findex = 0
            end if

            permindex = findex - 1

    end function permindex

    function nplus (nin, jin)
      integer(i4b), intent(in):: nin, jin
      integer(i4b):: nplus

      currentperm = perm(nin, :)
      currentperm(jin) = currentperm(jin) + 1

      nplus = permindex(currentperm)

      if (nplus == -1) then
        nplus = Nhier + 1 ! the Nhier + 1 'th density matrix contains only zeros
      end if

    end function nplus


    function nmin (nin, jin)
      integer(i4b), intent(in):: nin, jin
      integer(i4b):: nmin

      currentperm = perm(nin, :)

      nmin = -1
      if (currentperm(jin) > 0) then
        currentperm(jin) = currentperm(jin) - 1
        nmin = permindex(currentperm)
      end if

      if (nmin == -1) then
        nmin = Nhier + 1 ! the Nhier + 1 'th density matrix contains only zeros
      end if

    end function nmin


    function nu(j, m)
      integer(i4b), intent(in):: j, m
      complex(dpc):: nu
      if (m == 0) then
        nu = LLambda(j)
      else
        nu = 2*pi*m/beta(j)
      end if
    end function nu

    function cconst(j)
      integer(i4b), intent(in):: j
      complex(dpc):: cconst

      cconst = 2*lambda(j)/beta(j) - dcmplx(0.0,1.0) * lambda(j) * LLambda(j)

    end function cconst



    subroutine arend_initmult1
      ! initialize the propagator for a coherence |1><0|
      integer(i4b):: n, j
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)
      complex(dpc):: musum, nu1, jsum
      complex(dpc), allocatable:: identity(:,:)

      ALLOCATE(identity,(Nsys, Nsys))

      identity = 0.0_dp
      do n=1, Nsys
        identity(n, n) = 1.0_dp
      end do

      opLeft1 = 0.0
      opPlusLeft1 = 0.0
      opMinLeft1 = 0.0

      do n = 1, Nhier

        musum = 0
        do j=1, Nsys
          musum = musum + perm(n, j) * LLambda(j)
        end do


        opLeft1(n,:,:) = opLeft1(n,:,:) - iconst * HS

        opLeft1(n,:,:) = opLeft1(n,:,:) - musum * identity

        do j=1, Nsys
          ! first low temperature correction term, see Ishizaki PNAS 2009
          nu1 = 2*pi/beta(j)
          jsum = 2*(lambda(j) / beta(j)) * 2*LLambda(j)/(nu1*nu1 - LLambda(j)*LLambda(j))
          opLeft1(n,:,:) = opLeft1(n,:,:) - jsum * MATMUL(V(j,:,:), V(j,:,:))

          opPlusLeft1(n,j,:,:) =  opPlusLeft1(n,j,:,:) - iconst * V(j,:,:)

          opMinLeft1(n,j,:,:) =  opMinLeft1(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V(j,:,:)

          ! first low temperature correction term, see Ishizaki PNAS 2009
          opMinLeft1(n,j,:,:) = opMinLeft1(n,j,:,:) - iconst * perm(n,j)*jsum * LLambda(j) * V(j,:,:)

        end do

      end do

      DEALLOCATE(identity)
    end subroutine arend_initmult1





    subroutine arend_initmult2
      ! initialize the propagator for populations and coherences |1><1'|
      integer(i4b):: n, j
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)
      complex(dpc):: musum, nu1, jsum
      complex(dpc), allocatable:: identity(:,:)

      ALLOCATE(identity,(Nsys, Nsys))

      identity = 0.0_dp
      do n=1, Nsys
        identity(n, n) = 1.0_dp
      end do

      opLeft2 = 0.0
      opRight2 = 0.0
      opLRLeft2 = 0.0
      opLRRight2 = 0.0
      opPlusLeft2 = 0.0
      opPlusRight2 = 0.0
      opMinLeft2 = 0.0
      opMinRight2 = 0.0

      do n = 1, Nhier

        musum = 0
        do j=1, Nsys
          musum = musum + perm(n, j) * LLambda(j)
        end do


        opLeft2(n,:,:) = opLeft2(n,:,:) - iconst * HS
        opRight2(n,:,:) = opRight2(n,:,:) + iconst * HS

        opLeft2(n,:,:) = opLeft2(n,:,:) - musum * identity

        do j=1, Nsys
          ! first low temperature correction term, see Ishizaki PNAS 2009
          nu1 = 2*pi/beta(j)
          jsum = 2*(lambda(j) / beta(j)) * 2*LLambda(j)/(nu1*nu1 - LLambda(j)*LLambda(j))
          opLeft2(n,:,:) = opLeft2(n,:,:) - jsum * MATMUL(V(j,:,:), V(j,:,:))
          opRight2(n,:,:)= opRight2(n,:,:) -  jsum * MATMUL(V(j,:,:), V(j,:,:))
          opLRLeft2(n,j,:,:) =  opLRLeft2(n,j,:,:) + 2*jsum * V(j,:,:)
          opLRRight2(n,j,:,:) =  opLRRight2(n,j,:,:) + V(j,:,:)

          opPlusLeft2(n,j,:,:) =  opPlusLeft2(n,j,:,:) - iconst * V(j,:,:)
          opPlusRight2(n,j,:,:) =  opPlusRight2(n,j,:,:) + iconst * V(j,:,:)

          opMinLeft2(n,j,:,:) =  opMinLeft2(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V(j,:,:)
          opMinRight2(n,j,:,:) = opMinRight2(n,j,:,:) + iconst*perm(n, j)*conjg(cconst(j)) * V(j,:,:)

          ! first low temperature correction term, see Ishizaki PNAS 2009
          opMinLeft2(n,j,:,:) = opMinLeft2(n,j,:,:) - iconst * perm(n,j)*jsum * LLambda(j) * V(j,:,:)
          opMinRight2(n,j,:,:) = opMinRight2(n,j,:,:) + iconst * perm(n,j)*jsum * LLambda(j) * V(j,:,:)

        end do
      end do

      DEALLOCATE(identity)
    end subroutine arend_initmult2


    subroutine arend_initmult3
      ! initialize the propagator for 2-1 exciton coherences |2><1|
      integer:: n, j
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)
      complex(dpc):: musum, nu1, jsum
      complex(dpc), allocatable:: identity(:,:)

      ALLOCATE(identity,((Nsys*(Nsys-1))/2, (Nsys*(Nsys-1))/2))

      identity = 0.0_dp
      do n=1, (Nsys*(Nsys-1))/2
        identity(n, n) = 1.0_dp
      end do

      opLeft3 = 0.0
      opRight3 = 0.0
      opLRLeft3 = 0.0
      opLRRight3 = 0.0
      opPlusLeft3 = 0.0
      opPlusRight3 = 0.0
      opMinLeft3 = 0.0
      opMinRight3 = 0.0

      do n = 1, Nhier

        musum = 0
        do j=1, Nsys
          musum = musum + perm(n, j) * LLambda(j)
        end do


        opLeft3(n,:,:) = opLeft3(n,:,:) - iconst * HS2
        opRight3(n,:,:) = opRight3(n,:,:) + iconst * HS

        opLeft3(n,:,:) = opLeft3(n,:,:) - musum * identity

        do j=1, Nsys
          ! first low temperature correction term, see Ishizaki PNAS 2009
          nu1 = 2*pi/beta(j)
          jsum = 2*(lambda(j) / beta(j)) * 2*LLambda(j)/(nu1*nu1 - LLambda(j)*LLambda(j))
          opLeft3(n,:,:) = opLeft3(n,:,:) - jsum * MATMUL(V2(j,:,:), V2(j,:,:))
          opRight3(n,:,:)= opRight3(n,:,:) -  jsum * MATMUL(V(j,:,:), V(j,:,:))
          opLRLeft3(n,j,:,:) =  opLRLeft3(n,j,:,:) + 2*jsum * V2(j,:,:)
          opLRRight3(n,j,:,:) =  opLRRight3(n,j,:,:) + V(j,:,:)
          opPlusLeft3(n,j,:,:) =  opPlusLeft3(n,j,:,:) - iconst * V2(j,:,:)
          opPlusRight3(n,j,:,:) =  opPlusRight3(n,j,:,:) + iconst * V(j,:,:)
          opMinLeft3(n,j,:,:) =  opMinLeft3(n,j,:,:) - iconst*perm(n, j)*cconst(j) * V2(j,:,:)
          opMinRight3(n,j,:,:) = opMinRight3(n,j,:,:) + iconst*perm(n, j)*conjg(cconst(j)) * V(j,:,:)
          ! first low temperature correction term, see Ishizaki PNAS 2009
          opMinLeft3(n,j,:,:) = opMinLeft3(n,j,:,:) - iconst * perm(n,j)*jsum * LLambda(j) * V2(j,:,:)
          opMinRight3(n,j,:,:) = opMinRight3(n,j,:,:) + iconst * perm(n,j)*jsum * LLambda(j) * V(j,:,:)
        end do

      end do

      DEALLOCATE(identity)
    end subroutine arend_initmult3





    subroutine arend_Lmult1 (tt, rhoin, result)
      complex(dpc), intent(in)  :: rhoin(:,:)
      real(dp), intent(in)      :: tt ! this is completely unneccessary parameter to satisfy ode_rk4 function template
      complex(dpc), intent(out) :: result(:,:)
      integer:: n,j, nnp
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)

      result(:,:) = 0.0

      do n = 1, Nhier

        result(n,:) = result(n,:) + MATMUL(opLeft1(n,:,:), rhoin(n,:))

        do j=1, Nsys
          result(n,:) = result(n,:) + MATMUL(opPlusLeft1(n, j, :, :), rhoin(permplus(n,j),:))
          result(n,:) = result(n,:) + MATMUL(opMinLeft1(n,j,:,:), rhoin(permmin(n,j),:))
        end do

      end do
    end subroutine arend_Lmult1

    subroutine arend_propagate1reph (dt)
      real(dp), intent(in):: dt

      real(dp) :: t
      t = dt ! this is completely unneccessary parameter to satisfy ode_rk4 function template


      !!  second order
      !rho1 = conjg(rho1) ! Lmult1 propates the ket, want the bra for rephasing response.
      !call arend_Lmult1(t, rho1, prhox1)
      !rho1 = rho1 + dt * 0.25 * prhox1
      !prhodx1 = rho1 + 0.66666666666666666667 * dt * prhox1
      !call arend_Lmult1(t, prhodx1, prhox1) ! prhox1 -> prhodt1 for speed instead of memory
      !!rho1 = rho1 + dt * (0.25*prhox1 + 0.75*prhodt1) for speed instead of memory
      !rho1 = rho1 + dt * 0.75*prhox1
      !rho1 = conjg(rho1)

      rho1 = conjg(rho1) ! Lmult1 propates the ket, want the bra for rephasing response.
      call arend_Lmult1(t,rho1,prhodx1)
      call ode_rk4(rho1,prhodx1,t,dt,prhox1,arend_Lmult1)
      rho1 = prhox1
      rho1 = conjg(rho1)


    ! fourth order
    !  call arend_Lmult(rho1, prhodt1)
    !  prhox1 = rho1 + (dt/2.0) * prhodt1
    !  call arend_Lmult1(prhox1, prhodx1)
    !  prhox1 = rho1 + (dt/2.0) * prhodx1
    !  call arend_Lmult1(prhox1, prhodm1)
    !  prhox1 = rho1 + dt * prhodm1
    !  prhodm1 = prhodm1 + prhodx1
    !  call arend_Lmult1(prhox1, prhodx1)
    !  rho1 = rho1 + (dt/6.0)*(prhodt1 + prhodx1 + 2.0*prhodm1)
    end subroutine arend_propagate1reph


    subroutine arend_propagate1nonreph (dt)
      real(dp), intent(in):: dt

      real(dp) :: t
      t = dt ! this is completely unneccessary parameter to satisfy ode_rk4 function template

      !!  second order
      !call arend_Lmult1(t, rho1, prhox1)
      !rho1 = rho1 + dt * 0.25 * prhox1
      !prhodx1 = rho1 + 0.66666666666666666667 * dt * prhox1
      !call arend_Lmult1(t, prhodx1, prhox1)
      !rho1 = rho1 + dt * 0.75*prhox1

      call arend_Lmult1(t,rho1,prhodx1)
      call ode_rk4(rho1,prhodx1,t,dt,prhox1,arend_Lmult1)
      rho1 = prhox1
    end subroutine arend_propagate1nonreph


    subroutine arend_Lmult2 (tt, rhoin, result)
      complex(dpc), intent(in)  :: rhoin(:,:,:)
      real(dp), intent(in)      :: tt ! this is a dummy parameter to satisfy ode_rk4 function template
      complex(dpc), intent(out) :: result(:,:,:)
      integer:: n,j, nnp
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)

      result(:,:,:) = 0.0

      do n = 1, Nhier
        result(n,:,:) = result(n,:,:) + MATMUL(opLeft2(n,:,:), rhoin(n,:,:))
        result(n,:,:) = result(n,:,:) + MATMUL(rhoin(n,:,:), opRight2(n,:,:))

        do j=1, Nsys
          result(n,:,:) = result(n,:,:) + MATMUL(MATMUL(opLRLeft2(n,j,:,:), rhoin(n,:,:)), opLRRight2(n,j,:,:))
          result(n,:,:) = result(n,:,:) + MATMUL(opPlusLeft2(n,j, :,:), rhoin(permplus(n,j),:,:))
          result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permplus(n,j),:,:), opPlusRight2(n,j,:,:))
          result(n,:,:) = result(n,:,:) + MATMUL(opMinLeft2(n,j, :,:), rhoin(permmin(n,j),:,:))
          result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permmin(n,j),:,:), opMinRight2(n,j,:,:))
        end do

      if(light_hierarchy) then
        ! relaxation
        result(n,:(Nsys-1),:(Nsys-1)) = result(n,:(Nsys-1),:(Nsys-1)) - global_relax_rate * rhoin(n,:(Nsys-1),:(Nsys-1))

        do s=1, Nsys-1
          result(n,Nsys,Nsys) = result(n,Nsys,Nsys) + global_relax_rate * rhoin(n,s,s)
        end do
      end if

      end do
    end subroutine arend_Lmult2

    subroutine arend_propagate2 (dt)
      real(dp), intent(in) :: dt

      real(dp) :: t
      t = dt ! this is a dummy parameter to satisfy ode_rk4 function template

    !  !  second order
    !  call arend_Lmult2(t, rho2, prhox2)
    !  rho2 = rho2 + dt * 0.25 * prhox2
    !  prhodx2 = rho2 + 0.66666666666666666667 * dt * prhox2
    !  call arend_Lmult2(t, prhodx2, prhox2) ! prhox -> prhodt for speed instead of memory
    !  rho2 = rho2 + dt * 0.75*prhox2

      call arend_Lmult2(t,rho2,prhodx2)
      call ode_rk4(rho2,prhodx2,t,dt,prhox2,arend_Lmult2)
      rho2 = prhox2
    end subroutine arend_propagate2

    subroutine arend_LmultC (ttt, rhoin, result)
      complex(dpc), intent(in)  :: rhoin(:,0:,0:)
      real(dp), intent(in)      :: ttt ! time (fs)
      complex(dpc), intent(out) :: result(:,0:,0:)
      integer:: n,j, nnp
      complex(dpc), parameter   :: iconst = dcmplx(0.0, 1.0)
      real(dp)                  :: EE, tt

      result(:,:,:) = 0.0
      EE = bloch_strength
      tt = ttt

      if(bloch_term) then
        do n = 1, Nhier
          do s=1, Nsys
            result(n,0,0) = result(n,0,0) - iconst*cos(central_frequency*tt)*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s,0)*EE               &
                                          - conjg( iconst*cos(central_frequency*tt)*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s,0)*EE )

            result(n,s,0) = result(n,s,0) - iconst*cos(central_frequency*tt)*exp(iconst*rwa*tt)*mu(s, 1)*rhoin(n,0,0)*EE

          do s2=1, Nsys
            result(n,s,0) = result(n,s,0) + iconst*cos(central_frequency*tt)*exp(iconst*rwa*tt)*mu(s2, 1)*rhoin(n,s,s2)*EE

            result(n,s,s2) = result(n,s,s2) + iconst*cos(central_frequency*tt)*exp(-iconst*rwa*tt)*mu(s2, 1)*rhoin(n,s,0)*EE            &
                                            + conjg( iconst*cos(central_frequency*tt)*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s2,0)*EE )
          end do
          end do
        end do
      end if

      if(noise_term) then
        tt = min(max(0.5_dp, tt), dt*Ntimestept1*Ntimestept1In - 0.5)

        do n = 1, Nhier
          do s=1, Nsys
            result(n,0,0) = result(n,0,0) - iconst*oonoise(int(tt/dt+1))*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s,0)*EE               &
                                          - conjg( iconst*oonoise(int(tt/dt+1))*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s,0)*EE )

            result(n,s,0) = result(n,s,0) - iconst*oonoise(int(tt/dt+1))*exp(iconst*rwa*tt)*mu(s, 1)*rhoin(n,0,0)*EE

          do s2=1, Nsys
            result(n,s,0) = result(n,s,0) + iconst*oonoise(int(tt/dt+1))*exp(iconst*rwa*tt)*mu(s2, 1)*rhoin(n,s,s2)*EE

            result(n,s,s2) = result(n,s,s2) + iconst*oonoise(int(tt/dt+1))*exp(-iconst*rwa*tt)*mu(s2, 1)*rhoin(n,s,0)*EE            &
                                            + conjg( iconst*oonoise(int(tt/dt+1))*exp(-iconst*rwa*tt)*mu(s, 1)*rhoin(n,s2,0)*EE )
          end do
          end do
        end do
      end if

      ! Lmult1 part
      do n = 1, Nhier
        result(n,1:,0) = result(n,1:,0) + MATMUL(opLeft1(n,:,:), rhoin(n,1:,0))

        do j=1, Nsys
          result(n,1:,0) = result(n,1:,0) + MATMUL(opPlusLeft1(n, j, :, :), rhoin(permplus(n,j),1:,0))
          result(n,1:,0) = result(n,1:,0) + MATMUL(opMinLeft1(n,j,:,:), rhoin(permmin(n,j),1:,0))
        end do

!!        ! conjugate part of Lmult1
!!        do j=1, Nsys
!!          result(n,0,j) = conjg(result(n,j,0))
!!        end do
      end do

      ! Lmult2 part
      do n = 1, Nhier
        result(n,1:,1:) = result(n,1:,1:) + MATMUL(opLeft2(n,:,:), rhoin(n,1:,1:))
        result(n,1:,1:) = result(n,1:,1:) + MATMUL(rhoin(n,1:,1:), opRight2(n,:,:))

        do j=1, Nsys
          result(n,1:,1:) = result(n,1:,1:) + MATMUL(MATMUL(opLRLeft2(n,j,:,:), rhoin(n,1:,1:)), opLRRight2(n,j,:,:))
          result(n,1:,1:) = result(n,1:,1:) + MATMUL(opPlusLeft2(n,j, :,:), rhoin(permplus(n,j),1:,1:))
          result(n,1:,1:) = result(n,1:,1:) + MATMUL(rhoin(permplus(n,j),1:,1:), opPlusRight2(n,j,:,:))
          result(n,1:,1:) = result(n,1:,1:) + MATMUL(opMinLeft2(n,j, :,:), rhoin(permmin(n,j),1:,1:))
          result(n,1:,1:) = result(n,1:,1:) + MATMUL(rhoin(permmin(n,j),1:,1:), opMinRight2(n,j,:,:))
        end do

        ! relaxation
        result(n,1:,1:) = result(n,1:,1:) - global_relax_rate * rhoin(n,1:,1:)

        do s=1, Nsys
          result(n,0,0) = result(n,0,0) + global_relax_rate * rhoin(n,s,s)
        end do
      end do

      ! channel to the ground state here?
    end subroutine arend_LmultC

    subroutine arend_propagateC (dt, time)
      real(dp), intent(in) :: dt, time

      real(dp) :: t
      t = time

      call arend_LmultC(t,rhoC,prhodxC)
      call ode_rk4(rhoC,prhodxC,t,dt,prhoxC,arend_LmultC)
      rhoC = prhoxC
    end subroutine arend_propagateC

    subroutine arend_Lmult3 (tt, rhoin, result)
    complex(dpc), intent(in)  :: rhoin(:,:,:)
    real(dp), intent(in)      :: tt ! this is completely unneccessary parameter to satisfy ode_rk4 function template
    complex(dpc), intent(out) :: result(:,:,:)
    integer:: n,j, nnp
    complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)


    result(:,:,:) = 0.0


    do n = 1, Nhier

      result(n,:,:) = result(n,:,:) + MATMUL(opLeft3(n,:,:), rhoin(n,:,:))
      result(n,:,:) = result(n,:,:) + MATMUL(rhoin(n,:,:), opRight3(n,:,:))

      do j=1, Nsys
        result(n,:,:) = result(n,:,:) + MATMUL(MATMUL(opLRLeft3(n,j,:,:), rhoin(n,:,:)), opLRRight3(n,j,:,:))
        result(n,:,:) = result(n,:,:) + MATMUL(opPlusLeft3(n,j, :,:), rhoin(permplus(n,j),:,:))
        result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permplus(n,j),:,:), opPlusRight3(n,j,:,:))
        result(n,:,:) = result(n,:,:) + MATMUL(opMinLeft3(n,j, :,:), rhoin(permmin(n,j),:,:))
        result(n,:,:) = result(n,:,:) + MATMUL(rhoin(permmin(n,j),:,:), opMinRight3(n,j,:,:))
      end do

    end do

    end subroutine arend_Lmult3


    subroutine arend_propagate3 (dt)
      real(dp), intent(in):: dt

      real(dp) :: t
      t = dt ! this is completely unneccessary parameter to satisfy ode_rk4 function template

      !!  second order
      !call arend_Lmult3(t, rho3, prhox3)
      !rho3 = rho3 + dt * 0.25 * prhox3
      !prhodx3 = rho3 + 0.66666666666666666667 * dt * prhox3
      !call arend_Lmult3(t, prhodx3, prhox3) ! prhox -> prhodt for speed instead of memory
      !rho3 = rho3 + dt * 0.75*prhox3

      call arend_Lmult3(t,rho3,prhodx3)
      call ode_rk4(rho3,prhodx3,t,dt,prhox3,arend_Lmult3)
      rho3 = prhox3
    end subroutine arend_propagate3


    subroutine arend_propagate3s (dt)
      real(dp), intent(in):: dt

      real(dp) :: t
      t = dt ! this is completely unneccessary parameter to satisfy ode_rk4 function template

      !!  second order
      !call arend_Lmult1(t, rho3s, prhox3s)
      !rho3s = rho3s + dt * 0.25 * prhox3s
      !prhodx3s = rho3s + 0.66666666666666666667 * dt * prhox3s
      !call arend_Lmult1(t, prhodx3s, prhox3s)
      !rho3s = rho3s + dt * 0.75*prhox3s

      call arend_Lmult1(t,rho3s,prhodx3s)
      call ode_rk4(rho3s,prhodx3s,t,dt,prhox3s,arend_Lmult1)
      rho3s = prhox3s
    end subroutine arend_propagate3s



    subroutine init_generate_noise
      integer(i4b) :: b, b2, nb, clock
      real(dp), allocatable :: tmp(:,:)

      !real(dp), dimension(3,3) :: AAA, UUU
      !complex(dpc), dimension(2,2) :: BBB, VVV
      !real(dp), dimension(1)   :: ran
      !real(dp)                 :: sum, sum2

      !AAA(1,1) =   4
      !AAA(1,2) =  12
      !AAA(1,3) = -16
      !AAA(2,1) =  12
      !AAA(2,2) =  37
      !AAA(2,3) = -43
      !AAA(3,1) = -16
      !AAA(3,2) = -43
      !AAA(3,3) =  98

      !call cholesky(AAA,UUU)

      !write(*,*)
      !write(*,*) AAA
      !write(*,*)
      !write(*,*) UUU
      !write(*,*)

      !BBB(1,1) = 1
      !BBB(1,2) = cmplx(0,1)/2
      !BBB(2,1) = -cmplx(0,1)/2
      !BBB(2,2) = 1

      !call cholesky(BBB,VVV)

      !write(*,*)
      !write(*,*) BBB
      !write(*,*)
      !write(*,*) VVV
      !write(*,*)

      !sum = 0
      !sum2 = 0
      !do b=1,10000000
      !  call random_Normal(ran)
      !  !write(*,*) 'rrr', ran
      !  sum = sum + ran(1)
      !  sum2 = sum2 + ran(1)*ran(1)
      !end do
      !write (*,*)
      !write (*,*) sum/10000000, sum2/10000000
      !write (*,*)
      !stop

      CALL SYSTEM_CLOCK(COUNT=clock)
      call init_random(clock)

      nb = Ntimestept1*Ntimestept1in

      ALLOCATE(CholeskyCF,(nb, nb))
      ALLOCATE(tmp,(nb, nb))

      CholeskyCF = 0.0_dp
      tmp = 0.0_dp
      do b=1, nb
      do b2=1, nb
        if(b >= b2) then
          tmp(b,b2) = real(light_CF(dt*b,dt*b2))
        else
          tmp(b,b2) = real(light_CF(dt*b2,dt*b))
        end if
      end do
      end do

      call cholesky(tmp,CholeskyCF)

            !write(*,*) real(tmp(:5,:5))
            !write(*,*)
            !write(*,*) real(CholeskyCF(1,:))

      tmp = transpose(CholeskyCF)

      CholeskyCF = tmp

      DEALLOCATE(tmp)

    end subroutine init_generate_noise

    subroutine generate_noise(out)
      real(dpc), intent(out) :: out(:)

      real(dp)     :: oo(size(out))
      integer(i4b) :: b

      call random_Normal(oo)

      out = matmul(CholeskyCF,oo)
    end subroutine generate_noise















    !
    ! this part runs the phi, external program by Strumpfer and Schulten
    ! -------------------------------------------------------------------------------------
    !
    subroutine fill_evolution_superoperator_hierarchy_phi(type)
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
    end subroutine fill_evolution_superoperator_hierarchy_phi

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
        write(32,'(A)') 'LLambda:'
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

