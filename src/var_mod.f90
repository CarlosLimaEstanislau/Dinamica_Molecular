!----------------------------------------------------------------------!
! Module resposable to handle the variables of the simulation.         !
! Here is defined types of variables for more readibilit, being        !  
! sys: system variable                                                 !
! particles: particles variables                                       !
! params: parametres of the simulation as well thermostat related      !
! variables                                                            !
! Other thing this module does is to deal with the alocation and       !
! initiaion of the variables.                                          !
!----------------------------------------------------------------------!
module var_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: kappa = 0.0279_dp ! Screening lengh
    real(dp), parameter :: eps = 1000.0_dp     ! WCA strength
    real(dp), parameter :: dt = 0.001_dp       ! Integration step
    real(dp), parameter :: pi = acos(-1.0_dp)  ! pi constant

    type :: system_g
        real(dp) :: box, total_energy, temp_target
        integer  :: g
        real(dp) :: pot_DLVO, pot_WCA, total_potential
        real(dp), dimension(:,:), allocatable :: total_forces, force_DLVO, force_WCA !(3,n)
        integer :: num_particles
        real(dp):: rho, frac_particles,frac_charges
    end type system_g

    type :: particles
        real(dp), dimension(:,:) , allocatable:: positions, velocities !(3,n)
        real(dp), dimension(:), allocatable :: masses, charges, radius !(n)
        real(dp) :: potential
        real(dp):: Z, radius1, radius2
    end type particles 

    type :: parameters
        character(len=256):: filename
        character(len=256):: file_nml
        character(len=256):: outdir
        integer :: max_steps
        logical :: is_first
        logical :: interrupt
        logical :: exists 
        real(dp):: gamma
    end type parameters 

    contains

        subroutine init_params(params)
            implicit none
            type(parameters) :: params
            
            params%filename = 'init_conf.dat'
            params%file_nml = 'config.nml'
            params%is_first = .true.
            params%interrupt = .false.
            params%exists = .false.
            params%gamma = 2.0_dp 
        end subroutine

        subroutine init_basic(part, sys)
            implicit none
            
            type(particles), intent(inout):: part
            type(system_g), intent(in):: sys
            
            allocate(part%positions(3, sys%num_particles))
            allocate(part%velocities(3, sys%num_particles))
            part%positions = 0.0_dp
            part%velocities = 0.0_dp
        end subroutine

        subroutine init_particles(part, sys)
            implicit none
            type(particles), intent(inout):: part
            type(system_g), intent(in):: sys
            integer:: num_radius, num_charges

            call init_basic(part,sys)
            allocate(part%masses(sys%num_particles))
            allocate(part%charges(sys%num_particles))
            allocate(part%radius(sys%num_particles))
            
            part%charges = 0.0_dp
            part%radius  = 0.0_dp
            
            part%masses = 1.0_dp

            if (sys%frac_particles .le. 1.0_dp .and. sys%frac_particles .ge. 0.0_dp) then 
                num_radius  = max(1, nint(sys%frac_particles * sys%num_particles))
            else 
                error stop "Erro: Atribua um valor para fração de particulas entre 1 e 0!"
            end if
            
            if (sys%frac_charges .le. 1.0_dp .and. sys%frac_charges .ge. 0.0_dp) then 
                num_charges = max(1, nint(sys%frac_charges   * sys%num_particles))
            else 
                error stop "Erro: Atribua um valor para fração de cargas entre 1 e 0!"
            end if

            if(num_radius .eq. sys%num_particles) then
                part%radius(1:num_radius) = part%radius1
            else
                part%radius(1:num_radius) = part%radius1
                part%radius(num_radius+1:sys%num_particles) = part%radius2 
            end if
            
            if(num_charges .eq. sys%num_particles) then
                part%charges(1:num_charges) = part%Z
            else
                part%charges(1:num_charges) = part%Z
                part%charges(num_charges+1:sys%num_particles) = -part%Z
            end if   

        end subroutine

        subroutine end_basic(part)
            implicit none
            type(particles), intent(inout) :: part

            if(allocated(part%positions)) deallocate(part%positions)
            if(allocated(part%velocities)) deallocate(part%velocities)
        end subroutine

        subroutine end_particles(part)
            implicit none
            type(particles), intent(inout) :: part

            call end_basic(part)
            if(allocated(part%masses)) deallocate(part%masses)
            if(allocated(part%charges)) deallocate(part%charges)
            if(allocated(part%radius)) deallocate(part%radius)
        end subroutine 

        subroutine init_system(sys)
            implicit none
            type(system_g), intent(inout) :: sys

            sys%g = 3*sys%num_particles
            sys%total_energy = 0.0_dp
            sys%temp_target = 1.0_dp
            sys%pot_WCA = 0.0_dp
            sys%pot_DLVO = 0.0_dp
            sys%total_potential = 0.0_dp

            allocate(sys%total_forces(3, sys%num_particles))
            allocate(sys%force_DLVO(3, sys%num_particles))
            allocate(sys%force_WCA(3, sys%num_particles))
            sys%total_forces = 0.0_dp
            sys%force_DLVO = 0.0_dp
            sys%force_WCA = 0.0_dp
        end subroutine

        subroutine end_system(sys)
            implicit none
            type(system_g), intent(inout)::sys

            if(allocated(sys%total_forces)) deallocate(sys%total_forces)
            if(allocated(sys%force_DLVO)) deallocate(sys%force_DLVO)
            if(allocated(sys%force_WCA)) deallocate(sys%force_WCA)
        end subroutine
end module