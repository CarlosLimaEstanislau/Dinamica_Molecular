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
    use rng_mod, only: shuffle_array
    implicit none

    real(dp), parameter :: kappa = 2.79_dp      ! Screening lengh (micrometro^-1)
    real(dp), parameter :: bjerrum = 0.00072_dp ! Bjerrum lenght (micrometro)
    real(dp), parameter :: eps = 50.0_dp        ! WCA strength (KbT)
    real(dp), parameter :: dt = 0.001_dp        ! Integration step
    real(dp), parameter :: pi = acos(-1.0_dp)   ! pi constant
    real(dp), parameter :: a0 = 0.53_dp         ! Reference radii (micrometro)

    type :: system_g
        real(dp) :: box, total_energy, temp_target
        integer  :: g
        real(dp) :: pot_DLVO, pot_WCA, total_potential
        real(dp), dimension(:,:), allocatable :: total_forces, force_DLVO, force_WCA !(3,n)
        real(dp), allocatable :: r_save(:,:)
        integer, allocatable  :: list(:), point(:)
        real(dp) :: skin, r_cut_dlvo, r_list
        integer :: n, max_pairs
        integer :: num_particles
        real(dp):: rho, frac_particles,frac_charges
        real(dp):: kappa, bjerrum, max_radius
        real(dp):: gamma                                                         
    end type system_g

    type :: particles
        real(dp), dimension(:,:) , allocatable:: positions, velocities !(3,n)
        real(dp), dimension(:), allocatable :: masses, charges, radius !(n)
        real(dp) :: potential
        real(dp):: Z, radius1, radius2
    end type particles 

    type :: parameters
        character(len=256):: init_filename, part_filename, file_nml, outdir
        integer :: max_steps
        logical :: is_first
        logical :: interrupt
        logical :: exists 
    end type parameters 

    contains

        subroutine init_params(params)
            implicit none
            type(parameters), intent(out) :: params
            
            params%init_filename = 'init_conf.dat'
            params%part_filename = 'particles.part'
            params%file_nml = 'config.nml'
            params%is_first = .true.
            params%interrupt = .false.
            params%exists = .false.
        end subroutine

        subroutine alloc_particles(part, sys)
            implicit none
            type(particles), intent(inout):: part
            type(system_g), intent(in):: sys

            allocate(part%positions(3, sys%num_particles))
            allocate(part%velocities(3, sys%num_particles))
            allocate(part%masses(sys%num_particles))
            allocate(part%charges(sys%num_particles))
            allocate(part%radius(sys%num_particles))
            part%masses = 1.0_dp

        end subroutine

        subroutine set_particles(part, sys)
            implicit none
            type(particles), intent(inout):: part
            type(system_g), intent(in):: sys
            integer:: num_radius, num_charges

            part%positions = 0.0_dp
            part%velocities = 0.0_dp
            
            part%charges = 0.0_dp
            part%radius  = 0.0_dp
            
            part%masses = 1.0_dp

            if (sys%frac_particles .le. 1.0_dp .and. sys%frac_particles .ge. 0.0_dp) then 
                num_radius  = nint(sys%frac_particles * sys%num_particles)
            else 
                error stop "Erro: Atribua um valor para fração de particulas entre 0 e 1 !"
            end if
            
            if (sys%frac_charges .le. 1.0_dp .and. sys%frac_charges .ge. 0.0_dp) then 
                num_charges = nint(sys%frac_charges * sys%num_particles)
            else 
                error stop "Erro: Atribua um valor para fração de cargas entre 0 e 1!"
            end if

            part%radius = part%radius1/a0
            if (num_radius < sys%num_particles) then
                part%radius(num_radius+1:) = part%radius2/a0
            end if
            
            part%charges = part%Z
            if (num_charges < sys%num_particles) then
                part%charges(num_charges+1:) = -part%Z
            end if

            call shuffle_array(part%radius, part%charges)
            !call shuffle_array(part%radius)
            !call shuffle_array(part%charges)
        end subroutine set_particles

        subroutine dealloc_particles(part)
            implicit none
            type(particles), intent(inout) :: part

            if(allocated(part%positions))deallocate(part%positions)
            if(allocated(part%velocities))deallocate(part%velocities)
            if(allocated(part%masses)) deallocate(part%masses)
            if(allocated(part%charges)) deallocate(part%charges)
            if(allocated(part%radius)) deallocate(part%radius)
        end subroutine 

        subroutine init_verlet_list(sys, part)
            implicit none
            type(system_g), intent(inout) :: sys
            type(particles), intent(in)   :: part
            
            sys%max_radius = maxval(part%radius)
            sys%skin = 0.2_dp * sys%r_cut_dlvo
            sys%r_list = sys%r_cut_dlvo + sys%skin + 2.0_dp * sys%max_radius
            sys%max_pairs = sys%num_particles * (sys%num_particles - 1) / 2

            allocate(sys%list(sys%max_pairs))
            allocate(sys%point(sys%num_particles +1))
            allocate(sys%r_save(3,sys%num_particles ))

        end subroutine init_verlet_list

        subroutine init_system(sys)
            implicit none
            type(system_g), intent(inout) :: sys

            sys%g = 3*(sys%num_particles-1)
            sys%total_energy = 0.0_dp
            sys%temp_target = 1.0_dp
            sys%pot_WCA = 0.0_dp
            sys%pot_DLVO = 0.0_dp
            sys%total_potential = 0.0_dp
            sys%kappa = kappa*a0 
            sys%bjerrum = bjerrum/a0
            sys%r_cut_dlvo = 5.0_dp / sys%kappa

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