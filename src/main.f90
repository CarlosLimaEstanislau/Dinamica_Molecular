program simulation
    use, intrinsic:: iso_fortran_env, only: dp=>real64, output_unit
    use forces_mod
    use handler_mod
    use io_mod
    use ui_mod
    use var_mod
    implicit none

    !------------types-----------------
    type(parameters) :: params
    type(particles)  :: part
    type(system_g)   :: sys

    !------------variables-------------
    integer  :: step
    real(dp) :: kinetic, temp

    call init_params(params)
    call get_command_argument(1, params%outdir)

    call read_nml(sys, part, params)

    call init_system(sys)
    call init_particles(part, sys)

    call read_init_config(sys, part, params)
    call initialize(sys, params, part, dt)

    call DLVO(part%charges, part%positions, sys%box, part%radius, kappa, sys%force_DLVO, sys%pot_DLVO)
    call WCA(part%positions, sys%box, part%radius, eps, sys%force_WCA, sys%pot_WCA)

    sys%total_forces =  sys%force_WCA + sys%force_DLVO
    sys%total_potential =  sys%pot_WCA + sys%pot_DLVO 
    
    do step = 1, params%max_steps
        
        call u2_propagator(dt/2.0_dp, part%velocities, sys%total_forces)
        call langevin_step(dt/2.0_dp, part%velocities, part%masses, params%gamma, sys%temp_target)
        call u1_propagator(dt, part%positions, part%velocities, sys%box)
        
        call DLVO(part%charges, part%positions, sys%box, part%radius, kappa, sys%force_DLVO, sys%pot_DLVO)
        call WCA(part%positions, sys%box, part%radius, eps, sys%force_WCA, sys%pot_WCA)

        sys%total_forces =  sys%force_WCA + sys%force_DLVO
        sys%total_potential =  sys%pot_WCA + sys%pot_DLVO 
        
        call langevin_step(dt/2.0_dp, part%velocities, part%masses, params%gamma, sys%temp_target)
        call u2_propagator(dt/2.0_dp, part%velocities, sys%total_forces)
            
    	kinetic = 0.5_dp*sum(part%masses*sum(part%velocities**2, dim = 1))  
        sys%total_energy  = kinetic + sys%total_potential 
        temp = (2.0_dp*kinetic)/(real(sys%g, dp))
        
        if (mod(step, 100) == 0) then
            call write_pdb(sys, part, params%outdir)
            !call write_log(sys, part, params%outdir, step, temp, params%is_first)
            call write_out(step, sys%total_energy, temp)
            params%is_first = .false.
        end if
        
        call check_interruption(params%interrupt)
        
        if(step .eq. params%max_steps .or. params%interrupt) then
            call write_final_config(sys, part, params%outdir)
        end if
    end do
    
    call end_particles(part)
    call end_system(sys)

end program

    