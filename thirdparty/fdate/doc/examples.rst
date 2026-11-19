========
Examples
========

This section provides practical examples of using FortranDate in scientific and engineering applications.

Basic Time Series Analysis
==========================

This example shows how to calculate statistics over timestamped data:

.. code-block:: fortran

   program time_series_analysis
      use mod_datetime
      implicit none
      
      integer, parameter :: num_samples = 24
      type(t_datetime) :: timestamps(num_samples)
      real :: measurements(num_samples)
      real :: hourly_avg
      integer :: i, hour
      
      ! Generate hourly data for one day
      do i = 1, num_samples
         timestamps(i) = t_datetime(2025, 5, 16, i-1, 0, 0)
         
         ! Simulate temperature measurements (higher during day, lower at night)
         hour = timestamps(i)%hour()
         if (hour >= 6 .and. hour <= 18) then
            measurements(i) = 20.0 + 5.0 * sin((hour - 6.0) * 3.14159 / 12.0)
         else
            measurements(i) = 15.0
         end if
      end do
      
      ! Print the data
      print *, "Time Series Data:"
      print *, "=================="
      do i = 1, num_samples
         print *, timestamps(i)%format("%H:%M"), measurements(i)
      end do
      
      ! Calculate average for a specific time window (10:00 to 14:00)
      type(t_datetime) :: window_start, window_end
      real :: sum_values
      integer :: count_values
      
      window_start = t_datetime(2025, 5, 16, 10, 0, 0)
      window_end = t_datetime(2025, 5, 16, 14, 0, 0)
      
      sum_values = 0.0
      count_values = 0
      
      do i = 1, num_samples
         if (timestamps(i) >= window_start .and. timestamps(i) <= window_end) then
            sum_values = sum_values + measurements(i)
            count_values = count_values + 1
         end if
      end do
      
      if (count_values > 0) then
         hourly_avg = sum_values / count_values
         print *, "Average between 10:00 and 14:00:", hourly_avg
      end if
      
   end program time_series_analysis

Simulation with Timesteps
=========================

This example demonstrates managing a simulation with fixed timesteps:

.. code-block:: fortran

   program timestep_simulation
      use mod_datetime
      implicit none
      
      ! Simulation parameters
      type(t_datetime) :: sim_start, sim_end, current_time
      type(t_timespan) :: timestep, elapsed
      real :: result, total
      integer :: step_count
      
      ! Initialize simulation time bounds
      sim_start = t_datetime(2025, 5, 16, 0, 0, 0)
      sim_end = t_datetime(2025, 5, 16, 6, 0, 0)
      
      ! Set timestep to 30 minutes
      timestep = t_timespan(minutes=30)
      
      ! Initialize simulation
      current_time = sim_start
      step_count = 0
      total = 0.0
      
      print *, "Starting simulation at", sim_start%format("%Y-%m-%d %H:%M:%S")
      print *, "Timestep:", timestep%total_minutes(), "minutes"
      
      ! Main simulation loop
      do while (current_time <= sim_end)
         ! Calculate some result for this timestep (example: temperature model)
         result = calculate_temperature(current_time)
         
         ! Accumulate results
         total = total + result
         step_count = step_count + 1
         
         ! Output current state
         print *, "Time:", current_time%format("%H:%M"), "Temperature:", result
         
         ! Advance to next timestep
         current_time = current_time + timestep
      end do
      
      ! Calculate elapsed simulation time
      elapsed = sim_end - sim_start
      
      ! Output summary
      print *, "Simulation complete"
      print *, "Elapsed simulation time:", elapsed%total_hours(), "hours"
      print *, "Number of timesteps:", step_count
      print *, "Average temperature:", total / step_count
      
   contains
   
      ! Example function to calculate temperature at a given time
      ! Simple model: Temperature varies sinusoidally with time of day
      function calculate_temperature(dt) result(temp)
         type(t_datetime), intent(in) :: dt
         real :: temp
         real :: hour_frac
         
         ! Convert hour and minute to a decimal hour
         hour_frac = real(dt%hour()) + real(dt%minute()) / 60.0
         
         ! Temperature model (min at 3AM, max at 3PM)
         temp = 15.0 + 10.0 * sin((hour_frac - 3.0) * 3.14159 / 12.0)
      end function calculate_temperature
      
   end program timestep_simulation

Processing Timestamped Experiment Data
======================================

This example shows how to process experimental data with timestamps:

.. code-block:: fortran

   program experiment_analysis
      use mod_datetime
      implicit none
      
      ! Data structure for experimental measurement
      type :: measurement_type
         type(t_datetime) :: timestamp
         real :: value
         character(len=16) :: sample_id
      end type measurement_type
      
      ! Sample data (would typically be loaded from a file)
      type(measurement_type), dimension(5) :: measurements = (/ &
         measurement_type(t_datetime(2025, 5, 15, 9, 30, 0), 23.5, "Sample-A"), &
         measurement_type(t_datetime(2025, 5, 15, 10, 45, 0), 24.1, "Sample-B"), &
         measurement_type(t_datetime(2025, 5, 15, 13, 15, 0), 24.8, "Sample-C"), &
         measurement_type(t_datetime(2025, 5, 15, 15, 0, 0), 23.9, "Sample-D"), &
         measurement_type(t_datetime(2025, 5, 15, 16, 30, 0), 22.7, "Sample-E") &
      /)
      
      ! Analyzing time spans between measurements
      type(t_timespan) :: time_diff, total_duration
      real :: avg_time_gap, values_per_hour
      integer :: i
      
      ! Calculate total experiment duration
      total_duration = measurements(5)%timestamp - measurements(1)%timestamp
      
      ! Calculate average gap between measurements
      avg_time_gap = 0.0
      do i = 2, 5
         time_diff = measurements(i)%timestamp - measurements(i-1)%timestamp
         avg_time_gap = avg_time_gap + time_diff%total_minutes()
      end do
      avg_time_gap = avg_time_gap / 4.0  ! 4 intervals for 5 measurements
      
      ! Calculate measurement frequency
      values_per_hour = 60.0 / avg_time_gap
      
      ! Output analysis
      print *, "Experiment Analysis"
      print *, "=================="
      print *, "Start time:", measurements(1)%timestamp%format("%Y-%m-%d %H:%M:%S")
      print *, "End time:", measurements(5)%timestamp%format("%Y-%m-%d %H:%M:%S")
      print *, "Total duration:", total_duration%total_hours(), "hours", &
               total_duration%minutes(), "minutes"
      print *, "Average time between measurements:", avg_time_gap, "minutes"
      print *, "Measurement frequency:", values_per_hour, "per hour"
      
      ! Find measurement nearest to a specific time
      type(t_datetime) :: target_time
      type(t_timespan) :: min_diff, current_diff
      integer :: nearest_idx
      
      target_time = t_datetime(2025, 5, 15, 14, 0, 0)  ! 2:00 PM
      min_diff = t_timespan(hours=24)  ! Start with a large value
      nearest_idx = 1
      
      do i = 1, 5
         ! Get absolute time difference
         current_diff = measurements(i)%timestamp - target_time
         
         ! Check if this measurement is closer than the previous best
         if (abs(current_diff%total_minutes()) < abs(min_diff%total_minutes())) then
            min_diff = current_diff
            nearest_idx = i
         end if
      end do
      
      print *, "Measurement closest to 14:00:"
      print *, "Sample:", measurements(nearest_idx)%sample_id
      print *, "Time:", measurements(nearest_idx)%timestamp%format("%H:%M")
      print *, "Value:", measurements(nearest_idx)%value
      print *, "Difference:", abs(min_diff%total_minutes()), "minutes"
      
   end program experiment_analysis

Simulation Checkpointing
========================

This example demonstrates using the datetime library for simulation checkpointing:

.. code-block:: fortran

   program checkpoint_simulation
      use mod_datetime
      implicit none
      
      ! Simulation parameters
      type(t_datetime) :: sim_start, current_time, next_checkpoint
      type(t_timespan) :: timestep, checkpoint_interval, runtime
      integer :: step_count, checkpoint_count
      real :: simulation_state
      character(len=64) :: checkpoint_filename
      
      ! Initialize simulation
      sim_start = now()  ! Current time when simulation starts
      checkpoint_interval = t_timespan(minutes=10)  ! Checkpoint every 10 minutes
      timestep = t_timespan(seconds=1)  ! 1-second timesteps
      
      ! Setup initial values
      current_time = sim_start
      next_checkpoint = sim_start + checkpoint_interval
      step_count = 0
      checkpoint_count = 0
      simulation_state = 0.0
      
      print *, "Starting simulation at", sim_start%format("%Y-%m-%d %H:%M:%S")
      
      ! Main simulation loop - run for 30 minutes of simulated time
      do while (step_count < 1800)  ! 30 minutes = 1800 seconds
         ! Update simulation (simple example: accumulate value)
         simulation_state = simulation_state + 0.1
         
         ! Check if it's time for a checkpoint
         if (current_time >= next_checkpoint) then
            checkpoint_count = checkpoint_count + 1
            
            ! Create checkpoint filename with timestamp
            write(checkpoint_filename, '(A,I2.2,A)') &
               "checkpoint_", checkpoint_count, "_" // &
               trim(current_time%format("%Y%m%d_%H%M%S")) // ".dat"
            
            ! Save checkpoint (in a real simulation, would write to file)
            call save_checkpoint(checkpoint_filename, step_count, simulation_state)
            
            ! Schedule next checkpoint
            next_checkpoint = current_time + checkpoint_interval
         end if
         
         ! Advance simulation time
         current_time = current_time + timestep
         step_count = step_count + 1
         
         ! In a real simulation, we'd do more work here
         if (mod(step_count, 100) == 0) then
            ! Calculate and print runtime every 100 steps
            runtime = now() - sim_start
            print *, "Step", step_count, "Complete. Runtime:", &
                     runtime%total_seconds(), "seconds"
         end if
      end do
      
      ! Calculate total runtime
      runtime = now() - sim_start
      
      ! Output summary
      print *, "Simulation complete"
      print *, "Total steps:", step_count
      print *, "Checkpoints created:", checkpoint_count
      print *, "Final simulation state:", simulation_state
      print *, "Total runtime:", runtime%total_seconds(), "seconds"
      
   contains
   
      ! Simulate saving a checkpoint file
      subroutine save_checkpoint(filename, step, state)
         character(len=*), intent(in) :: filename
         integer, intent(in) :: step
         real, intent(in) :: state
         
         print *, "Creating checkpoint:", trim(filename)
         print *, "   Step:", step
         print *, "   State:", state
      end subroutine save_checkpoint
      
   end program checkpoint_simulation
