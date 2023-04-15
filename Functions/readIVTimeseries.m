function [ return_status, iv_timeseries ] = readIVTimeseries( folder_path, filename_cell_array )
%readIVTimeseries reads the voltage/current data from a list of measurement
%files corresponding to the same measuring conditions and returns a list of
%I = I(V) curves, each corresponding to a voltage cycle

    fprintf('Runnning readIVTimeseries\n');

    voltage_sign = +1;
    current_sign = -1; %(IV-Curve convention: invert current axis);

    iv_timeseries.N_files = length(filename_cell_array);
    fprintf(['Number of files: ',num2str(iv_timeseries.N_files),'\n']);

    fprintf(['Folder path: ',folder_path,'\n']);
    fprintf('File names:\n');
    for nf = 1:iv_timeseries.N_files
        fprintf([filename_cell_array{nf},'\n']);
    end

    NC_tot_array = nan(iv_timeseries.N_files,1);
    V_ADC_step_array = nan(iv_timeseries.N_files,1);
    I_ADC_step_array = nan(iv_timeseries.N_files,1);
    time_step_array = nan(iv_timeseries.N_files,1);
    time_steps_per_cycle_array = nan(iv_timeseries.N_files,1);

    %Collect all the relevant information to form the I_cycles and V_cycles
    %frome each measurement file
    for nf = 1:iv_timeseries.N_files

        full_file_path = strcat(folder_path,filename_cell_array{nf});
        %Check metadata
        metadata = h5info(full_file_path);
        %Extract timestep from metadata [s]
        time_step = metadata.Groups.Attributes(2,1).Value(1);
        time_step_array(nf) = time_step;
        %Extract date and time from metadata
        iv_timeseries.date_array(nf) = metadata.Datasets.Attributes(1,1).Value;
        iv_timeseries.time_array(nf) = metadata.Datasets.Attributes(2,1).Value;

        %Read data
        data = h5read(full_file_path,'/dataset');
        %DAQ delay in s
        %DAQ_delay = 20e-3;
        %Transpose
        data = data';
        %Extract voltage
        voltage = data(:,1);
        %Optional inversion of voltage signal (depends on the data source)
        voltage = voltage_sign.*voltage;
        %Extract current
        current = data(:,2);
        %Optional inversion of current signal (depends on the data source)
        current = current_sign.*current;
        %Extract uncertainty coming from ADC discretization;
        V_ADC_step_array(nf) = signalToADCStep(voltage);
        I_ADC_step_array(nf) = signalToADCStep(current);

        %Calculate sampling frequency [Hz]
        sampling_freq = 1/time_step;

        %Split the timeseries signal into chunks of equal duration
        %(each corresponding to one voltage cycle)
        [V_cycles,iv_timeseries.initial_time_steps_offset_array(nf)] =  timeseriesSplitSinefit(voltage,voltage,sampling_freq);
        I_cycles =  timeseriesSplitSinefit(voltage,current,sampling_freq);

        %Total number of cycles
        NC_tot_array(nf) = size(V_cycles,2);
        %Time steps per cycle
        time_steps_per_cycle_array(nf) = size(V_cycles,1);



        %Produce the time axis
        %time_axis = time_step*ts_per_cycle*conversion_factor*(linspace(1,NC_tot,NC_tot) - 1/2);
        %time_axis = time_axis + DAQ_delay*conversion_factor;


        %total_time_axis = (1/sampling_freq)*conversion_factor*linspace(0,length(voltage)-1,length(voltage));
        %total_time_axis = total_time_axis + DAQ_delay*conversion_factor;

        iv_cycles_struct_array(nf).cycle_length = size(V_cycles,1);
        hl = round(iv_cycles_struct_array(nf).cycle_length/2);

        iv_cycles_struct_array(nf).V_cycles_rise = V_cycles(1:hl,:);
        iv_cycles_struct_array(nf).I_cycles_rise = I_cycles(1:hl,:);

        iv_cycles_struct_array(nf).V_cycles_fall = V_cycles(hl+1:end,:);
        iv_cycles_struct_array(nf).I_cycles_fall = I_cycles(hl+1:end,:);

    end

    %CHECK IF ALL MEASUREMENT FILES YELD THE SAME NC_tot, V_ADC_step,
    %I_ADC_step, time_step
    if all(NC_tot_array == NC_tot_array(1))
        fprintf('All measurement files yeld the same number of cycles\n');
        iv_timeseries.NC_tot = NC_tot_array(1);
    else
        fprintf('WARNING!: Measurement files yeld different numbers of cycles: First value is selected\n');
        iv_timeseries.NC_tot = NC_tot_array(1);
    end

    if all(V_ADC_step_array == V_ADC_step_array(1))
        fprintf('All measurement files yeld the same value for Voltage ADC step\n');
        iv_timeseries.V_ADC_step = V_ADC_step_array(1);
    else
        fprintf('WARNING!: Measurement files yeld different value for Voltage ADC step: First value is selected\n');
        iv_timeseries.V_ADC_step = V_ADC_step_array(1);
    end

    if all(I_ADC_step_array == I_ADC_step_array(1))
        fprintf('All measurement files yeld the same value for Current ADC step\n');
        iv_timeseries.I_ADC_step = I_ADC_step_array(1);
    else
        fprintf('WARNING!: Measurement files yeld different value for Current ADC step: First value is selected\n');
        iv_timeseries.I_ADC_step = I_ADC_step_array(1);
    end

    if all(time_step_array == time_step_array(1))
        fprintf('All measurement files yeld the same value for the time step\n');
        iv_timeseries.time_step_in_seconds = time_step_array(1);
        iv_timeseries.sampling_freq = 1/time_step_array(1);
        iv_timeseries.time_steps_per_cycle = time_steps_per_cycle_array(1);
        iv_timeseries.time_uncertainty = round(iv_timeseries.time_steps_per_cycle/4)*time_step;
    else
        fprintf('WARNING!: Measurement files yeld different value for the time step: First value is selected\n');
        iv_timeseries.time_step_in_seconds = time_step_array(1);
        iv_timeseries.sampling_freq = 1/time_step_array(1);
        iv_timeseries.time_steps_per_cycle = time_steps_per_cycle_array(1);
        iv_timeseries.time_uncertainty = round(iv_timeseries.time_steps_per_cycle/4)*time_step;
    end
    if all(iv_timeseries.initial_time_steps_offset_array == iv_timeseries.initial_time_steps_offset_array(1))
        fprintf('All measurement files yeld the same value for the initial time step offset\n');
        iv_timeseries.initial_time_steps_offset = iv_timeseries.initial_time_steps_offset_array(1);
    else
        fprintf('WARNING!: Measurement files yeld different value for the initial time step offset: MODE value is selected\n');
        iv_timeseries.initial_time_steps_offset = mode(iv_timeseries.initial_time_steps_offset_array);
    end

    %Total number of time steps
    iv_timeseries.total_time_steps = iv_timeseries.initial_time_steps_offset + iv_timeseries.time_steps_per_cycle*iv_timeseries.NC_tot;
    fprintf('Total number of time steps extracted\n');
    %Total time axis for rising flank
    iv_timeseries.total_time_axis_rise = iv_timeseries.time_step_in_seconds*(iv_timeseries.initial_time_steps_offset + round(iv_timeseries.time_steps_per_cycle/4) + (0:iv_timeseries.NC_tot-1)*iv_timeseries.time_steps_per_cycle);
    fprintf('Total time axis for rising flank extracted\n');
    %Total time axis for falling flank
    iv_timeseries.total_time_axis_fall = iv_timeseries.time_step_in_seconds*(iv_timeseries.initial_time_steps_offset + round(iv_timeseries.time_steps_per_cycle*(3/4)) + (0:iv_timeseries.NC_tot-1)*iv_timeseries.time_steps_per_cycle);
    fprintf('Total time axis for falling flank extracted\n');
    
    %Stack together the datapoints and data uncertainties which correspond to
    %the same cycle
    iv_cycles_stack = struct('V_cycles_rise',[],'I_cycles_rise',[],'V_cycles_fall',[],'I_cycles_fall',[],'V_ADC_step_array',[],'I_ADC_step_array',[]);
    for nf = 1:iv_timeseries.N_files
        iv_cycles_stack.V_cycles_rise = cat(1,iv_cycles_stack.V_cycles_rise,iv_cycles_struct_array(nf).V_cycles_rise);
        iv_cycles_stack.I_cycles_rise = cat(1,iv_cycles_stack.I_cycles_rise,iv_cycles_struct_array(nf).I_cycles_rise);
        iv_cycles_stack.V_cycles_fall = cat(1,iv_cycles_stack.V_cycles_fall,iv_cycles_struct_array(nf).V_cycles_fall);
        iv_cycles_stack.I_cycles_fall = cat(1,iv_cycles_stack.I_cycles_fall,iv_cycles_struct_array(nf).I_cycles_fall);
    end
    fprintf('Datapoints have been stacked together\n');

    %Sort the V(t) and I(t) signals in every cycle following the V-order
    for nc = 1:iv_timeseries.NC_tot
        [iv_cycles_stack.V_cycles_rise(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_rise(:,nc));
        iv_cycles_stack.I_cycles_rise(:,nc) = iv_cycles_stack.I_cycles_rise(sort_ind,nc);

        [iv_cycles_stack.V_cycles_fall(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_fall(:,nc));
        iv_cycles_stack.I_cycles_fall(:,nc) = iv_cycles_stack.I_cycles_fall(sort_ind,nc);
    end
    %MAYBE IT CAN BE DONE FASTER...
    % [iv_cycles_stack.V_cycles_rise,sort_ind] = sort(iv_cycles_stack.V_cycles_rise,1);
    % iv_cycles_stack.I_cycles_rise = iv_cycles_stack.I_cycles_rise(sort_ind);
    fprintf('V(t) and I(t) signals sorted according to V-order\n');


    for nc = 1:iv_timeseries.NC_tot

        rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_rise(:,nc),iv_cycles_stack.I_cycles_rise(:,nc));
        iv_timeseries.rise(nc).V_axis = rec_curve_out.x;
        iv_timeseries.rise(nc).I_of_V = rec_curve_out.y_of_x;

        rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_fall(:,nc),iv_cycles_stack.I_cycles_fall(:,nc));
        iv_timeseries.fall(nc).V_axis = rec_curve_out.x;
        iv_timeseries.fall(nc).I_of_V = rec_curve_out.y_of_x;

    end
    fprintf('IV-timeseries ready\n');

    return_status = 0;

end

