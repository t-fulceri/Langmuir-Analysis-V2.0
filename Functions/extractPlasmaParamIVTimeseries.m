function [return_status, plasma_param_timeseries] = extractPlasmaParamIVTimeseries(iv_timeseries)
    
    fprintf('Running extractPlasmaParamIVTimeseries with following parameters:\n\n');
    return_status = 1;
    
    %FUNCTION PARAMETERS
    %Frame length fraction for Savitzky-Golay filter
    frame_length_fraction = 1/4
    %Polynomial grade for Savitzky-Golay filter
    n_sgolay = 5
    %Plasma detector: max/abs(min) threshold
    max_over_abs_min_th = 4;
    %Plasma detector: max(abs)/delta_current threshold
    max_abs_over_delta_current_th = 4;
    %Multiplication factor for the up-shifting of logI
    min_multiply_factor = 1.1
    %Polynomial grade for the approximation of I fixed to 5
    %Polynomial grade for the approximation of logI fixed to 5
    %Maximum Temperature allowed [eV]
    T_e_max = 30
    %Minimum Temperature allowed [eV]
    T_e_min = 1
    %Minimum absolute magnitude of Ion Saturation Current allowed [A]
    I_sat_i_abs_min = 0.001;
    %Voltage where to pick the ion saturation current [V]
    V_i_sat_leftward_offset = 100

    
    NC_tot = length(iv_timeseries.rise);
    delta_I = iv_timeseries.I_ADC_step;
    
    %Output data structure initialization
    for nc = 1:NC_tot
        %Output parameters from rising edges
        plasma_param_timeseries.rise.I_of_V_sgolay{nc} = [];
        plasma_param_timeseries.rise.fitresult_cell_array{nc} = [];
        plasma_param_timeseries.rise.gof_cell_array{nc} = [];
        plasma_param_timeseries.rise.I_of_V_fitted{nc} = [];
        plasma_param_timeseries.rise.I_i{nc} = [];

        plasma_param_timeseries.rise.T_e(nc) = nan;
        plasma_param_timeseries.rise.T_e_upper(nc) = nan;
        plasma_param_timeseries.rise.T_e_lower(nc) = nan;

        plasma_param_timeseries.rise.V_float(nc) = nan;
        plasma_param_timeseries.rise.V_plasma(nc) = nan;
        plasma_param_timeseries.rise.I_sat_i(nc) = nan;

        %Output parameters from falling edges
        plasma_param_timeseries.fall.I_of_V_sgolay{nc} = [];
        plasma_param_timeseries.fall.fitresult_cell_array{nc} = [];
        plasma_param_timeseries.fall.gof_cell_array{nc} = [];
        plasma_param_timeseries.fall.I_of_V_fitted{nc} = [];
        plasma_param_timeseries.fall.I_i{nc} = [];

        plasma_param_timeseries.fall.T_e(nc) = nan;
        plasma_param_timeseries.fall.T_e_upper(nc) = nan;
        plasma_param_timeseries.fall.T_e_lower(nc) = nan;

        plasma_param_timeseries.fall.V_float(nc) = nan;
        plasma_param_timeseries.fall.V_plasma(nc) = nan;
        plasma_param_timeseries.fall.I_sat_i(nc) = nan;
        
        %Rise-fall-averaged output parameters
        plasma_param_timeseries.mean.T_e(nc) = nan;

        plasma_param_timeseries.mean.V_float(nc) = nan;
        plasma_param_timeseries.mean.V_plasma(nc) = nan;
        plasma_param_timeseries.mean.I_sat_i(nc) = nan;
    end
    
    %Main cycle
    %For every voltage cycle
    for  nc = 1:NC_tot
        fprintf(['\n\nAnalyzing cycle number ',num2str(nc),'/',num2str(NC_tot),'\n\n']);
        
        %For every rising and faling edge
        for rise_fall_select = 1:2
            
            switch rise_fall_select
                case 1 %Rise
                    fprintf('(Rising edge)\n');
                    V_axis = iv_timeseries.rise(nc).V_axis;
                    I_of_V = iv_timeseries.rise(nc).I_of_V;
                case 2 %Fall
                    fprintf('(Falling edge)\n');
                    V_axis = iv_timeseries.fall(nc).V_axis;
                    I_of_V = iv_timeseries.fall(nc).I_of_V;
            end

            L = length(V_axis);
            frame_length = round(frame_length_fraction*L);
            if mod(frame_length,2) == 0
                frame_length = frame_length + 1;
            end

            %Apply Savitzky–Golay filter
            I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
            fprintf('Savitzky–Golay filter applied\n');

            %Plasma detector:
            if (max(I_of_V_sgolay)/abs(min(I_of_V_sgolay)) < max_over_abs_min_th && max(abs(I_of_V_sgolay))/delta_I < max_abs_over_delta_current_th)
                fprintf('No plasma detected: abort parameter extraction\n');
                %Exit from the current iteration in 1:NC_tot
                %and go to the next
                continue;
            else
                fprintf('Plasma detected\n');
            end
            
            %Polynomial approximations of I(V), I'(V), I''(V)
            n_poly = 5;
            p_I = polyfit(V_axis,I_of_V_sgolay,n_poly);
            p_DI = polyder(p_I);
            p_D2I = polyder(p_DI);
            
            %Shift the I = I(V) curve upward to avoid negative values in the logarithm
            if min(I_of_V_sgolay) < 0
                shift = -min_multiply_factor*min(I_of_V_sgolay);
                fprintf(['I(V) shifted upward, I_shift = ',num2str(shift),' A\n']);
            else
                shift = 0;
            end
            
            %Polynomial approximation of log(I+I_shift) and log(I+I_shift)'
            logI = log(I_of_V_sgolay + shift);
            fprintf('Logarithm of I=I(V) + I_shift taken\n');
            n_poly_log = 5;
            p_logI = polyfit(V_axis,logI,n_poly_log);
            fprintf(['Logarithm of I=I(V) + I_shift approximated with polynomial of degree ',num2str(n_poly_log),'\n']);
            p_DlogI = polyder(p_logI);
            fprintf('Polynomial derivative of log(I+I_shfit) taken\n');

            %Extract floating potential
            [ ~, ind_V_float ] = min(abs(I_of_V_sgolay));
            if ismember(ind_V_float,1:length(V_axis))
                V_float = V_axis(ind_V_float);
                fprintf(['Floating potential extracted: V_float = ',num2str(V_float),' V\n']);
            else
                fprintf('Floating potential is out of range: plasma parameter extraction aborted\n');
                continue
            end
            
            %Extract plasma potential (a root of D2I)
            r = roots(p_D2I);
            r = r(r > min(V_axis) & r < max(V_axis) & r >= V_float & r >= 0);
            if isempty(r)
                fprintf('I''''(V) has no roots greater than V_float: plasma parameter extraction aborted\n');
                continue
            end
            V_plasma = min(r);
            fprintf(['Plasma potential extracted: V_plasma = ',num2str(V_plasma),' V\n']);
            [~,V_plasma_ind] = min(abs(V_axis-V_plasma));
            
            %Look for peaks of DLogI
            DlogI = polyval(p_DlogI,V_axis);
            pks = findpeaks(DlogI);
            if ~isempty(pks)
                mp = pks(1);
            else
                fprintf('Derivative of log(I+I_shift) has no peaks: parameter extraction aborted\n');
                %Exit from the current iteration in 1:NC_tot
                %and go to the next
                continue;
            end

            %NEW FIT RANGE: from beginning of V_axis to plasma potential
            fit_range = 1:V_plasma_ind;
            if isempty(fit_range)
                fprintf('Fit range does not contain any indices \n');
                %Exit from the current iteration in 1:NC_tot
                %and go to the next
                continue;
            else
                fprintf('Fit range indentified\n');
            end
            
            V_tf = V_axis(fit_range);
            I_tf = I_of_V_sgolay(fit_range);
            
            fprintf('Starting fitting procedure...\n');
            [xData, yData] = prepareCurveData( V_tf, I_tf );

            %FIT PARAMETERS (A LOT)
            ft = fittype( 'I_0 + a*V + b*exp(c*V)', 'independent', 'V', 'dependent', 'I' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [-Inf 0 0.001 0.01];
            opts.MaxFunEvals = 600;
            opts.MaxIter = 600;
            opts.Robust = 'Bisquare';
            opts.StartPoint = [0 0.0001 0.1 mp];
            opts.TolFun = 0.0001;
            opts.TolX = 0.0001;
            opts.Upper = [0 Inf 1 Inf];

            % Fit model to data.
            try
                [fitresult, gof] = fit( xData, yData, ft, opts );
            catch ex
                fprintf('Fit procedure failed: plasma parameter extraction aborted');
                continue
            end
            cv = coeffvalues(fitresult);
            ci = confint(fitresult);
            fprintf('Fitting procedure successful!\n');

            I_0 = cv(1);
            a = cv(2);
            b = cv(3);
            c = cv(4);

            %Plasma parameters:
            
            %Electron Temperature
            T_fit = 1/c;
            fprintf(['Electron temperature extracted: T_e = ',num2str(T_fit),' eV\n']);
            T_upper = 1/ci(1,4);
            if T_upper < T_fit
                T_upper = T_fit;
            end
            T_lower = 1/ci(2,4);
            if T_lower > T_fit
                T_lower = T_fit;
            end              


            %Ion saturation current
            I_sat_i_fit = I_0 + a*(V_float - V_i_sat_leftward_offset);
            fprintf(['Ion saturation current extracted: I_sat_i_fit = ',num2str(I_sat_i_fit),' A\n']);

            if T_fit > T_e_min && T_fit < T_e_max
                fprintf('Electron temperature is in range\n');
            else
                fprintf('Temperature out of range: plasma parameters discarded\n');
                %Exit from the current iteration in 1:NC_tot
                %and go to the next
                continue;
            end

            if abs(I_sat_i_fit) > I_sat_i_abs_min
                fprintf('Ion Saturation Current is in range\n');
            else
                fprintf('Ion Saturation Current is out of range: plasma parameters discarded\n');
                %Exit from the current iteration in 1:NC_tot
                %and go to the next
                continue;
            end
            
            %If arrived to this point, it means that parameter extraction
            %went good during this iteration: assign the extracted parameters to the output
            %structure
            fprintf('Parameter extraction successful, assigning parameters values to output structure...\n');
            switch rise_fall_select
                case 1 %Rise
                    plasma_param_timeseries.rise.I_of_V_sgolay{nc} = I_of_V_sgolay;
                    plasma_param_timeseries.rise.fitresult_cell_array{nc} = fitresult;
                    plasma_param_timeseries.rise.gof_cell_array{nc} = gof;
                    plasma_param_timeseries.rise.I_of_V_fitted{nc} = I_0 + a.*V_axis + b.*exp(c.*V_axis);
                    plasma_param_timeseries.rise.I_i{nc} = I_0 + a.*V_axis;
                    
                    plasma_param_timeseries.rise.T_e(nc) = T_fit;
                    plasma_param_timeseries.rise.T_e_upper(nc) = T_upper;
                    plasma_param_timeseries.rise.T_e_lower(nc) = T_lower;

                    plasma_param_timeseries.rise.V_float(nc) = V_float;
                    plasma_param_timeseries.rise.V_plasma(nc) = V_plasma;
                    plasma_param_timeseries.rise.I_sat_i(nc) = I_sat_i_fit;
                case 2 %Fall
                    plasma_param_timeseries.fall.I_of_V_sgolay{nc} = I_of_V_sgolay;
                    plasma_param_timeseries.fall.fitresult_cell_array{nc} = fitresult;
                    plasma_param_timeseries.fall.gof_cell_array{nc} = gof;
                    plasma_param_timeseries.fall.I_of_V_fitted{nc} = I_0 + a.*V_axis + b.*exp(c.*V_axis);
                    plasma_param_timeseries.fall.I_i{nc} = I_0 + a.*V_axis;

                    plasma_param_timeseries.fall.T_e(nc) = T_fit;
                    plasma_param_timeseries.fall.T_e_upper(nc) = T_upper;
                    plasma_param_timeseries.fall.T_e_lower(nc) = T_lower;
                    
                    plasma_param_timeseries.fall.V_float(nc) = V_float;
                    plasma_param_timeseries.fall.V_plasma(nc) = V_plasma;
                    plasma_param_timeseries.fall.I_sat_i(nc) = I_sat_i_fit;
            end
            fprintf('Done\n');
        end
    end
    
    %Produce the rise-fall-averaged values for the output parameters
    
    V_float_timeseries_rise = plasma_param_timeseries.rise.V_float;
    V_float_timeseries_fall = plasma_param_timeseries.fall.V_float;
    plasma_param_timeseries.mean.V_float = mean([V_float_timeseries_rise;V_float_timeseries_fall],1);
        
    V_plasma_timeseries_rise = plasma_param_timeseries.rise.V_plasma;
    V_plasma_timeseries_fall = plasma_param_timeseries.fall.V_plasma;
    plasma_param_timeseries.mean.V_plasma = mean([V_plasma_timeseries_rise;V_plasma_timeseries_fall],1);
    
    I_sat_i_timeseries_rise = plasma_param_timeseries.rise.I_sat_i;
    I_sat_i_timeseries_fall = plasma_param_timeseries.fall.I_sat_i;
    plasma_param_timeseries.mean.I_sat_i = mean([I_sat_i_timeseries_rise;I_sat_i_timeseries_fall],1);
    
    T_e_timeseries_rise = plasma_param_timeseries.rise.T_e;
    T_e_timeseries_fall = plasma_param_timeseries.fall.T_e;
    plasma_param_timeseries.mean.T_e = mean([T_e_timeseries_rise;T_e_timeseries_fall],1);
    
 
    
    %If arrived untill here it means that everything went good:
    %set the return status properly
    fprintf('\n\nextractPlasmaParamIVTimeseries executed successfully\n');
    return_status = 0;

end

