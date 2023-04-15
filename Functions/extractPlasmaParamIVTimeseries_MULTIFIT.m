function [return_status, plasma_param_timeseries] = extractPlasmaParamIVTimeseries(iv_timeseries)
    
    fprintf('Running extractPlasmaParamIVTimeseries with following parameters:\n');
    
    %FUNCTION PARAMETERS
    %Polynomial grade for Savitzky-Golay filter
    n_sgolay = 5
    %Margin of data excluded from anlysis
    margin_percent = 0.1
    %Multiplication factor for the up-shifting of logI
    min_multiply_factor = 1.1
    %Polynomial grade for the approximation of logI
    n_poly = 10
    %Margins around the candidate focal point of fitting procedure
    delta_left  = 1000
    delta_right = 10
    %Maximum Temperature allowed [eV]
    T_e_max = 30
    %Minimum Temperature allowed [eV]
    T_e_min = 1
    %Voltage where to pick the ion saturation current [V]
    V_ion_saturation = -200
    
    NC_tot = length(iv_timeseries.rise);
    for rise_fall_select = 1:2
        for nc = 1:NC_tot
            fprintf(['Analyzing cycle number ',num2str(nc),'\n']);
            switch rise_fall_select
                case 1 %Rise
                    fprintf('(Rising flank)\n');
                    V_axis = iv_timeseries.rise(nc).V_axis;
                    I_of_V = iv_timeseries.rise(nc).I_of_V;
                case 2 %Fall
                    fprintf('(Falling flank)\n');
                    V_axis = iv_timeseries.fall(nc).V_axis;
                    I_of_V = iv_timeseries.fall(nc).I_of_V;
            end


            %FUNCTION PARAMETER
            L = length(V_axis);
            frame_length = round(L/4);
            if mod(frame_length,2) == 0
                frame_length = frame_length + 1;
            end

            %Apply Savitzky–Golay filter
            I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
            fprintf('Savitzky–Golay filter applied\n');

            %Focus on a meaningful range of indices (exclude margins)

            margin = round(margin_percent*length(V_axis));
            margin_range = margin:length(V_axis)-margin;

            if min(I_of_V_sgolay) < 0
                shift = -min_multiply_factor*min(I_of_V_sgolay);
                fprintf(['I(V) shifted upward, I_shift = ',num2str(shift),' A\n']);
            else
                shift = 0;
            end

            logI = log(I_of_V_sgolay + shift);
            fprintf('Logarithm of I=I(V) + I_shift taken\n');

            warning('off','all');
            p_logI = polyfit(V_axis,logI,n_poly);
            warning('on','all')
            fprintf(['Logarithm of I=I(V) + I_shift approximated with polynomial of grade ',num2str(n_poly),'\n']);

            logI_poly = polyval(p_logI,V_axis);

            DlogI = diff(logI_poly)./diff(V_axis);
            fprintf('Derivative of log(I+I_shfit) taken\n');

            %Extract floating potential and d
            [ d, ind_V_float ] = min(abs(I_of_V_sgolay));
            if ismember(ind_V_float,1:length(V_axis))
                V_float = V_axis(ind_V_float);
            else
                fprintf('Floating potential is out of range\n');
                V_float = nan;
            end
            
            %HERE WE PRODUCE A FIT FOR EVERY peak of DLogI

            % %NEWEST FIT-RANGE FINDER: look for the peak of Dlog nearest to V = 0
            [~,locs] = findpeaks(DlogI);
            locs = locs(locs > margin_range(1) & locs < margin_range(end));
            pks = DlogI(locs);

            [~,mp_ind_in_locs] = min(abs(V_axis(locs)));
            mp = pks(mp_ind_in_locs);


            %Apply the fitting procedure around the mp point

            %fit_range = intersect(find(DlogI > mp*0.8),locs(mp_ind_in_locs)-delta_left:locs(mp_ind_in_locs)+delta_right);
            data_range = locs(mp_ind_in_locs)-delta_left:locs(mp_ind_in_locs)+delta_right;
            fit_range = intersect(margin_range,data_range);
            if isempty(fit_range)
                fprintf('Fit range does not contain any indices \n');
                return
            end

            V_tf = V_axis(fit_range);
            I_tf = I_of_V_sgolay(fit_range);


            [xData, yData] = prepareCurveData( V_tf, I_tf );

            %FIT PARAMETERS (A LOT)
            ft = fittype( 'I_0 + a*V + b*exp(c*V)', 'independent', 'V', 'dependent', 'I' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [-Inf 0 0.001 0.01];
            opts.MaxFunEvals = 1000;
            opts.MaxIter = 1000;
            opts.Robust = 'Bisquare';
            opts.StartPoint = [0 0.0001 0.1 mp];
            opts.TolFun = 0.0001;
            opts.TolX = 0.0001;
            opts.Upper = [0 Inf 1 Inf];

            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            cv = coeffvalues(fitresult);
            ci = confint(fitresult);

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

            if T_fit > T_e_min && T_fit < T_e_max
                fprintf('Electron temperature is in range\n');
                %Plasma potential
                V_plasma = V_float + T_fit
                %Ion saturation current
                I_sat_i_fit = I_0 + a*V_ion_saturation
            else
                fprintf('Temperature out of range: plasma parameters discarded\n');

                T_fit = nan;
                T_upper = nan;
                T_lower = nan;

                V_plasma = nan;
                I_sat_i_fit = nan;
            end
            
            switch rise_fall_select
                case 1 %Rise
                    plasma_param_timeseries.rise.fitresult_cell_array{nc} = fitresult;

                    plasma_param_timeseries.rise.T_e(nc) = T_fit;
                    plasma_param_timeseries.rise.T_e_upper(nc) = T_upper;
                    plasma_param_timeseries.rise.T_e_lower(nc) = T_lower;


                    plasma_param_timeseries.rise.V_float(nc) = V_float;
                    plasma_param_timeseries.rise.V_plasma(nc) = V_plasma;
                    plasma_param_timeseries.rise.I_sat_i(nc) = I_sat_i_fit;
                case 2 %Fall
                    plasma_param_timeseries.fall.fitresult_cell_array{nc} = fitresult;

                    plasma_param_timeseries.fall.T_e(nc) = T_fit;
                    plasma_param_timeseries.fall.T_e_upper(nc) = T_upper;
                    plasma_param_timeseries.fall.T_e_lower(nc) = T_lower;


                    plasma_param_timeseries.fall.V_float(nc) = V_float;
                    plasma_param_timeseries.fall.V_plasma(nc) = V_plasma;
                    plasma_param_timeseries.fall.I_sat_i(nc) = I_sat_i_fit;
            end

        end
    end
    
    return_status = 0;

end