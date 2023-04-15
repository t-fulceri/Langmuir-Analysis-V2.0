function [ signal_cycles, total_output_offset ] = timeseriesSplitSinefit(reference, signal, sampling_freq)
%timeseries_split_sinefit Summary of this function goes here
%   Detailed explanation goes here

%Keep note of important values
L = length(signal);
time_step = 1/sampling_freq;
t = time_step*linspace(0,L-1,L);

%Some starting parameters
[ main_freq, main_phase ] = extractMainFreqAndPhase( reference, sampling_freq );

pulsation_start = main_freq*2*pi*time_step;

[xData, yData] = prepareCurveData( [], reference );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [max(reference) pulsation_start 2.5];




% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
cv = coeffvalues(fitresult);
A = cv(1);
omega = cv(2)/time_step;
phi = cv(3);

fitted_oscillation = A*sin(omega*t + phi);

%Calculate cycle length
cycle_length = round(sampling_freq*2*pi/omega);

%Find first anti-peak (which marks the beginning of a cycle)    
[neg_peaks, neg_locs] = findpeaks(-fitted_oscillation);
offset = neg_locs(1) - 1;
%Cut untill the beginning of the first complete cycle
if offset ~= 0 && offset < cycle_length
    signal = signal(offset+1:end);
end

total_output_offset = offset;
%Calculate total number of cycles
N = floor(length(signal)./cycle_length);
%Round signal length to the nearest multiple of cycle_length
signal = signal(1:cycle_length*N);

%Perform the splitting
%(One long periodic signal ---> N cycles)
%Everything ordered in a matrix
signal_cycles = reshape(signal,cycle_length,N);

end

