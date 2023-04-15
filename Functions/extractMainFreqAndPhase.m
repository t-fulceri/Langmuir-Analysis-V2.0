function [ main_freq, main_phase ] = extractMainFreqAndPhase( signal, freq_sampling )
%EXTRACT MAIN FREQUENCY AND PHASE Applies FFT to a signal and returns the
%frequency and phase of the most important harmonic contribution
%   Tiziano Fulceri 2018-01-29
    
    % Sampling frequency
    switch nargin
        case 1
            Fs = 1e7;  %DEFAULT VALUE (10 MHz) 
        case 2
            Fs = freq_sampling;
    end
    
    % Length of signal
    L = length(signal);             

    %Prepare the array of frequencies
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y

    % Use only positive frequencies as reference
    freq = Fs/2*linspace(0,1,NFFT/2+1);

    %Perform the Fast Fourier Transform using the array of frequencies
    fft_signal = fft(signal,NFFT)/L;
    
    %Consider only positive frequencies
    fft_positive = fft_signal(1:NFFT/2+1);
    
    %For every +/- frequency pair, sum their amplitudes
    double_abs_fft_signal = 2*abs(fft_positive);
    
    %Look for maximum amplitude in the spectrum
    [M,M_ind] = max(double_abs_fft_signal);
    %Find corresponding frequency
    main_freq = freq(M_ind);
    %Find corresponding phase
    main_phase = angle(fft_positive(M_ind));
    
%     %PLOTTING (JUST FOR DEBUG)
%     % Plot single-sided amplitude distribution in frequency space
%     f_fft_signal = figure;
%     movegui(f_fft_signal,'north');
%     ax_fft_signal = axes;
%     double_abs_fft_signal_normalized = double_abs_fft_signal./M;
%     loglog(ax_fft_signal,freq,double_abs_fft_signal_normalized,'.-');
%     title(ax_fft_signal,'Normalized amplitude distribution (AKA spectrum)');
%     xlabel(ax_fft_signal,'Frequency [Hz]');
%     ylabel(ax_fft_signal,'Normalized Amplitude');  
%     
%     % Plot single-sided phase distribution in frequency space
%     f_phase = figure;
%     movegui(f_phase,'northeast');
%     ax_phase = axes;
%     semilogx(ax_phase,freq,angle(fft_positive),'.-');
%     title(ax_phase,'Phase distribution');
%     xlabel(ax_phase,'Frequency [Hz]');
%     ylabel(ax_phase,'Phase');
%     
%     % Plot original signal plus extracted main component
%     f_supimp = figure;
%     movegui(f_supimp,'center');
%     ax_supimp = axes;
%     hold on;
%     plot(signal);
%     time_step = (1/Fs);
%     t = time_step*linspace(0,L-1,L);
%     main_component = M*cos(2*pi*main_freq*t + main_phase);
%     plot(main_component);
%     hold off;
%     title(ax_supimp,'Original signal superimposed with main Fourier component');
%     xlabel(ax_supimp,'Time [timesteps]');
%     ylabel(ax_supimp,'Phase');

end

