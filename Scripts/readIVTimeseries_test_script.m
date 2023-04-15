clear all;
close all;

addpath('..\Functions\');

%folder_path = '\\pc-e5-ws-2\D\tzf\Plasma_Profile_050kHz_4V_5_rep_CentralZone_REPRISE\';
folder_path = '\\pc-e5-ws-2\D\tzf\2018_04_19_Plasma_Profile_050kHz_4V_5_rep_CompleteProfile_X-Drive_OFF\'
filename_cell_array = cell(1,1);
filename_cell_array{1} = 'meas_0305.h5';
% filename_cell_array{2} = 'meas_0306.h5';
% filename_cell_array{3} = 'meas_0307.h5';
% filename_cell_array{4} = 'meas_0308.h5';
% filename_cell_array{5} = 'meas_0309.h5';

[ return_status, iv_timeseries ] = readIVTimeseries( folder_path, filename_cell_array )

nc_selected = [6,9,12,18,21,24];

time_axis = 1e6.*(iv_timeseries.total_time_axis_rise + iv_timeseries.total_time_axis_fall)/2;

n_sub_plot = 1;
for nc = nc_selected
    
    V_axis_rise = iv_timeseries.rise(nc).V_axis;
    I_of_V_rise = iv_timeseries.rise(nc).I_of_V;

    V_axis_fall = iv_timeseries.fall(nc).V_axis;
    I_of_V_fall = iv_timeseries.fall(nc).I_of_V;
    
    delta_I = iv_timeseries.I_ADC_step;
    subplot(2,3,n_sub_plot);
    set(gcf,'Renderer','painters');
    errorbar(V_axis_rise,I_of_V_rise,ones(length(V_axis_rise),1).*delta_I,'or');
    hold on;
    errorbar(V_axis_fall,I_of_V_fall,ones(length(V_axis_fall),1).*delta_I,'ok');
    hold off;
    title(['IV_characteristic @ time = ',num2str(time_axis(nc)),' µs']);
    xlim([-80,+80]);
    ylim([-0.1,3]);
    xlabel('Bias voltage [V]');
    ylabel('Probe current [A]');
    legend('Data points from RISING edge','Datapoints from FALLING edge');
    
    n_sub_plot = n_sub_plot + 1;

end