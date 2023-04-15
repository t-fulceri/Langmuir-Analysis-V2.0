clear;

addpath('..\Functions\');

folder_path = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_02_20_Plasma_Profile_050kHz_4V_5_rep_CentralZone_REPRISE\';
filename_cell_array = cell(1,1);
filename_cell_array{1} = 'meas_0300.h5';
% filename_cell_array{2} = 'meas_0496.h5';
% filename_cell_array{3} = 'meas_0497.h5';
% filename_cell_array{4} = 'meas_0498.h5';
% filename_cell_array{5} = 'meas_0499.h5';

[ return_status, iv_timeseries ] = readIVTimeseries( folder_path, filename_cell_array );

[return_status, plasma_param_timeseries] = extractPlasmaParamIVTimeseries(iv_timeseries);




close all;


% set(gcf,'Renderer','painters');
% axes_IV = axes;
% 
% for nc = 1:iv_timeseries.NC_tot
% 
%     V_axis_rise = iv_timeseries.rise(nc).V_axis;
%     I_of_V_rise = iv_timeseries.rise(nc).I_of_V;
%     I_of_V_sgolay_rise = plasma_param_timeseries.rise.I_of_V_sgolay{nc};
% 
%     V_axis_fall = iv_timeseries.fall(nc).V_axis;
%     I_of_V_fall = iv_timeseries.fall(nc).I_of_V;
%     I_of_V_sgolay_fall = plasma_param_timeseries.fall.I_of_V_sgolay{nc};
% 
% 
%     %MUST FIX: DIMENSIONS DO NOT AGREE IN GENERAL
%     %histeresis_measure = I_of_V_sgolay_rise - I_of_V_sgolay_fall
%     cla(axes_IV);
%     hold(axes_IV,'on');
%     
%     plot(axes_IV,V_axis_rise,I_of_V_rise,'or');
%     plot(axes_IV,V_axis_rise,I_of_V_sgolay_rise,'--r','LineWidth',2);
%     
%     plot(axes_IV,V_axis_fall,I_of_V_fall,'ok');
%     plot(axes_IV,V_axis_fall,I_of_V_sgolay_fall,'--k','LineWidth',2);
%     
%     hold(axes_IV,'off');
%     
%     legend('IV Datapoints (Rising Flank)','Savitzky-Golay Filtered Datapoints (Rising Flank)','IV Datapoints (Falling Flank)','Savitzky-Golay Filtered Datapoints (Falling Flank)');
%     
%     title(['IV-characteristic @ cycle number ',num2str(nc),'/',num2str(iv_timeseries.NC_tot)]);
%     ylabel('Current [I]');
%     xlabel('Voltage [V]');
%     
%     xlim([-80 80]);
%     ylim([-0.5 3]);
%     drawnow;
%     pause(2.0);
% end

% set(gcf,'Renderer','painters');
% hold on
% %Temperature timeseries extracted from rising flanks
% T_e_timeseries_rise = plasma_param_timeseries.rise.T_e;
% T_e_timeseries_lower_rise = plasma_param_timeseries.rise.T_e_lower;
% T_e_timeseries_upper_rise = plasma_param_timeseries.rise.T_e_upper;
% 
% time_axis_rise = iv_timeseries.total_time_axis_rise;
% time_axis_rise = time_axis_rise*1e6;
% delta_T_e_negative_rise = T_e_timeseries_rise-T_e_timeseries_lower_rise;
% delta_T_e_positive_rise = T_e_timeseries_upper_rise-T_e_timeseries_rise;
% 
% plot(time_axis_rise,T_e_timeseries_rise,'r*-');
% %errorbar(time_axis_rise,T_e_timeseries_rise,delta_T_e_negative_rise,delta_T_e_positive_rise,'r*-');
% 
% %Temperature timeseries extracted from falling flanks
% T_e_timeseries_fall = plasma_param_timeseries.fall.T_e;
% T_e_timeseries_lower_fall = plasma_param_timeseries.fall.T_e_lower;
% T_e_timeseries_upper_fall = plasma_param_timeseries.fall.T_e_upper;
% 
% time_axis_fall = iv_timeseries.total_time_axis_fall;
% time_axis_fall = time_axis_fall*1e6;
% delta_T_e_negative_fall = T_e_timeseries_fall-T_e_timeseries_lower_fall;
% delta_T_e_positive_fall = T_e_timeseries_upper_fall-T_e_timeseries_fall;
% 
% plot(time_axis_fall,T_e_timeseries_fall,'k*-');
% %errorbar(time_axis_fall,T_e_timeseries_fall,delta_T_e_negative_fall,delta_T_e_positive_fall,'k*-');
% 
% time_axis_mean = mean([time_axis_rise;time_axis_fall],1);
% T_e_timeseries_mean = mean([T_e_timeseries_rise;T_e_timeseries_fall],1);
% 
% plot(time_axis_mean,T_e_timeseries_mean,'b*-');
% 
% title('Electron Temperature Evolution Over Time');
% xlabel('Time [µs]');
% ylabel('T_e [eV]');
% 
% % legend('T_e extracted from rising edges','','T_e extracted from falling edges','','T_e mean from rising and falling edges');
% legend('T_e extracted from rising edges','T_e extracted from falling edges','T_e mean from rising and falling edges');
% hold off

nc = 8;
rise_fall_select = 1;

delta_I = iv_timeseries.I_ADC_step;

switch rise_fall_select
    case 1
        V_axis = iv_timeseries.rise(nc).V_axis;
        I_of_V = iv_timeseries.rise(nc).I_of_V;
        I_of_V_sgolay = plasma_param_timeseries.rise.I_of_V_sgolay{nc};
        time = iv_timeseries.total_time_axis_rise(nc);
        V_plasma = plasma_param_timeseries.rise.V_plasma(nc);
        I_of_V_fitted = plasma_param_timeseries.rise.I_of_V_fitted{nc};
        V_float = plasma_param_timeseries.rise.V_float(nc);
        I_i = plasma_param_timeseries.rise.I_i{nc};
        T_e = plasma_param_timeseries.rise.T_e(nc);
        I_sat_i = plasma_param_timeseries.rise.I_sat_i(nc);
    case 2
        V_axis = iv_timeseries.fall(nc).V_axis;
        I_of_V = iv_timeseries.fall(nc).I_of_V;
        I_of_V_sgolay = plasma_param_timeseries.fall.I_of_V_sgolay{nc};
        time = iv_timeseries.total_time_axis_fall(nc);
        V_plasma = plasma_param_timeseries.fall.V_plasma(nc);
        I_of_V_fitted = plasma_param_timeseries.fall.I_of_V_fitted{nc};
        V_float = plasma_param_timeseries.fall.V_float(nc);
        I_i = plasma_param_timeseries.fall.I_i{nc};
        T_e = plasma_param_timeseries.fall.T_e(nc);
        I_sat_i = plasma_param_timeseries.fall.I_sat_i(nc);
end

time = time*1e6;

if return_status == 0
    figure
    set(gcf,'Renderer','painters');
    ax_IV = axes;
    
    vec = ones(length(V_axis),1);
    errorbar(ax_IV,V_axis,I_of_V,delta_I.*vec,'ko');
    hold(ax_IV,'on');
    plot(ax_IV,V_axis,I_of_V_sgolay,'b-','LineWidth',2);
    set(ax_IV,'fontsize',28);
    switch rise_fall_select
        case 1
            title(ax_IV,['IV-characteristic at time = ',num2str(time),' µs (Rising edge)']);
        case 2
            title(ax_IV,['IV-characteristic at time = ',num2str(time),' µs (Falling edge)']);
    end
    xlabel(ax_IV,'Voltage [V]');
    ylabel(ax_IV,'Current [A]');
    xlim_temp = xlim;
    [~,V_plasma_ind] = min(abs(V_axis-V_plasma));
    plot_range = 1:V_plasma_ind;
    plot(ax_IV,V_axis(plot_range),I_of_V_fitted(plot_range),'r-','LineWidth',2);
    plot(ax_IV,V_axis,I_i,'r--','LineWidth',2);
    hold(ax_IV,'off');

    %Draw a vertical line at floating potential position
    axes(ax_IV);
    x_bound = [V_float V_float];
    y_bound = get(ax_IV,'YLim');
    line(x_bound,y_bound,'Color',[0 1 0],'LineWidth',2);

    %Draw a vertical line at plasma potential position
    axes(ax_IV);
    x_bound = [V_plasma V_plasma];
    y_bound = get(ax_IV,'YLim');
    line(x_bound,y_bound,'Color',[1 0 1],'LineWidth',2);

%     h_leg = legend(ax_IV,'Data points','Savitzky-Golay filtered data','Fit with model I(V) = I_0 + aV + be^{cV}','Ion current I_i(V) = I_0 + aV','Floating potential','Plasma potential');
%     set(h_leg,'FontSize',10);

    %cla(ax_IV);
%     axes(ax_IV);
%     descr = {'Plasma parameters from fit:';
%         ['V_{plasma} = ',num2str(V_plasma),' V'];
%         ['V_{float} = ',num2str(V_float),' V'];
%         ['T_e = ',num2str(T_e),' eV'];
%         ['I_{sat,i} = ',num2str(I_sat_i*1000),' mA']};
%     text(min(V_axis),max(I_of_V)*0.75,descr,'FontSize',12);

    % figure
    % set(gcf,'Renderer','painters');
    % movegui('northeast');
    % 
    % ax_logI = subplot(3,1,1);
    % plot(ax_logI,V_axis,logI);
    % title('log(I+I_{shift})');
    % xlim(xlim_temp);
    % 
    % ax_DlogI = subplot(3,1,2);
    % plot(ax_DlogI,V_axis,DlogI);
    % title('Dlog(I+I_{shift})');
    % 
    % ax_D2logI = subplot(3,1,3);
    % plot(ax_D2logI,V_axis,D2logI);
    % title('D2log(I+I_{shift})');
    % 
    % figure
    % set(gcf,'Renderer','painters');
    % movegui('southeast');
    % 
    % ax_I = subplot(3,1,1);
    % plot(ax_I,V_axis,I);
    % title('I = I(V)');
    % xlim(xlim_temp);
    % 
    % ax_DI = subplot(3,1,2);
    % plot(ax_DI,V_axis,DI);
    % title('I'' = I''(V)');
    % 
    % ax_D2I = subplot(3,1,3);
    % plot(ax_D2I,V_axis,D2I);
    % title('I'''' = I''''(V)');
    % x_bound = get(ax_IV,'XLim'); 
    % y_bound = [0 0];
    % line(x_bound,y_bound,'Color',[1 0 0]);
end