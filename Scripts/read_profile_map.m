clear all
close all

addpath('..\Functions\');

%Access grid scan directory

folder_path = '\\pc-e5-ws-2\D\tzf\2018_05_28_Plasma_Profile_050kHz_4V_5_rep_CompleteProfile_X-Drive_On_1000V\'
%Set the number of repetition per position
rep_per_point = 5;

map_file = 'map.map'
voltage_sign = +1;
current_sign = -1;
full_file_path = strcat(folder_path,map_file);

%Raster step
raster_step = 10;

%Read the map file
B = tdfread(full_file_path);
ind = B.x0x25_idx;
x_pos = B.x;
y_pos = B.y;
ind = ind(:,1);
x_pos = x_pos(:,1);
y_pos = y_pos(:,1);
%Number of points

%CORRECTION
x_pos_mod = 10*floor(x_pos/10);
x_pos = x_pos_mod;
%----------

NX = (max(x_pos) - min(x_pos))/raster_step + 1;
NY = (max(y_pos) - min(y_pos))/raster_step + 1;

%Organize data in multidimensional array
%Measurement number 2Dspace*1Drepetition array
meas_number_STR = NaN(NX,NY,rep_per_point);

for k_x = 1:NX
    for k_y = 1:NY
        for k_r = 1:rep_per_point
            %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
            %file_names_STR(k_x,k_y,repetition) = file_name;
            k = k_x + (k_y-1)*NX;
            meas_number_STR(k_x,k_y,k_r) = rep_per_point*(k-1)+(k_r-1);
        end
    end
end

%Select sub-grid to process
NX_begin = 1;
NX_end = NX;
nX = NX_end - NX_begin + 1;

NY_begin = 1;
NY_end = NY;
nY = NY_end - NY_begin + 1;

rep_begin = 1;
rep_end = 5;
nr = rep_end - rep_begin +1;

%FLoating potential 2Dspace*1Dtime*1Drepetition array
V_float_XYTR_rise = NaN(nX,nY,100,nr);
V_float_XYTR_fall = NaN(nX,nY,100,nr);
%Plasma potential 2Dspace*1Dtime*1Drepetition array
V_plasma_XYTR_rise = NaN(nX,nY,100,nr);
V_plasma_XYTR_fall = NaN(nX,nY,100,nr);
%Electron temperature 2Dspace*1Dtime*1Drepetition array
T_e_XYTR_rise = NaN(nX,nY,100,nr);
T_e_XYTR_fall = NaN(nX,nY,100,nr);
%Ion saturation current 2Dspace*1Dtime*1Drepetition array
I_sat_i_XYTR_rise = NaN(nX,nY,100,nr);
I_sat_i_XYTR_fall = NaN(nX,nY,100,nr);

for k_x = 1:nX
    for k_y = 1:nY
        for k_r = 1:nr
            file_number = meas_number_STR(k_x+NX_begin-1,k_y+NY_begin-1,k_r + rep_begin-1);
            filename = strcat('meas_',num2str(file_number,'%04.f'),'.h5')
            filename_cell_array = {filename};
            try
                [ return_status_read, iv_timeseries ] = readIVTimeseries( folder_path, filename_cell_array );
            catch err_read
                fprintf('IV Timeseries Read failed: jumping to next interation');
                continue
            end
            try
                [ return_status_extract, plasma_param_timeseries ] = extractPlasmaParamIVTimeseries(iv_timeseries);
            catch err_extract
                fprintf('Plasma Parameter extraction failed: jumping to next interation');
                continue
            end
        
            V_float_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_float;
            V_float_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_float;
            
            V_plasma_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_plasma;
            V_plasma_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_plasma; 
            
            T_e_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.T_e;
            T_e_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.T_e;
            
            I_sat_i_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.I_sat_i;
            I_sat_i_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.I_sat_i;
        end
    end
end

save_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_05_28_X_Drive_ON_1000V_Langmuir\';

V_float_XYT_rise = nanmean(V_float_XYTR_rise,4);
save([save_folder,'V_float_XYT_rise.mat'],'V_float_XYT_rise');
V_float_XYT_fall = nanmean(V_float_XYTR_fall,4);
save([save_folder,'V_float_XYT_fall.mat'],'V_float_XYT_fall');

V_plasma_XYT_rise = nanmean(V_plasma_XYTR_rise,4);
save([save_folder,'V_plasma_XYT_rise.mat'],'V_plasma_XYT_rise');
V_plasma_XYT_fall = nanmean(V_plasma_XYTR_fall,4);
save([save_folder,'V_plasma_XYT_fall.mat'],'V_plasma_XYT_fall');

T_e_XYT_rise = nanmean(T_e_XYTR_rise,4);
save([save_folder,'T_e_XYT_rise.mat'],'T_e_XYT_rise');
T_e_XYT_fall = nanmean(T_e_XYTR_fall,4);
save([save_folder,'T_e_XYT_fall.mat'],'T_e_XYT_fall');

I_sat_i_XYT_rise = nanmean(I_sat_i_XYTR_rise,4);
save([save_folder,'I_sat_i_XYT_rise.mat'],'I_sat_i_XYT_rise');
I_sat_i_XYT_fall = nanmean(I_sat_i_XYTR_fall,4);
save([save_folder,'I_sat_i_XYT_fall.mat'],'I_sat_i_XYT_fall');
