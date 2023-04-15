close all;
%WARNING: NEEDS AN iv_timesereis to have the time_axis

%Load profiles

profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_04_12\02_Averaging_after_fit\';
load([profile_folder, 'T_e_XYT_rise']);
load([profile_folder, 'T_e_XYT_fall']);
load([profile_folder, 'I_sat_i_XYT_rise']);
load([profile_folder, 'I_sat_i_XYT_fall']);
load([profile_folder, 'V_plasma_XYT_rise']);
load([profile_folder, 'V_plasma_XYT_fall']);

[size_x,size_y,size_t] = size(T_e_XYT_rise);

T_e_XYT = nan(size_x,size_y,size_t);
I_sat_i_XYT = nan(size_x,size_y,size_t);
V_plasma_XYT = nan(size_x,size_y,size_t);

for ix = 1:size_x
    for iy = 1:size_y
        for it = 1:size_t
            T_e_XYT(ix,iy,it) = (T_e_XYT_rise(ix,iy,it) + T_e_XYT_fall(ix,iy,it))/2;
            
            %REVERSE SIGN FOR ION SATURATION CURRENT!
            I_sat_i_XYT(ix,iy,it) = -1e3*(I_sat_i_XYT_rise(ix,iy,it) + I_sat_i_XYT_fall(ix,iy,it))/2;
            
            V_plasma_XYT(ix,iy,it) = real((V_plasma_XYT_rise(ix,iy,it) + V_plasma_XYT_fall(ix,iy,it))/2);
        end
    end
end

%Reconstruct surface from gridded data
F = scatteredInterpolant((1:size_x)',(1:size_y)',(1:size_t)',T_e_XYT)
