close all;
%WARNING: NEEDS AN iv_timesereis to have the time_axis

%Load profiles

%profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Erlangen_2018\Matlab_profiles\'
% load([profile_folder, 'n_e_min_XYT']);
% load([profile_folder, 'T_e_min_XYT']);
% T_scale = [3 10];
% n_scale = [1 5]*1e19;

%profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_04_12\02_Averaging_after_fit\';
% x_range = 98:10:198;
% y_range = 65:10:165;

profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_05_18_X_Drive_OFF_Langmuir\';
x_range = 0:10:220;
y_range = 0:10:170;

%Load the timeseries data structure
load([profile_folder, 'iv_timeseries']);

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

T_scale = [0 20];
I_sat_i_scale = [0 90];
V_plasma_scale = [0 30];

% %Smooth the data a little bit...
space_smoothing = [3,3,1];
time_smoothing = [1,1,3];
spacetime_smoothing = 3;
% T_e_XYT = smooth3(T_e_XYT,'box',space_smoothing);
% I_sat_i_XYT = smooth3(I_sat_i_XYT,'box',space_smoothing);
V_plasma_XYT = smooth3(V_plasma_XYT,'box',space_smoothing);

n_isopotential = 30;

%snapshots = [6 ,12, 18];
%snapshots = [8];
%snapshots = [4 ,8, 12, 16];
snapshots = [8 ,12, 16, 20];
Ns_tot = length(snapshots);

f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');

%Load the E_gun CAD image
picture_dir = 'C:\Users\tzf\Desktop\VINETA II Lab\Presentations for Olaf\2018_03_20_PPT\';
file_name = 'E_gun_CAD.bmp';
complete_picture_path = [picture_dir,file_name];
E_gun_CAD_image = imread(complete_picture_path);
E_gun_CAD_image = E_gun_CAD_image(203:203+278,203:203+278);
scale_factor = 100/size(E_gun_CAD_image,1);
E_gun_CAD_image = imresize(E_gun_CAD_image,scale_factor);
alpha_data = ~E_gun_CAD_image;

for ns = 1:Ns_tot
    
    subplot(3,Ns_tot,ns);
    hold on;

    
    [C,h] = contourf(x_range,y_range,(T_e_XYT(:,:,snapshots(ns)))');
    %[C,h] = contourf(x_range,y_range,T_e_XYT(:,:,snapshots(ns)));
    set(h,'LineColor','none')
    
%     hp = imshow(E_gun_CAD_image,'InitialMagnification','fit');
%     set(hp, 'AlphaData', alpha_data);

    %Not exactly the same but ok for testing
    time_axis = 1e6.*(iv_timeseries.total_time_axis_rise);
    
    title(strcat('Electron Temperature profile @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,hot);
    caxis(T_scale);
    h = colorbar(gca);
    ylabel(h, 'T_e [eV]');
    
    hold off;
    
    subplot(3,Ns_tot,ns + Ns_tot);
    [~,h] = contourf(x_range,y_range,(I_sat_i_XYT(:,:,snapshots(ns)))');
    set(h,'LineColor','none')
    title(strcat('Ion saturation current @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,parula);
    caxis(I_sat_i_scale)
    h = colorbar(gca);
    ylabel(h, 'I_{sat,i} [mA]');
    
    hold off;
    
    subplot(3,Ns_tot,ns + 2*Ns_tot);
    [~,~] = contourf(x_range,y_range,(V_plasma_XYT(:,:,snapshots(ns)))',n_isopotential);
    %set(h,'LineColor','none')
    title(strcat('Plasma potential @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,cool);
    %caxis(V_plasma_scale)
    h = colorbar(gca);
    ylabel(h, 'V_{plasma} [V]');
    
end

f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
%movegui(gcf,'center');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

frame_rep = 10;
%Record video
frames = cell(length(time_axis)*frame_rep,1);

for ts = 1:length(time_axis)

%     subplot(2,2,[1,3]);
%     plot((out(1).t)*1e6 , out(1).data);
%     title('Electron Gun 1 Cathode Current');
%     xlabel(time_label_string);
%     ylabel('I_{Cathode} [A]');
%     %time_limits = xlim;
%     xlim(time_limits);
% 
%     SP= time_axis(ts); %your point goes here 
%     line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);


    subplot(2,3,4);
    [~,h] = contourf(x_range,y_range,(T_e_XYT(:,:,ts))');
    set(h,'LineColor','none')
    title(strcat('Electron Temperature @ t = ',num2str(time_axis(ts)),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,hot);
    caxis(T_scale);
    h = colorbar(gca);
    ylabel(h, 'T_e [eV]');

    subplot(2,3,5);
    [~,h] = contourf(x_range,y_range,(I_sat_i_XYT(:,:,ts))');
    set(h,'LineColor','none')
    title(strcat('Ion saturation current @ t = ',num2str(time_axis(ts)),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,parula);
    caxis(I_sat_i_scale)
    h = colorbar(gca);
    ylabel(h, 'I_{sat,i} [mA]');

    subplot(2,3,6);
    [~,h] = contourf(x_range,y_range,(V_plasma_XYT(:,:,ts))',n_isopotential);
    %set(h,'LineColor','none')
    title(strcat('Plasma potential @ t = ',num2str(time_axis(ts)),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 max(y_range)/max(x_range) 1]);
    colormap(gca,cool);
    %caxis(V_plasma_scale)
    h = colorbar(gca);
    ylabel(h, 'V_{plasma} [V]');

    drawnow;
   for frame_copy_index = 1:frame_rep
       frame_index = (ts-1)*frame_rep + frame_copy_index;
       frames{frame_index} = getframe(f);
   end

end

%Save video to file
frames = cell2mat(frames);
v = VideoWriter([profile_folder,'profiles_video.avi'],'Uncompressed AVI');

open(v);
writeVideo(v,frames);
close(v);