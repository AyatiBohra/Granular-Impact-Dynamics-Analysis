clear functions
clear all
clc

format long;
C = {'k','b','r','g','m','c','y','k','b','r','g','m','c','y','k','b','r','g','m','c','y'};


%% Folder Path
newFolder = "C:\Users\ayati\Downloads\impact_data_r_10_15\impact_data_r_10_15";
cd(newFolder)

files = dir('*.xyz_all');
names = {files.name};

% Sort files numerically
nums = cellfun(@(x) str2double(regexp(x,'\d+','match','once')), names);
[~,order] = sort(nums);
files = files(order);

nof = length(files);

%% Time settings
freq = 500;
delT = 3.2*10^-6 * freq;

%% Preallocate
xi  = zeros(1,nof);
yi  = zeros(1,nof);
vxi = zeros(1,nof);
vyi = zeros(1,nof);
fyi = zeros(1,nof);
time = zeros(1,nof);

%   READ DATA

for ik = 1:nof

 % Read from 10th row
I1 = readmatrix(files(ik).name, 'FileType','text', 'NumHeaderLines', 9);
   
    % Extract values
xi(ik)  = I1(end,3);
yi(ik)  = I1(end,4);

vxi(ik) = I1(end,10);
vyi(ik) = I1(end,11);

fyi(ik) = I1(end,8);   
time(ik) = (ik-1)*delT;
   
end



%   CLEAN CONTINUOUS PLOT

[MM, II]=max(fyi);

%Defining grid size
grid_size = 250;

x_min_zoom = -7 * 0.1015; 
x_max_zoom =  7 * 0.1015; 

% Vertical: See 2 radii above (wake) and 6 radii below (compaction bulb)
y_min_zoom = -5 * 0.1015; 
y_max_zoom =  5 * 0.1015;


[X, Y] = meshgrid( ...
    linspace(x_min_zoom, x_max_zoom, 250), ... % Increased to 400 for sharpness
    linspace(y_min_zoom, y_max_zoom, 250));
start_frame = II + 1;
end_frame   = II + 8;

fprintf('Start frame: %d | File: %s | Time: %.6f s\n', ...
    start_frame, files(start_frame).name, time(start_frame));

fprintf('End frame: %d | File: %s | Time: %.6f s\n', ...
    end_frame, files(end_frame).name, time(end_frame));
Vx_sum = 0;
Vy_sum = 0;
count = 0;

for jj=start_frame:end_frame
   jm=jj;

 I2=readmatrix(files(jm).name, 'FileType','text', 'NumHeaderLines', 9);
 yintrudcen=I2(length(I2),4);
 xintrudcen=I2(length(I2),3);
 intrudrad=(I2(length(I2),6))/2;
 x_wrt_intrud=(I2(1:(length(I2)-1),3)-xintrudcen);
 y_wrt_intrud=(I2(1:(length(I2)-1),4)-yintrudcen);

 x_lab= I2(1:(length(I2)-1),3);
 y_lab=  I2(1:(length(I2)-1),4);
 
 vyintrud=-I2(length(I2),11);
 vxintrud=I2(length(I2),10);

 v_i = sqrt(vxintrud^2 + vyintrud^2);
 vx_wrt_intrud= (I2(1:(length(I2)-1),10))/v_i;
 vy_wrt_intrud=-(I2(1:(length(I2)-1),11))/v_i;

 vx_lab= I2(1:(length(I2)-1),10);
 vy_lab=I2(1:(length(I2)-1),11);
 
 rho_particle = 1280;
 type=I2(1:(length(I2)-1),2);
 mass = zeros(size(type));
 dia=I2(1:(length(I2)-1),6);
 thickness = 0.003;
 grain_indices = (type == 1 | type == 2);
 mass(grain_indices) = rho_particle * (pi * (dia(grain_indices).^2) / 6* thickness);
 



 
  diapmin = min(dia);
  sigma=1*diapmin;

  [Vx_grid_lab, Vy_grid_lab, mass_grid] = coarseGrainVelocityfast_same_dx_dy_lab(x_wrt_intrud, y_wrt_intrud, vx_wrt_intrud, vy_wrt_intrud,mass, dia, X ,Y, sigma);
   Vx_sum = Vx_sum + Vx_grid_lab;
   Vy_sum = Vy_sum + Vy_grid_lab;

   count = count + 1;
    
   
end  

Vx_plot_lab = Vx_sum / count;
Vy_plot_lab = Vy_sum / count;
r_grid = sqrt((X ).^2 + (Y ).^2);
mask_intruder = r_grid <= intrudrad;


Vx_plot_lab(mask_intruder) = 0;
Vy_plot_lab(mask_intruder) = 1;

figure

fields = {
    '$v_x / v_i$', Vx_plot_lab;
    '$v_y / v_i$', Vy_plot_lab;
};
cmap = turbo;

% --- Intruder position in LAB frame ---
theta=linspace(0,2*pi,300);
xcir=intrudrad*cos(theta);
ycir=intrudrad*sin(theta);

centx = xintrudcen;
centy = yintrudcen;
rr = intrudrad;

for i = 1:2
    ax = subplot(1,2,i);

    data = fields{i,2};

    % -------- CLEAN --------
    data(~isfinite(data)) = NaN;
    
    data = imgaussfilt(data, 0.8);
    data(mask_intruder) = fields{i,2}(mask_intruder);
    % -------- PLOT --------
   % plot (IMPORTANT CHANGE)
    imagesc(X(1,:), Y(:,1), data);
    
    set(gca,'YDir','normal');



    colormap(ax, turbo);
    colorbar(ax);
    caxis([-0.5 0.5]);

    axis equal tight
    xlim([-0.2 0.2])
    ylim([-0.18 0.15])

    title(fields{i,1}, 'Interpreter','latex');

    hold on

    % only boundary (no white fill)
    plot(xcir, ycir, 'k', 'LineWidth',2);
end

