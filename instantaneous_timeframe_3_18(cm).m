clear functions
clear all
clc

format long;
C = {'k','b','r','g','m','c','y','k','b','r','g','m','c','y','k','b','r','g','m','c','y'};


%% Folder Path
newFolder = "C:\Users\ayati\Downloads\impact_data_r_3_18\impact_data_r_3_18";
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
figure(1)
plot(time,fyi,'color','b','marker','d')
hold on 
xlabel('Analysis Time (s)'); ylabel('Force Fy');
title('Impact Force Profile');
set(gca, 'Color', 'w');
% Display of max impact force
peak_file = files(II).name;
fprintf('Peak force frame: %d\n', II);
disp(['Max Fy: ', num2str(MM)])
disp(['Peak time        : ', num2str(time(II)), ' s'])
disp(['Peak force file: ', peak_file])
% DISPLAY OF TIME BEFORE AND AT IMPACT

time_before_impact_abs = time(1);
time_at_impact_abs     = time(2) ;
    
fprintf('--- Impact Timing Information ---\n');
fprintf('File 1 (Before): %s | Time: %.6f s\n', files(1).name, time_before_impact_abs);
fprintf('File 2 (Impact): %s  | Time: %.6f s\n', files(2).name, time_at_impact_abs);
fprintf('Time Step between files: %.6f s\n', time_at_impact_abs - time_before_impact_abs);
fprintf('---------------------------------\n');


%Defining grid size
grid_size = 250;

x_min_zoom = -7 * 0.0318; 
x_max_zoom =  7 * 0.0318; 

% Vertical: See 2 radii above (wake) and 6 radii below (compaction bulb)
y_min_zoom = -5 * 0.0318; 
y_max_zoom =  5 * 0.0318;


[X, Y] = meshgrid( ...
    linspace(x_min_zoom, x_max_zoom, 250), ... % Increased to 400 for sharpness
    linspace(y_min_zoom, y_max_zoom, 250));


for jj=II+1
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

 vx_wrt_intrud= (I2(1:(length(I2)-1),10)-vxintrud);
 vy_wrt_intrud=(I2(1:(length(I2)-1),11));

 vx_lab= I2(1:(length(I2)-1),10);
 vy_lab=I2(1:(length(I2)-1),11);
 
 rho_particle = 1280;
 type=I2(1:(length(I2)-1),2);
 mass = zeros(size(type));
 dia=I2(1:(length(I2)-1),6);
 thickness = 0.003;
 grain_indices = (type == 1 | type == 2);
 mass(grain_indices) = rho_particle * (pi * (dia(grain_indices).^2) / 6* thickness);
 

 v_i = sqrt(vxintrud^2 + vyintrud^2);

 
  diapmin = min(dia);
  sigma=1*diapmin;

  [Vx_grid_lab, Vy_grid_lab, mass_grid] = coarseGrainVelocityfast_same_dx_dy_lab(x_wrt_intrud, y_wrt_intrud, vx_wrt_intrud, vy_wrt_intrud,mass, dia, X ,Y, sigma);
  
    
   
end  


r_grid = sqrt((X ).^2 + (Y ).^2);
mask_intruder = r_grid <= intrudrad;

% RADIAL / ANGULAR
% =========================
dx = X ;
dy = Y ;
%r_grid(r_grid==0)=1e-8;

erx = dx./r_grid; ery = dy./r_grid;

Vx_plot_lab = Vx_grid_lab/(v_i) ;
Vy_plot_lab = -Vy_grid_lab/(v_i) ;

Vx_plot_lab(mask_intruder) = 0;
Vy_plot_lab(mask_intruder) = 1;

vr_plot = (Vx_plot_lab.*erx + Vy_plot_lab.*ery);
vtheta_plot = (Vx_plot_lab.*ery - Vy_plot_lab.*erx);



mask = r_grid <= intrudrad;
vr_plot(mask) = vyintrud / v_i;   % constant (±1)
vtheta_plot(mask) = 0;                     % no rotation insided;


%vr_intrud = vx_intrud .* (X./r_grid) + vy_intrud .* (Y./r_grid);
%vr_plot = vr;
%vr_plot(mask_intruder) = vr_intrud(mask_intruder);
%vtheta_intrud = vx_intrud .* (Y./r_grid) - vy_intrud .* (X./r_grid);
%vtheta_plot = vtheta;
%vtheta_plot(mask_intruder) = vtheta_intrud(mask_intruder);


figure(2)

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
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])

    title(fields{i,1}, 'Interpreter','latex');

    hold on

    % only boundary (no white fill)
    plot(xcir, ycir, 'k', 'LineWidth',2);
end

%% =========================
% PLOT ALL
% =========================
figure

fields = {
    
    '$v_r/v_i$', vr_plot;
    '$v_\theta/v_i$', vtheta_plot;
};

theta = linspace(0,2*pi,300);
xc = intrudrad*cos(theta);
yc =intrudrad*sin(theta);

for i=1:2
    subplot(1,2,i)

    data = fields{i,2};
    data = fillmissing(data,'linear');
    data = imgaussfilt(data,0.8);
    data(mask) = fields{i,2}(mask);

    imagesc(X(1,:),Y(:,1),data)
    %set(gca,'YDir','reverse')

    colormap(turbo)
    colorbar
    caxis([-0.5 0.5])

    axis equal tight
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    hold on
    plot(xc,yc,'k','LineWidth',2)

    title(fields{i,1},'Interpreter','latex')
end



%  CLARK-STYLE PIV


% smoothing 
Vx_s = imgaussfilt(Vx_grid_lab, 1);
Vy_s = imgaussfilt(Vy_grid_lab, 1);

% REMOVE BACKGROUND FLOW 

far_field = (abs(X) > 3*intrudrad) & (Y < -2*intrudrad);   % away from intruder
Vy_bg = nanmean(Vy_s(far_field));
Vx_bg = nanmean(Vx_s(far_field));

Vy_s = Vy_s - Vy_bg;
Vx_s = Vx_s - Vx_bg;

% --- Magnitude ---
Vmag = sqrt(Vx_s.^2 + Vy_s.^2);

% --- Mask intruder ---
r_grid = sqrt(X.^2 + Y.^2);
mask_intruder = r_grid <= intrudrad;

Vx_s(mask_intruder) = NaN;
Vy_s(mask_intruder) = NaN;

% PLOT
% =========================
figure;
hold on;

skip = 3;

Xq = X(1:skip:end, 1:skip:end);
Yq = Y(1:skip:end, 1:skip:end);
Vxq = Vx_s(1:skip:end, 1:skip:end);
Vyq = Vy_s(1:skip:end, 1:skip:end);

quiver(Xq, Yq, Vxq, Vyq, 1.5, 'b', 'LineWidth', 0.8);
theta = linspace(0,2*pi,300);
plot(intrudrad*cos(theta), intrudrad*sin(theta), 'k', 'LineWidth', 2);

axis equal;
xlim([-3*intrudrad 3*intrudrad]);
ylim([-4*intrudrad 1*intrudrad]);

set(gca,'YDir','normal');

xlabel('x');
ylabel('y');
title('PIV Flow Field');

box on;
set(gca, 'Color', 'w');

fprintf('processing frame: %d\n', jm);
fprintf('File name: %s\n', files(jm).name);
fprintf('Time: %.6f seconds\n', time(jm));