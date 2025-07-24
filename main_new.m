clear all; close all; clc; 
addpath(fullfile(cd,'Advection'));
addpath(fullfile(cd,'LagrangianDeform'));
addpath(fullfile(cd,'miscFtns'));

% cell index
for i=1:4
home = getenv('HOME');
base_dir = fullfile(home, 'NB_Oscillations', 'meshed_data_poisson', sprintf('cell%d', i));
dat = fullfile(base_dir, 'MAT_inputs', 'DEC_inputs.mat');
load(dat); 
out_dir = fullfile(home, 'NB_Oscillations', "DEC_analysis", sprintf("cell_%d",i),"Lagrangian_Flows");
mkdir(out_dir)

cpu_num = 4;
Nt = size(timeArr,2);
dt = 1;
ct_i = 1;
ct_f = Nt;
regulFac = 0;

%% Initialize parallel pool

% Make parallel 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)                                          % if parpool is not open
    parpool('local',cpu_num);
elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
    delete(gcp)
    parpool('local',cpu_num);
end

%% Advect the grid backward 
clc; close;

tau0_Advect = 0; tauf_Advect = timeArr(ct_f)-timeArr(ct_i); 
xqBck = x{ct_f}; yqBck = y{ct_f}; zqBck = z{ct_f}; % Initial condition for backward advect 
xqFor = x{ct_i}; yqFor = y{ct_i}; zqFor = z{ct_i}; % Initial condition for forward advect

tf_Absolute = timeArr(ct_f);
tauArr_Advect = tau0_Advect:dt:tauf_Advect;
tauSave_Advect = tauArr_Advect(1:1:end); 


[xt_Advect,yt_Advect,zt_Advect] = AdvctBck(cpu_num,tf_Absolute,tauSave_Advect, ...
    tauArr_Advect,xqBck,yqBck,zqBck,timeArr,v,x,y,z,TrianT);


% % Advect the grid forward
t0_Advect_for = timeArr(ct_i); tf_Advect_for = timeArr(ct_f);
tArr_Advect_for = t0_Advect_for:dt:tf_Advect_for;
tSave_Advect_for = tArr_Advect_for(1:1:end); 


[xt_Advect_for,yt_Advect_for,zt_Advect_for] = AdvctFor(cpu_num,tSave_Advect_for, ...
    tArr_Advect_for,xqFor,yqFor,zqFor,timeArr,v,x,y,z,TrianT);



 %% backwards time
 max_time_dist = length(timeArr);

for time_dist=2:max_time_dist
x0 = squeeze(xt_Advect(1,:));
y0 = squeeze(yt_Advect(1,:));
z0 = squeeze(zt_Advect(1,:));
xf = squeeze(xt_Advect(time_dist,:));
yf = squeeze(yt_Advect(time_dist,:));
zf = squeeze(zt_Advect(time_dist,:));
tf = timeArr(end-time_dist+1);
ti = timeArr(end); % Lagrangian time interval of analysis 

Tri_pts_Uni = TrianT{end}; % The triangulation of the mesh at tf
xfData = squeeze(x{end-time_dist+1});
yfData = squeeze(y{end-time_dist+1});
zfData = squeeze(z{end-time_dist+1});clc
Trif_Data = squeeze(TrianT{end-time_dist+1});
[lambdaField,lambdaIsoField,maxDefEigVec] = lagDefCompute(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data,tf,ti,regulFac);

scals = {lambdaIsoField, lambdaField};
names = {"lambda_iso","lambda_length"};
% write_vtk_scalar_mesh(fullfile(out_dir,sprintf('lambda_output_t=%04d.vtk',time_dist)), x0, y0, z0, Tri_pts_Uni, lambdaIsoField, 'lambdaIso');
write_vtk_scalar_mesh(fullfile(out_dir,sprintf('attractors_backward_lambda_output_t=%04d.vtk',time_dist)), ...
    [x0.',y0.',z0.'], Tri_pts_Uni, scals, names)

particle_id = (1:length(xf))';
write_vtk_scalar_mesh(fullfile(out_dir,sprintf('backward_tracers_output=%04d.vtk',time_dist)), ...
    [xf.',yf.',zf.'], Tri_pts_Uni, {particle_id}, {"partcile_id"})

end
% 
% figure;
% p = trisurf(Tri_pts_Uni,x0,y0,z0,lambdaIsoField,'FaceAlpha',1,'Edgecolor','none');
% p.SpecularStrength = 0.3; 
% p.AmbientStrength = 0.4;
% hold on;  
% c = colorbar('southoutside');
% axis equal; 
% shading interp;  
% axis off;
% c.Color = 'w';
% rotate3d on

%%
max_time_dist = length(timeArr);

for time_dist=2:max_time_dist

x0 = squeeze(xt_Advect_for(1,:)); 
y0 = squeeze(yt_Advect_for(1,:)); 
z0 = squeeze(zt_Advect_for(1,:));
xf = squeeze(xt_Advect_for(time_dist,:)); 
yf = squeeze(yt_Advect_for(time_dist,:)); 
zf = squeeze(zt_Advect_for(time_dist,:));
tf = timeArr(time_dist); 
ti = timeArr(1); % Lagrangian time interval of analysis

Tri_pts_Uni = TrianT{1}; 
xfData = squeeze(x{time_dist}); 
yfData = squeeze(y{time_dist}); 
zfData = squeeze(z{time_dist}); 
Trif_Data = squeeze(TrianT{time_dist});
[lambdaField,lambdaIsoField,maxDefEigVec] = lagDefCompute(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data,tf,ti,regulFac);


scals = {lambdaIsoField, lambdaField};
names = {"lambda_iso","lambda_length"};
write_vtk_scalar_mesh(fullfile(out_dir,sprintf('repellors_forward_lambda_output_t=%04d.vtk',time_dist)), ...
    [x0.',y0.',z0.'], Tri_pts_Uni, scals, names)

 particle_id = (1:length(xf))';
write_vtk_scalar_mesh(fullfile(out_dir,sprintf('forward_tracers_output=%04d.vtk',time_dist)), ...
    [xf.',yf.',zf.'], Tri_pts_Uni, {particle_id}, {"partcile_id"})

end
% 
% figure;
% p = trisurf(Tri_pts_Uni,x0,y0,z0,lambdaField,'FaceAlpha',1,'Edgecolor','none');
% p.SpecularStrength = 0.3; 
% p.AmbientStrength = 0.4;
% hold on;  
% c = colorbar('southoutside');
% axis equal; 
% shading interp;  
% axis off;
% c.Color = 'w';
% rotate3d on
end