%— script: build_DEC_inputs.m —%

clear all; close all; clc;
% cell index
i = 3;

% build dirs
home        = getenv('HOME');
base_dir    = fullfile(home, 'NB_Oscillations', 'meshed_data_poisson', sprintf('cell%d', i));
params_file = fullfile(home, 'NB_Oscillations', 'segmented_data', sprintf('data_%d', i), sprintf('params_%d.json', i));
meshes_dir  = fullfile(base_dir, 'filtMeshes');
vel_dir     = fullfile(base_dir, 'filtMeshVels');
out_dir     = fullfile(base_dir, 'MAT_inputs');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% load start/end frames
params = jsondecode( fileread(params_file) );
st     = params.startFrame;
et     = params.endTime;

% number of time‑steps
Nt = (et - 1) - st;  
timeArr = 1:Nt;
% preallocate cell arrays
x      = cell(Nt,1);
y      = cell(Nt,1);
z      = cell(Nt,1);
TrianT = cell(Nt,1);
v      = cell(3, Nt);

% loop over each timestep k = 1..Nt
for k = 1:Nt
    t = st + (k-1);
    
    % --- read mesh ---
    mesh_fn = fullfile(meshes_dir, sprintf('filtMeshT=%04d.ply', t));
    mesh     = readSurfaceMesh(mesh_fn);  % requires a PLY reader on the path
    
    % extract vertex coords
    V = mesh.Vertices;  % (Np×3)
    x{k} = V(:,1);
    y{k} = V(:,2);
    z{k} = V(:,3);
    
    % extract faces and convert to 1‑based indexing
    TrianT{k} = mesh.Faces;
    
    % --- read velocities ---
    vel_fn = fullfile(vel_dir, sprintf('VsT=%04d.txt', t));
    VS     = dlmread(vel_fn);  % (Np×3)
    v{1,k} = VS(:,1);
    v{2,k} = VS(:,2);
    v{3,k} = VS(:,3);

end

% save all four into one .mat
out_mat = fullfile(out_dir, 'DEC_inputs.mat');
save(out_mat, 'x', 'y', 'z', 'TrianT', 'v', 'timeArr', '-v7.3');
fprintf('✓ saved DEC inputs to %s\n', out_mat);
