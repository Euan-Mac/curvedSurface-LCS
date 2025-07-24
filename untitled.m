function write_vtk_points(file_path, x, y, z, scalars, scalar_names)
    % Write VTK PolyData file with points (no cells), and optional scalars

    fid = fopen(file_path, 'w');

    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Point data\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');

    n_points = length(x);
    fprintf(fid, 'POINTS %d float\n', n_points);
    fprintf(fid, '%f %f %f\n', [x(:) y(:) z(:)]');

    fprintf(fid, 'VERTICES %d %d\n', n_points, 2*n_points);
    for i = 0:n_points-1
        fprintf(fid, '1 %d\n', i);
    end

    if ~isempty(scalars)
        fprintf(fid, 'POINT_DATA %d\n', n_points);
        for i = 1:length(scalars)
            fprintf(fid, 'SCALARS %s float 1\n', scalar_names{i});
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', scalars{i});
        end
    end

    fclose(fid);
end