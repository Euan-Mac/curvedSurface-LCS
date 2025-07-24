function write_vtk_scalar_mesh(filename, vertices, faces, scalarFields, scalarNames)
    % Check if scalarFields and scalarNames are valid
    if length(scalarFields) ~= length(scalarNames)
        error('The number of scalar fields must match the number of scalar names.');
    end

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Unable to open file for writing: %s', filename);
    end
    
    % Write the VTK file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Surface mesh with multiple scalar fields\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');
    
    % Write the points section
    numVertices = size(vertices, 1);
    fprintf(fid, 'POINTS %d float\n', numVertices);
    for i = 1:numVertices
        fprintf(fid, '%f %f %f\n', vertices(i, :));
    end
    
    % Write the polygons section
    numFaces = size(faces, 1);
    fprintf(fid, 'POLYGONS %d %d\n', numFaces, numFaces * 4);
    for i = 1:numFaces
        fprintf(fid, '3 %d %d %d\n', faces(i, :)-1); % VTK uses 0-based indexing
    end
    
    % Write the scalar fields for vertices
    fprintf(fid, 'POINT_DATA %d\n', numVertices);

    % Loop over the scalar fields
    for j = 1:length(scalarFields)
        scalarField = scalarFields{j};
        scalarName = scalarNames{j};

        if length(scalarField) ~= numVertices
            error('Scalar field "%s" does not match the number of vertices.', scalarName);
        end
        
        % Write the scalar field section
        fprintf(fid, 'SCALARS %s float 1\n', scalarName);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i = 1:numVertices
            fprintf(fid, '%f\n', scalarField(i));
        end
    end
    
    % Close the file
    fclose(fid);
end
