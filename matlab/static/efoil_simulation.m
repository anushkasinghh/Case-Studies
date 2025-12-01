common_var;  % import common_var

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.2;  % mesh size (m)

% material
DEFAULT_MATERIAL_NAME = 'aluminum';


% __main__
model = CreatePDEModel(STL_FILENAME, DEFAULT_MATERIAL_NAME);
model = GenerateMesh(model, DEFAULT_MESH_ELEMENT_SIZE);
model = SetBoundaryConditions(model, WING_SPAN);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);

VizualizeSolution(model, solution);


function model = CreatePDEModel(filename, material_name)
    model = femodel(AnalysisType='structuralStatic', Geometry=filename);
    model.MaterialProperties = materialProperties(Material=material_name);
    figure('Name', 'Imported geometry'), pdegplot(model, 'FaceLabels', 'on'), axis equal
    title('Imported geometry')
end


function model = GenerateMesh(model, mesh_element_size)
    model = generateMesh(model, 'Hmax', mesh_element_size, 'GeometricOrder', 'linear');
    figure('Name', 'Mesh'), pdemesh(model), title('FEM'); axis equal
end


function model = SetBoundaryConditions(model, wing_span)
    num_faces = model.Geometry.NumFaces;
    vertex_coords = model.Geometry.Mesh.Nodes.';

    % get centroids for each face from mesh nodes
    face_centers = nan(num_faces, 3);
    for i = 1:num_faces
        face_vertex_ids = findNodes(model.Geometry.Mesh, 'region', 'Face', i);
        if isempty(face_vertex_ids)
            continue;
        end
        pts = vertex_coords(face_vertex_ids, :);
        face_centers(i, :) = mean(pts, 1);
    end
    
    % find root face (fixed relatively from the mast)
    y_min = min(vertex_coords(:, 2));
    root_face_ids = find(abs(face_centers(:, 2) - y_min) < 5e-3 * wing_span);
    fprintf('Applying structural constraint to faces: %s\n', mat2str(root_face_ids));
    model.FaceBC(root_face_ids) = faceBC(Constraint='fixed');
 
    % calculate hydrostatic pressure on submerged faces (NOT SURE)
    % add frontal force
    % calculate boyant force
end

function VizualizeSolution(model, solution)
    figure('Name', 'Mesh and deformed shape', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

    p1 = subplot(1, 2, 1);
    pdeplot3D(model.Geometry.Mesh)
    title('Undeformed mesh'); axis equal

    d = solution.Displacement;
    scale_factor = 100; % visualization scale

    p2 = subplot(1, 2, 2);
    pdeplot3D(model.Geometry.Mesh, 'ColorMapData', sqrt(d.ux.^2 + d.uy.^2 + d.uz.^2), 'Deformation', scale_factor * [d.ux d.uy d.uz]);
    title(['Deformed (scaled \times' num2str(scale_factor) ')']); axis equal
    colorbar
end
