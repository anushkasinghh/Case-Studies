common_var;  % import common_var

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.05;  % mesh size (m)

% material
POISSONS_RATIO = 0.27;
YOUNGS_MODULUS = 100e9; % (Pa)
MASS_DENSITY = 1800;    % (kg/m^3)

% current simulation values
TIP_FORCE = 100000; % force applied to the tip of the wing (N)
CONV_PTS_N = 1000;


% __main__
do_sub_visualization = false;
do_print_values = false;
do_main_visualization = true;
mesh_grid = [0.2 0.15 0.1 0.07 0.05 0.03 0.02 0.015 0.01 0.007 0.005];
mesh_labels = ["0.2" "0.15" "0.1" "0.07" "0.05" "0.03" "0.02" "0.015" "0.01" "0.007" "0.005"];
displacements = [];
vm_stresses = [];
for mesh_size = mesh_grid
    fprintf('Mesh size: %.4f\n', mesh_size);
    model = CreatePDEModel(STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
    model = GenerateMesh(model, mesh_size, do_sub_visualization);
    model = SetConditions(model, WING_SPAN, TIP_FORCE);
    
    fprintf('Solving linear elasticity equation\n');
    solution = solve(model);
    [cur_displacements, cur_vm_stresses] = VizualizeSolution(model, solution, WING_SPAN, mesh_size, CONV_PTS_N, CHORD_LEN, WATER_LEVEL, do_print_values, false);
    displacements = [displacements; cur_displacements(1, :)];
    vm_stresses = [vm_stresses; cur_vm_stresses(1, :)];
end

VisualizeConvergenceGrid(mesh_labels, displacements, vm_stresses);

% default stress graph
fprintf('Default stress graph: Mesh size: %.4f\n', DEFAULT_MESH_ELEMENT_SIZE);
model = CreatePDEModel(STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
model = GenerateMesh(model, DEFAULT_MESH_ELEMENT_SIZE, do_sub_visualization);
model = SetConditions(model, WING_SPAN, TIP_FORCE);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);
VizualizeSolution(model, solution, WING_SPAN, DEFAULT_MESH_ELEMENT_SIZE, CONV_PTS_N, CHORD_LEN, WATER_LEVEL, do_print_values, true);

% end __main__

function model = CreatePDEModel(filename, poissons_ratio, youngs_modulus, mass_density, do_visualization)
    model = femodel(AnalysisType='structuralStatic', Geometry=filename);
    model.MaterialProperties = materialProperties(PoissonsRatio=poissons_ratio, YoungsModulus=youngs_modulus, MassDensity=mass_density);
    if do_visualization
        figure('Name', 'Imported geometry'), pdegplot(model, 'FaceLabels', 'on'), axis equal
        title('Imported geometry')
    end
end


function model = GenerateMesh(model, mesh_element_size, do_visualization)
    model = generateMesh(model, 'Hmax', mesh_element_size, 'GeometricOrder', 'linear');
    if do_visualization
        figure('Name', 'Mesh'), pdemesh(model), title('FEM'); axis equal
    end
end


function model = SetConditions(model, wing_span, tip_force)
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
 
    % apply force to the tip of the wing
    y_max = max(vertex_coords(:, 2));
    tip_face_ids = find(abs(face_centers(:, 2) - y_max) < 5e-3 * wing_span);
    fprintf('Applying forces to the faces: %s\n', mat2str(tip_face_ids));
    model.FaceLoad(tip_face_ids) = faceLoad(SurfaceTraction=[tip_force 0 0]);

    % calculate hydrostatic pressure on submerged faces (NOT SURE)
    % add frontal force
    % calculate boyant force
end


function [picked_displacements, picked_vm_stress] = VizualizeSolution(model, solution, wing_span, mesh_size, pick_pts_n, pick_x_coord, pick_z_coord, do_print_values, do_visualization)
    if do_visualization
        figure('Name', sprintf('Mesh and deformed shape: mesh_size - %.4f', mesh_size));
    
        p1 = subplot(1, 3, 1);
        pdeplot3D(model.Geometry.Mesh)
        title('Undeformed mesh'); axis equal
    end

    d = solution.Displacement;
    d_len = sqrt(d.ux.^2 + d.uy.^2 + d.uz.^2);
    scale_factor = 100; % visualization scale
    if do_print_values
        fprintf('Displacement: [%f, %f, %f], len = %f\n', d.ux, d.uy, d.uz, d_len);
    end

    if do_visualization
        p2 = subplot(1, 3, 2);
        pdeplot3D(model.Geometry.Mesh, 'ColorMapData', d_len);
        title('Displacement'); axis equal
        colorbar
    end

    % Compute von Mises stress and plot
    % s = solution.Stress;
    vm_stress = solution.VonMisesStress; % sqrt(s.sxx.^2 - s.sxx.*s.syy + s.syy.^2 + 3 * s.sxy.^2);
    % if do_print_values
    %     fprintf('Stress: [%f, %f, %f], len = %f\n', s.sxx, s.sxy, s.syy, stress);
    % end

    if do_visualization
        p3 = subplot(1, 3, 3);
        pdeplot3D(model.Geometry.Mesh, 'ColorMapData', vm_stress);
        title('Approx. von Mises stress'); axis equal
        colorbar
    end

    picked_displacements = zeros(pick_pts_n);
    picked_vm_stress = zeros(pick_pts_n);
    y_vals = linspace(0, wing_span, pick_pts_n);
    for i = 1:pick_pts_n
        v = [pick_x_coord; y_vals(i); pick_z_coord];
        node_id = findNodes(solution.Mesh, 'nearest', v);
        picked_displacements(1, i) = d_len(node_id, 1);
        picked_vm_stress(1, i) = vm_stress(node_id, 1);
    end
end


function VisualizeConvergenceGrid(x_labels, displacements, vm_stresses)
    x_labels_t = x_labels.';
    n = size(displacements, 1);
    displacement_diffs = [];
    vm_stress_diffs = [];
    for i = 1:n-1
        displacement_diffs = [displacement_diffs; norm(displacements(i) - displacements(n), 2)];
        vm_stress_diffs = [vm_stress_diffs; norm(vm_stresses(i) - vm_stresses(n), 2)];
    end

    figure('Name', sprintf('Mesh convergence chart'));

    p1 = subplot(1, 2, 1);
    plot(1:n-1, displacement_diffs.', '-o');
    xticks(1:n-1);
    xticklabels(x_labels_t(1:n-1));
    title('Convergence by displacement');
    ylabel('L2-norm diff by displacement');

    p2 = subplot(1, 2, 2);
    plot(1:n-1, vm_stress_diffs.', '-o');
    xticks(1:n-1);
    xticklabels(x_labels_t(1:n-1));
    title('Convergence by von Mises stress');
    ylabel('L2-norm diff by von Mises stress');
end