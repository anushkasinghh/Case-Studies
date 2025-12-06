common_var;  % import common_var

EPS = 1e-9;

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.05;  % mesh size (m)

% material carbon fiber
POISSONS_RATIO = 0.27;
YOUNGS_MODULUS = 100e9; % (Pa)
MASS_DENSITY = 1800;    % (kg/m^3)

% current simulation values
TIP_FORCE = 1000; % force applied to the tip of the wing (N)


% __main__
do_sub_visualization = false;
do_print_values = false;
do_main_visualization = true;
mesh_grid = [0.1 0.07 0.05 0.04 0.03 0.02 0.015 0.01 0.0085 0.007];
mesh_labels = ["0.1" "0.07" "0.05" "0.04" "0.03" "0.02" "0.015" "0.01" "0.0085" "0.0075"];
solutions = [];
for mesh_size = mesh_grid
    fprintf('Mesh size: %.4f\n', mesh_size);
    model = CreatePDEModel(STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
    model = GenerateMesh(model, mesh_size, do_sub_visualization);
    model = SetConditions(model, WING_SPAN, TIP_FORCE);
    
    fprintf('Solving linear elasticity equation\n');
    solution = solve(model);
    VizualizeSolution(model, solution, mesh_size, do_print_values, false);
    solutions = [solutions; solution];
end

VisualizeConvergenceGrid(mesh_labels, solutions, EPS);

% default stress graph
fprintf('Default stress graph: Mesh size: %.4f\n', DEFAULT_MESH_ELEMENT_SIZE);
model = CreatePDEModel(STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
model = GenerateMesh(model, DEFAULT_MESH_ELEMENT_SIZE, do_sub_visualization);
model = SetConditions(model, WING_SPAN, TIP_FORCE);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);
VizualizeSolution(model, solution, DEFAULT_MESH_ELEMENT_SIZE, do_print_values, true);

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
 
    % remove gravity
    model.FaceLoad = faceLoad(Gravity=[0 0 0]);

    % apply force to the tip of the wing
    y_max = max(vertex_coords(:, 2));
    tip_face_ids = find(abs(face_centers(:, 2) - y_max) < 5e-3 * wing_span);
    fprintf('Applying forces to the faces: %s\n', mat2str(tip_face_ids));
    model.FaceLoad(tip_face_ids) = faceLoad(SurfaceTraction=[tip_force 0 0]);

    % calculate hydrostatic pressure on submerged faces (NOT SURE)
    % add frontal force
    % calculate boyant force
end


function VizualizeSolution(model, solution, mesh_size, do_print_values, do_visualization)
    if do_visualization
        figure('Name', sprintf('Mesh and deformed shape: mesh_size - %.4f', mesh_size));
    
        p1 = subplot(1, 3, 1);
        pdeplot3D(model.Geometry.Mesh)
        title('Undeformed mesh'); axis equal
    end

    displacement = solution.Displacement;
    d_len = sqrt(displacement.ux.^2 + displacement.uy.^2 + displacement.uz.^2);
    if do_print_values
        fprintf('Displacement: [%f, %f, %f], len = %f\n', displacement.ux, displacement.uy, displacement.uz, d_len);
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
end


function [by_displacement, by_vm_stress] = GetDiffValues(sc, sf, eps)  % coarse, fine grids
    cmesh_nodes = sc.Mesh.Nodes.';
    n = size(cmesh_nodes, 1);
    by_displacement = 0;
    by_vm_stress = 0;
    for cnode_id = 1:n
        i = findNodes(sf.Mesh, 'nearest', cmesh_nodes(cnode_id, :).');

        % add displacement
        cd = sc.Displacement;
        fd = sf.Displacement;
        mg = (cd.Magnitude(cnode_id) + fd.Magnitude(i)) / 2;
        if abs(mg) > eps
            sqdx = (cd.ux(cnode_id) - fd.ux(i)).^2;
            sqdy = (cd.uy(cnode_id) - fd.uy(i)).^2;
            sqdz = (cd.uz(cnode_id) - fd.uz(i)).^2;
            by_displacement = by_displacement + sqrt(sqdx + sqdy + sqdz) / mg;
        end

        % add stress
        cvms = sc.VonMisesStress(cnode_id);
        fvms = sf.VonMisesStress(i);
        st = (abs(cvms) + abs(fvms)) / 2;
        if abs(st) > eps
            by_vm_stress = by_vm_stress + abs(cvms - fvms) / st;
        end
    end

    by_displacement = by_displacement / n;
    by_vm_stress = by_vm_stress / n;
end


function VisualizeConvergenceGrid(x_labels, solutions, eps)
    x_labels_t = x_labels.';
    n = size(solutions, 1);
    displacement_diffs = [];
    vm_stress_diffs = [];
    for i = 1:n-1
        [d, s] = GetDiffValues(solutions(i), solutions(i+1), eps);
        displacement_diffs = [displacement_diffs; d];
        vm_stress_diffs = [vm_stress_diffs; s];
    end
    displacement_diffs = [displacement_diffs; displacement_diffs(n-1)];
    vm_stress_diffs = [vm_stress_diffs; vm_stress_diffs(n-1)];

    fprintf('Displacements: %s\n', mat2str(displacement_diffs));
    fprintf('VM Stress: %s\n', mat2str(vm_stress_diffs));

    figure('Name', sprintf('Mesh convergence chart'));

    p1 = subplot(1, 2, 1);
    plot(1:n, displacement_diffs.', '-o');
    xticks(1:n);
    xticklabels(x_labels_t(1:n));
    title('Convergence by displacement');
    ylabel('Average L1-norm difference by displacement');

    p2 = subplot(1, 2, 2);
    plot(1:n, vm_stress_diffs.', '-o');
    xticks(1:n);
    xticklabels(x_labels_t(1:n));
    title('Convergence by von Mises stress');
    ylabel('Average L1-norm difference by von Mises stress');
end