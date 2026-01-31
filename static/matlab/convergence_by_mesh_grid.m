common_var;  % import common_var

EPS = 1e-9;

% current simulation values
TIP_FORCE = [100; 0; 0];    % force applied to the tip of the wing (N)


% __main__
do_sub_visualization = false;
do_print_values = false;
mesh_grid = [0.1 0.07 0.05 0.04 0.03 0.02 0.015 0.01 0.0085 0.0075];
mesh_labels = ["0.1" "0.07" "0.05" "0.04" "0.03" "0.02" "0.015" "0.01" "0.0085" "0.0075"];
geometric_order = "linear";
solutions = [];

for mesh_size = mesh_grid
    fprintf('Mesh size: %.4f; Geometric order: %s\n', mesh_size, geometric_order);
    model = efoil_simulation.CreatePDEModel(DEFAULT_STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
    model = efoil_simulation.GenerateMesh(model, mesh_size, 1.0, geometric_order, do_sub_visualization);
    model = efoil_simulation.SetDefaultTestConditions(model, WING_SPAN, TIP_FORCE);
   
    fprintf('Solving linear elasticity equation\n');
    solution = solve(model);
    efoil_simulation.VizualizeSolution(model, solution, mesh_size, do_print_values, false);
    solutions = [solutions; solution];
end

VisualizeConvergenceGrid(mesh_labels, solutions, EPS);

% end __main__


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
        st = abs(cvms) + abs(fvms);
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