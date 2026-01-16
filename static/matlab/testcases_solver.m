common_var;             % import common_var

TIP_FORCE = [1000000 0 0];  % force applied to the tip of the wing (N)
root_pos = [0, 0, -1.0]; % 1m under water

% __main__

naca_shapes = ["0012" "2412" "0015" "1234"];

do_sub_visualization = false;
do_print_values = false;
do_visualization = false;

hmin_ratio = 1.0;            % not adaptive meshing
geometric_order = 'linear';  % mesh gen with linear geom order

mesh_grid = [0.02 0.015 0.01 0.009 0.007];
mesh_labels = ["0.02" "0.015" "0.01" "0.009" "0.007"];
for naca_shape = naca_shapes
    % gen shape
    [vertices, faces] = naca_4digits_gen.GetWingVerticesAndFaces(naca_shape, N_CHORD, CHORD_LEN, @GetDefaultXRoot, WING_SPAN, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, root_pos, false, 0);
    stl_filename = GetSTLFilename(naca_shape);
    naca_4digits_gen.DumpSTL(vertices, faces, stl_filename);

    % solve
    fprintf('Naca shape: %s\n', naca_shape);
    for i = 1:size(mesh_grid, 2);
        mesh_label = mesh_labels(i);
        mesh_size = mesh_grid(i);
        fprintf('Mesh size: %.4f\n', mesh_size);
        model = efoil_simulation.CreatePDEModel(stl_filename, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
        model = efoil_simulation.GenerateMesh(...
            model, mesh_size, hmin_ratio, geometric_order, do_visualization);
        model = efoil_simulation.SetDefaultTestConditions(model, WING_SPAN, TIP_FORCE);
       
        fprintf('Solving linear elasticity equation\n');
        solution = solve(model);
        results_table = GetResultsTable(solution);
        csv_filename = GetCSVFilename(naca_shape, mesh_label);
        writetable(results_table, csv_filename);
    end
    fprintf("-------------------\n\n");
end

% end __main__

function results_table = GetResultsTable(solution)
    x = solution.Mesh.Nodes(1, :).';
    y = solution.Mesh.Nodes(2, :).';
    z = solution.Mesh.Nodes(3, :).';
    sigma_x = solution.Stress.sxx;
    sigma_y = solution.Stress.syy;
    sigma_z = solution.Stress.szz;
    u = solution.Displacement.Magnitude;
    % eps ??
    columns = {'X', 'Y', 'Z', 'sigma_x', 'sigma_y', 'sigma_z', 'u'};
    results_table = table(x, y, z, sigma_x, sigma_y, sigma_z, u, 'VariableNames', columns);
end


function filename = GetSTLFilename(naca_shape)
    filename = sprintf("testcases/naca_%s/naca_%s_test_shape.stl", naca_shape, naca_shape);
end


function filename = GetCSVFilename(naca_shape, mesh_label)
    filename = sprintf("testcases/naca_%s/MATLAB_%s.csv", naca_shape, mesh_label);
end


function x = GetDefaultXRoot(~, ~, ~, ~)
    x = 0;
end