common_var;  % import common_var

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.05;  % mesh size (m)

% current simulation values
TIP_FORCE = [100; 0; 0];  % force applied to the tip of the wing (N)
wing_span = 1.0;          % 1m
chord_len = 0.3;          % 0.3m


% __main__

do_sub_visualization = false;
do_print_values = false;
do_visualization = true;

[vertices, faces] = naca_4digits_gen.GetWingVerticesAndFaces(NACA0012, N_CHORD, chord_len, @GetDefaultXRoot, wing_span, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, ROOT_POS, false, 0);
naca_4digits_gen.DumpSTL(vertices, faces, DEFAULT_STL_FILENAME);
naca_4digits_gen.VisualizeShape(DEFAULT_STL_FILENAME);

fprintf('Default stress graph: Mesh size: %.4f\n', DEFAULT_MESH_ELEMENT_SIZE);
model = efoil_simulation.CreatePDEModel(DEFAULT_STL_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
model = efoil_simulation.GenerateDefaultMesh(model, DEFAULT_MESH_ELEMENT_SIZE, do_sub_visualization);
model = efoil_simulation.SetDefaultTestConditions(model, wing_span, TIP_FORCE);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);
efoil_simulation.VizualizeSolution(model, solution, DEFAULT_MESH_ELEMENT_SIZE, do_print_values, do_visualization);

% end __main__


function x = GetDefaultXRoot(~, ~, ~, ~)
    x = 0;
end