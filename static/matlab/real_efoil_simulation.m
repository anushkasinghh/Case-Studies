common_var;  % import common_var

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.05;  % mesh size (m)

% current simulation values
TIP_FORCE = [10000 -10000 0];  % force applied to the tip of the wing (N)

% __main__

do_sub_visualization = false;
do_print_values = false;
do_visualization = true;
hmin_ratio = 0.5;
geometric_order = 'linear';

fprintf('Default stress graph: Mesh size: %.4f\n', DEFAULT_MESH_ELEMENT_SIZE);
model = efoil_simulation.CreatePDEModel(WING_WITH_MAST_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
model = efoil_simulation.GenerateMesh(model, DEFAULT_MESH_ELEMENT_SIZE, hmin_ratio, geometric_order, do_visualization);
model = efoil_simulation.SetRealBC(model, WING_SPAN, TIP_FORCE);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);
efoil_simulation.VizualizeSolution(model, solution, DEFAULT_MESH_ELEMENT_SIZE, do_print_values, do_visualization);

% end __main__