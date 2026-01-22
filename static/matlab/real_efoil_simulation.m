common_var;  % import common_var

% combined filename from ANSYS
COMBINED_WING_FILENAME = 'WWM.stl';

% mesh params
DEFAULT_MESH_ELEMENT_SIZE = 0.007;  % mesh size (m)

% current simulation values
TIP_FORCE = [100; -100; 0];  % force applied to the tip of the wing (N)

% __main__

do_sub_visualization = true;
do_print_values = false;
do_visualization = true;
hmin_ratio = 0.15;
geometric_order = 'linear';

fprintf('Default stress graph: Mesh size: %.4f\n', DEFAULT_MESH_ELEMENT_SIZE);
model = efoil_simulation.CreatePDEModel(COMBINED_WING_FILENAME, POISSONS_RATIO, YOUNGS_MODULUS, MASS_DENSITY, do_sub_visualization);
model = efoil_simulation.GenerateMesh(model, DEFAULT_MESH_ELEMENT_SIZE, hmin_ratio, geometric_order, do_sub_visualization);
model = efoil_simulation.SetRealBC(model, WING_SPAN, TIP_FORCE);

fprintf('Solving linear elasticity equation\n');
solution = solve(model);
efoil_simulation.VizualizeSolution(model, solution, DEFAULT_MESH_ELEMENT_SIZE, do_print_values, do_visualization);

% end __main__