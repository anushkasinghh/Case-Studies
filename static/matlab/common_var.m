DEFAULT_STL_FILENAME = "default_wing.stl";
WING_WITH_MAST_FILENAME = "wing_with_mast_not_merged.stl";

% efoil params
NACA0012 = '0012';
N_CHORD = 40;                       % points along airfoil surface (per section)
CHORD_LEN = 0.15;                    % chord length (m)
WIDTH_SECTIONS = 50;                % number of spanwise sections
WING_SPAN = 0.4;                     % wing span by y-coord from 0 to WING_SPAN
REAR_WING_SPAN = WING_SPAN / 2;
REAR_WING_CHORD_LEN = CHORD_LEN / 2;
CHORD_TO_MAST_RATIO = 0.5;
MAST_WIDTH = 0.01;                   % 1cm
ROD_LEN = 0.5;                       % 0.5m
ROD_HEIGHT = 0.02;                   % 2cm

% efoil position
INIT_ROTATION_ANGLE = 0;   % initial rotation of the overall wing
ROTATE_ANGLE = 0;          % rotation around the axis Oy-plane
WATER_LEVEL = -0.5;              % z coordinate of free surface (m) 
ROOT_POS = [0, 0, WATER_LEVEL];  % root position: x - by length, y - by width, z - centerline height
MAST_ROOT_POS = [(ROOT_POS(1) + CHORD_LEN) ROOT_POS(2) (ROOT_POS(3) + ROD_HEIGHT / 2)];
REAR_WING_ROOT_POS = [(ROOT_POS(1) + ROD_LEN - REAR_WING_CHORD_LEN) ROOT_POS(2) ROOT_POS(3)];

% material carbon fiber
POISSONS_RATIO = 0.27;
YOUNGS_MODULUS = 1e11; % (Pa)
MASS_DENSITY = 1800;    % (kg/m^3)
