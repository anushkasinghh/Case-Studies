STL_FILENAME = 'naca_wing_temp.stl';

% efoil params
NACA0012 = '0012';
N_CHORD = 40;              % points along airfoil surface (per section)
CHORD_LEN = 0.3;           % chord length (m)
WIDTH_SECTIONS = 25;       % number of spanwise sections
WING_SPAN = 1.0;           % wing span by y-coord from 0 to WING_SPAN

% efoil position
INIT_ROTATION_ANGLE = 0;   % initial rotation of the overall wing
ROTATE_ANGLE = 0;          % rotation around the axis Oy-plane
WATER_LEVEL = -1.0;              % z coordinate of free surface (m)
ROOT_POS = [0, 0, WATER_LEVEL];  % root position: x - by length, y - by width, z - centerline height

% common physical params
WATER_DENSITY = 1000;      % kg / m^3
GRAVITY = 9.81;            % gravity m/s^2
