common_var;  % import common_var


% __main__
coords_2d = Build2DNaca4DigitsFoilCoordinatesXZ(NACA0012, N_CHORD, CHORD_LEN);
[vertices, faces] = ExtrudeAlongWidthAndCreateSurfaceTriangles(WING_SPAN, WIDTH_SECTIONS, coords_2d, INIT_ROTATION_ANGLE, ROTATE_ANGLE, ROOT_POS);
DumpSTL(vertices, faces, STL_FILENAME);


function coords = Build2DNaca4DigitsFoilCoordinatesXZ(naca, n_chord, chord_len)
    % parse values
    m = str2double(naca(1)) / 100;         % camber as percentage of the chord
    p = str2double(naca(2)) / 10;          % distance of maximum camber from the foil leading edge in tenths of the chord
    t = str2double(naca(3:4)) / 100;       % maximum thickness of the airfoil as percent of the chord
    
    % x distribution (cosine spacing for smooth leading edge)
    n_points = round(n_chord / 2);
    a = linspace(0, pi, n_points);
    x_direct = (1 - cos(a)) / 2;
    x_reversed = fliplr(x_direct(2:end-1));
    x = [x_direct, x_reversed];            % full circle

    % yt is the half thickness at a given value of x (centerline to surface),
    % x is the position along the chord from 0 to 1.00 (0 to 100%)
    yt = 5 * t * (0.2969 * sqrt(x) - 0.1260 * (x) - 0.3516 * (x).^2 + 0.2843 * (x).^3 - 0.1015 * (x).^4);
    
    % build Equation for a cambered 4-digit NACA airfoil
    yc = zeros(size(x));
    yc_diff = zeros(size(x));
    for i = 1:length(x)
        xi = x(i);
        if 0 <= xi && xi < p
            yc(i) = (m / p^2) * (2 * p * xi - xi.^2);
            yc_diff(i) = (m / p^2) * (2 * p - 2 * xi);
        else
            yc(i) = (m / (1 - p).^2) * ((1 - 2 * p) + 2 * p * xi - xi.^2);
            yc_diff(i) = (m / (1 - p).^2) * (2 * p - 2 * xi);
        end
    end

    theta = atan(yc_diff);
    
    % build the upper and lower airfoil surface
    xu = x - yt .* sin(theta);
    xl = x + yt .* sin(theta);
    yu = yc + yt .* cos(theta);
    yl = yc - yt .* cos(theta);
    
    % combine into single surface loop (upper forward, lower back)
    x_full = [xu(1:n_points), xl(n_points+1:end)];
    y_full = [yu(1:n_points), yl(n_points+1:end)];
    coords = [x_full(:), y_full(:)] * chord_len;
end


function [vertices, faces] = ExtrudeAlongWidthAndCreateSurfaceTriangles(width, width_sections, coords_2d, init_rotatation_angle, rotate_angle, root_pos)
    % build vertices
    R_init = [
        cos(init_rotatation_angle)  0 sin(init_rotatation_angle);
        0                           1                          0;
        -sin(init_rotatation_angle) 0 cos(init_rotatation_angle);
    ];
    n_points = size(coords_2d, 1);
    y_sections = linspace(0, width, width_sections);
    vertices = [];
    for j = 1:width_sections
        yj = y_sections(j);
        phi = rotate_angle * (yj / width);
        R = [cos(phi) 0 sin(phi);
             0        1        0;
            -sin(phi) 0 cos(phi)];
        section_points = zeros(n_points, 3);
        for i = 1:n_points
            v = [coords_2d(i, 1); 0; coords_2d(i, 2)];
            v = R_init * R * v;   % rotate along y-axis
            section_points(i, :) = [v(1) + root_pos(1), yj + root_pos(2), v(3) + root_pos(3)];
        end
        vertices = [vertices; section_points];
    end
    
    % build faces
    faces = [];
    for j = 1:(width_sections - 1)
        offset1 = (j - 1) * n_points;
        offset2 = j * n_points;
        for i = 1:n_points
            i1 = i;
            i2 = mod(i, n_points) + 1;

            p1 = offset1 + i1;
            p2 = offset1 + i2;
            p3 = offset2 + i1;
            p4 = offset2 + i2;

            faces = [faces; p1 p3 p4];
            faces = [faces; p1 p4 p2];
        end
    end
    
    % combine
    root_idx = 1:n_points;
    tip_idx = (width_sections - 1) * n_points + (1:n_points);
    
    centroid_root = mean(vertices(root_idx,:), 1);
    vertices = [vertices; centroid_root];
    cr = size(vertices, 1);
    for i = 1:n_points
        i1 = root_idx(i);
        i2 = root_idx(mod(i, n_points) + 1);
        faces = [faces; i1 i2 cr];
    end
    
    centroid_tip = mean(vertices(tip_idx, :), 1);
    vertices = [vertices; centroid_tip];
    ct = size(vertices, 1);
    for i = 1:n_points
        i1 = tip_idx(i);
        i2 = tip_idx(mod(i, n_points) + 1);
        faces = [faces; i2 i1 ct];
    end
end


function WriteStlFile(filename, vertices, faces)
    f = fopen(filename, 'w');
    fprintf(f,'solid wing\n');
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        n = cross(v2 - v1, v3 - v1);
        if norm(n) < eps
            n = [0 0 0];
        else
            n = n / norm(n);
        end
        fprintf(f, ' facet normal %e %e %e\n', n);
        fprintf(f, '  outer loop\n');
        fprintf(f, '   vertex %e %e %e\n', v1);
        fprintf(f, '   vertex %e %e %e\n', v2);
        fprintf(f, '   vertex %e %e %e\n', v3);
        fprintf(f, '  endloop\n');
        fprintf(f, ' endfacet\n');
    end
    fprintf(f,'endsolid wing\n');
    fclose(f);
end


function DumpSTL(vertices, faces, filename)
    WriteStlFile(filename, vertices, faces);
    fprintf('STL file dumped: %s\n', filename);
end
