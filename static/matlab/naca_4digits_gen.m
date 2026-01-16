classdef naca_4digits_gen
    methods(Static)
        function coords = Build2DNaca4DigitsFoilCoordinatesXZ(naca_str, n_chord, chord_len, x_root)
            naca = convertStringsToChars(naca_str);
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
            x_full = [xu(1:n_points), xl(n_points+1:end)] * chord_len;
            y_full = [yu(1:n_points), yl(n_points+1:end)] * chord_len;
            coords = [x_full(:) + x_root, y_full(:)];
        end
        
        
        function [vertices, faces] = GetWingVerticesAndFaces(naca, n_chord, chord_len, get_cur_x_root_func, width, width_sections, init_rotatation_angle, rotate_angle, root_pos, two_sided, init_side_y_coord)
            % build vertices
            % get_cur_x_root_func - function of x coordinate of naca
            %   profile depending on y, chord_len, wing_span (has to go through (0, 0)

            R_init = [
                cos(init_rotatation_angle)  0 sin(init_rotatation_angle);
                0                           1                          0;
                -sin(init_rotatation_angle) 0 cos(init_rotatation_angle);
            ];

            y_sections = [];
            if two_sided
                y_sections = linspace(-width, width, width_sections);
            else
                y_sections = linspace(init_side_y_coord, width, width_sections);
            end
            vertices = [];
            n_points = 0; % default will be rewritten
            for j = 1:width_sections
                yj = y_sections(j);
                x_root = get_cur_x_root_func(yj, chord_len, abs(width), init_side_y_coord);
                cur_chord_len = chord_len - x_root;
                coords_2d = naca_4digits_gen.Build2DNaca4DigitsFoilCoordinatesXZ(naca, n_chord, cur_chord_len, x_root);
                n_points = size(coords_2d, 1);
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


        function [vertices, faces] = GetRecangularVerticesAndFaces(bottom_left_p, top_right_p)
            vertices = [...
                bottom_left_p(1) bottom_left_p(2) top_right_p(3);...
                bottom_left_p(1) top_right_p(2) top_right_p(3);...
                top_right_p(1) bottom_left_p(2) top_right_p(3);...
                top_right_p;...

                bottom_left_p;...
                bottom_left_p(1) top_right_p(2) bottom_left_p(3);...
                top_right_p(1) bottom_left_p(2) bottom_left_p(3);...
                top_right_p(1) top_right_p(2) bottom_left_p(3);...
            ];
            faces = [...
                1 3 2;...
                2 3 4;...
                8 2 4;...
                6 2 8;...
                8 7 6;...
                5 6 7;...
                1 5 3;...
                7 3 5;...
                1 2 5;...
                6 5 2;...
                8 4 3;...
                7 8 3;...
            ];
        end

        function [vertices, faces] = GetRecangularMastVerticesAndFaces(chord_len, chord_to_mast_ratio, mast_width, root_pos)
            x_len = chord_to_mast_ratio * chord_len;
            y_width = mast_width / 2;
            [vertices, faces] = naca_4digits_gen.GetRecangularVerticesAndFaces(...
                [root_pos(1) (root_pos(2) - y_width) root_pos(3)],...
                [(root_pos(1) + x_len) (root_pos(2) + y_width) 0]);
        end


        function [vertices, faces] = GetRecangularRodVerticesAndFaces(x_len, y_width, z_height, root_pos)
            y_width = y_width / 2;
            z_height = z_height / 2;
            [vertices, faces] = naca_4digits_gen.GetRecangularVerticesAndFaces(...
                [root_pos(1) (root_pos(2) - y_width) (root_pos(3) - z_height)],...
                [(root_pos(1) + x_len) (root_pos(2) + y_width) (root_pos(3) + z_height)]);
        end


        function DumpSTL(vertices, faces, filename)
            stlwrite(triangulation(faces, vertices), filename);
        end

        
        function VisualizeShape(filename)
            model = femodel(AnalysisType='structuralStatic', Geometry=filename);
            figure('Name', 'Imported geometry'), pdegplot(model, 'FaceLabels', 'on'), axis equal
        end
    end
end