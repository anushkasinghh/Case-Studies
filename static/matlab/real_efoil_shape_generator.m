common_var; % import common_var

LEFT_WING_FILENAME = "tmp/left_wing.stl";
RIGHT_WING_FILENAME = "tmp/right_wing.stl";
REAR_LEFT_WING_FILENAME = "tmp/rear_left_wing.stl";
REAR_RIGHT_WING_FILENAME = "tmp/rear_right_wing.stl";
REAR_WING_FILENAME = "tmp/rear_right_wing.stl";
ROD_FILENAME = "tmp/rod.stl";
MAST_FILENAME = "tmp/mast.stl";


% __main__

% left wing
[left_wing_vertices, left_wing_faces] = naca_4digits_gen.GetWingVerticesAndFaces(NACA0012, N_CHORD, CHORD_LEN, @GetParabolicRoot, -WING_SPAN, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, ROOT_POS, false, -MAST_WIDTH / 2.0);
naca_4digits_gen.DumpSTL(left_wing_vertices, left_wing_faces, LEFT_WING_FILENAME);
naca_4digits_gen.VisualizeShape(LEFT_WING_FILENAME);

% right wing
[right_wing_vertices, right_wing_faces] = naca_4digits_gen.GetWingVerticesAndFaces(NACA0012, N_CHORD, CHORD_LEN, @GetParabolicRoot, WING_SPAN, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, ROOT_POS, false, MAST_WIDTH / 2.0);
naca_4digits_gen.DumpSTL(right_wing_vertices, right_wing_faces, RIGHT_WING_FILENAME);
naca_4digits_gen.VisualizeShape(RIGHT_WING_FILENAME);

% rear left wing
[rear_left_wing_vertices, rear_left_wing_faces] = naca_4digits_gen.GetWingVerticesAndFaces(NACA0012, N_CHORD, REAR_WING_CHORD_LEN, @GetParabolicRoot, -REAR_WING_SPAN, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, REAR_WING_ROOT_POS, false, -MAST_WIDTH / 2.0);
naca_4digits_gen.DumpSTL(rear_left_wing_vertices, rear_left_wing_faces, REAR_LEFT_WING_FILENAME);
naca_4digits_gen.VisualizeShape(REAR_LEFT_WING_FILENAME);

% rear right wing
[rear_right_wing_vertices, rear_right_wing_faces] = naca_4digits_gen.GetWingVerticesAndFaces(NACA0012, N_CHORD, REAR_WING_CHORD_LEN, @GetParabolicRoot, REAR_WING_SPAN, WIDTH_SECTIONS, INIT_ROTATION_ANGLE, ROTATE_ANGLE, REAR_WING_ROOT_POS, false, MAST_WIDTH / 2.0);
naca_4digits_gen.DumpSTL(rear_right_wing_vertices, rear_right_wing_faces, REAR_RIGHT_WING_FILENAME);
naca_4digits_gen.VisualizeShape(REAR_RIGHT_WING_FILENAME);

% rod
[rod_vertices, rod_faces] = naca_4digits_gen.GetRecangularRodVerticesAndFaces(ROD_LEN, MAST_WIDTH, ROD_HEIGHT, ROOT_POS);
naca_4digits_gen.DumpSTL(rod_vertices, rod_faces, ROD_FILENAME);
naca_4digits_gen.VisualizeShape(ROD_FILENAME);

% mast
[mast_vertices, mast_faces] = naca_4digits_gen.GetRecangularMastVerticesAndFaces(CHORD_LEN, CHORD_TO_MAST_RATIO, MAST_WIDTH, MAST_ROOT_POS);
naca_4digits_gen.DumpSTL(mast_vertices, mast_faces, MAST_FILENAME);
naca_4digits_gen.VisualizeShape(MAST_FILENAME);

% merge
CombineAndDumpSTL([...
    LEFT_WING_FILENAME,...
    RIGHT_WING_FILENAME,...
    ROD_FILENAME,...
    MAST_FILENAME,...
    REAR_LEFT_WING_FILENAME,...
    REAR_WING_FILENAME,...
], WING_WITH_MAST_FILENAME);
naca_4digits_gen.VisualizeShape(WING_WITH_MAST_FILENAME);

% end __main__

function x = GetParabolicRoot(y, chord_len, wing_span, ~)
    if y >= 0
        x = GetParabolicRootRightSide(y, chord_len, wing_span);
    else
        x = GetParabolicRootLeftSide(y, chord_len, wing_span);
    end
end


function x = GetParabolicRootRightSide(y, chord_len, wing_span)
    x = chord_len * (1 - sqrt(1 - y / wing_span));
end

function x = GetParabolicRootLeftSide(y, chord_len, wing_span)
    x = chord_len * (1 - sqrt(1 + y / wing_span));
end

function CombineAndDumpSTL(filenames, combined_filename)
    % taken from page https://www.mathworks.com/matlabcentral/answers/2131056-combine-two-stl-into-one
    Tc{1} = stlread(filenames(1));
    for i = 2:size(filenames, 2)
        Tc{i} = stlread(filenames(i));
    end
    [faces vertices] = tricat(Tc); % attached
    stlwrite(triangulation(faces, vertices), combined_filename);
end
