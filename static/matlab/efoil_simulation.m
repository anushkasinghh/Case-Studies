classdef efoil_simulation
    methods(Static)
        function model = CreatePDEModel(filename, poissons_ratio, youngs_modulus, mass_density, do_visualization)
            model = femodel(AnalysisType='structuralStatic', Geometry=filename);
            model.MaterialProperties = materialProperties(PoissonsRatio=poissons_ratio, YoungsModulus=youngs_modulus, MassDensity=mass_density);
            if do_visualization
                figure('Name', 'Imported geometry'), pdegplot(model, 'FaceLabels', 'on'), axis equal
                title('Imported geometry')
            end
        end
        
        
        function model = GenerateDefaultMesh(model, mesh_element_size, do_visualization)
            model = generateMesh(model, 'Hmax', mesh_element_size, 'GeometricOrder', 'linear');
            if do_visualization
                figure('Name', 'Mesh'), pdemesh(model), title('FEM'); axis equal
            end
        end

        function model = GenerateMesh(model, mesh_element_size, hmin_ratio, geometric_order, do_visualization)
            model = generateMesh(...
                model, 'Hmax', mesh_element_size, 'Hmin', mesh_element_size .* hmin_ratio, 'GeometricOrder', geometric_order);
            if do_visualization
                figure('Name', 'Mesh'), pdemesh(model), title('FEM'); axis equal
            end
        end
        
        
        function VizualizeSolution(model, solution, mesh_size, do_print_values, do_visualization)
            if do_visualization
                figure('Name', sprintf('Mesh and deformed shape: mesh_size - %.4f', mesh_size));
            
                pdeplot3D(model.Geometry.Mesh)
                title('Undeformed mesh'); axis equal
            end
        
            displacement = solution.Displacement;
            d_len = sqrt(displacement.ux.^2 + displacement.uy.^2 + displacement.uz.^2);
            if do_print_values
                fprintf('Displacement: [%f, %f, %f], len = %f\n', displacement.ux, displacement.uy, displacement.uz, d_len);
            end
        
            if do_visualization
                figure('Name', sprintf('Mesh and deformed shape: mesh_size - %.4f', mesh_size));
                pdeplot3D(model.Geometry.Mesh, 'ColorMapData', d_len);
                title('Displacement'); axis equal
                colorbar
            end
        
            % Compute von Mises stress and plot
            vm_stress = solution.VonMisesStress; % sqrt(s.sxx.^2 - s.sxx.*s.syy + s.syy.^2 + 3 * s.sxy.^2);
        
            if do_visualization
                figure('Name', sprintf('Mesh and deformed shape: mesh_size - %.4f', mesh_size));
                pdeplot3D(model.Geometry.Mesh, 'ColorMapData', vm_stress);
                title('Approx. von Mises stress'); axis equal
                colorbar
            end    
        end


        function [vertex_coords, face_centers] = GetVertexCoordAndFaceCenters(model)
            num_faces = model.Geometry.NumFaces;
            vertex_coords = model.Geometry.Mesh.Nodes.';
        
            % get centroids for each face from mesh nodes
            face_centers = nan(num_faces, 3);
            for i = 1:num_faces
                face_vertex_ids = findNodes(model.Geometry.Mesh, 'region', 'Face', i);
                if isempty(face_vertex_ids)
                    continue;
                end
                pts = vertex_coords(face_vertex_ids, :);
                face_centers(i, :) = mean(pts, 1);
            end
        end


        function model = ApplyTipForce(model, wing_span, tip_force)
            [~, face_centers] = efoil_simulation.GetVertexCoordAndFaceCenters(model);

            % apply force to the tip of the wing
            y_max = max(face_centers(:, 2));
            tip_face_ids = find(abs(face_centers(:, 2) - y_max) < 5e-3 * wing_span);
            fprintf('Applying forces to the faces: %s\n', mat2str(tip_face_ids));
            model.FaceLoad(tip_face_ids) = faceLoad(SurfaceTraction=tip_force);
        end


        function model = SetDefaultTestConditions(model, wing_span, tip_force)
            [~, face_centers] = efoil_simulation.GetVertexCoordAndFaceCenters(model);
            
            % find root face (fixed relatively from the mast)
            y_min = min(face_centers(:, 2));
            root_face_ids = find(abs(face_centers(:, 2) - y_min) < 5e-3 * wing_span);
            fprintf('Applying structural constraint to faces: %s\n', mat2str(root_face_ids));
            model.FaceBC(root_face_ids) = faceBC(Constraint='fixed');
         
            % remove gravity
            model.FaceLoad = faceLoad(Gravity=[0 0 0]);
        
            % apply force to the tip of the wing
            model = efoil_simulation.ApplyTipForce(model, wing_span, tip_force);
        end


        function model = SetRealBC(model, wing_span, tip_force)
            [~, face_centers] = efoil_simulation.GetVertexCoordAndFaceCenters(model);
            
            % find mast root face (fixed relatively from the efoil board)
            z_max = min(face_centers(:, 3));
            root_face_ids = find(abs(face_centers(:, 3) - z_max) < 5e-3 * wing_span);
            fprintf('Applying structural constraint to faces: %s\n', mat2str(root_face_ids));
            model.FaceBC(root_face_ids) = faceBC(Constraint='fixed');
         
            % remove gravity
            model.FaceLoad = faceLoad(Gravity=[0 0 0]);
        
            % apply force to the tip of the wing
            model = efoil_simulation.ApplyTipForce(model, wing_span, tip_force);
        end
    end
end