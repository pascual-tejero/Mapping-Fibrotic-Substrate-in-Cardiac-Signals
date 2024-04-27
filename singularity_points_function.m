function [singularity_points,phase] = singularity_points_function(values, mesh)

%     singularity_points = [];
%     addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/EGMDataProcessing/FilteringAndProcessing')
%     values = filterSignalsClinically(values,333,0.1,166,2);
%     sigma = 1.5;
%     values = gaussFiltfilt(values, sigma); % Gauss filter
    
    values_n = normalize(values','zscore','robust');
    H = hilbert(values_n);
    phase = atan2(imag(H),real(H))';
    
    singularity_points = zeros(size(values,1),size(values,2));

    cells = mesh.cells;   
%     points = mesh.points;

%     TR = triangulation(double(mesh.cells),double(mesh.points));
    
    for j = 1:size(singularity_points,2) 
%         for k = 1:size(phase,1)
        for k = 1:size(cells,1)
            vector_cells = cells(k,:);
            sum_phase = sum(phase(vector_cells,j));
            
%             V = vertexAttachments(TR,k);
%             cells(V{:},:)
%             pts = unique(cells(V{:},:));
%             pts(pts==k) = [];
%             sum_phase = sum(phase(pts,j));

            if ((sum_phase >= (2*pi)-(2*pi)*0.001) && (sum_phase <= (2*pi)+(2*pi)*0.001))
                singularity_points(vector_cells,j) =  1;
            end
            
            if ((sum_phase >= (-2*pi)+(-2*pi)*0.001) && (sum_phase <= (-2*pi)-(-2*pi)*0.001))
                singularity_points(vector_cells,j) =  1;
            end

        end
    end
%     comb = vtkAppendPolyData({mesh, electrodes});
%     qtrip(comb.points, comb.cells, phase', 1041, [0 90 90], 0.5, [min(min(phase)) max(max(phase))], 10);
%     qtrip(comb.points, comb.cells, singularity_points', 1041, [0 90 90], 1, [min(min(singularity_points)) max(max(singularity_points))], 10);

    
    %figure; plot(phase(singularity_points(cont-1,4),:));
end