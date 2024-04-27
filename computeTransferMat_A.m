function A = computeTransferMat_A(tissue,electrodes,bbox)
    %% Add path
    addpath('/Volumes/bordeaux/IBT/matlab/common/vtkToolbox/MATLAB')
    addpath(genpath('functions'));
    
    %% Read meshes
    tiss_mesh = PrepareTriangleMesh(tissue.points, tissue.cells);
    elec_mesh = PrepareTriangleMesh(electrodes.points, electrodes.cells);
    bbox_mesh = PrepareTriangleMesh(bbox.points, bbox.cells);
    elecIds = extractElecPointIds(electrodes);
    
    %% compute transfer matrix

    c_in = 0.14; %0.05? % intracellular conductivity of source domain. Division between 
    %ci/co of the first element

    % bulk conductivities
    ci = [0.14 1 1]; % inside of compartments [0.14 1e3]
    co = [1 1 0]; % outside of compartments [1 1]
    
    meshes = {tiss_mesh, elec_mesh,bbox_mesh}; % source surface has to come first
    electrodesSurface = 2;

    zerolevel = -1;%Volume conductor infinite (monodomain)
    [Ltmv,inds] = TransferMatrix_TMV_Linear(meshes, ci, co, c_in, zerolevel);
    
    A_allPoints = Ltmv(inds{electrodesSurface},:);

    A = NaN(size(elecIds,1), size(A_allPoints,2));
    for i = 1:size(elecIds,1)
        A(i,:) = mean(A_allPoints(elecIds(i,:),:), 1);
    end
    
    
end

