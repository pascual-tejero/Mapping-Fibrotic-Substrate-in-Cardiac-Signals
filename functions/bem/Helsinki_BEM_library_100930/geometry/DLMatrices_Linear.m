function omega=DLMatrices_Linear(meshes,verbose)
% function omega=DLMatrices_Linear(meshes,verbose)
% Builds linear double layer matrices for a set of boundary surfaces.
% meshes= {mesh_1,mesh_2,...,mesh_n}
%
% verbose (optional): 0 for no output, 1 for some, 2 for a lot.
% Eq. (33)
nsurf=length(meshes);
if nsurf==1 && ~iscell(meshes)
    temp{1}=meshes;
    meshes=temp;    
end
if nargin==1
    verbose=0;
    verbose2=0;
end
if verbose==2
    verbose2=1;
else 
    verbose2=0;
end

if verbose
    fprintf('Building %d D matrices...',nsurf*nsurf);
end
for I=1:nsurf,
    for J=1:nsurf
        if verbose2
            fprintf('\n%d of %d; ',nsurf*(I-1)+J,nsurf*nsurf);
        end
        if I==J 
            omega{I,I}=DLMatrix_Linear_SameSurface(meshes{I},verbose2);
        else
            omega{I,J}=DLMatrix_Linear_DiffSurface(meshes{I},meshes{J},verbose2);
        end
    end
end
if verbose
    fprintf('Done.\n');
end
