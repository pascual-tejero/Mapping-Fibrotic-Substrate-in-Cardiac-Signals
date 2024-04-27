function res=Interpolate_CentroidsToNodes_Coarse(mesh,phi)
% function res=Interpolate_CentroidsToNodes_Coarse(mesh,phi)
% Interpolates values from mesh centroids to mesh nodes by in a very simple
% way: just takes mean of the values in the centroids of the
% triangles that are connected to the node.
% Use only for visualization etc, where high accuracy is not needed.
res=zeros(mesh.nop,1);
for I=1:mesh.nop,
    temp=mesh.ntri(I,1:mesh.ntri_n(I));
    res(I)=mean(phi(temp));
end