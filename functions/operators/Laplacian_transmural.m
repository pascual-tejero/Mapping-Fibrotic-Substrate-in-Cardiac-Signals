function L = Laplacian_transmural(vtk_sur, numSubDiv)
% Relies on cotmatrix from the gptoolbox by Alec Jacobson

if nargin < 2
    numSubDiv = 1;
end

volToSur = 1:size(vtk_sur.points,1);

if numSubDiv
    vtk_sur = vtkLinearSubdivisionFilter(vtk_sur, numSubDiv);
end
vtk_vol = tetrahedralizeTriangleMesh(vtk_sur);

L_vol = cotmatrix(double(vtk_vol.points), double(vtk_vol.cells));
M_vol = full(LaplaceInterpolation(L_vol, volToSur));

L = L_vol(volToSur,:) * M_vol;
L(abs(L)<1e-7) = 0;
L = sparse(L);

end