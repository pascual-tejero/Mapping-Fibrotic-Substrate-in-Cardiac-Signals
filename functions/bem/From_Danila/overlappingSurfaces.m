function [isOverlapping,ind1,ind2] = overlappingSurfaces(mesh1, mesh2)
% for BEM calculations of blood and heart potentials, but might be used for
% other purposes as well

nop1 = mesh1.nop;
nop2 = mesh2.nop;

eps_tol = 1e-4;
isOverlapping = 0;

ind1 = [];
ind2 = [];

for i=1:nop1
    p = mesh1.p(i,:);
    dif = sqrt(sum((repmat(p,nop2,1) - mesh2.p).^2,2));
    intersect_point = find(dif<eps_tol);
    if ~isempty(intersect_point)
        ind1 = [ind1; i];
        ind2 = [ind2; intersect_point];
        isOverlapping = 1;
    end
end

end