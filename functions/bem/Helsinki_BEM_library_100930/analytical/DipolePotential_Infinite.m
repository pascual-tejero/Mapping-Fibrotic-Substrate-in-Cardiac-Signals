function phi=DipolePotential_Infinite(points,dpos,dmoment)
% function phi=DipolePotential_Infinite(fieldpoints,dposition,dmoment)
% Electric potential for a current dipole in infinite medium of unit conductivity 
% fieldpoints = N * 3
% dposition = 1 * 3
% dmoment = 1 * 3
% Eq. 5
K=1/(4*pi);
nop=length(points);
R=points-ones(nop,1)*dpos;
absR = sqrt(sum(R.^2,2));
absRm3=1./(absR.^3);
source=ones(nop,1)*dmoment;
phi= K * sum(source.*R,2).*absRm3;    




