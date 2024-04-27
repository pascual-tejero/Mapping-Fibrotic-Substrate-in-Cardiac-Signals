function B=DipoleB_Infinite(points,dpos,dmoment,mu)
% function B=DipoleB_Infinite(points,dpos,dmoment,mu)
% Calculates magnetic flux intensity of a current dipole in an infinite
% medium
% fieldpoints = N * 3
% dpos = 1 * 3
% dmoment = 1 * 3
% mu optional; if not given, assume mu0;
% Eq. 7
% Update Feb 2010
if nargin==4 && ~isempty(mu)
    K=mu/(4*pi);
else
    K=1e-7;%mu0/(4pi)
end
nop=size(points,1);
R=points-ones(nop,1)*dpos;
absR = sqrt(sum(R.^2,2));
absRm3=1./(absR.^3)*[1 1 1];
source=ones(nop,1)*dmoment;
B = K * cross(source,R).*absRm3;    
function R=cross(R1,R2)
	R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
	R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
	R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
