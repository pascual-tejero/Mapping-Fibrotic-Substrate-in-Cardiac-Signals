function res=ProjectPointToSurface(rf,n,r0)
%function res=ProjectPointToSurface(rf,n,r0)
%res is a projection of rf on the surface defined by unit normal n and point r0 
%rf: one point, 1 x 3
%n:  N x 3
%r0: N x 3

rfm=ones(size(n,1),1)*rf;
tau = dots(n,(rfm-r0));
res=rfm-[tau tau tau].*n;

function dot=dots(R1,R2)
	dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
