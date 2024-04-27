function omega = SolidAngles_Epi(mesh)
% function omega = SolidAngles_Epi(elements,nodes)
% nodes = nop * 3 -matrix, 
% elements = noe * 3 -matrix
% This function calculates solid angles spanned by triangles according to
% Eq. 8 of  ref. vanOosterom83
% Eq. 25 in the article
% August 2013

elem = mesh.e;
nodes = mesh.p;
nop = mesh.nop;
omega = zeros(nop);
e1=elem(:,1);
e2=elem(:,2);
e3=elem(:,3);

for i=1:nop
    eyepoint = nodes(i,:);
    p=nodes-ones(nop,1)*eyepoint;
    absp=sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);

    absp1=absp(e1);
    absp2=absp(e2);
    absp3=absp(e3);

    p1=p(e1,:);
    p2=p(e2,:);
    p3=p(e3,:);

    dot12=dots(p1, p2);
    dot23=dots(p2, p3);
    dot13=dots(p1, p3);

    nom=dots(p1,cross(p2,p3));
    den=absp1.*absp2.*absp3 + dot12.*absp3 + dot13.*absp2 + dot23.*absp1;
    omegas = 2*atan2(nom,den); % was '-' sign
    
    omega(i,i) = sum(omegas)/(4*pi);
end

function dot=dots(R1,R2)
	dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
function R=cross(R1,R2)
	R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
	R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
	R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
