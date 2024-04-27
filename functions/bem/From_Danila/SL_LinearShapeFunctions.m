function shapes=SL_LinearShapeFunctions(fieldpoint,sourcemesh)
% function shapes=SL_LinearShapeFunctions_Test(fieldpoint,elements,nodes,normals,areas)
%
%  from paper (Rajesh Srivastava):
% "Efficient evaluation of integrals in three- dimensional 
% boundary element method using linear shape functions over plane triangular elements"
% Update 10 September 2013
eps = 1e-12;

non=sourcemesh.nop;
not=sourcemesh.noe;
triangles=sourcemesh.e;
nodes=sourcemesh.p;
normals=sourcemesh.un;
areas=sourcemesh.a;

shapes=zeros(not,3);

field_vec = ones(non,1)*fieldpoint;
dif = nodes - field_vec;
ind_i = find((dif(:,1).^2 + dif(:,2).^2 + dif(:,3).^2)<=eps);


sing_tri = [];
if ~isempty(ind_i)
    sing_tri = sourcemesh.ntri(ind_i,1:sourcemesh.ntri_n(ind_i));
    for i=1:length(sing_tri)
       I = sourcemesh.ntri_s(ind_i,i);
       switch I
        case 1
            J = 2;
            K = 3;
        case 2
            J = 3;
            K = 1;
        case 3
            J = 1;
            K = 2;
       end
       %sing_tri(i)
       pi = nodes(triangles(sing_tri(i),I),:);
       pj = nodes(triangles(sing_tri(i),J),:);
       pk = nodes(triangles(sing_tri(i),K),:);
       Rij = norm(pi-pj);
       Rki = norm(pk-pi);
       Rjk = norm(pj-pk);
       Ae = areas(sing_tri(i));
       cos_theta = dot(pi-pk,pj-pk)/(Rki*Rjk);
       
       Gi = Ae*log((Rij + Rki + Rjk)/(Rij + Rki - Rjk))/Rjk;
       Gj = Ae*(Rij - Rki + Rki*cos_theta*log((Rij+Rki+Rjk)/(Rij+Rki-Rjk)))/Rjk^2;
       Gk = Gi - Gj;
       
       shapes(sing_tri(i),I) = shapes(sing_tri(i),I) + Gi;
       shapes(sing_tri(i),J) = shapes(sing_tri(i),J) + Gj;
       shapes(sing_tri(i),K) = shapes(sing_tri(i),K) + Gk;
       
    end
end


tri = setdiff([1:sourcemesh.noe]',sing_tri);

p1 = nodes(triangles(tri,1),:);
p2 = nodes(triangles(tri,2),:);
p3 = nodes(triangles(tri,3),:);

p31 = p1 - p3;
R31 = sqrt(p31(:,1).*p31(:,1) + p31(:,2).*p31(:,2) + p31(:,3).*p31(:,3));

ar = 2*areas(tri)./R31;
ar2 = 2*areas(tri)./(R31.*R31);

% determine Gausian points: 7(11) scheme
nq = 11;
noe = length(tri);

g(1) = 0;% g(1) = 0;
g(2) = -0.2695431559523450;% g(2) = -0.4058451513773972;
g(3) = -g(2);% g(3) = -g(2);
g(4) = -0.5190961292068118;% g(4) = -0.7415311855993945;
g(5) = -g(4);% g(5) = -g(4);
g(6) = -0.7301520055740494;% g(6) = -0.9491079123427585;
g(7) = -g(6);% g(7) = -g(6);
g(8) = -0.8870625997680953;
g(9) = -g(8);
g(10) = -0.9782286581460570;
g(11) = -g(10);
% 
w(1) = 0.2729250867779006;% w(1) = 0.4179591836734694;
w(2) = 0.2628045445102467;% w(2) = 0.3818300505051189;
w(3) = w(2);% w(3) = w(2);
w(4) = 0.2331937645919905;% w(4) = 0.2797053914892766;
w(5) = w(4);% w(5) = w(4);
w(6) = 0.1862902109277343;% w(6) = 0.1294849661688697;
w(7) = w(6);% w(7) = w(6);
w(8) = 0.1255803694649046;
w(9) = w(8);
w(10) = 0.0556685671161737;
w(11) = w(10);

w = 0.5*w;
g = 0.5*(g+1);

g = g';
w = w';

P_vec = p2 - p1;
Q_vec = p2 - p3;

%tmp = zeros(nq,noe,3);
% P_p = zeros(nq,noe,3);
% Q_p = zeros(nq,noe,3);

% Rp = zeros(noe,nq);
% Rq = zeros(noe,nq);
% Rr = zeros(noe,nq);

G1 = zeros(noe,1);
G2 = zeros(noe,1);
G3 = zeros(noe,1);


field_vec = ones(noe,1)*fieldpoint;

for i=1:nq
    P_p = p1 + P_vec*g(i);
    Q_p = p3 + Q_vec*g(i);
    
    op = P_p - field_vec;
    Rp = sqrt(op(:,1).*op(:,1) + op(:,2).*op(:,2) + op(:,3).*op(:,3));
    
    oq = Q_p - field_vec;
    Rq = sqrt(oq(:,1).*oq(:,1) + oq(:,2).*oq(:,2) + oq(:,3).*oq(:,3));
    
    qp = P_p - Q_p;
    Rr = sqrt(qp(:,1).*qp(:,1) + qp(:,2).*qp(:,2) + qp(:,3).*qp(:,3));
    
    cos_theta = dots(qp,-oq)./(Rr.*Rq);
    
    num = Rp + Rq + Rr;
    den = Rp + Rq - Rr;
    l = log(num./den);
    int1 = Rp - Rq + Rq.*cos_theta.*l;
    int2 = g(i)*l;
    int3 = l;
    
    G1 = G1 + w(i)*int1;
    G2 = G2 + w(i)*int2;
    G3 = G3 + w(i)*int3;
end

G1 = ar2.*G1;
G2 = ar.*G2;
G3 = ar.*G3 - G1 - G2;

G1 = -G1;
G2 = -G2;
G3 = -G3;

shapes(tri,1) = -(shapes(tri,1) + G1);
shapes(tri,2) = -(shapes(tri,2) + G2);
shapes(tri,3) = -(shapes(tri,3) + G3);



function R=cross(R1,R2)
R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);

function dot=dots(R1,R2)
dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
