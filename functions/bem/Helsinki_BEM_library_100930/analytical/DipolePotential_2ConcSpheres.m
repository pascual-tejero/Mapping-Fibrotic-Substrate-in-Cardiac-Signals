function phi=DipolePotential_2ConcSpheres(momentvec,R1,R2,s1,s2,s3,points)
% function phi=DipolePotential_2ConcSpheres(momentvec,R1,R2,s1,s2,s3,points)
%Potential of a current dipole placed in the midpoint of a volume conductor 
%consisting of two concentric spheres
%dmom: dipole moment
%R1, R2: radii of the spheres
%s1: conductivity inside R1
%s2: conductivity between the spheres (R1<r<R2)
%s3: conductivity outside the spheres
%fieldpoints: N x 3 matrix
den=2*R1^3*(s1-s2)*(s2-s3)+R2^3*(s1+2*s2)*(s2+2*s3);
A=2*(R1^3*(2*s1+s2)*(s2-s3)+R2^3*(s1-s2)*(s2+2*s3))/(R1^3*den);
B=6*s1*(s2-s3)/den;
C=3*R2^3*s1*(s2+2*s3)/den;
D=9*R2^3*s1*s2/den;

noe=size(points,1);
phi=zeros(noe,1);
momval=sqrt(sum(momentvec.^2));
source=ones(noe,1)*momentvec;
absR2=points(:,1).^2+points(:,2).^2+points(:,3).^2;
absR=sqrt(absR2);
cosines=sum(source.*points,2)./absR/momval;
Q=momval/(4*pi*s1);
indexin1=find(absR<=R1);
indexin2=find(R1<absR & absR <= R2);
indexout=find(absR>R2);
A=A*Q;
B=B*Q;
C=C*Q;
D=D*Q;
phi(indexin1)=(Q./absR2(indexin1)+A.*absR(indexin1)).*cosines(indexin1);
phi(indexin2)=(B.*absR(indexin2)+C./absR2(indexin2)).*cosines(indexin2);
phi(indexout)=D.*cosines(indexout)./absR2(indexout);
