function [omegaB,inds]=TransferMatrix_B_Constant(coils,meshes,ci,co)
% function [omegaB,inds]=TransferMatrix_B_Constant(fieldpoints,meshes,ci,co)
% Builds omega matrix (see Eq. 6) for calculation of B, when Phi is
% calculated with constant elements.
%
% fieldpoints: points, where B is to be calculated: N x 3 matrix
%
% meshes= {mesh1,mesh2,mesh3}
% Triangle meshes for boundary surfaces
%
% ci=[ci1;ci2;ci3], co=[co1;co2;co3]
% Conductivities inside and outside each boundary surface. Length of these
% vectors must be the same as number of boundary surfaces
%
% Any number of meshes can be given.
% inds = indices for nodes of each surface in L
%
% Use of this routine demands some skill from the user. 
nsurf=length(meshes);
if length(ci)~=nsurf || length(co)~=nsurf 
    error('Number of meshes and size of the conductivity matrix do not match')
end
if nsurf==1 && ~iscell(meshes)
    temp{1}=meshes;
    meshes=temp;
end

sind(1)=1;
noe=0;
for I=1:nsurf,
    noe=noe+meshes{I}.noe;
    eind(I)=sind(I)+meshes{I}.noe-1;
    if I<nsurf
        sind(I+1)=eind(I)+1;
    end
    inds{I}=sind(I):eind(I);
end
mu0over4pi=1e-7;
noc=size(coils,1);
omegaB=zeros(3,noe,noc);
for I=1:nsurf,
    omegaBtemp=BMatrix_Constant(coils,meshes{I},1);%Eq. 6
    omegaBtemp=omegaBtemp*mu0over4pi*(ci(I)-co(I));
    omegaB(:,sind(I):eind(I),:)=omegaBtemp;
end
