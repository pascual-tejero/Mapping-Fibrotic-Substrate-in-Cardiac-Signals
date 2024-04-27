function [L,inds,b,c]=TransferMatrix_Phi_Linear(meshes,ci,co,zerolevel,b_in,omegas)
% function [L,inds,b,c]=TransferMatrix_Phi_Linear(meshes,ci,co,zerolevel, b_in, omegas)
% Builds linear transfer matrix for Phi in case of many boundary surfaces.
% meshes= {mesh1,mesh2,mesh3}
% Triangle meshes for boundary surfaces
%
% ci=[ci1;ci2;ci3], co=[co1;co2;co3]
% Conductivities inside and outside each boundary surface. Length of these
% vectors must be the same as number of boundary surfaces
%
% zerolevel
% Define surface, where the zero level of the potential is set. give 0 for "sum of
% all nodes", or I for "sum of all nodes on surface I".
% if the volume conductor is infinite, zero level does not have to be
% defined; in that case, set zerolevel to -1 or leave this field empty. 
%
% b_in
% Defines, are the constant terms of b (Eq. 16) multiplied inside the transfer matrix;
% 1 = yes, 0 = no; default value is 1.
%
% omegas (optional)
% cell array of double layer matrices as provided by DLMatrices_Linear
%
% Any number of meshes can be given.
%
% inds: indices for nodes of each surface in L
% b:    constant terms from Eq. (16)
% c:    see Eq. (17)
%
% Use of this routine demands some skill from the user. 
% June 30 2008
nsurf=length(meshes);
if length(ci)~=nsurf || length(co)~=nsurf 
    error('Number of meshes and size of the conductivity matrix do not match')
end
if nsurf==1 && ~iscell(meshes)
    temp{1}=meshes;
    meshes=temp;
end

c=zeros(nsurf);
sind(1)=1;%index for the first node of surface I
nop=0;%number of nodes in the whole model
for I=1:nsurf,
    nop=nop+meshes{I}.nop;
    eind(I)=sind(I)+meshes{I}.nop-1;%index for the last node of surface I
    if I<nsurf
        sind(I+1)=eind(I)+1; %start node for surface I+1
    end
    inds{I}=[sind(I):eind(I)]';%indices for points on surface I
end

D=zeros(nop);

if nargin<6 || isempty(omegas)
    flag_calcomegas=1;
else
    flag_calcomegas=0;
end
for I=1:nsurf,
    b(I)=2/(ci(I)+co(I));%Eq. 16
    for J=1:nsurf
        c(I,J)=2*(ci(J)-co(J))/(ci(I)+co(I));%Eq. 17
        if flag_calcomegas
            if I==J && c(I,J)~=0
                omegas{I,J}=DLMatrix_Linear_SameSurface(meshes{I},1);%Eq. 33
            elseif I~=J && c(I,J)~=0
                omegas{I,J}=DLMatrix_Linear_DiffSurface(meshes{I},meshes{J},1);%Eq. 33
            else
                omegas{I,J}=0;
            end
        end
        D(sind(I):eind(I),sind(J):eind(J))=c(I,J)*omegas{I,J};% Eq. 26/27
    end
end
Imatrix=eye(size(D));
Linv=Imatrix+D;%Eq. 28
clear Imatrix D
if nargin>=4 && ~isempty(zerolevel) && zerolevel>-1
    %Set zero level; see ref. Fischer02
    w=zeros(1,nop);
    if zerolevel==0
        w(:)=1/nop;
    else
        w(sind(zerolevel):eind(zerolevel))=1/meshes{zerolevel}.nop;
    end
    defmatrix=ones(nop,1)*w;
    Linv=Linv+defmatrix;
    clear defmatrix
end
L=inv(Linv);

%Multiply b-terms inside the matrix

if nargin<5 || isempty(b_in) 
    b_in=1;
end
if b_in==0
    return
else
    for I=1:nsurf,
        L(:,sind(I):eind(I))=L(:,sind(I):eind(I))*b(I);
    end
end