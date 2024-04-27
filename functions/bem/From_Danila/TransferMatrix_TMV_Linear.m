function [L,inds,b,c]=TransferMatrix_TMV_Linear(meshes, ci,co, c_in, zerolevel,omegas)
% function [L,inds,b,c]=TransferMatrix_TMV_Linear(meshes,h_mesh, ci,co, c_in, zerolevel, omegas)
% Builds linear transfer matrix for Phi in case of many boundary surfaces.
% meshes= {mesh1,mesh2,mesh3}
% Triangle meshes for boundary surfaces
%
% ci=[ci1;ci2;ci3], co=[co1;co2;co3]
% Conductivities inside and outside each boundary surface. Length of these
% vectors must be the same as number of boundary surfaces
%
% assumptions are: heart mesh is the 1st, torso - last
% h_mesh - index of the heart mesh = 1 by default
% t_mesh - index of the torso mesh = nsurf by default
%
% c_in - intracellular conductivity of the heart
%
% zerolevel
% Define surface, where the zero level of the potential is set. give 0 for "sum of
% all nodes", or I for "sum of all nodes on surface I".
% if the volume conductor is infinite, zero level does not have to be
% defined; in that case, set zerolevel to -1 or leave this field empty. 
%
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
% August 2013
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
for I=1:nsurf
    nop=nop+meshes{I}.nop;
    eind(I)=sind(I)+meshes{I}.nop-1;%index for the last node of surface I
    if I<nsurf
        sind(I+1)=eind(I)+1; %start node for surface I+1
    end
    inds{I}=[sind(I):eind(I)]';%indices for points on surface I
end


h_mesh = 1;

D=zeros(nop);
B = zeros(nop,meshes{h_mesh}.nop);

b = zeros(nsurf,1);

if nargin<6 || isempty(omegas)
    flag_calcomegas=1;
else
    flag_calcomegas=0;
end


for I=1:nsurf
    b(I)=-2*c_in/(ci(I)+co(I));%Eq. 16
    for J=1:nsurf
        c(I,J)=2*(ci(J)-co(J))/(ci(I)+co(I));%Eq. 17
        if flag_calcomegas
            if I==J && c(I,J)~=0
                omegas{I,J}=DLMatrix_Linear_SameSurface(meshes{I},1);%Eq. 33
            elseif I==J && J==h_mesh
                omegas{I,J}=DLMatrix_Linear_SameSurface(meshes{I},1);%Eq. 33
            elseif I~=J && c(I,J)~=0
                omegas{I,J}=DLMatrix_Linear_DiffSurface(meshes{I},meshes{J},1);%Eq. 33
            elseif I~=J && J==h_mesh
                omegas{I,J}=DLMatrix_Linear_DiffSurface(meshes{I},meshes{J},1);%Eq. 33
            else
                omegas{I,J}=0;
            end
        end
        D(sind(I):eind(I),sind(J):eind(J))=c(I,J)*omegas{I,J};% Eq. 26/27
        if (h_mesh==J)
            B(sind(I):eind(I),:) = b(I)*omegas{I,J};
        end
    end
end

Imatrix=eye(size(D));
Linv=Imatrix+D;%Eq. 28
clear Imatrix D

if nargin>=5 && ~isempty(zerolevel) && all(zerolevel>-1)
    %Set zero level; see ref. Fischer02
    w=zeros(1,nop);
    
    if size(zerolevel,2) == 1 && zerolevel == 0
        %sum of all surface points
        w(:) = 1/nop;
        
    elseif size(zerolevel,2) == 1 && zerolevel(1) > 0
        %'grounded' surface
        w(sind(zerolevel):eind(zerolevel)) = 1/meshes{zerolevel}.nop;
        
    elseif size(zerolevel,2) > 1 && zerolevel(1) > 0
        %multiple reference electrodes
        w(sind(zerolevel(1))+zerolevel(2:end)-1) = 1/(size(zerolevel,2)-1);
        
    else
        error('Problem with the zerolevel.');
    end
    
    defmatrix=ones(nop,1)*w;
    Linv=Linv+defmatrix;
    clear defmatrix
end

L=Linv\B;
