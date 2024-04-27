function [L,inds]=TransferMatrix_Epi_Linear(meshes,ci,co,omegas)
% function [L,inds]=TransferMatrix_Epi_Linear(meshes,ci,co,omegas)
% Builds linear transfer matrix for epicardial potentials in case of many boundary surfaces.
% meshes= {mesh1,mesh2,mesh3}
% Triangle meshes for boundary surfaces
%
% ci=[ci1;ci2;ci3], co=[co1;co2;co3]
% Conductivities inside and outside each boundary surface. Length of these
% vectors must be the same as number of boundary surfaces
%
% omegas (optional)
% cell array of double layer matrices as provided by DLMatrices_Linear
%
% Any number of meshes can be given.
%
% inds: indices for nodes of each surface in L
%
% Use of this routine demands some skill from the user:
% the first mesh is for the heart, the last one - for torso
% August 28 2013
nsurf=length(meshes);
if length(ci)~=nsurf || length(co)~=nsurf 
    error('Number of meshes and size of the conductivity matrix do not match')
end
if nsurf==1 && ~iscell(meshes)
    temp{1}=meshes;
    meshes=temp;
end


h_mesh = 1;%heart mesh is the first one
t_mesh = nsurf;%torso mesh is the last one
h_size = meshes{h_mesh}.nop;
t_size = meshes{t_mesh}.nop;

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

T=zeros(nop);
T_der = zeros(nop);
B = zeros(nop,meshes{h_mesh}.nop);

if nargin<4 || isempty(omegas)
    flag_calcomegas=1;
else
    flag_calcomegas=0;
end
for I=1:nsurf
    G{I} = SLMatrix_Linear(meshes{I},meshes{h_mesh},1);
    T(sind(I):eind(I),1:h_size) = G{I};
    T_der(sind(I):eind(I),1:h_size) = 0;
    for J=1:nsurf
        if J~=t_mesh && J~=h_mesh
            c(J) = ci(J)/ci(t_mesh) - 1;
            c_plus(J) = ci(J)/ci(t_mesh) + 1;
        else
            c(J) = 1;
        end
        if flag_calcomegas
            if I==J && c(J)~=0
                omegas{I,J}=DLMatrix_Linear_SameSurface(meshes{I},1);%Eq. 33               
                solid = SolidAngles_Epi(meshes{J});
                solid = 0.5*eye(length(solid));
            elseif I~=J && c(J)~=0
                omegas{I,J}=DLMatrix_Linear_DiffSurface(meshes{I},meshes{J},1);%Eq. 33
            else
                omegas{I,J}=0;
                solid = 0;
            end
        end
        
        if I~=J && J==h_mesh
            B(sind(I):eind(I),:) = omegas{I,J};
        end
        if I~=J && J~=h_mesh && J~=t_mesh
            T(sind(I):eind(I),sind(J):eind(J)) = c(J)*omegas{I,J};
        end
        if I==J && I~=h_mesh && I~=t_mesh
            T(sind(I):eind(I),sind(I):eind(I)) = c(J)*(omegas{I,J} + solid) + eye(eind(I)-sind(I)+1);
        end
        if I==J && I==h_mesh
            B(sind(I):eind(I),:) = omegas{I,J} + solid - eye(h_size);
        end
        if I==J && I==t_mesh
            T(sind(I):eind(I),sind(I):eind(I)) = omegas{I,J} + solid;
            T_der(sind(I):eind(I),sind(I):eind(I)) = 0;
        end
        if I~=J && J==t_mesh
            T(sind(I):eind(I),sind(J):eind(J)) = omegas{I,J};
            T_der(sind(I):eind(I),sind(J):eind(J)) = 0;
        end
    end
end

L=T\B;
