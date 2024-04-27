function Bvol=B_VolumeCurrents(phi,omegaB)
% function Bvol=B_VolumeCurrents(phi,omegaB,inds)
% Calculates B due to volume currents
% phi is calculated with matrix L resulting from function
% "TransferMatrix_Phi_Constant" or "TransferMatrix_Phi_Linear" 
% omegaB is formed with function "TransferMatrix_B_Constant" or "TransferMatrix_B_Linear"
nofp=size(omegaB,3);
Bvol=zeros(nofp,3);
for I=1:nofp,
    Bvol(I,:)=(omegaB(:,:,I)*phi)';%Integral of Eq. 6
end