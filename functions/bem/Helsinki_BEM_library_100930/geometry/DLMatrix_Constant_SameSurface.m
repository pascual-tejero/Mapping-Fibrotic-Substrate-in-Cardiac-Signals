function angles=DLMatrix_Constant_SameSurface(mesh,verbose)
% function angles=DLMatrix_Constant_SameSurface(mesh,verbose)
% Calculates, in which solid angle element J is seen from midpoint of element I of the same mesh.
% The result is normalized to unity.
% verbose: 1 for some output
% Eq. 25
angles=zeros(mesh.noe);
if nargin==1
    verbose=0;
end

if verbose
    disp('Creating DL matrix');
    dispstep=round(mesh.noe/10);
end
for I=1:mesh.noe
    if verbose & ~mod(I,dispstep)
        disp(sprintf('%2d percents done',round(100*I/mesh.noe)));
    end
    
    angles(I,:)=SolidAngles(mesh.e, mesh.p,mesh.mp(I,:))';
    angles(I,I)=0;
end
angles = angles/(4*pi);
