function angles=DLMatrix_Constant_DiffSurface(fieldmesh,sourcemesh,verbose)
% function angles=DLMatrix_Constant_DiffSurface(fieldmesh,sourcemesh,verbose)
% Calculates, in which solid angle element J is seen from fieldpoint I.
% The result is normalized to unity.
% verbose: 1 for some output
% Eq. 25
if nargin==2
    verbose=0;
end

angles=zeros(fieldmesh.noe,sourcemesh.noe);
if verbose
    disp('Creating DL matrix');
    dispstep=round(fieldmesh.noe/10);
end
for I=1:fieldmesh.noe
    if verbose & ~mod(I,dispstep)
        disp(sprintf('%2d percents done',round(100*I/fieldmesh.noe)));
    end
    angles(I,:)=SolidAngles(sourcemesh.e, sourcemesh.p,fieldmesh.mp(I,:))';
end
angles = angles/(4*pi);
