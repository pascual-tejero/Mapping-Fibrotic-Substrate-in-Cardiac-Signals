function angles=DLMatrix_Constant_FieldPoints(fp,sourcemesh,verbose)
% function angles=DLMatrix_Constant_FieldPoints(fp,sourcemesh,verbose)
% Calculates, in which solid angle elements of the sourcemesh are seen from fieldpoints.
% The result is normalized to unity.
% fieldpoints: N x 3
% verbose: 1 for some output
% Eq. 25
if nargin==2
    verbose=0;
end

nof=size(fp,1);
angles=zeros(nof,sourcemesh.noe);
if verbose
    disp('Creating DL matrix');
    dispstep=round(nof/10);
end
for I=1:nof
    if verbose & ~mod(I,dispstep)
        disp(sprintf('%2d percents done',round(100*I/nof)));
    end
    angles(I,:)=SolidAngles(sourcemesh.e, sourcemesh.p,fp(I,:))';
end
angles = angles/(4*pi);
