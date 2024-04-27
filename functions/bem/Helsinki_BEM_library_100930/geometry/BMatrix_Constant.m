function res=BMatrix_Constant(fieldpoints,sourcemesh,verbose)
% function res=BMatrix_Constant(fieldpoints,sourcemesh,verbose)
% fieldpoints: N x 3
% verbose (optional): 0 for no output, 1 for some, 2 for a lot of.
% Eq. 6
if nargin==2
    verbose=0;
    verbose2=0;
end
if verbose==2
    verbose2=1;
else 
    verbose2=0;
end

nop=size(fieldpoints,1);
res=zeros(3,sourcemesh.noe,nop);
if verbose
    fprintf('Creating basis functions for calculation of B...');
    dispstep=round(nop/10);
end
for I=1:nop,
    if verbose2 && ~mod(I,dispstep)
        fprintf('\n%2d percents done',round(100*I/nop));
    end
    res(:,:,I)=B_ConstantShapeIntegrals(fieldpoints(I,:),sourcemesh.e,sourcemesh.p)';  
end
if verbose
    fprintf('Done.');
end
