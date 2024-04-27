function omega=SLMatrix_Linear(fieldmesh,sourcemesh,verbose)
% function omega=SLMatrix_Linear(fieldmesh,sourcemesh,verbose)
% verbose: 1 for some output
% Calculates single-layer integrals
if nargin==2
    verbose=0;
end

non=sourcemesh.nop;
not=sourcemesh.noe;
triangles=sourcemesh.e;
nodes=sourcemesh.p;
normals=sourcemesh.un;
areas=sourcemesh.a;

p1=triangles(:,1);
p2=triangles(:,2);
p3=triangles(:,3);

fieldpoints=fieldmesh.p;
nof=size(fieldpoints,1);

omega=zeros(nof,non);
shapes1=zeros(nof,not);
shapes2=zeros(nof,not);
shapes3=zeros(nof,not);

dq = parallel.pool.DataQueue;
if verbose
    disp('Creating SL matrix - calculating shape functions');
    progressStep=floor(nof/10);
    progressCount=0;
    afterEach(dq, @printProgress);
end
function printProgress(~)
    progressCount=progressCount+1;
    if ~mod(progressCount,progressStep)
        fprintf('%2d percents done\n', 10*progressCount/progressStep);
    end
end

parfor I=1:nof
    shapes = SL_LinearShapeFunctions(fieldpoints(I,:),sourcemesh);
    shapes1(I,:)=shapes(:,1)';
    shapes2(I,:)=shapes(:,2)';
    shapes3(I,:)=shapes(:,3)';
    send(dq, 0);
end

for J=1:not
    omega(:,p1(J))=omega(:,p1(J))+shapes1(:,J);
    omega(:,p2(J))=omega(:,p2(J))+shapes2(:,J);
    omega(:,p3(J))=omega(:,p3(J))+shapes3(:,J);
end
omega=omega/(4*pi);

end