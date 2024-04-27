function omega=DLMatrix_Linear_DiffSurface(fieldmesh,sourcemesh,verbose)
% function omega=DLMatrix_Linear_DiffSurface(fieldmesh,sourcemesh,verbose)
% verbose: 1 for some output
% Eq. 33
if nargin==2
    verbose=0;
end

non=sourcemesh.nop;
not=sourcemesh.noe;
triangles=sourcemesh.e;
nodes=sourcemesh.p;
normals=sourcemesh.n;
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
    disp('Creating DL matrix - calculating shape functions');
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
    shapes=DL_LinearShapeFunctions_DiffSurface(fieldpoints(I,:),triangles,nodes,normals,areas);
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

[~,ind1,ind2] = overlappingSurfaces(fieldmesh,sourcemesh);
for I=1:length(ind1)
    temp=sum(omega(ind1(I),:));
    omega(ind1(I),ind2(I))= -.5-temp + .5;
end

end