function omega=DLMatrix_Linear_FieldPoints(fp,sourcemesh,verbose)
% function omega=DLMatrix_Linear_FieldPoints(fp,sourcemesh,verbose)
% fieldpoints: N x 3
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

nof=size(fp,1);

p1=triangles(:,1);
p2=triangles(:,2);
p3=triangles(:,3);

shapes1=zeros(nof,not);
shapes2=zeros(nof,not);
shapes3=zeros(nof,not);

omega=zeros(nof,non);
if verbose
    disp('Creating DL matrix - calculating shape functions');
    dispstep=round(nof/10);
end
for I=1:nof,
    if verbose && ~mod(I,dispstep)
        disp(sprintf('%2d percents done',round(100*I/nof)));
    end
    shapes=DL_LinearShapeFunctions_DiffSurface(fp(I,:),triangles,nodes,normals,areas);
    shapes1(I,:)=shapes(:,1)';
    shapes2(I,:)=shapes(:,2)';
    shapes3(I,:)=shapes(:,3)';
    
end    
if verbose
    disp('Building matrix');
end
for J=1:not,
    omega(:,p1(J))=omega(:,p1(J))+shapes1(:,J);
    omega(:,p2(J))=omega(:,p2(J))+shapes2(:,J);
    omega(:,p3(J))=omega(:,p3(J))+shapes3(:,J);
end
omega=omega/(4*pi);
