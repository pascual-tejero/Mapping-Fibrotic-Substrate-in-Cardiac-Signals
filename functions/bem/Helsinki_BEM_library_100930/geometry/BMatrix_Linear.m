function omega=BMatrix_Linear(fieldpoints,sourcemesh,verbose)
% function omega=BMatrix_Linear(fieldpoints,sourcemesh,verbose)
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

non=sourcemesh.nop;
not=sourcemesh.noe;
triangles=sourcemesh.e;

p1=triangles(:,1);
p2=triangles(:,2);
p3=triangles(:,3);

nof=size(fieldpoints,1);

omeganewx=zeros(non,nof);
omeganewy=zeros(non,nof);
omeganewz=zeros(non,nof);

if verbose
    fprintf('Creating basis functions for calculation of B...');
    dispstep=round(nof/10);
end

S1x=zeros(nof,not);
S2x=zeros(nof,not);
S3x=zeros(nof,not);

S1y=zeros(nof,not);
S2y=zeros(nof,not);
S3y=zeros(nof,not);

S1z=zeros(nof,not);
S2z=zeros(nof,not);
S3z=zeros(nof,not);

for I=1:nof,
    if verbose2 && ~mod(I,dispstep)
        fprintf('\n%2d percents done',round(100*I/nof));
    end
    S=B_LinearShapeIntegrals(fieldpoints(I,:),sourcemesh);
    S1x(I,:)=S(:,1,1);
    S1y(I,:)=S(:,1,2);
    S1z(I,:)=S(:,1,3);
    S2x(I,:)=S(:,2,1);
    S2y(I,:)=S(:,2,2);
    S2z(I,:)=S(:,2,3);
    S3x(I,:)=S(:,3,1);
    S3y(I,:)=S(:,3,2);
    S3z(I,:)=S(:,3,3);
end
for J=1:not,
    omeganewx(p1(J),:)=omeganewx(p1(J),:)+S1x(:,J)';
    omeganewx(p2(J),:)=omeganewx(p2(J),:)+S2x(:,J)';
    omeganewx(p3(J),:)=omeganewx(p3(J),:)+S3x(:,J)';
    omeganewy(p1(J),:)=omeganewy(p1(J),:)+S1y(:,J)';
    omeganewy(p2(J),:)=omeganewy(p2(J),:)+S2y(:,J)';
    omeganewy(p3(J),:)=omeganewy(p3(J),:)+S3y(:,J)';
    omeganewz(p1(J),:)=omeganewz(p1(J),:)+S1z(:,J)';
    omeganewz(p2(J),:)=omeganewz(p2(J),:)+S2z(:,J)';
    omeganewz(p3(J),:)=omeganewz(p3(J),:)+S3z(:,J)';
end
omega=zeros(3,non,nof);
omega(1,:,:)=omeganewx;
omega(2,:,:)=omeganewy;
omega(3,:,:)=omeganewz;
if verbose
    fprintf('Done.');
end

    
 