function R=CrossProduct(R1,R2)
% function R=CrossProduct(R1,R2)
% optimized cross product routine; compare to Matlab function cross.m

R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
