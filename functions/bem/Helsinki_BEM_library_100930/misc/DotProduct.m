function dot=DotProduct(R1,R2)
% function dot=DotProduct(R1,R2)
% optimized dot product routine; compare to Matlab function dot.m
dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
