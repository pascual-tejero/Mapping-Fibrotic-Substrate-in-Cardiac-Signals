function R=Correlation(set1,set2)
% function R=Correlation(set1,set2)
% Correlation or directional cosine between two n-dimensional vectors
% Eq. 35
d1=set1-mean(set1);
d2=set2-mean(set2);
l1=sqrt(sum(d1.*d1));
l2=sqrt(sum(d2.*d2));
R=sum(d1.*d2)/(l1*l2);

