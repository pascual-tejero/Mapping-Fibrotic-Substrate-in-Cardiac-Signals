function scale=SetMapScale(datamax,ncols)
% function scale=SetMapScale(datamax,ncols)
% Computes a "reasonable" value for color scaling. 
% datamax: maximum value of the data / reference value for scaling
% ncols:   number of colors in the colormap
step=2*datamax/ncols;
steplog=log10(step);
gain=10^(floor(steplog));
steptemp=step/gain;
stepvec=[1 2 2.5 5 10];
[foo,ind]=min(abs(steptemp-stepvec));
step=stepvec(ind)*gain;
scale=step*ncols/2;