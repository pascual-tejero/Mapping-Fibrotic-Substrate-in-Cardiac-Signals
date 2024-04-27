function h=PlotMesh_IsoSurf(mesh,data,scale,cmap,colorbarflag)
% function h=PlotMesh_IsoSurf(mesh,data,scale,colormap,colorbarflag)
% scale, ncolors, colorbarflag optional
% scale: if scalar, then the scale will be [-scale scale],
%        can also be a vector [scalemin scalemax]
% colormap: a single integer M or a N x 3 array
%           if M, then jet colormap with M colors is used (default: 10)
%           if N x 3: array of RGB values, e.g. an output of a Matlab colormap
%           routine
% h: plot handle
N=nargin;

%initialization
pichandle=gcf;
set(pichandle,'DoubleBuffer','on');
set(pichandle,'RendererMode','manual');
set(pichandle,'Renderer','zbuffer');

N=nargin;
if N<4 |isempty(cmap)
    colormap(jet(10));
elseif N>3 & length(cmap)==1
    colormap(jet(cmap));
else
    colormap(cmap);
end

if N<3 | isempty(scale)
    scale=max(abs(data));
    caxis([-scale scale]);
elseif length(scale)==1
    caxis([-scale scale]);
else
    caxis([scale(1) scale(end)]);
end
    

%camera
cameratoolbar('SetCoordSys','y')
cameratoolbar('SetMode','orbit') 
camzoom(1);
axis equal;
axis off;

%the basic drawing routine
h = patch('faces',mesh.e,'vertices',mesh.p,'facevertexcdata',data,'facecolor','interp');

if N==5 & colorbarflag>0
    colorbar;
end
