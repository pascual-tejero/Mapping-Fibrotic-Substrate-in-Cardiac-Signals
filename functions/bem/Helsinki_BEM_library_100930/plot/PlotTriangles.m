function h=PlotTriangles(mesh,eleind,facecolor,facealpha,edgecolor,edgealpha)
% function h=PlotTriangles(mesh,triangle_ind,facecolor,facealpha,edgecolor,edgealpha)
% triangle_ind: indices of triangles to be plotted
% facecolor, edgecolor: optional, [r g b]
% facealpha, edgealpha: optional, between 0 and 1
% h: plot handle
n_arg=nargin;
if n_arg<3 | isempty(facecolor)
    facecolor=[.7 .7 .7];
end
if n_arg<4 | isempty(facealpha)
    facealpha=.5;
end
if n_arg<5 | isempty(edgecolor)
    edgecolor=[0 0 0];
end
if n_arg<5 | isempty(edgecolor)
    edgecolor=[0 0 0];
end
if n_arg<6 | isempty(edgealpha)
    edgealpha=1;
end

h = patch('faces',mesh.e(eleind,:),'vertices',[mesh.p(:,1),mesh.p(:,2),mesh.p(:,3)],'facecolor',facecolor,...
    'edgecolor',edgecolor,'facealpha',facealpha,'edgealpha',edgealpha);
