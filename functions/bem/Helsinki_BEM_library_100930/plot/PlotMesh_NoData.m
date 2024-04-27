function h=PlotMesh_NoData(mesh,facecolor,facealpha,edgecolor,edgealpha)
%function h=PlotMesh_NoData(mesh,facecolor,facealpha,edgecolor,edgealpha)
%facecolor and edgecolor: RGB values between [0 0 0] and [1 1 1] or color
%symbols understood by Matlab, e.g., 'red', 'blue' (optional)
%facealpha and edgealpha: scalars between 0 and 1. (optional)
%h: plot handle
%If alpha is used, colors have to be given as well.
%If alpha parameters are set, this routine needs OpenGL.
%for further information see help for the trimesh function.

nar=nargin;
if nar==1
    h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor','black');
elseif nar==2 | (nar ==3 & isempty(facealpha))
    h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor','black','FaceColor',facecolor);
elseif nar==3
    set(gcf,'Renderer','OpenGL');
    h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor','black','FaceColor',facecolor,'FaceAlpha',facealpha);
elseif nar==4
    if isempty(facealpha)
        h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor);
    else    
        set(gcf,'Renderer','OpenGL');
        h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor,'FaceAlpha',facealpha);
    end
else
    if isempty(facealpha) & isempty(edgealpha)
        h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor);
    else
        set(gcf,'Renderer','OpenGL');
        if isempty(facealpha)
            h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor,'EdgeAlpha',edgealpha);
        elseif isempty(edgealpha)
            h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor,'FaceAlpha',facealpha);
        else
            h=trimesh(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'EdgeColor',edgecolor,'FaceColor',facecolor,'FaceAlpha',facealpha,'EdgeAlpha',edgealpha);
        end
    end    
end
cameratoolbar('SetCoordSys','y')
cameratoolbar('SetMode','orbit') 
camzoom(1);
axis equal;
axis off;
