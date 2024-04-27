function h=PlotCoordinateAxes(linewidth,xstyle,ystyle,zstyle)
% function h=PlotCoordinateAxes(linewidth,xstyle,ystyle,zstyle)
% all arguments optional
% for xstyle etc see the help of plot3.m
% h: plot handles
h=gca;
xlim=get(h,'Xlim')*1.2;
ylim=get(h,'Ylim')*1.2;
zlim=get(h,'Zlim')*1.2;
ze=[0 0];

holdstate=ishold;
hold on;
if nargin==4
    h1=plot3([0 xlim(2)],ze, ze,xstyle);
    h2=plot3(ze, [0 ylim(2)],ze,ystyle);
    h3=plot3(ze,ze, [0 zlim(2)],zstyle);
else
    h1=plot3([0 xlim(2)],ze, ze,'r');
    h2=plot3(ze, [0 ylim(2)],ze,'g');
    h3=plot3(ze,ze, [0 zlim(2)],'b');
end    
if nargin>0 & ~isempty(linewidth)
    set(h1,'LineWidth',linewidth);
    set(h2,'LineWidth',linewidth);
    set(h3,'LineWidth',linewidth);
end
h=[h1;h2;h3];

if ~holdstate
    hold off;
end