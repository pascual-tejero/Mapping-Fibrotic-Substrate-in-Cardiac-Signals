function qtrip(points, cells, data, port, angles, step, limits, fps)

if nargin < 8, fps = 10; end
if nargin < 7, limits = round([min(data(:)) max(data(:))]); end
if nargin < 6, step = round(diff(limits)/10); end
if nargin < 5, angles = [0 0 0]; end
if nargin < 4, port = 1041; end
if nargin < 3, error('Not enough input arguments.'); end

% flip orientation of faces
tmp = cells(:,1);
cells(:,1) = cells(:,2);
cells(:,2) = tmp;

pnet('closeall');
global qtriplot_port;
global qtriplot_con;
qtriplot_con = -1;
qtriplot_port = port;
try
    qtriplot('reset')
catch
    qtriplot_path = sprintf('%s/qtriplot/qtriplot.app/Contents/MacOS/qtriplot', fileparts(mfilename('fullpath')));
    system(sprintf('screen -d -m %s -p %i', qtriplot_path, port));
    pause(1);
    qtriplot('reset')
end
qtriplot('bgdcolor white')
qtriplot(points, cells)
qtriplot(sprintf('angle %i,%i,%i', round(angles)))
qtriplot(data)
qtriplot(sprintf('step %.3f', step))
qtriplot(sprintf('funscale %.3f', limits))
qtriplot('antialias true')
qtriplot(sprintf('movie %i', round(fps)))

end
