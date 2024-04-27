function qtriplot(varargin);
%
% QTRIPLOT - Function to start and send data to the qtriplot application.
%
% -- Usage -------------------------------------------------------
%
% qtriplot(surfStruct [,'name'])
% 
% Load the surface described by SurfStruct into the qtriplot application.
% SurfStruct should be a structure containing the fields .pnt and .tri
% (or .node and .face), describing the vertices and triangles respectively.
% If supplied, 'name' will be the name of the geometry in the qtriplot
% application. Otherwise a default name wil be supplied: 'def_1' for the 
% first object, etc.
% Note: if surfStruct does not contain a .tri field, qtriplot will
% plot the points in the .pnt field as dots.
%
% qtriplot(ver, tri, [,'name'])
% 
% Load the surface described by the vertex array ver and the triangle array tri
% into the qtriplot application.
% If supplied, 'name' will be the name of the geometry in the qtriplot
% application. Otherwise a default name wil be supplied: 'def_1' for the 
% first object, etc.
% Note: if tri==[], qtriplot will plot ver as dots.%
%
% In the following cases, 'name' (if supplied) refers to the geometry in 
% qtriplot to which the object described by the first parameter(s) will be 
% assigned. If 'name' is not supplied,  the object will be assigned to 
% the surface most recently defined.
%
% qtriplot(fun, [,'name'])
% 
% Load the matrix fun as a function structure into the qtriplot application.
% Note: if fun is a 1-column, integer array, it is considered to be an
% electrode array rather than a function (see next).
%
% qtriplot(int32(el), [,'name'])
% 
% Load the vertex-electrode description el into the qtriplot 
% application (see qtriplot help for details).
%
% qtriplot(index, lm, [,'name'])
% 
% Load the lm-electrode description defined by index and lm into the 
% qtriplot application (see qtriplot help for details).
%
% qtriplot(elStruct, [,'name'])
% 
% Load the electrode description elStruct into the qtriplot application.
% elStruct must contain a .index field, and may contain a .lm field.
% if it contains no .lm field, it represents a vertex-electrode descrition,
% otherwise it contains a lm-electrode description (see qtriplot help for 
% details).
%
% qtriplot('mask', int32(mask), [,'name'])
% 
% Load the mask description mask into the qtriplot 
% application (see qtriplot help for details).
%
% qtriplot('command')
% 
% Send 'command' to the qtriplot application. 
% Qtriplot will try to execute 'command' as a qtriplot command.
%
% -- Background -------------------------------------------------------
%
% The qtriplot application is a versatile program to plot triangulated 
% surfaces and functions on those surfaces 
% (see http://www.ecgsim.org/qtriplot).
%
% The function qtriplot first determines wether the qtriplot application is
% running. If not, it will start the application, setting up tcp communication.
% The default tcp port is 1041. This can be overridden by setting
% the variable qtriplot_port.
%
% The program qtriplot must be in the PATH environment variable in order for
% this script to be able to start it. As an alternative, you can define the 
% global variable qtriplot_path to point to the location of the program.
% Note: on MacOS this must point to the directory where the actual program 
% resides within the application folder:  qtriplot.app/Contents/MacOS
%
% See the manual of the qtriplot application for information how to
% operate that application.

% a note on byte order:
% all integers and double must be sent in big-endian, which is the network
% standard.
% Note that the manual of the tcp-udp-ip toolbox, from which the function pnet
% is used, uses the incorrect (i.e. swapped) definition of little-endian and
% big-endian.


% tcp-port to communicate with qtriplot application
global qtriplot_port;
global qtriplot_con;

if size(qtriplot_port)==0
  qtriplot_port=1042; % default value for port; should be unassigned
  %qtriplot_port=119; % default value for port; should be unassigned
end
port=qtriplot_port;
% misschien hier een pause van 0.1 om te voorkomen dat qtriplot stymied raakt
% dat komt doordat ik vroeger bij elke call van qtriplot.m een
% socket opende en weer sloot. Als het niet snel genoeg werd
% afgewerkt, kwamen de commando's niet in dezelfde volgorde aan als
% ze verstuurd werden. Daarom open ik nu eenmaling een socket, en
% hou die open.

% check wether there is already an connection

if (size(qtriplot_con)==0)
  qtriplot_con=-1;
end
con=qtriplot_con;
if (con>=0) % check wether still open
  try
%    fprintf('trying\n');
    writeCon(con, '%');               % have to write someting
    if (pnet(con, 'status')==0)  % otherwise 'status' will be
%      fprintf('status==0\n');    % nonzero if qtriplot stopped
      pnet('closeall');
      con=-1;
%    else
% unfortunately, pnet(con, 'write') does not return anything
% status returns nonzero
%      fprintf('attempting to write\n');
%      pnet(con, 'write', int32(length('%')))
%      pnet(con, 'write', '%')
    end
  catch err
%    fprintf('catching\n');
    con=-1;
    pnet('closeall');
  end
end
%fprintf('ready\n');
%con
  
% try to connect to qtriplot application via port

% watch it: Matlab has its own version of the Qt-frameworks in
% /Applications/MATLAB_R2011a.app/bin/maci64, that are propably
% incompatible with the qtriplot version. You can overwrite those
% im maci64 by the qtriplot ones


if (con<0)
%   fprintf('Opening a new connection\n');
  con=pnet('tcpconnect', 'localhost', port);
end

% if unsuccessful, start qtriplot application
if (con<0)
  fprintf('Starting program qtriplot...\n');
  global qtriplot_path;
  if size(qtriplot_path)==0
    res=system(['qtriplot -p ' mat2str(port) ' &']);
  else
    res=system([qtriplot_path '/qtriplot -p ' mat2str(port) ' &']);
  end
  if res~=0
    error('Cannot start program qtriplot. ');
  end
  fprintf('Establishing tcp connection with program qtriplot...\n');
  ntry=20;
  for i=1:ntry
    con=pnet('tcpconnect', 'localhost', port);
    if con>=0
      break;
    end
    pause(0.1);
  end

  if i>=ntry
    error(['Cannot establish connection to qtriplot on port ' mat2str(port)]);
  end
  fprintf('Connection established.\n');
end

qtriplot_con=con;

% Send datagram to qtriplot application.
% datagram consists of uint4, indiciting size of datablock in bytes,
% followed by datablock.
% Datablock always starts with a \0-ended string
% if string starts with |, the datablock contains matlab data,
% otherwise the string is send to qtriplot, that will accept it as a command.

% Determine number and type of arguments, and from that what action to take

if length(varargin)==0
  return;
end

type=0;

% a single string is a command.
if length(varargin)==1
  if ischar(varargin{1})
    type = 1;			% qtriplot command
  end
end

% if 2 matrices, and optionally a string,
% the matrices are the vertex and index arrays of a triangulate surface.
% the name will be used in qtriplot as the name for the object.
if (length(varargin)==2) || ((length(varargin)==3) && (ischar(varargin{3})))
  if (isnumeric(varargin{1})) && (isnumeric(varargin{2}))
    [dummy p3] = size(varargin{1});
    [dummy t3] = size(varargin{2});
    if (p3==3) && ( (t3==3) || (t3==0) )
      type = 11;		% triangled surface in 2 arrays
    end
  end
end

% if a structure (containing at least .pnt field), and optionally a string,
% it is a triangulate surface.
% the name will be used in qtriplot as the name for the object.
if (length(varargin)==1) || ((length(varargin)==2) && (ischar(varargin{2})))
  if (isstruct(varargin{1})) && (isfield(varargin{1},'pnt'))
    [dummy p3] = size(varargin{1}.pnt);
    if isfield(varargin{1}, 'tri')
      [dummy t3] = size(varargin{1}.tri);
      if dummy==0
        t3=3;   % no triangles; prevent error message
      end
    else
      t3=3;
    end
    if (p3==3) && (t3==3)
      type = 12;		% triangled surface in structure (pnt/tri)
    end
  end
end

% if a structure (containing at least .node field), and optionally a string,
% it is a triangulate surface.
% the name will be used in qtriplot as the name for the object.
if (length(varargin)==1) || ((length(varargin)==2) && (ischar(varargin{2})))
  if (isstruct(varargin{1})) && (isfield(varargin{1},'node'))
    [dummy p3] = size(varargin{1}.node);
    if isfield(varargin{1}, 'face')
      [dummy t3] = size(varargin{1}.face);
      if dummy==0
        t3=3;   % no triangles; prevent error message
      end
    else
      t3=3;
    end
    if (p3==3) && (t3==3)
      type = 13;		% triangled surface in structure (node/face)
    end
  end
end

% a single matrix, optionally followed by a string,
% is a function matrix and the name of the corresponding geometry in qtriplot.
% except when the array is of integer type, then it is an vertex electrode
% array.
if (length(varargin)==1) || ((length(varargin)==2) && ischar(varargin{2}))
  if (isnumeric(varargin{1}))
    type = 21;			% function on triangulated surface
    if (isa(varargin{1}, 'int16')) ...
       || (isa(varargin{1}, 'uint16')) ...
       || (isa(varargin{1}, 'int32')) ...
       || (isa(varargin{1}, 'uint32'))
      [n m]=size(varargin{1});
      if (m==1) || (n==1)
        type = 31;		% vertex electrode array
      end
    end
  end
end

% two matrices, the first containing 1 column and the second 2, having the same
% number of rows, andoptionally followed by a string,
% are an electrode description and the name of the corresponding geometry in 
% qtriplot.
% the first array are the indices of the triangles inside which an electode is
% positioned, and the second one contains the lambda-mu parameters describing
% the position within the triangle (see qtriplot help for more details).
if (length(varargin)==2) || ((length(varargin)==3) && ischar(varargin{3}))
  if isnumeric(varargin{1}) && isnumeric(varargin{2})
    [n m1]=size(varargin{1});
    [n m2]=size(varargin{2});
    if m1==1 && m2==2
      type = 32;		% lm electrode array
    end
  end
end

% if a structure (containing at least .index field), and optionally a string,
% it is an electrode description.
% the name will be used in qtriplot as the name for the object.
if (length(varargin)==1) || ((length(varargin)==2) && (ischar(varargin{2})))
  if (isstruct(varargin{1})) && (isfield(varargin{1},'index'))
    [nel i1] = size(varargin{1}.index);
    if isfield(varargin{1}, 'lm')
      [nl l2] = size(varargin{1}.lm);
    else
      l2=2;
      nl==nel;
    end
    if (i1==1) && (l2==2) && (nl==nel)
      type = 33;		% electrode description in structure
    end
  end
end


% if the string 'mask', an array and optionally a string,
% it is a mask description.
if (ischar(varargin{1}))
  if (strcmp(varargin{1},'mask'))
    if (length(varargin)==2) || ((length(varargin)==3) && (ischar(varargin{3})))
      if (isa(varargin{2}, 'int16')) ...
         || (isa(varargin{2}, 'uint16')) ...
         || (isa(varargin{2}, 'int32')) ...
         || (isa(varargin{2}, 'uint32'))
        [n m]=size(varargin{2});
        if (m==1) || (n==1)
          type = 41;		% mask array
        end
      end
    end
  end
end


% if there are more than 1 parameters, and the last one is a string,
% it is the name to use for the qtriplot object.
% retreive that name.
if (length(varargin)>1) && (ischar(varargin{length(varargin)}))
  name=varargin{length(varargin)};
else
  name='';
end

% qtriplot command.
if (type==1)
  writeCon(con, varargin{1});

% function matrix
elseif (type==21)
  [nr nc]=size(varargin{1});
  command=['|fun ' name char(0)];
  while (mod(length(command),8)~=0)
    command=[command char(0)];
  end
  tot=length(command)+8+8*nr*nc;
  nr=int32(nr);
  nc=int32(nc);
  tot=int32(tot);
  pnet(con, 'write', tot);
  pnet(con, 'write', command);
  pnet(con, 'write', nr);
  pnet(con, 'write', nc);
  f=double(varargin{1});
  pnet(con, 'write', f);

% triangulated surface
elseif (type==11 || type==12 || type==13)
  if (type==11)
    ver=varargin{1};
    tri=varargin{2};
  elseif (type==12)
    ver=varargin{1}.pnt;
    if isfield(varargin{1},'tri')
      tri=varargin{1}.tri;
    else
      tri=zeros(0);
    end
  elseif (type==13)
    ver=varargin{1}.node;
    if isfield(varargin{1},'face')
      tri=varargin{1}.face;
    else
      tri=zeros(0);
    end
  end
  [nver dim]=size(ver);
  if (dim~=3)
    error(['Vertex array must have 3 columns']);
  end
  if isempty(tri)
    ntri=0;
  else
    [ntri dim]=size(tri);
    if (dim~=3)
      error(['Index array must have 3 columns']);
    end
  end
  command=['|tri ' name char(0)];
  while (mod(length(command),8)~=0)
    command=[command char(0)];
  end
  tot=length(command)+8+8*nver*3+4*ntri*3;
  nver=int32(nver);
  ntri=int32(ntri);
  tot=int32(tot);
  pnet(con, 'write', tot);
  pnet(con, 'write', command);
  pnet(con, 'write', nver);
  pnet(con, 'write', ntri);
  ver=double(ver);
  pnet(con, 'write', ver);
  if ~isempty(tri)
    tri=int32(tri);
    pnet(con, 'write', tri);
  end

% electrode array
elseif (type==31 || type==32 || type==33)
  lm=zeros(0);
  if (type==31 || type == 32)
    index=varargin{1};
    [n m] = size(index);
    if (m~=1)
      index=index';
    end
    if type==32
      lm=varargin{2};
      if (m~=1)
        lm=lm';
      end
    end
  else
    index=varargin{1}.index;
    if isfield(varargin{1}, 'lm')
      lm=varargin{1}.lm;
    end
  end
  [nel dim]=size(index);
  if (dim~=1)
    error(['Electrode index array must have 1 column']);
  end
  if isempty(lm)
    lmtype=0;
  else
    lmtype=1;
    [nlm dim]=size(lm);
    if (dim~=2)
      error(['Lambda-mu array must have 2 columns']);
    end
    if (nlm~=nel)
      error(['Index and Lambda-mu must have same number of rows']);
    end
  end
  command=['|ele ' name char(0)];
  while (mod(length(command),8)~=0)
    command=[command char(0)];
  end
  tot=length(command)+8+4*nel;
  if ~isempty(lm)
    tot = tot+2*8*nel;
  end
  nel=int32(nel);
  lmtype=int32(lmtype);
  tot=int32(tot);
  pnet(con, 'write', tot);
  pnet(con, 'write', command);
  pnet(con, 'write', nel);
  pnet(con, 'write', lmtype);
  if ~isempty(lm)
    lm=double(lm);
    pnet(con, 'write', lm);
  end
  index=int32(index);
  pnet(con, 'write', index);

% mask array
elseif (type==41)
  mask=varargin{2};
  [n m] = size(mask);
  if (m~=1)
    mask=mask';
  end
  [nmask dim]=size(mask);
  if (dim~=1)
    error(['Mask array must have 1 column']);
  end
  command=['|mask ' name char(0)];
  while (mod(length(command),8)~=0)
    command=[command char(0)];
  end
  tot=length(command)+4+4*nmask;
  nmask=int32(nmask);
  tot=int32(tot);
  pnet(con, 'write', tot);
  pnet(con, 'write', command);
  pnet(con, 'write', nmask);
  mask=int32(mask);
  pnet(con, 'write', mask);

else
  error('Invalid parameters in qtriplot');
end

%pnet(con, 'close');
return;

%------------------
function writeCon(con, x)
pnet(con, 'write', int32(length(x)));
pnet(con, 'write', x);
return;