function F = sd_interp3(varargin)
% sd_interp3 INTERP3 function with some error checking removed
% Take care that input is monotonic increasing meshgrid with equal spacing

% sd_custom
arg1= varargin{1}; %x
arg2= varargin{2}; %y
arg3= varargin{3}; %z
arg4= varargin{4}; %v
arg5= varargin{5}; %xi
arg6= varargin{6}; %yi
arg7= varargin{7}; %zi
method= varargin{8}; %not used
ExtrapVal= varargin{9}; 


%LINEAR 3-D trilinear data interpolation.
%   VI = LINEAR(X,Y,Z,V,XI,YI,ZI) uses trilinear interpolation
%   to find VI, the values of the underlying 3-D function in V
%   at the points in arrays XI, YI and ZI.  Arrays X, Y, and
%   Z specify the points at which the data V is given.  X, Y,
%   and Z can also be vectors specifying the abscissae for the
%   matrix V as for MESHGRID. In both cases, X, Y, and Z must be
%   equally spaced and monotonic.
%
%   Values of NaN are returned in VI for values of XI, YI and ZI
%   that are outside of the range of X, Y, and Z.
%
%   If XI, YI, and ZI are vectors, LINEAR returns vector VI
%   containing the interpolated values at the corresponding
%   points (XI,YI,ZI).
%
%   VI = LINEAR(V,XI,YI,ZI) assumes X = 1:N, Y = 1:M, and
%   Z = 1:P where [M,N,P] = SIZE(V).
%
%   VI = LINEAR(V,NTIMES) returns the matrix VI expanded by
%   interleaving bilinear interpolates between every element,
%   working recursively for NTIMES.  LINEAR(V) is the same
%   as LINEAR(V,1).
%
%   This function needs about 4 times SIZE(XI) memory to be
%   available.
%
%   See also INTERP3.

%{
  arg1 = 1:50; 
  arg2 = 1:25; 
  arg3 = 1:25; 
  arg4 = flow;
arg5= [21 21]; %xi
arg6= [1 25] ; %yi
arg7= [11 NaN]; %zi
  ExtrapVal = NaN;

%   arg3(1)=NaN; arg3(8)=NaN;
  yi = sd_interp3(arg1,arg2,arg3,arg4,arg5,arg6,arg7,0,NaN) 
%   plot(arg1,arg2,'o',arg3,yi,'x')



%}
  [nrows,ncols,npages] = size(arg4);
  mx = numel(arg1); my = numel(arg2); mz = numel(arg3);
%   if ~isequal([my mx mz],size(arg4)) && ...
%      ~isequal(size(arg1),size(arg2),size(arg3),size(arg4))
%     error(message('MATLAB:interp3:linear:XYZLengthMismatchV'));
%   end
  s = 1 + (arg5-arg1(1))/(arg1(mx)-arg1(1))*(ncols-1);
  t = 1 + (arg6-arg2(1))/(arg2(my)-arg2(1))*(nrows-1);
  w = 1 + (arg7-arg3(1))/(arg3(mz)-arg3(1))*(npages-1);
  
% end

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols) | isnan(s));
if ~isempty(sout), s(sout) = ones(size(sout)); end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows)| isnan(t));
if ~isempty(tout), t(tout) = ones(size(tout)); end

% Check for out of range values of w and set to 1
wout = find((w<1)|(w>npages)| isnan(w));
if ~isempty(wout), w(wout) = ones(size(wout)); end

% Matrix element indexing
nw = nrows*ncols;
ndx = floor(t)+floor(s-1)*nrows+floor(w-1)*nw;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Compute intepolation parameters, check for boundary value.
if isempty(w), d = w; else d = find(w==npages); end
w(:) = (w - floor(w));
if ~isempty(d), w(d) = w(d)+1; ndx(d) = ndx(d)-nw; end

% Now interpolate.

  F =  (( arg4(ndx).*(1-t) + arg4(ndx+1).*t ).*(1-s) + ...
        ( arg4(ndx+nrows).*(1-t) + arg4(ndx+(nrows+1)).*t ).*s).*(1-w) +...
       (( arg4(ndx+nw).*(1-t) + arg4(ndx+1+nw).*t ).*(1-s) + ...
        ( arg4(ndx+nrows+nw).*(1-t) + arg4(ndx+(nrows+1+nw)).*t ).*s).*w;


% Now set out of range values to ExtrapVal.
if ~isempty(sout), F(sout) = ExtrapVal; end
if ~isempty(tout), F(tout) = ExtrapVal; end
if ~isempty(wout), F(wout) = ExtrapVal; end
end %sd custom
