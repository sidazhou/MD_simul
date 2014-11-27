function F = sd_interp1(varargin)
% sd_interp1 based on sd_interp3, by removing the 2 extra dimension
% yi = sd_interp1(x,y,xi)
% all input need to be column vector


% error(nargchk(1,9,nargin,'struct')); % allowing for an ExtrapVal

% sd_custom
arg1= varargin{1};
arg2= varargin{2}; %was 4
arg3= varargin{3}; %was 5
ExtrapVal= varargin{4}; % was 8

%{
  arg1 = 0:10; arg2 = sin(arg1); arg3 = 0:.25:10; ExtrapVal = NaN;
  arg3(1)=NaN; arg3(8)=NaN;
  yi = sd_interp1(arg1,arg2,arg3,NaN); plot(arg1,arg2,'o',arg3,yi,'x')

%}

  [nrows,ncols] = size(arg2);
%    mx = numel(arg1); 

arg3=arg3-arg1(1);
s = rem(arg3,arg1(2)-arg1(1))/(arg1(2)-arg1(1));

ndx = floor( arg3/(arg1(2)-arg1(1)) )+1;
%   s = 1 + (arg3-arg1(1))/(arg1(mx)-arg1(1))*(nrows-1);
% Check for out of range values of s and set to 1
% sout = find((s<1)|(s>=nrows));
% if ~isempty(sout), s(sout) = ones(size(sout)); end

% Matrix element indexing
% ndx = floor(s);
ind = (ndx<1)|(ndx>=nrows) | isnan(ndx);

% Compute intepolation parameters, check for boundary value.
% if isempty(s), d = s; else d = find(s==nrows); end
% if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end
% 
% s(:) = (s - floor(s));

% Now interpolate.
F(~ind) = arg2(ndx(~ind)).*(1-s(~ind))  +  arg2(ndx(~ind)+1) .*(s(~ind));  
F(ind) = ExtrapVal;
% Now set out of range values to ExtrapVal.
% if ~isempty(sout), F(sout) = ExtrapVal; end


end %sd custom