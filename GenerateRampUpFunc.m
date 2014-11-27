function []=GenerateRampUpFunc(dt,NstepRampUp,NstepB4RampDown)
if ~exist('dt','var')
    dt=10;
end

if ~exist('NstepRampUp','var')
    NstepRampUp = ceil(53.5e-6/ 1e-9 / dt);
end

if ~exist('NstepB4RampDown','var')
    NstepB4RampDown = ceil(73.5e-6/ 1e-9 / dt);
end

x=0:NstepRampUp;   % max range we want to calculate, NstepRampUp can be substituted with NstepB4RampDown

% the shape of result highly depends on how many steps in x we have, since 
% are using element, x(2) to calculate the
% constants c and m, and value x(2) changes
% depending on the step size

y = ones(size(x,2),1);
x0 = NstepRampUp;  % ram to full in

x_tmp=x(x<x0); %for getting the start and the end point of the ramp
% %% solving based on prior knowledge (should really be a+b*log(x+c), and not a+b*log(x)
% % simulaneous eq: MX=Y
% % c*log(x0)+m = 1
% % c*log(0)+m = 0
% 
% M = [[log(x_tmp(2)) 1 ];[log(x_tmp(end)) 1]]; % offset from index=1, to avoid log(0) = -Inf
% 
% Y = [0 1];  % and function at start point equals to 0%
%             % require that function at end point equals to 100%
%             
% 
% X = M\Y'; % same as inv(M)*Y', solving simultaneous equation to get X = [c m]
% 
% 
%  y(x<x0) = log(x(x<x0))* X(1)+X(2);
% %  plot(x,y,'x'); hold on; plot(x,y);
% 
% RpUpTimeArray = x;
% RpUpPercentArray = y;
% 
% %% 70us rising time, fitted to experiment, coeff given by yang
% xGiven=0:1e-6:70e-6; 
% a = 0.15648; b=0.0126; c=4.826e-6;
% yGiven = a+b.*log(xGiven+c);
% yGiven = yGiven./max(yGiven);
% 
% xGiven = xGiven/1e-9/50;  % converting the unit from SI to my unit
% y = interp1(xGiven,yGiven,x); % interpolating at integer times
% 
% % plot(x,y,'.')
% 
% % check that  RpUpPercentArray() - 1 < eps; or there will be trouble
% 
% RpUpTimeArray = x;
% RpUpPercentArray = y;
% 
% RpUpPercentArray(isinf(RpUpPercentArray)) = 0; % removing infinity
% 
% %manual cleaning, dont ask why
% RpUpPercentArray(end)=[];
% RpUpTimeArray(end)=[];
% RpUpPercentArray(end)=1;
% RpUpTimeArray(end)=1;
% 
% save(['RampUpFunc' num2str(dt) '.mat'], 'RpUpTimeArray', 'RpUpPercentArray');

%% 73.5us total time, fitted to experiment, coeff given by yang
% note: NEED TO MATCH THIS NUMBER 53.5us, to main.m, or NaN errors will occur
if NstepRampUp ~= ceil(53.5e-6/ 1e-9 / dt)
    disp('error, check GenerateRampUpFunc.m line 72')
end
xGiven=0:0.1e-6:53.5e-6;  %MANUAL
a = 0.13547; b=0.01075; c=4.144e-6;
yGiven = a+b.*log(xGiven+c);
yGiven = yGiven./max(yGiven);

xGiven = xGiven/1e-9/dt;  % converting the unit from SI to my unit
y = interp1(xGiven,yGiven,x); % interpolating at integer times

% plot(x,y,'.')

%% output
RpUpTimeArray = x;
RpUpPercentArray = y;
%manual cleaning, dont ask why
RpUpPercentArray(isinf(RpUpPercentArray)) = 0; % removing infinity
RpUpPercentArray = RpUpPercentArray(1:NstepRampUp);
RpUpTimeArray = RpUpTimeArray(1:NstepRampUp);
RpUpPercentArray=RpUpPercentArray(1:NstepRampUp);
RpUpPercentArray(1) = 0;
RpUpPercentArray(NstepRampUp+1)=1; % check that  RpUpPercentArray() - 1 < eps; or there will be trouble

if(any(isnan(RpUpPercentArray))), disp('warning, NaN in RpUpPercentArray'); end

save(['RampUpFunc' num2str(dt) '.mat'], 'RpUpTimeArray', 'RpUpPercentArray');
end % end function
