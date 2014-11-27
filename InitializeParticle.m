function [xy,Vxy]=InitializeParticle(Npt,nozzleDistance,nozzleOpeningTime,targetVel,cR)
 %% Getting a initial gaussian distribution 
s = RandStream('mcg16807','Seed',0);  % fixing the random number generator
RandStream.setDefaultStream(s);

stddev = (nozzleOpeningTime*targetVel)/6; %eg full width is 60us with 360m/s velocity) = 6*stddev

xy=zeros(3,Npt);        %initialization
Vxy=zeros(3,Npt);

sdmean =  -nozzleDistance - 3*stddev;        %in meters       % the mean is          (-nozzleDistance - width of distribution/2) 

% for i=1:Npt
 xy(1,:) = (randn(Npt,1) * cR/3 + 0);           % mean at 0mm, with 99% within cR
 xy(2,:) = (randn(Npt,1) * cR/3 + 0);           % mean at 0mm, with 99% within cR
 xy(3,:) = (randn(Npt,1) * stddev + sdmean); 
 
 Vxy(1,:) = (randn(Npt,1) * 0.015*targetVel + 0)* 1e-9;                % velocity in x is 1%  
 Vxy(2,:) = (randn(Npt,1) * 0.015*targetVel + 0)* 1e-9;                % velocity in y is 1% 
 Vxy(3,:) = (randn(Npt,1) * 0.03 + 1)* targetVel * 1e-9; % velocity in z is = 5% * targetVel, this value is determined by trial and error



%% Set the target molecule index always to Npt/2, other bits of code are written
% around Npt/2
xy(1,round(Npt/2)) = 0   ;
xy(2,round(Npt/2)) = 0   ;
xy(3,round(Npt/2)) = sdmean   ;

Vxy(1,round(Npt/2)) = 0   ;  
Vxy(2,round(Npt/2)) = 0   ;
Vxy(3,round(Npt/2)) = targetVel * 1e-9   ;

% xy
% Vxy
end

