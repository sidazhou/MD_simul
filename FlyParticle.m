function [targetPtArrivalTime,maxtime,cSeq]=FlyParticle(firstRunFlag,sdParam,currentFolder,maxtime,cR,Npt,NstepRampUp,NstepRampDown,mass,DetectorDistance,DetectorRadius,raleyLength,nozzleOpeningTime,sdPhase,dt,Ncoil,cSpacing,xy,Vxy,cDim,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,cSeq)
%% Necessary Initializations
load(['RampUpFunc' num2str(dt) '.mat']); % getting RpUpTimeArray and RpUpPercentArray

if(ispc()) %if PC
load 'C:\Users\sidz\Desktop\ZeemanShiftO2.mat'
else %on westgrid, which returns ''
load '~/ZeemanShiftO2.mat'
end

cFlag = zeros(1,Ncoil);                         % Memory pre-allocation
fpw = fopen([currentFolder '/sdPlotThis' num2str(sdParam) '.txt'],'wt');
 targetPtArrivalTime = [];

fvcCounter=1; % initiation, velocity after each coil counter
fpfvc = fopen([currentFolder '/velocityAfterEachStage.txt'],'wt'); % velocity after each coil

% fpp = fopen('./TESTcFlag.txt','wt'); %TEST cFlag
% visualization: 
% for i=1:15
%     plot(TESTcFlag(:,i)+2*i); hold on
% end

ptArrivedDetectorFlag = zeros(1,Npt);
ptIndex = 1:Npt; ptIndex = ptIndex';  % for writing sdDetectorClickTimes*.txt

%% MAIN LOOP
for t = 0 : dt :maxtime                       % current time
%% Update position and velocity, storing into the ptBoxNold, then updating the new ptBoxN
% fprintf(fpp,repmat('%+e\t',1,15),cFlag); %TEST cFlag
% fprintf(fpp,'\n');                       %TEST cFlag


 cSeq=GetcSeq(t,Ncoil,Npt,sdPhase,cSpacing,cDim,xy,cSeq,currentFolder);
 cFlag = GetFlag(t,Ncoil,NstepRampUp,NstepRampDown,cFlag,cSeq,RpUpPercentArray);                    % Gets the new cFlag, ie get which coils are flagged to be turned on and which are turned off

% check if it hits outer wall, there are two different indications
ind1=sqrt(xy(1,:).^2+ xy(2,:).^2) > cR  &  xy(3,:) >= 0    &   xy(3,:) < cDim(3)*Ncoil + sum(cSpacing(1:Ncoil-1));
ind2=isnan(xy(1,:)) | isnan(xy(2,:))| isnan(xy(3,:))| isnan(Vxy(1,:))| isnan(Vxy(2,:))| isnan(Vxy(3,:));
ind0 = (ind1 | ind2);

xy(1,ind0)= cR + eps;
xy(2,ind0)= cR + eps;
xy(3,ind0)= -eps;
Vxy(1,ind0)=0;
Vxy(2,ind0)=0;
Vxy(3,ind0)=0;

%     t
%   xy
%  Vxy
%    xy(:,round(Npt/2))
%  Vxy(3,round(Npt/2))/1e-9
%   dbstop in FlyParticle at 54 if t==271000
%{
 dbstop in FlyParticle at 56 if Vxy(3,round(Npt/2))/1e-9<=0 
     for k = 1:Ncoil
    ind = (xy(3,round(Npt/2)) >= cDim(3)*(k-1) + sum(cSpacing(1:k-1))) & (xy(3,round(Npt/2)) < cDim(3)*k + sum(cSpacing(1:k-1))); 
     if(ind==1) disp(k); end; 
    end
%}

[xy(:,~ind0),Vxy(:,~ind0)] = rk4Update(sdParam,cR,dt,mass,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,xy(:,~ind0),Vxy(:,~ind0),BfieldSummary,eigenValuesDerivativeSummary);

% recording velocity after each coil of the target particle, base on the position of the target particle
if (firstRunFlag==0 && fvcCounter <= Ncoil && sdParam == 9 || firstRunFlag==1 && fvcCounter <= Ncoil )  % MANUAL
    if (xy(3,round(Npt/2)) >= cDim(3)*(fvcCounter-1) + sum(cSpacing(1:fvcCounter-1))) && (xy(3,round(Npt/2)) < cDim(3)*fvcCounter + sum(cSpacing(1:fvcCounter-1)));   
    fprintf(fpfvc,'%d\t%.2f\n',fvcCounter,Vxy(3,round(Npt/2))/1e-9);
    fvcCounter = fvcCounter + 1;
    end
end
    
    
% post check for particles that are moved, if they have any NaN values, which screws up the file output
ind3=isnan(xy(1,~ind0)) | isnan(xy(2,~ind0)) | isnan(xy(3,~ind0))| isnan(Vxy(1,~ind0))| isnan(Vxy(2,~ind0))| isnan(Vxy(3,~ind0));
ind4=sqrt(xy(1,~ind0).^2+ xy(2,~ind0).^2) > cR  &  xy(3,~ind0) >= 0    &   xy(3,~ind0) < cDim(3)*Ncoil + sum(cSpacing(1:Ncoil-1));
ind5=(ind3 | ind4);

xy(1,ind5)= cR + eps;
xy(2,ind5)= cR + eps;
xy(3,ind5)= -eps;
Vxy(1,ind5)=0;
Vxy(2,ind5)=0;
Vxy(3,ind5)=0;


%% Recording sdDetectorClickTimes.txt
    ind_notArrived_N_justArrival = ptArrivedDetectorFlag==0 & sqrt((xy(3,:)-(cDim(3)*Ncoil + sum(cSpacing(1:Ncoil-1)) +DetectorDistance)   ).^2 + xy(2,:).^2) <=DetectorRadius    &  abs(xy(1,:))   <=  raleyLength;
% %
% fprintf(fpt,'%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\n',n,t*1e-9,xy(1,n),xy(2,n),xy(3,n),Vxy(1,n)/1e-9,Vxy(2,n)/1e-9,Vxy(3,n)/1e-9); % not in use %
    tmp_writeThis = [ptIndex(ind_notArrived_N_justArrival) t*ones(size(xy(:,ind_notArrived_N_justArrival),2),1)*1e-9 xy(:,ind_notArrived_N_justArrival)' Vxy(:,ind_notArrived_N_justArrival)'/1e-9]; %#ok<NASGU>
    save([currentFolder '/sdDetectorClickTimes' num2str(sdParam) '.txt'],'tmp_writeThis','-ascii','-append'); %extract the time at which pt crosses the detector
    ptArrivedDetectorFlag(ind_notArrived_N_justArrival)=1; % mark the ones that has arrived
    clear tmp_writeThis;

% sdout.txt format;
% ptN t1 x1 y1 Vx1 Vy1
% .. .. .. ... ...
% ptN t1 xn yn Vxn Vyn
%
% ptN t2 x1 y1 Vx1 Vy1 
% .. .. .. ... ...
% ptN t2 xn yn Vxn Vyn

        %set maxtime
    if (sqrt((xy(3,round(Npt/2))-(cDim(3)*Ncoil + sum(cSpacing(1:Ncoil-1)) +DetectorDistance)   ).^2 + xy(2,round(Npt/2)).^2) <=DetectorRadius    &&    abs(xy(1,round(Npt/2)))   <=  raleyLength  &&  isempty(targetPtArrivalTime)==1 )   
      targetPtArrivalTime = t;
      targetPtArrivalTime-nozzleOpeningTime/2  %this is arrival time if target is at nozzle opening at t=0 (eg instead of packet front at nozle opening)
      Vxy(3,round(Npt/2)) /1e-9 
      
      if(~exist([currentFolder './sdoutMaxTime.txt'],'file'))
        sdAdditionalTime = max([0.3e-3 nozzleOpeningTime]); %additional time after pt reaches detector, 0.3ms is arbitrary, lokos good
        maxtime=t+(sdAdditionalTime /1e-9);  %record for 0.7ms after target pt exits
      end
    end
        
        %stopping the simulation
        if (t >= maxtime)
        break;
        end
        
%% Recording the trajectory sdPlotThis.txt       

% if (exist('sdoutMaxTime.txt','file'))  %only record when   % NOTE, THIS LINE DOESNT WORK ON THE CLUSTER, REASON UNKNOWN       
    timeRecordInterval=round(maxtime/10); %determine how many points to plot, lets say 10
    timeRecordInterval=timeRecordInterval-mod(timeRecordInterval,dt);  %rounding to integer arb units

    if (    round(mod(t , timeRecordInterval)     )==0 )   %extract 100 points for plotting
        for n = 1 : Npt
        fprintf(fpw,'%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\t%+e\n',n,t*1e-9,xy(1,n),xy(2,n),xy(3,n),Vxy(1,n)/1e-9,Vxy(2,n)/1e-9,Vxy(3,n)/1e-9);  %extraxt 1000 points for plotting
        end
    end
% end




end   %END OF MAIN LOOP



fclose(fpw); 

fclose(fpfvc); disp('fclose successful.')
end

function [xy,Vxy]=rk4Update(sdParam,cR,dt,mass,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,xy,Vxy,BfieldSummary,eigenValuesDerivativeSummary)

% RK4 for position and velocity simultaneously (from Natsunari)
k1x = Vxy;
k1v = GetForce(sdParam,cR,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,xy,BfieldSummary,eigenValuesDerivativeSummary)/mass;

k2x = Vxy + 0.5 * k1v * dt;
k2v = GetForce(sdParam,cR,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,(xy+0.5*k1x*dt),BfieldSummary,eigenValuesDerivativeSummary)/mass;

k3x = Vxy + 0.5 * k2v * dt;
k3v = GetForce(sdParam,cR,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,(xy+0.5*k2x*dt),BfieldSummary,eigenValuesDerivativeSummary)/mass;

k4x = Vxy + k3v * dt;
k4v = GetForce(sdParam,cR,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,(xy+k3x*dt),BfieldSummary,eigenValuesDerivativeSummary)/mass;

xy = xy + 1/6 * dt * ( k1x + 2*k2x + 2*k3x + k4x );
Vxy = Vxy + 1/6 * dt * ( k1v + 2*k2v + 2*k3v + k4v );


end




