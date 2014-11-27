function [DetectorRadius,raleyLength,cDim,cR,sdPhase,cSpacing]=main(firstRunFlag,sdParam,currentFolder,Npt)
%% True Main program
    if ~exist('Npt','var')
        Npt =3; %for testing purposes
    end
    
    if (exist('./sdoutMaxTime.txt','file'))
        tmp=importdata('./sdoutMaxTime.txt');  
        maxtime=tmp(1)/1e-9;              % 1st element is maxtime, 2nd element is time when pt crosses detector
        clear tmp;
    else
        maxtime = round(10e-3 /1e-9); %absolute maxtime usually wont be reached;
    end
    
% reduced mass of particles
mass = 32; % enter in au  % total mass for O2
% defined such that 1 au mass = 1.660538782e-27 kg
dt = 100;  % time step, (!)UPDATE RampUpFunc when dt~=50, using GenerateRampUpFunc()
nozzleDistance   = 0.1048 - 0.001; %-1mm is correction to Ncoil*cH   % enter in meters
DetectorDistance = 0.2817 - 0.001 ;  % enter in meters %-1mm is correction to Ncoil*cH
DetectorRadius   = 2e-3 ; % enter in meters, no conversion unit of *100 %was 1e-3
raleyLength = 3e-3; %laser shape % was 1.5e-3
sdDelay = 0; %artificial delay %enter in seconds
targetVel = 430; %enter in m/s
nozzleOpeningTime = 60e-6; %enter in seconds
sdPhase = 45;                       % 90 ... 0 degree
numberofcoils =15;

cSeq=[];
cSeq=cSeq-mod(cSeq,dt);  %t jumps from 3,6,9,etc without this it will never stop on cSeq=1000, for example     

cSpacing = ones(1,numberofcoils-1)*eps;  %enter in meters     % coil spacing, 5mm, (which is already included instrinsicly 
tempcheck = diff(cSpacing);
if ~all(tempcheck == tempcheck(1)), disp('warning, cSpacing need to be equal, see 11 in GetForce.m'); end;

Ncoil = size(cSpacing,2)+1; % number of coils

% How fast ramping up and ramping down
if (firstRunFlag==0)
NstepRampUp   = ceil(53.5e-6/ 1e-9 / dt); 
elseif (firstRunFlag==1)
NstepRampUp   = ceil(1e-9/ 1e-9 / dt); % instantaneous for determining cseq 
end

NstepRampDown = ceil(8e-6 / 1e-9 /dt); % exponential decay of 7us, approximated by 7us linear decay
NstepB4RampDown = ceil(73.5e-6/ 1e-9 / dt); %not used currently, need to MANUALly change line 60
% uncomment if we want to update this function
GenerateRampUpFunc(dt,NstepRampUp,NstepB4RampDown) %update this if dt changes

%%%%%%%%%%%%%%%%%%%%
[cR,cDim,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB]=LoadField();
%%%%%%%%%%%%%%%%%%%%

deceleratorLegnth = Ncoil*cDim(3); %for reference %in meters
totalLength = nozzleDistance + deceleratorLegnth + DetectorDistance;
totalLength

%%%%%%%%%%%%%%%%%%%%
[xy,Vxy]=InitializeParticle(Npt,nozzleDistance,nozzleOpeningTime,targetVel,cR);
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
[targetPtArrivalTime,maxtime,cSeq]=FlyParticle(firstRunFlag,sdParam,currentFolder,maxtime,cR,Npt,NstepRampUp,NstepRampDown,mass,DetectorDistance,DetectorRadius,raleyLength,nozzleOpeningTime,sdPhase,dt,Ncoil,cSpacing,xy,Vxy,cDim,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,cSeq);
%%%%%%%%%%%%%%%%%%%%

if ~(exist('coilSequence.txt','file'))
    cSeq(1,:) = cSeq(2,:) - 73.5e-6/1e-9;   % Turning ON the coils 40 microseconds before it should be turned OFF
    cSeq=cSeq-mod(cSeq,dt);   % rounding to integer au
    cSeq(1,cSeq(1,:) < 0) = 0;  % no negative timing
    
    fpa=fopen([currentFolder '/coilSequence.txt'],'w');
    for count2 = 1:size(cSeq,1)
        for count1 = 1:Ncoil
            try 
                fprintf(fpa,'%g ',cSeq(count2,count1)+ sdDelay/1e-9-mod(sdDelay/1e-9,dt));       % errors if we dont use all the coils, as in case of H
                                                              % see GetcSeq.m line 6
            catch
                fprintf(fpa,'NaN ');
            end
        end
    fprintf(fpa,'\n\n');
    end
    fclose(fpa);
end %endif

fpm = fopen('sdoutMaxTime.txt','w');
fprintf(fpm,'%+e\t%+e\n',maxtime*1e-9,targetPtArrivalTime*1e-9);
fclose(fpm);

cSeq
disp('main DONE!!');

end % end main



