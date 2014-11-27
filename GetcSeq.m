%changing the phase from 90degree to 45degree, also will change cSeq
%initially to be -1, so that now each coil is switched individually

function [cSeq]=GetcSeq(t,Ncoil,Npt,sdPhase,cSpacing,cDim,xy,cSeq,currentFolder)

persistent cSeqCounter;
if (isempty(cSeqCounter))   % initializes  state
    cSeqCounter=ones(Ncoil,1);  
end
persistent cSeqONOFFCounter;
if (isempty(cSeqONOFFCounter))   % initializes  state
    cSeqONOFFCounter=zeros(Ncoil,1);  
end

sdtemp = sdPhase/90 * 0.5;    % 0.5 ... 0 multiple of coil length

% Npt/2 is
% ...
% .x.
% ...
% that middle one in phase space, which is our target molecule

%% scheme1 (OLD VERSION): all ON in start, turn off when target arrives at center of each coil 
% for k = 1:Ncoil
% % if the target pt is ant any coil's first half, and if the coils are off, then turn it ON
% 
% %if ( abs(cSeqONOFFCounter(k)) <= eps)
% if (cSeqONOFFCounter(k) == 0 )
%    %if (xy(1,round(Npt/2)) >= cH*(k-1) + sum(cSpacing(1:k-1))) && (xy(1,round(Npt/2)) < cH*(k-0.5) + sum(cSpacing(1:k-1))) %turning ON when pt enters coil 
%    if ~(xy(1,round(Npt/2)) >= cH*(k-1+sdtemp) + sum(cSpacing(1:k-1))) && (xy(1,round(Npt/2)) < cH*(k-0.0) + sum(cSpacing(1:k-1))) %ON from the very start
%    cSeq(cSeqCounter(k),k)=t;  
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=1;
% %    break; %bug to add break?
%     end
% 
% %elseif (abs(cSeqONOFFCounter(k)-1) <=eps) 
% elseif (cSeqONOFFCounter(k) == 1 )
% % if the target pt is ant any coil's last half, and if the coils are off, then turn it OFF
%     if (xy(1,round(Npt/2)) >= cH*(k-1+sdtemp) + sum(cSpacing(1:k-1))) && (xy(1,round(Npt/2)) < cH*(k-0.0) + sum(cSpacing(1:k-1)))  %turning OFF
%    cSeq(cSeqCounter(k),k)=t;   % was   cSeq(cSeqCounter,:)=t+dt;
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=0;
% %    break; %bug to add break?
%     end
% 
% end %endif
% 
% end %endfor

%% scheme2 (BETTER WRITTEN): all ON, switch OFF like in scheme1 and turning back ON after a delay 

% for k = 1:Ncoil                         %turning all ON initially
% if (cSeqONOFFCounter(k) == 0 && t == 0)
%    
%     cSeq(cSeqCounter(k),k)=t+dt;
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=1;
% 
% end
% end
%     
%     
%     
% for k = 1:Ncoil                      % when to turn OFF
% if (cSeqONOFFCounter(k) == 1  && xy(1,round(Npt/2)) >= cH*(k-1+sdtemp) + sum(cSpacing(1:k-1))) && (xy(1,round(Npt/2)) < cH*(k-0.0) + sum(cSpacing(1:k-1)) ...
%     && cSeqCounter(k) == 2) %enforce to be on 2nd line
%    
%     cSeq(cSeqCounter(k),k)=t;
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=0;           % turned OFF
%     
% end
% end    
% 
% 
% 
% 
% for k = 1:Ncoil                         %turning all ON again
%     
%     
% % enforcing onto the 3rd line of cSeq
% if (cSeqONOFFCounter(k) == 0    &&   cSeqCounter(k) == 3  )  % HARDCODED  %after a coil has beeng turned
%                                                              % ON and OFF exactly once then do following 
% 
%    cSeq(cSeqCounter(k),k)=cSeq(cSeqCounter(k)-1,k)+10000;     % switch ON after 10 000 timestep
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=1;               %turned ON
% 
% 
% end
% end

%% scheme3:
% if (exist('./coilSequence.txt','file'))       %if this file exist
%     
%     if (isempty(cSeq)==1) %if cSeq doesnt exist, then import, otherwise do nothing
% 
%     cSeq=importdata('./coilSequence.txt');
%     cSeq(1,:) = cSeq(2,:) - 80e-6/1e-9;   % Turning ON the coils 40 microseconds before it should be turned OFF
%     cSeq=cSeq-mod(cSeq,dt);   % rounding to integer au
% 
%     end
% 
% else       %when file dont exist, then we build
% 
% 
% 
% sdmode = 1; %switching 1 coils at a time
% 
% for k = 1:Ncoil                         %turning all ON initially
% if (cSeqONOFFCounter(k) == 0 && t == 0)
%    
%     cSeq(cSeqCounter(k),k)=t;
%    
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=1;
% 
% end
% end
%     
%     
%     
% for masterCounter = 1:sdmode:Ncoil    % 1, 4, 7 etc...                       % when to turn OFF
% k = masterCounter+sdmode-1;           % make an offset to do 3,6,9 etc ...
% if k > Ncoil
%     continue;            % do nothing if out of bound
% end
% 
% if (cSeqONOFFCounter(k) == 1  && xy(1,round(Npt/2)) >= cH*(k-1+sdtemp) + sum(cSpacing(1:k-1))) && (xy(1,round(Npt/2)) < cH*(k-0.0) + sum(cSpacing(1:k-1)) ...
%     && cSeqCounter(k) == 2) %enforce to be on 2nd line
%    
%    for tempk =  k :-1: k-sdmode+1       %set all previous coils at the same time
%    cSeq(cSeqCounter(tempk),tempk)=t;
%    cSeqCounter(tempk)=cSeqCounter(tempk)+1;
%    cSeqONOFFCounter(tempk)=0;           % turned OFF
%    end %endfor
%    
% end %endif
% 
% end %endfor
%     
% 
% 
% end %end if exist


%% scheme3 simplified
if (exist('coilSequence.txt','file'))       %if this file exist
    
    if (isempty(cSeq)==1) %if cSeq doesnt exist, then import, otherwise do nothing
    cSeq=importdata('coilSequence.txt');
    end
else       %when file dont exist, then we build

    for k = 1:Ncoil                         
        if (cSeqONOFFCounter(k) == 0 && t == 0) %turning all ON initially
           cSeq(cSeqCounter(k),k)=t;
           cSeqCounter(k)=cSeqCounter(k)+1;
           cSeqONOFFCounter(k)=1;
        end
    end

    
    
for k = 1:Ncoil    % 1, 4, 7 etc...                       % when to turn OFF
    if k > Ncoil
        continue;            % do nothing if out of bound
    end
    
    
if (cSeqONOFFCounter(k) == 1  && xy(3,round(Npt/2)) >= cDim(3)*(k-1+sdtemp) + sum(cSpacing(1:k-1))) && (xy(3,round(Npt/2)) < cDim(3)*(k-0.0) + sum(cSpacing(1:k-1)) ...
    && cSeqCounter(k) == 2) %enforce to be on 2nd line
   % IF 5th coil AND 21th coil is broken:
%     if ((k == 5 || k == 21 ) &&  cSeqONOFFCounter(k) == 1 && cSeqCounter(k) == 2  )             % do nothing if coils broken
%    cSeq(cSeqCounter(k),k)=0;
%    cSeqCounter(k)=cSeqCounter(k)+1;
%    cSeqONOFFCounter(k)=0;
%     continue;            % do nothing if coils broken
%     end

   cSeq(cSeqCounter(k),k)=t;
   cSeqCounter(k)=cSeqCounter(k)+1;
   cSeqONOFFCounter(k)=0;           % turned OFF
   
   
end %endif

end %endfor
    


end %end if exist








end %endfunction