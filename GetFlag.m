% Successfully implemented the the coil switching input format, but still with instantaneous switching
% eps used to fix error in GetFlag_ver3.m
function [cFlag] = GetFlag(t,Ncoil,NstepRampUp,NstepRampDown,cFlag,cSeq,RpUpPercentArray)
% RampUpFcn   = linspace(0,1,NstepRampUp);        % linearly increasing from 0% to 100%    Can put a custom fcn here
RampUpFcn = RpUpPercentArray; % log increasing, instead of linear
RampDownFcn = linspace(0,1,NstepRampDown+1);      % linearly decreasing from 100% to 0%

persistent rampCounter;           % ERROER if write persistent rampCounter holdcFlag
persistent holdcFlag;             

if (isempty(rampCounter))   % initializes  state
    rampCounter=zeros(1,Ncoil);  
end
if (isempty(holdcFlag))    % initializes  state
    holdcFlag=zeros(1,Ncoil);
end



for k = 1:Ncoil                                 % k denotes which coil we are considering

    %%% reset counter from to 0  when reached NstepRampUp if ramping up OR when reached NstepRampDown if  ramping down
        
        if ( rampCounter(k) == NstepRampUp+1 && abs(holdcFlag(k))  <=  eps   ||   rampCounter(k) == NstepRampDown+1 && abs(holdcFlag(k)-1) <= eps )           %%% NOTE (rampCounter(k) == 1) YIELDS FALSE
            rampCounter(k) = 0;                   %%% REASON was machine calc error rampCounter(k) was 1.0000000000000000, but not 1, hence eps need to be invoked
         end
        %%%  
       
        if(k > size(cSeq,2))        % if some coils, towards the end are not defined, 
                                    % then dont do anthing
            continue;
        end
        
% if no match for the t we are at now, find(x,1) the one means we only need to 
% consider first element to determine if its empty of not     
% AND we are not in ramping mode
    if (isempty(find(cSeq(:,k)==t, 1)) && rampCounter(k) == 0)   
                                                            
                                                            
        continue;                         % skipping, ie no modification to cFlag for this coil
    
    else      % then we initiate ramping
        


            
            if (rampCounter(k)  == 0)
            holdcFlag(k) = cFlag(k);    %hold the value in begining of ramping
            end
            
            
            if (abs(holdcFlag(k))  <=  eps)    %RampUp
                
                
                rampCounter(k) = rampCounter(k)+1; 
                cFlag(k) = 0+RampUpFcn(rampCounter(k));



            elseif (abs(holdcFlag(k)-1) <= eps)    %RampDown
                
                
                rampCounter(k) = rampCounter(k)+1;
                cFlag(k) = 1-RampDownFcn(rampCounter(k));

            end
            %%% end Ramp
        
        
       

    end

        
        %%%% Ramp    



end
