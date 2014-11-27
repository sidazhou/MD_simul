%% Plotting
% clear all; close all; clc;
disp('starting plotfileio.m');
%% Auto-detect Npt
try
fid=fopen('./sdPlotThis.txt','r');
temp=fgetl(fid); %skip first line, which should have index1
count=1; %initializing
while(1)
temp=fgetl(fid);
temp=str2num(temp);

if (temp(1)==1)     %if index is 1 AGAIN
Npt=count; break;   % ASSIGNING Npt
end
count = count+1;
end
clear temp count; fclose(fid);

catch err
Npt = 100; % if dont succeed, then Npt is set to 1000 by default
clear temp count;
end

%% Plot tragectory


fp = fopen('./sdPlotThis.txt','r');
fig = figure; 
axes; 


indicator=1;
a=zeros(Npt,8);

while(indicator)
    
try               %MOST IDIOTIC WAY TO DETERMINE EOF, THAT IS, TO BREAK WHEN ERROR READING THE FILE, OMFG

for i = 1:Npt
% n,t*1e-9,xy(1,n),xy(2,n),xy(3,n),Vxy(1,n)/1e-9,Vxy(2,n)/1e-9,Vxy(3,n)/1e-9
            a(i,:)= fscanf(fp,'%e',8);
end 
       plot3(a(:,3),a(:,4),a(:,5),'.'); 
        xlim([-1.5e-3 1.5e-3]); %MANUAL
        ylim([-1.5e-3 1.5e-3]);
        zlim([-0.2 1]);
        
        pause(0.5)
catch err          
   indicator=0; %WHOLE INDICATOR THING IS TO STOP THE HANDLES FROM DELETING THEMSELVES, 
end

end




