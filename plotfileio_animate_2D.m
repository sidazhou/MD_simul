clear all; close all; clc;
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
Npt = 1000; % if dont succeed, then Npt is set to 1000 by default
clear temp count; fclose(fid);
end
%%%%%%%%%%%%%%%%%%
fid=fopen('./sdPlotThis.txt','r');
fseek(fid,-86,'eof'); %set to last line
tmp= fscanf(fid,'%e',6);
maxtime = tmp(2); % autoset maxtime
clear tmp; fclose(fid);

%loading field
if ~exist('cH','var')
if(ispc()) %if PC
load 'C:\Users\sidz\Desktop\3mm_600A_5.2T_3D.mat'
else %on westgrid, which returns ''
load '~/3mm_600A_5.2T_3D.mat'
end
end

fp = fopen('./sdPlotThis.txt','r');
% X,Y,Z are from 3mm_600A_5.2T_3D.mat
cH = range(Z);   % region height            
if range(X) == range(Y)  % if this is round, then we can define a Radius
cR = range(X)/2 ; %=range(Y)/2
else
disp('error, cR not defined');
cR = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = importdata('coilSequence.txt');
numberofcoils = size(c,2);
cSpacing = ones(1,numberofcoils-1)*eps;      %MANUAL from main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Plotting
fig = figure('units','normalized','OuterPosition',[0.005, 0.005, 0.99, 0.99]); 
fig2 = subplot(4,3,[2 6]); %phase space
fig1 = subplot(4,3,[7 9]); %position space
fig3 = subplot(4,3,[10 12]); %Coils ON/OFF
fig4 = subplot(4,3,4);   % pseudo TOF

uicontrol(... % Button for updating selected plot
    'Parent', fig, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.00 0.95 0.05 0.05],...
    'String','Pause Simulation',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Callback', 'pause()');

%initiating coil on/off
subplot(fig3); 
for ii=1:numberofcoils  
hLine(ii)=line([   ((ii-1)*cH+sum(cSpacing(1:ii-1)))        (ii*cH+sum(cSpacing(1:ii-1)))],[0 0]);  %initializing for fig3

if mod(ii,2)==0, set(hLine(ii),'color',[0 0 0.5]); end % alternating color

set(hLine(ii),'LineWidth',45);
set(hLine(ii),'Visible','off');
end
%%%%%%%%%%%%%
indCSeq=ones(numberofcoils,1);  % 1st element of array of cSeq
flagONOFF = zeros(numberofcoils,1);  %everything is off
%%%%%%%%%%%%%

a = zeros(Npt,8);
indicator=1;




%% START OF MAIN LOOP
while(indicator)
    
try               %MOST IDIOTIC WAY TO DETERMINE EOF, THAT IS, TO BREAK WHEN ERROR READING THE FILE, OMFG
        for i = 1:Npt
% n,t*1e-9,xy(1,n),xy(2,n),xy(3,n),Vxy(1,n)/1e-9,Vxy(2,n)/1e-9,Vxy(3,n)/1e-9
            a(i,:)= fscanf(fp,'%e',8);
        end
        Ind= ~(a(:,5) <= 0 & sqrt(a(:,3).^2 + a(:,4).^2) >= cR); %all molecules which are not hitting the wall

        
catch err          
   indicator=0; %WHOLE INDICATOR THING IS TO STOP THE HANDLES FROM DELETING THEMSELVES,
   disp('sdERROR');
end

%% subplot position space
subplot(fig1); hold on
axis([0 (cH*numberofcoils+sum(cSpacing(1:numberofcoils-1))) -cR cR]);

if ~exist('xEdge','var')
xEdgeN=200; %Manual
xEdge=linspace(0,cH*numberofcoils+sum(cSpacing(1:numberofcoils-1)),xEdgeN);
yEdgeN=50; %Manual
yEdge=linspace(-cR,cR,yEdgeN);
end

[temp,xEdgep,yEdgep]=hist2(a(Ind,5),a(Ind,4),xEdge,yEdge); %putting all pts into 2D histogram

h1=pcolor(xEdgep,yEdgep,temp); %plotting
set(h1,'LineStyle','none'); 
colormap(flipud(gray)); %flipped gray colormap
caxis([0 Npt/xEdgeN/yEdgeN*8]);  %manual *8 is from trial and error
colorbar;
clear temp;

xlabel('Z (m)'); ylabel('Y (m)');

% h5=plot(a(a(:,1)==round(Npt/2),5),a(a(:,1)==round(Npt/2),4),'r.'); hold on;
% set(h5,'MarkerSize',15); % doesnt work


%% subplot phase space
subplot(fig2); hold on

if ~exist('xEdge2','var')
% get axis from scatter plot
temp=figure('visible','off');
plot(a(Ind,5)-mean(a(Ind,5)),a(Ind,8),'k.'); 
sdaxislimFor2=axis(gca);
close(temp);

% edges
xEdge2N=100;  %Manual
xEdge2=linspace(sdaxislimFor2(1),sdaxislimFor2(2),xEdge2N);
yEdge2N=100;  %Manual
yEdge2=linspace(sdaxislimFor2(3),sdaxislimFor2(4),yEdge2N);

end

axis(sdaxislimFor2); %set axis for 2d hist
[temp2,xEdge2p,yEdge2p]=hist2(a(Ind,5)-mean(a(Ind,5)),a(Ind,8) ,xEdge2,yEdge2);

h2=pcolor(xEdge2p,yEdge2p,temp2);
set(h2,'LineStyle','none'); 
colormap(flipud(bone)); %flipped gray colormap
caxis([0 Npt/xEdge2N/yEdge2N*12]);  %manual *8 is from trial and error
colorbar;
clear temp2;

 xlabel('Z (m)'); ylabel('Vz (m/s)')

% h2=plot(a(:,5)-mean(a(:,5)),a(:,8),'k.'); xlabel('Z (m)'); ylabel('Vz (m/s)');
% h4=plot(a(a(:,1)==round(Npt/2),5)-mean(a(:,5)),a(a(:,1)==round(Npt/2),8),'r.');

% h2=plot(a(:,4),a(:,5),'k.'); xlabel('y (m)'); ylabel('Vx (m/s)');hold on;
% h4=plot(a(round(Npt/2),4),a(round(Npt/2),5),'r.'); xlim([0 3e-3]);hold on;

%% pseudo TOF

subplot(fig4); hold on

temp4=histc(a(Ind,5)-mean(a(Ind,5)),xEdge2);
h4=plot(xEdge2,temp4);
 xlabel('Z (m)'); ylabel('intensity (au)')

%% time
h3=annotation(fig,'textbox',...
    [0.0755 0.8178 0.0926 0.0408],...
    'String',{'time(s)' a(1,2) 'time(au)' a(1,2)/1e-9},... %time
    'FitBoxToText','on');
% 
%% subplot coil on/off

    subplot(fig3); 
    axis([0 (cH*numberofcoils+sum(cSpacing(1:numberofcoils-1))) -cR cR]);
    
    for ii=1:numberofcoils
    if indCSeq(ii)> size(c,1) 
        continue;
    else
        if(a(1,2)/1e-9 >= c(indCSeq(ii),ii))


            if(flagONOFF(ii)==0)
            flagONOFF(ii)=1;
            set(hLine(ii),'Visible','on');
            elseif(flagONOFF(ii)==1)
            flagONOFF(ii)=0;
            set(hLine(ii),'Visible','off');
            end


            indCSeq(ii) = indCSeq(ii)+1; %move array pointer to next element
        end
    end
    end
%%
if (indicator~=0)

    
 refresh; drawnow;  
%  pause(.2);

delete(h4);
delete(h3);
delete(h2);
delete(h1);


end


end

fclose(fp);
disp(err);
disp('plotfileio DONE!!!');

