function []=plotfileio_plots(sdParam,currentFolder,DetectorRadius,raleyLength,cDim,cR,Npt,sdPhase,cSpacing)
%% Initializations
%{
 sdParam = 1;
 DetectorDistance=0.2600;
 cDim(3) = 0.012; 
 cR = 0.0015;
 Npt = 500000;
 sdPhase=60; 
 raleyLength = 3e-3; 
 DetectorRadius=1e-3;
 currentFolder = pwd;
%}
cH = cDim(3);

folderNumber = GetfolderNumber(currentFolder); %subfunctions
% [lastBlock firstBlock Z Vz Zinitial Vzinitial NptFirst NptFinal] = GetLastFirstBlock(currentFolder,Npt); %subfunctions

%%
%% *a.jpg
%% 
%% Plot Histogram, 
a = importdata([currentFolder '/sdDetectorClickTimes' num2str(sdParam) '.txt']);

a = a(abs(a(:,4)) <=DetectorRadius & abs(a(:,3))<=raleyLength,:); % removing stuff thats outside the detector rad, in x, think about the direction of travel

sdBins = 100;
 
[sdHist sdHistPosition] = hist(a(:,2),sdBins);

%%%%%%%%%%%%% convolution of signal
% sdBinWidth = (max(a(:,2)) - min(a(:,2)) )/ sdBins;   % how many sec each bin correspond to
% temp1 = randn(1000,1); % generate 1000 normally distributed numbers
% sdResponseFunc = hist(temp1,round(0.25e-3/sdBinWidth)*2 ); % *2 because FWHM looks better like 0.2ms
% clear temp1;
% sdConvolutedSignal = conv(sdHist,sdResponseFunc);
% plot(sdHistPosition,sdConvolutedSignal(1:sdBins)); % 100 bins
% 
% 
% hold on;
% plot(a(a(:,1)==round(Npt/2),2),0,'r.','MarkerSize',15); % find target molecules position
% xlabel('seconds');
% ylabel('clicks');
% title(strcat('When pt enters detector(fig',folderNumber,'a)'));
% 
% annotation(gcf,'textbox',[0.1903 0.6452 0.2954 0.2214],...
%     'String',{'Npt=',Npt,'Phase=',sdPhase},...
%     'FitBoxToText','on');
% 
% saveas(gcf, strcat(folderNumber,'aConv'), 'jpg')

%%%%%%%%%%%%%%%%%%%%%% non convolved version
figure;
plot(sdHistPosition,sdHist)
hold on;
plot(a(a(:,1)==round(Npt/2),2),0,'r.','MarkerSize',15); % find target molecules position
xlabel('seconds');
ylabel('clicks');
title(strcat('When pt enters detector(fig',folderNumber,'a)'));

annotation(gcf,'textbox',[0.1903 0.6452 0.2954 0.2214],...
    'String',{'Npt=',Npt,'Phase=',sdPhase},...
    'FitBoxToText','on');

saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'a'], 'jpg')
saveas(gcf, [currentFolder '/'  num2str(folderNumber) num2str(sdParam) 'a'], 'fig')


%%
%% *b.jpg
%%
%% Displaying final phase space distribution
% figure1=figure;
% 
% plot(lastBlock(:,5)-mean(lastBlock(:,5)),lastBlock(:,8),'k.'); xlabel('Z (m)'); ylabel('Vz (m/s)');
% 
% title(strcat('final phase space distribution(fig',folderNumber,'b)'));hold on;
% plot(lastBlock(lastBlock(:,1)==round(Npt/2),5)-mean(lastBlock(:,5)),lastBlock(lastBlock(:,1)==round(Npt/2),8),'r.'); hold on;  %weird error on reading the lastblock, indexes are not in right place  

%% scanning the packet around the target molecule, and getting the density
% 
% % limits of the search, eg
% scanRangeZ=0.0025; %+-0.2m
% scanRangeVz=20; %+-25m/s
% 
% index=[]; %init
% 
% 
% % searching which pt are within limit of the search
% for n=1:NptFinal
% if (lastBlock(n,5) >= Z-scanRangeZ && lastBlock(n,5) <= Z+scanRangeZ && ...    %x
%     lastBlock(n,8) >= Vz-scanRangeVz)   %Vx % 223was VX+scanRangeVX
% 
% index = [index lastBlock(n,1)];  %contains all pt index which are within the limits
% 
% end
% end
% 
% % Getting tmp_packet, which is bottom+left+right wall of the box, but not top wall 
% count = 1;
% for n=1:NptFinal
%    
%    if (any(lastBlock(n,1)==index) == 1)    % all the pt which are in the "packet"
%      tmp_packet(count,:) =lastBlock(n,:);
%    count = count + 1;
%    end
%    
% end
% 
% 
% for n=size(index,2) : -1 : 1    %counting backwards as we remove the array named index 
% if(tmp_packet(n,8) > Vz+Vz-min(tmp_packet(:,5)))            %Getting the top wall for the box
%     index(n) = [];  %contains all pt index which are within the limits
% end
% end
% 
% clear count;
% clear tmp_packet;
% %  next, we narrow down the box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% searching for the packet and its initial position
% count = 1;
% for n=1:NptFinal
%    
%    if (any(lastBlock(n,1)==index) == 1)    % all the pt which are in the "packet"
%      packet(count,:) =lastBlock(n,:);
%    count = count + 1;
%    end
%    
% end
% clear count;
% 
% 
% count = 1;
% % for n=1:NptFirst
% for n=1:Npt
%    if (any(firstBlock(n,1)==index) == 1)    % all the pt which are in the "packet"
%      packetInit(count,:) =firstBlock(n,:);
%    count = count + 1;
%    end
%    
% end
% clear count;
% 
% %%
% try
%  boxWidth = scanRangeZ*2;  %x
% boxHeight = scanRangeVz*2;%  boxHeight = max(packet(:,8)) - min(packet(:,8));  %Vx
% 
% rectangle('Position',[Z-scanRangeZ-mean(lastBlock(:,5)),Vz-scanRangeVz,boxWidth,boxHeight],...
%     'EdgeColor','b'); %drawing the box
% 
% 
% % Density =  size(packet,1)/boxWidth  ;  %DONE!! the apprx phase space density
% % was (boxWidth * boxHeight)
% 
% % Create textbox
% annotation(figure1,'textbox',[0.1903 0.6452 0.2954 0.2214],...
%     'String',{'boxWidth',boxWidth,'Density (# of pt/boxwidth)',Density,'Final velocity of taget pt (m/s)',Vz},...
%     'FitBoxToText','on');
% catch
%     disp('could not draw a box');
% annotation(figure1,'textbox',[0.1903 0.6452 0.2954 0.2214],...
%     'String',{'Final velocity of taget pt (m/s)',Vz},...
%     'FitBoxToText','on');
% 
% end

% saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'b'], 'jpg')
% saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'b'], 'fig')











%%
%% *c.jpg
%%

% 
% %% Drawing initial coordinates of the pt that are inside the decelerated packet 
% sdPhase = sdPhase*2*pi/360; %converts deg to rad
% % actual plot
% figure;
% 
% % plot(   (firstBlock(:,3)-Xinitial)/(.5*(spacing+cH))+sdPhase   ,   firstBlock(:,5)-VXinitial   ,'k.');  hold on; %plotting initial distribution, in rads
% % plot(   (packetInit(:,3)-Xinitial)/(.5*(spacing+cH))+sdPhase   ,   packetInit(:,5)-VXinitial   ,'c.','MarkerSize',16);  hold on; %plotting initial distribution, in rads
% % plot(0+sdPhase,0,'r.');
% % xlabel('phase (radians)'); 
% % set(gca,'XTick',-3/2*pi :pi/2: 3/2*pi)
% % set(gca,'XTickLabel',{'-3/2 pi','-pi','-pi/2','0','pi/2','pi','3/2 pi'})
% 
% plot(   (firstBlock(:,5)-Zinitial)   ,   firstBlock(:,8)-Vzinitial   ,'k.');  hold on; %plotting initial distribution, in rads
% plot(   (packetInit(:,5)-Zinitial)   ,   packetInit(:,8)-Vzinitial   ,'m.');  hold on; %plotting initial distribution, in rads
% plot(0,0,'r.');
% xlabel('Z (m)'); 
% 
% % first, center the packet, such that target molecule is at [0 0]
% % thenrelative position is divide by (.5*(spacing+cH)), the distance corresponding to 90 degree 
% % then this is offset by sdPhase
% 
% ylabel('Vz difference from target molecule (m/s)');
% title(strcat('Initial phase space plot of decelerated packet(fig',folderNumber,'c)'));
% % xlim([-3/2*pi 3/2*pi]);
% % ylim([-90 90]);
% annotation(gcf,'textbox',[0.1903 0.6452 0.2954 0.2214],...
%     'String',{'% of pt didnt hit the wall',NptFinal/NptFirst},...
%     'FitBoxToText','on');
% 
% 
% 
% saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'c'], 'jpg')
% saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'c'], 'fig')
% 
% %select what to display on top
% figure(figure1);
% 
% 



%%
%% FUNCTIONS
%%
%% Auto-detect folders number
    function finalStr=GetfolderNumber(currentFolder)

temp = max(findstr(currentFolder,'/')); % the pointer position of the last folder
sdStr = currentFolder(temp+1:size(currentFolder,2));  % extract 2 digits after the \
clear temp currentFolder;

%if .*From00 (assuming from00 is always at the end)
% remove all char after 'F'

[mat idx] = regexp(sdStr, '[Ff]rom', 'match', 'start');
 if ~isempty(idx)  
     sdStr = sdStr(1:idx-1);
 end
clear mat idx;
%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [mat idx] = regexp(sdStr, '\d', 'match', 'start');
 
 %while the last number is a number (ie str2num doesnt returns []), 
%  then look at previous char, see if its a number or a char
 
count=0;
try   %is number is at begining, indexing problem, and errors
 while ( ~isempty( str2num(sdStr( max(idx)-count )  )))   
     count = count + 1;  % count then gives how long the number was, either 00 or 000 etc.
 end
 
catch
disp('nice try..')
end

finalStr = sdStr( max(idx)-count+1 : max(idx));
% % if finalstr is empty, it means there was no numerics in the original string
if isempty(finalStr)
    finalStr = sdStr;
    
end

    end

%% Getting last time block and the final target molecule properties
function [lastBlock firstBlock Z Vz Zinitial Vzinitial NptFirst NptFinal] = GetLastFirstBlock(currentFolder,Npt)

fid=fopen([currentFolder '/sdPlotThis' num2str(sdParam) '.txt'],'r');

%% determing how much bytes is one line:

tmpStart=ftell(fid);
fscanf(fid,'%e',8);
tmpStop=ftell(fid);

if(ispc())
sdLineLength = tmpStop - tmpStart + 2; %LF+CR = 2 characters, when using 'wt' flag in fprintf
else
sdLineLength = tmpStop - tmpStart + 1; %LF+CR = 1 character in linux, dont know why
end

frewind(fid);
clear tmpStart;
clear tmpStop;

%%


lastBlock = zeros(Npt,8); %initialization %last time block

% Getting arrival time of target pt
tmp=importdata([currentFolder '/sdoutMaxTime.txt']);
targetPtArrivalTime=tmp(2); % arrival time of target pt
clear tmp;

while(1)
tmp1= fscanf(fid,'%e',8);
fseek(fid, 2, 'cof'); % move LF+CR 
% ftell(fid)
fseek(fid, sdLineLength*(Npt-1), 'cof'); %skip rest of the block for the same time

    if (tmp1(2)>=targetPtArrivalTime)      %find the starting
        break;
    end
end

% fseek(fid,-sdLineLength,'cof'); %moves back one line


% fwthis = fopen([currentFolder '/LastFirstBlock' num2str(sdParam) '.txt'],'w'); 

% for n=1:Npt
% lastBlock(n,:) = fscanf(fid,'%e',8);
% fprintf(fwthis,'%e ',lastBlock(n,:)); %write LastFirstBlock.txt
% fprintf(fwthis,'\n');
% end
% fprintf(fwthis,'\n'); %separate last from first block


% n t x y z Vx Vy Vz
Z =lastBlock(lastBlock(:,1)==round(Npt/2),5); %target molecule's final z in m
Vz=lastBlock(lastBlock(:,1)==round(Npt/2),8); %target molecule's final Vz in m/s

[lastBlock NptFinal]=RmPtRanIntoWall(lastBlock);

% Getting first time block
frewind(fid);  % set file pointer to start of file
firstBlock = zeros(Npt,8);

for n=1:Npt
firstBlock(n,:) = fscanf(fid,'%e',8);
% fprintf(fwthis,'%e ',firstBlock(n,:));  %write LastFirstBlock.txt
% fprintf(fwthis,'\n');
end
% fclose(fwthis);

Zinitial =firstBlock(firstBlock(:,1)==round(Npt/2),5);
Vzinitial=firstBlock(firstBlock(:,1)==round(Npt/2),8);

% [firstBlock NptFirst]=RmPtRanIntoWall(firstBlock); 
NptFirst=Npt;  %TESTING, should be correct

fclose(fid);
    end


% %% auto-detect the phase
%     function sdPhase = GetsdPhase()
% tmpText=fileread('./main.m');
% tmpIndex=findstr(tmpText,'sdPhase');
% tmpIndex= tmpIndex(1);
% sdPhase = tmpText(tmpIndex+9:tmpIndex+11); % manual
% clear tmpText tmpIndex;
% 
%     end

%% Remove particle that hit the wall
    function [blockName NptNEW]= RmPtRanIntoWall(blockName)
        temp=0;
        NptNEW=Npt;
        for i = Npt:-1:1   %removing particles that hit the wall of coil, see Flyparticle.m
         if( blockName(i,5) <= 0 && sqrt(blockName(i,3)^2 + blockName(i,4)^2) >= cR)
         blockName(i,:) = [];
         temp=temp+1;
         end
        end
        NptNEW = NptNEW - temp;
    end




end % end main


