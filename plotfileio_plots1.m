function []=plotfileio_plots1(sdParam,currentFolder,DetectorRadius,raleyLength,cDim,cR,Npt,sdPhase,cSpacing)
%% Initializations
%{
 sdParam = 6;
 DetectorDistance=0.2600;
 cDim(3) = 0.012; 
 cR = 0.0015;
 Npt = 10;
 sdPhase=60; 
 raleyLength = 3e-3;
 DetectorRadius=1e-3;
 currentFolder = pwd;
%}
cH = cDim(3);

folderNumber = GetfolderNumber(currentFolder); %subfunctions

%%
%% *d.jpg
%% 
maxtime=importdata([currentFolder '/sdoutMaxTime.txt']); % reading [maxtime targetptArrivingtime]
maxtime = maxtime(1); %getting maxtime from [maxtime targetptArrivingtime]
maxtime = maxtime * 0.9; % lets plot 90% of the maxtime, everything should be pass the detector by now

currentBlock = zeros(Npt,8);
XPosSummary = zeros(51,110); % 51 from sdEdge, 110 i think we record only 100snapshots max
tIndex=1;

sdEdges = -5e-3:0.2e-3:5e-3; % from -5mm to +5mm, in x position, 1cm total, size=51 

fid=fopen([currentFolder '/sdPlotThis' num2str(sdParam) '.txt'],'r');

while(1)    
    for n=1:Npt
    currentBlock(n,:) = fscanf(fid,'%e',8);
    end
    
    % removing particles in the wall
    ind = currentBlock(:,5) <= 0 & currentBlock(:,3) >= cR & currentBlock(:,4) >= cR;
    currentBlock(ind,:) = []; 
    
    
    temp = histc(currentBlock(:,3),sdEdges); % hist the x position
    XPosSummary(:,tIndex) = temp;
    TargetPtZpos(tIndex) = currentBlock(currentBlock(:,1)==round(Npt/2),5);
    tIndex = tIndex + 1;
    
    if currentBlock(1,2) >= maxtime
        break;
    end
end

% plotting
xruler=sdEdges;
yruler=TargetPtZpos;

figure;
plothandle = pcolor(yruler,xruler,XPosSummary(1:size(xruler,2),1:size(yruler,2))); % x and y mixed up, dont ask why
set(plothandle,'LineStyle','none'); 
colormap(flipud(gray)); %flipped gray colormap
colorbar;
% caxis([0 2000]);

xlabel('z(m)');
ylabel('x(m)');


saveas(gcf, [currentFolder '/' num2str(folderNumber) num2str(sdParam) 'd'], 'jpg')
saveas(gcf, [currentFolder '/'  num2str(folderNumber) num2str(sdParam) 'd'], 'fig')


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



end % end main


