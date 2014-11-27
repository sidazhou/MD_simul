function []=REALmain(currentFolderLocal,currentFolderRemote)
                           % clear all; clears persistent variables, but also clears breakpoints
                           % clear; clears variables, but dont clear breakpoint (good)
                           % but also dont clear persistent variables(bad)
                           % (BUT APPARANTLY clear all does clear persistent variables 
tic % recording running time
% string manipulate currentFolderLocal to get a nice tree


if(ispc())
    currentFolder = pwd;
else
    currentFolder=currentFolderRemote;

eval(['!mkdir -p "' currentFolder '"']);  % make dir
cd(currentFolder);                % cd into that remote directory
end

if(ispc()) %if PC
Npt = 1;    % on PC, 1 pt    
else
Npt = 100000;  % running on cluster 100000 particle
end


firstRunFlag=0;
parfor sdParam=1:9
[DetectorRadius,raleyLength,cDim,cR,sdPhase,cSpacing]=main(firstRunFlag,sdParam,currentFolder,Npt);
plotfileio_plots(sdParam,currentFolder,DetectorRadius,raleyLength,cDim,cR,Npt,sdPhase,cSpacing);
% plotfileio_plots1(sdParam,currentFolder,DetectorRadius,raleyLength,cDim,cR,Npt,sdPhase,cSpacing);

delete([currentFolder '/sdPlotThis' num2str(sdParam) '.txt'])
disp('sdPlotThis.txt deleted');
disp(['plotfileio_plots DONE!! for sdParam= ' num2str(sdParam)]);

end

%% reprocessing of coilSequence.txt 
a=importdata('coilSequence.txt');
b=a';
c = b(:,2) * 1e-9;

% Auto-detect folders number
% temp = max(findstr(currentFolder,'/')); % the pointer position of the last folder
% folderNumber = currentFolder(temp+1:temp+2);  % extract 2 digits after the \
% clear temp;

 fp = fopen([currentFolder '/cSeq.txt'],'wt');

 for i = 1: size(c)
   fprintf(fp,'%.9f\r\n',c(i));
 end
 fclose(fp);

totalExecutionTime = toc; 
save([currentFolder '/totalExecutionTime.txt'],'totalExecutionTime', '-ascii');

 
 
end %end REALmain