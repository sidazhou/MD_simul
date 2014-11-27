function []=sdsubPar()
%% getting the absolute paths on local machine and on remote machine
% get currentFolderLocal
tic
currentFolderLocal = pwd;
% get currentFolderRemote
  [~,matchend]= regexp(currentFolderLocal,'global.home.pcyzsz1');
  sdPath=currentFolderLocal(matchend+1:end);
  sdPath( sdPath=='\') = '/';  %flipping the folder divider
currentFolderRemote = ['/global/home/pcyzsz1'  sdPath] ;  % currentFolderRemote

% saving path information for retrieving files and cleaning up
save('misc_Info.mat','currentFolderLocal','currentFolderRemote');
%% get coilSequence.txt
if ~exist('coilSequence.txt','file')
Npt = 1;
sdParam = 9;     %MANUAL
firstRunFlag=1;
main(firstRunFlag,sdParam,currentFolderLocal,Npt);       % to get cSeq, and to get maxtime for the real calculation
delete([currentFolderLocal '/sdPlotThis*.txt']); % removing stuff from 
delete([currentFolderLocal '/sdDetectorClickTimes*.txt']); % removing stuff from 
end

toc
%% create parallel sched object
sched = GETCLUSTER();

%{  
%For debugging
sched = GETCLUSTER();
pjob=findJob(sched);
ptask=findTask(pjob)
ptask(1).ErrorMessage
destroy(findJob(sched))
%}

% create the matlab job
pjob=createMatlabPoolJob(sched);

% set the number of CPUs to use
set(pjob, 'MaximumNumberOfWorkers', 10)  % Need number of worker+1, +1 is for master who coordinates the tasks, special to matlabpooljobs
set(pjob, 'MinimumNumberOfWorkers', 10)

%building fileList for set(pjob, 'FileDependencies', fileList)
fileList_tmp= dir('*.m');
fileList = '';
for i = 1:size(fileList_tmp,1)
fileList{i}=fileList_tmp(i).name;
end
fileList = [fileList 'coilSequence.txt']; %add coilSequence.txt to file dependency

set(pjob, 'FileDependencies', fileList);

ptask=createTask(pjob, @REALmain, 0, {currentFolderLocal,currentFolderRemote}); % 0 output arguments, and 2 input arguments

%{ 
%for debugging
get(pjob, 'FileDependencies')
get(pjob, 'JobData')
get(pjob, 'PathDependencies')
get(pjob, 'Tasks')
%}

%submit PBS job
    submit(pjob);

%     waitForState(pjob);
%     results=getAllOutputArguments(pjob)
% 
%     load('misc_Info.mat');
%     eval(['!winscp.exe /console "pcyzsz1@seawolf1.westgrid.ca" /command "synchronize both -filemask=|sdDetectorClickTimes?.txt ""' currentFolderLocal '"" ""' currentFolderRemote '"" " "exit"']);

%     temp = rmdir([currentFolderLocal '/Job*'], 's'); % temp is to supress output error when directory dont exist                                                        
%     clear temp;                                                                                                                                                                                             
%     delete([currentFolderLocal '/Job*.mat']); % 
%     delete([currentFolderLocal '/matlab_metadata.mat']); %                                                                                                                                                                                                

end

function [sched] =  GETCLUSTER()
WestgridID='pcyzsz1';
Email='sidazhou@gmail.com';
Nprocs='10';
Wtime='3:00:00';  % per process
Memory='2gb';  % total memory required by this job

% submitArguments=strcat(' -l procs=',Nprocs,',mem=',Memory,',walltime=',Wtime,',software=MDCE:',Nprocs,' -m bea -M ',Email);
submitArguments=strcat(' -l procs=',Nprocs,',mem=',Memory,',walltime=',Wtime,',software=MDCE:',Nprocs);
VER=version('-release');
switch VER
    case '2011a'
        remoteMatlabRoot='/global/software/matlab-2011a';
    otherwise
        fprintf(' Matlab version %s is not supported\n',VER);
        return;
end
clusterHost='orcinus.westgrid.ca';
remoteDataLocation=strcat('/global/scratch/',WestgridID);
sched = findResource('scheduler','type','generic');
set(sched,'ClusterSize',str2num(Nprocs));
set(sched,'ClusterOsType', 'unix');
set(sched,'HasSharedFilesystem',0);
set(sched,'ClusterMatlabRoot',remoteMatlabRoot);
set(sched,'GetJobStateFcn',@getJobStateFcn);
set(sched,'DestroyJobFcn',@destroyJobFcn);
set(sched,'SubmitFcn',{@distributedSubmitFcn,clusterHost,remoteDataLocation,submitArguments});
set(sched,'ParallelSubmitFcn',{@parallelSubmitFcn,clusterHost,remoteDataLocation,submitArguments});
currentFolder=pwd;
% set(sched, 'DataLocation', [currentFolder '/']);
%get(sched);
end