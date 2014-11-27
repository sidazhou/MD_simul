function [cR,cDim,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB]=LoadField()
%% Importing data 3D

if(ispc()) %if PC
load 'C:\Users\sidz\Desktop\3mm_600A_4.83T_3D_SmoothEdgeNaN.mat'
else %on westgrid, which returns ''
load '~/3mm_600A_4.83T_3D_SmoothEdgeNaN.mat'
end



cDim(1) =  max(X) - min(X) ;  % X
cDim(2) =  max(Y) - min(Y) ;  % Y
cDim(3) =  max(Z) - min(Z) ;  % Z

if cDim(1) == cDim(1)  % if this is round, then we can define a Radius
cR = cDim(1)/2 ; 
else
disp('error, cR not defined');
cR = [];
end

indexedB=(4/4.83)*(600/600)*Bnorm; % or =sqrt(Br.^2 + Bz.^2)

% X and Y dimension may need to the swapped (due to matlab column indexing), not sure where tho.
XX=XX; %#ok<ASGSL,NODEF>
YY=YY; %#ok<ASGSL,NODEF>
ZZ=ZZ-min(min(min(ZZ)));  %#ok<NODEF>  %make Z axis from 0...2*range, instead of -range ... +range

gTN = size(indexedB); %get number of grid points

gSize(1) = cDim(1) / (gTN(1)-1); % X  % grid size
gSize(2) = cDim(2) / (gTN(2)-1); % Y
gSize(3) = cDim(3) / (gTN(3)-1); % Z

[indexeddBbdX indexeddBbdY indexeddBbdZ]= gradient(indexedB,gSize(1),gSize(2),gSize(3));
end %end function


%% visualization
% data=dBzdZ;
% 
% hold on;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% 
% sdzlim = size(data,3);
% 
% w=vol3d('cdata',data(:,:, 1:2),'texture','2D');  % make figure
% zlim([0 sdzlim]);
% view(3);
% colorbar;
% %  alphamap('rampdown');
% alphamap('decrease');
% 
% for i=1:100:sdzlim
% % q=slice(data,[],[],i);  set(q,'LineStyle','none');
% 
% if ~(i>sdzlim-1)
%     w=vol3d('cdata',data(:,:, 1:i+1),'texture','2D'); 
% end
% 
% pause(0.1)
% 
% % delete(q);
% if ~(i>sdzlim-2)
%     delete(w.handles); %for vol3d
% end
% 
% end
% % 
% 

