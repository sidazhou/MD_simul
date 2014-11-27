function [Fxy] = GetForce(sdParam,cR,Npt,Ncoil,cSpacing,cDim,cFlag,XX,YY,ZZ,indexeddBbdX,indexeddBbdY,indexeddBbdZ,indexedB,xy,BfieldSummary,eigenValuesDerivativeSummary)

Fxy = zeros(3,size(xy,2));                             % Memory pre-allocation
dBbdxyz = zeros(3,size(xy,2));  
dUbdB = zeros(3,size(xy,2));  


 % Which coil the pt is in
 % cheking each interval from the start, if pt is in it
    
    ind = xy(3,:) >= (cDim(3)+cSpacing(1)) * Ncoil  |  xy(3,:) < 0;
    ptBoxN = ceil( xy(3,:) / (cDim(3)+cSpacing(1)) ) ;
    ptBoxN(ind) = 0;

    
% if outside coil region, force=0
ind1 = ptBoxN(1,:)==0;     
  Fxy(1,ind1)=0; 
  Fxy(2,ind1)=0;
  Fxy(3,ind1)=0;
  
ind0=isnan(xy(1,:)) | isnan(xy(2,:))| isnan(xy(3,:));

%% interpolate more to the potential at that point

% Following calculates F = - dB/dx * dU/dB
ind2 = ~(ind1 | ind0); % for all rest of the molecules

% x needs the be the position relative to the coil, ie mod
dBbdxyz(1,ind2)= sd_interp3(XX,YY,ZZ,indexeddBbdX,xy(1,ind2),xy(2,ind2),mod(   (xy(3,ind2) - sum(  cSpacing(1:ptBoxN(1,ind2)-1))  ),cDim(3)),'linear',NaN);      % gets dBbdx in SI units
dBbdxyz(2,ind2)= sd_interp3(XX,YY,ZZ,indexeddBbdY,xy(1,ind2),xy(2,ind2),mod(   (xy(3,ind2) - sum(  cSpacing(1:ptBoxN(1,ind2)-1))  ),cDim(3)),'linear',NaN);
dBbdxyz(3,ind2)= sd_interp3(XX,YY,ZZ,indexeddBbdZ,xy(1,ind2),xy(2,ind2),mod(   (xy(3,ind2) - sum(  cSpacing(1:ptBoxN(1,ind2)-1))  ),cDim(3)),'linear',NaN);

    % fit result for the first curve (from bottom), from zeemanO2*.txt
     % f(x) = p(1)*x^3+p(2)*x^2+p(3)*x+p(4);
     
% % FitCoeff=[-0.46741 * 29.9792458e9 ; 0.46741 * 29.9792458e9    ]; %29.9792458e9 makes cm^-1 -> Hz
% FitCoeff=[ 372233660.118598,	-5296246572.37845,	606668984.469696,	7288453640.32503;
% 4238276.28432441,	-81280350.3442507,	-27142976111.1704,	69858025906.0871;
% 238869032.652676,	-2900392051.29761,	-14645481436.2115,	70164835174.3623;
% 119694466.783500,	-280036725.241385,	-5232909900.32078,	71086299214.4007;
% 228503516.181931,	-2759517831.01523,	12397289136.3992,	70176402650.9427;
% 1564243008.72053,	-15514215245.3966,	46462210676.0573,	65839132988.8684;
% -1838916910.21184,	18398776005.7756,	-32409401435.0558,	129862844250.366;
% -492675913.109559,	5565260678.31869,	4623175656.62356,	124968591580.161;
% -237763773.352120,	2750257651.99618,	15147501458.2087,	125805664998.323;];
% 
% % higher level polynomial
% FitCoeff=[[     7.3533e+9,7.8472e+7,-4.4181e+9,-1.515e+8,1.2743e+8,-1.081e+7];
% [	  6.9854e+10,-2.7128e+10,-1.0516e+8,1.154e+7,-9.6985e+5,41887];
% [	  6.9819e+10,-1.322e+10,-4.2597e+9,6.9767e+8,-5.6946e+7,1.6035e+6];
% [	  6.9826e+10,5.1869e+8,-6.5377e+9,2.7239e+9,-4.5857e+8,2.8529e+7];
% [	  6.9824e+10,1.3865e+10,-4.1811e+9,7.1727e+8,-6.4253e+7,2.2007e+6];
% [	  6.9858e+10,2.7098e+10,-5.1817e+7,-5.2108e+7,1.0834e+7,-1.7185e+6];
% [	   1.2616e+11,-1.4357e+10,4.2102e+9,-6.9313e+8,5.6564e+7,-1.5891e+6];
% [	   1.2616e+11,-6.0146e+8,1.0932e+10,-2.578e+9,3.3255e+8,-1.7928e+7];
% [	   1.2615e+11,1.3688e+10,4.1718e+9,-7.5096e+8,7.2325e+7,-3.2916e+6]];
% 
% % p=FitCoeff(sdParam,:);% 1st state (lowest in energy)
B = sd_interp3(XX,YY,ZZ,indexedB,xy(1,ind2),xy(2,ind2),mod(   (xy(3,ind2) - sum(  cSpacing(1:ptBoxN(1,ind2)-1))  ),cDim(3)),'linear',NaN);      % gets dBbdx in SI units

% dUbdB(1,ind2) = (3*p(1)*B.^2 + 2*p(2)*B + p(3) ) * 6.62606957e-34;  % 6.62606957e-34 = plancks constant, turns Hz into Joules
% dUbdB(2,ind2) = (3*p(1)*B.^2 + 2*p(2)*B + p(3) ) * 6.62606957e-34;  % result is Joule/Tesla
% dUbdB(3,ind2) = (3*p(1)*B.^2 + 2*p(2)*B + p(3) ) * 6.62606957e-34;

% higher level polynomial
% dUbdB(1,ind2) =(p(2) +2*p(3)*B +3*p(4)*B.^2 +4*p(5)*B.^3 +5*p(6)*B.^4 ) * 6.62606957e-34;
% dUbdB(2,ind2) =(p(2) +2*p(3)*B +3*p(4)*B.^2 +4*p(5)*B.^3 +5*p(6)*B.^4 ) * 6.62606957e-34;
% dUbdB(3,ind2) =(p(2) +2*p(3)*B +3*p(4)*B.^2 +4*p(5)*B.^3 +5*p(6)*B.^4 ) * 6.62606957e-34;

B = B .* cFlag(ptBoxN(1,ind2)); % scaling the field
% using interpolation
dUbdB(1,ind2) = sd_interp1(BfieldSummary,eigenValuesDerivativeSummary(:,sdParam),B',NaN); % faster version of interp1
dUbdB(2,ind2) = dUbdB(1,ind2); % saving time for not needing to interp1 again.
dUbdB(3,ind2) = dUbdB(1,ind2); % saving time for not needing to interp1 again.


% Force = - dB/dx * dU/dB
Fxy(1,ind2) = - ( cFlag(ptBoxN(1,ind2)) .* dBbdxyz(1,ind2) .* dUbdB(1,ind2) ) * 6.022141290116742e8;  % *6.022141290116742e+08 turns kg*m/s^2 into au
Fxy(2,ind2) = - ( cFlag(ptBoxN(1,ind2)) .* dBbdxyz(2,ind2) .* dUbdB(2,ind2) ) * 6.022141290116742e8;  % defined such that 1 au mass = 1.660538921e-27 kg
Fxy(3,ind2) = - ( cFlag(ptBoxN(1,ind2)) .* dBbdxyz(3,ind2) .* dUbdB(3,ind2) ) * 6.022141290116742e8;  % defined such that 1 au mass = 1.660538921e-27 kg


end  %END GetForce.m

