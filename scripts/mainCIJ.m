%% This is the main script that should produce all the pictures for the paper and moves them
%
close all
clearvars
addpath([pwd '/tools']);
addpath([pwd '/maple']);
addpath([pwd '/core']);
system('mkdir FIG');

%% define parameters
tNumPar = 10; % total number of parameters that are considered for inversion
mycmap = myColormap; % standard seismic colormap blue-white-red
isoKnownFlag = 0; % 1 - exclude isotropic parameters from inversion

SVDthreshArr = [0.1, 3*10^-1]; %
CijFlag = 1; % 1 = (Cij, rho) parametrization instead of whatever

kappa = 1/sqrt(3); %vp/vs ratio
VsFlag = 1; % 0 -- pseudo-acoustic approximation is used,
% 1 -- full set of elastic parmaeters with C_44, C_55, C_66
denFlag=1;  % 0 -- exclude density from inversion
% 1 -- density is included
KzMin = 0*2/sqrt(2); % minimum normalized wavenumber for inversion
KzMax = 2; % maximum normalized wavenumber (effective angle concept)
dKz = 0.01; % sampling in wavenumber - smaller numbers refine pictures but reduce performance
phiMax = pi; % maximum azimuth available
dPhi = pi/72; % sampling in azimuth
iFig = 11;


%% automatic initiation
nKz = size(KzMin:dKz:KzMax,2); % number of wavenumbers for illumination
nPhi = size(0:dPhi:phiMax,2); % number of azimuths

TsensTotal.par = zeros(tNumPar,1);
TsensTotal.ij = zeros(36,1);

path_pattern_save = '../latex/Fig/patterns/';

% collect everything into a structure
parSET = v2struct;


% %% radiation patterns in 3D view
% % produce patterns for PP scattering
%
% % wave type incident wave - scattered wave
WT = 'PP';
% mkdir(path_pattern_save);
% for i=1:6
%     for j=i:6
%         drawPatternsIJ(WT,i,j,path_pattern_save);
%     end
% end


%% spectral pattern for C_55 P-P wave scattering

figure(55);
Cij = zeros(6);
Cij(5,5) = 1;
Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, 0);
imagesc(0:5:180*phiMax/pi,KzMin:dKz:KzMax,Tsens');
title('Spectral sensitivity to C_5_5');
set(gca,'xtick',[0 180])
set(gca,'xticklabel',[0 180])
set(gca,'FontSize',25)

ylabel('K_z/k_0');
xlabel('Azimuth(^o)');
caxis ([-1 1])
colormap(mycmap);
axis xy tight

fig2 = gcf;
fig2.PaperPosition = [0 0 10 7];
print(strcat('FIG/',WT,'_C_55'),'-depsc','-r0')


%% PP waves
parSET.WTCellArray = {'PP'};

% Cij parameterization
resFunc(parSET);

% new parameterization
parSET.CijFlag = 0;
resFunc(parSET);

system('rm -rf ../latex/Fig/PP_Full')
system('mv FIG ../latex/Fig/PP_Full')


%% P-SV waves only
system('mkdir FIG');
close all
% Cij parameterization

parSET.WTCellArray = {'PSV'};
parSET.CijFlag = 1;
parSET.VsFlag = 1;
parSET.denFlag=1;

resFunc(parSET);

parSET.CijFlag = 0;
resFunc(parSET);

system('rm -rf ../latex/Fig/PSV')
system('mv FIG ../latex/Fig/PSV')
system('mkdir FIG');

%% P-SH waves only
system('mkdir FIG');
close all
% Cij parameterization

parSET.WTCellArray = {'PSH'};
parSET.CijFlag = 1;

resFunc(parSET);

parSET.CijFlag = 0;
resFunc(parSET);

system('rm -rf ../latex/Fig/PSH')
system('mv FIG ../latex/Fig/PSH')
system('mkdir FIG');


%% P-P,SV waves together
close all

system('mkdir FIG');
parSET.WTCellArray = {'PP','PSV'};

parSET.CijFlag = 1;
resFunc(parSET);

%% P-P,SV,SH waves together
close all

system('mkdir FIG');
parSET.WTCellArray = {'PP','PSV','PSH'};

parSET.CijFlag = 1;
resFunc(parSET);

system('rm -rf ../latex/Fig/PP_PSV_PSH')
system('mv FIG ../latex/Fig/PP_PSV_PSH')

%% SV-SV waves
close all

system('mkdir FIG');
parSET.WTCellArray = {'SVSV'};

parSET.CijFlag = 1;
resFunc(parSET);

parSET.CijFlag = 0;
resFunc(parSET);


system('rm -rf ../latex/Fig/SVSV')
system('mv FIG ../latex/Fig/SVSV')

%% SH-SH waves
close all

system('mkdir FIG');
parSET.WTCellArray = {'SHSH'};

parSET.CijFlag = 1;
resFunc(parSET);

parSET.CijFlag = 0;
resFunc(parSET);

system('rm -rf ../latex/Fig/SHSH')
system('mv FIG ../latex/Fig/SHSH')

%% SV-SH waves
close all

system('mkdir FIG');
parSET.WTCellArray = {'SVSH'};

parSET.CijFlag = 1;
resFunc(parSET);

parSET.CijFlag = 0;
resFunc(parSET);

system('rm -rf ../latex/Fig/SVSH')
system('mv FIG ../latex/Fig/SVSH')


