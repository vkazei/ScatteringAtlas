%% This is the main script that should produce all the pictures for the paper and moves them
close all
clearvars
addpath([pwd '/tools']);
addpath([pwd '/maple']);
system('mkdir FIG');

%% radiation patterns in 3D view
% produce patterns for PP scattering
 WT = 'PP';

for i=1:6
    for j=i:6
        drawPatternsIJ(WT,i,j);
    end
end

%%
% define parameters
tNumPar = 10; % total number of parameters that are considered for inversion
mycmap = myColormap; % standard seismic colormap blue-white-red
isoKnownFlag = 0; % 1 - exclude isotropic parameters from inversion

SVDthreshArr = [0.1, 3*10^-1]; % 
CijFlag = 1; % 1 = (Cij, rho) parametrization instead of whatever

kappa = 1/sqrt(3); %vp/vs ratio
VsFlag = 1;
denFlag=1;
KzMin = 0*2/sqrt(2);
KzMax = 2;
dKz = 0.01;
phiMax = pi;
dPhi = pi/72;
nKz = size(KzMin:dKz:KzMax,2);
nPhi = size(0:dPhi:phiMax,2);
TsensTotal.par = zeros(tNumPar,1);
TsensTotal.ij = zeros(36,1);

%parInCijTensor=loadFromMapleDen('parInCij.mat',tNumPar-1,CijFlag,VsFlag,denFlag,isoKnownFlag);
%parInCijTensor = eye(tNumPar);

iFig = 11;

WTCellArray = {'PP'};
%WTCellArray = {'PP','PSV','SVSV','SHSH','SVSH'};

% collect everything into a structure
parSET = v2struct;

% cases of full illumination 

% P-P waves only

%% spectral pattern for C_55

%%  
figure(55);
Cij = zeros(6);
Cij(5,5) = 1;
    Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, 0);
    %subaxis(1,1,1,'Spacing',0.25,'Margin',0.3);    
    imagesc(0:5:180*phiMax/pi,KzMin:dKz:KzMax,Tsens');
    title('Spectral sensitivity to C_5_5');
    set(gca,'xtick',[0 180])
    set(gca,'xticklabel',[0 180])
    %set(gca,'ytick',[])
    set(gca,'FontSize',25)
    
    
    ylabel('K_z/k_0');
    xlabel('Azimuth(^o)');
    caxis ([-1 1])
    colormap(mycmap);
    %colorbar southoutside
    axis xy tight
    
    fig2 = gcf;
    fig2.PaperPosition = [0 0 10 7];
    print(strcat('FIG/',WT,'_C_55'),'-depsc','-r0')
    

%% PP waves


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


