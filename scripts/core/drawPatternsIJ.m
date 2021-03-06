%% plots radiation pattern  for a perturbation of C_ij parameter
% inputs:
% WT - incident - scattered waves e.g. PP, PSV, SVSH,
% i, j define the perturbed cell in the stiffness matrix
% outputs:
% none
% 

function drawPatternsIJ(WT, i, j, path_pattern_save)

% initiate Cij matrix with zeros
Cij = zeros(6);

% perturb ij cell
Cij(i,j) = 1;
Cij(j,i) = 1;

% types of patterns to plot, see drawPattern for specifications
patternArray = {'diffZ', 'diffX', 'diffXZ', 'diffXYZ', 'trans', 'reflZ'};%

% plot radiation patterns of each type

path_save_pattern_Cij = [path_pattern_save 'C_' num2str(i) num2str(j) '/'];
mkdir(path_save_pattern_Cij);

for patType = patternArray
    drawPattern(WT,Cij,0,patType{1});
    fig1 = gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperSize = [20 10];
    fig1.PaperPosition = [15 5 19 19];
    
    print_N_note([path_save_pattern_Cij,WT,'_',patType{1}]);
    
%     fprintf('saving figure %s \n',[path_save_pattern_Cij,WT,'_',patType{1}])
%     
%     print(fig1,[path_save_pattern_Cij,WT,'_',patType{1}],'-depsc','-r0')
    
end
end