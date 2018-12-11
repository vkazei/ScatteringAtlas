%% resolution for Cijs
% this function evaluates spectral sensitivities for all Cij parameters and
% plots resolution patterns
%
% it takes as input the type of parameterization
% parSET is defined in mainCIJ.m

function resFunc(parSET)

% make folder for pictures if not made yet
mkdir(parSET.path_pattern_save);

%%read parameters from structure
v2struct(parSET);

fprintf('\n Starting scattering resolution analysis for');
disp(parSET.WTCellArray);

% set of parameters to be useed for description of the model
if CijFlag == 0
    % parameterization of Oh-Alkhalifah
    parNameArray = {'\rho', 'V_p', 'V_s', '\epsilon_1', '\epsilon_d', '\eta_1', '\eta_d', '\delta_3', '\gamma_1', '\gamma_d'};
else
    % parameterization with elastic constants
    parNameArray = {'\rho', 'C_{11}', 'C_{22}', 'C_{33}', 'C_{12}', 'C_{13}', 'C_{23}', 'C_{44}', 'C_{55}', 'C_{66}'};
end

disp('Parameterization');
disp(parNameArray);

%partial derivatives computed in Maple are loaded
% number of Cij parameters
par_maple.NCijPar = tNumPar-1;
%CijFlag = 1 - parameterize by Cij 
par_maple.CijFlag = CijFlag;
%drop out Vs related parameters = 0
par_maple.VsFlag = VsFlag;
%include density = 1
par_maple.denFlag = denFlag;
% fix isotropic parameters = 1
par_maple.isoKnownFlag = isoKnownFlag;
% save pictures to the same location
par_maple.path_pattern_save = parSET.path_pattern_save;

parInCijTensor=loadFromMapleDen('parInCij.mat', par_maple);

iFig = 11;

%WTCellArray = {'PP','PSV','SVSV','SHSH','SVSH'};
%%
for cellWT = WTCellArray
    % item is a 1x1 cell array, not the actual string contents.
    % item{1} is the string contents.
    %disp(item{1})
    WT = cellWT{1};
    
    %%
    % for normalization
    maxAbs=0;
    fig1=figure(1);
    
    % sensitivity array
    TsensAll.(WT) = zeros(tNumPar,nPhi,nKz);
    
    disp('evaluating the sensitivities for orthorhombic parameters');
    for iPar=1:tNumPar
        % partial derivatives for parameter number iPar
        Cij(:,:) = parInCijTensor(iPar,:,:);
        % drawing diffraction-based radiation patterns incidence from top
        %drawPatternTA(WT,Cij,denFlag,'diffZ');
        
        % drawing reflection-based radiation patterns horizontal reflaector
        
        % evaluate the sensitivities - this is inefficient, but we don't care
        % for now
        % it would be better to compute the sensitivities only once and then
        % compose sensitivities for parameters from those of Cij
        Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, denFlag);
        denFlag = 0;
        % correction for density
        
        subaxis(1,tNumPar,iPar,'Spacing',0.005,'Margin',0.1);
        imagesc(0:5:90,KzMin:dKz:KzMax,Tsens');
        colormap(mycmap);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        axis xy
        axis on
        maxAbs = max(max(max(abs(Tsens))),maxAbs);
        TsensAll.(WT)(iPar,:,:)=Tsens;
        % title is the respective parameter
        title(parNameArray{iPar});
    end;
    
    %% find similar patterns and put them into the frames of the same color
    colorArr = {'g','c','m','y','r','b','k','w'};
    cNum=1;
    for iPar=1:tNumPar
        subaxis(1,tNumPar,iPar,'Spacing',0.005,'Margin',0.1);
        caxis([-maxAbs maxAbs])
        set(gca,'FontSize',20)
        Tsens1(:,:) = TsensAll.(WT)(iPar,:,:);
        for jPar = iPar+1:tNumPar
            Tsens2(:,:) = TsensAll.(WT)(jPar,:,:);
            % check the angle between patterns as linear spaces and that
            % they are non-zero
            if (subspace(Tsens1(:),Tsens2(:))<pi/18) ...
                    && norm(Tsens1)*norm(Tsens2) > 10^-5
                
                colorVec = colorArr{cNum};
                
                subaxis(1,tNumPar,iPar,'Spacing',0.005,'Margin',0.1);
                ax = gca;
                ax.XColor = colorVec;
                ax.YColor = colorVec;
                ax.LineWidth = 5;
                
                subaxis(1,tNumPar,jPar,'Spacing',0.005,'Margin',0.1);
                ax = gca;
                ax.XColor = colorVec;
                ax.YColor = colorVec;
                ax.LineWidth = 5;
                drawnow
                cNum = cNum+1;
            end
        end
    end
    
    fig1 = gcf;
    fig1.PaperPosition = [0 0 10 1];
    print(strcat(path_pattern_save,WT,'nPar',num2str(CijFlag)),'-depsc','-r0')
    %%
    
    %perturbing Cij parameters
    
    Cij = zeros(6);
    TsensCijAll.(WT) = zeros(6,6,nPhi,nKz);
    fig2 = figure(2);
    disp('evaluating sensitivities for all 21 Cij parameters');
    for i = 1:6
        for j = i:6
            Cij = zeros(6);
            Cij(i,j) = 1;
            Cij(j,i) = 1;
            Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, denFlag);
            TsensCijAll.(WT)(i,j,:,:) = Tsens;
            TsensCijAll.(WT)(j,i,:,:) = Tsens;
            subaxis(6,6,(i-1)*6+j,'Spacing',0.01,'Margin',0.01);
            imagesc(0:5:90,KzMin:dKz:KzMax,Tsens'/max(10^-5+abs(Tsens(:))));
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            caxis ([-1 1])
            colormap(mycmap);
            axis xy
            ylim([0 KzMax])
        end
    end
    hold on
    %% find similar patterns for the Cij parameters
    cNum = 1;
    % determining angles
    for i1=1:6
        for j1=i1:6
            for i2=i1:6
                for j2=i2:6
                    Tsens1(:,:) = TsensCijAll.(WT)(i1,j1,:,:);
                    Tsens2(:,:) = TsensCijAll.(WT)(i2,j2,:,:);
                    if (subspace(Tsens1(:),Tsens2(:))<pi/18) ...
                            && ~isequal([i1 j1],[i2 j2]) && ~isequal([i1 j1],[j2 i2]) ...
                            && norm(Tsens1)*norm(Tsens2) > 10^-10
                        
                        colorVec = colorArr{cNum};
                        
                        subaxis(6,6,(i1-1)*6+j1,'Spacing',0.01,'Margin',0.01);
                        ax = gca;
                        ax.XColor = colorVec;
                        ax.YColor = colorVec;
                        ax.LineWidth = 4;
                        
                        subaxis(6,6,(i2-1)*6+j2,'Spacing',0.01,'Margin',0.01);
                        ax = gca;
                        ax.XColor = colorVec;
                        ax.YColor = colorVec;
                        ax.LineWidth = 4;
                        drawnow
                        cNum = cNum+1;
                    end
                end
            end
        end
    end
    %%
    
    
    fig2 = gcf;
    fig2.PaperPositionMode = 'auto';
    print_N_note([path_pattern_save,WT,'Cij'])
    
    %% density picture separate
    
    Cij = zeros(6);
    % setting denFlag to 1
    Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, 1);
    subaxis(6,6,26,'Spacing',0.001,'Margin',0.01);
    imagesc(0:5:180*phiMax/pi,KzMin:dKz:KzMax,Tsens');
    title Density
    set(gca,'xtick',[0 180])
    set(gca,'xticklabel',[0 180])
    set(gca,'ytick',[0 2*kappa 2])
    set(gca,'yticklabel',{'0','2\omega/V_p', '2\omega/V_s'})
    %set(gca,'YAxisLocation', 'right')
    set(gca,'FontSize',15)
    set(gca,'TickLength',[0.1 0.1])
    
    ylabel('K_z');
    xlabel('Azimuth(^o)');
    caxis ([-1 1])
    colormap(mycmap);
    %colorbar southoutside
    axis xy
    
    fig2 = gcf;
    fig2.PaperPositionMode = 'auto';
    print_N_note([path_pattern_save,WT,'density'])
    
    
%     %% SVD analysis for "time domain"
%           
%     fig2=figure(iFig);
    superSens.(WT) = reshape(TsensAll.(WT),tNumPar,nPhi*nKz);
     [~,S.(WT),V.(WT)] = svd(superSens.(WT)','econ');
TsensTotal.par = [TsensTotal.par superSens.(WT)];
%     
%     subplot(1,2,1);
%     
%     imagesc(log10(S.(WT)(1:min(tNumPar,nKz*nPhi),1:min(tNumPar,nKz*nPhi))/max(max(S.(WT)))));
%     caxis ([-2 0])
%     title(strcat('New param, log_{10}(singular values), ',WT));
%     colormap('parula');
%     colorbar
%     
%     subplot(1,2,2);
%     
%     imagesc(V.(WT));
%     caxis ([-1 1])
%     title(strcat('Singular vectors, ',WT));
%     colormap(mycmap);
%     colorbar
%     iFig = iFig+1;
    %%
    
    
    for SVDthresh = SVDthreshArr
        figure(351);
        %SS1 = S.(WT)(1:10,1:10)/(S.(WT)(1:10,1:10)+eye(10)*10^-1*S.(WT)(1,1))
        thresH =  SVDthresh* S.(WT)(1,1);
        SS1 = (S.(WT)>thresH);
        SS1 = SS1(1:10,1:10);
        %imagesc((V.(WT))'*V.(WT));
        caxis ([-1 1])
        
        imagesc(V.(WT)*SS1*(V.(WT))');
        caxis ([-1 1])
        set(gca,'YTickLabel',parNameArray)
        set(gca,'XTickLabel',parNameArray)
        title(strcat('Hierarchical param. resolution matrix'));
        colormap(mycmap);
        colorbar
        set(gca,'YTick',1:10)
        set(gca,'XTick',1:10)
        if CijFlag==1
            title(strcat('C_{ij}, \rho param. resolution matrix'));
        end
        set(gca,'FontSize',20)
        axis equal
        axis tight
        grid on
        colorbar
        
        fig2 = gcf;
        fig2.PaperPosition = [0 0 10 10];
        print(strcat(path_pattern_save,WT,'_resMatrix',num2str(CijFlag),'sTol',num2str(SVDthresh*100)),'-depsc2','-r0');
    end
    
    %%
    superSensCij.(WT) = reshape(TsensCijAll.(WT),36,nPhi*nKz);
    %[Uij.(WT),Sij.(WT),Vij.(WT)] = svd(superSensCij.(WT)');
    [~,Sij.(WT),Vij.(WT)] = svd(superSensCij.(WT)','econ');
    TsensTotal.ij = [TsensTotal.ij superSensCij.(WT)];
     
    iFig = iFig+1;
  
    
end

%% resolution of all data types together

[~, STotalPar, VTotalPar] = svd((TsensTotal.par)','econ');
[~, STotalCij, VTotalCij] = svd((TsensTotal.ij)','econ');

for SVDthresh = SVDthreshArr
    figure(351);
    %SS1 = S.(WT)(1:10,1:10)/(S.(WT)(1:10,1:10)+eye(10)*10^-1*S.(WT)(1,1))
    thresH =  SVDthresh* STotalPar(1,1);
    SS1 = (STotalPar>thresH);
    SS1 = SS1(1:10,1:10);
    %imagesc((V.(WT))'*V.(WT));
    caxis ([-1 1])
    
    imagesc(VTotalPar*SS1*(VTotalPar)');
    caxis ([-1 1])
    set(gca,'YTickLabel',parNameArray)
    set(gca,'XTickLabel',parNameArray)
    title(strcat('Hierarchical param. total resolution matrix'));
    colormap(mycmap);
    colorbar
    set(gca,'YTick',1:10)
    set(gca,'XTick',1:10)
    if CijFlag==1
        title(strcat('C_{ij}, \rho param. resolution matrix'));
    end
    set(gca,'FontSize',20)
    axis equal
    axis tight
    grid on
    colorbar
    
    fig2 = gcf;
    fig2.PaperPosition = [0 0 10 10];
    print(strcat(path_pattern_save,'resMatrix',num2str(SVDthresh*100),'Total',num2str(CijFlag)),'-depsc2','-r0');
    
end

%%
figure(iFig);
iFig = iFig+1;
%subplot(2,1,1);

STdiag = diag((STotalPar(1:min(tNumPar,nKz*nPhi),1:min(tNumPar,nKz*nPhi))/...
    max(max(STotalPar)))+10^-30);

semilogy(STdiag,'b','LineWidth',4);
%caxis ([-3 0])
title(strcat('Singular values'));
%title(strcat('Singular values ',WT,' scattering'));
set(gca,'FontSize',20)
axis tight
axis([1 tNumPar 10^-3 1])
grid on
hold on

if CijFlag
    cijSVD = STdiag;
    save(strcat('cijSVD',WT,'.mat'),'cijSVD');
else
    load(strcat('cijSVD',WT,'.mat'));
    semilogy(cijSVD,'r','LineWidth',4);
    
end

fig2 = gcf;
fig2.PaperPosition = [0 0 10 10];
print_N_note([path_pattern_save,'TotalSingVal']);

%%
figure(iFig);
iFig = iFig+1;

imagesc(VTotalPar);
caxis ([-1 1])
title(strcat('Singular vectors'));
set(gca,'YTickLabel',parNameArray)
colormap(mycmap);
colorbar
set(gca,'YTick',1:10)
set(gca,'XTick',1:10)
if CijFlag==1
    set(gca,'YTickLabel',parNameArray)
    title(strcat('C_{ij}, \rho param., singular vectors'));
end
set(gca,'FontSize',20)
axis equal
axis tight
grid on
colorbar

fig2 = gcf;
fig2.PaperPosition = [0 0 10 10];

print_N_note([path_pattern_save,'TotalSingVec',num2str(CijFlag)]);

end

