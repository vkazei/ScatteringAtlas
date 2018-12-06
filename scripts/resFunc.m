%% simple resolution for Cijs

function res = resFunc(parSET)

%read parameters from structure
v2struct(parSET);

if CijFlag == 0
    parNameArray = {'\rho', 'V_p', 'V_s', '\epsilon_1', '\epsilon_d', '\eta_1', '\eta_d', '\delta_3', '\gamma_1', '\gamma_d'};
else
    parNameArray = {'\rho', 'C_{11}', 'C_{22}', 'C_{33}', 'C_{12}', 'C_{13}', 'C_{23}', 'C_{44}', 'C_{55}', 'C_{66}'};
end;

%partial derivatives computed in Maple
parInCijTensor=loadFromMapleDen('parInCij.mat',tNumPar-1,CijFlag,VsFlag,denFlag,isoKnownFlag);

iFig = 11;

%WTCellArray = {'PP','PSV','SVSV','SHSH','SVSH'};
%%
for cellWT = WTCellArray
    % item is a 1x1 cell array, not the actual string contents.
    % item{1} is the string contents.
    %disp(item{1})
    WT = cellWT{1}
    
    %%
    % for normalization
    maxAbs=0;
    fig1=figure(1);
    
    % sensitivity array
    TsensAll.(WT) = zeros(tNumPar,nPhi,nKz);
    
    % evaluating the sensitivities
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
    
    %% drawing specials
    colorArr = {'g','c','m','y','r','b','k','w'};
    cNum=1;
    for iPar=1:tNumPar
        subaxis(1,tNumPar,iPar,'Spacing',0.005,'Margin',0.1);
        caxis([-maxAbs maxAbs])
        set(gca,'FontSize',20)
        Tsens1(:,:) = TsensAll.(WT)(iPar,:,:);
        for jPar = iPar+1:tNumPar
            Tsens2(:,:) = TsensAll.(WT)(jPar,:,:);
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
    print(strcat('FIG/',WT,'nPar',num2str(CijFlag)),'-depsc','-r0')
    %%
    
    %perturbing Cij parameters
    
    Cij = zeros(6);
    TsensCijAll.(WT) = zeros(6,6,nPhi,nKz);
    fig2 = figure(2);
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
%%   
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
    print(strcat('FIG/',WT,'Cij'),'-depsc','-r0')
    
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
    print(strcat('FIG/',WT,'density'),'-depsc','-r0')
    
    
    
    %% SVD analysis for "time domain"
    
    
    
    fig2=figure(iFig)
    superSens.(WT) = reshape(TsensAll.(WT),tNumPar,nPhi*nKz);
    [~,S.(WT),V.(WT)] = svd(superSens.(WT)','econ');
    TsensTotal.par = [TsensTotal.par superSens.(WT)];
    
    subplot(1,2,1);
    
    imagesc(log10(S.(WT)(1:min(tNumPar,nKz*nPhi),1:min(tNumPar,nKz*nPhi))/max(max(S.(WT)))));
    caxis ([-2 0])
    title(strcat('New param, log_{10}(singular values), ',WT));
    colormap('parula');
    colorbar
    
    subplot(1,2,2);
    
    imagesc(V.(WT));
    caxis ([-1 1])
    title(strcat('Singular vectors, ',WT));
    colormap(mycmap);
    colorbar
    iFig = iFig+1;
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
        print(strcat('FIG/',WT,'_resMatrix',num2str(CijFlag),'sTol',num2str(SVDthresh*100)),'-depsc2','-r0');
    end;
    
    %%
    superSensCij.(WT) = reshape(TsensCijAll.(WT),36,nPhi*nKz);
    %[Uij.(WT),Sij.(WT),Vij.(WT)] = svd(superSensCij.(WT)');
    [~,Sij.(WT),Vij.(WT)] = svd(superSensCij.(WT)','econ');
    TsensTotal.ij = [TsensTotal.ij superSensCij.(WT)];
    %
    % subplot(2,2,3);
    %
    % imagesc(log10(Sij.(WT)(1:21,1:21)/max(max(Sij.(WT)))));
    % caxis ([-15 0])
    % title(strcat('C_{ij}, log_{10}(singular values), ',WT));
    % colormap('parula');
    % colorbar
    %
    %
    % subplot(2,2,4);
    %
    % imagesc(Vij.(WT)(1:21,1:21));
    % caxis ([-1 1])
    % title(strcat('C_{ij}, singular vectors for ',WT));
    % colormap(mycmap);
    % colorbar
    
    iFig = iFig+1;
    %
    % fig2 = gcf;
    % %fig2.PaperPositionMode = 'auto';
    % fig2.PaperPosition = [0 0 24 8];
    % print(strcat('FIG/',WT,'_SVD4'),'-dpng','-r0')
    % %surf(Tsens);
    % %hold on
    
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
    print(strcat('FIG/resMatrix',num2str(SVDthresh*100),'Total',num2str(CijFlag)),'-depsc2','-r0');
    
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
    
end;
% SVDthreshPlot = SVDthreshArr(1)*ones(size(cijSVD));
% semilogy(SVDthreshPlot,'g','LineWidth',3);
% SVDthreshPlot = SVDthreshArr(2)*ones(size(cijSVD));
% semilogy(SVDthreshPlot,'k','LineWidth',3);
fig2 = gcf;
fig2.PaperPosition = [0 0 10 10];
%print(strcat('FIG/TotalSingVal',num2str(CijFlag)),'-depsc2','-r0')

print('FIG/TotalSingVal','-depsc2','-r0')
%%
%subplot(2,1,2);
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
print(strcat('FIG/TotalSingVec',num2str(CijFlag)),'-depsc2','-r0')

% subplot(2,2,3);
%
% imagesc(log10(STotalCij(1:21,1:21)/max(max(STotalCij))));
% caxis ([-13 0])
% title(strcat('C_{ij}, log_{10}(singular values)'));
% colormap('parula');
% colorbar
%
%
% subplot(2,2,4);
%
% imagesc(VTotalCij(1:21,1:21));
% caxis ([-1 1])
% title(strcat('C_{ij}, singular vectors together'));
% colormap(mycmap);
% colorbar
%
% fig1 = gcf;
% fig1.PaperPositionMode = 'auto';
% print(strcat('FIG/SVD_AllTogether'),'-dpng','-r0')

for SVDthresh = SVDthreshArr
    
    nParRes(TsensAll.(WT),SVDthresh);
    %nParResAngle(TsensAll.(WT),SVDthresh,parSET);
    if strcmp(WT,'PSV') || strcmp(WT,'PSV')
        nParResAnglePS(TsensAll.(WT),SVDthresh,parSET);
    elseif strcmp(WT,'SVSV') || strcmp(WT,'SVSH') || strcmp(WT,'SHSH')
        blablbllbljkl = 100000000
    elseif strcmp(WT,'PP')
        nParResAngle(TsensAll.(WT),SVDthresh,parSET);
    end;
    
end;

%nParRes(TsensAll.PP,0.01);
% nParRes2(TsensAll,0.1);
% nParRes2(TsensAll,0.01);
% nParRes3(TsensAll,0.1);
% nParRes3(TsensAll,0.01);

res = 0;
end

