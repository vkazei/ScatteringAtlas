% this function computes the number of singular values larger than tol of Tsens matrix-Tensor as a function of minKz maxKz 

function numParRes = nParResAnglePS(Tsens, sTol, parSET)


%% we presume that Tsens has the sizes as zeros(9,nPhi,nKz);

nKz = size(Tsens,3);
nPhi = size(Tsens,2);
nParam = size(Tsens,1);
kappa = parSET.kappa;


ijTrig = 1;

myEps = 10^-30;
thetaMax = deg2rad(90);
thetaMin = deg2rad(0);
iTheta05=thetaMax:-pi/180:thetaMin;
nPar = zeros(size(iTheta05,2), nPhi);
nParHigh = nPar;

maxSingVal = zeros(nKz, nPhi);
for iTh=1:size(iTheta05,2)
    i = max(1,round((cos(iTheta05(iTh))*kappa+...
        (1-(kappa^2)*sin(thetaMin)^2)^0.5)*nKz/2))
    
    iMax = round((cos(thetaMin)*kappa+...
        (1-(kappa^2)*sin(thetaMin)^2)^0.5)*nKz/2)
    
    iMin = max(1,round((cos(thetaMax)*kappa+...
        (1-(kappa^2)*sin(thetaMax)^2)^0.5)*nKz/2))
    
    i,nKz
    for j=1:nPhi
        % we crop out the Kz values between iKz and jKz
        
        % high wavenumbers
        reshIJKzHigh = reshape(Tsens(:,j:nPhi,i:iMax),nParam,(nPhi+1-j)*(iMax+1-i));
        % low wavenumbers
        reshIJKz = reshape(Tsens(:,j:nPhi,iMin:i),nParam,(nPhi+1-j)*(i-iMin+1));
        
        % full azimuth -- does not work yet
        %reshIJKzFullA = reshape(Tsens(:,1:nPhi,iMin:iMax),nParam,nPhi*(i-iMin+1)); 
        
        
        S=svd(reshIJKz);
        SHigh=svd(reshIJKzHigh);
        if ijTrig == 1
            MAXS = abs(S(1,1));
            ijTrig = 0;
        end
        normFact = 1;%sqrt(MAXS/S(1,1));%sqrt((nPhi+1-j)*(nKz+1-i)/(nPhi*nKz));
        nPar(iTh,j) = nnz(abs(S(:))>sTol*(S(1,1)*normFact+sTol*myEps));
        nParHigh(iTh,j) = nnz(abs(SHigh(:))>sTol*(SHigh(1,1)*normFact+sTol*myEps));
        maxSingVal(iTh,j) = S(1,1);        
    end;
end;

figure(555);
A = colormap(jet);
imagesc(90:-0.5:0,90:-5:0,maxSingVal);
colorbar
axis xy
figure(333)
imagesc(180:-0.5:0,rad2deg(iTheta05),nPar);
axis ij
colorbar
shading interp
fig2 = gcf;
fig2.PaperPosition = [0 0 5 4];
xlabel('Max azimuth(^o)');
%ylabel('Min normalized wavenumber (2(v_s/v_p)cos\theta)'); 
ylabel('\theta_{min}/2'); 
set(gca,'FontSize',20)
set(gca,'XTick',[0:45:180])
set(gca,'YTick',rad2deg(iTheta05(end:-10:1)))
colormap(A(4:6:64,:))
caxis([0 10]);
cbh=colorbar('v');
set(cbh,'YTick',[0:1:10])
print(strcat('FIG/nParPS',num2str(parSET.CijFlag),'_sTol',num2str(sTol),'.eps'),'-depsc2','-r0');

%%

figure(334)
imagesc(180:-0.5:0,rad2deg(iTheta05),nParHigh);
axis ij
colorbar
shading interp
fig2 = gcf;
fig2.PaperPosition = [0 0 5 4];
xlabel('Max azimuth(^o)');
%ylabel('Min normalized wavenumber (2(v_s/v_p)cos\theta)'); 
ylabel('\theta_{max}/2'); 
set(gca,'FontSize',20)
set(gca,'XTick',[0:45:180])
set(gca,'YTick',rad2deg(iTheta05(end:-10:1)))
colormap(A(4:6:64,:))
caxis([0 10]);
cbh=colorbar('v');
set(cbh,'YTick',[0:1:10])
print(strcat('FIG/nParHighPS',num2str(parSET.CijFlag),'_sTol',num2str(sTol),'.eps'),'-depsc2','-r0');


numParRes = 0;
end


     