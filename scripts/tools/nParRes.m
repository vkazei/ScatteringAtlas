% this function computes the number of singular values larger than tol of Tsens matrix-Tensor as a function of minKz maxKz 

function numParRes = nParRes(Tsens, sTol)

%we presume that Tsens has the sizes as zeros(9,nPhi,nKz);

nKz = size(Tsens,3);
nPhi = size(Tsens,2);
nParam = size(Tsens,1);

nPar = zeros(nKz, nPhi);
maxSingVal = zeros(nKz, nPhi);

ijTrig = 1;

myEps = 10^-30;

for i=1:nKz
    i
    for j=1:nPhi
        % we crop out the Kz values between iKz and jKz
        reshIJKz = reshape(Tsens(:,j:nPhi,i:nKz),nParam,(nPhi+1-j)*(nKz+1-i));
        S=svd(reshIJKz,'econ');
        if ijTrig == 1
            MAXS = abs(S(1,1));
            ijTrig = 0;
        end
        normFact = 1;%sqrt(MAXS/S(1,1));%sqrt((nPhi+1-j)*(nKz+1-i)/(nPhi*nKz));
        nPar(i,j) = nnz(abs(S(:))>sTol*(S(1,1)*normFact+sTol*myEps));
        maxSingVal(i,j) = S(1,1);        
    end;
end;
%%
figure();
A = colormap(jet);
imagesc(90:-0.5:0,0:0.2:2,maxSingVal);
colorbar
axis xy
figure(777)
imagesc(180:-0.5:0,0:0.2:2,nPar);
axis xy
colorbar
shading interp
fig2 = gcf;

fig2.PaperPosition = [0 0 5 4];
xlabel('Max azimuth(^o)');
%ylabel('Min normalized wavenumber (2(v_s/v_p)cos\theta)'); 
ylabel('2cos(\theta_{max}/2)'); 
set(gca,'FontSize',20)
set(gca,'XTick',[0:45:180])
colormap(A(4:6:64,:))
caxis([0 10]);
cbh=colorbar('v');
set(cbh,'YTick',[0:1:10])
print(strcat('FIG/nPar_sTol',num2str(sTol),'.eps'),'-depsc2','-r0');



end


     