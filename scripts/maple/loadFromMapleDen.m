%load the tensor of coefficients for parameters from .mat file generated in
%Maple (see Cij.mw to modify)

function parInCijTensor = loadFromMapleDen(fileName, par_maple)

v2struct(par_maple);

if CijFlag == 0
    parNameArray = {'\rho', 'V_p', 'V_s', '\epsilon_1', '\epsilon_d', '\eta_1', '\eta_d', '\delta_3', '\gamma_1', '\gamma_d'};
else
    parNameArray = {'\rho', 'C_{11}', 'C_{22}', 'C_{33}', 'C_{12}', 'C_{13}', 'C_{23}', 'C_{44}', 'C_{55}', 'C_{66}'};
end

CijNameArray = {'\rho', 'C_{11}', 'C_{22}', 'C_{33}', 'C_{12}', 'C_{13}', 'C_{23}', 'C_{44}', 'C_{55}', 'C_{66}'};

mycmap = myColormap;
load(fileName)
%NCijPar=9;

if CijFlag==1
    parInCij = eye(10);
end

figure(100);
subaxis(2,1,1);
% we can normalize later
%parInCij = normc(parInCij);
imagesc(parInCij);
set(gca,'xticklabel',parNameArray)
set(gca,'yticklabel',CijNameArray)
title('Full parameterization');

if VsFlag == 0
    parInCij(8:10,:)=0; % C_44,55,66 are not updated
    if CijFlag~=1
        parInCij(:,9:10)=0; % gamma_1, gamma_d are not updated
        parInCij(:,3)=0; %   Vs is not updated
    end
end
% for the C_ij parametrization
% Cij related perturbations
if denFlag==0
    parInCij(:,1) = 0;
end
if isoKnownFlag
    parInCij(:,1:3) = 0;
    parInCij(1:3,:) = 0;
end
subaxis(2,1,2);
imagesc(parInCij);
set(gca,'xticklabel',parNameArray)
set(gca,'yticklabel',CijNameArray)
title('After muting');
fig1 = gcf;
fig1.PaperPosition = [0 0 10 20];
print_N_note([path_pattern_save,'parDerivCij1D',num2str(CijFlag)]);

parInCijMatrix = zeros(NCijPar+1, NCijPar);
%first parameter is rho there
parInCijMatrix(:,:) = parInCij(2:NCijPar+1,1:NCijPar+1)';
%imagesc(parInCijMatrix)
% map each column of the matrix into the ij notation
iParVector = zeros(1,NCijPar);
parInCijTensor = zeros(NCijPar+1,6,6);
figure(166);
for iPar=1:NCijPar+1
    iParVector(:) = parInCijMatrix(iPar,:);
    %subplot(1,NCijPar+1,iPar);
    subaxis(1,NCijPar+1,iPar,'Spacing',0.005,'Margin',0.1);
    plot(iParVector);
    iParMatrix = CijFromOrthV(iParVector);
    imagesc(iParMatrix);
    colormap(mycmap);
    shading faceted
    axis equal tight ij
    grid on
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    %title(parNameArray{iPar});
    caxis([-2 2]);
    parInCijTensor(iPar,:,:) = iParMatrix(:,:);
end
fig1 = gcf;
fig1.PaperPosition = [0 0 10 1];
print_N_note([path_pattern_save,'parDerivCij']);

end