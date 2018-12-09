% this script is for benchmarking the radiation patterns for Cij
% it draws the radiation pattern of given type in 3D
% the path is currently hardcoded to "patterns/"

function drawPattern(WT,Cij,denFlag,patternType)

% PARAMETERS
%% Cij is a 6x6 perturbation of Voigt matrix
% for example C_11 perturbation results in 
% Cij =
% 
%      1     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%      0     0     0     0     0     0
%

%WT = 'SHSH';

%% denFlag == 1, makes a perturbation in density
%  denFlag == 0, no perturbation

%% patternType encodes the type of scattering that is of interest
% 'diffZ' - for the vertical incidence (along 3rd axis)
% 'reflZ' - for reflection from a horizontal reflector
%  other types of scattering can be deduced from the code chunk below
%
%    switch patternType
%             case 'diffZ'
%                 sv = [0 0 1];
%             case 'diffX'
%                 sv = [1 0 0];
%             case 'diffY'
%                 sv = [0 1 0];
%             case 'reflZ'
%                 if strcmp(WT,'PSH') || strcmp(WT,'PSV')
%                     sxy = norm(gv(1:2));
%                     if k >= sxy
%                         sz = sign(gv(3))*sqrt(k^2 - sxy^2);
%                         sv = [-gv(1), -gv(2), sz];                        
%                     else
%                         sv = [0; 0; 0];                        
%                     end
%                 else
%                     sv = [-x(i,j);-y(i,j);z(i,j)];
%                 end
%             case 'reflXZ'
%                 if strcmp(WT,'PSH') || strcmp(WT,'PSV')
%                     sxy = norm(gv(1:2));
%                     if k >= sxy
%                         sz = sign(gv(3))*sqrt(k^2 - sxy^2);
%                         sv = [-gv(1), -gv(2), sz];                        
%                     else
%                         sv = [0; 0; 0];                        
%                     end
%                 else
%                     sv = [-x(i,j);-y(i,j);z(i,j)];
%                 end
%                 sv = rotx(45)*sv;
%                 gv = rotx(45)*gv;
%             case 'trans'
%                 sv = -gv;
%             case 'diffXYZ'
%                 sv = (1/sqrt(3))*[1;1;1];
%             case 'diffXZ'
%                 sv = [1/(2^0.5) 0 1/(2^0.5)];
%             case 'obliq'
%                 %sv = [1/2 0  sqrt(3)/2];
%                 sv = [0.5; 0.5; 0.5^0.5];

n=21;
%
%
% Cij = zeros(6);
% 
% lambda or vp perturbation
% Cij(1:3,1:3)=1;
% 
% ith parameter perturbation 
% Cij(:,:) = parInCijTensor(2,:,:);
% 
%mu pertubation
% for i=1:6
%     Cij(i,i)=1;
% end
k=1/sqrt(3); % Poisson ratio = 0.25
%
[x,y,z] = sphere(n-1);
mesh(x,y,z);
AP = zeros(n);
for i=1:n
    for j=1:n
        gv = [x(i,j);y(i,j);z(i,j)];
        switch patternType
            case 'diffZ'
                sv = [0 0 1];
            case 'diffX'
                sv = [1 0 0];
            case 'diffY'
                sv = [0 1 0];
            case 'reflZ'
                if strcmp(WT,'PSH') || strcmp(WT,'PSV')
                    sxy = norm(gv(1:2));
                    if k >= sxy
                        sz = sign(gv(3))*sqrt(k^2 - sxy^2);
                        sv = [-gv(1), -gv(2), sz];                        
                    else
                        sv = [0; 0; 0];                        
                    end
                else
                    sv = [-x(i,j);-y(i,j);z(i,j)];
                end
            case 'reflXZ'
                if strcmp(WT,'PSH') || strcmp(WT,'PSV')
                    sxy = norm(gv(1:2));
                    if k >= sxy
                        sz = sign(gv(3))*sqrt(k^2 - sxy^2);
                        sv = [-gv(1), -gv(2), sz];                        
                    else
                        sv = [0; 0; 0];                        
                    end
                else
                    sv = [-x(i,j);-y(i,j);z(i,j)];
                end
                sv = rotx(45)*sv;
                gv = rotx(45)*gv;
            case 'trans'
                sv = -gv;
            case 'diffXYZ'
                sv = (1/sqrt(3))*[1;1;1];
            case 'diffXZ'
                sv = [1/(2^0.5) 0 1/(2^0.5)];
            case 'obliq'
                %sv = [1/2 0  sqrt(3)/2];
                sv = [0.5; 0.5; 0.5^0.5];
        end
        
        
        cijkl = MS_cij2cijkl(Cij);
        AP(i,j) = RsgCijklWT(sv,gv,cijkl,WT);
        
        
        if denFlag
            AP(i,j) = AP(i,j) + RsgDenWT(sv,gv,WT);
        end
        
                    
        
        %[x(i,j),y(i,j),z(i,j)] = eP';% * AP(i,j);
    end;
end;


%%
AP = fillmissing(AP,'nearest');
eP = cat(3,x,y,z);
%eSH = cat(3,-y,x,0)
AP = AP/max(abs(AP(:))+10^-10);
%eP' * M * eP
figure(15);
subaxis(1,1,1,'Spacing',0.000,'Margin',-0,'MarginLeft',0.05,'MarginBottom',0.05,'Padding',0.01);
ePA = eP;
for i=1:3 
  eP(:,:,i) = eP (:,:,i).* abs(AP);
end;

m=surf(eP(:,:,1),eP(:,:,2),eP(:,:,3),AP,'EdgeColor','black');
%tri = delaunay(x(:),y(:),z(:));
%m = trisurf(tri,eP(:,:,1),eP(:,:,2),eP(:,:,3),AP,'EdgeColor','none');

%set(m,'facecolor','none')
set(m,'facealpha',0.6);
labelSize = 50;
xlabel('\bf{g_x}','FontSize',labelSize);
ylabel('\bf{g_y}','FontSize',labelSize);
zlabel('\bf{g_z}','FontSize',labelSize);

%% drawing

load('mycmap.colormap','mycmap','-mat');
colormap(mycmap);


%axis ij
%axis off
l = light('Position',[0 10 10]);
set(gca,'CameraPosition',[3 3 2])
set(gca,'CameraTarget',[0 0 0])
lighting phong
shading interp
%colorbar EastOutside
%first view
%camproj('perspective')
%next view
%camproj('orthographic')


arw.length = 1.15;
arw.stemWidth = 0.01;
arw.tipWidth = 0.02;
arw.faceAlpha = 0.5;

h = mArrow3([0 0 0],[0 0 arw.length],'color','black',...
    'stemWidth',arw.stemWidth, 'tipWidth',arw.tipWidth,'facealpha',arw.faceAlpha);
h = mArrow3([0 0 0],[arw.length 0 0],'color','black',...
    'stemWidth',arw.stemWidth, 'tipWidth',arw.tipWidth,'facealpha',arw.faceAlpha);
h = mArrow3([0 0 0],[0 arw.length 0],'color','black',...
    'stemWidth',arw.stemWidth, 'tipWidth',arw.tipWidth,'facealpha',arw.faceAlpha);

set(gcf,'color','w');
caxis([-max(abs(caxis)) max(abs(caxis))])
axis equal tight square

xlim([-1.15 1.15])
ylim([-1.15 1.15])
zlim([-1.15 1.15])

%% axis tight

end