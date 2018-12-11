% this script computes the resolution pattern - the sensitivity to vertical
% wavenumbers for a given wavetype and Cij perturbation


function Tsens = simpleRes(Cij, WT, KzMin, KzMax, dKz, phiMax, dPhi, denFlag)

% k = vs/vp ratio
k=1/sqrt(3); % Poisson ratio = 0.25

%the goal is to plot the radiation pattern in (Kv, phi) coordinates
[TKz, Tphi] = meshgrid(KzMin:dKz:KzMax, 0:dPhi:phiMax);
nKz = size(TKz,2);
nPhi = size(Tphi,1);

Tsens = zeros(nPhi,nKz);

for iKz=1:nKz
    Kz = KzMin+dKz * (iKz-1);
    for iPhi=1:nPhi
        phi = dPhi * (iPhi-1);
        if strcmp(WT,'PSH') || strcmp(WT,'PSV')
            if (k+1>=Kz) && (1-k<=Kz)
                kRoot = sqrt(2*k^2*Kz^2-k^4+2*k^2-Kz^4+2*Kz^2-1);
                sxy = kRoot/(2*k*Kz);
                sz = (k^2+Kz^2-1)/(2*k*Kz);
                sv = [sxy*cos(phi); sxy*sin(phi); sz];
                gxy = kRoot/(2*Kz);
                gz = (Kz^2-k^2+1)/(2*Kz);
                gv = [gxy*cos(phi); gxy*sin(phi); gz];
            else
                sv = [0; 0; 0];
                gv = [0; 0; 0];
            end
        elseif strcmp(WT,'SHSH') || strcmp(WT,'SVSV') || strcmp(WT,'SVSH')
            % the scale of the diagram is defined by the S waves
            sqKz2 = sqrt(1-(Kz/2)^2);
            sv = [sqKz2*cos(phi); sqKz2*sin(phi); Kz/2];
            gv = [-sv(1); -sv(2); Kz/2];
        elseif strcmp(WT,'PP')
            if (Kz > 2*k)
                sv = [0; 0; 0];
                gv = [0; 0; 0];
            else
                sqKzP2 = sqrt(1-(Kz/(2*k))^2);
                sv = [sqKzP2*cos(phi); sqKzP2*sin(phi); Kz/(2*k)];
                gv = [-sv(1); -sv(2); Kz/(2*k)];
            end
        end
        
        if denFlag==1
            Tsens(iPhi,iKz) = RsgDenWT(sv,gv,WT);
            %denFlag=0;
        end
        cijkl = MS_cij2cijkl(Cij);
        Tsens(iPhi,iKz) = Tsens(iPhi,iKz)+RsgCijklWT(sv,gv,cijkl,WT);
        
    end
end
end






