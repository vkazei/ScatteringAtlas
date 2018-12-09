% compute WT (wavetype) radiation pattern on cijkl as function of sv gv for 3-D
% WT = 'PP' 'SHSH' 'SVSV'
function Rp = RsgCijklWT(sv,gv,cijkl,WT)
% sv = [0;0;1];
% gv = [0;1;0];
% 
% Cij = diag([1 1 1 1 1 1]); % mu perturbation
% 
% cijkl = MS_cij2cijkl(Cij);
kappa3 = 1;%3^1.5; % this is the cube of v_p/v_s ratio
kappa6 = 1;%3^3;
ez = [0;0;1];
regCoeff = 10^(-30);
Rp = 0;
rx90 = 0;%rotx(90);

%% compute polarization vectors
svSH = V_cross(sv,ez);
svSH = (svSH + regCoeff*rx90*sv)/(norm(svSH)+regCoeff);
svSV = V_cross(sv,svSH);
gvSH = V_cross(gv,ez);
gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
gvSV = V_cross(gv,gvSH);


for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                if strcmp(WT,'PP')
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gv(k)*gv(l);
                elseif strcmp(WT,'SHSH')
                    Rp = Rp + sv(i)*svSH(j)*cijkl(i,j,k,l)*gv(k)*gvSH(l)*kappa6;                    
                elseif strcmp(WT,'SVSV')                    
                    Rp = Rp + sv(i)*svSV(j)*cijkl(i,j,k,l)*gvSV(k)*gv(l)*kappa6;
                elseif strcmp(WT,'PSH')
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gvSH(k)*gv(l)*kappa3;
                elseif strcmp(WT,'PSV')
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gvSV(k)*gv(l)*kappa3;
                elseif strcmp(WT,'SVSH')                                     
                    Rp = Rp + sv(i)*svSV(j)*cijkl(i,j,k,l)*gvSH(k)*gv(l)*kappa6;
                end
            end
        end
    end
end

end

function [c] = V_cross(a, b)
    % Can be quicker than regular cross

    c = [a(2).*b(3)-a(3).*b(2)
         a(3).*b(1)-a(1).*b(3)
         a(1).*b(2)-a(2).*b(1)];
    
return
end