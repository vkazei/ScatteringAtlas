% compute WT (wavetype) radiation pattern on cijkl as function of sv gv for 3-D
% WT = 'PP' 'SHSH' 'SVSV'
function Rp = RsgDenWT(sv,gv,WT)
% sv = [0;0;1];
% gv = [0;1;0];
% 
% Cij = diag([1 1 1 1 1 1]); % mu perturbation
% 
% cijkl = MS_cij2cijkl(Cij);
ez = [0;0;1];
regCoeff = 10^(-30);
Rp = 0;
rx90 = 0;%rotx(90);

                if strcmp(WT,'PP')
                    Rp = Rp + sv'*gv;
                elseif strcmp(WT,'SHSH')
                    svSH = V_cross(sv,ez);
                    svSH = (svSH + regCoeff*rx90*sv)/(norm(svSH)+regCoeff);
                    gvSH = V_cross(gv,ez);
                    gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
                    Rp = Rp + svSH'*gvSH;
                elseif strcmp(WT,'SVSV')
                    svSH = V_cross(sv,ez);
                    svSH = (svSH + regCoeff*rx90*sv)/(norm(svSH)+regCoeff);
                    svSV = V_cross(sv,svSH);
                    gvSH = V_cross(gv,ez);
                    gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
                    gvSV = V_cross(gv,gvSH);
                    Rp = Rp + svSV'*gvSV;
                elseif strcmp(WT,'PSH')
                    gvSH = V_cross(gv,ez);
                    gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
                    Rp = Rp + sv'*gvSH;
                elseif strcmp(WT,'PSV')
                    gvSH = V_cross(gv,ez);
                    gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
                    gvSV = V_cross(gv,gvSH);
                    Rp = Rp + sv'*gvSV;
                elseif strcmp(WT,'SVSH')
                    svSH = V_cross(sv,ez);
                    svSH = (svSH + regCoeff*rx90*sv)/(norm(svSH)+regCoeff);
                    svSV = V_cross(sv,svSH);
                    gvSH = V_cross(gv,ez);
                    gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);                 
                    Rp = Rp + svSV'*gvSH;
                end
                
%Rp

end