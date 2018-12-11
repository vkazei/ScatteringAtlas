% compute WT (wavetype) radiation pattern on cijkl as function of sv gv for 3-D
% WT = 'PP' 'SHSH' 'SVSV'
function Rp = RsgDenWT(sv,gv,WT)

% example with P-P scattering:
%
% WT = 'PP';
% vertical incidence
% sv = [0;0;1];
% scattering along the second (y) axis
% gv = [0;1;0];
%

%% compute polarization vectors

% regularization of the denominator (SV, SH are not well defined near
% vertical propagation direction
ez = [0;0;1];
regCoeff = 10^(-30);

% could also be used for for setting SV and SH near vertical direction
Rp = 0;
rx90 = 0;%rotx(90);

svSH = V_cross(sv,ez);
svSH = (svSH + regCoeff*rx90*sv)/(norm(svSH)+regCoeff);
svSV = V_cross(sv,svSH);
gvSH = V_cross(gv,ez);
gvSH = (gvSH + regCoeff*rx90*gv)/(norm(gvSH)+regCoeff);
gvSV = V_cross(gv,gvSH);

%%

if strcmp(WT,'PP')
    Rp = Rp + sv'*gv;
elseif strcmp(WT,'SHSH')
    Rp = Rp + svSH'*gvSH;
elseif strcmp(WT,'SVSV')
    Rp = Rp + svSV'*gvSV;
elseif strcmp(WT,'PSH')
    Rp = Rp + sv'*gvSH;
elseif strcmp(WT,'PSV')
    Rp = Rp + sv'*gvSV;
elseif strcmp(WT,'SVSH')
    Rp = Rp + svSV'*gvSH;
end

end