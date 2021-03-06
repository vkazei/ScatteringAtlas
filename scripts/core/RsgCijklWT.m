% compute WT (wavetype) radiation pattern on cijkl as function of sv gv for 3-D
% WT = 'PP' 'SHSH' 'SVSV'
% returns Rp as a number as a function of
% sv, gv - incident and scattered directions
% cijkl - perturbation of parameters
% WT - type of scattering experiment

function Rp = RsgCijklWT(sv,gv,cijkl,WT)

% example with P-P scattering:
%
% WT = 'PP';
% vertical incidence
% sv = [0;0;1];
% scattering along the second (y) axis
% gv = [0;1;0];%
% Cij = diag([1 1 1 1 1 1]); % mu perturbation
% mapping it to the tensor notation
% % cijkl = MS_cij2cijkl(Cij);

% normalization of scattering coefficients is off kappaI = 1
kappa3 = 1;%3^1.5; % this is the cube of v_p/v_s ratio
kappa6 = 1;%3^3;

% vertical vector for SV, SH polarization computation
ez = [0;0;1];

%% compute polarization vectors

% regularization of the denominator (SV, SH are not well defined near
% vertical propagation direction
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

%% evaluating the radiation patterns as scattering functions
% R = \sv\sp : cv : \gv\gp
% where
% \sv and \sp - are propagation and polarization vectors for incident
% wavefield, \gv and \gp - same for scattered


if strcmp(WT,'PP')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gv(k)*gv(l);
                end
            end
        end
    end
elseif strcmp(WT,'SHSH')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*svSH(j)*cijkl(i,j,k,l)*gv(k)*gvSH(l)*kappa6;
                end
            end
        end
    end
elseif strcmp(WT,'SVSV')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*svSV(j)*cijkl(i,j,k,l)*gvSV(k)*gv(l)*kappa6;
                end
            end
        end
    end
elseif strcmp(WT,'PSH')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gvSH(k)*gv(l)*kappa3;
                end
            end
        end
    end
elseif strcmp(WT,'PSV')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*sv(j)*cijkl(i,j,k,l)*gvSV(k)*gv(l)*kappa3;
                end
            end
        end
    end
elseif strcmp(WT,'SVSH')
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rp = Rp + sv(i)*svSV(j)*cijkl(i,j,k,l)*gvSH(k)*gv(l)*kappa6;
                end
            end
        end
    end
end

end

function [c] = V_cross(a, b)
% might be quicker than regular cross

c = [a(2).*b(3)-a(3).*b(2)
    a(3).*b(1)-a(1).*b(3)
    a(1).*b(2)-a(2).*b(1)];

return
end