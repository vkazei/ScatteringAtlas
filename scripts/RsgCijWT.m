% this function evaluates the Cij radiation pattern for WT wave type scattering
% PP SHSH and SVSV are supported
% the input is two indexes (i,j) and two vectors (sv,gv)

function Rij = RsgCijWT(i,j,WT,sv,gv)

Cij = zeros(6);
Cij(i,j) = 1;
Cij(j,i) = 1;
%we need to add the symmetry

% then we convert into the tensor shape
cijkl = MS_cij2cijkl(Cij);

Rij = RsgCijklWT(sv,gv,cijkl,WT);

end


