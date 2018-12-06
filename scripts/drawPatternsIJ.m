function drawPatternsIJ(WT, i, j)

Cij = zeros(6);
Cij(i,j) = 1;
Cij(j,i) = 1;

mkdir patterns

patternArray = {'diffZ', 'diffX', 'diffXZ', 'diffXYZ', 'trans', 'reflZ'};%

for patType = patternArray
    drawPattern(WT,Cij,0,patType{1});
end
% drawPattern(WT,Cij,0,'diffXYZ');
% drawPattern(WT,Cij,0,'diffX');
% drawPattern(WT,Cij,0,'trans');
% drawPattern(WT,Cij,0,'reflZ');
% drawPattern(WT,Cij,0,'reflXZ');

mkdir ../latex/Fig/patterns/
system(['mkdir ../latex/Fig/patterns/C_',num2str(i),num2str(j)])
system(['mv patterns/* ../latex/Fig/patterns/C_',num2str(i),num2str(j)])

end