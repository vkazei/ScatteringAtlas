% this function creates stiffness matrix from the Cij for orthorhombic
% given as a vector C11,C22,C33,C12,C13,C23,C44,C55,C66

function Cij = CijFromOrthV(CijOrth)

    Cij=zeros(6);
    for i=1:3
        Cij(i,i)=CijOrth(i);
        Cij(i+3,i+3)=CijOrth(6+i);
    end
    Cij(1,2)=CijOrth(4);
    Cij(2,1)=CijOrth(4);
    Cij(1,3)=CijOrth(5);
    Cij(3,1)=CijOrth(5);
    Cij(2,3)=CijOrth(6);
    Cij(3,2)=CijOrth(6);
end