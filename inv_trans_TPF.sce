
///////////////////
// inversion Abel transformation
// See works: 
// C. J. Dasch,... 
// P. S. Kolhe and A. K. Agrawal,... 
// H. Chehouani and M. El Fagrich,...
// for Scilab
// 02-2019
// Jury Barinov
//
// (two-point formula, TPF)
function [yout]=inv_trans_TPF(yi)
    n1=length(yi); // poins on R
    deff("[a]=Aij(i,j)","a=sqrt((j)^2-(i-1)^2)-sqrt((j-1)^2-(i-1)^2)");
    deff("[b]=Bij(i,j)","b=log((j+sqrt((j)^2-(i-1)^2))./((j-1)+sqrt((j-1)^2-(i-1)^2)))");
    //
    Dij=0;
    for i=1:n1
        D=0;
        for j=i:n1 
            if j>i & j<>2 then
                Dij=(1/%pi);
                Dij=Dij.*(Aij(i,j)-Aij(i,j-1)-(j).*Bij(i,j)+(j-2).*Bij(i,j-1));

            elseif j>i & j==2 then
                Dij=(1/%pi);
                Dij=Dij.*(Aij(i,j)-(j).*Bij(i,j)-1);//  
            end
            if j==i & i<>1 then
                Dij=(1/%pi);
                Dij=Dij.*(Aij(i,j)-(j).*Bij(i,j));
            end
            if (i==j & i==1) | j<i then
                Dij=0;
            end               
            D=D+Dij.*yi(j);
        end
        b(i)=D;
    end;
    yout=b;
endfunction
