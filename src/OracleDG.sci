function [F,G,ind]=OracleDG(lambda,ind)
    a=Ar'*pr+Ad'*lambda;
    radicande=a./r;
    q#=-sign(radicande).*sqrt(abs(radicande));
    
    F=0;
    G=zeros(md,1);
    
    if ind==2 | ind==4 then
        F=-((1/3)*q#'*(r.*q#.*abs(q#))+pr'*(Ar*q#)+lambda'*(Ad*q#-fd));
    end
    
    if ind==3 | ind==4 then
        G=-Ad*q#+fd;
    end
endfunction
