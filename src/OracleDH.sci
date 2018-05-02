function [F,G,H,ind]=OracleDH(lambda,ind)
    a=Ar'*pr+Ad'*lambda;
    radicande=a./r;
    q#=-sign(radicande).*sqrt(abs(radicande));
    
    F=0;
    G=zeros(md,1);
    H=zeros(md,md);
    
    if ind==2 | ind==4 | ind==7 then
        F=-((1/3)*q#'*(r.*q#.*abs(q#))+pr'*(Ar*q#)+lambda'*(Ad*q#-fd));
    end
    
    if ind==3 | ind==4 | ind==6 |ind==7 then
        G=-Ad*q#+fd;
    end
    
    if ind==5 | ind==6 | ind==7 then
        H=(1/2)*Ad*diag(1./sqrt(abs(a.*r)))*Ad';
    end
endfunction
