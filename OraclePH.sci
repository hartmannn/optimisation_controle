function [F,G,H,ind]=OraclePH(qc,ind)
    F=0;
    s=n-md;
    G=zeros(s,1);
    H=zeros(s,s);
    
    if ind==2 | ind==4 | ind==7 then
        Pg1=q0+B*qc;
        Pd1=r.*Pg1.*abs(Pg1);
        
        Pd2=Ar*Pg1;
        
        F=(1/3)*Pg1'*Pd1+pr'*Pd2;
    end
    
    if ind==3 | ind==4 | ind==6 |ind==7 then
        termeDroit = B'*Ar'*pr;        
        termeGauche = B' * (r .* (q0 + B*qc) .* abs(q0 + B*qc));
        G= termeDroit + termeGauche
    end
    
    if ind==5 | ind==6 | ind==7 then
        H=B'*diag(2*r.*abs(q0+B*qc))*B;
    end
endfunction
