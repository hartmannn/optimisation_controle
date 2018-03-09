function [F,G,ind]=OraclePG(qc,ind)
    F=0;
    s=n-md;
    G=zeros(s,1);
    
    Pg1=q0+B*qc;
    
    if ind==2 | ind==4 then
        Pd1=r.*Pg1.*abs(Pg1);
        
        Pd2=Ar*Pg1;
        
        F=(1/3)*Pg1'*Pd1+pr'*Pd2;
    end
    
    if ind==3 | ind==4 then
        termeDroit = B'*Ar'*pr;        
        termeGauche = B' * (r .* Pg1 .* abs(Pg1));
        G= termeDroit + termeGauche
    end
endfunction
