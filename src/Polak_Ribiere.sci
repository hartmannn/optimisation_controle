function [fopt,xopt,gopt]=Polak_Ribiere(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//                Algorithme de Polak-Ribière                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de l algorithme de Polak-Ribière";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"1";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   
// - initialisation des variables :
   dim_x=size(xini)(1);
   Gprec=zeros(dim_x,1);
   D=zeros(dim_x,1);
   
   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [F,G] = Oracle(x,ind);
      
//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente
      beta_=((G-Gprec)'*G)/(Gprec'*Gprec);
      D = -G+beta_*D;  
      
//    - on garde en mémoire la valeur du gradient
      Gprec=G;  

//    - calcul de la longueur du pas de gradient

      [alpha,ok]=Wolfe(alphai,x,D,Oracle);

//    - mise a jour des variables

      x = x + (alpha*D);
      
//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de Polak_Ribiere')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction

