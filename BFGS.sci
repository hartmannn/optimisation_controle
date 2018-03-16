exec("Wolfe_Skel.sci");
function [fopt,xopt,gopt]=BFGS(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//                         Algorithme de BFGS                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de l algorithme de BFGS";
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
   s=n-md;
   W=eye(s,s);
   Gprec=zeros(s,1);
   xprec=zeros(s,1);
   
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
      delta_x=x-xprec;
      delta_G=G-Gprec;
      facteur=1/(delta_G'*delta_x);
      
      W=(eye(s,s)-facteur*(delta_x*delta_G'))*W...
        *(eye(s,s)-facteur*(delta_G*delta_x'))...
        +facteur*(delta_x*delta_x');
      
      D = -W*G;
      
//    - on garde en mémoire les valeurs du point et du gradient
      Gprec=G;
      xprec=x;

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
   disp('Fin de la methode de gradient a pas fixe')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction


