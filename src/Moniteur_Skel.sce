///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
    clear;
// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   //stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_R.sce');
   exec('Structures_R.sce');
   
   // Affichage des resultats
   titrgr="";
   exec('Visualg.sci');
   
   // Verification  des resultats

   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   //exec('OraclePG.sci');
   //exec('Optim_Scilab.sci');
   //titrgr = "Fonction optim de Scilab sur le probleme primal";

    exec('OracleDG.sci');
    //exec('Wolfe_Skel.Sci');
    exec('Gradient_F.sci');

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal
   dim_p=n-md;
   dim_d=md;
   xini = 0.1 * rand(dim_d,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------

   // Exemple : la fonction "optim" de Scilab
   //
   //[fopt,xopt,gopt] = Optim_Scilab(OraclePG,xini);

   [fopt,xopt,gopt]=Gradient_F(OracleDG,xini);

// --------------------------
// Verification des resultats
// --------------------------

   [q,z,f,p] = HydrauliqueD(xopt);

   Verification(q,z,f,p);

//
