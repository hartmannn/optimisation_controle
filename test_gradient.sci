clear;
cd /Users/armanddenis/Desktop/IMI_S4_local/TP_OC/Scilab
exec("Probleme_R.sce");
exec("Structures_R.sce");
exec("OraclePG.sci");
exec("Visualg.sci");
exec("Gradient_F.sci");

titrgr="Gradient test";
qcini=zeros(n-md,1);

//[fopt,xopt,gopt]=Gradient_F(OraclePG,qcini)

//[F,G,H,ind]=OraclePH(qcini,7)

exec("Wolfe_Skel.sci");
exec("Gradient_V.sci");

[fopt,xopt,gopt]=Gradient_V(OraclePG,qcini)
