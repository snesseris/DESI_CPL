#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(* Definitions *)


(* ::Text:: *)
(* Set  the  directory *)


(* data path *)
root\[Ellipsis]path="/home/snesseris/DESI_runs/CPLAN/";
SetDirectory[root\[Ellipsis]path];
datapath=root\[Ellipsis]path;

Off[General::"obspkg",NIntegrate::"inumr",NIntegrate::ncvb];

Print["Setting up the physics...."];


(* ::Text:: *)
(* Various  scales  and  data *)


(* ::Text:: *)
(* Various  params *)


Tcmb=2.7255;
c=299792.458;
cH0=2997.92458;
KtoeV=8.621738*10^-5;
cmtoeV=1/8065.544005;
\[Rho]cr\[Ellipsis]over\[Ellipsis]h2=8.098*10^-11 eV^4;


(* ::Text:: *)
(* DE  equation  of  state  w (a)  and  energy  density  f (a) *)


w[a_,w0_,wa_,wb_,wc_]:=w0+wa(1-a)+wb (1-a)^2+wc (1-a)^3
f[a_,w0_,wa_,wb_,wc_]:=a^(-3 (1+w0+wa+wb+wc)) E^(1/2 (-1+a) (6 wa-3 (-3+a) wb+(11+a (-7+2 a)) wc))
(*
 Tests: 
-1-1/3a D[Log[f[a,w0,wa]],a]-w[a,w0,wa]//Simplify
Exp[-3Integrate[(1+w[aa,w0,wa])/aa,{aa,1,a},GenerateConditions\[Rule]False]]//PowerExpand//FullSimplify
 *)


(* H(z) for the w0waCDM model *) 
Neff=3.044;
\[CapitalOmega]\[Gamma]h2=Pi^2/30  g\[Gamma] T^4/\[Rho]cr\[Ellipsis]over\[Ellipsis]h2/.T->Tcmb K/.g\[Gamma]->2/.GeV->10^9 eV/.K->KtoeV eV;
\[CapitalOmega]\[Gamma][h_]:=\[CapitalOmega]\[Gamma]h2/h^2;
\[CapitalOmega]\[Nu][h_,Neff_]:=Neff*7/8 (4/11)^(4/3) \[CapitalOmega]\[Gamma][h]
\[CapitalOmega]r[h_,Neff_]:=\[CapitalOmega]\[Gamma][h]+\[CapitalOmega]\[Nu][h,Neff];
H[a_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=100h Sqrt[om a^-3+\[CapitalOmega]r[h,Neff] a^-4+(1-om-\[CapitalOmega]r[h,Neff])f[a,w0,wa,wb,wc]]


(* ::Text:: *)
(* Sound  speed, see  Eq . (3.66)  in  the  old  version  of  Dodelson' s  book *)


cs[a_?NumberQ,obh2_?NumberQ]:=c Sqrt[1/(3(1+3/4 obh2/\[CapitalOmega]\[Gamma]h2 a))]
(* Old expression from Wang's paper on CMB shift params *)
cs\[Ellipsis]old[a_?NumberQ,obh2_?NumberQ]:=c/Sqrt[3(1+(31500 obh2 (Tcmb/2.7)^-4)a) ];
(* Redshift at decoupling (see Eqs. (66)-(68) in Komatsu et al.2009) *)
zcmb\[Ellipsis]EH[omh2_?NumberQ,obh2_?NumberQ]:=1048(1+0.00124(obh2)^-0.738)(1+((0.0783(obh2)^-0.238)/(1+39.5(obh2)^0.763))(omh2)^(0.560/(1+21.1(obh2)^1.81)))
(* Drag redshift (see Eq. (4) of Hu astro-ph/9709112) *)
zdrag\[Ellipsis]EH[omh2_?NumberQ,obh2_?NumberQ]:=1291 omh2^0.251/(1+0.659 omh2^0.828) (1+(0.313(omh2)^-0.419 (1+0.607omh2^0.674))(obh2)^(0.238omh2^0.223))
(* Redshift at decoupling Eq. (A4) from 2106.00428 *)
zcmb\[Ellipsis]GA[\[Omega]m_?NumberQ,\[Omega]b_?NumberQ]:=(391.672\!\(\*SuperscriptBox[\(\[Omega]m\), \(\[Minus]0.372296\)]\)+937.422\!\(\*SuperscriptBox[\(\[Omega]b\), \(\[Minus]0.97966\)]\))/(\!\(\*SuperscriptBox[\(\[Omega]m\), \(\[Minus]0.0192951\)]\) \!\(\*SuperscriptBox[\(\[Omega]b\), \(\[Minus]0.93681\)]\))+\!\(\*SuperscriptBox[\(\[Omega]m\), \(\[Minus]0.731631\)]\)
(* Redshift at decoupling Eq. (A2) from 2106.00428 *)
zdrag\[Ellipsis]GA[\[Omega]m_?NumberQ,\[Omega]b_?NumberQ]:=1/\[Omega]m^0.714129 (1+428.169 \[Omega]b^0.256459 \[Omega]m^0.616388+925.56 \[Omega]m^0.751615)


(* ::Text:: *)
(* Proper  (not  comoving)  angular  diameter  distance  with  curvature, see  Komatsu  et  al . 0803.0547  Eq . (2)  zcmb[om, obh2] *)


Clear[DLsol,DL,dL]
DLsol[om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=(DLsol[om,w0,wa,wb,wc,h]=NDSolve[{D[dL[zz]/(1+zz),zz]==1/(H[1/(1+zz),om,w0,wa,wb,wc,h]/(100h)),dL[0]==0},dL,{zz,0,1500},MaxSteps->Infinity])
DL[z_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=(c/(100h) dL[z]/.DLsol[om,w0,wa,wb,wc,h])[[1]]//Chop


(* ::Text:: *)
(* Angular  diameter  distance, see  Komatsu  et  al . 0803.0547  *)


DA[z_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=1/(1+z)^2 DL[z,om,w0,wa,wb,wc,h]
(* rs approximation from Eisenstein and Hu, ApJ,496,605, Eq. (26) *)
rs\[Ellipsis]EH[omh2_,obh2_]:=44.5Log[9.83/omh2]/(\[Sqrt](1+10(obh2)^(3/4)));
(* Comoving sound horizon at drag redshift *)
rs\[Ellipsis]num[ze_,om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=NIntegrate[cs[A,obh2]/(A^2H[A,om,w0,wa,wb,wc,h]),{A,0,1/(1+ze[om h^2,obh2])}]
(* Sound horizon from 2106.00428 *)
rsd\[Ellipsis]GA[\[Omega]m_,\[Omega]b_]:=Module[{a1=0.00257366,a2=0.05032,a3=0.013,a4=0.7720642,a5=0.24346362,a6=0.00641072,a7=0.5350899,a8=32.7525,a9=0.315473},1/(a1 \[Omega]b^a2+a3 \[Omega]b^a4 \[Omega]m^a5+a6 \[Omega]m^a7 )-a8/\[Omega]m^a9]
(* Sound horizon from 2212.04522 *)
rsd\[Ellipsis]HGM[\[Omega]m_,\[Omega]b_]:=147.05 (\[Omega]m/0.1432)^-0.23 (3.044/3.04)^-0.1 (\[Omega]b/0.02236)^-0.13
(* Dilation scale *)
Dv[z_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=(((1+z)DA[z,om,w0,wa,wb,wc,h] )^2  ((c*z)/H[1/(1+z),om,w0,wa,wb,wc,h]))^(1/3)
(* BAO Subscript[d, z] ratio *)
dz[z_?NumberQ,om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=rs\[Ellipsis]num[zdrag\[Ellipsis]GA,om,obh2,w0,wa,wb,wc,h]/Dv[z,om,w0,wa,wb,wc,h]
(* DH "distance" *)
DH[z_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=c/H[1/(1+z),om,w0,wa,wb,wc,h]
(* Scaled distance at recombination *)
R[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Sqrt[om (100h)^2]DA[zcmb\[Ellipsis]GA[om h^2,obh2],om,w0,wa,wb,wc,h](1+zcmb\[Ellipsis]GA[om h^2,obh2])/c
(* Angular scale of sound horizon at recombination *)
la[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=\[Pi] (DA[zcmb\[Ellipsis]GA[om h^2,obh2],om,w0,wa,wb,wc,h](1+zcmb\[Ellipsis]GA[om h^2,obh2]))/(rs\[Ellipsis]num[zcmb\[Ellipsis]GA,om,obh2,w0,wa,wb,wc,h])

(* DM "distance" *)
DM[z_?NumberQ,om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=(1+z)DA[z,om,w0,wa,wb,wc,h]

(* DM "distance" ratio over rs(zd) *)
DM\[Ellipsis]over\[Ellipsis]rs[z_?NumberQ,om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=((1+z)DA[z,om,w0,wa,wb,wc,h])/rs\[Ellipsis]num[zdrag\[Ellipsis]GA,om,obh2,w0,wa,wb,wc,h]


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(* CMB  and  BAO  \[Chi]^2  terms *)


(* ::Text:: *)
(* CMB : Planck18 *)


(* ::Text:: *)
(**)


Print["Reading the CMB data...."];
(* Central values and 1\[Sigma] errors for {R, Subscript[l, a], Subscript[\[CapitalOmega], b]h^2} from Planck 2018 for non-flat universe *)
(* See 1811.07425, page 7, Eq.(32) for details *)
datacmb={1.74451,301.76918,0.022483};
\[Sigma]CMB={0.093,0.0051,0.00016};
(* Inverse Covariance Matrix using  for {R, Subscript[l, a], Subscript[\[CapitalOmega], b]h^2} *)
covcmb=10^-8 \!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"2556.7782", "23212.222", 
RowBox[{"-", "57.345815"}]},
{"23212.222", "830122.02", 
RowBox[{"-", "628.56261"}]},
{
RowBox[{"-", "57.345815"}], 
RowBox[{"-", "628.56261"}], "2.5300094"}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
invcov=Inverse[covcmb];
(* \[Chi]^2 for R *)
vec[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:={R[om,obh2,w0,wa,wb,wc,h]-datacmb[[1]],la[om,obh2,w0,wa,wb,wc,h]-datacmb[[2]],obh2-datacmb[[3]]};
chi2\[Ellipsis]CMB[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=vec[om,obh2,w0,wa,wb,wc,h] . invcov . vec[om,obh2,w0,wa,wb,wc,h]


(* ::Text:: *)
(* BAO : DESI - Y1 *)


Print["Reading the DESI data...."]
(* The errors of the data, see Table 1 in 2404.03002 *)
BGS\[Ellipsis]\[Sigma]=0.15;
LRG\[Ellipsis]05\[Ellipsis]\[Sigma]={0.25,0.61,-0.445};
LRG\[Ellipsis]07\[Ellipsis]\[Sigma]={0.32,0.60,-0.420};
LRG\[Ellipsis]ELG\[Ellipsis]09\[Ellipsis]\[Sigma]={0.28,0.35,-0.389};
ELG\[Ellipsis]13\[Ellipsis]\[Sigma]={0.69,0.42,-0.444};
QSO\[Ellipsis]\[Sigma]=0.67;
Lya\[Ellipsis]QSO\[Ellipsis]233\[Ellipsis]\[Sigma]={0.94,0.17,-0.477};
(* The values *)
data\[Ellipsis]BAO={{0.30,7.93},{0.51,13.62},{0.51,20.98},{0.71,16.85},{0.71,20.08},{0.93,21.71},{0.93,17.88},{1.32,27.79},{1.32,13.82},{1.49,26.07},{2.33,39.71},{2.33,8.52}};


(* This makes the covariance/Fisher matrices *)
make\[Ellipsis]cij[vec_List]:={{vec[[1]]^2,vec[[1]] vec[[2]]vec[[3]]},{vec[[1]] vec[[2]]vec[[3]],vec[[2]]^2}}
make\[Ellipsis]cij[vec_]:=vec^2
Cij\[Ellipsis]DESI=ConstantArray[0,{12,12}];
Cij\[Ellipsis]DESI[[1,1]]=make\[Ellipsis]cij@BGS\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[{2,3},{2,3}]]=make\[Ellipsis]cij@LRG\[Ellipsis]05\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[{4,5},{4,5}]]=make\[Ellipsis]cij@LRG\[Ellipsis]07\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[{6,7},{6,7}]]=make\[Ellipsis]cij@LRG\[Ellipsis]ELG\[Ellipsis]09\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[{8,9},{8,9}]]=make\[Ellipsis]cij@ELG\[Ellipsis]13\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[10,10]]=make\[Ellipsis]cij@QSO\[Ellipsis]\[Sigma];
Cij\[Ellipsis]DESI[[{11,12},{11,12}]]=make\[Ellipsis]cij@Lya\[Ellipsis]QSO\[Ellipsis]233\[Ellipsis]\[Sigma];
Fij\[Ellipsis]DESI=Inverse[Cij\[Ellipsis]DESI];


(* BBN prior *)
chi2\[Ellipsis]BBN\[Ellipsis]DESI[obh2_?NumberQ]:=((obh2-0.02218)/0.00055)^2

(* The chi^2 *)
chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]Gil[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{V\[Ellipsis]data,V\[Ellipsis]theory,\[CapitalDelta]V},V\[Ellipsis]theory={Dv[data\[Ellipsis]BAO[[1,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[2,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[3,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[4,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[5,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[6,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[7,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[8,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[9,1]],om,w0,wa,wb,wc,h],Dv[data\[Ellipsis]BAO[[10,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[11,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[12,1]],om,w0,wa,wb,wc,h]}/rsd\[Ellipsis]HGM[om h^2,obh2];
V\[Ellipsis]data=data\[Ellipsis]BAO[[All,2]];
\[CapitalDelta]V=V\[Ellipsis]theory-V\[Ellipsis]data;
\[CapitalDelta]V . Fij\[Ellipsis]DESI . \[CapitalDelta]V]

chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]GA[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{V\[Ellipsis]data,V\[Ellipsis]theory,\[CapitalDelta]V},V\[Ellipsis]theory={Dv[data\[Ellipsis]BAO[[1,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[2,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[3,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[4,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[5,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[6,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[7,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[8,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[9,1]],om,w0,wa,wb,wc,h],Dv[data\[Ellipsis]BAO[[10,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[11,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[12,1]],om,w0,wa,wb,wc,h]}/rsd\[Ellipsis]GA[om h^2,obh2];
V\[Ellipsis]data=data\[Ellipsis]BAO[[All,2]];
\[CapitalDelta]V=V\[Ellipsis]theory-V\[Ellipsis]data;
\[CapitalDelta]V . Fij\[Ellipsis]DESI . \[CapitalDelta]V]

(* The chi^2 *)
chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]rdh[om_?NumberQ,rdh_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{V\[Ellipsis]data,V\[Ellipsis]theory,\[CapitalDelta]V},V\[Ellipsis]theory=h {Dv[data\[Ellipsis]BAO[[1,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[2,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[3,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[4,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[5,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[6,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[7,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[8,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[9,1]],om,w0,wa,wb,wc,h],Dv[data\[Ellipsis]BAO[[10,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[11,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[12,1]],om,w0,wa,wb,wc,h]}/rdh;
V\[Ellipsis]data=data\[Ellipsis]BAO[[All,2]];
\[CapitalDelta]V=V\[Ellipsis]theory-V\[Ellipsis]data;
DownValues[DLsol]=Take[DownValues[DLsol],-1];
\[CapitalDelta]V . Fij\[Ellipsis]DESI . \[CapitalDelta]V]

chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]num[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{V\[Ellipsis]data,V\[Ellipsis]theory,\[CapitalDelta]V},V\[Ellipsis]theory={Dv[data\[Ellipsis]BAO[[1,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[2,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[3,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[4,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[5,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[6,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[7,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[8,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[9,1]],om,w0,wa,wb,wc,h],Dv[data\[Ellipsis]BAO[[10,1]],om,w0,wa,wb,wc,h],DM[data\[Ellipsis]BAO[[11,1]],om,w0,wa,wb,wc,h],DH[data\[Ellipsis]BAO[[12,1]],om,w0,wa,wb,wc,h]}/rs\[Ellipsis]num[zdrag\[Ellipsis]GA,om,obh2,w0,wa,wb,wc,h];
V\[Ellipsis]data=data\[Ellipsis]BAO[[All,2]];
\[CapitalDelta]V=V\[Ellipsis]theory-V\[Ellipsis]data;
\[CapitalDelta]V . Fij\[Ellipsis]DESI . \[CapitalDelta]V]


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(* Pantheon + SHOES  SnIa  \[Chi]^2  term *)


(* ::Text:: *)
(* Read  the  SnIa  data *)


Print["Reading the Pantheon+ data...."]
path=".//data//data_Pantheon_Plus_SH0ES//";
data\[Ellipsis]Pantheon\[Ellipsis]plus=Import[path<>"Pantheon+SH0ES.dat","Table"];
ndat\[Ellipsis]Pantheon\[Ellipsis]plus=Length[data\[Ellipsis]Pantheon\[Ellipsis]plus]-1;
file\[Ellipsis]covs="Pantheon+SH0ES_STAT+SYS.zip";
script\[Ellipsis]mcut=10;
z\[Ellipsis]plus\[Ellipsis]min=0.01;
(*z\[Ellipsis]plus\[Ellipsis]SHOES\[Ellipsis]min=0.023;*)
z\[Ellipsis]plus\[Ellipsis]SHOES\[Ellipsis]min=0.01;


(* ::Text:: *)
(* The  redshifts  and  magnitudes  data *)


z\[Ellipsis]cmb=data\[Ellipsis]Pantheon\[Ellipsis]plus[[2;;-1,3]];
z\[Ellipsis]hel=data\[Ellipsis]Pantheon\[Ellipsis]plus[[2;;-1,7]];
Mb\[Ellipsis]corr=data\[Ellipsis]Pantheon\[Ellipsis]plus[[2;;-1,9]];
is\[Ellipsis]calibrator=data\[Ellipsis]Pantheon\[Ellipsis]plus[[2;;-1, 14]];
ceph\[Ellipsis]D=data\[Ellipsis]Pantheon\[Ellipsis]plus[[2;;-1,13]];
(* The "raw" covariance *)
Dij0\[Ellipsis]Stat\[Ellipsis]Sys = Import[path<>file\[Ellipsis]covs,{"ZIP","Pantheon+SH0ES_STAT+SYS.cov","List"}];
Cij\[Ellipsis]tot=Partition[Dij0\[Ellipsis]Stat\[Ellipsis]Sys[[2;;-1]],{Dij0\[Ellipsis]Stat\[Ellipsis]Sys[[1]]}];
Dij0\[Ellipsis]Stat\[Ellipsis]Sys=.


(* ::Input::Initialization:: *)
(* Pantheon+ chi^2 *)


(* ::Input::Initialization:: *)
(* We can easily find the position of the point\[Sigma] above z_min=0.01 *)
conditions\[Ellipsis]Pantheon\[Ellipsis]Plus=If[#>z\[Ellipsis]plus\[Ellipsis]min,1,0]&/@z\[Ellipsis]cmb;
pos\[Ellipsis]PP=Flatten[Position[conditions\[Ellipsis]Pantheon\[Ellipsis]Plus,val_/;val==1]];
(* The marginalized covariance *)
InvCovTotal\[Ellipsis]Pplus=Inverse[Cij\[Ellipsis]tot[[pos\[Ellipsis]PP,pos\[Ellipsis]PP]]];
Id\[Ellipsis]plus=ConstantArray[1,Length[pos\[Ellipsis]PP]];


(* ::Input::Initialization:: *)
chi2\[Ellipsis]Pantheon\[Ellipsis]plus\[Ellipsis]Marge[om_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ] := Module[{dist, \[CapitalDelta]m, A, B, Ee},
dist = Table[
If[z\[Ellipsis]cmb[[i]]>z\[Ellipsis]plus\[Ellipsis]min,5 Log10[(1+z\[Ellipsis]hel[[i]])/(1+z\[Ellipsis]cmb[[i]]) DL[z\[Ellipsis]cmb[[i]],om,w0,wa,wb,wc,h]]+25,Nothing],{i, 1, ndat\[Ellipsis]Pantheon\[Ellipsis]plus}];
\[CapitalDelta]m=Mb\[Ellipsis]corr[[pos\[Ellipsis]PP]]-dist;
A=\[CapitalDelta]m . InvCovTotal\[Ellipsis]Pplus . \[CapitalDelta]m;
B=\[CapitalDelta]m . InvCovTotal\[Ellipsis]Pplus . Id\[Ellipsis]plus;
Ee=Id\[Ellipsis]plus . InvCovTotal\[Ellipsis]Pplus . Id\[Ellipsis]plus;
A+Log[Ee/(2Pi)]-B^2/Ee]


(* ::Input::Initialization:: *)
(* end *)


(* ::Input:: *)
(*(* end *)*)


(* ::Text:: *)
(* Combine  all  data  (SnIa, CMB, BAO) *)


(* ::Text:: *)
(* The  number  of  data  points *)


ndat\[Ellipsis]BAO\[Ellipsis]DESI=12;
ndat\[Ellipsis]CMB=Length[datacmb];
ndat\[Ellipsis]tot\[Ellipsis]DESI=ndat\[Ellipsis]Pantheon\[Ellipsis]plus+ndat\[Ellipsis]CMB+ndat\[Ellipsis]BAO\[Ellipsis]DESI;

(* The total chi^2 for the new data *)
chi2\[Ellipsis]total[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{chi20},chi20=chi2\[Ellipsis]CMB[om,obh2,w0,wa,wb,wc,h]+chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]Gil[om,obh2,w0,wa,wb,wc,h]+chi2\[Ellipsis]Pantheon\[Ellipsis]plus\[Ellipsis]Marge[om,w0,wa,wb,wc,h];DownValues[DLsol]=DownValues[DLsol][[-2;;-1]];
If[Im[chi20]!=0,10^20,chi20]]

(* The individual parts of the chi^2 for the old data *)
chi2\[Ellipsis]total\[Ellipsis]vec[om_?NumberQ,obh2_?NumberQ,w0_?NumberQ,wa_?NumberQ,wb_?NumberQ,wc_?NumberQ,h_?NumberQ]:=Module[{chi20},chi20={chi2\[Ellipsis]CMB[om,obh2,w0,wa,wb,wc,h],chi2\[Ellipsis]BAO\[Ellipsis]DESI\[Ellipsis]Gil[om,obh2,w0,wa,wb,wc,h],chi2\[Ellipsis]Pantheon\[Ellipsis]plus\[Ellipsis]Marge[om,w0,wa,wb,wc,h]};DownValues[DLsol]=DownValues[DLsol][[-2;;-1]];
chi20]


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(* MCMC  stuff *)


(* Some code borrowed from *)
"https://mathematica.stackexchange.com/questions/6877/do-i-have-to-code-each-case-of-this-grid-full-of-plots-separately";
Options[plotGrid]={ImagePadding->40};
plotGrid[l_List,w_,h_,opts:OptionsPattern[]]:=Module[{nx,ny,sidePadding=OptionValue[plotGrid,ImagePadding],topPadding=0,widths,heights,dimensions,positions,frameOptions=FilterRules[{opts},FilterRules[Options[Graphics],Except[{ImagePadding,Frame,FrameTicks}]]]},{ny,nx}=Dimensions[l];
widths=(w-2 sidePadding)/nx Table[1,{nx}];
widths[[1]]=widths[[1]]+sidePadding;
widths[[-1]]=widths[[-1]]+sidePadding;
heights=(h-2 sidePadding)/ny Table[1,{ny}];
heights[[1]]=heights[[1]]+sidePadding;
heights[[-1]]=heights[[-1]]+sidePadding;
positions=Transpose@Partition[Tuples[Prepend[Accumulate[Most[#]],0]&/@{widths,heights}],ny];
Graphics[Table[Inset[Show[l[[ny-j+1,i]],ImagePadding->{{If[i==1,sidePadding,0],If[i==nx,sidePadding,0]},{If[j==1,sidePadding,0],If[j==ny,sidePadding,topPadding]}},AspectRatio->Full],positions[[j,i]],{Left,Bottom},{widths[[i]],heights[[j]]}],{i,1,nx},{j,1,ny}],PlotRange->{{0,w},{0,h}},ImageSize->{w,h},Evaluate@Apply[Sequence,frameOptions]]]
plempty[{left_,right_}]:=Plot[0,{x,-1,1},Frame->True,Axes->False,PlotStyle->White,FrameStyle->{{left,right},{White,White}},FrameTicks->None];
pleint=plempty[{Automatic,White}];
pleout=plempty[{White,White}];


(* ::Text:: *)
(* This  is  the  jump  function // Proposes  the  next  point *)


jumpfunc[point_,s_]:=point+RandomVariate[NormalDistribution[0,s],Length[point]]
jumpfuncmat[point_,cijtest_]:=RandomVariate[MultinormalDistribution[point,cijtest]]
(* \[CapitalDelta]\[Chi]^2 for 2 params, up to 3\[Sigma] *)
dx2=Table[2 InverseGammaRegularized[Nn/2,1-Erf[x\[Sigma]/Sqrt[2]]]/.Nn->2//N,{x\[Sigma],1,3}];
mysymmetrize[Cij_]:=(Cij+Transpose[Cij])/2;
(* Contours+color from the MCMC // Uses ConvexHull // Advanced++ *)
Needs["ComputationalGeometry`"];
Off[CompiledFunction::cflist];
Confid2D[data_List,color_,inter_]:=Module[{cv,cv1,cv2,qq},cv=ConvexHull[data,AllPoints->True];cv1=data[[#]]&/@cv;cv2=Join[cv1,{cv1[[1]],cv1[[2]]}];qq=Show[Graphics[{color,FilledCurve[BSplineCurve[cv2,SplineDegree->inter,SplineKnots->"Unclamped"]],Antialiasing->True}],Frame->True];qq]
normalize[bins_,counts_]:=counts/Max[counts]
(* Convergence diagnostics *)
Rstatistic[chain_,param_,burnin_]:=Module[{n,m,W,B,Rminus1},
m=chain//Length;
n=Table[chain[[i,burnin;;-1,1,param]]//Length,{i,1,m}]//Mean;
W=Mean[Table[Variance[chain[[j,burnin;;-1,1,param]]],{j,1,m}]];
B=n Variance[Table[Mean[chain[[j,burnin;;-1,1,param]]],{j,1,m}]];
Rminus1=Sqrt[(1-1/n)+1/n B/W]-1]
(* Utils *)
time[dt1_]:=Which[dt1<60,ToString[dt1]<>" secs",60<=dt1<1000,ToString[dt1/60.]<>" mins",1000<=dt1,ToString[dt1/3600.]<>" hrs"]


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(* Cache *)


obh2pl=0.02225;
ompl=0.3153;
s8pl=0.8111;
hpl=67.36/100;
obh2\[Ellipsis]BBN=0.02218;


(* ::Text:: *)
(* end *)


(* ::Text:: *)
(*(* Minimizations and MCMC *)*)


(* ::Text:: *)
(*(* Minimize all data *)*)


Print["Doing some tests...."]
chi2\[Ellipsis]total[0.315,0.0222,-1,0,0,0,0.673]//Print
chi2\[Ellipsis]total[0.315,0.0222,-1,0.1,0,0,0.673]//Print
chi2\[Ellipsis]total[0.315,0.0222,-1,0.1,0.1,0,0.673]//Print
chi2\[Ellipsis]total[0.315,0.0222,-1,0.1,0.1,0.1,0.673]//Print


(* ::Text:: *)
(*(* MCMC w0-wa-wb *)*)


(* Choose the params to vary *)
chi2[om_, w0_, wa_, wb_, h_] := chi2\[Ellipsis]total[om, obh2\[Ellipsis]BBN, w0, wa, wb, 0, h]


(* Set the initial priors for {Subscript[\[CapitalOmega], m],Subscript[w, 0]} *)
params = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(m\)]\)", "\!\(\*SubscriptBox[\(w\), \(0\)]\)", "\!\(\*SubscriptBox[\(w\), \(a\)]\)", "\!\(\*SubscriptBox[\(w\), \(b\)]\)", "h"};
priors = {{0.2, 0.4}, {-2, 0}, {-10, 10}, {-30,20}, {0.55, 0.85}};
nparams = priors // Length;
smallchainname = "smallchain_w0wawb_v1";
bigchainname1 = "bigchain_w0wawb_v1";


(* ::Text:: *)
(*(* Do a small MCMC to get a test cov. matrix *)*)


(* Do a small MCMC to get a test cov. matrix *)
Print["Starting a small MCMC...."]
step=0;
t1=AbsoluteTime[];
totalsteps=5000;
initpoint={0.308,-0.824,-0.795,-0.5,0.6805};
jumpsize={0.002,0.004,0.1,0.1,0.0015};
currentpoint=initpoint;
chi2current=chi2@@currentpoint;
smallchain={};
While[step<totalsteps+1,
( 
proposedpoint=jumpfuncmat[currentpoint,DiagonalMatrix[jumpsize^2]];
If[And@@Table[IntervalMemberQ[Interval[priors[[i]]],proposedpoint[[i]]],{i,1,Length[priors]}],
chi2prop=chi2@@proposedpoint;
jumpprob=Exp[-(chi2prop-chi2current)/2];
If[RandomReal[]<=Min[{1,jumpprob}],
currentpoint=proposedpoint;chi2current=chi2prop;AppendTo[smallchain,{currentpoint,chi2current}];,currentpoint=currentpoint;],currentpoint=currentpoint;];
);
step++]
Print["Time needed (secs): ",AbsoluteTime[]-t1]


(* Save the chain in zip format to save space *)
Export[smallchainname<>".zip",{smallchainname<>".m"->smallchain}];


(* Read in the chain // for postprocessing etc *)
smallchain=Import[smallchainname<>".zip",smallchainname<>".m"];


(* Best-fit from the MCMC *)
chi2min=Min[smallchain[[All,2]]]
bf=smallchain[[Position[smallchain,chi2min][[1,1]]]][[1]];
Print["The best-fit chi2 from small chain: ",{chi2min,bf}];


(* Cov. Matrix directly from the chain *)
chlen=Length[smallchain];
burnin=Max[Round[0.1*chlen],1];
Cijchain0=Covariance[smallchain[[All,1,All]][[burnin;;-1]]]//mysymmetrize
PositiveDefiniteMatrixQ[Cijchain0]//Print
SymmetricMatrixQ[Cijchain0]//Print
Cijchain0//Diagonal//Sqrt//Print


(* ::Text:: *)
(*(* end *)*)


(* ::Text:: *)
(*(* This is the main MCMC *)*)


(* Do the MCMC *)
step=0;
t1=AbsoluteTime[];
totalsteps=700000;
initpoint=bf;
currentpoint=initpoint;
chi2current=chi2@@currentpoint;
bigchain={};
While[step<totalsteps+1,
( 
proposedpoint=jumpfuncmat[currentpoint,Cijchain0];
If[And@@Table[IntervalMemberQ[Interval[priors[[i]]],proposedpoint[[i]]],{i,1,Length[priors]}],
chi2prop=chi2@@proposedpoint;
jumpprob=Exp[-(chi2prop-chi2current)/2];
If[RandomReal[]<=Min[{1,jumpprob}],
currentpoint=proposedpoint;chi2current=chi2prop;AppendTo[bigchain,{currentpoint,chi2current}];,currentpoint=currentpoint;],currentpoint=currentpoint;];
);
step++]
Print["Time needed (secs): ",AbsoluteTime[]-t1]


(* Save the chain in zip format to save space *)
Export[bigchainname1<>".zip",{bigchainname1<>".m"->bigchain}];


(* Best-fit from the MCMC *)
chi2min=Min[bigchain[[All,2]]]
bf=bigchain[[Position[bigchain,chi2min][[1,1]]]][[1]]
Print["The best-fit chi2 from big chain: ",{chi2min,bf}];


(* Mean values from the MCMC *)
bigchain[[All,1]]//Mean//Print
bigchain[[All,1]]//StandardDeviation//Print


Print["The end"];
