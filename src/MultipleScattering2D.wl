(* ::Package:: *)

BeginPackage["MultipleScattering2D`"]

ClearAll @@ Names["MultipleScattering2D`*"];

Needs["NumericalCalculus`"];

ExportListWave::usage="ExportListWave[{Xs,rs,n,{Max\[Omega],NN,\[Delta]}},listWave]";
ImportListWave::usage="ImportListWave[n,NN] returns {Header, listWave}, where Header is of the form {Xs,rs0,n0,{Max\[Omega],NN1,\[Delta]}}";

PlotWaves::usage="PlotWaves[waves], where waves is an output from AcousticScattering or AcousticImpulse, returns a list of plots of the wave at different instances of time.";
ReplaceMiddleZero::usage="Given a list rng, ReplaceMiddleZero[rng] returns rng but with middle element replaced with (Last[rng]-First[rng])/1000;";
CombineWaves::usage="CombineWaves[Xs,listWavesOnxtx\[Theta]xr,options] where Xs is a list of the position of the centre of the waves given in listWavesOnxtx\[Theta]xr, and listWavesOnxtx\[Theta]xr is a list of listWaveFromCoefficients. The call returns {dimension,data}, where dimensions = {{minx,maxx},{miny,maxy}} and data can be used as ListPlot@/data, where each element in data refers to one instant in time.";
SingleScattererCoefficients::usage="SingleScattererCoefficients[rng\[Omega],rngn, Scatterer Radius:0.2, Distance of incidence wave from scatterer centre:1.8]";
SeveralScattererCoefficients::usage="SingleScattererCoefficients[rng\[Omega],rngn, Scatterer Radius:0.2, Distance of incidence wave from scatterer centre:1.8]";
ScatteringCoefficients::usage="ScatteringCoefficients[n0_Integer,rs] returns a function TCoeffF[wE,X,k] which gives the scattering coefficients of a cyclinder that has been excited by the wave wE[{x,y},k] ";

FrequencyFromCoefficients::usage="FrequencyFromCoefficients[k_,Xs_,scatteringCoefficients_,listeners_:{{0,0}}] calculates the total wave, in frequency space, at given listener positions. ";
ListenersOutsideScatterers::usage= "ListenersOutsideScatterers[radius_,Xs_,rngX_,rngY_]";

ConvolutionTest::usage="ConvolutionTest[b,r0,rng\[Omega]Fourier] returns a plot of the convolution of b with the 2D Green's function.";
WaveF\[Theta]FromCoefficients::usage="WaveF\[Theta]FromCoefficients[listaOnx\[Omega],rng\[Omega],t0,rngr,Options] returns the result as function of \[Theta], which for each \[Theta] gives a list {{\!\(\*SubscriptBox[\(r\), \(0\)]\),\!\(\*SubscriptBox[\(w\), \(0\)]\)},...,{\!\(\*SubscriptBox[\(r\), \(n\)]\),\!\(\*SubscriptBox[\(w\), \(n\)]\)}}.";
listWaveFromCoefficients::usage=" listWaveFromCoefficients[listaOnx\[Omega],rng\[Omega],rngt,rng\[Theta],rngr,options] returns a pure list of the form listWOtx\[Theta]xr, where each element W ={\!\(\*SubscriptBox[\(t\), \(i\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \(i\)]\),\!\(\*SubscriptBox[\(r\), \(i\)]\)}}." ;

AcousticScattering::usage="AcousticScattering[CylinderRadius,{xImpulse,yImpulse},{t_min, t_max ,\[Delta]t}||t,{ Rmax,\[Delta]R,\[Delta]\[Theta]},options] returns a list { Wavetmin,...,Wavetmax}. Each Wavet = {Wavet\[Theta]0,..,Wavet\[Pi]}, each Wavet\[Theta]={{t,\[Theta], R0, wS[t,\[Theta],r0]}, ..., {t,\[Theta], Rmax, wS[t,\[Theta],Rmax]}}. ";
AcousticImpulse::usage="AcousticImpulse[t,options] returns a list {{r1,wS[r1]},...,{rn,wS[rn]}}.";





Begin["`Private`"]

ExportListWave[{Xs_,rs_,n_,{Max\[Omega]_,NN_,\[Delta]_}},listWave_]:=Module[{strDirectory},
   (*strDirectory = ToString[ Now[[1,3]]]<>"|"<>ToString[ Now[[1,2]]]<>"-"<>ToString[#[[1,1]]]<>":"<>ToString[#[[1,2]]] &@Now[[2]];*)
   strDirectory=ToString[Length[Xs]]<>" waves n="<>ToString[n];
   CreateDirectory[NotebookDirectory[]<>strDirectory];
   Export[NotebookDirectory[]<>strDirectory<>"/"<>ToString[NN]<>"-Info.txt",{Xs,rs,n,{Max\[Omega],NN,\[Delta]}}];
   Export[NotebookDirectory[]<>strDirectory<>"/"<>ToString[NN]<>"-listWave.mx",listWave,"MX"];
];

ImportListWave[NWaves_,n_,NN_]:=Module[{strDirectory,Header,listW},
   strDirectory=ToString[NWaves]<>" waves n="<>ToString[n];
   Header =ToExpression/@Import[NotebookDirectory[]<>strDirectory<>"/"<>ToString[NN]<>"-Info.txt","Data"];
   listW= Import[NotebookDirectory[]<>strDirectory<>"/"<>ToString[NN]<>"-listWave.mx","MX"];
 Return[{Header,listW}]
];


PlotWaves[listwOtx\[Theta]xr_,RS_:0,options__:{}]:=Module[{Amp,red,blue,mid,plotLegend,bound01,CoolColor},

  bound01 = If[#>1,1, If[#<0,0,#]]&;
  red={1,0.,0.0}; blue ={.0,0.2,1.}; mid={0.5,.9,0.5};
  
  If[Length@Dimensions[listwOtx\[Theta]xr]==4,
    Amp=Max[Transpose[Flatten[Abs@Re@listwOtx\[Theta]xr,2]][[4]]];  
    Print["Max Amp: ",Amp];
    CoolColor[z_] :=RGBColor@@bound01/@( blue+ (red-blue) (z +Amp )/(2Amp)+mid Exp[-2(z/Amp)^2] ) ;
    plotLegend= Placed[BarLegend[{CoolColor[#]&,{-Amp,Amp}},LegendLabel->"Wave Amplitude t="<>ToString[#]],Below] &;
    Return[
      ListDensityPlot[{#[[3]] Cos[#[[2]]],#[[3]] Sin[#[[2]]],#[[4]] }&/@Flatten[#,1],PlotRange->All,PlotLegends->plotLegend[#[[1,1,1]]],(*AspectRatio->1/2,*)ColorFunctionScaling->False,ColorFunction->CoolColor,PlotRange->All,RegionFunction->Function[{x,y,z},x^2+y^2>RS^2]]&/@(Re@listwOtx\[Theta]xr)
    ];
  ,
    If[Length@Dimensions[listwOtx\[Theta]xr]==2,
        Amp=Max[Transpose[Re@listwOtx\[Theta]xr][[2]]];  
        Print["Max Amp: ",Amp];
        CoolColor[z_] :=RGBColor@@bound01/@( blue+ (red-blue) (z +Amp )/(2Amp)+ 0.7 mid Exp[-2(z/Amp)^2] ) ;
        plotLegend=Placed[BarLegend[{CoolColor[#]&,{-Amp,Amp}},LegendLabel->"Wave Amplitude"],Right];
      Return[ListPlot[Re@listwOtx\[Theta]xr,PlotRange->All,PlotLegends->plotLegend,AspectRatio->1/2,ColorFunctionScaling->False,ColorFunction->CoolColor,PlotRange->All,Joined->True]];    
    ];
  ];

];

ReplaceMiddleZero[rng_,factor_:0.001]:=If[EvenQ[Length@rng]==True,
                             Print["The list passed to ReplaceMiddleZero has to have an odd number of elements"]; Return[rng] 
                           , 
                             Return@ReplacePart[rng,(Length[rng]-1)/2+1->(Last[rng]-First[rng])factor]
                          ];

SingleScattererCoefficients[rng\[Omega]Fourier_,rngn_:Range[-2,2], RS_:0.2,RI_:1.8]:=(
   (*If[First[rng\[Omega]]!=- Last[rng\[Omega]], Print["The frequency range rng\[Omega] has to satisfy First[rng\[Omega]] == - Last[rng\[Omega]]"]; Return[{{0}}]  ];*)
   (*change the frequency \[Omega]=0 to \[Omega] = \[Delta]\[Omega]/20. *)   
   Return[ -(I/4)Outer[  (BesselJ[-1+#1,#2 RS]-BesselJ[1+#1,#2 RS])/(HankelH1[-1+#1,#2 RS]-HankelH1[1+#1,#2 RS]) HankelH1[#1,#2 RI] & ,rngn, rng\[Omega]Fourier] ]; 
);

CombineWaves[Xs_, listWavesOnxtx\[Theta]xr_, options__:{}]:=Module[{Amp,dx, minx,maxx,miny,maxy,data,rngx,rngy,rngLt,
Nwaves = Length@Xs, Widths=(1-10^-3)Cos[\[Pi]/4]#[[1,1,-1,3]]&/@listWavesOnxtx\[Theta]xr,WaveOtXwtFxXy,plotLegend},

   WaveOtXwtFxXy=Table[Map[Interpolation[#,InterpolationOrder->1]&,
         Map[{Xs[[n]]+ {#[[3]] Cos[#[[2]]],#[[3]] Sin[#[[2]]]},#[[4]] }&/@Flatten[#,1]&, listWavesOnxtx\[Theta]xr[[n]]]]
      ,{n,Nwaves}];
   WaveOtXwtFxXy=Transpose[WaveOtXwtFxXy,{2,1}];

(*Maximum dimension where all waves have data*)
   minx =Max[Xs\[Transpose][[1]]-Widths];maxx =Min[Xs\[Transpose][[1]]+Widths];
   miny =Max[Xs\[Transpose][[2]]-Widths];maxy =Min[Xs\[Transpose][[2]]+Widths];
   If[NumericQ["MeshSize"/.options], dx = Min[(maxx-minx)/25,(maxy-miny)/25, "MeshSize"/.options ]; Print["Mesh Size:", dx],  dx = Min[(maxx-minx)/25,(maxy-miny)/25 ] ];
   rngx=Range[minx,maxx,dx];
   rngy = Range[miny,maxy,dx];
   rngLt =Range[ Length@listWavesOnxtx\[Theta]xr[[1]]];

   data =  Flatten[ Outer[{#1,#2,Sum[ WaveOtXwtFxXy[[#3,n]][#1,#2] ,{n,Nwaves}]}&,rngx,rngy,rngLt] ,{{3},{1,2}}];

   If[("Plot"/.options)=!=False, 
      Amp= Max[Transpose[Flatten[Abs@Re@listWavesOnxtx\[Theta]xr[[1]],2]][[4]] ]; 
      If[  NumericQ[pAmp =("PlotAmplitude"/.options)], Amp= pAmp Amp; Print["Amp: ",Amp]];
      red={1,0.,0.0}; blue ={.0,0.2,1.}; mid={0.5,.9,0.5};
      bound01 = If[#>1,1, If[#<0,0,#]]&;
      CoolColor[z_] :=RGBColor@@bound01/@( blue+ (red-blue) (z +Amp )/(2Amp)+mid Exp[-2(z/Amp)^2] );
      plotLegend= Placed[BarLegend[{CoolColor[#]&,{-Amp,Amp}},LegendLabel->"Wave Amplitude t="<>ToString[Round[10 #]/10 //N ]],Below] &;
      Return[ListDensityPlot[data[[#]] ,PlotRange->All,PlotLegends->plotLegend[listWavesOnxtx\[Theta]xr[[1,#]][[1,1,1]]],AspectRatio->(maxy-miny)/(maxx-minx),ColorFunctionScaling->False,ColorFunction->CoolColor,PlotRange->All]&/@rngLt];
   (* Return[{{{minx,maxx},{miny,maxy}},data}];*)
   ,
    Return[{{{minx,maxx},{miny,maxy}},data}];
   ];
];

ListenersOutsideScatterers[radius_,Xs_,rngX_,rngY_]:= 
Reap[Outer[
  If[Min[Norm/@Transpose[Xs\[Transpose]-{#1,#2}]]>radius ,Sow[{#1,#2}]]&
,rngX,rngY]][[2,1]];

FrequencyFromCoefficients[k_,Xs_,scatteringCoefficients_,listeners_:{{0,0}}]:=
Module[{coeffs =scatteringCoefficients[k], N0 = (Length[scatteringCoefficients[0.1][[1]]]-1)/2},

  Return@Map[
     Sum[Sum[coeffs[[j]][[n+N0+1]]HankelH1[n,k Norm[#-Xs[[j]]]] E^(I n ArcTan@@(#-Xs[[j]])),{n,-N0,N0}],{j,Length[Xs]}]&
  ,listeners]
];


ScatteringCoefficients[n0_Integer,rs_,Options__:{}]:=Module[{p,Dp, CC,TCoeffF,Xi, rng\[Theta], rngk,k,r,maxk, options =Flatten[ List@Options],wEs,wEN,w0k,ErrFkx\[Theta],rngErr},
   (*If[Options!= {}, options =Flatten[ List@Options]; Print@options];*)
   CC[m_,j_,n_]:= Sum[(-1)^((-n)/2 +j-m1) Binomial[m,m1]Binomial[2 j-m,-n/2+j-m1],{m1, Max[0,-n/2-j+m],Min[m,-n/2+j] }];
   p[n_,f_,X_,k_,rsOrder_]:=Sum[ (D[f[{x,y},k],{x,m},{y,2j-m}] /.{x-> X[[1]],y-> X[[2]]}) (CC[m,j,n](I^(m-2 j )) )/(m!(2 j -m)!) (rs/2)^(2 j),{j,Abs[n]/2,(rsOrder +1)/2},{m,0,2 j}]; 
   If[ ("BoundaryConditions"/.options)=== "Neumann",
      Dp[n_,f_,X_,k_,rsOrder_]:=Sum[ (D[f[{x,y},k],{x,m},{y,2j-m}] /.{x-> X[[1]],y-> X[[2]]}) (CC[m,j,n](I^(m-2 j )) )/(m!(2 j -m)!) j (rs/2)^(2 j-1),{j,Abs[n]/2,(rsOrder +1)/2},{m,0,2 j}]; 
      TCoeffF[f_,X_,k_]:= - (Dp[#,f,X,k,n0]/D[HankelH1[#,k r],r]/.r->rs)&/@Range[-n0,n0],
      TCoeffF[f_,X_,k_]:= - (p[#,f,X,k,n0]/HankelH1[#,k rs])&/@Range[-n0,n0]
   ];
   If[("PrintChecks"/.options)==True || NumericQ["CheckDistance"/.options],
      If[NumericQ["Maxk"/.options],maxk="Maxk"/.options, maxk=1.2/rs];
      If[NumericQ["CheckDistance"/.options],Xi ={0, "CheckDistance"/.options} , Xi= {5 rs, 5 rs} ];      
      wEs= HankelH1[0, #2 Norm[{#[[1]],#[[2]]}-Xi]/.Abs->Identity ]+HankelH1[1, #2 Norm[{#[[1]],#[[2]]}-Xi]/.Abs->Identity ]& ;

      rng\[Theta]={0,2 \[Pi],2\[Pi]/4.9};
      rngk=Range[rs, maxk,(maxk-rs)/21 ];
      wEN[\[Theta]_,k_]=  Sum[ p[n,wEs,{0,0},k,n0]E^(I n \[Theta]),{n,-n0,n0}];
      w0k= Abs[Plus@@(wEs[rs{Cos[#],Sin[#] },maxk/2.]&/@rng\[Theta] )];
      Print["|Wi[k]\!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\): ",w0k];
      ErrFkx\[Theta] = Abs[wEs[rs{Cos[#2],Sin[#2] },#1]-wEN[#2, #1]]&;
      (*Print["p[n,wEs,{0,0},k,n0,rs]: ", p[1,wEs,{0,0},0.5,n0]];*)

      rngErr=(Plus@@#/w0k)&/@Outer[ErrFkx\[Theta],rngk,rng\[Theta]];
      Print@ListPlot[{rngk,rngErr}\[Transpose],PlotRange->All,Joined->True,InterpolationOrder->2,Filling->Bottom,
                  AxesLabel->{"k","|W[k]-  \!\(\*SubscriptBox[\(W\), \(n\)]\)[k] \*SuperscriptBox[\(\[ExponentialE]\), \(\[ImaginaryI] n \[Theta]\)] \!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\)/ |W[k]\!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\) "},
                  PlotLabel->"W= \!\(\*SuperscriptBox[\(H\), \((0)\)]\)(k rs)+\!\(\*SuperscriptBox[\(H\), \((1)\)]\)(k rs) with "<>ToString[n0]<>" Fourier components, rs: "<>ToString[rs],
                  PlotLegends-> "Distance from source = "<>ToString[Norm@Xi]<>"- rs"
            ];
   ]; 
   Return[TCoeffF]
];

SeveralScattererCoefficients[N0_:3, Xs_:{{0.5,0.5},{-0.5,0.5}}, RS_:0.1,Options_:{}]:=
Module[{x,y,a,c,k,Maxk,Xi,parameters, parameterNames, options =Flatten[ List@Options], 
OutWave, rngParticles = Range[Length@Xs], OutWaves,SourceWave, kstep, kTon, CoeffEqsFkOn, vars, TCoeffsF, m1, argNSolveFk, coeffsFk,subvarsFk},
   Protect[x,y,a,c,k,rngParticles];
   $Assumptions=$Assumptions~Join~{x\[Element]Reals,y\[Element]Reals,k\[Element] Reals};
   
   (*We check if options sets any of the parameters, if not we assign some default values.*)
   parameters = {{Maxk,10.},{Xi,{0,0}}};
   parameterNames = {"Maxk","SourcePosition"};
   Evaluate[Transpose[parameters][[1]]] =If[Head[N@#[[1]]]===Head[N@#[[2]]],#[[1]],#[[2]]]&/@Transpose[{ReplaceAll[parameterNames,options],Transpose[parameters][[2]]}];

   (*Set the type of incident source wave WaveI and outgoing waves listWaveE*)
   OutWave[Xj_,aj_,n_]:= Module[{f},f[X_,k_]:= Plus@@Array[(aj[#] (HankelH1[#,k Norm[X-Xj] /.Abs->Identity]E^(I # ArcTan@@(X-Xj)) ))&,2n+1,-Abs@n];Return@f];
   OutWaves = OutWave[Xs[[#]],a[#],N0]&/@rngParticles; 
   If[NumericQ[("SourceWave"/.options)[{1,1},1]], 
     SourceWave[{x_,y_},k_] = ("SourceWave"/.options)[{x,y}-Xi,k],
     SourceWave[{x_,y_},k_] = I/4 HankelH1[0,k Norm[{x,y}-Xi] /.Abs->Identity]
   ];
   Protect[SourceWave,OutWaves];

   (*Here we decide how many Hankel functions to solve for in the scattered fields. 
  Note that Neumann always uses at least Subscript[H, 0] and Subscript[H, 1]. Bothe Dirchlett and Neumann 
   are setup to use the lowest amount of Hankels for 0< k\[LessEqual] RS. *)
   If[ ("BoundaryCondition"/.options) === "Neumann" ,m1=1,m1=0];
   If[("PrintChecks"/.options)== True, Print["Using a "<>If[m1===1,"Neumann","Dirchlett"]<>" boundary condition."] ];
   With[{m0=m1},
     kstep = ConstantArray[0,N0+2 -m0];
     kstep[[N0+2 - m0]]= Maxk;
     Map[(If[#==1,kstep[[#]]=0, kstep[[#]]=(RS^(2/(N0-m0))/Maxk^(2/(N0-m0)+1)kstep[[#+1]]^(#+1))^(1/#)];
         With[{k0=kstep[[#]],k1=kstep[[#+1]]},
           kTon[k_?NumericQ]:= #-1+m0/;k0<Abs[k]<=k1])&
      ,Range[N0+1 -m0,1,-1]];
     If[("PrintChecks"/.options)== True, 
      Print["To include up to \!\(\*SubsuperscriptBox[\(H\), \(n\), \((1)\)]\)(k r) in the scatterred field depends on k:"]; 
      MapThread[Print["up too "<>ToString[#1]<>" if "<>ToString[#2]]&,{Array[ "H["<>ToString[#-1+m0]<>"]"&, Length[kstep]-1],   Thread[Less[Drop[kstep,-1],"k",Drop[kstep,1]   ]]}];
     ];
   ];m1=.;
   
   vars=Array[a[#],2 N0+1,-Abs@N0]&/@rngParticles;
   TCoeffsF = ScatteringCoefficients[N0,RS,Sequence@@options,"CheckDistance"-> Union[Flatten@Outer[N@Norm[#1-#2]&, Xs, Xs]][[2]]];
   CoeffEqsFkOn[k_] = Map[ TCoeffsF[SourceWave,Xs[[#]],k]-vars[[#]]+ Plus@@First@Outer[TCoeffsF[#2,Xs[[#1]],k]&,{#},Drop[OutWaves,{#}] ] &,rngParticles];

   With[{n=#},
      argNSolveFk[k_?((kTon[Abs@#]==n )&)]=  
         {Thread[Flatten[#[[N0+1-n;;N0+1+n]]&/@CoeffEqsFkOn[k],1]==0]/.Thread[Rule[Flatten@Drop[vars,None,{N0+1 -n,N0+1 +n}],0]]
         ,Flatten@Take[vars,All,{N0+1-n,N0+1+n }]};
   ]&/@Range[0,N0];

   subvarsFk[k_]:=Flatten[{NSolve@@argNSolveFk[k], Thread[Rule[Flatten@Drop[vars,None,{N0+1-kTon[Abs@k],N0+1+kTon[Abs@k] }],0]]}];
   coeffsFk[k_]:=Module[{subNvars},
                      If[k==0, Print["Frequency k =0 has been replaced with k=", k=10^(-4$MachinePrecision)] ];
                      subNvars =subvarsFk[k];
                      If[!NumericQ[Total[vars/.subNvars,3]], 
                         Print["For k=",k,", either k> ",Maxk,", or possibly Ill-posed system "];
                         Return[ vars 0 ] 
                      ];
                      Return[vars/.subNvars]
                ];

  If[("PrintChecks"/.options)==True,
     Module[{rngk,rng\[Theta], SourceWaveNorms, BoundaryNorms},
        rngk =Range[kstep[[1]] + kstep[[2]]/5,Last[kstep],(Last[kstep]-kstep[[1]] - kstep[[2]]/5)/(2Length[kstep] +1)];
        rng\[Theta] =Range[0.,2.\[Pi], \[Pi]/(N0+2)];
        SourceWaveNorms= Abs@Apply[Plus,Outer[SourceWave[Xs[[#2]]+RS{Cos[#1],Sin[#1]},#3] &, rng\[Theta], rngParticles,rngk]];
        BoundaryNorms= Abs@Apply[Plus,Outer[SourceWave[Xs[[#2]]+RS{Cos[#1],Sin[#1]},#3] + Apply[Plus,Through[OutWaves[Xs[[#2]]+RS{Cos[#1],Sin[#1]},#3]]/.subvarsFk[#3]  ]& , rng\[Theta], rngParticles,rngk ]];
        PlotStyleBlended = Table[{Blend[{Blue,Red},x],Thick},{x,0,1,1/(Length[Xs]-1)}];

    Print@ListPlot[ {rngk,BoundaryNorms[[#]]/SourceWaveNorms[[#]]}\[Transpose]&/@rngParticles ,
                   Joined-> True,InterpolationOrder-> 2, PlotRange-> All,
                   AxesLabel-> {"Freq. k","Error"}, PlotLabel-> "Check ||Boundary[k]|| / ||IncidentWave[k]|| for all scatterers", 
                   PlotLegends->Evaluate[("X="<>ToString[#])&/@Xs], Filling->Axis];
    ];
  ];

Return[coeffsFk];
];

ConvolutionTest[b_,r0_,rng\[Omega]Fourier_]:= Module[{g,g\[Omega],listbO\[Omega],listwOt,NN=Length[rng\[Omega]Fourier], rngt,wFt,lpt,pt,\[Delta],T},
(*Here we check our discrete Fourier transform for the convolution of the Impulse with a Hankel function *)
         g\[Omega] = I/4 HankelH1[0,# r0] &;
         g =1/(2\[Pi]) HeavisideTheta[#-r0](#^2-r0^2)^(-.5) &; 
         T = (2\[Pi])/(rng\[Omega]Fourier[[2]]-rng\[Omega]Fourier[[1]]);
         \[Delta]=First@rng\[Omega]Fourier; 
         rngt = Range[0,T (NN-1)/NN ,T/NN];
         
         listbO\[Omega] =T /Sqrt[NN] Fourier[(b/@rngt)E^(I rngt \[Delta])];
         listwOt= Sqrt[NN]/T InverseFourier[listbO\[Omega] g\[Omega][rng\[Omega]Fourier]] E^(-I rngt \[Delta]); 
         wFt[t_]:= 1/(2\[Pi]) NIntegrate[ b[t-\[Tau]](\[Tau]^2-r0^2)^(-.5),{\[Tau],r0,Max[r0,t]}];       
         lpt=ListPlot[{rngt,Re@listwOt}\[Transpose],PlotStyle->{Orange,PointSize[Large]},PlotLegends->{"Frequency sampled \!\(\*SubsuperscriptBox[\(w\), \(N\), \(n\)]\)"}];
         pt =Plot[Re@wFt[t],{t,0,T},PlotRange-> All, PlotLabel->"Convolution of Impluse b(t) with 2D Green's function.",PlotLegends->{"analytic w[t]"},AxesLabel->{"t"}];
(*Print["wFt[t]: ", wFt[t], ",  wFt[1.]:",wFt[1.]]; *)        
    Return@Show[pt,lpt];
];

listWaveFromCoefficients[listaOnx\[Omega]_,rng\[Omega]Fourier_,rng\[Theta]_:Range[0.,\[Pi],0.2],rngR_:{},Options__:{}]:=
Module[{options= If[ Flatten[{Options}]=={},{},Flatten@ List@Options],b,\[Delta],bT,T,bN,rngt,rngtSample,Lrngt,rngr,rngn,maxn, Maxt,nt,
listbO\[Omega],listbN\[Omega]O\[Omega],bolPrintChecks,parameters,parameterNames, NN,t,wIN,\[Delta]\[Omega],listIntegrateSimp,listwOnxr,listwOtx\[Theta]xr},
        bolPrintChecks = TrueQ["PrintChecks"/.options];
 
    (*Select parameters from options or give default value*)
        parameters = ConstantArray[0,4]; 
        parameters[[1]] = E^(- 1/(1.00000001-(2#/bT -1)^2)) HeavisideTheta[#-.00000001]HeavisideTheta[bT-#]  &;
        parameters[[2]] =  2.; parameters[[3]] = 15; parameters[[4]] =  8.;
        parameterNames = {"Impulse","ImpulsePeriod","MaxTimeSamples","MaxTime"};
        parameters=If[NumericQ[#[[1]]]|| NumericQ[#[[1]]@0.],#[[1]],#[[2]]]&/@Transpose[{ReplaceAll[ parameterNames,options],parameters}];     
         b= parameters[[1]];bT=parameters[[2]]; nt=Round[parameters[[3]]]; Maxt = parameters[[4]];
    (*End parameters setup *)  

       \[Delta]\[Omega]=rng\[Omega]Fourier[[2]]-rng\[Omega]Fourier[[1]]; T=(2\[Pi])/\[Delta]\[Omega];
       \[Delta]=First@rng\[Omega]Fourier; If[\[Delta]==0, Print["Error: The frequencyy \[Omega]=0 does not work for Hankel functions.."]; Return[{}]];
       (*rng\[Omega]Fourier=rng\[Omega]~Join~Reverse@Drop[-rng\[Omega]+2\[Delta],1];*)
       NN=Length[rng\[Omega]Fourier];

       If[rngR=={}, rngr=Range[0.1, T+0.2, T 4/NN ] ,rngr=rngR];
       If[T< bT, Print["The frequency range is too coarse to cover the Impluse. The max time period covered is at most T=", T, ", while the impulse period is bT=",bT ] ];
       rngt = Range[0,T (NN-1)/NN ,T/NN ];
       (*listbO\[Omega]= T/Sqrt[2 NN+1] (#[[2]]~Join~#[[1]])&@TakeDrop[Fourier[b/@rngt//Chop], NN+1]; *)

       If[bolPrintChecks,
           Print@ConvolutionTest[b,T 0.5,rng\[Omega]Fourier];
           listbO\[Omega] =T /Sqrt[NN] Fourier[(b/@rngt) E^(I rngt \[Delta]) ];
           listbN\[Omega]O\[Omega] = NIntegrate[b[t]E^(I # t),{t,0.,T}]&/@rng\[Omega]Fourier;           
           Print@ListPlot[{{rng\[Omega]Fourier,Abs@listbN\[Omega]O\[Omega]}\[Transpose],{rng\[Omega]Fourier,Abs@listbO\[Omega]}\[Transpose]},PlotRange->All,PlotLabel->"Frequency of Impluse",PlotStyle->{{PointSize[Large]},{Thickness[0.002],Orange}},PlotLegends->{"b\[Omega] Fourier Transform","b\[Omega] Discrete Fourier"},AxesLabel->{"\[Omega]"}];
               
      ];
      listbO\[Omega] =T /Sqrt[NN] Fourier[(b/@rngt)E^(I rngt \[Delta])]; 
      maxn=(Length[listaOnx\[Omega]]-1 )/2; rngn=Range[ -maxn,maxn];
      listIntegrateSimp=(\[Delta]\[Omega] listbO\[Omega] # )&/@listaOnx\[Omega]; 

(*Print["L@rngt: ", Length@rngt, ", NN: ",NN];
Print["L@first@listaOnx\[Omega]: ", Length@First@listaOnx\[Omega]];*)
(*Return[{}];*)
(*Print["maxn: ", maxn,",  Length@rng\[Omega]:", Length@rng\[Omega]," ,  Length@listbO\[Omega]:", Length@listbO\[Omega],",  Length@listaOnx\[Omega][[1]]:", Length@(listaOnx\[Omega][[1]])];*)
      

     Lrngt=Min[Length[rngt], Round[Maxt Length[rngt]/ Max[rngt]]];  
      If[IntegerQ[nt],
         If[Length[rngt] >nt,rngtSample= rngt[[Round[#]]]&/@ Range[1,Lrngt,Lrngt/nt], rngtSample= rngt[[1;;Lrngt]] ];
      ,
        rngtSample= rngt[[1;;Lrngt]];
         
      ];
(*      Print["rngtSample:", rngtSample//Evaluate];
Print["nt:", nt]; Print["IntegerQ[nt]:", IntegerQ[nt]];
Print["listIntegrateSimp:", listIntegrateSimp];
Print["Outer[(E^(-I rng\[Omega]Fourier #2)) . HankelH1[#1,#3 rng\[Omega]Fourier ] &, rngn,rngtSample,rngr];:", Outer[E^(-I rng\[Omega]Fourier #2) . HankelH1[#1,#3 rng\[Omega]Fourier ] &, rngn,rngtSample,rngr]];
*)
      listwOnxr=1/(2.\[Pi])Outer[(E^(-I rng\[Omega]Fourier #2)listIntegrateSimp[[#1+maxn+1]]) . HankelH1[#1,#3 rng\[Omega]Fourier ] &, rngn,rngtSample,rngr];
(*Print["listwOnxr:", listwOnxr];*)
      listwOtx\[Theta]xr = Transpose[(Exp[I rngn #] .listwOnxr)&/@rng\[Theta],{2,1,3}];
      listwOnxr=.;
    Return[Outer[{rngtSample[[#1]],rng\[Theta][[#2]],rngr[[#3]],listwOtx\[Theta]xr[[#1,#2,#3]]}& ,Range@Length@rngtSample,Range@Length@rng\[Theta],Range@Length@rngr ]]
];

WaveF\[Theta]FromCoefficients[listaOnx\[Omega]_,rng\[Omega]_,t0_Real:1.,rngR_:{},Options__:{}]:=
Module[{options= If[ Flatten[{Options}]=={},{},Flatten@ List@Options],b,bT,T,bN,rngt,rngr,rngn,maxn,  listbO\[Omega],listbN\[Omega]O\[Omega],bolPrintChecks,parameters,parameterNames, NN,t,wIN,\[Delta]\[Omega],listIntegrateSimp,listIntegrate,listwOnxr,listwxrF\[Theta]},
       If[rngR=={}, rngr=Range[0.1,Max[t0]+0.2,Min[{(t0+0.2)/10,0.2 }]  ],rngr=rngR];
   (*    Print["options: \n", options];  *)
       bolPrintChecks = TrueQ["PrintChecks"/.options];

(*Select parameters from options or give default value*)
       parameters = { E^(- 1/(1.00000001-(#-1)^2)) HeavisideTheta[#]HeavisideTheta[2-#]  &,2.};
       parameterNames = {"Impulse","ImpulsePeriod"};
       parameters=If[NumericQ[#[[1]]]|| NumericQ[#[[1]]@0.],#[[1]],#[[2]]]&/@Transpose[{ReplaceAll[ parameterNames,options],parameters}];     
       b= parameters[[1]]; bT=parameters[[2]];
       \[Delta]\[Omega]=rng\[Omega][[2]]-rng\[Omega][[1]];T = (2\[Pi])/\[Delta]\[Omega];
       NN=(Length[rng\[Omega]]-1)/2;

      If[T< bT, Print["The frequency range is too coarse to cover the Impluse. The max time period covered is at most T=", T, ", while the impulse period is bT=",bT ] ];
      rngt = Range[0,T (2NN)/(2NN+1) ,T 1/(2NN+1) ];
      listbO\[Omega]= T/Sqrt[2 NN+1] (#[[2]]~Join~#[[1]])&@TakeDrop[Fourier[b/@rngt//Chop], NN+1]; 
      If[bolPrintChecks,
         listbN\[Omega]O\[Omega]=NIntegrate[b[t]E^(I # t),{t,0.,T}]&/@rng\[Omega];           
         Print@ListPlot[{{rng\[Omega],Re@listbN\[Omega]O\[Omega]}\[Transpose],{rng\[Omega],Re@listbO\[Omega]}\[Transpose]},PlotRange->All,PlotLabel->"Frequency of Impluse",PlotStyle->{{PointSize[Large]},{Thickness[0.002],Orange}},PlotLegends->{"b\[Omega] Fourier Transform","b\[Omega] Discrete Fourier"},AxesLabel->{"\[Omega]"}];
         bN[t_] = 1./(2\[Pi]) listbO\[Omega].Exp[-I t rng\[Omega] ] \[Delta]\[Omega]; 
         Print@Plot[{b[t],bN[t]//Re},{t,0,T},PlotRange->All,PlotStyle->{{Thickness[0.005]},{Thickness[0.003],Orange}},PlotLabel->"Impluse b(t) and frequency sampled impulse bN(t)",PlotLegends->{"b[t]","bN[t]"},AxesLabel->{"t"}];
       ];
listIntegrateSimp =\[Delta]\[Omega] E^(-I rng\[Omega] t0)  listbO\[Omega];
listIntegrateSimp=(listIntegrateSimp # )&/@listaOnx\[Omega];    
maxn=(Length[listaOnx\[Omega]]-1 )/2;
rngn=Range[ -maxn,maxn];
listwOnxr=1/(2.\[Pi])Outer[listIntegrateSimp[[#1+maxn+1]] . HankelH1[#1,#2( ReplaceMiddleZero@rng\[Omega] )] &, rngn,rngr];
    With[{rngn2 =rngn,rngr2=rngr,listwOnxr2= listwOnxr},
        listwxrF\[Theta][\[FormalT]_]:=   {rngr2,Exp[I rngn2 \[FormalT]].listwOnxr2}\[Transpose];
    ];
  Return[listwxrF\[Theta]]
];

AcousticScattering[RS_Real,XI_List,TimeMesh_,SpatialMesh_List:{0,0.2,0.1},Options__:{}]:=Module[{options= If[ Flatten[{Options}]=={},{}, List@Options],b, bT,bN,bolPrintChecks,T,RI= Norm[XI],maxR,\[Delta]r,\[Delta]\[Theta],\[CapitalOmega],T0,maxT,\[Delta]t,N\[Omega],Na,a,n,k,parameters,parameterNames},
(*Setup parameters and options*)  
  If[SpatialMesh[[1]]==0 || !NumericQ[SpatialMesh[[1]]],maxR=Sqrt[5](RS+RI),maxR=SpatialMesh[[1]]];
  \[Delta]r=SpatialMesh[[2]];\[Delta]\[Theta]=SpatialMesh[[3]];
  If[NumericQ[TimeMesh], T0=TimeMesh;maxT=T0;\[Delta]t=1;, If[Head[TimeMesh]==List && Length[TimeMesh]== 3, T0=TimeMesh[[1]]; maxT=TimeMesh[[2]]; \[Delta]t=TimeMesh[[3]];] ];

  bolPrintChecks = TrueQ["PrintChecks"/.options];
  
(*Select parameters from options or give default value*)
  parameters = {2.,20,5};
  parameterNames = {"ImpulsePeriod","FrequencyModes","NAngularModes"};
  parameters=If[NumericQ[#[[1]]],#[[1]],#[[2]]]&/@ Transpose[{ReplaceAll[ parameterNames,options],parameters}];

  If[ ("Impulse"/.options)==="Impulse",  
        b[\[Tau]_]= E^(- 1/(1-(\[Tau]-1)^2)) HeavisideTheta[\[Tau]]HeavisideTheta[2-\[Tau]];
        bT=2;
      , b =("Impulse"/.options);
        bT = parameters[[1]];
  ];

(*Approximate the impulse b with it's truncated Fourier series with N modes*)  
  N\[Omega]=parameters[[2]](*Number of freq modes*); T=maxT+bT (*Total time considered*); \[CapitalOmega] =2\[Pi] N\[Omega]/(T+1)//N;
  \[Delta]\[Omega]=\[CapitalOmega]/N\[Omega];rng\[Omega] =Range[-\[CapitalOmega]-\[Delta]\[Omega]/4,\[CapitalOmega],\[Delta]\[Omega]];

  listb\[Omega]=NIntegrate[b[t]E^(I # t),{t,0.,bT}]&/@rng\[Omega];
  If[bolPrintChecks,
    bN[t_] = 1./(2\[Pi]) listb\[Omega].Exp[-I t rng\[Omega] ] \[Delta]\[Omega];
    Print@Plot[{b[t],bN[t]//Re},{t,0,T},PlotRange->All,PlotStyle->{{Thickness[0.003]},{Orange}},PlotLabel->"Impluse b(t) and frequency sampled impulse bN(t)",PlotLegends->{"b[t]","bN[t]"},AxesLabel->{"t"}];
  ];
  If[("Boundary"/.options)==="Dirichlet",
   a[n_,k_]= -(I/4) BesselJ[n,RS k]/HankelH1[n,RS k] HankelH1[n,k RI],
   a[n_,k_]= -(I/4) D[BesselJ[n,RS k],k ]/D[HankelH1[n,RS k],k] HankelH1[n,k RI];
  ];
  
(*Test convergence of the \[Sum] a[n,k,r1]*)
(*The number of coefficients chosen should also reflect how many coefficients the expansion of H^(1)_0(k |XS - XI|) with Graff Addition theory needs. Typically 10 is enough and the more coeficients the closer XS is to XI. *)
  Na=parameters[[3]]; (*number of angular modes: \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = \(-Na\)\), \(Na\)]\ \(an\ HankelH1[n, k\ r]\ 
\*SuperscriptBox[\(\[ExponentialE]\), \(\(\[ImaginaryI]\)\(\ \)\(n\)\(\ \)\(\[Theta]\)\(\ \)\)]\)\)*) 
  If[bolPrintChecks, 
    s[m_]:= Sum[a[n,\[CapitalOmega]],{n,-m,m}];
    Print@ListPlot [(Abs[s[#]]&/@Range[Na+5]),PlotRange->All,PlotLabel->"Test convergence of the \!\(\*SubsuperscriptBox[\(\[Sum]\), \(-N\), \(N\)]\) \!\(\*SubscriptBox[\(a\), \(n\)]\)[k]\!\(\*SubsuperscriptBox[\(H\), \(n\), \((1)\)]\)(k r)",AxesLabel->{"N","\!\(\*SubsuperscriptBox[\(\[Sum]\), \(-N\), \(N\)]\) \!\(\*SubscriptBox[\(a\), \(n\)]\)[maxk]"}];
  ]; 

  listb\[Omega]Simspsons= \[Delta]\[Omega]/3 Array[listb\[Omega][[#]]If[EvenQ[#],2,4]&,Length[rng\[Omega]]]; listb\[Omega]Simspsons[[1]]=\[Delta]\[Omega]/3 listb\[Omega][[1]];listb\[Omega]Simspsons[[-1]]=\[Delta]\[Omega]/3 listb\[Omega][[-1]];
  rngr =Range[RS,maxR,\[Delta]r]; (*don't include r=0 due to singularities*)
  rng\[Theta] =Range[0.,\[Pi],\[Delta]\[Theta]]; 
  rngt =Range[T0,maxT,\[Delta]t]; (*best not to include T0=0 due to singularities*)

  gS2[r_,k_,\[Theta]_]= Sum[a[n,k]HankelH1[n,k r] E^(I n \[Theta]),{n,-Na,Na}];


  wSN[t_Real,\[Theta]_Real?NumericQ,r_Real?NumericQ]:=1/(2.\[Pi]) Re[listb\[Omega]Simspsons.(E^(-I # t) gS2[r,#,\[Theta]]&/@rng\[Omega])];
  
  Return[
     Map[ {Sequence@@#,If[#[[3]]-2 \[Delta]r< #[[1]]-Norm[RS{Cos[#[[2]]],Sin[#[[2]]]}-XI], wSN@@#,0]} &, Outer[List,rngt,rng\[Theta],rngr],{3}]
  ];

];

AcousticImpulse[t0_Real,{maxR_:0,\[Delta]r_:0.2},Options__:{}]:=Module[{options= If[ Flatten[{Options}]=={},{}, List@Options],b, bT,bN,bolPrintChecks,T,rngr,parameters,parameterNames,N\[Omega],rng\[Omega],listb\[Omega],t,listb\[Omega]Simspsons},
  If[maxR==0, maxR=t0+2 \[Delta]r];
  bolPrintChecks = TrueQ["PrintChecks"/.options];

(*Select parameters from options or give default value*)
  parameters = {2.,20};
  parameterNames = {"ImpulsePeriod","FrequencyModes"};
  parameters=If[NumericQ[#[[1]]],#[[1]],#[[2]]]&/@ Transpose[{ReplaceAll[ parameterNames,options],parameters}];

  If[ ("Impulse"/.options)==="Impulse",
        b = E^(- 1/(1.00000001-(2#/bT -1)^2)) HeavisideTheta[#-.00000001]HeavisideTheta[bT-#]  &;
        bT=2;
      , b =("Impulse"/.options);
        bT =parameters[[1]];
  ];

  N\[Omega]=parameters[[2]](*Number of freq modes*); T=t0+bT (*Total time considered*); \[CapitalOmega] =2\[Pi] N\[Omega]/(T+2)//N;
  \[Delta]\[Omega]=\[CapitalOmega]/N\[Omega];rng\[Omega] =Range[-\[CapitalOmega]-\[Delta]\[Omega]/4,\[CapitalOmega],\[Delta]\[Omega]];
Print["{\[Delta]\[Omega],\[CapitalOmega]}: ", {\[Delta]\[Omega], \[CapitalOmega]}];

  listb\[Omega]=NIntegrate[b[t]E^(I # t),{t,0.,bT}]&/@rng\[Omega];
  If[bolPrintChecks,
    bN[t_]=1/(2\[Pi]) listb\[Omega].Exp[-I t rng\[Omega] ] \[Delta]\[Omega];
    Print@Plot[{b[t],bN[t]//Re},{t,0,T}, PlotRange->All,PlotStyle->{{Thickness[0.003]},{Orange,Dashed}},PlotLabel->"Impluse b(t) and frequency sampled impulse bN(t)",PlotLegends->{"b[t]","bN[t]"},AxesLabel->{"t"}];
  ];

  (*listb\[Omega]= 1/3 Array[listb\[Omega][[#]]If[EvenQ[#],2,4]&,Length[rng\[Omega]]]; listb\[Omega][[1]]=1/3 listb\[Omega][[1]];listb\[Omega][[-1]]=1/3 listb\[Omega][[-1]];*)
  wIN[r_Real]:= 1/(2.\[Pi]) Re[\[Delta]\[Omega] listb\[Omega]. (E^(-I # t0) I/4 HankelH1[0,r #]&/@rng\[Omega])];

  rngr =Range[\[Delta]r,maxR,\[Delta]r];
  
  Return[{#,wIN[#]}&/@rngr];

];


End[];
EndPackage[];

