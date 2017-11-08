(* ::Package:: *)

(* ::Code::Bold:: *)
BeginPackage["MultipleScattering2D`"];

ClearAll @@ Names["MultipleScattering2D`*"];

Needs["NumericalCalculus`"];

ExportListFrequency::usage="";
ImportListFrequency::usage="";
ImportBoundaryFrequency::usage="ImportBoundaryFrequency[stres,lenXs,HighestHankelH1Order,options] returns Outer[{#1, #2, totalWave[#1,#2] }&, rngPositions,rngFrequency,1,1] "
ExportListWaves::usage="ExportListWaves[{ListScattererPositions, ScattererRadius, HighestHankelH1Order,{MaxFrequency,MaxFrequencySamples}},listWaves,opts__:{}]";
ImportListWaves::usage="ImportListWaves[NWaves,HighestHankelH1Order,MaxFrequencySamples]";
ExportListCylindricalWave::usage="ExportListCylindricalWave[{Xs,rs,n,{Max\[Omega],NN}},listWave]";
ImportListCylindricalWave::usage="ImportListCylindricalWave[n,NN] returns {Header, listWave}, where Header is of the form {Xs,rs0,n0,{Max\[Omega],NN1,\[Delta]}}";

ListenersOutsideScatterers::usage= "ListenersOutsideScatterers[radius_,Xs_,rngX_,rngY_]";

DrawScatterers::usage= "DrawScatterers[Xs, listrs:0.1,opts:OptionsPattern[]]";
SetParametersAndReturnOptions::usage="SetParametersAndReturnOptions[parameters,options] will take parameters of the form {{\!\(\*SubscriptBox[\(String\), \(1\)]\), \!\(\*SubscriptBox[\(Symbol\), \(1\)]\), \!\(\*SubscriptBox[\(DefaultValue\), \(1\)]\)},{\!\(\*SubscriptBox[\(String\), \(2\)]\), \!\(\*SubscriptBox[\(Symbol\), \(2\)]\), \!\(\*SubscriptBox[\(DefaultValue\), \(2\)]\)},...} and will set \!\(\*SubscriptBox[\(Symbol\), \(i\)]\) to the value from the rule for \!\(\*SubscriptBox[\(String\), \(i\)]\) in options, or will set to \!\(\*SubscriptBox[\(DefaultValue\), \(i\)]\) and return an updated options. Note N is applied to all DefaultValues.";
GenerateParticles::usage="GenerateParticles[NParticles,rs,{{x1,x2},{y1,y2}}] returns a list of positions {x,y} for NParticles with radius rs in the rectangle {{x1,x2},{y1,y2}}";
StatsFromParticles::usage="StatsFromParticles[rs_Real,Xs_List] return {e, s}, where e is the expected distance between the particles and s is the standard deviation of the distances.";

ListPlotSequence::usage="PlotWave[listxyWaveOt_,rngt_:{},Options__:{}] returns a list of ListDensityPlot of each element of listxyWaveOt. Options for ListDensityPlot may be included in Options.";
PlotCylindricalWave::usage="PlotCylindricalWave[CylindricalWave, rs_:0, options__:{} ] where CylindricalWave is given from CylindricalWaveFromCoefficients or each element of ListCylindricalWaveDueToImpulse, and Plotwaves returns a list of plots of the wave at different instances of time.";

ListWavesDueToImpulse::usauge"ListWavesDueToImpulse generates a wave in time and even saves it.";
CombineWaves::usage="CombineWaves[Xs,listCompleteWavesOnxtx\[Theta]xr,options] where Xs is a list of the position of the centre of the waves given in listWavesOnxtx\[Theta]xr, and listCompleteWavesOnxtx\[Theta]xr is a list of CylindricalWaveFromCoefficients. The call returns {dimension,data}, where dimensions = {{minx,maxx},{miny,maxy}} and data can be used as ListPlot@/data, where each element in data refers to one instant in time.";
rng\[Omega]FourierOffset::usage="rng\[Omega]FourierOffset[MaxFrequency, NumFrequencySamples] returns the frequency range used in Fourier";
ListCylindricalWaveDueToImpulse::usage="ListCylindricalWaveDueToImpulse[Xs, rs:0.1, n0:3, opts__:{}] returns a list where each element represent the scatterered wave of one scatterer in the form listWOtx\[Theta]xr, where each W ={\!\(\*SubscriptBox[\(WaveAmp\), \(i\)]\),\!\(\*SubscriptBox[\(t\), \(i\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \(i\)]\),\!\(\*SubscriptBox[\(r\), \(i\)]\)}}.";

FrequencyFromScatterers::usage="FrequencyFromScatterers[NN , X, Xs, rs_:0.1, rngk_:Range[0.1,10.1,1], Options__:{}] returns the frequency response of the system measure at X={x,y}, where the position of the scatterers are given by the list Xs, their radius rs, the order of the hankel functions used for the scatterers is NN and frequency range rngk";
AttenuationFromListFrequency::usage="AttenuationFromListFrequency[X_, listFreq_,options:OptionsPattern[]] returns a list with each element being the attenuation against frequency measured at each position in X.  "

TotalWave::usage="TotalWave[Xs_List?(Depth[#]\[Equal]2&),listCoeffsFk_Symbol,opts:OptionsPattern[]] returns a function wave, where wave[{x,y},k] gives the total wave in the frequency at {x,y}.";
OutWave::usage="OutWave[{x,y},a,n0] returns general scatterered wave from cylinder. "

CylindricalWaveFromCoefficients::usage=" CylindricalWaveFromCoefficients[listaOnx\[Omega],rng\[Omega],rngt,rng\[Theta],rngr,options] returns Outer[{#1,#2,#3, Wave[#1,#2,#3]}& ,rngt,rng\[Theta],rngr], in other words at depth 3 each element is of the form W ={\!\(\*SubscriptBox[\(WaveAmp\), \(i\)]\), \!\(\*SubscriptBox[\(t\), \(i\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \(i\)]\),\!\(\*SubscriptBox[\(r\), \(i\)]\)}}." ;
MultipleScatteringCoefficients::usage="MultipleScatteringCoefficients[ xs_:{{0.5,0.5},{-0.5,0.5}}, RS_:0.1, N0_:3,Options_:{}] returns a function ListScatteringCoefficients[k] that for each frequency k returns a list of the scattering coefficients of each cylinder with centre given by xs and takes into account all multiple scattering. See ?ScatteringCoefficientOperator for more detail.";
ScatteringCoefficientOperator::usage="ScatteringCoefficientOperator[n0_Integer,rs] returns a function TCoeffF[wE,X,k] which gives a list {\!\(\*SubscriptBox[\(a\), \(-N\)]\),\!\(\*SubscriptBox[\(a\), \(\(-N\) + 1\)]\),...,\!\(\*SubscriptBox[\(a\), \(N - 1\)]\),\!\(\*SubscriptBox[\(a\), \(N\)]\)} where \!\(\*UnderoverscriptBox[\(\[Sum]\), \(n = \(-N\)\), \(N\)]\)\!\(\*SubscriptBox[\(a\), \(n\)]\)\!\(\*SubsuperscriptBox[\(H\), \(n\), \((1)\)]\)[k r]\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\[ImaginaryI]\\\ k\\\ n\)]\) is the scattered wave from a cyclinder with centre X that has been excited by the wave wE[{x,y},k].";

SingleScattererCoefficientsFromImpulse::usage="SingleScattererCoefficients[rng\[Omega],MaxOrderOfHankel, ScattererRadius:0.1, Distance of incidence wave from scatterer centre:1.8]";
ConvolutionTest::usage="ConvolutionTest[b,r0,rng\[Omega]Fourier] returns a plot of the convolution of b with the 2D Green's function.";


(* ::Code::Bold:: *)
Begin["`Private`"]
Unprotect[failTag];ClearAll[failTag];Protect[failTag];

SetDirectory[NotebookDirectory[]];

ImpulseDirectory[LXs_,n_]:= Directory[]<>"/"<>ToString[LXs]<>"-ImpulseWaves_n="<>ToString[n];
FrequencyDirectory[LXs_,n_]:= Directory[]<>"/"<>ToString[LXs]<>"-WaveDispersion_n="<>ToString[n];
RulesQ=Catch[Map[If[Head@# =!= Rule,Throw[False]]& ,Flatten@List@##]; Throw[True]  ]&;
protectedQ=Quiet[Or@@Map[MatchQ[#,Protected]&, Attributes[#]]]& ;
variableQ=Quiet@ListQ@NSolve[{},#]&;
(*UnprotectVariableQ = (variableQ [#] && !protectedQ[#] &);*)

labelData="volfrac:";
strVol[volfrac_]:=ToString[volfrac]<>"%";

ExportListFrequency[{X_List,Xs_List,rs_Real,n_Integer},listFreq_List, listBoundaryCheck_List:{},opts:(__?RulesQ):{}]:=
Module[{strDirectory, optDirectoryExtension, options =Flatten[ List@opts],plotSequence,label,listlabels, strvol = strVol[StatsFromParticles[rs,Xs][[1]]] },
  label = labelData<>strvol;

   If[(("DirectoryExtension"/.options)=!= "DirectoryExtension") && StringQ["DirectoryExtension"/.options],
     strDirectory=FrequencyDirectory[Length[Xs],Round@n]<>("DirectoryExtension"/.options);
     optDirectoryExtension="DirectoryExtension"->("DirectoryExtension"/.options) 
   ,
     strDirectory=FrequencyDirectory[Length[Xs],Round@n];
     optDirectoryExtension =Sequence[];
   ];

   If[Quiet[CreateDirectory[strDirectory]]=== $Failed,
      Print["Creating the directory: ",strDirectory,", has failed. No data was saved"]; Return[{}],   
      Print["Saving wave data in the directoy:"<>strDirectory<>"/"];
   ];
  (*listBoundaryCheck is of the form Outer[{#1,#2,totalWave[#1,#2]}&,rngX,rngk,1,1] *) 
   If[listBoundaryCheck!={}, Export[strDirectory<>"/"<>label<>"-BoundaryCheck.mx",listBoundaryCheck,"MX"]];
   If[ (listlabels=Quiet[ ToExpression/@Import[strDirectory<>"/fileLabels.txt","List"]])=== $Failed
     ,Export[strDirectory<>"/fileLabels.txt",{strvol}] 
     ,Export[strDirectory<>"/fileLabels.txt",Union[listlabels~Join~{strvol}]] 
   ];

   Export[strDirectory<>"/"<>label<>"-Info.txt",{rs,X,Xs}];
   Export[strDirectory<>"/"<>label<>"-listFrequency.mx",listFreq,"MX"];
   Export[strDirectory<>"/"<>label<>"-options.mx",options,"MX"];

   Print["ImportListFrequency returns {Header,optsFreq,listFrequency}"];
   Print["where Header = {rs,listX,listXs} and each element of listFrequency is the frequency reponse at each X \[Element] listX."];
   If[listBoundaryCheck!={}, Print["To load the boundary data call:"];
      Print[Inactive[ImportBoundaryFrequency][StatsFromParticles[rs,Xs][[1]],Length[Xs],n,optDirectoryExtension]];
   ];
  
  Return[Inactive[ImportListFrequency][StatsFromParticles[rs,Xs][[1]],Length[Xs],n,optDirectoryExtension]]
];

ImportBoundaryFrequency[Nvol_,lenXs_Integer,n_Integer,optDirectoryExtension:OptionsPattern[]]:= 
Module[{strDirectory,options = Flatten[List@optDirectoryExtension],label = labelData<>strVol[Round[10Nvol]/10.]},
  If[(("DirectoryExtension"/.options)=!= "DirectoryExtension") && StringQ["DirectoryExtension"/.options],
     strDirectory=FrequencyDirectory[lenXs,Round@n]<>("DirectoryExtension"/.options);
   ,
     strDirectory=FrequencyDirectory[lenXs,Round@n];
   ];
  Print["ImportBoundaryFrequency returns  Outer[{#1, #2, totalWave[#1,#2] }&, rngPositions,rngFrequency,1,1]."];
  Return[Import[strDirectory<>"/"<>label<>"-BoundaryCheck.mx"]]
]; 

ImportListFrequency[Nvol_,lenXs_Integer,n_,opts:(__?RulesQ):{}]:=
Module[{strDirectory,Header,listF, options =Flatten[ List@opts],label = labelData<>strVol[Round[10Nvol]/10.]},
   If[(("DirectoryExtension"/.options)=!= "DirectoryExtension") && StringQ["DirectoryExtension"/.options],
     strDirectory=FrequencyDirectory[lenXs,Round@n]<>("DirectoryExtension"/.options);
   ,
     strDirectory=FrequencyDirectory[lenXs,Round@n];
   ];
   listF= Import[strDirectory<>"/"<>label<>"-listFrequency.mx","MX"];
   If[listF== $Failed, Throw[{},"Failed to Import file:"<>strDirectory<>"/"<>label<>"-listFrequency.mx"] ];
   Header =ToExpression/@Import[strDirectory<>"/"<>label<>"-Info.txt","Data"];
   options=Import[strDirectory<>"/"<>label<>"-options.mx","MX"];
 Return[{Header,options,listF}]
];

ExportListWaves[{Xs_,rs_,n_,{Max\[Omega]_,NN_}},listWave_,opts:(__?RulesQ):{}]:=Module[{strDirectory, options =Flatten[List@opts],plotSequence},
   strDirectory=ImpulseDirectory[Length[Xs],n];
   If[Quiet[DirectoryQ[strDirectory]===False && CreateDirectory[strDirectory]=== $Failed],
      Print["Creating the directory: ",strDirectory,", has failed. No data was saved"]; 
      Return[]   
   ];
   Print["Saving data in the directoy:"<>strDirectory<>"/"];
   Export[strDirectory<>"/"<>ToString[Round@NN]<>"-Info.txt",{Xs,rs,n,{Max\[Omega],Round@NN}}];
   Export[strDirectory<>"/"<>ToString[Round@NN]<>"-listWave.mx",listWave,"MX"];
   Export[strDirectory<>"/"<>ToString[Round@NN]<>"-options.mx",options,"MX"];
   Print["To load this wave data call: \n {Header,options,listW}= ImportListWaves["<>ToString[Length[Xs]]<>
","<>ToString[n]<>","<>ToString[Round@NN]<>"]"];
   If["SaveGIF"/.options,
      plotSequence= ListPlotSequence[listWave,  If[ Head["RangeTime"/.options]=== List, "RangeTime"/.options, {}],options];
      Export[strDirectory<>"/"<>ToString[Round@NN]<>"-Wave_"<>ToString[Length@Xs]<>"-Scatterers.GIF",Show[#,DrawScatterers[Xs]]&/@plotSequence];
      Print["The GIF has beed saved as: `",strDirectory<>"/"<>ToString[Round@NN]<>"-Wave-"<>ToString[Length@Xs]<>"Scatterers.GIF`"]
   ];
  Return[strDirectory]
];

ImportListWaves[NWaves_,n_,NN_]:=Module[{strDirectory,Header,listW,options},
   strDirectory=ImpulseDirectory[NWaves,n];
   Header =ToExpression/@Import[strDirectory<>"/"<>ToString[Round@NN]<>"-Info.txt","Data"];
   listW= Import[strDirectory<>"/"<>ToString[Round@NN]<>"-listWave.mx","MX"];
   options=Import[strDirectory<>"/"<>ToString[Round@NN]<>"-options.mx","MX"];
 Return[{Header,options,listW}]
];

ExportListCylindricalWave[{Xs_,rs_,n_,{Max\[Omega]_,NN_}},listWave_]:=Module[{strDirectory, msg},
   strDirectory=ImpulseDirectory[Length[Xs],n];
   If[Quiet[CreateDirectory[strDirectory]]=== $Failed, 
      Print["Creating the directory: ",strDirectory,", has failed. No data was saved"]; Return[],   
      Print["Saving wave data in the directoy:"<>strDirectory<>"/"];
   ];
   Export[strDirectory<>"/"<>ToString[NN]<>"-Info.txt",{Xs,rs,n,{Max\[Omega],NN}}];
   msg=Export[strDirectory<>"/"<>ToString[NN]<>"-listWave.mx",listWave,"MX"];
   Print["To load this wave data call: \n ImportListCylindricalWave[",Length[Xs],",",n,",",NN "]"];
];

ImportListCylindricalWave[NWaves_,n_,NN_]:=Module[{strDirectory,Header,listW},
   strDirectory=ImpulseDirectory[NWaves,n];
   Header =ToExpression/@Import[strDirectory<>"/"<>ToString[NN]<>"-Info.txt","Data"];
   listW= Import[strDirectory<>"/"<>ToString[NN]<>"-listWave.mx","MX"];
 Return[{Header,listW}]
];

ListenersOutsideScatterers[radius_,Xs_,rngX_,rngY_]:= 
Reap[Outer[
  If[Min[Norm/@Transpose[Xs\[Transpose]-{#1,#2}]]>radius ,Sow[{#1,#2}]]&
,rngX,rngY]][[2,1]];


(* ::Code::Bold:: *)
Options[DrawScatterers]= {"ReceiverPositions"-> {},(*"ImpulsePosition"-> {},*)"Axes"-> True};
DrawScatterers[Xs_, listrs_:0.1, opts:OptionsPattern[]]:= 
Module[{particles,  options =Flatten[ List@opts], Xreceivers,receivers},
    Xreceivers = If[Length@OptionValue["ReceiverPositions"]==2 && Depth@OptionValue["ReceiverPositions"]==2,{OptionValue["ReceiverPositions"]},OptionValue["ReceiverPositions"] ];
    particles=Thread[{listrs,Xs}];
    receivers = Flatten@Map[{Orange,Disk[#, 2 Max@listrs,{ArcTan@@(-#)+Pi/4, ArcTan@@(-#) +2 Pi -Pi/4 }]}&,Xreceivers];
    (*{RGBColor[.2,.6,.5],Disk[OptionValue["ImpulsePosition"] , 2 Max@listrs,{3Pi/2-Pi/4,3Pi/2+Pi/4}]}*)	
    Return@Graphics[{Gray}~Join~(Disk[#[[2]],#[[1]]]&/@particles)~Join~receivers,Axes->OptionValue["Axes"],AxesLabel-> {"x","y"}];
];

SetAttributes[SetParametersAndReturnOptions,HoldFirst];
SetParametersAndReturnOptions[parameters_,options_,opts_:{}]:=Module[{tmpoptions},
   Evaluate[parameters[[All,2]]] =If[ Head[N@#[[1]]]===Head[N@#[[2]]],N@#[[1]],N@#[[2]]]&/@Transpose[{ReplaceAll[parameters[[All,1]],options],parameters[[All,3]]}];
   tmpoptions=DeleteDuplicates[ N[options~Join~Thread@Rule[parameters[[All,1]],N@parameters[[All,2]]] ]];
   If[Head@Unevaluated@opts==Symbol && Head[opts]== List,opts=tmpoptions];
  Return@tmpoptions
];
StatsFromParticles[rs_Real,Xs_List]:= Module[{distances,expectedDistance, NParticles= Length@Xs,volfrac,stdDeviation},
  distances=DeleteCases[Flatten@Outer[N@Norm[#1-#2]&, Xs, Xs,1,1],0.]-2rs;
  If[(distances=!={})&&(N[distances] =!= N[0]),
    expectedDistance =Total[distances/(2NParticles (NParticles-1))];
    stdDeviation =Sqrt[Total[(expectedDistance-distances)^2]/(2NParticles (NParticles-1))];
    volfrac= 100 NParticles Pi rs^2 /((Max[#[[2]]]-Min[#[[2]]])(Max[#[[1]]]-Min[#[[1]]])&@(Xs\[Transpose]));
  ,
    expectedDistance=0;
    stdDeviation=0;
    volfrac= 1.;
  ];
  volfrac= Round[10 volfrac]/10.;
  Return[{volfrac,expectedDistance,stdDeviation}]
];

Options[GenerateParticles]= {"ParticlePositions"-> {},"Norm"-> Norm};
GenerateParticles[NParticles_Integer,rs_Real,{{x1_,x2_},{y1_,y2_}},OptionsPattern[]]:=
  Module[{NNewP = NParticles, rngDelete, Xs=OptionValue["ParticlePositions"], norm=OptionValue["Norm"]},  
    Do[
      Xs=Xs~Join~Transpose[{RandomReal[{x1,x2},NNewP],RandomReal[{y1,y2},NNewP]}]; 
      rngDelete= Quiet@DeleteDuplicates@First@Transpose@DeleteDuplicates[ Position[Outer[norm[#1-#2]&,Xs,Xs,1,1] , _?(0.< #<2rs&)], (#1[[{2,1}]]== #2)& ];
      Xs=Delete[Xs,{#}&/@rngDelete];
      If[(NNewP=Length@rngDelete)== 0, Return[]];
    ,{100}];
    If[Length[Xs]<NParticles, Throw[{{},0,0},"There was not enough space to put "<>ToString[NParticles]<> " particles"]];
  
  Return[Xs]
 ];

PlotCylindricalWave[listwOtx\[Theta]xr_,RS_:0,opts:(__?RulesQ):{}]:=Module[{Amp,red,blue,mid,plotLegend,bound01,CoolColor, options =Flatten[ List@opts],
patOptions,optsListDensityPlot},

  bound01 = If[#>1,1, If[#<0,0,#]]&;
  red={1,0.,0.0}; blue ={.0,0.2,1.}; mid={0.5,.9,0.5};
  
  If[Length@Dimensions[listwOtx\[Theta]xr]==4,
    If[options!= {},
      patOptions= Alternatives@@Transpose[ Options[ListDensityPlot]/.Rule->List/.RuleDelayed->List][[1]];
      optsListDensityPlot= Extract[options,Position[Transpose[options/.Rule-> List/.RuleDelayed-> List][[1]],patOptions ]],
      optsListDensityPlot=Sequence[]
    ];
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

ListPlotSequence[listxyWaveOt_,rngt_:{},Opts__:{}]:=Module[{options =Flatten[ List@Opts], Amp,red,blue,mid,darkred,darkblue,
plotLegend,bound01,CoolColor, patOptions,optsListDensityPlot,boundz,f},

   If[NumericQ["PlotAmplitude"/.options], Amp="PlotAmplitude"/.options, Amp = Max[Transpose[Flatten[listxyWaveOt,1]][[3]]] ]
   If[NumericQ["PlotAmplitudeRatio"/.options], Amp= ("PlotAmplitudeRatio"/.options) Amp, Amp= 0.7 Amp ];
   If[options!= {},
     patOptions= Alternatives@@Transpose[ Options[ListDensityPlot]/.Rule->List/.RuleDelayed->List][[1]];
     optsListDensityPlot= Extract[options,Position[Transpose[options/.Rule-> List/.RuleDelayed-> List][[1]],patOptions ]],
     optsListDensityPlot=Sequence[]
   ];
   (*bound01 = If[#>1,1, If[#<0,0,#]]&;*)
   darkred={.2,0.0,0.0}; red={1,0.,0.0};  mid={1.,1.,1.}; blue ={.0,0.2,1.};darkblue={0.,0.0,0.2};
   boundz = If[#>1,1, If[#<-1,-1,#]]&;
   f= Interpolation[{{-1,darkblue},{- 0.6 ,blue},{0,mid},{0.6,red},{1,darkred}},InterpolationOrder-> 1];
   CoolColor[z_] :=RGBColor@@f[boundz[z/Amp]];
   (*CoolColor[z_] :=RGBColor@@bound01/@( blue+ (red-blue) (z +Amp )/(2Amp)+mid Exp[-2(z/Amp)^2] );*)
   If[Length[rngt]=== Length[listxyWaveOt],
      plotLegend= Placed[BarLegend[{CoolColor[#]&,{-Amp,Amp}}, LegendLabel->"Wave Amplitude at t="<>ToString[Round[10 rngt[[#]]]/10 //N ]],Below] &,
      plotLegend= Placed[BarLegend[{CoolColor[#]&,{-Amp,Amp}}, LegendLabel->"Wave Amplitude"],Below] &
   ];
  Return[ListDensityPlot[listxyWaveOt[[#]],
      PlotRange->All,PlotLegends->plotLegend[#], optsListDensityPlot ,ColorFunctionScaling->False,
      ColorFunction->CoolColor,PlotRange->All]&/@Range[Length[listxyWaveOt]]
  ]  
];

FrequencyFromScatterers[X_List, Xs_List, rs_:0.1, n_, rngk_List:(Range[0,4,1]+0.1), opts:(__?RulesQ):{}]:=
Module[{options =Flatten[ List@opts],listCoeffFk, parameters,  Wi, freq, bcX, bck,totalWave,listFreq}, 
  parameters = {{"SourceWave", Wi, Exp[I #2 #1[[1]]]&}};
  options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts];
  listCoeffFk=MultipleScatteringCoefficients[Xs,rs,n,Unevaluated@options];
  If[Head@Unevaluated@opts == Symbol && Head@opts == List,opts=options];

  totalWave= TotalWave[Xs,listCoeffFk,options];
  listFreq = Outer[{#2,totalWave[#1,#2]}&, If[Depth[X] ==2, {X}, X], rngk,1,1];
  If["SaveData"/.options,
    bcX= ({0,rs}+#)&/@Union@Flatten[Extract[Xs,Position[Transpose[Xs][[2]],#]]&/@ {Min@#,Max@#}&@(Transpose[Xs][[2]]),1];
    bck={First@rngk, Median@rngk,Last@rngk}//N;
    Print@ExportListFrequency[{X,Xs,rs,n},listFreq, Outer[{#1,#2, totalWave[#1,#2]}&,bcX,bck,1,1], options];
  ];

Return@listFreq
];

AttenuationFromListFrequency[X_, listFreq_,options:OptionsPattern[]]:=
Module[{SourceWave,listSourceWaves,data},
  If[Length@X=!= Length@listFreq, Print["Length of X (the receiver position list) is different
   from listFreq (the list of frequency measured at the receivers)"] ; 
    Return[{}];
  ];
  If[NumericQ[("SourceWave"/.Flatten@List@options)[X[[1]],0.2]],
    SourceWave= ("SourceWave"/.Flatten@List@options), 
    SourceWave=Exp[I #2 #1[[2]]]&;
  ];
  listSourceWaves=Thread[SourceWave[#[[1]],#[[2]]\[Transpose][[1]]],List,-1]&/@Thread[{X,listFreq}];
Return[{#[[2]]\[Transpose][[1]],Re[(#[[2]]\[Transpose][[2]])/(#[[1]])-1 ]}\[Transpose]&/@Transpose[{listSourceWaves,listFreq}]]
]; 


(* ::Code::Bold:: *)
rng\[Omega]FourierOffset[Maxk_,NN_,c\[Delta]_:0.3]:=Module[{\[Delta]\[Omega], \[Delta]},
  \[Delta]\[Omega] = Maxk /(NN+c\[Delta]); (*\[DoubleLongLeftArrow] NN \[Delta]\[Omega] \[Equal] Max\[Omega] - c\[Delta] \[Delta]\[Omega] *) 
  \[Delta] = c\[Delta] \[Delta]\[Omega]; (*As you can see c\[Delta] determines how much of \[Delta]\[Omega] is used to offset the frequency*)
  (*T= 2\[Pi] NN /Maxk; \[Delta]t=T/(2NN+1) ; *)
  Print["For the frequency given, we have mesh t \[Element] Range[0, ",2\[Pi] NN /Maxk,",",\[Pi] /Maxk 2NN/(2NN+1),"]" ];
  Return[Range[0,Maxk-\[Delta],\[Delta]\[Omega]]~Join~Range[-Maxk+\[Delta], -\[Delta]\[Omega], \[Delta]\[Omega]] +\[Delta]]
];


(* ::Code::Bold:: *)
ListWavesDueToImpulse::fail="ListWavesDueToImpulse failed because `1` ";

ListWavesDueToImpulse[Xs_List, rs_:0.1, n0_:3, opts:(__?RulesQ):{}]:=Module[{options =Flatten[ List@opts], 
x0,x1,y0,y1,maxR, parameters, listWaves,Xi, Maxk, c\[Delta]=0.3, listCylindricalWaves, PAmp, meshSize, plotSequence},
   parameters = {{"MaxFrequency", Maxk, 1.1/rs},{"ImpulsePosition",Xi,{0.,0.}},{"PlotAmplitudeRatio",PAmp,0.8},{"MeshSize",meshSize,rs/2}};
   options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts];

  If[!TrueQ@NumericQ[maxR="MaxRadius"/.options],  maxR=2rs+ 1.6 Max[Outer[Norm[#1-#2]&,#,#,1,1]]&[{Xi}~Join~Xs] ];
  Catch[
    listCylindricalWaves= ListCylindricalWaveDueToImpulse[Xs, rs(*Scatterer radius*), n0(*Order of HankelH1 functions *), Unevaluated@options];
    {{{x0,x1},{y0,y1}},listWaves}= CombineWaves[{"ImpulsePosition"/.options}~Join~Xs,listCylindricalWaves, Unevaluated@options];
  
    If[("SaveGIF"/.options)||("SaveData"/.options),
      If[!NumericQ[AspectRatio/.options], AppendTo[options, AspectRatio->(y1-y0)/(x1-x0)]];
      ExportListWaves[{Xs,rs,n0,{Maxk,"MaxFrequencySamples"/.options}},listWaves,Unevaluated@options, "RangeTime"-> First@Transpose[listCylindricalWaves[[1]][[All,1,1]]]];
    ];
   Return[{{{x0,x1},{y0,y1}},listWaves}]
   ,_failTag, (Message[ListWavesDueToImpulse::fail,Style[First@#2,Red]];Return[#1])&
  ];
];

CombineWaves[Xs_List,listCylindricalWaveOtx\[Theta]xr_List, opts:(__?RulesQ):{}]:=Module[{Amp,dx, minx,maxx,miny,maxy,
listxyWaveOt,rngx,rngy,rngLt,Nwaves = Length@Xs, Widths=(1-10^-3)Cos[\[Pi]/4]#[[1,1,-1,3]]&/@listCylindricalWaveOtx\[Theta]xr
,WaveOtXwtFxXy,plotLegend, options =Flatten[ List@opts]},

   WaveOtXwtFxXy=Table[Map[Interpolation[#,InterpolationOrder->1]&,
         Map[{Xs[[n]]+ {#[[3]] Cos[#[[2]]],#[[3]] Sin[#[[2]]]},#[[4]] }&/@Flatten[#,1]&, listCylindricalWaveOtx\[Theta]xr[[n]]] ]
      ,{n,Nwaves}];
   WaveOtXwtFxXy=Transpose[WaveOtXwtFxXy,{2,1}];

(*Maximum dimension where all waves have data*)
   minx =Max[Xs\[Transpose][[1]]-Widths];maxx =Min[Xs\[Transpose][[1]]+Widths];
   miny =Max[Xs\[Transpose][[2]]-Widths];maxy =Min[Xs\[Transpose][[2]]+Widths];

   If[NumericQ["MeshSize"/.options], dx = Min[(maxx-minx)/25,(maxy-miny)/25, "MeshSize"/.options ]; Print["Mesh Size:", dx],  dx = Min[(maxx-minx)/25,(maxy-miny)/25 ] ];
   rngx=Range[minx,maxx,dx];
   rngy=Range[miny,maxy,dx];
   rngLt=Range[ Length@listCylindricalWaveOtx\[Theta]xr[[1]]];
  listxyWaveOt =  Flatten[ Outer[{#1,#2,Sum[ WaveOtXwtFxXy[[#3,n]][#1,#2] ,{n,Nwaves}]}&,rngx,rngy,rngLt] ,{{3},{1,2}}];
  Return[{{{minx,maxx},{miny,maxy}}, listxyWaveOt}];
];


(* ::Code::Bold:: *)
(*Error Handling*)
ListCylindricalWaveDueToImpulse[Xs_List, rs_Real:0.1, n0_Integer:3, opts:(__?RulesQ):{}]:= Throw[$Failed,
failTag["ListCylindricalWaveDueToImpulse was called with option \"SourceWave\" different to \[ImaginaryI] 0.25 HankelH1[0, k Norm[{x,y}"<>If[List===Head["ImpulsePosition"/.Flatten@List[opts]],
"- (\"ImpulsePosition\"/.options)",""]<>"], which is what all the impulse functions are setup to work with (to be later convoluted with the option \"ImpulseFunction\")."]]/;
( NumericQ@N[("SourceWave"/.Flatten@List[opts])[{1.3,1.2},1.1]] && N[("SourceWave"/.Flatten@List[opts])[{1.3,1.2},1.1]]=!= N[I 0.25 HankelH1[0, 1.1 Norm[{1.3,1.2}-If[List===Head["ImpulsePosition"/.Flatten@List[opts]],("ImpulsePosition"/.Flatten@List[opts]),{0,0}] ] ] ]) ; 
    
(*Main function call*)
ListCylindricalWaveDueToImpulse[Xs_List, rs_Real:0.1, n0_Integer:3, opts:(__?RulesQ):{}]:=Module[{parameters, options =Flatten[ List@opts], 
Maxk,Xi,NN, rng\[Omega]Fourier,listCoeffFk,CoeffsOscatxnx\[Omega],rngr, rng\[Theta], b,bT, listCylindricalWaves, Header, MaxR},
  parameters = {{"MaxFrequency", Maxk,10.},{"ImpulsePosition",Xi,{0,0}},{"MaxFrequencySamples",NN,45} };
  options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts];
  parameters ={{"MaxRadius",MaxR,2rs+ 1.6 Max[Outer[Norm[#1-#2]&,#,#,1,1]]&[{Xi}~Join~Xs]}};
  With[{tmpbT = 2\[Pi] NN/(Maxk(2NN+1))(*T/(2NN+1)*)},
    parameters =parameters~Join~{{"ImpulseFunction",b, If[0<=#<=tmpbT,1,0]&},{"ImpulsePeriod",bT,tmpbT}};   
  ];
  (*tmpbT=1; parameters =parameters~Join~{{"ImpulsePeriod",bT,2tmpbT},{"ImpulseFunction",b,(2 E^(- 1/(1.00000001-(2#/(2tmpbT) -1)^2)) HeavisideTheta[#-.00000001]HeavisideTheta[2tmpbT-#]  &)}};*)
  options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts];
  rng\[Omega]Fourier = rng\[Omega]FourierOffset[Maxk,NN];

(*Test the choosen impulse in a convolution*)
  If["PrintChecks"/.options, Print@ConvolutionTest[b, \[Pi] NN/Maxk (*r*) ,rng\[Omega]Fourier]];
  
(*Scattering coefficients for the scatterers*)
  With[{tmp=MultipleScatteringCoefficients[Xs,rs,n0,Unevaluated@options]},
  (*Add the source wave to the coefficients, where SourceWave will be given by the options in MultipleScatteringCoefficients,
if it was not already specified by the user*)
    listCoeffFk[\[FormalK]_]:= {ReplacePart[ConstantArray[0,2n0+1],n0+1->I 0.25]}~Join~tmp[\[FormalK]];
  ];
  If["PrintMethodOptions"/.options, Print["The chosen options up too calculating Cylindrical Waves: \n",options] ];

  CoeffsOscatxnx\[Omega]= Transpose/@Transpose[listCoeffFk/@rng\[Omega]Fourier];

  (*Choose the range for the time, radius and angle in cylindrical coordinates*)
  rngr=Range[rs,MaxR,rs];
  rng\[Theta]=Range[0,2\[Pi]+.1,.1];

  listCylindricalWaves= Re@CylindricalWaveFromCoefficients[#,rng\[Omega]Fourier,rng\[Theta],rngr,Unevaluated@options]&/@CoeffsOscatxnx\[Omega];
  If[TrueQ["SaveEachWave"/.options],
    Header={Xs,rs,n0,{Maxk,NN}};
    ExportListCylindricalWave[Header,listCylindricalWaves];
  ];
Return[listCylindricalWaves]
];


(* ::Code::Bold:: *)
CylindricalWaveFromCoefficients[listaOnx\[Omega]_,rng\[Omega]Fourier_,rng\[Theta]_:Range[0.,\[Pi],0.2],rngR_:{},opts:(__?RulesQ):{}]:=
Module[{options= If[ Flatten[{opts}]=={},{},Flatten@ List@opts],b,\[Delta],bT,T,bN,rngt,rngtSample,Lrngt,rngr,rngn,maxn, Maxt,nt,
listbO\[Omega],listbN\[Omega]O\[Omega],bolPrintChecks,parameters,parameterNames, NN,t,wIN,\[Delta]\[Omega],listIntegrateSimp,listwOnxr,listwOtx\[Theta]xr},

   (*Select parameters from options or give default value*)
   parameters = {{"PrintChecks", bolPrintChecks,False},{"ImpulseFunction",b , (E^(- 1/(1.00000001-(2#/bT -1)^2)) HeavisideTheta[#-.00000001]HeavisideTheta[bT-#]  &)}
,{"ImpulsePeriod",bT,2.},{"MaxTimeSamples",nt,.5},{"MaxTime",Maxt,6}};
   options = SetParametersAndReturnOptions[parameters,options, Unevaluated@opts];
  (*End parameters setup *)  
   
       \[Delta]\[Omega]=rng\[Omega]Fourier[[2]]-rng\[Omega]Fourier[[1]]; T=(2\[Pi])/\[Delta]\[Omega];
       \[Delta]=First@rng\[Omega]Fourier; If[\[Delta]==0, Print["Error: The frequencyy \[Omega]=0 does not work for Hankel functions.."]; Return[{}]];
       NN= (Length[rng\[Omega]Fourier]-1)/2;

       If[rngR=={}, rngr=Range[0.1, T+0.2, T 4/NN ] ,rngr=rngR];
       If[T< bT, Print["The frequency range is too coarse to cover the Impluse. The max time period covered is at most T=", T, ", while the impulse period is bT=",bT ] ];
       rngt = Range[0,T (2NN)/(2NN+1) ,T/(2NN+1) ];
     
       listbO\[Omega] =T /Sqrt[NN] Fourier[(b/@rngt)E^(I rngt \[Delta])]; 
       maxn=(Length[listaOnx\[Omega]]-1 )/2; rngn=Range[ -maxn,maxn];
       listIntegrateSimp=(\[Delta]\[Omega] listbO\[Omega] # )&/@listaOnx\[Omega]; 

       Lrngt=Min[Length[rngt], Round[Maxt Length[rngt]/ Max[rngt]]];  
       If[IntegerQ[nt],
         If[Lrngt >nt ,rngtSample= rngt[[Round[#]]]&/@ Range[1,Lrngt,Lrngt/nt], rngtSample= rngt[[1;;Lrngt]] ];
       ,
         rngtSample= rngt[[1;;Lrngt]];   
       ];
       listwOnxr=1/(2.\[Pi])Outer[(E^(-I rng\[Omega]Fourier #2)listIntegrateSimp[[#1+maxn+1]]) . HankelH1[#1,#3 rng\[Omega]Fourier ] &, rngn,rngtSample,rngr];
       listwOtx\[Theta]xr = Transpose[(Exp[I rngn #] .listwOnxr)&/@rng\[Theta],{2,1,3}];
       listwOnxr=.;
    Return[Outer[{rngtSample[[#1]],rng\[Theta][[#2]],rngr[[#3]],listwOtx\[Theta]xr[[#1,#2,#3]]}& ,Range@Length@rngtSample,Range@Length@rng\[Theta],Range@Length@rngr ]]
];





(* ::Code::Bold:: *)
TotalWave[Xs_List?(ArrayDepth[#]==2&), listCoeffsFk_,opts:OptionsPattern[]]:= Module[{wave, SourceWave = ("SourceWave"/.Flatten@List@opts), n = (Length[listCoeffsFk[1.][[1]]]-1)/2},
(*,k_?(Head[#]\[Equal] Real||Head[#]==Integer&)*)
  wave[{x_,y_},k_]:=SourceWave[{x,y},k] + Plus@@Map[OutWave[#[[1]],#[[2]],n][{x,y},k]&, Transpose[{Xs,listCoeffsFk[k]}]];
  Return[wave]
]/;("SourceWave"/.Flatten@List[opts])=!= "SourceWave";(*/;(NumericQ@N[("SourceWave"/.Flatten@List[opts])[{1.3,1.2},1.1]]); *)

OutWave[Xj_List?(ArrayDepth[#]==1&),Coeffs_List?(ArrayDepth[#]==1&),n_Integer]:= Module[{f},
    f[X_,k_]:= Coeffs.Array[ HankelH1[#,k Norm[X-Xj] /.Abs->Identity]E^(I # ArcTan@@(X-Xj)) &,2n+1,-Abs@n]; 
  Return@f];

OutWave[Xj_List?(ArrayDepth[#]==1&),aj_?variableQ, n_Integer]:= Module[{f},
    f[X_,k_]:= Plus@@Array[(aj[#] (HankelH1[#,k Norm[X-Xj] /.Abs->Identity]E^(I # ArcTan@@(X-Xj)) ))&,2n+1,-Abs@n]; 
  Return@f];

MultipleScatteringCoefficients[Xs_List:{{0.5,0.5},{-0.5,0.5}}, RS_:0.1, N0_:3, opts:(__?RulesQ):{}]:=
Module[{x,y,a,k,Maxk,parameters, parameterNames, options =Flatten[ List@opts], rngParticles = Range[Length@Xs], 
OutWaves,SourceWave, CoeffEqsFkOn, vars, TCoeffsF, m1, argNSolveFk, coeffsFk, checkDis,CoeffsFromScatteredWave,CoeffsFromScatteredWavePre},
   Protect[x,y,a, k,rngParticles];
  
   (*We check if options sets any of the parameters, if not we assign some default values.*)
   parameters = {{"MaxFrequency", Maxk, 10.}}; 
   options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts]; 
   If[NumericQ[("SourceWave"/.options)[{1,1},1]], 
     SourceWave[{x_,y_},k_] = ("SourceWave"/.options)[{x,y},k]
   , 
     parameters = {{"SourceWave", SourceWave, I/4 HankelH1[0,#2 Norm[{#1[[1]],#1[[2]]}] /.Abs->Identity]&}} ;   
     options = SetParametersAndReturnOptions[parameters,options,Unevaluated@opts]; 
   ]; 

   If[Length[Xs]>1,checkDis=Union[Flatten@Outer[N@Norm[#1-#2]&, Xs, Xs,1,1]][[2]]];
   (*"CheckDistance"->Union[Flatten@Outer[N@Norm[#1-#2]&, Xs~Join~{Xi}, Xs~Join~{Xi},1,1]][[2]];*)
   TCoeffsF = ScatteringCoefficientOperator[N0,RS,Sequence@@options,"CheckDistance"->checkDis];
(*This block does as much precalculations as it can for a the contribution to the scattering coefficient due to another scattered wave*)
   Clear@@Evaluate[("h"<>ToString[#])&/@Range[0,N0]];
   Block[{b,\[Kappa],xi,xs,dx,yi,ys,dy,dist,ksqrdist,e\[Theta],m,subPreCalculate,hvars=  ToExpression["h"<>ToString[#]]&/@Range[0,N0]},
       subPreCalculate[{dx_,dy_,dist_,ksqrdist_,e\[Theta]_},Thread[Pattern[Evaluate[ hvars],Blank[]]],b_,\[Kappa]_]= {(xi-xs)^2+(yi-ys)^2-> dist,
      (-xi+xs)^2+(-yi+ys)^2-> dist,dist^(1/2)-> ksqrdist /\[Kappa],xi-> dx + xs,yi-> dy + ys,Power[E,Times[ Complex[0, m_],ArcTan[-dx,-dy]] ]-> e\[Theta]^m }
      ~Join~Thread[Rule[HankelH1[Range[0,N0],ksqrdist],hvars]];
       CoeffsFromScatteredWavePre[{dx_,dy_,dist_,ksqrdist_,e\[Theta]_},Thread[Pattern[Evaluate[ hvars],Blank[]]],b_,\[Kappa]_]= Chop@Together[
      TCoeffsF[OutWave[{xi,yi},b,N0],{xs,ys},\[Kappa]]//.subPreCalculate[{dx,dy,dist,ksqrdist,e\[Theta]},hvars,b,\[Kappa]]]//.subPreCalculate[{dx,dy,dist,ksqrdist,e\[Theta]},hvars,b,\[Kappa]];
       CoeffsFromScatteredWave[{xi_,yi_},{xs_,ys_},b_,\[Kappa]_]:=CoeffsFromScatteredWavePre[{xi-xs,yi-ys,(xi-xs)^2+(yi-ys)^2,\[Kappa] ((xi-xs)^2+(yi-ys)^2)^(1/2),
      Exp[I ArcTan[-xi+xs,-yi+ys]]},HankelH1[Range[0,N0],\[Kappa] ((xi-xs)^2+(yi-ys)^2)^(1/2)],b,\[Kappa]];
   ];

   Clear@@Evaluate[("h"<>ToString[#])&/@Range[0,N0]];

   Protect[SourceWave,CoeffsFromScatteredWave];

   vars=Array[a[#],2 N0+1,-Abs@N0]&/@rngParticles;
   CoeffEqsFkOn[k_]= Chop@Together@Map[TCoeffsF[SourceWave,Xs[[#]],k]-vars[[#]]+Plus@@First@Outer[CoeffsFromScatteredWave[Xs[[#2]],Xs[[#1]],a[#2],k ]&
      ,{#},Drop[rngParticles,{#}] ,1] &,rngParticles];

   argNSolveFk[k_]={Thread[Flatten[CoeffEqsFkOn[k],1]==0], Flatten[vars]};
   coeffsFk[k_]:=
      Module[{subNvars},
        If[k==0, Print["Frequency k =0 has been replaced with k=", k=10^(-4$MachinePrecision)] ];
        subNvars =Flatten[NSolve@@argNSolveFk[k]];
        If[!NumericQ[Total[vars/.subNvars,3] || Length[subNvars]!= Length[Flatten@vars] ], 
           Print["Non-numeric/unique solution for MultipleScatteringCoefficients found for k=",k,"(< ",Maxk,"?), solution found: \n", subNvars];
           Return[ vars 0 ] 
        ];
        Return[vars/.subNvars]
      ];
  
  If[("PrintChecks"/.options)==True,
     Module[{rngk,rng\[Theta], SourceWaveNorms, BoundaryNorms,BoundaryNorms2, rng3Particles, PlotStyleBlended, totalWave= TotalWave[Xs,coeffsFk,"SourceWave"-> SourceWave]},
        rngk =Range[RS, Maxk, (Maxk-RS)/(2N0+1)];
        rng\[Theta] =Range[0.,2.\[Pi], 2\[Pi]/(N0+2)];
        rng3Particles = If[Length@rngParticles >3, rngParticles[[{1,Round[Length[rngParticles]/3],Length[rngParticles]}]] , rngParticles];
       If[("BoundaryCondition"/.options)=== "Neumann",
           Needs["NumericalCalculus`"];
          SourceWaveNorms= Abs@Apply[Plus,Outer[ND[SourceWave[Xs[[#2]]+r{Cos[#1],Sin[#1]},#3],r,RS] &, rng\[Theta], rng3Particles,rngk]];
           BoundaryNorms= Abs@Apply[Plus,Outer[ND[totalWave[Xs[[#2]]+r{Cos[#1],Sin[#1]},#3],r,RS] & , rng\[Theta], rng3Particles,rngk ]];
        ,
           SourceWaveNorms= Abs@Apply[Plus,Outer[SourceWave[Xs[[#2]]+RS{Cos[#1],Sin[#1]},#3] &, rng\[Theta], rng3Particles,rngk]];
           BoundaryNorms= Abs@Apply[Plus,Outer[totalWave[Xs[[#2]]+RS{Cos[#1],Sin[#1]},#3]& , rng\[Theta], rng3Particles,rngk ]]; 

       ];
        PlotStyleBlended = Table[{Blend[{Blue,Red},x],Thick},{x,0,1, 1/(Length[Xs]- If[Length[Xs]>1,1,0])}];
        Print@ListPlot[Transpose[{rngk,#}]&/@(BoundaryNorms/SourceWaveNorms),
                   Joined-> True,InterpolationOrder-> 2, PlotRange-> {0,0.01},
                   AxesLabel-> {"Freq. k","Error"}, PlotLabel-> "Check ||Boundary[k]|| / ||IncidentWave[k]|| for up to 3 scatterers", 
                   PlotLegends->Evaluate[("X="<>ToString[#])&/@Xs], Filling->Axis
         ];
    ];
  ];
Return[coeffsFk];
];

ScatteringCoefficientOperator[n0_Integer,rs_,opts:(__?RulesQ):{}]:=Module[{p,Dp, CC,TCoeffF,Xi, rng\[Theta], rngk,k,r,maxk, options =Flatten[ List@opts],
wEs,wEN,w0k,ErrFkx\[Theta],rngErr,\[Theta],plotLabel},
   CC[m_,j_,n_]:= Sum[(-1)^((-n)/2 +j-m1) Binomial[m,m1]Binomial[2 j-m,-n/2+j-m1]//N,{m1, Max[0,-n/2-j+m],Min[m,-n/2+j] }];
   If[ ("BoundaryCondition"/.options)=== "Neumann",
      Dp[n_,f_,X_,k_,rsOrder_]:=Sum[ (D[f[{x,y},k],{x,m},{y,2j-m}] /.{x-> X[[1]],y-> X[[2]]}) (CC[m,j,n](I^(m-2 j )) )/(m!(2 j -m)!) j (rs/2)^(2j-1)
         ,{j, If[n==0, 1, Abs[n]/2],(rsOrder +1)/2},{m,0,2 j}]; 
      TCoeffF[f_,X_,k_]:= - (2/k Dp[#,f,X,k,n0]/(HankelH1[-1+#,k rs]-HankelH1[1+#,k rs]))&/@Range[-n0,n0]
   ,
      p[n_,f_,X_,k_,rsOrder_]:=Sum[(D[f[{x,y},k],{x,m},{y,2j-m}] /.{x-> X[[1]],y-> X[[2]]}) (CC[m,j,n](I^(m-2 j )) )/(m!(2 j -m)!) (rs/2)^(2j)
         ,{j,Abs[n]/2,(rsOrder +1)/2},{m,0,2j}]; 
      TCoeffF[f_,X_,k_]:= - (p[#,f,X,k,n0]/HankelH1[#, k rs])&/@Range[-n0,n0]
   ]; 

   If[("PrintChecks"/.options)=!=False && ("BoundaryCondition"/.options)=== "Dirichlet",
      If[NumericQ["MaxFrequency"/.options],maxk="MaxFrequency"/.options, maxk=1.2/rs];
      If[NumericQ["CheckDistance"/.options],Xi ={0, "CheckDistance"/.options} , Xi= {5 rs, 5 rs} ];      
      If[("BoundaryCondition"/.options)=== "Neumann",
        wEs= HankelH1[0, #2 Norm[{#[[1]],#[[2]]}-Xi]/.Abs->Identity ]+HankelH1[1, #2 Norm[{#[[1]],#[[2]]}-Xi]/.Abs->Identity ]&,
        wEs= HankelH1[0, #2 Norm[{#[[1]],#[[2]]}-Xi]/.Abs->Identity ]&
      ];

      rng\[Theta]={0,2 \[Pi],2\[Pi]/4.9};
      rngk=Range[rs, maxk,(maxk-rs)/21 ];
      wEN[\[Theta]_,k_]=  Sum[ p[n,wEs,{0,0},k,n0]E^(I n \[Theta]),{n,-n0,n0}];
      w0k= Abs[Plus@@(wEs[rs{Cos[#],Sin[#] },maxk/2.]&/@rng\[Theta] )];
      Print["|Wi[k]\!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\): ",w0k];
      ErrFkx\[Theta] = Abs[wEs[rs{Cos[#2],Sin[#2] },#1]-wEN[#2, #1]]&;
      (*Print["p[n,wEs,{0,0},k,n0,rs]: ", p[1,wEs,{0,0},0.5,n0]];*)

      rngErr=(Plus@@#/w0k)&/@Outer[ErrFkx\[Theta],rngk,rng\[Theta]];
      plotLabel ="W= \!\(\*SuperscriptBox[\(H\), \((0)\)]\)(k rs)"<>If[("BoundaryCondition"/.options)=== "Neumann","+\!\(\*SuperscriptBox[\(H\), \((1)\)]\)(k rs)",""]<> "with "<>ToString[n0]<>" Fourier components, rs: "<>ToString[rs];
      Print@ListPlot[{rngk,rngErr}\[Transpose],PlotRange->All,Joined->True,InterpolationOrder->2,Filling->Bottom,
                  AxesLabel->{"k","|W[k]-  \!\(\*SubscriptBox[\(W\), \(n\)]\)[k] \*SuperscriptBox[\(\[ExponentialE]\), \(\[ImaginaryI] n \[Theta]\)] \!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\)/ |W[k]\!\(\*SubscriptBox[\(|\), \(\[Theta]\)]\) "},
                  PlotLabel->plotLabel,
                  PlotLegends-> "Min Distance="<>ToString[Norm@Xi]<>"- rs"
            ];
   ]; 
   Return[TCoeffF]
];

SingleScattererCoefficientsFromImpulse[rng\[Omega]Fourier_,n_:2, RS_:0.1 ,RI_:1.8]:=Module[{rngn= Range[-n,n]},
   (*If[First[rng\[Omega]]!=- Last[rng\[Omega]], Print["The frequency range rng\[Omega] has to satisfy First[rng\[Omega]] == - Last[rng\[Omega]]"]; Return[{{0}}]  ];*)
   (*change the frequency \[Omega]=0 to \[Omega] = \[Delta]\[Omega]/20. *)   
   Return[ -(I/4)Outer[  (BesselJ[-1+#1,#2 RS]-BesselJ[1+#1,#2 RS])/(HankelH1[-1+#1,#2 RS]-HankelH1[1+#1,#2 RS]) HankelH1[#1,#2 RI] & ,rngn, rng\[Omega]Fourier] ]; 
];

ConvolutionTest[b_,r0_,rng\[Omega]Fourier_]:= Module[{g,g\[Omega],listbO\[Omega],listwOt,listwOt2,listgOt,rngt,wFt,lpt,pt, T,
NN=(Length[rng\[Omega]Fourier]-1)/2, \[Delta]=First@rng\[Omega]Fourier},
(*Here we check our discrete Fourier transform for the convolution of the Impulse with a Hankel function *)
         g\[Omega] = I/4 HankelH1[0,# r0] &;
         g =1/(2\[Pi]) HeavisideTheta[#-r0](#^2-r0^2)^(-.5) &; 
         T= 2\[Pi] NN /Max[rng\[Omega]Fourier];
         rngt = Range[0,T 2NN/(2NN+1) ,T/(2NN+1)];         
         listbO\[Omega] = T/Sqrt[2NN+1] Fourier[(b/@rngt)E^(I rngt \[Delta])];
         listgOt = (Sqrt[2NN+1]/T) Re[ InverseFourier[ g\[Omega][rng\[Omega]Fourier]]E^(-I rngt \[Delta])];

         listwOt2= T/(2NN+1)(listgOt.(b/@(#-rngt)))&/@rngt;
         listwOt= Sqrt[2NN+1]/T InverseFourier[listbO\[Omega] g\[Omega][rng\[Omega]Fourier]] E^(-I rngt \[Delta]); 
         wFt[t_]:=  NIntegrate[g[\[Tau]]b[t-\[Tau]],{\[Tau],r0,Max[r0,t]}];       
         lpt=ListPlot[{{rngt,Re@listwOt}\[Transpose](*,{rngt,Re@listwOt2}\[Transpose]*)},PlotStyle->{{Orange,PointSize[Large]},{Red,PointSize[Large]}},PlotLegends->{"Frequency sampled \!\(\*SubsuperscriptBox[\(w\), \(N\), \(n\)]\)","Time sampled \!\(\*SubsuperscriptBox[\(w\), \(N\), \(n\)]\)"},PlotRange-> All];     
         pt =Plot[Re[wFt[t]],{t,0,T},PlotRange-> All, PlotLabel->"Convolution of Impluse b(t) with 2D Green's function.",PlotLegends->{"analytic w[t]"},AxesLabel->{"t"}];    

  Return@Show[pt,lpt,PlotRange->All];
];

End[];
EndPackage[];

