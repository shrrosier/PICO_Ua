
function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar)

% to plot the SSA terms:
CtrlVar.MUA.MassMatrix=true; 

CtrlVar.Experiment='RegressionTest';
CtrlVar.doplots=0; CtrlVar.doRemeshPlots=0;
UserVar.OutputsDir = 'ResultsFiles';

xd=600e3; xu=0e3 ; yl=600e3 ; yr=-600e3; % 12km wide
MeshBoundaryCoordinates=flipud([xu yr ; xd yr ; xd yl ; xu yl]);

save MeshBoundaryCoordinates.mat MeshBoundaryCoordinates;

%% Types of runs
CtrlVar.TimeDependentRun=0;
CtrlVar.time=0 ;
CtrlVar.dt=0.01;
CtrlVar.ResetTime = 1;
CtrlVar.TotalNumberOfForwardRunSteps=1;
CtrlVar.AdaptiveTimeStepping=1 ;
CtrlVar.TotalTime=0;
CtrlVar.ThicknessConstraints=0;
CtrlVar.FlowApproximation='SSTREAM' ;  % 'SSTREAM'|'SSHEET'|'Hybrid'

%make sure that width of grounding line is well below 1m which is the boundary layer width that we are interested in 
CtrlVar.kH=10; % 0.001 m width

%CtrlVar.Implicituvh=0;      CtrlVar.TG3=0;
%CtrlVar.uvhTimeSteppingMethod='theta';  % theta | tg3 | supg
%CtrlVar.uvhTimeSteppingMethod='supg';  % theta | tg3 | supg
%CtrlVar.SUPG.beta0=0.5 ; CtrlVar.SUPG.beta1=0.0 ;
%CtrlVar.theta=0.5;

%CtrlVar.uvhTimeSteppingMethod='tg3';  CtrlVar.TG3=1 ; % theta | tg3 | supg
%CtrlVar.uvhTimeSteppingMethod='shocks';


%CtrlVar.SpeedZero=1e-10;
%% Solver
CtrlVar.NLtol=1e-15; % this is the square of the error, i.e. not root-mean-square error
CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.InfoLevel=10;
CtrlVar.LineSeachAllowedToUseExtrapolation=1;

%% Restart
CtrlVar.Restart=0; 
CtrlVar.AdaptMesh=0;

CtrlVar.WriteRestartFile=1;
CtrlVar.NameOfRestartFiletoRead='./ResultsFiles/Ua2D_Restartfile.mat';
CtrlVar.NameOfRestartFiletoWrite='./ResultsFiles/Ua2D_Restartfile.mat';

%% Mesh generation and remeshing parameters

CtrlVar.ReadInitialMesh=0;    % if true then read FE mesh (coordinates, connectivity) directly from a .mat file
% unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='NewMeshFile.mat';

CtrlVar.TriNodes=3 ;
CtrlVar.MeshSize=7e3; %
CtrlVar.MeshSizeMin=7e3;% 
CtrlVar.MeshSizeMax=CtrlVar.MeshSize;
CtrlVar.MaxNumberOfElements=10000; %

%% for adaptive meshing
CtrlVar.MeshGenerator='mesh2d';  % possible values: {mesh2d|gmsh}

CtrlVar.AdaptMeshInitial=1  ; % remesh in first run-step irrespecitivy of the value of AdaptMeshInterval
CtrlVar.AdaptMeshRunStepInterval = 1; % Number of run-steps between mesh adaptation
CtrlVar.AdaptMeshMaxIterations=4;  % Number of adapt mesh iterations within each run-step.
CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan=0;

CtrlVar.InfoLevelAdaptiveMeshing=1;
CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';


CtrlVar.MeshAdapt.GLrange=[1500 500 ; 1000 250 ; 500 100 ; 200 50; 100 20 ]; %
                                            
%% plotting
CtrlVar.doplots=0;
CtrlVar.PlotLabels=0 ; CtrlVar.PlotMesh=0; CtrlVar.PlotBCs=0;
CtrlVar.PlotXYscale=1000;     % used to scale x and y axis of some of the figures, only used for plotting purposes



end
