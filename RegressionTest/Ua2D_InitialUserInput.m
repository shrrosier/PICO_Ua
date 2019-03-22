
function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar)


CtrlVar.Experiment='RegressionTest';
CtrlVar.doplots=0; CtrlVar.doRemeshPlots=0;
UserVar.OutputsDir = 'ResultsFiles';

% Quadratic area with a round ice shelf its the center
xd=600e3; xu=-600e3 ; yl=600e3 ; yr=-600e3; 
MeshBoundaryCoordinates=flipud([xu yr ; xd yr ; xd yl ; xu yl]);
xd_inner = 1e3; xu_inner = -1e3; yl_inner = 1e3; yr_inner = -1e3;
MeshBoundaryCoordinates = [MeshBoundaryCoordinates; NaN NaN; flipud([xu_inner yr_inner ; xd_inner yr_inner ; xd_inner yl_inner ; xu_inner yl_inner]);]

save MeshBoundaryCoordinates.mat MeshBoundaryCoordinates;

%% 
CtrlVar.TimeDependentRun=0;
CtrlVar.time=0 ;
CtrlVar.dt=0.01;
CtrlVar.ResetTime = 1;
CtrlVar.TotalNumberOfForwardRunSteps=1;
CtrlVar.AdaptiveTimeStepping=1 ;
CtrlVar.TotalTime=0;
CtrlVar.ThicknessConstraints=0;
CtrlVar.FlowApproximation='SSTREAM' ;  % 'SSTREAM'|'SSHEET'|'Hybrid'

CtrlVar.kH=10; 

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

CtrlVar.ReadInitialMesh=0;    
CtrlVar.ReadInitialMeshFileName='NewMeshFile.mat';

CtrlVar.TriNodes=3 ;
CtrlVar.MeshSize=10e3; %
CtrlVar.MeshSizeMin=10e3;% 
CtrlVar.MeshSizeMax=CtrlVar.MeshSize;
CtrlVar.MaxNumberOfElements=1e6; %

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
