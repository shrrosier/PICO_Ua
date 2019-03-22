
function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FTBD)


% Defines model geometry

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
alpha=0.;

gamma=0.01;
S=x*0 ; 
B=-500 + zeros(MUA.Nnodes,1) ;
Circle = ((x/1000).^2 + (y/1000).^2)<=(500)^2; %circle around (0,0), 500km radius
B(Circle) = -1000;
b=B;

s=b+1000;


end
