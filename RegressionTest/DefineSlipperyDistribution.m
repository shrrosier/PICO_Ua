function  [UserVar,C,m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

m=3; 

%  u=C*tau^m
%  tau=20 kPa ; u=1 m/d -> C=(1/20)^(1/m)

C=(1/20)^(1/m)+MUA.coordinates(:,1)*0;

C=C*1e-10;

end
