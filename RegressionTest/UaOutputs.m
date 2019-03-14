
function  UserVar=UaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);

v2struct(F);

time=CtrlVar.time; 


plots='-ubvb-e-save-mua-dhdt(x)-h(x)-';
%plots='-udvd-ubvb-ub(x)-sbSB(x)-txzb(x)-';
%plots='-ub(x)-h(x)-sbSB(x)-';
%plots='-h(x)-ub(x)-dhdt(x)-';

TRI=[];
x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);

if ~isempty(strfind(plots,'-save-'))

    % save data in files with running names
    % check if folder 'ResultsFiles' exists, if not create

    if strcmp(CtrlVar.UaOutputsInfostring,'First call ') && exist('ResultsFiles','dir')~=7 ;
        mkdir('./ResultsFiles') ;
    end
    
    if strcmp(CtrlVar.UaOutputsInfostring,'Last call')==0

        FileName=['./ResultsFiles/',sprintf('%07i',CtrlVar.UaOutputsCounter),'-',CtrlVar.Experiment];
        
        fprintf(' Saving data in %s \n',FileName)
        save(FileName,'CtrlVar','MUA','time','F','GF')
        
    end
end

% % only do plots at end of run
% if ~strcmp(CtrlVar.UaOutputsInfostring,'Last call') ; return ; end

[~,I]=sort(x) ;

if ~isempty(strfind(plots,'-txzb(x)-'))
    
    [txzb,tyzb]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,AGlen,n,C,m,GF,s,b,ub,vb);
    
    figure ;  plot(x/CtrlVar.PlotXYscale,txzb) ; title('txzb(x)')
    
end


if ~isempty(strfind(plots,'-ub(x)-'))
    figure
    plot(x(I)/CtrlVar.PlotXYscale,ub(I)) ;
    title(sprintf('u_b(x) at t=%-g ',time)) ; xlabel('x') ; ylabel('u_b')
    drawnow
end


if ~isempty(strfind(plots,'-dhdt(x)-'))
    figure
    plot(x(I)/CtrlVar.PlotXYscale,dhdt(I)) ;
    title(sprintf('dhdt(x) at t=%-g ',time)) ; xlabel('x') ; ylabel('dh/dt')
    %drawnow
    %export_fig(strcat('./Outputs/', 'dhdt_at_time', num2str(time)), '-r400');
end


if ~isempty(strfind(plots,'-h(x)-'))
    figure;
    plotyy(x(I)/CtrlVar.PlotXYscale,h(I),x(I)/CtrlVar.PlotXYscale,GF.node(I)) ;
    
    if CtrlVar.Implicituvh
        title(sprintf('fully-implicit h(x) at t=%-g (%s)',time,CtrlVar.uvhTimeSteppingMethod)) ;
    else
        title(sprintf('semi-implicit h(x) at t=%-g (TG3=%i)',time,CtrlVar.TG3)) ;
    end
    xlabel('x') ; ylabel('h')
    drawnow
end

if ~isempty(strfind(plots,'-ud(x)-'))
    figure
   plot(x/CtrlVar.PlotXYscale,ud) ;
    title(sprintf('u_d(x) at t=%-g ',time)) ; xlabel('x') ; ylabel('u_d')
end


if ~isempty(strfind(plots,'-sbSB(x)-'))
    figure
    
    plot(x(I)/CtrlVar.PlotXYscale,S(I),'k--') ; hold on
    plot(x(I)/CtrlVar.PlotXYscale,B(I),'k') ; 
    plot(x(I)/CtrlVar.PlotXYscale,b(I),'b') ; 
    plot(x(I)/CtrlVar.PlotXYscale,s(I),'b') ;
    
    title(sprintf('sbSB(x) at t=%-g ',time)) ; xlabel('x') ; ylabel('z')
    drawnow
end


if ~isempty(strfind(plots,'-sbB-'))
    figure(5)
    hold off
    if isempty(TRI) ;  TRI = delaunay(x,y); end
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ; hold on
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'EdgeColor','none') ;
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
    view(50,20); lightangle(-45,30) ; lighting phong ;
    xlabel('y') ; ylabel('x') ;
    colorbar ; title(colorbar,'(m)')
    hold on
    
    title(sprintf('sbB at t=%#5.1g ',time))
    axis equal ; tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale]); axis tight
    hold off
end


if ~isempty(strfind(plots,'-ubvb-'))
    % plotting horizontal velocities
    figure
    N=1;

    QuiverColorGHG(x(1:N:end),y(1:N:end),ub(1:N:end),vb(1:N:end),CtrlVar);
    hold on
    title(sprintf('(ub,vb) t=%-g ',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
    axis equal tight
    %export_fig(strcat('./Outputs/', 'ubvb_at_time', num2str(time)), '-r400');
    
end

if ~isempty(strfind(plots,'-udvd-'))
    % plotting horizontal velocities
    figure
    N=1;

    QuiverColorGHG(x(1:N:end),y(1:N:end),ud(1:N:end),vd(1:N:end),CtrlVar);
    hold on
    title(sprintf('(ud,vd) t=%-g ',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
    axis equal tight
    
end

if ~isempty(strfind(plots,'-e-'))
    % plotting effectiv strain rates
    
    % first get effective strain rates, e :
    [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
    % all these variables are are element variables defined on integration points
    % therfore if plotting on nodes, must first project these onto nodes
    eNod=ProjectFintOntoNodes(MUA,e);
    
    figure
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,eNod,CtrlVar)    ;
    title(sprintf('e t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)')
    %export_fig(strcat('./Outputs/', 'strainrates_at_time', num2str(time)), '-r400');
    
end

if ~isempty(strfind(plots,'-ub-'))
    
    figure
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,ub,CtrlVar)    ;
    title(sprintf('ub t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)')
    
end

if ~isempty(strfind(plots,'-mua-'))

    figure;
    PlotMuaMesh(CtrlVar,MUA)
    hold on 
  
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r','LineWidth',2);
    title(sprintf('t=%g',time))
    hold off
    %export_fig(strcat('./Outputs/', 'MUA_at_time', num2str(time)), '-r400');
end




end
