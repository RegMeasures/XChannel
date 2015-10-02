function [XsFigure] = PlotXS(CellCentrePos,CellEdgePos,IniBedLevel,BedLevel,WL,WetCells,VAvVel,Tau_S,qsS_flow_kg,Dg_m,SubDg_m,T)
% Plot cross-section

% Create new figure in the middle of the screen with a wide aspect ratio
scrsz = get(groot,'ScreenSize');
XsFigure.FigureH = figure('Position',[scrsz(3)/2 - 500, scrsz(4)/2-400, 1000, 800]);


%% Cross-section shape plot (top panel)
XsFigure.XsAxesH = subplot(4,1,[1:2]);

% Plot initial bed level
stairs(CellEdgePos, [IniBedLevel; IniBedLevel(end)],'k--');
hold on

% Plot current bed level and velocity
[XsFigure.BedVelPlotH, XsFigure.BedLineH, XsFigure.VelLineH] = plotyy(CellEdgePos, [BedLevel; BedLevel(end)], CellCentrePos, VAvVel, 'stairs', 'plot' );

% Plot water level
WlForPlot = WL * ones(size(CellCentrePos));
WlForPlot(~WetCells) = NaN;
XsFigure.WlLineH = stairs(CellEdgePos, [WlForPlot; WlForPlot(end)],'b-','LineWidth',1);

% Labels and formatting etc
set(XsFigure.XsAxesH,'XTickLabel','')
ylabel(XsFigure.BedVelPlotH(1),'Elevation (m)');
ylabel(XsFigure.BedVelPlotH(2),'Velocity (m/s)');
%legend('Initial bed','Current bed','Water level',...
%       'Location','SouthOutside','Orientation','horizontal');
%XScale = get(XsAxesH,'XLim');
XScale = CellEdgePos([1,end]);
set(XsFigure.BedVelPlotH(1),'XLim', XScale,'YColor','k');
set(XsFigure.BedVelPlotH(2),'XLim', XScale);
set(XsFigure.BedLineH,'Color','k','LineWidth',1)
XsFigure.XsPlotTitleH = title(['Model time = ', num2str(T), ' s']);

hold off

%% Shear stress and transport plot (2nd panel)
subplot(4,1,3);

% Plot stress and transport
[XsFigure.ShearTransPlotH, XsFigure.ShearLineH, XsFigure.TransLineH] = plotyy(CellCentrePos, Tau_S, CellCentrePos, qsS_flow_kg);

% Labels and formatting
set(XsFigure.ShearTransPlotH(1),'XLim',XScale,'XTickLabel','','YLimMode','auto','YTickMode','auto','box','off')
set(XsFigure.ShearTransPlotH(2),'XLim',XScale,'XTickLabel','','YLimMode','auto','YTickMode','auto')
ylabel(XsFigure.ShearTransPlotH(1),'Shear stress (kgm^-^1s^-^2)');
ylabel(XsFigure.ShearTransPlotH(2),'Transport rate (kg/s/m)');

%% Grain size plot (3rd panel)
subplot(4,1,4)

% Plot substrate initial condition
plot(CellCentrePos, SubDg_m*1000,'b--');
hold on

% calculate armour index
ArmourIndex = (Dg_m)./(SubDg_m);

% Plot active layer grain size & armour index

[XsFigure.DgPlotH, XsFigure.DgLineH, XsFigure.ArmourLineH] = plotyy(CellCentrePos, Dg_m*1000, CellCentrePos, ArmourIndex, 'semilogy', 'plot');

% Labels and formatting
%set(XsFigure.DgAxesH(1),'XLim',XScale)
xlabel('Distance across section (m)');
set(XsFigure.DgPlotH(1),'XLim',XScale,'YLimMode','auto','YTickMode','auto','box','off')
set(XsFigure.DgPlotH(2),'XLim',XScale,'YLimMode','auto','YTickMode','auto','box','off')
ylabel(XsFigure.DgPlotH(1),'D_g (mm)');
ylabel(XsFigure.DgPlotH(2),'Armour index');
hold off

end
