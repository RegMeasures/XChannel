function [XsFigureH, BedPlotH, WlPlotH, XsPlotTitleH] = PlotXS(CellCentrePos,CellEdgePos,IniBedLevel,BedLevel,WL,DryFlc,T)
% Plot cross-section

% Create new figure in the middle of the screen with a wide aspect ratio
scrsz = get(groot,'ScreenSize');
XsFigureH = figure('Position',[scrsz(3)/2 - 500, scrsz(4)/2-250, 1000, 500]);
hold on

% Plot initial bed level
stairs(CellEdgePos, [IniBedLevel; IniBedLevel(end)],'k-');

% Plot current bed level
BedPlotH = stairs(CellEdgePos, [BedLevel; BedLevel(end)],'g-');

% Plot water level
WlForPlot = WL * ((WL - BedLevel) > DryFlc);
WlForPlot(WlForPlot==0) = NaN;
WlPlotH = plot(CellCentrePos, WlForPlot,'b-');

% Label axes
xlabel('Distance across channel (m)');
ylabel('Elevation (m)');

% Add legend and title
legend('Initial bed','Current bed','Water level',...
       'Location','SouthOutside','Orientation','horizontal');
XsPlotTitleH = title(['Model time = ', num2str(T), ' s']);
   
hold off
end
