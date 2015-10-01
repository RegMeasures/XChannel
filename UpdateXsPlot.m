function UpdateXsPlot(BedPlotH, WlPlotH, XsPlotTitleH, BedLevel,WL,DryFlc,T)
% Update bed level and water level on cross-section plot
% See also PlotXS


% Update bed level
set(BedPlotH, 'YData', [BedLevel; BedLevel(end)])

% Update water level
WlForPlot = WL * ((WL - BedLevel) > DryFlc);
WlForPlot(WlForPlot==0) = NaN;
set(WlPlotH, 'YData', WlForPlot);

% title(['Model time = ', num2str(T), ' s'])
set(XsPlotTitleH, 'String', ['Model time = ', num2str(T), ' s'])   

drawnow

end
