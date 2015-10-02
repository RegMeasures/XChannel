function UpdateXsPlot(XsFigure, BedLevel,WL,WetCells,VAvVel,Tau_S,qsS_flow_kg,Dg_m,SubDg_m,T)
% Update bed level and water level on cross-section plot
% See also PlotXS

% Update bed level
set(XsFigure.BedLineH, 'YData', [BedLevel; BedLevel(end)])

% Update velocity
set(XsFigure.VelLineH, 'YData', VAvVel)

% Update water level
WlForPlot = WL * ones(size(VAvVel));
WlForPlot(~WetCells) = NaN;
set(XsFigure.WlLineH, 'YData', [WlForPlot; WlForPlot(end)]);

% Update Shear stress
set(XsFigure.ShearLineH, 'YData', Tau_S)


% Update Transport rate
set(XsFigure.TransLineH, 'YData', qsS_flow_kg)

% Update grain size
set(XsFigure.DgLineH, 'YData', 1000*Dg_m)

% Update armour index
ArmourIndex = (Dg_m)./(SubDg_m);
set(XsFigure.ArmourLineH, 'YData', ArmourIndex)

% Update title
set(XsFigure.XsPlotTitleH, 'String', ['Model time = ', num2str(T), ' s'])   

drawnow

end
