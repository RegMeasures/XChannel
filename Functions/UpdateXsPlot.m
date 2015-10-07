function UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T)
% Update bed level and water level on cross-section plot
% See also PlotXS

% Update bed level
set(XsFigure.BedLineH, 'YData', [Cell.Z; Cell.Z(end)])

% Update bank erosion
IsBank = Edge.IsBank~=0;
BankIDPlotX = [Edge.N(IsBank)',NaN; 
               Edge.N(IsBank)',NaN;
               NaN(1,sum(IsBank)+1)];
BankIDPlotX = BankIDPlotX(:);
BankIDPlotY = [Cell.Z(IsBank(1:end-1))',NaN; 
               Cell.Z(IsBank(2:end))',NaN;
               NaN(1,sum(IsBank)+1)];
BankIDPlotY = BankIDPlotY(:);
set(XsFigure.BankIDLineH, 'XData', BankIDPlotX, 'YData', BankIDPlotY)

StencilPlotX = [Cell.N(Bank.Top)',NaN;
                Edge.N(IsBank)',NaN;
                Cell.N(Bank.Bottom)',NaN;
                NaN(1,Bank.NBanks+1)];
StencilPlotX = StencilPlotX(:);
StencilPlotY = [Cell.Z(Bank.Top)',NaN;
                (Cell.Z(IsBank(1:end-1)) + Cell.Z(IsBank(2:end)))' ./2,NaN;
                Cell.Z(Bank.Bottom)',NaN;
                NaN(1,Bank.NBanks+1)];
StencilPlotY = StencilPlotY(:);
set(XsFigure.BankStencilH, 'XData', StencilPlotX, 'YData', StencilPlotY)

% Update velocity
set(XsFigure.VelLineH, 'YData', Cell.U)

% Update water level
WlForPlot = WL * ones(size(Cell.U));
WlForPlot(~Cell.Wet) = NaN;
set(XsFigure.WlLineH, 'YData', [WlForPlot; WlForPlot(end)]);

% Update Shear stress
set(XsFigure.ShearLineH, 'YData', Cell.Tau_S)


% Update Transport rate
set(XsFigure.TransLineH, 'YData', Cell.qsS_flow_kg)

% Update grain size
set(XsFigure.DgLineH, 'YData', 1000*Cell.Dg_m)

% Update armour index
ArmourIndex = (Cell.Dg_m)./(Cell.SubDg_m);
set(XsFigure.ArmourLineH, 'YData', ArmourIndex)

% Update title
set(XsFigure.XsPlotTitleH, 'String', ['Model time = ', num2str(T), ' s'])   

drawnow

end
