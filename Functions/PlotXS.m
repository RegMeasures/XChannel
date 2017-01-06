function [XsFigure] = PlotXS(Cell, Edge, Bank, WL, T, Flow, PlotSed)
%PLOTXS   Plot cross-section shape, velocity, WL, transport, etc.
%Create updateable figure with 3 axes: top plot shows initial, current and
%final cross-section shape; middle plot shows shear stress and transport
%rate; and bottom panel (optional) shows grain size. Generated figure is
%subsequently updated/animated in XChannel model by the UpdateXsPlot
%function.
%
%   [XsFigure] = PLOTXS(Cell, Edge, Bank, WL, T, Flow, PlotSed)
%
%   Inputs:
%      Cell    = Struct of cell center properties initialised by
%                InitialiseVariables and set in earlier steps of XChannel
%      Edge    = Struct of cell edge properties initialised by
%                InitialiseVariables and set in earlier steps of XChannel
%      Bank    = Struct of eroding bank properties initialised by
%                InitialiseVariables and set in bank erosion steps of 
%                XChannel
%      WL      = Water level in cross-section [m]
%      T       = Model time [s]
%      Flow    = Flow rate [m3/s]
%      PlotSed = optional flag to plot sediment size (default = true)
%
%   Outputs:
%      XsFigure = Structure array containing handles to created figure,
%                 axes, lines, etc.
%
%   See also: UPDATEXSPLOT, XCHANNELMODEL, INITIALISEVARIABLES.

if ~exist('PlotSed', 'var')
    PlotSed = true;
end

% Create new figure in the middle of the screen with a wide aspect ratio
scrsz = get(groot, 'ScreenSize');

if PlotSed
    XsFigure.FigureH = ...
        figure('Position', [scrsz(3)/2 - 500, scrsz(4)/2-400, 1000, 800]);
else
    XsFigure.FigureH = ...
        figure('Position', [scrsz(3)/2 - 500, scrsz(4)/2-400, 1000, 600]);
end

%% Cross-section shape plot (top panel)
if PlotSed
    XsFigure.XsAxesH = subplot(4,1,(1:2));
else
    XsFigure.XsAxesH = subplot(3,1,(1:2));
end


% Plot initial bed level
stairs(Edge.N, [Cell.Zinitial; Cell.Zinitial(end)],'k--');
hold on

% Plot final bed elevation
stairs(Edge.N, [Cell.Zfinal; Cell.Zfinal(end)],'c--');

% Plot current bed level and velocity
[XsFigure.BedVelPlotH, XsFigure.BedLineH, XsFigure.VelLineH] = ...
    plotyy(Edge.N, [Cell.Z; Cell.Z(end)], Cell.N, Cell.U, ...
                    'stairs', 'plot' );

% Plot water level
WlForPlot = WL * ones(size(Cell.N));
WlForPlot(~Cell.Wet) = NaN;
XsFigure.WlLineH = stairs(Edge.N, [WlForPlot; WlForPlot(end)], ...
                          'b-', 'LineWidth', 1);

% Identified bank markers
IsBank = Edge.IsBank~=0;
BankIDPlotX = [Edge.N(IsBank)', NaN; 
               Edge.N(IsBank)', NaN;
               NaN(1,sum(IsBank)+1)];
BankIDPlotX = BankIDPlotX(:);
BankIDPlotY = [Cell.Z(IsBank(1:end-1))', NaN; 
               Cell.Z(IsBank(2:end))', NaN;
               NaN(1,sum(IsBank)+1)];
BankIDPlotY = BankIDPlotY(:);
XsFigure.BankIDLineH = plot(BankIDPlotX, BankIDPlotY, 'r-','LineWidth',2);
StencilPlotX = [Cell.N(Bank.Top)', NaN;
                Edge.N(IsBank)', NaN;
                Cell.N(Bank.Bottom)', NaN;
                NaN(1,Bank.NBanks+1)];
StencilPlotX = StencilPlotX(:);
StencilPlotY = [Cell.Z(Bank.Top)', NaN;
                (Cell.Z(IsBank(1:end-1)) + Cell.Z(IsBank(2:end)))' ./2,NaN;
                Cell.Z(Bank.Bottom)', NaN;
                NaN(1,Bank.NBanks+1)];
StencilPlotY = StencilPlotY(:);
XsFigure.BankStencilH = plot(StencilPlotX, StencilPlotY, 'ro:');

% Labels and formatting etc
%set(XsFigure.XsAxesH,'XTickLabel','')
ylabel(XsFigure.BedVelPlotH(1), 'Elevation (m)');
ylabel(XsFigure.BedVelPlotH(2), 'Velocity (m/s)');
%legend('Initial bed','Current bed','Water level',...
%       'Location','SouthOutside','Orientation','horizontal');
XScale = Edge.N([1,end]);
set(XsFigure.BedVelPlotH(1),'XLim', XScale, 'YColor', 'k', 'XGrid', 'on');
set(XsFigure.BedVelPlotH(2),'XLim', XScale, 'YLim', [0,3], ...
    'YTickMode', 'auto');% 'YLimMode','auto');
set(XsFigure.BedLineH,'Color','k','LineWidth',1)
XsFigure.XsPlotTitleH = ...
    title(sprintf('Model time = %.0f s, Flow = %.1f, WL = %.2f', ...
                  T, Flow, WL));

xlabel('Distance across section (m)');

hold off

%% Shear stress and transport plot (2nd panel)
if PlotSed
    subplot(4,1,3);
else
    subplot(3,1,3);
end

% Plot stress and transport
[XsFigure.ShearTransPlotH, XsFigure.ShearLineH, XsFigure.TransLineH] = ...
    plotyy(Cell.N, Cell.Tau_S, Cell.N, Cell.qsS_flow_kg);

% Labels and formatting
set(XsFigure.ShearTransPlotH(1), 'XLim', XScale, 'XTickLabel', '', ...
    'YLimMode', 'auto', 'YTickMode', 'auto', 'box', 'off', 'XGrid', 'on')
set(XsFigure.ShearTransPlotH(2), 'XLim', XScale, 'XTickLabel', '', ...
    'YLimMode', 'auto', 'YTickMode', 'auto')
ylabel(XsFigure.ShearTransPlotH(1), '\tau_S (N/m^2/m)');
ylabel(XsFigure.ShearTransPlotH(2), 'q_s_,_S (kg/s/m)');
if PlotSed
    xlabel('Distance across section (m)');
end

%% Grain size plot (3rd panel)
if PlotSed
    subplot(4,1,4)

    % Plot substrate initial condition
    plot(Cell.N, Cell.SubDg_m*1000,'b--');
    hold on

    % calculate armour index
    ArmourIndex = (Cell.Dg_m)./(Cell.SubDg_m);

    % Plot active layer grain size & armour index

    [XsFigure.DgPlotH, XsFigure.DgLineH, XsFigure.ArmourLineH] = ...
        plotyy(Cell.N, Cell.Dg_m*1000, Cell.N, ArmourIndex, ...
               'semilogy', 'plot');

    % Labels and formatting
    %set(XsFigure.DgAxesH(1),'XLim',XScale)
    xlabel('Distance across section (m)');
    set(XsFigure.DgPlotH(1), 'XLim', XScale, 'YLimMode', 'auto', ...
        'YTickMode', 'auto', 'box', 'off', 'XGrid', 'on')
    set(XsFigure.DgPlotH(2), 'XLim', XScale, 'YLimMode', 'auto', ...
        'YTickMode', 'auto', 'box', 'off')
    ylabel(XsFigure.DgPlotH(1), 'D_g (mm)');
    ylabel(XsFigure.DgPlotH(2), 'Armour index');
    hold off
end

end
