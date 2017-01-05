function [Delta_i_bank] = BankFlux(Options, Cell, Frac, dT, Bank)
%BANKFLUX   Fractional volumetric sediment flux due to bank erosion   
%BANKFLUX is called by XChannelModel to calculate fractional sediment
%flux between cells due to bank erosion. IndentifyBanks, BankStencil and 
%TriggerBanks are used to identify active banks and the cells associated 
%with their top and bottom prior to calling BANKFLUX.
%
%   [Delta_i_bank] = BANKFLUX(Options, Cell, Edge, Frac, dT, Bank)
%   
%   Inputs:
%      Options    = Bank flux calculation options as read in to 
%                   Inputs.Bank.Flux struct by ReadModelInputs
%      Cell       = Struct of cell center properties initialised by
%                   InitialiseVariables and set in earlier steps of 
%                   XChannelModel
%      Frac       = Struct of sediment fraction properties created by
%                   InitialiseVariables
%      dT         = timestep in seconds
%      Bank       = Struct of bank properties created by BankStencil
%
%   Outputs:
%      Delta_i_bank = NCells x NFracs matrix containing fractional flux
%                     rate [m3/m/s] into and out of each cell due to bank 
%                     erosion.
%
%   See also: XCHANNELMODEL, IDENTIFYBANKS, BANKSTENCIL, TRIGGERBANKS,
%   INITIALISEVARIABLES, READMODELINPUTS

Delta_i_bank = zeros(Cell.NCells,Frac.NFracs);

for BankNo = 1:Bank.NBanks
    % if bank is active then calculate flux
    if Bank.Active(BankNo)

        Top = Bank.Top(BankNo);
        Bottom = Bank.Bottom(BankNo);

        %% Calculate total flux from top to toe
        % FluxRate = rate of sediment flux due to bank erosion [m3/m/s]
        switch Options.Approach
            case 0
                % No bank erosion flux
                FluxRate = 0;
            case 1
                % Excess slope bank erosion flux
                CellSeperation = abs(Cell.N(Top) - Cell.N(Bottom));
                ExcessHeight = max((Cell.Z(Top) - Cell.Z(Bottom)) - ...
                                   Options.Repose * CellSeperation, 0);
                FluxRate = Options.SlipRatio * ExcessHeight/2 * ...
                           sum(Cell.Width([Top,Bottom]))/2 / dT;
            case 2
                % Flux based on bank toe erosion rate (THETSD approach)
                ToeErosionRate = -sum(Cell.Delta_i_flow(Bottom,:) + ...
                                      Cell.Delta_i_slope(Bottom,:), 2);
                ToeErosionRate = max(ToeErosionRate,0);
                FluxRate = Options.ThetSD * ToeErosionRate;
            case 3
                % Flux based on bank toe transport rate
                ToeTransportRate = Cell.qsTot_flow(Bottom);
                FluxRate = Options.QsBeRatio * ToeTransportRate;
            case 4
                % Flux based on toe transport rate * slope
                Slope = (Cell.Z(Top) - Cell.Z(Bottom)) / ...
                        abs(Cell.N(Top) - Cell.N(Bottom));
                ToeTransportRate = Cell.qsTot_flow(Bottom);
                FluxRate = Options.BErodibility * ToeTransportRate * Slope;
        end

        %% Calculate fractional flux in and out of each cell
        if Options.StencilMix
            % Mix sediment through each cell between bank top and bottom
            if Top>Bottom
                Delta_i_bank(Bottom+1:Top,:) = ...
                    Delta_i_bank(Bottom+1:Top,:) - ...
                    FluxRate * Cell.Fi(Bottom+1:Top,:);
                Delta_i_bank(Bottom:Top-1,:) = ...
                    Delta_i_bank(Bottom:Top-1,:) + ...
                    FluxRate * Cell.Fi(Bottom+1:Top,:);
            else
                Delta_i_bank(Top:Bottom-1,:) = ...
                    Delta_i_bank(Top:Bottom-1,:) - ...
                    FluxRate * Cell.Fi(Top:Bottom-1,:);
                Delta_i_bank(Top+1:Bottom,:) = ...
                    Delta_i_bank(Top+1:Bottom,:) + ...
                    FluxRate * Cell.Fi(Top:Bottom-1,:);
            end
        else
            % Transport sediment directly from bank top to toe
            FluxRate_i = FluxRate * Cell.Fi(Top,:);
            Delta_i_bank(Top,:) = Delta_i_bank(Top,:) - FluxRate_i;
            Delta_i_bank(Bottom,:) = Delta_i_bank(Bottom,:) + FluxRate_i;
        end
    end
end

end

