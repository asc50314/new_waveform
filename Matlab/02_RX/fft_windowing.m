function [SymbolFD Frame_RX]=fft_windowing(Mode,Param,Frame_RX)


switch Mode.Trans
  case 'OFDM'
    %-----------------------------------------------------
    % Removing CP
    %-----------------------------------------------------
    SymbolTD = Frame_RX(1:Param.CPLength+Param.FFTSize);
    SymbolTD = SymbolTD(Param.CPLength+1:end);
    Frame_RX(1:Param.CPLength+Param.FFTSize) = [];
  case 'WOLA'
    switch Mode.OLOverhead
      case '0'
        %-----------------------------------------------------
        % Removing CP
        %-----------------------------------------------------
        SymbolTD = Frame_RX(1:Param.CPLength+Param.RollOffPeriod+Param.FFTSize);
        SymbolTD = SymbolTD(Param.CPLength+1:end);
        %-----------------------------------------------------
        % RX WOLA
        %-----------------------------------------------------
        % Weighting
        SymbolTD(1:Param.RollOffPeriod) = SymbolTD(1:Param.RollOffPeriod) .* Param.WeightFunc;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) .* fliplr(Param.WeightFunc);
        WeightingA = SymbolTD(1:Param.RollOffPeriod);
        WeightingB = SymbolTD(end-Param.RollOffPeriod+1:end);
        % Overlaping and adding
        SymbolTD(1:Param.RollOffPeriod)         = SymbolTD(1:Param.RollOffPeriod) + WeightingB;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) + WeightingA;
        % Truncate to fft window
        SymbolTD = SymbolTD(1+Param.RollOffPeriod/2:end-Param.RollOffPeriod/2);
        Frame_RX(1:Param.CPLength+Param.FFTSize) = [];
      case 'ROP/2'
        %-----------------------------------------------------
        % Removing CP
        %-----------------------------------------------------
        SymbolTD = Frame_RX(1:Param.CPLength+Param.RollOffPeriod+Param.FFTSize);
        SymbolTD = SymbolTD(Param.CPLength+1:end);
        %-----------------------------------------------------
        % RX WOLA
        %-----------------------------------------------------
        % Weighting
        SymbolTD(1:Param.RollOffPeriod) = SymbolTD(1:Param.RollOffPeriod) .* Param.WeightFunc;
        % SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) .* fliplr(Param.WeightFunc);
        WeightingA = SymbolTD(1:Param.RollOffPeriod);
        WeightingB = SymbolTD(end-Param.RollOffPeriod+1:end);
        % Overlaping and adding
        SymbolTD(1:Param.RollOffPeriod)         = SymbolTD(1:Param.RollOffPeriod) + WeightingB;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) + WeightingA;
        % Truncate to fft window
        SymbolTD = SymbolTD(1+Param.RollOffPeriod/2:end-Param.RollOffPeriod/2);
        Frame_RX(1:Param.CPLength+Param.FFTSize+Param.RollOffPeriod/2) = [];
      case 'ROP'
        %-----------------------------------------------------
        % Removing CP
        %-----------------------------------------------------
        SymbolTD = Frame_RX(1:Param.CPLength+2*Param.RollOffPeriod+Param.FFTSize);
        SymbolTD = SymbolTD(Param.CPLength+Param.RollOffPeriod/2+1:end-Param.RollOffPeriod/2);
        %-----------------------------------------------------
        % RX WOLA
        %-----------------------------------------------------
        % Weighting
        SymbolTD(1:Param.RollOffPeriod) = SymbolTD(1:Param.RollOffPeriod) .* Param.WeightFunc;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) .* fliplr(Param.WeightFunc);
        WeightingA = SymbolTD(1:Param.RollOffPeriod);
        WeightingB = SymbolTD(end-Param.RollOffPeriod+1:end);
        % Overlaping and adding
        SymbolTD(1:Param.RollOffPeriod)         = SymbolTD(1:Param.RollOffPeriod) + WeightingB;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) + WeightingA;
        % Truncate to fft window
        SymbolTD = SymbolTD(1+Param.RollOffPeriod/2:end-Param.RollOffPeriod/2);
        Frame_RX(1:Param.CPLength+Param.FFTSize+Param.RollOffPeriod) = [];
    end
end

SymbolFD = fftshift(fft(SymbolTD));