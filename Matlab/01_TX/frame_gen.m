function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-----------------------------
% Symbol aid structure parameter setting
%-----------------------------
switch Mode.Trans
  case 'WOLA'
    %-- Weighting function(Refer "OFDM Versus FBMC" p.95 FILTERING)
    ProtoFl = ones(1,Param.FFTSize+Param.CPLength);
    WeightFunc = sin(pi*(0:Param.RollOffPeriod-1)/Param.RollOffPeriod)/sum(sin(pi*(0:Param.RollOffPeriod-1)/Param.RollOffPeriod));
    WeightFunc = conv(ProtoFl,WeightFunc);
    WeightFunc = WeightFunc(1:Param.RollOffPeriod);
    clear ProtoFl;
end

%-----------------------------
% Frame Generating
%-----------------------------
Gray16QAMmap = [  -3-3i -1-3i 1-3i 3-3i ...
                  -3-i  -1-i  1-i  3-i  ...
                  -3+3i -1+3i 1+3i 3+3i ...
                  -3+i  -1+i  1+i  3+i  ];
GrayQPSKmap = [  -3-3i 3-3i ...
                  -3+3i 3+3i];
%Frame initializing
switch Mode.Trans
  case 'OFDM'
    Frame.Frame_TX = [];
  case 'WOLA'
    Frame.Frame_TX = zeros(1,Param.RollOffPeriod/2);
end
%Original Data Bit Stream (16QAM)
Frame.Data_Bitstream = randi([0,1],[2,Param.SymbolNum*Param.ToneNum]);
%-----------------------------
% Symbol Generating
%-----------------------------
for symbol_count = 1:Param.SymbolNum
  %-----------------------------
  % Constellation Mapping (16QAM)
  %-----------------------------
  Symbol16QAM = zeros(1,Param.FFTSize);
  for data_count = Param.BandStart:Param.BandStart+Param.ToneNum-1  % random sequence mapping to 16QAM
    GrayIndex = num2str([Frame.Data_Bitstream(:,data_count-Param.BandStart+1)]);
    GrayIndex = bin2dec(GrayIndex .');
    Symbol16QAM(data_count) = GrayQPSKmap(GrayIndex+1)/10^0.5;
  end
  %===============test========================
  % Symbol16QAM = zeros(1,Param.FFTSize);
  % Symbol16QAM(2) = 1;
  %===============test========================
  %-----------------------------
  % Time-domain Symbol Generating
  %-----------------------------
  %iFFT
  SymbolTD = ifft(Symbol16QAM);
  %Symbol Up-sample
  if(Param.SymbolUpSample > 1)
    for data_count = 1:Param.FFTSize
      SymbolTDup((data_count-1)*Param.SymbolUpSample+1:(data_count)*Param.SymbolUpSample) = [SymbolTD(data_count) zeros(1,Param.SymbolUpSample-1)];
    end
    %RC Interpolation
    SymbolTDupInterpo = cconv(SymbolTDup, Param.SymbolInterpoFunc, length(SymbolTDup));
    SymbolTDupInterpo = circshift(SymbolTDupInterpo,[0 -Param.SymbolUpSample*Param.SymbolInterpoDelay]);
    SymbolTD = SymbolTDupInterpo;
    clear SymbolTDupInterpo;
  end
  %-----------------------------
  % Adding aid structure for specific transmission
  %-----------------------------
  switch Mode.Trans
    case 'OFDM'
      %-- Add CP
      SymbolTD = [SymbolTD(end-Param.CPLength+1:end) SymbolTD];
      %-----------------------------
      % Real frame
      %-----------------------------
      % Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
      %-----------------------------
      % Symbol frame
      %-----------------------------
      Frame.Frame_TX(end+1,:) = SymbolTD;
    case 'WOLA' %(Refer "OFDM Versus FBMC" p.95 FILTERING)
      %Add CPrefix and CPostfix
      SymbolTD = [SymbolTD(end-Param.CPLength-Param.RollOffPeriod/2+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod/2)];
      %CP weighting
      SymbolTD(1:Param.RollOffPeriod) = SymbolTD(1:Param.RollOffPeriod) .* WeightFunc;
      SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) .* fliplr(WeightFunc);
      %-----------------------------
      % Real frame
      %-----------------------------
      % %Symbol overlap adding
      % Frame.Frame_TX(end-Param.RollOffPeriod/2+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod/2+1:end) + SymbolTD(1:Param.RollOffPeriod/2);
      % %Attach rest of the symbol to frame
      % Frame.Frame_TX(end+1:end+Param.CPLength+Param.FFTSize+Param.RollOffPeriod/2) = SymbolTD(Param.RollOffPeriod/2+1:end);
      %-----------------------------
      % Symbol frame
      %-----------------------------
      %Symbol overlap adding
      Frame.Frame_TX(end,1:Param.RollOffPeriod/2) = Frame.Frame_TX(end-Param.RollOffPeriod/2+1:end) + SymbolTD(1:Param.RollOffPeriod/2);
      %Attach rest of the symbol to frame
      Frame.Frame_TX(end,Param.RollOffPeriod/2+1:Param.RollOffPeriod/2+Param.CPLength+Param.FFTSize+Param.RollOffPeriod/2) = SymbolTD(Param.RollOffPeriod/2+1:end);
      Frame.Frame_TX(end+1,:) = 0;
  end
end

if(Mode.Trans eq 'WOLD')
  Frame.Frame_TX(end,:) = [];
end
