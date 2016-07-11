function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-- Symbol aid structure parameter setting
switch Mode.Trans
  case 'WOLA'
    %-- Weighting function
    ProtoFl = ones(1,Param.FFTSize+Param.CPLength);
    WeightFunc = sin(pi*(0:Param.RollOffPeriod-1)/Param.RollOffPeriod)/sum(sin(pi*(0:Param.RollOffPeriod-1)/Param.RollOffPeriod));
    WeightFunc = conv(ProtoFl,WeightFunc);
    WeightFunc = WeightFunc(1:Param.RollOffPeriod);
    clear ProtoFl;
end

for run_count = 1:Param.run
  %Frame initializing
  switch Mode.Trans
    case 'OFDM'
      Frame.Frame_TX = [];
    case 'WOLA'
      Frame.Frame_TX = zeros(1,Param.RollOffPeriod/2);
  end
  %-- Original Data Bit Stream (16QAM)
  Frame.Data_Bitstream = randi([0,1],[4,Param.SymbolNum*Param.ToneNum]);
  Symbol16QAM = zeros(1,Param.FFTSize);
  for symbol_count = 1:Param.SymbolNum
    %-- Constellation Mapping (16QAM)
    Gray16QAMmap = [  -3-3i -1-3i 1-3i 3-3i ...
                      -3-i  -1-i  1-i  3-i  ...
                      -3+3i -1+3i 1+3i 3+3i ...
                      -3+i  -1+i  1+i  3+i  ];
    Symbol16QAM = zeros(1,Param.FFTSize);
    for data_count = Param.BandStart:Param.BandStart+Param.ToneNum-1  % random sequence mapping to 16QAM
      GrayIndex = num2str([Frame.Data_Bitstream(:,data_count+(symbol_count-1)*Param.FFTSize)]);
      GrayIndex = bin2dec(GrayIndex .');
      Symbol16QAM(data_count) = Gray16QAMmap(GrayIndex+2)/10^0.5;
    end
    %-- Time-domain Symbol Generating
    %-- iFFT
    SymbolTD = ifft(Symbol16QAM);

    if(Param.SymbolUpSample > 1)
      %-- Symbol Up-sample
      for data_count = 1:Param.FFTSize
        SymbolTDup((data_count-1)*Param.SymbolUpSample+1:(data_count)*Param.SymbolUpSample) = [SymbolTD(data_count) zeros(1,Param.SymbolUpSample-1)];
      end
      %-- RC Interpolation
      SymbolTDupInterpo = cconv(SymbolTDup, Param.SymbolInterpoFunc, length(SymbolTDup));
      SymbolTDupInterpo = circshift(SymbolTDupInterpo,[0 -Param.SymbolUpSample*Param.SymbolInterpoDelay]);
      SymbolTD = SymbolTDupInterpo;
      clear SymbolTDupInterpo;
    end
    %-- Adding aid structure for specific transmission
    switch Mode.Trans
      case 'OFDM'
        %-- Add CP
        SymbolTD = [SymbolTD(end-Param.CPLength+1:end) SymbolTD]; 
        % Frame.Frame_TX((symbol_count-1)*(Param.CPLength+Param.FFTSize)+1:symbol_count*(Param.CPLength+Param.FFTSize)) = SymbolTD;
        Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
      case 'WOLA'
        %-- Add CP
        SymbolTD = [SymbolTD(end-Param.CPLength-Param.RollOffPeriod/2+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod/2)];
        SymbolTD(1:Param.RollOffPeriod) = SymbolTD(1:Param.RollOffPeriod) .* WeightFunc;
        SymbolTD(end-Param.RollOffPeriod+1:end) = SymbolTD(end-Param.RollOffPeriod+1:end) .* fliplr(WeightFunc);
        % Frame.Frame_TX((symbol_count-1)*(Param.CPLength+Param.FFTSize)+1:symbol_count*(Param.CPLength+Param.FFTSize)) = SymbolTD;
        %Overlap adding
        Frame.Frame_TX(end-Param.RollOffPeriod/2+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod/2+1:end) + SymbolTD(1:Param.RollOffPeriod/2);
        %Attach to frame
        Frame.Frame_TX(end-Param.CPLength-Param.FFTSize-Param.RollOffPeriod/2+1:end) = SymbolTD(Param.RollOffPeriod/2+1:end);
    end
  end
end
