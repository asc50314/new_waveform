function [Frame]=frame_gen_SC(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-----------------------------
% Symbol aid structure parameter setting
%-----------------------------
switch Mode.Trans
  case 'SC-WOLA'
    % Weighting function(Refer "OFDM Versus FBMC" p.95 FILTERING)
    ProtoFl = ones(1,(Param.FFTSize+Param.CPLength)*Param.OverSample);
    WeightFunc = sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample))/sum(sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample)));
    WeightFunc = conv(ProtoFl,WeightFunc);
    WeightFunc = WeightFunc(1:Param.RollOffPeriod*Param.OverSample);
    clear ProtoFl; 
end

%-----------------------------
% Original Data Bit Stream Generating
%-----------------------------
Gray16QAMmap = [  -3-3i -1-3i 1-3i 3-3i ...
                  -3-i  -1-i  1-i  3-i  ...
                  -3+3i -1+3i 1+3i 3+3i ...
                  -3+i  -1+i  1+i  3+i  ];
GrayQPSKmap = [   -3-3i 3-3i ...
                  -3+3i 3+3i];
%Original Data Bit Stream (16QAM)
switch Mode.Mapping
  case 'QPSK'
    Frame.Data_Bitstream = randi([0,1],[2,Param.SymbolNum*Param.FFTSize]);
  case '16QAM'
    Frame.Data_Bitstream = randi([0,1],[4,Param.SymbolNum*Param.FFTSize]);
end
%-----------------------------
% Symbol Generating
%-----------------------------
%Frame initializing
switch Mode.Trans
  case 'SC'
    Frame.Frame_TX = [];
  case 'SC-WOLA'
    Frame.Frame_TX = zeros(1,Param.RollOffPeriod*Param.OverSample);
end

for symbol_count = 1:Param.SymbolNum
  %-----------------------------
  % Constellation Mapping
  %-----------------------------
  SymbolTD = zeros(1,Param.FFTSize);
  for data_i = 1:Param.FFTSize  % random sequence mapping to 16QAM
    GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.FFTSize+data_i)]);
    GrayIndex = bin2dec(GrayIndex .');
    switch Mode.Mapping
      case 'QPSK'
        SymbolTD(data_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
      case '16QAM'
        SymbolTD(data_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
    end
  end
  %-----------------------------
  % Time-domain Symbol Generating
  %-----------------------------
  switch Param.OverSampleType
    case 'FFT'
      SymbolFD = fftshift(fft(SymbolTD));
      SymbolTD = ifft(ifftshift([ zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2) ...
                                  SymbolFD ...
                                  zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2)]));
    case {'SRRC', 'RC'}
      if(Param.OverSample > 1)
        SymbolTD = upsample(SymbolTD.',Param.OverSample).';
        SymbolTD = cconv(SymbolTD,Param.PulseShapeFunc,Param.FFTSize*Param.OverSample);
        SymbolTD = circshift(SymbolTD,[1,-Param.OverSample*4]);
      end
  end
  
  %-----------------------------
  % Adding aid structure for specific transmission
  %-----------------------------
  switch Mode.Trans
    case 'SC'
      % Add CP to current symbol
      SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample+1:end) SymbolTD];
      % Attach current symbol to Frame_TX
      Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
    case 'SC-WOLA' %(Refer "OFDM Versus FBMC" p.95 FILTERING)
      %Add CPrefix and CPostfix to current symbol
      SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample-Param.RollOffPeriod*Param.OverSample/2+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod*Param.OverSample/2)];
      %CP weighting
      SymbolTD(1:Param.RollOffPeriod*Param.OverSample) = SymbolTD(1:Param.RollOffPeriod*Param.OverSample) .* WeightFunc;
      SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) = SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) .* fliplr(WeightFunc);
      %Symbol overlap adding
      Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) + SymbolTD(1:Param.RollOffPeriod*Param.OverSample);
      %Attach rest of the symbol to Frame_TX
      Frame.Frame_TX(end+1:end+Param.CPLength*Param.OverSample+Param.FFTSize*Param.OverSample) = SymbolTD(Param.RollOffPeriod*Param.OverSample+1:end);
  end
end
