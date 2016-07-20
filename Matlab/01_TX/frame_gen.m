function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-----------------------------
% Symbol aid structure parameter setting
%-----------------------------
switch Mode.Trans
  case 'WOLA'
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
    Frame.Data_Bitstream = randi([0,1],[2,Param.SymbolNum*Param.ToneNum]);
  case '16QAM'
    Frame.Data_Bitstream = randi([0,1],[4,Param.SymbolNum*Param.ToneNum]);
end
%-----------------------------
% Symbol Generating
%-----------------------------
%Frame initializing
switch Mode.Trans
  case 'OFDM'
    Frame.Frame_TX = [];
  case 'WOLA'
    Frame.Frame_TX = zeros(1,Param.RollOffPeriod*Param.OverSample);
end

for symbol_count = 1:Param.SymbolNum
  %-----------------------------
  % Constellation Mapping
  %-----------------------------
  SymbolFD = zeros(1,Param.FFTSize);
  data_count = 1;
  if Param.DCTerm == 0
    for tone_i = [Param.FFTSize/2 - ceil(Param.ToneNum/2)+1 : Param.FFTSize/2 ...
        Param.FFTSize/2+2 : Param.FFTSize/2+2 + floor(Param.ToneNum/2)-1]  % random sequence mapping to 16QAM
        GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
        GrayIndex = bin2dec(GrayIndex .');
        switch Mode.Mapping
          case 'QPSK'
            SymbolFD(tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
          case '16QAM'
            SymbolFD(tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
        end
        data_count = data_count + 1;
    end
  else
    for tone_i = [Param.FFTSize/2 - ceil(Param.ToneNum/2)+2 : Param.FFTSize/2 + floor(Param.ToneNum/2)+1]  % random sequence mapping to 16QAM
        GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
        GrayIndex = bin2dec(GrayIndex .');
        switch Mode.Mapping
          case 'QPSK'
            SymbolFD(tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
          case '16QAM'
            SymbolFD(tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
        end
        data_count = data_count + 1;
    end
  end
  %===============test========================
  % Symbol16QAM = zeros(1,Param.FFTSize);
  % Symbol16QAM(2) = 1;
  %===============test========================
  %-----------------------------
  % Time-domain Symbol Generating
  %-----------------------------
  switch Param.OverSampleType
    case 'FFT'
      SymbolTD = ifft(ifftshift([ zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2) ...
                                  SymbolFD ...
                                  zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2)]));
    case 'SRRC'
      SymbolTD = ifft(ifftshift(SymbolFD));
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
    case 'OFDM'
      % Add CP to current symbol
      SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample+1:end) SymbolTD];
      % Attach current symbol to Frame_TX
      Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
    case 'WOLA' %(Refer "OFDM Versus FBMC" p.95 FILTERING)
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
