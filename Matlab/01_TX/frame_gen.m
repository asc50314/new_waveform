function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
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
    switch Mode.OLOverhead
      case '0'
        Frame.Frame_TX = zeros(1,Param.RollOffPeriod*Param.OverSample);
      case 'ROP/2'
        Frame.Frame_TX = zeros(1,Param.RollOffPeriod/2*Param.OverSample);
      case 'ROP'
        Frame.Frame_TX = zeros(1,Param.RollOffPeriod*Param.OverSample);
    end
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
  Frame.SymbolFD(symbol_count,:) = SymbolFD;

  % Scatter plot
  % figure (1)
  % hold on
  % scatter(real(SymbolFD(symbol_count,:)),imag(SymbolFD(symbol_count,:)),'b')

  %-----------------------------
  % Time-domain Symbol Generating
  %-----------------------------
  switch Param.OverSampleType
    case 'FFT'
      SymbolTD = ifft(ifftshift([ zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2) ...
                                  SymbolFD ...
                                  zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2)]));
    case {'SRRC', 'RC'}
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
      switch Mode.OLOverhead
        case '0'
          %Add CPrefix and CPostfix to current symbol
          SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample-Param.RollOffPeriod*Param.OverSample/2+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod*Param.OverSample/2)];
          %CP weighting
          SymbolTD(1:Param.RollOffPeriod*Param.OverSample) = SymbolTD(1:Param.RollOffPeriod*Param.OverSample) .* Param.WeightFunc;
          SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) = SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) .* fliplr(Param.WeightFunc);
          %Symbol overlap adding
          Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) + SymbolTD(1:Param.RollOffPeriod*Param.OverSample);
          %Attach rest of the symbol to Frame_TX
          Frame.Frame_TX(end+1:end+(Param.CPLength+Param.FFTSize)*Param.OverSample) = SymbolTD(Param.RollOffPeriod*Param.OverSample+1:end);
        case 'ROP/2'
          %Add CPrefix and CPostfix to current symbol
          SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample-Param.RollOffPeriod*Param.OverSample/2+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod*Param.OverSample/2)];
          %CP weighting
          SymbolTD(1:Param.RollOffPeriod*Param.OverSample) = SymbolTD(1:Param.RollOffPeriod*Param.OverSample) .* Param.WeightFunc;
          SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) = SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) .* fliplr(Param.WeightFunc);
          %Symbol overlap adding
          Frame.Frame_TX(end-Param.RollOffPeriod/2*Param.OverSample+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod/2*Param.OverSample+1:end) + SymbolTD(1:Param.RollOffPeriod/2*Param.OverSample);
          %Attach rest of the symbol to Frame_TX
          Frame.Frame_TX(end+1:end+(Param.CPLength+Param.FFTSize+Param.RollOffPeriod/2)*Param.OverSample) = SymbolTD(Param.RollOffPeriod/2*Param.OverSample+1:end);
        case 'ROP'
          %Add CPrefix and CPostfix to current symbol
          SymbolTD = [SymbolTD(end-Param.CPLength*Param.OverSample-Param.RollOffPeriod*Param.OverSample+1:end) SymbolTD SymbolTD(1:Param.RollOffPeriod*Param.OverSample)];
          %CP weighting
          SymbolTD(1:Param.RollOffPeriod*Param.OverSample) = SymbolTD(1:Param.RollOffPeriod*Param.OverSample) .* Param.WeightFunc;
          SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) = SymbolTD(end-Param.RollOffPeriod*Param.OverSample+1:end) .* fliplr(Param.WeightFunc);
          %Symbol overlap adding
          Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) = Frame.Frame_TX(end-Param.RollOffPeriod*Param.OverSample+1:end) + SymbolTD(1:Param.RollOffPeriod*Param.OverSample);
          %Attach rest of the symbol to Frame_TX
          Frame.Frame_TX(end+1:end+(Param.CPLength+Param.FFTSize+Param.RollOffPeriod)*Param.OverSample) = SymbolTD(Param.RollOffPeriod*Param.OverSample+1:end);
      end

      
  end
end
