function [Frame]=frame_gen_UFMC(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-----------------------------
% UFMC TX Filter
%-----------------------------
% TXflt = chebwin(Param.TXfltTap,Param.TXfltSideAttenu).';
TXflt = fir1(Param.TXfltTap*Param.OverSample-1,(Param.RBsize/2)/(512*Param.OverSample),chebwin(Param.TXfltTap*Param.OverSample,Param.TXfltSideAttenu));

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
Frame.Frame_TX = [];
for symbol_count = 1:Param.SymbolNum
  %-----------------------------
  % Constellation Mapping
  % RBs of Left Part
  %-----------------------------
  data_count = 1;
  RB_count = 1;
  for RB = [-ceil(Param.RBnum/2):-1]
    SymbolFD = zeros(1,Param.FFTSize);
    % Random sequence mapping to 16QAM for each tone in 1 RB
    for tone_i = [Param.FFTSize/2+1+(RB+1)*Param.RBsize-Param.RBsize : ...
                  Param.FFTSize/2+1+(RB+1)*Param.RBsize-1]
      GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
      GrayIndex = bin2dec(GrayIndex .');
      switch Mode.Mapping
        case 'QPSK'
          SymbolFD(RB_count,tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
        case '16QAM'
          SymbolFD(RB_count,tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
      end
      data_count = data_count + 1;
    end
    %-----------------------------
    % Time-domain Symbol Generating for 1 RB
    %-----------------------------
    switch Param.OverSampleType
      case 'FFT'
        % dbstop 61
        SymbolTDtemp = ifft(ifftshift([ zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2) ...
                                    SymbolFD(RB_count,:) ...
                                    zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2)]));
      case {'SRRC', 'RC'}
        SymbolTDtemp = ifft(ifftshift(SymbolFD(RB_count,:)));
        if(Param.OverSample > 1)
          SymbolTDtemp = upsample(SymbolTDtemp.',Param.OverSample).';
          SymbolTDtemp = cconv(SymbolTDtemp,Param.PulseShapeFunc,Param.FFTSize*Param.OverSample);
          SymbolTDtemp = circshift(SymbolTDtemp,[1,-Param.OverSample*4]);
        end
    end
    % SymbolTDtemp = ifft(ifftshift(SymbolFD(RB_count,:)));
    % Pass modulated TX filter
    SymbolTD(RB_count,:) = conv(SymbolTDtemp, ...
      TXflt.*exp(1i*2*pi*mean([tone_i-Param.FFTSize/2-1 tone_i-Param.FFTSize/2-1-Param.RBsize+1])*(0:Param.TXfltTap*Param.OverSample-1)/(Param.FFTSize*Param.OverSample)));
    RB_count = RB_count + 1;
  end
  %-----------------------------
  % Constellation Mapping
  % RBs of Right Part
  %-----------------------------
  for RB = [1:floor(Param.RBnum/2)]
    SymbolFD = zeros(1,Param.FFTSize);
    % Random sequence mapping to 16QAM for each tone in 1 RB
    for tone_i = [Param.FFTSize/2+1+(RB-1)*Param.RBsize+1 : ...
                  Param.FFTSize/2+1+(RB-1)*Param.RBsize+Param.RBsize]
      GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
      GrayIndex = bin2dec(GrayIndex .');
      switch Mode.Mapping
        case 'QPSK'
          SymbolFD(RB_count,tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
        case '16QAM'
          SymbolFD(RB_count,tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
      end
      data_count = data_count + 1;
    end
    %-----------------------------
    % Time-domain Symbol Generating for 1 RB
    %-----------------------------
    switch Param.OverSampleType
      case 'FFT'
        SymbolTDtemp = ifft(ifftshift([ zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2) ...
                                    SymbolFD(RB_count,:) ...
                                    zeros(1,Param.OverSample*Param.FFTSize/2-Param.FFTSize/2)]));
      case {'SRRC', 'RC'}
        SymbolTDtemp = ifft(ifftshift(SymbolFD(RB_count,:)));
        if(Param.OverSample > 1)
          SymbolTDtemp = upsample(SymbolTDtemp.',Param.OverSample).';
          SymbolTDtemp = cconv(SymbolTDtemp,Param.PulseShapeFunc,Param.FFTSize*Param.OverSample);
          SymbolTDtemp = circshift(SymbolTDtemp,[1,-Param.OverSample*4]);
        end
    end
    % SymbolTDtemp = ifft(ifftshift(SymbolFD(RB_count,:)));
    % Pass modulated TX filter
    SymbolTD(RB_count,:) = conv(SymbolTDtemp, ...
      TXflt.*exp(1i*2*pi*mean([tone_i-Param.FFTSize/2-1 tone_i-Param.FFTSize/2-1-Param.RBsize+1])*(0:Param.TXfltTap*Param.OverSample-1)/(Param.FFTSize*Param.OverSample)));
    RB_count = RB_count + 1;
  end
  % Sum all the SymbolTD from all the RBs
  SymbolTD = sum(SymbolTD,1);
  % Attach current symbol to Frame
  Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
end
