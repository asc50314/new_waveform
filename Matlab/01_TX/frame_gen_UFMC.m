function [Frame]=frame_gen_UFMC(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-----------------------------
% UFMC TX Filter
%-----------------------------
TXflt = chebwin(Param.TXfltTap,Param.TXfltSideAttenu).';    

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
  %-----------------------------
  data_count = 1;
  RB_count = 1;
  for RB = [-ceil(Param.RBnum/2):-1]
    SymbolFD = zeros(1,Param.FFTSize);
    for tone_i = [Param.FFTSize/2+1+(RB+1)*Param.RBsize-Param.RBsize : ...
                  Param.FFTSize/2+1+(RB+1)*Param.RBsize-1]  % random sequence mapping to 16QAM
      GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
      GrayIndex = bin2dec(GrayIndex .');
      switch Mode.Mapping
        case 'QPSK'
          SymbolFD(data_count,tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
        case '16QAM'
          SymbolFD(data_count,tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
      end
      data_count = data_count + 1;
    end
    %-----------------------------
    % Time-domain Symbol Generating
    %-----------------------------
    SymbolTDtemp(RB_count,:) = ifft(ifftshift(SymbolFD(RB_count,:)));
    SymbolTD(RB_count,:) = conv(SymbolTDtemp(RB_count,:), ...
      TXflt.*exp(1i*2*pi*mean([tone_i-Param.FFTSize-1 tone_i-Param.FFTSize-1+Param.RBsize-1])*(0:Param.TXfltTap-1)/Param.FFTSize));
    RB_count = RB_count + 1;
  end
  for RB = [1:floor(Param.RBnum/2)]
    SymbolFD = zeros(1,Param.FFTSize);
    for tone_i = [Param.FFTSize/2+1+(RB+1)*Param.RBsize+1 : ...
                  Param.FFTSize/2+1+(RB+1)*Param.RBsize+Param.RBsize]  % random sequence mapping to 16QAM
      GrayIndex = num2str([Frame.Data_Bitstream(:,(symbol_count-1)*Param.ToneNum+data_count)]);
      GrayIndex = bin2dec(GrayIndex .');
      switch Mode.Mapping
        case 'QPSK'
          SymbolFD(RB,tone_i) = GrayQPSKmap(GrayIndex+1)/10^0.5;
        case '16QAM'
          SymbolFD(RB,tone_i) = Gray16QAMmap(GrayIndex+1)/10^0.5;
      end
      data_count = data_count + 1;
    end
    %-----------------------------
    % Time-domain Symbol Generating
    %-----------------------------
    SymbolTDtemp(RB,:) = ifft(ifftshift(SymbolFD(RB,:)));
    SymbolTD(RB,:) = conv(SymbolTDtemp(RB,:), ...
      TXflt.*exp(1i*2*pi*mean([tone_i-Param.FFTSize-1 tone_i-Param.FFTSize-1+Param.RBsize-1])*(0:Param.TXfltTap-1)/Param.FFTSize));
  end
  SymbolTD = sum(SymbolTD,1);
  Frame.Frame_TX(end+1:end+length(SymbolTD)) = SymbolTD;
end
