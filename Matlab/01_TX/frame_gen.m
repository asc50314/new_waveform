function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-- Original Data Bit Stream (16QAM)
Frame.Data_Bitstream = randi([0,1],[4,Param.SymbolNum*Param.FFTSize]);

Symbol16QAM = zeros(1,Param.FFTSize);

for symbol_count = 1:Param.SymbolNum
  %-- Constellation Mapping (16QAM)
  gray_16QAM_map = [-3-3i -1-3i 1-3i 3-3i ...
                    -3-i  -1-i  1-i  3-i ...
                    -3+3i -1+3i 1+3i 3+3i ...
                    -3+i  -1+i  1+i  3+i];
  for data_count = 1:Param.FFTSize  % random sequence mapping to 16QAM
    index = num2str([Frame.Data_Bitstream(:,data_count+(symbol_count-1)*Param.FFTSize)]);
    index = bin2dec(index.');
    Symbol16QAM(data_count) = gray_16QAM_map(index+1)/10^0.5;
  end
  %-- iFFT
  SymbolTD = ifft(Symbol16QAM);

  %-----------------------------------------------------
  % Symbol Generating
  %-----------------------------------------------------
  %-- ‰æùÊ??≥Ëº∏Ê®°Â?Ê∑ªÂ?ËºîÂä©?Ñframe ?∂Ê?Ôºå‰?Â¶ÇCP-OFDM Ë¶ÅÂ?CP
  switch Mode.Trans
    case 'OFDM'
        %-- Add CP
        SymbolTD = [SymbolTD(end-Param.CPLength+1:end) SymbolTD]; 
        Frame.Frame_TX((symbol_count-1)*(Param.CPLength+Param.FFTSize)+1:symbol_count*(Param.CPLength+Param.FFTSize)) = SymbolTD;

      otherwise
          
  end
      
      
 end
  
end
