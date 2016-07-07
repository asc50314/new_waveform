function [Frame]=frame_gen(Mode,Param)

%-----------------------------------------------------
% Code for PHY Payload Field
%-----------------------------------------------------
%-- Aid structure parameter setting



%-- Original Data Bit Stream (16QAM)
<<<<<<< HEAD
for run_count = 1:Param.run
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

    %-----------------------------------------------------
    % Symbol Generating
    %-----------------------------------------------------
    %-- iFFT
    SymbolTD = ifft(Symbol16QAM);
    %-- Symbol up-sample
    for data_count = 1:Param.FFTSize
      SymbolTDup((data_count-1)*Param.SymbolUpSample+1:(data_count)*Param.SymbolUpSample) = [SymbolTD(data_count) 0 0 0];
    end
    %-- RC interpolation
    SymbolTDup = conv(SymbolTDup, Parm.PlusShap);
    SymbolTDup = SymbolTDup(17:end-16);
    %-- Adding aid structure for specific transmission 
    switch Mode.Trans
      case 'OFDM'
        %-- Add CP
        SymbolTD = [SymbolTD(end-Param.CPLength+1:end) SymbolTD]; 
        Frame.Frame_TX((symbol_count-1)*(Param.CPLength+Param.FFTSize)+1:symbol_count*(Param.CPLength+Param.FFTSize)) = SymbolTD;
      case 'WOLA'
        %-- Add CP
        SymbolTD = [SymbolTD(end-Param.CPLength+1:end) SymbolTD]; 
        Frame.Frame_TX((symbol_count-1)*(Param.CPLength+Param.FFTSize)+1:symbol_count*(Param.CPLength+Param.FFTSize)) = SymbolTD;
      otherwise

    end     
  end
end
