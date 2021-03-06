clear all
clc
close

cd('..\')
load([pwd '\01_TX\output\Frame.mat']);
cd ('02_RX')

Temp.SqErr = [];
Temp.FDVar = [];
Temp.Frame_RX_FD = [];
Temp.Frame_TX_FD = [];
if Param.DCTerm == 0
  Temp.tone_i = [Param.FFTSize/2 - ceil(Param.ToneNum/2)+1 : Param.FFTSize/2 ...
    Param.FFTSize/2+2 : Param.FFTSize/2+2 + floor(Param.ToneNum/2)-1];
else
  Temp.tone_i = [Param.FFTSize/2 - ceil(Param.ToneNum/2)+2 : Param.FFTSize/2 + floor(Param.ToneNum/2)+1];
end

for run_i = 1:Param.run
  [Frame(run_i).Frame_RX BDSuccess] = boundary_detect(Mode,Param,Frame(run_i).Frame_RX);
  if(BDSuccess == 1)
    Temp.SqErr(end+1) = 0;
    Temp.FDVar(end+1) = 0;
    for symbol_i = 1:Param.SymbolNum
      [Temp.SymbolFD Frame(run_i).Frame_RX] = fft_windowing(Mode, Param, Frame(run_i).Frame_RX);
      Temp.Frame_RX_FD = [Temp.Frame_RX_FD Temp.SymbolFD(Temp.tone_i)];
      if Param.DCTerm == 0
        Temp.SqErr(end) = Temp.SqErr(end) + sum(abs(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i) - Temp.SymbolFD(Temp.tone_i)).^2);
        Temp.FDVar(end) = Temp.FDVar(end) + sum(abs(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i)).^2);
      else
        Temp.SqErr(end) = Temp.SqErr(end) + sum(abs(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i) - Temp.SymbolFD(Temp.tone_i)).^2);
        Temp.FDVar(end) = Temp.FDVar(end) + sum(abs(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i)).^2);
      end
      Temp.Frame_TX_FD = [Temp.Frame_TX_FD Frame(run_i).SymbolFD(symbol_i,Temp.tone_i)];
    end
  end
end

Temp.MSE = sum(Temp.SqErr);
Temp.SignalPower = sum(Temp.FDVar);
SINR = 10*log10(Temp.SignalPower/Temp.MSE)


scatterplot(Temp.Frame_TX_FD)
axis([-1.5 1.5 -1.5 1.5])
title('TX constellation with 0 overhead WOLA')
% scatter(real(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i)),imag(Frame(run_i).SymbolFD(symbol_i,Temp.tone_i)),'b')

scatterplot(Temp.Frame_RX_FD)
axis([-1.5 1.5 -1.5 1.5])
title('RX constellation with 0 overhead WOLA')
% hold on
% scatter(real(Temp.SymbolFD(Temp.tone_i)),imag(Temp.SymbolFD(Temp.tone_i)),'r')
