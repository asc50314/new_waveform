%**************************************************************************
%--------------------------------------------------------------------------
%   TX pattern generating
%   2016/7/6 by Friedrich Lee
%--------------------------------------------------------------------------
%**************************************************************************
clear all; close all; clc;
%--------------------------------------------------------------------------
% Mode Parameters
%--------------------------------------------------------------------------
Mode.Trans    = 'OFDM'; % OFDM/WOLA/FBMC/UFMC
Mode.Mapping  = 'QPSK'; % QPSK/16QAM
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
%-----------------------------
% Parameters Setting
%-----------------------------
Param.run             = 100;
Param.sample_rate     = 1200;
Param.SymbolNum       = 1;
Param.FFTSize         = 1024;
Param.CPratio         = 0.1;
Param.ToneNum         = 600;

Param.UpSampleDAC  = 1;
if(Param.UpSampleDAC > 1)
  Param.DACInterpoFunc   = rcosine(1, Param.UpSampleDAC, 'fir', 0.25, 3);
end

%For WOLA
Param.RollOffRatio    = 0.0781;

%--------------------------------------------------------------------------
% Auto Generated Parameters 
%--------------------------------------------------------------------------
switch Mode.Trans
  case 'OFDM'
    Param.CPLength = round(Param.FFTSize/(1-Param.CPratio)*Param.CPratio);
  case 'WOLA'
    Param.CPLength = round(Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.CPratio);
    Param.RollOffPeriod = round(Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.RollOffRatio);
end

%For PSD plot
Param.PlotUpSample    = 4;
Param.PlotRightBand   = 55;
Param.PlotLeftBand    = 45;

%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for case_mode = 1:2
    if(case_mode == 1)
        Mode.Trans = 'OFDM';
    else
        Mode.Trans = 'WOLA';
    end
     FrameFD = [];
    for run_count = 1:Param.run
        Frame = frame_gen(Mode,Param);
%         Frame.Frame_TX(end+1:end+Param.STFTsize) = zeros(1,Param.STFTsize);   
%         TempFrameFD = zeros(1,Param.STFTsize);
%         for move_count = 1:Param.FFTMove:length(Frame.Frame_TX)-Param.STFTsize
%         TempFrameFD = TempFrameFD + abs(fft(Frame.Frame_TX(move_count:move_count+Param.STFTsize-1))).^2;
            % TempFrameFD = fftshift(TempFrameFD);
%         end!
        TempFrameFD = fftshift(fft(Frame.Frame_TX.',Param.PlotUpSample*Param.FFTSize)).';
        PSD = TempFrameFD.*conj(TempFrameFD);
%         TempFrameFD = TempFrameFD./move_count;
        FrameFD = [FrameFD;PSD]; 
    end

    %--------------------------------------------------------------------------
    % Frame PSD
    %--------------------------------------------------------------------------
   
    AvgFrame = mean(FrameFD);
   AvgFrame = AvgFrame/max(AvgFrame);
    
%     set(gca,'ytick',[-40 -30 -20 -10 0]);
    if(case_mode == 1)
        plot([-(Param.PlotLeftBand+Param.PlotRightBand)/2:1/Param.PlotUpSample:(Param.PlotLeftBand+Param.PlotRightBand)/2],...
        10*log10(AvgFrame((length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample)...
                         :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample))),'k');
        hold on      
    else
        plot([-(Param.PlotLeftBand+Param.PlotRightBand)/2:1/Param.PlotUpSample:(Param.PlotLeftBand+Param.PlotRightBand)/2],...
        10*log10(AvgFrame((length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample)...
                         :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample))),'g');
        ylabel('dB');    
        xlabel('Normalized freq [1/T]');
    end
    
end

legend('CP-OFDM','WOLA');
grid on
axis([-inf inf -100 0]);
