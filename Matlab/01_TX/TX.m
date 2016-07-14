%**************************************************************************
%--------------------------------------------------------------------------
%   TX pattern generating
%   2016/7/6 by Wu Chang Ting & Friedrich Lee
%--------------------------------------------------------------------------
%**************************************************************************
clear all; close all; clc;
%--------------------------------------------------------------------------
% Mode Parameters
%--------------------------------------------------------------------------
Mode.Trans    = 'WOLA'; % OFDM/WOLA/FBMC/UFMC
Mode.Mapping  = 'QPSK'; % QPSK/16QAM
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
%-----------------------------
% Parameters Setting
%-----------------------------
Param.run             = 10;
Param.sample_rate     = 1200;
Param.SymbolNum       = 20;
Param.FFTSize         = 1024;
Param.CPratio         = 0.1;
Param.ToneNum         = 600;
Param.CarrierSp       = 0.015;

Param.UpSampleDAC  = 16;
if(Param.UpSampleDAC > 1)
  Param.DACInterpoFunc   = rcosine(1, Param.UpSampleDAC, 'fir', 0.2, 4);
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
    Param.RollOffPeriod = round((Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.RollOffRatio)/2)*2;
end

%For PSD plot
Param.PlotUpSample    = 4;
Param.PlotRightBand   = 1600;
Param.PlotLeftBand    = 1600;

%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for case_mode = 1:2
    if(case_mode == 1)
        Mode.Trans = 'OFDM';
    else
        Mode.Trans = 'WOLA';
        Param.CPLength = round(Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.CPratio);
        Param.RollOffPeriod = round((Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.RollOffRatio)/2)*2;
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
        if(Param.UpSampleDAC > 1)
          Frame.Frame_TX = upsample(Frame.Frame_TX.',Param.UpSampleDAC).';
          Frame.Frame_TX = conv(Frame.Frame_TX,Param.DACInterpoFunc);
        end
        PSD_FFTSize = ceil(length(Frame.Frame_TX)/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC)*Param.FFTSize*Param.PlotUpSample*Param.UpSampleDAC;
        PSD = fftshift(fft(Frame.Frame_TX.',PSD_FFTSize)).';
        PSD = downsample(PSD,PSD_FFTSize/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC);
        PSD = PSD.*conj(PSD);
        FrameFD = [FrameFD;PSD];
    end

    %--------------------------------------------------------------------------
    % Frame PSD
    %--------------------------------------------------------------------------
   
    AvgFrame(case_mode,:) = mean(FrameFD,1);
    clear FrameFD;
    AvgFrame(case_mode,:) = AvgFrame(case_mode,:)./max(AvgFrame(case_mode,:));
    
%     set(gca,'ytick',[-40 -30 -20 -10 0]);
    if(case_mode == 1)
        plot([-(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2:Param.CarrierSp/Param.PlotUpSample:(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2],...
        10*log10(AvgFrame(1,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                         :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'k');
        hold on      
    else
        plot([-(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2:Param.CarrierSp/Param.PlotUpSample:(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2],...
        10*log10(AvgFrame(2,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                         :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'g');
        ylabel('dB');    
        xlabel('freq [MHz]');
    end
    
end

legend('CP-OFDM','WOLA');
grid on
axis([-inf inf -120 0]);
