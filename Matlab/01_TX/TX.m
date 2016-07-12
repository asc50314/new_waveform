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
Mode.Trans = 'OFDM'; % OFDM/WOLA/FBMC/UFMC
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
Param.run             = 1000;
Param.sample_rate     = 1200;
Param.SymbolNum       = 60;
Param.FFTSize         = 1024;
Param.CPLength        = round(Param.FFTSize*0.1); %CP number of sample after symbol up-sample
Param.BandStart       = 508;
Param.ToneNum         = 12;

Param.SymbolUpSample  = 1;
Param.SymbolInterpoDelay  = 4;
if(Param.SymbolUpSample > 1)
  Param.SymbolInterpoFunc   = rcosine(1, Param.SymbolUpSample, 'fir', 0.25, Param.SymbolInterpoDelay);
end

Param.UpSampleDAC  = 1;
if(Param.UpSampleDAC > 1)
  Param.DACInterpoFunc   = rcosine(1, Param.UpSampleDAC, 'fir', 0.25, 3);
end
%For WOLA
Param.RollOffPeriod   = round(Param.FFTSize*0.0781/2)*2; %Number of sample after symbol up-sample

%For PSD plot
% Param.STFTupSample    = 4;
% Param.STFTsize        = Param.FFTSize*Param.STFTupSample;
% Param.FFTCentral      = (Param.BandStart+Param.ToneNum/2)*Param.STFTupSample;
% Param.FFTBand         = 100*Param.STFTupSample;
% Param.FFTMove         = Param.FFTSize*Param.STFTupSample/2;
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
        TempFrameFD = fft(Frame.Frame_TX.',4*Param.FFTSize).';
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
        plot([-50:0.25:50],10*log10(AvgFrame((Param.BandStart-46)*4:(Param.BandStart+54)*4)),'k');
        hold on      
    else
       plot([-50:0.25:50],10*log10(AvgFrame((Param.BandStart-46)*4:(Param.BandStart+54)*4)),'g');
        ylabel('dB');    
        xlabel('Normalized freq [1/T]');
    end
    
end

legend('CP-OFDM','WOLA');
grid on
