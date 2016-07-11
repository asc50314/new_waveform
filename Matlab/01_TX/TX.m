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
Param.run             = 10;
Param.sample_rate     = 1200;
Param.SymbolNum       = 10;
Param.FFTSize         = 1024;
Param.CPLength        = round(Param.FFTSize*0.1); %CP number of sample after symbol up-sample
Param.BandStart       = 256;
Param.ToneNum         = 12;

Param.SymbolUpSample  = 1;
Param.SymbolInterpoDelay  = 4;
if(Param.SymbolUpSample > 1)
  Param.SymbolInterpoFunc   = rcosine(1, Param.SymbolUpSample, 'fir', 0.25, Param.SymbolInterpoDelay);
end

Param.UpSampleDAC     = 8;
Param.DAC_LPF         = rcosine(1, Param.UpSampleDAC, 'fir', 0.078, 3);

%For WOLA
Param.RollOffPeriod   = round(Param.FFTSize*0.0781/2)*2; %Number of sample after symbol up-sample

%For PSD plot
Param.FFTCentral      = Param.BandStart+Param.ToneNum/2;
Param.FFTBand         = 100;    
Param.FFTMove         = Param.FFTBand/2;
%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for run_count = 1:Param.run
    Frame(run_count) = frame_gen(Mode,Param);
end

%--------------------------------------------------------------------------
% Frame PSD
%--------------------------------------------------------------------------
FrameFD = [];
for run_count = 1:Param.run
    for move_count = 1:Param.FFTMove:length(Frame(run_count).Frame_TX)-Param.FFTBand
        TempFrameFD = fft(Frame(run_count).Frame_TX(move_count:move_count+Param.FFTBand-1));
        TempFrameFD = fftshift(TempFrameFD);
        FrameFD = [FrameFD;TempFrameFD]; 
    end
    
end
AvgFrame = mean(FrameFD);
plot([-50:49],10*log10(abs(AvgFrame)));



