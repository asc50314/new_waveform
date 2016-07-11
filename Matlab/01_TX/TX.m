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
Mode.Trans = 'WOLA'; % OFDM/WOLA/FBMC/UFMC
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
Param.run             = 20;
Param.sample_rate     = 1200;
Param.SymbolNum       = 30;
Param.FFTSize         = 128;
Param.CPLength        = round(Param.FFTSize*0.1); %CP number of sample after symbol up-sample
Param.BandStart       = 64;
Param.ToneNum         = 12;

Param.SymbolUpSample  = 1;
Param.SymbolInterpoDelay  = 4;
if(Param.SymbolUpSample > 1)
  Param.SymbolInterpoFunc   = rcosine(1, Param.SymbolUpSample, 'fir', 0.25, Param.SymbolInterpoDelay);
end

Param.UpSampleDAC     = 1;
if(Param.UpSampleDAC > 1)
  Param.DAC_LPF         = rcosine(1, Param.UpSampleDAC, 'fir', 0.078, 3);
end

%For WOLA
Param.RollOffPeriod   = round(Param.FFTSize*0.0781/2)*2; %Number of sample after symbol up-sample

%For PSD plot
Param.STFTupSample    = 4;
Param.STFTsize        = Param.UpSampleDAC*Param.FFTSize*Param.STFTupSample;
Param.FFTCentral      = (Param.BandStart+Param.ToneNum/2)*Param.STFTupSample*Param.STFTupSample;
Param.FFTBand         = 100*Param.UpSampleDAC*Param.STFTupSample;
Param.FFTMove         = Param.UpSampleDAC*Param.FFTSize*Param.STFTupSample/2;
%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for run_count = 1:Param.run
    Frame(run_count) = frame_gen(Mode,Param);
    Frame(run_count).Frame_TX(end+1:end+Param.STFTsize) = zeros(1,Param.STFTsize);
end

%--------------------------------------------------------------------------
% Frame PSD
%--------------------------------------------------------------------------
FrameFD = [];
for run_count = 1:Param.run
    TempFrameFD = zeros(1,Param.STFTsize);
    for move_count = 1:Param.FFTMove:length(Frame(run_count).Frame_TX)-Param.STFTsize
        TempFrameFD = TempFrameFD + abs(fft(Frame(run_count).Frame_TX(move_count:move_count+Param.STFTsize-1)));
        % TempFrameFD = fftshift(TempFrameFD);
    end
    TempFrameFD = TempFrameFD./move_count;
    FrameFD = [FrameFD;TempFrameFD]; 
end
AvgFrame = mean(FrameFD);
plot(10*log10(AvgFrame));



