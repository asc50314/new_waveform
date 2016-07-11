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
Param.run             = 1;
Param.sample_rate     = 1200;
Param.SymbolNum       = 1;
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

%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for run_count = 1:Param.run
    Frame(run_count) = frame_gen(Mode,Param);
end

% %--------------------------------------------------------------------------
% % Frame PSD
% %--------------------------------------------------------------------------
% Frame_TX = zeros(1,780);
% for run_count = 1:Param.run
%     Frame_TX = Frame_TX+Frame(run_count).Frame_TX;
% end
% Frame_TX = Frame_TX/Param.run;

% % fd_Frame = fftshift(fft(Frame_TX));
% % plot(abs(fd_Frame))
% a =  zeros(1,Param.UpSampleDAC*length(Frame_TX));
% a(1:Param.UpSampleDAC:end) = Frame_TX;

% a = conv(h, a);

% psd_Frame = fftshift(fft(a));
% % figure(2)
% % plot(abs(psd_Frame))
% psd_Frame = abs(psd_Frame).^2;
% psd_Frame = psd_Frame./max(psd_Frame);
% f_spacing = Param.sample_rate/Param.FFTSize;
% freq = -f_spacing/2:f_spacing/(length(psd_Frame)-1):f_spacing/2;
% % figure(3)

% plot(freq,10*log10(psd_Frame));    
% grid on
