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
Param.run = 1;
Param.sample_rate = 1200;
Param.SymbolNum = 60;
Param.FFTSize = 12;
Param.CPLength = 1;
Param.UpSampleDAC = 8;
h  =rcosine(1, Param.UpSampleDAC, 'fir', 0.078, 3);
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
Frame_TX = zeros(1,780);
for run_count = 1:Param.run
    Frame_TX = Frame_TX+Frame(run_count).Frame_TX;
end
Frame_TX = Frame_TX/Param.run;

% fd_Frame = fftshift(fft(Frame_TX));
% plot(abs(fd_Frame))
a =  zeros(1,Param.UpSampleDAC*length(Frame_TX));
a(1:Param.UpSampleDAC:end) = Frame_TX;

a = conv(h, a);

psd_Frame = fftshift(fft(a));
% figure(2)
% plot(abs(psd_Frame))
psd_Frame = abs(psd_Frame).^2;
psd_Frame = psd_Frame./max(psd_Frame);
f_spacing = Param.sample_rate/Param.FFTSize;
freq = -f_spacing/2:f_spacing/(length(psd_Frame)-1):f_spacing/2;
% figure(3)

plot(freq,10*log10(psd_Frame));    
grid on
