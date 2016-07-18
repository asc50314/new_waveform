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
Param.run             = 100;
Param.sample_rate     = 1200;
Param.SymbolNum       = 60;
Param.FFTSize         = 1024;
Param.CPratio         = 0.1;
Param.ToneNum         = 12;
Param.CarrierSp       = 0.015;

%For WOLA
Param.RollOffRatio    = 0.0781;

%For UFMC
Param.RBsize          = 12;
Param.TXfltTap        = 102;
Param.TXfltSideAttenu = 60; %(dB)

%Symbol oversample
Param.OverSample      = 1;
if(Param.OverSample > 1)
  Param.PulseShapeFunc    = SRRCFlt(Param.OverSample, 0.2, 4);
end

%DAC up-sample
Param.UpSampleDAC  = 1;
if(Param.UpSampleDAC > 1)
  Param.DACInterpoFunc   = SRRCFlt(Param.UpSampleDAC*Param.OverSample, 0.4, 6);
end
Param.DCTerm          = 1;  % 0: DC = 0 ; 1: DC != 0
Param.ClipThreshold   = inf; % [dB], inf for no clipping

%For PSD plot
Param.PlotUpSample    = 4;
Param.PlotRightBand   = 50;
Param.PlotLeftBand    = 50;
Param.AxisModel       = 'SC'; % CF(analog freq)/DF(discrete freq)/SC(subcarrier)
Param.SpectrumMask    = 0; % 0:no spectrum mask ; 1: with spectrum mask

%--------------------------------------------------------------------------
% Auto Generated Parameters 
%--------------------------------------------------------------------------
switch Mode.Trans
  case 'OFDM'
    % Param.CPLength = round(Param.FFTSize/(1-Param.CPratio)*Param.CPratio);
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
  case 'WOLA'
    % Param.CPLength = round(Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.CPratio);
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
    % Param.RollOffPeriod = round((Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.RollOffRatio)/2)*2;
    Param.RollOffPeriod = round((Param.FFTSize*Param.RollOffRatio)/2)*2;
end

%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for case_mode = 1:2
  if(case_mode == 1)
    Param.ToneNum         = 12;
  else
    Param.ToneNum         = 12;
  end
  for clip_mode = 1:3
    if(case_mode == 1)
      Mode.Trans = 'OFDM';
      Param.ClipThreshold   = inf;
    else
      Mode.Trans = 'WOLA';
      % Param.CPLength = round(Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.CPratio);
      % Param.RollOffPeriod = round((Param.FFTSize/(1-Param.CPratio-Param.RollOffRatio)*Param.RollOffRatio)/2)*2;
      
      switch clip_mode
        case 1
          Param.ClipThreshold   = 6;
        case 2
          Param.ClipThreshold   = 8;
        case 3
          Param.ClipThreshold   = inf;
      end
    end

    FrameFD = [];
    for run_count = 1:Param.run
      Frame = frame_gen(Mode,Param);

      if(Param.OverSample > 1)
        Frame.Frame_TX = upsample(Frame.Frame_TX.',Param.OverSample).';
        Frame.Frame_TX = conv(Frame.Frame_TX,Param.PulseShapeFunc);
      end

      if(Param.UpSampleDAC > 1)
        Frame.Frame_TX = upsample(Frame.Frame_TX.',Param.UpSampleDAC).';
        Frame.Frame_TX = conv(Frame.Frame_TX,Param.DACInterpoFunc);
      end
      
      %----------------------   power checking      ---------------------------
      %---------------------   to do the clipping simulation ------------------
      if clip_mode ~= 3
        FramePower = Frame.Frame_TX.*conj(Frame.Frame_TX);
        varFrame = mean(FramePower);
        for i = 1:length(FramePower)
          if(FramePower(i) > varFrame * (10^(Param.ClipThreshold/10)))
            Frame.Frame_TX(i) = Frame.Frame_TX(i)*sqrt((varFrame * (10^(Param.ClipThreshold/10)))/FramePower(i));
          end
        end
      end
      %-----------------------------------------------------------------------
      
      PSD_FFTSize = ceil(length(Frame.Frame_TX)/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample)*Param.FFTSize*Param.PlotUpSample*Param.UpSampleDAC*Param.OverSample;
      PSD = fftshift(fft(Frame.Frame_TX.',PSD_FFTSize)).';
      PSD = downsample(PSD,PSD_FFTSize/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample);
      PSD = PSD.*conj(PSD);
      FrameFD = [FrameFD;PSD];
    end
    %--------------------------------------------------------------------------
    % Frame PSD
    %--------------------------------------------------------------------------
    AvgFrame(case_mode,:) = mean(FrameFD,1);
    clear FrameFD;
    AvgFrame(case_mode,:) = AvgFrame(case_mode,:)./max(AvgFrame(case_mode,:));
    
    switch Param.AxisModel
      case 'SC'
        plot_axis = [-(Param.PlotLeftBand+Param.PlotRightBand)/2:1/Param.PlotUpSample:(Param.PlotLeftBand+Param.PlotRightBand)/2];
      case 'CF'
        plot_axis = [-(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2:Param.CarrierSp/Param.PlotUpSample:...
                   (Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2];
      case 'DF'
    end
    if(case_mode == 1)
      plot(plot_axis,10*log10(AvgFrame(1,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                                   :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'k');
      hold on      
      break
    elseif(case_mode == 2 && clip_mode == 1)
      plot(plot_axis,10*log10(AvgFrame(2,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                                     :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'r');
    elseif(case_mode == 2 && clip_mode == 2)
      plot(plot_axis,10*log10(AvgFrame(2,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                                     :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'b');
    elseif(case_mode == 2 && clip_mode == 3)
      plot(plot_axis,10*log10(AvgFrame(2,(length(AvgFrame)/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
                                     :(length(AvgFrame)/2+Param.PlotRightBand*Param.PlotUpSample+1))),'g');                             
    end    
  end
end

axis([-inf inf -100 0]);
grid on
switch Param.AxisModel
  case 'SC'
      xlabel('Normalozed freq [1/T]');    
  case 'CF'
      xlabel('freq [MHz]');
  case 'DF'
end
ylabel('dB');  
% legend('CP-OFDM','WOLA');
legend('CP-OFDM: no clipping','WOLA: clip at 6 dB','WOLA: clip at 8 dB','WOLA: no clipping');

%--------------------------------------------------------------------------
% Spectrum Mask
%--------------------------------------------------------------------------
if(Param.SpectrumMask == 1)
  hold on
  mask_x = [-20:0.1:20]; % MHz
  mask = zeros(1,length(mask_x));
  mask(1:50) = -99;   % 10~15 MHz
  mask(51:90) = -74;  % 6~10 MHz
  mask(91:100) = -61; % 5~6 MHz
  mask(101:122) = -48; % 2.8~5 MHz
  mask(123:125) = -38; % 2.5~2.8 MHz
  mask(126:140) = -28; % 1~2.5 MHz
  mask(141:150) = -18; % 0~1 MHz
  mask(151:200) = 0;   % in band
  mask = [mask(1:200) 0 fliplr(mask(1:200))];
  plot(mask_x,mask)
end