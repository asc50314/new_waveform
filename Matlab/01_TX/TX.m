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
Mode.Trans    = 'SC'; % OFDM/WOLA/FBMC/UFMC/SC/SC-WOLA
Mode.Mapping  = 'QPSK'; % QPSK/16QAM
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
%-----------------------------
% Parameters Setting
%-----------------------------
Param.run             = 10;
Param.SapmleRate      = 15.36;
Param.SymbolNum       = 2;
Param.FFTSize         = 128;
Param.CPratio         = 0;
Param.ToneNum         = 600;
Param.CarrierSp       = Param.SapmleRate/Param.FFTSize;

%For WOLA
Param.RollOffRatio    = 0.01;

%For UFMC
Param.RBsize          = 12;
Param.TXfltTap        = 102;
Param.TXfltSideAttenu = 60; %(dB)

%Symbol oversample
Param.OverSampleType  = 'FFT';  % FFT/SRRC/RC
Param.OverSample      = 8;
if(Param.OverSample > 1)
  switch Param.OverSampleType
    case 'SRRC'
      Param.PulseShapeFunc = SRRCFlt(Param.OverSample/2*3, 0.1, 12);
    case 'RC'
      Param.PulseShapeFunc = rcosine(1,Param.OverSample/2*3,'fir/sqrt', 0.3, 12);
  end
end

%DAC up-sample
Param.UpSampleDAC  = 1;
if(Param.UpSampleDAC > 1)
  Param.DACInterpoFunc   = SRRCFlt(Param.UpSampleDAC*Param.OverSample, 0.3, 6);
end
Param.DCTerm          = 0;  % 0: DC = 0 ; 1: DC != 0
Param.ClipThreshold   = inf; % [dB], inf for no clipping

%For PSD plot
Param.PlotUpSample    = 8;
Param.PlotRightBand   = 511;
Param.PlotLeftBand    = 512;
Param.AxisModel       = 'CF'; % CF(analog freq)/DF(discrete freq)/SC(subcarrier)
Param.SpectrumMask    = 0; % 0:no spectrum mask ; 1: with spectrum mask

%--------------------------------------------------------------------------
% Auto Generated Parameters 
%--------------------------------------------------------------------------
switch Mode.Trans
  case {'OFDM', 'SC'}
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
  case {'WOLA', 'SC-WOLA'}
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
    Param.RollOffPeriod = round((Param.FFTSize*Param.RollOffRatio)/2)*2;
end
%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
% [Param] = param_setting(Mode);
for case_mode = 1:2
  if(case_mode == 1)
    Param.ToneNum         = 600;
  else
    Param.ToneNum         = 666;
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
    Param.RollOffPeriod = round((Param.FFTSize*Param.RollOffRatio)/2)*2;
  end
  for clip_mode = 1:3
    if(case_mode == 1)
      Mode.Trans = 'SC';
      Param.ClipThreshold   = inf;
    else
      Mode.Trans = 'SC-WOLA';  
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
      switch Mode.Trans
        case {'SC', 'SC-WOLA'}
          Frame = frame_gen_SC(Mode,Param);
        otherwise
          Frame = frame_gen(Mode,Param);
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
legend('CP-SC: no clipping','SC-WOLA: clip at 6 dB','SC-WOLA: clip at 8 dB','SC-WOLA: no clipping');

%--------------------------------------------------------------------------
% Spectrum Mask
%--------------------------------------------------------------------------
if(Param.SpectrumMask == 1)
  hold on
  mask_x = [-20 -15 -15 -11 -11 -10 -10 -7.8 -7.8 -7.5 -7.5 -6 -6 -5 -5]; % MHz
  mask = zeros(1,length(mask_x));
  mask(1:2) = -99;   % 10~15 MHz
  mask(3:4) = -74;  % 6~10 MHz
  mask(5:6) = -61; % 5~6 MHz
  mask(7:8) = -48; % 2.8~5 MHz
  mask(9:10) = -38; % 2.5~2.8 MHz
  mask(11:12) = -28; % 1~2.5 MHz
  mask(13:14) = -18; % 0~1 MHz
  mask(15) = 0;   % in band
  mask = [mask fliplr(mask)];
  mask_x = [mask_x -fliplr(mask_x)];
  plot(mask_x,mask)
end