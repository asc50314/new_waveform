%**************************************************************************
%--------------------------------------------------------------------------
%   TX pattern generating
%   2016/7/6 by Wu Chang Ting & Friedrich Lee
%--------------------------------------------------------------------------
%**************************************************************************
clear all; clc; close all;
%--------------------------------------------------------------------------
% Mode Parameters
%--------------------------------------------------------------------------
Mode.Trans    = 'WOLA'; % ['OFDM' 'WOLA' 'FBMC' 'UFMC']
Mode.Mapping  = 'QPSK'; % ['QPSK' '16QAM']
% For WOLA
Mode.OLOverhead = '0'; % ['0' 'ROP/2' 'ROP']
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
%-----------------------------
% Parameters Setting
%-----------------------------
Param.run             = 60;
Param.sample_rate     = 1200;
Param.SymbolNum       = 1;
Param.FFTSize         = 1024;
Param.CPratio         = 0.1;
Param.ToneNum         = 24;
Param.CarrierSp       = 0.015; % MHz

%For WOLA
Param.RollOffRatio    = 0.0781;

%For UFMC
Param.RBsize          = 12;
Param.RBnum           = 2;
Param.TXfltTap        = 102;
Param.TXfltSideAttenu = 60; %(dB)

%Symbol oversample
Param.OverSampleType  = 'FFT';  % ['FFT' 'SRRC' 'RC']
Param.OverSample      = 1;
if(Param.OverSample > 1)
  switch Param.OverSampleType
    case 'SRRC'
      Param.PulseShapeFunc = SRRCFlt(Param.OverSample, 0.1, 4);
    case 'RC'
      Param.PulseShapeFunc = rcosine(1,Param.OverSample,'fir', 0.1, 4);
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
Param.PlotUpSample    = 4;
Param.PlotRightBand   = 450;
Param.PlotLeftBand    = 450;
Param.PlotOOB         = [5 50];
Param.AxisModel       = 'SC'; % CF(analog freq)/DF(discrete freq)/SC(subcarrier)
Param.SpectrumMask    = 0; % 0:no spectrum mask ; 1: with spectrum mask
%-----------------------------
% Auto Generated Parameters 
%-----------------------------
switch Mode.Trans
  case 'OFDM'
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
  case 'WOLA'
    Param.CPLength = round(Param.FFTSize*Param.CPratio);
    Param.RollOffPeriod = round((Param.FFTSize*Param.RollOffRatio)/2)*2;
    % WOLA weighting function
    ProtoFl = ones(1,(Param.FFTSize+Param.CPLength)*Param.OverSample);
    Param.WeightFunc = sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample))/sum(sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample)));
    Param.WeightFunc = conv(ProtoFl,Param.WeightFunc);
    Param.WeightFunc = Param.WeightFunc(1:Param.RollOffPeriod*Param.OverSample);
    clear ProtoFl; 
  case 'UFMC'
    Param.ToneNum = Param.RBsize*Param.RBnum;
end
%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
Param.ClipThreshold = [8 6 4];
ROP = [4:4:80];

for rop_i = 1:length(ROP)
  % Param.RollOffRatio = ROP(rop_i);
  % Param.CPLength = round(Param.FFTSize*Param.CPratio);
  Param.RollOffPeriod = ROP(rop_i);
  % WOLA weighting function
  ProtoFl = ones(1,(Param.FFTSize+Param.CPLength)*Param.OverSample);
  Param.WeightFunc = sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample))/sum(sin(pi*(0:Param.RollOffPeriod*Param.OverSample-1)/(Param.RollOffPeriod*Param.OverSample)));
  Param.WeightFunc = conv(ProtoFl,Param.WeightFunc);
  Param.WeightFunc = Param.WeightFunc(1:Param.RollOffPeriod*Param.OverSample);
  clear ProtoFl; 

  for clip_i = 1:length(Param.ClipThreshold)
    FrameFD = [];
    rng(1);
    for run_count = 1:Param.run
      switch Mode.Trans
        case 'UFMC'
          Frame(rop_i,clip_i,run_count) = frame_gen_UFMC(Mode,Param);
        otherwise
          Frame(run_count) = frame_gen(Mode,Param);
      end

      if(Param.UpSampleDAC > 1)
        Frame(run_count).Frame_TX = upsample(Frame(run_count).Frame_TX.',Param.UpSampleDAC).';
        Frame(run_count).Frame_TX = conv(Frame(run_count).Frame_TX,Param.DACInterpoFunc);
      end
      %-----------------------------
      % Clipping
      %-----------------------------
      if Param.ClipThreshold(clip_i) ~= inf
        FramePower = Frame(run_count).Frame_TX.*conj(Frame(run_count).Frame_TX);
        varFrame = mean(FramePower);
        for i = 1:length(FramePower)
          if(FramePower(i) > varFrame * (10^(Param.ClipThreshold(clip_i)/10)))
            Frame(run_count).Frame_TX(i) = Frame(run_count).Frame_TX(i)*sqrt((varFrame * (10^(Param.ClipThreshold(clip_i)/10)))/FramePower(i));
          end
        end
      end
      %-----------------------------
      % PSD FFT to Get FrameFD
      %-----------------------------
      PSD_FFTSize = ceil(length(Frame(run_count).Frame_TX)/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample)*Param.FFTSize*Param.PlotUpSample*Param.UpSampleDAC*Param.OverSample;
      FrameFDTemp = fftshift(fft(Frame(run_count).Frame_TX.',PSD_FFTSize)).';
      FrameFDTemp = downsample(FrameFDTemp,PSD_FFTSize/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample);
      FrameFDTemp = FrameFDTemp.*conj(FrameFDTemp);
      FrameFD = [FrameFD;FrameFDTemp];
    end
    %-----------------------------
    % Averaging FrameFD to Get PSD
    %-----------------------------
    PSD(rop_i,:,clip_i) = mean(FrameFD,1);
    clear FrameFD;
    PSD(rop_i,:,clip_i) = PSD(rop_i,:,clip_i)./max(PSD(rop_i,:,clip_i));  
  end
end

switch Param.AxisModel
  case 'SC'
    plot_axis = [-(Param.PlotLeftBand+Param.PlotRightBand)/2 : ...
      1/Param.PlotUpSample: ...
      (Param.PlotLeftBand+Param.PlotRightBand)/2];
  case 'CF'
    plot_axis = [-(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2:Param.CarrierSp/Param.PlotUpSample:...
               (Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2];
end

for clip_i = 1:length(Param.ClipThreshold)
  figure
  % hold on
  % plot(plot_axis,10*log10(PSD(1,(length(PSD(1,:,1))/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
  %   :(length(PSD(1,:,1))/2+Param.PlotRightBand*Param.PlotUpSample+1),1)),'g');


  % plot(plot_axis,10*log10(PSDClipOOB(clip_i-1,:)),'k');
  % figure(4)
  % hold on
  % axis([-inf inf -70 -20]);
  % grid on
  % title('OOB Leakage When Clipped at 5dB')
  % xlabel('OOB Normalozed Freq. [1/T]')
  % ylabel('OOB Additive Leakage [dB]')
  % legend('OFDM','WOLA','UFMC')
  % figure
  surf(plot_axis.',ROP.',10*log10(PSD(:, ...
    (length(PSD(1,:,1))/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
    :(length(PSD(1,:,1))/2+Param.PlotRightBand*Param.PlotUpSample+1),clip_i)))
  hold on
  contour3(plot_axis.',ROP.',10*log10(PSD(:, ...
    (length(PSD(1,:,1))/2-Param.PlotLeftBand*Param.PlotUpSample+1)...
    :(length(PSD(1,:,1))/2+Param.PlotRightBand*Param.PlotUpSample+1),clip_i)),'k')
  % axis([-inf inf -inf inf -80 0])
  grid on
  title('PSD under different ROP')
  xlabel('Normalozed freq [1/T]');
  ylabel('Roll-off-period (Samples)');  
  zlabel('PSD (dB)');
  shading interp
end



% for run_count = 1:Param.run
%   Frame(run_count).Frame_RX = Frame(run_count).Frame_TX;
% end
% Frame = rmfield(Frame,'Frame_TX');
% %--------------------------------------------------------------------------
% % Save Frame
% %--------------------------------------------------------------------------
% mkdir 'output';
% save('output/Frame(clip_i,run_count).mat','Frame','Param','Mode');