%**************************************************************************
%--------------------------------------------------------------------------
%   TX pattern generating
%   2016/7/6 by Wu Chang Ting & Friedrich Lee
%--------------------------------------------------------------------------
%**************************************************************************
clear all; clc;
%--------------------------------------------------------------------------
% Mode Parameters
%--------------------------------------------------------------------------
Mode.Trans    = 'WOLA'; % ['OFDM' 'WOLA' 'FBMC' 'UFMC']
Mode.Mapping  = 'QPSK'; % ['QPSK' '16QAM']
% For WOLA
Mode.OLOverhead = 'ROP'; % ['0' 'ROP/2' 'ROP']
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
Param.ToneNum         = 12;
Param.CarrierSp       = 0.015; % MHz

%For WOLA
Param.RollOffRatio    = 0.0781;

%For UFMC
Param.RBsize          = 12;
Param.TXfltTap        = 102;
Param.TXfltSideAttenu = 60; %(dB)

%Symbol oversample
Param.OverSampleType  = 'FFT';  % ['FFT' 'SRRC' 'RC']
Param.OverSample      = 4;
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
Param.DCTerm          = 1;  % 0: DC = 0 ; 1: DC != 0
Param.ClipThreshold   = inf; % [dB], inf for no clipping

%For PSD plot
Param.PlotUpSample    = 4;
Param.PlotRightBand   = 50;
Param.PlotLeftBand    = 50;
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
end
%--------------------------------------------------------------------------
% Frame Generating
%--------------------------------------------------------------------------
clip_i = 1;
Param.ClipThreshold = [inf 4:0.1:8.5];
for clip_i = 1:length(Param.ClipThreshold)
  FrameFD = [];
  for run_count = 1:Param.run
    Frame(clip_i,run_count) = frame_gen(Mode,Param);
    if(Param.UpSampleDAC > 1)
      Frame(clip_i,run_count).Frame_TX = upsample(Frame(clip_i,run_count).Frame_TX.',Param.UpSampleDAC).';
      Frame(clip_i,run_count).Frame_TX = conv(Frame(clip_i,run_count).Frame_TX,Param.DACInterpoFunc);
    end
    %-----------------------------
    % Clipping
    %-----------------------------
    if Param.ClipThreshold(clip_i) ~= inf
      FramePower = Frame(clip_i,run_count).Frame_TX.*conj(Frame(clip_i,run_count).Frame_TX);
      varFrame = mean(FramePower);
      for i = 1:length(FramePower)
        if(FramePower(i) > varFrame * (10^(Param.ClipThreshold(clip_i)/10)))
          Frame(clip_i,run_count).Frame_TX(i) = Frame(clip_i,run_count).Frame_TX(i)*sqrt((varFrame * (10^(Param.ClipThreshold(clip_i)/10)))/FramePower(i));
        end
      end
    end
    %-----------------------------
    % PSD FFT to Get FrameFD
    %-----------------------------
    PSD_FFTSize = ceil(length(Frame(clip_i,run_count).Frame_TX)/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample)*Param.FFTSize*Param.PlotUpSample*Param.UpSampleDAC*Param.OverSample;
    FrameFDTemp = fftshift(fft(Frame(clip_i,run_count).Frame_TX.',PSD_FFTSize)).';
    FrameFDTemp = downsample(FrameFDTemp,PSD_FFTSize/Param.FFTSize/Param.PlotUpSample/Param.UpSampleDAC/Param.OverSample);
    FrameFDTemp = FrameFDTemp.*conj(FrameFDTemp);
    FrameFD = [FrameFD;FrameFDTemp];
  end
  %-----------------------------
  % Averaging FrameFD to Get PSD
  %-----------------------------
  PSD(clip_i,:) = mean(FrameFD,1);
  clear FrameFD;
  PSD(clip_i,:) = PSD(clip_i,:)./max(PSD(clip_i,:));  
end

switch Param.AxisModel
  case 'SC'
    plot_axis = [Param.ToneNum/2+Param.PlotOOB(1) : ...
      1/Param.PlotUpSample : ...
      Param.ToneNum/2+Param.PlotOOB(2)];
  case 'CF'
    plot_axis = [-(Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2:Param.CarrierSp/Param.PlotUpSample:...
               (Param.PlotLeftBand+Param.PlotRightBand)*Param.CarrierSp/2];
end

PSD = PSD(:,length(PSD(clip_i,:))/2+(Param.ToneNum/2+Param.PlotOOB(1))*Param.PlotUpSample+1 : ...
  length(PSD(clip_i,:))/2+(Param.ToneNum/2+Param.PlotOOB(2))*Param.PlotUpSample+1);

for clip_i = 2:length(Param.ClipThreshold)
  %-----------------------------
  % Find Clipping OOB Leakage
  %-----------------------------
  PSDClipOOB(clip_i-1,:) = PSD(clip_i,:) - PSD(1,:);
end
PSDClipOOB(find(PSDClipOOB<0)) = 0;
% for clip_i = 2:length(Param.ClipThreshold)
%   figure
%   plot(plot_axis,10*log10(PSDClipOOB(clip_i-1,:)),'k');
%   axis([-inf inf -100 0]);
%   grid on
% end
figure(1)
surf(plot_axis,Param.ClipThreshold(2:end),10.*log10(PSDClipOOB))
axis([-inf inf -inf inf -80 -20])
grid on
xlabel('Normalozed freq [1/T]');
ylabel('Clipping Threshold (dB)');  
zlabel('Out-of-band Leakage (dB)');


% for run_count = 1:Param.run
%   Frame(run_count).Frame_RX = Frame(run_count).Frame_TX;
% end
% Frame = rmfield(Frame,'Frame_TX');
% %--------------------------------------------------------------------------
% % Save Frame
% %--------------------------------------------------------------------------
% mkdir 'output';
% save('output/Frame(clip_i,run_count).mat','Frame','Param','Mode');