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
Mode.Trans    = 'WOLA'; % ['OFDM' 'WOLA' 'FBMC' 'UFMC']
Mode.Mapping  = 'QPSK'; % ['QPSK' '16QAM']
% For WOLA
Mode.OLOverhead = 'ROP/2'; % ['0' 'ROP/2' 'ROP']
%--------------------------------------------------------------------------
% Execution Parameters
%--------------------------------------------------------------------------
%-----------------------------
% Parameters Setting
%-----------------------------
Param.run             = 2;
Param.sample_rate     = 1200;
Param.SymbolNum       = 20;
Param.FFTSize         = 1024;
Param.CPratio         = 0.1;
Param.ToneNum         = 600;
Param.CarrierSp       = 0.015;

%For WOLA
Param.RollOffRatio    = 0.0781;

%For UFMC
Param.RBsize          = 12;
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
% [Param] = param_setting(Mode);
for run_count = 1:Param.run
  Frame(run_count) = frame_gen(Mode,Param);
  if(Param.UpSampleDAC > 1)
    Frame(run_count).Frame_TX = upsample(Frame(run_count).Frame_TX.',Param.UpSampleDAC).';
    Frame(run_count).Frame_TX = conv(Frame(run_count).Frame_TX,Param.DACInterpoFunc);
  end
  %-----------------------------
  % Clipping
  %-----------------------------
  if Param.ClipThreshold ~= inf
    FramePower = Frame(run_count).Frame_TX.*conj(Frame(run_count).Frame_TX);
    varFrame = mean(FramePower);
    for i = 1:length(FramePower)
      if(FramePower(i) > varFrame * (10^(Param.ClipThreshold/10)))
        Frame(run_count).Frame_TX(i) = Frame(run_count).Frame_TX(i)*sqrt((varFrame * (10^(Param.ClipThreshold/10)))/FramePower(i));
      end
    end
  end
end
for run_count = 1:Param.run
  Frame(run_count).Frame_RX = Frame(run_count).Frame_TX;
end
Frame = rmfield(Frame,'Frame_TX');
%--------------------------------------------------------------------------
% Save Frame
%--------------------------------------------------------------------------
mkdir 'output';
save('output/Frame.mat','Frame','Param','Mode');