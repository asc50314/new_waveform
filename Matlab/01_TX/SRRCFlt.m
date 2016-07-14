function SRRC = SRRCFlt(UpSample,RollOffFactor,Delay,Plot)
Plot(2) = 0;

t = [-Delay:1/UpSample:Delay];

%%  SRRC
% t = 0
SRRC1 = (1-RollOffFactor+4*RollOffFactor/pi);
% t = 1/4/RollOffFactor
SRRC2 = RollOffFactor/sqrt(2)*((1+2/pi)*sin(pi/4/RollOffFactor)+(1-2/pi)*cos(pi/4/RollOffFactor));
SRRC  = (sin(pi.*t.*(1-RollOffFactor)) + 4.*RollOffFactor.*t.*cos(pi.*t.*(1+RollOffFactor))) ./ ...
  (pi.*t.*(1-(4.*RollOffFactor.*t).^2));

for t_i = 1:length(t)
  if(t(t_i) == 0)
    SRRC(t_i) = SRRC1;
  end
  if(t(t_i) == 1/4/RollOffFactor)
    SRRC(t_i) = SRRC2;
  end
  if(t(t_i) == -1/4/RollOffFactor)
    SRRC(t_i) = SRRC2;
  end
end

if(Plot(1) == 1)
  H_SRRC = fft(SRRC);
  figure(2)
  subplot(2,1,1); stem(SRRC);
  subplot(2,1,2); stem(abs(H_SRRC));
end