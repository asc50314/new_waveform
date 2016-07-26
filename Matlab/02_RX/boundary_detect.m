function [Frame_RX BDSuccess]=boundary_detect(Mode,Param,Frame_RX)

BDSuccess = 0;

switch Mode.Trans
  case 'OFDM'
  case 'WOLa'
end
BDSuccess = 1;