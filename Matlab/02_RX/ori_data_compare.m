function [SqErr]=ori_data_compare(Mode, Param, SymbolFD_Ori, SymbolFD_RX)

SqErr = sum(abs(SymbolFD_Ori - SymbolFD_RX).^2);