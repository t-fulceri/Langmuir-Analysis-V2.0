function [ step ] = signalToADCStep( signal )
%SIGNAL_TO_STEP Returns the minimum step of a signal coming from an ADC

D = diff(signal);
%Find the  step
step = min(abs(D(D ~= 0)));

end

