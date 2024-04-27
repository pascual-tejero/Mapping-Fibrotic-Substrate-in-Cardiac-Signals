function sig = gaussFiltfilt(sig, sigma)

bt = sqrt(log(2))/(2*pi*sigma);
span = round(6*sigma/2)*2;
gaussCoeffs = gaussdesign(bt, span, 1);
pad = ceil((3*span-size(sig,2))/2)+1;
if pad > 0
    sig = [repmat(sig(:,1),1,pad) sig repmat(sig(:,end),1,pad)];
else
    pad = 0;
end

% sig = filtfilt(gaussCoeffs, 1, sig')';
sig = fftfiltfilt(gaussCoeffs, sig')'; % equivalent, but faster than filtfilt

sig = sig(:,1+pad:end-pad);

end