function c = fftConv(a, b)

a = a(:)';
b = b(:)';
    
c = ifft(fft([a, zeros(1, numel(b))]) .* fft([b, zeros(1, numel(a))]));

if any(c(numel(a):end) > 1e-10)
    warning('Elements of convolved vector are not zero')
end
c = c(1:numel(a))';
   