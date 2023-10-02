data = load('tide_data');
fftdata = fft(data(2,:));
N = length(data);

subplot(2,3,1);
plot(data(2,:));
title('Original Data');

subplot(2,3,2);
plot(fftdata/N,'*');
title(sprintf('Original, # points = %d',N));


subplot(2,3,3);
% NN=40; plot(0:NN-1,abs(fftdata(1:NN))/N);
plot(-N/2:(N/2-1),abs(fftshift(fftdata))/N);
title('ABS(Fourier Transform)');


x = linspace(0,N,1000);
cutoff = 0.2;
maxfft = max(abs(fftdata(2:end)));
filter = cutoff*maxfft;
compressed = zeros(1,length(fftdata));
filtered_points = 1;
for n = 1:length(fftdata)
    if (abs(fftdata(n))>filter)
        compressed(n) = fftdata(n);
        filtered_points = filtered_points +1;
    end
end
    
compressed = zeros(1,length(fftdata));
compressed(1:cutoff) = fftdata(1:cutoff);
compressed(end:-1:(end-cutoff+1)) = fftdata(end:-1:(end-cutoff+1));

subplot(2,3,4);
plot(real(ifft(compressed)));
title('Reconstructed Data');

subplot(2,3,5);
plot(compressed/N,'*');
title(sprintf('Filtered, # points = %d',filtered_points));

subplot(2,3,6);
% NN=40; plot(0:NN-1,abs(fftdata(1:NN))/NN);
plot(-N/2:(N/2-1),abs(fftshift(compressed))/N,'*');
title('ABS(Fourier Transform)');
height = max(abs(compressed))/N;
text(0.2*N,height,sprintf('Compression: %2.1f%%',100*(1-filtered_points/N)));


% omega = 2*pi/N;
% A0 = real(fftdata(1));
% A11 = real(fftdata(12)); B11 = imag(fftdata(12)); 
% A21 = real(fftdata(21)); B21 = imag(fftdata(21)); 
% A32 = real(fftdata(32)); B32 = imag(fftdata(32)); 
% reconstruct = (A0 + ...
%                   + A11*cos(omega*11*x) + B11*cos(omega*11*x) ...
%                   + A21*cos(omega*21*x) + B21*cos(omega*21*x) ...
%                   + A32*cos(omega*32*x) + B32*cos(omega*32*x))/N;
% plot(x,  reconstruct)

