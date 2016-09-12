% From Mike Cohen :http://www.mikexcohen.com/lectures.html
% modified Jan Klee 7.4.2015 !!Remember to concatinate trials for Speed!!

%Troubleshooting
% data=rand(1,10000);
% frequencies= 1:100;
% samplingRate=2000;

function [Power] = TFA(data,frequencies,cycles,samplingRate)

%% Create placeholder for dataMatrix
Power=zeros(length(frequencies),length(data));
%% some parameters taken out of the loop for speed
srate = samplingRate; % in hz
time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
nData = length(data);
nKern = length(time);
nConv = nData + nKern - 1;
% FFT of data
dataX = fft(data,nConv);

for i=length(frequencies):-1:1;    

 %% create a Morlet wavelet
% srate = samplingRate; % in hz
% time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
frex  = frequencies(i); % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp(1i*2*pi*frex.*time);

% create Gaussian window
s = cycles / (2*pi*frex); % this is the standard deviation of the gaussian
gaus_win  = exp((-time.^2)./(2*s^2));

% now create cimplex Morlet wavelet
cmw  = sine_wave .* gaus_win;

%% define convolution parameters

% nData = length(data);
% nKern = length(cmw);
% nConv = nData + nKern - 1;

%% FFTs

% note that the "N" parameter is the length of convolution, NOT the length
% of the original signals! Super-important!

% % FFT of data
% dataX = fft(data,nConv);

% FFT of wavelet, and amplitude-normalize in the frequency domain
cmwX = fft(cmw,nConv);
cmwX = cmwX ./ max(cmwX);

% now for convolution... and inverse Fourier transform
% conv_res = dataX.*cmwX;

conv_res_timedomain = ifft(dataX.*cmwX);


%% now back to the time domain

% cut 1/2 of the length of the wavelet from the beginning and from the end
half_wav = floor( length(cmw)/2 )+1;

% % take inverse Fourier transform
% conv_res_timedomain = ifft(conv_res);

conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav);

Power(i,:)=abs(conv_res_timedomain).^2;

%% plots 
% figure(2), clf
% plot(data,'k')
% hold on
% plot(abs(conv_res_timedomain).^2,'r','linew',2)

end
