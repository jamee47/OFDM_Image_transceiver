clc, clearvars, close all;

%% Simlation Parameters
n_fft = 64; %IFFT/FFT size
n_cpe = 16; %cyclic prefix extension
snr = 30;   %snr(dB) 
n_taps = 8; %number of channel taps
ch_est_method = 'LS'; %channel estimation methd
save_file = 0;

%calculate modulation order from modulation method
mod_method = 'BPSK'; %modulation method
mod_order = 1;

%% input data do binary stream
im = imread('../data/bsmraau.jpg');
im_bin = dec2bin(im(:))';
im_bin = im_bin(:);

%% binary stream to symbols
%parse binary stream into mod_order bit symbols
%pads input signal to appropriate length
sym_rem = mod( mod_order-mod( length( im_bin), mod_order), mod_order);
padding = repmat ( '0', sym_rem, 1);
im_bin_padded = [im_bin; padding];
cons_data = reshape( im_bin_padded, mod_order, length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

%% symbol modulation
%BPSK

mod_ind = 2^(mod_order-1);
n = 0:pi/mod_ind:2*pi-pi/mod_ind;
in_phase = cos(n);
quadrature = sin(n);
symbol_book = (in_phase + quadrature*1i)';


X = symbol_book(cons_sym_id+1);

%% use IFFT to move to time domain
%pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X), n_fft), n_fft);
X_padded = [X; zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft, length(X_padded)/n_fft);
x = ifft(X_blocks);

%add cyclic prefix extension and shift from parallel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);

%% add AWGN
%calculate data power
data_pwr = mean(abs(x_s.^2));

%add noise to channel
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s)) + normrnd(0, sqrt(noise_pwr/2), size(x_s))*1i;
x_s_noise = x_s + noise;

%measure SNR
snr_meas = 10*log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));

%% apply fading channel
g = exp(-(0:n_taps-1));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise, g, 'same');

%% use fft to move to frequency domain
%remove cyclic prefix extension and shift from serial to parallel
x_p = reshape(x_s_noise_fading, n_fft+n_cpe, length(x_s_noise_fading)/(n_fft+n_cpe));
x_p_cpr = x_p(n_cpe + 1:end,:);

%move to frequency domain
X_hat_blocks = fft(x_p_cpr);

%% estimate channel
if n_taps > 1
 switch(ch_est_method)
 case 'none'
 case 'LS'
  G = X_hat_blocks(:,1)./X_blocks(:,1);
  X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));
 end
end

%% symbol demodulation
%remove fft padding
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

%recover data from modulated symbols
rec_syms = knnsearch([real(symbol_book) imag(symbol_book)], [real(X_hat) imag(X_hat)]) - 1;

%parse to binary stream & remove symbol padding
rec_syms_cons = dec2bin(rec_syms);
rec_im_bin = reshape(rec_syms_cons', numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);
ber = sum(abs(rec_im_bin - im_bin))/length(im_bin);

%% recover image
rec_im = reshape(rec_im_bin,8, numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

%% generate plots
%transmit constellation
figure;
subplot(2,2,1);
plot(X, 'x', 'linewidth', 2, 'markersize', 10);
axis([-1.5 1.5 -2 2]);
xlabel('In Phase');
ylabel('Quadrature');
if n_taps > 1
 title(sprintf('Transmitted Constellation Diagram\n\\rm%s Modulation\nMultipath Channel Taps: %d', mod_method, n_taps));
else
 title(sprintf('Transmitted Constellation Diagram\n\\rm%s Modulation', mod_method));
end
grid on;

%recovered constellation
subplot(2,2,2);
plot(X_hat(1:500:end),'x','markersize',3);
axis([-1.5 1.5 -2 2]);
xlabel('In Phase');
ylabel('Quadrature');
if n_taps > 1
 title(sprintf('Received Constellation Diagram\n\\rmMeasured SNR: %.2f db\nChannel Estimation: %s', snr_meas, ch_est_method));
else
 title(sprintf('bfReceived Constellation Diagram\n\\rmMeasured SNR: %.2f dB', snr_meas));
end
grid on;

%original image

subplot(2,2,3);
imshow(im);
title('\bfTransmitted Image');

%recovered image
subplot(2,2,4);
imshow(rec_im);
<<<<<<< HEAD
title(sprintf('\\bfReceived Image\n\\rmBER: %.2g', ber));
=======
title(sprintf('\\bfReceived Image\n\\rmBER: %.2g', ber));
>>>>>>> 5758a8df132e5513d99fe31a0f4497763711c32b
