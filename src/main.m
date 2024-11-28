%clearing the workplace
clc, clearvars, close all;

%Simlation Parameters
mod_method = 'BPSK'; %modulation method
n_fft = 64; %IFFT/FFT size
n_cpe = 16; %cyclic prefix extension
snr = 20;   %snr(dB) 
n_taps = 8; %number of channel taps
ch_est_method = 'LS'; %channel estimation methd
save_file = 0;

%calculation modulation order from modulation method
mod_methods = {'BPSK','QPSK','8PSK','16QAM','32QAM','64QAM'};
mod_order = find(ismember(mod_methods,mod_method));

%Input data to binary stream
im = imread('../data/image.jpeg');
im_bin = dec2bin(im(:))';
im_bin = im_bin(:);

%binary stream to symbols
sym_rem = mod(mod_order - mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded,mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

%symbol modulation
%BPSK
if mod_order==1
    mod_ind = 2^(mod_order-1);
    n = 0: pi/mod_ind :2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    symbol_book = (in_phase+quadrature*1i)';
end

%Modulate data according to symbol book
X = symbol_book(cons_sym_id+1);

%use IFFT to move to time domain
%pad input signal to appropiate length
fft_rem = mod(n_fft - mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded,n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);

%Add cyclic prefix
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);

%AWGN
data_pwr = mean(abs(x_s.^2)); %calculate data power
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s + noise;
snr_meas = 10*log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));

%apply fading channel
g = exp(-(0:n_taps-1));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise,g,'same');


%use fft to move to frequency domain
%removing the cyclic prefix
x_p = reshape(x_s_noise_fading,n_fft+n_cpe, length(x_s_noise_fading)/(n_fft+n_cpe));
x_p_cpr = x_p(n_cpe+1:end, :);

%move to frequency domain
X_hat_blocks = fft(x_p_cpr);

%estimate channel
if n_taps>1
    switch(ch_est_method)
        case 'none'
        case 'LS'
            G = X_hat_blocks(:,1)./X_blocks(:,1);
            X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));
    end
end

%symbol demodulation
%remove fft padding
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

%recover data
rec_syms = knnsearch([real(symbol_book) imag(symbol_book)],[real(X_hat) imag(X_hat)]) - 1;

%parse to binary stream and remove symbol padding
rec_syms_cons = dec2bin(rec_syms);
rec_im_bin = reshape(rec_syms_cons',numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);
ber = sum(abs(rec_im_bin-im_bin))/length(im_bin);

%recover image 
rec_im = reshape(rec_im_bin,8,numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

figure;
imshow(im);
title('transmitted image');

figure;
imshow(rec_im);
title('Received image');



