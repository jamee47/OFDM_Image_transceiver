% Load the grayscale image
clc, clear, close all;
img = imread('image.jpeg');  % Replace 'image.jpeg' with your grayscale image file
if size(img, 3) == 3  % Check if the image is in color (RGB), convert to grayscale if necessary
    img = rgb2gray(img);
end
figure; 
imshow(img); 
title('Original Grayscale Image');

% Convert to double precision for processing
img = double(img);

% Apply IFFT and FFT to reconstruct the image (IFFT first, then FFT)
ifft_img = ifft2(img);    % IFFT
fft_img = fft2(ifft_img); % FFT
reconstructed_img = real(ifft2(fft_img));   % IFFT of the FFT

% Show the intermediate IFFT and FFT results for the grayscale image
figure;
subplot(1,3,1);
imshow(real(ifft_img), []); title('IFFT of Grayscale Image');
subplot(1,3,2);
imshow(log(abs(fft_img) + 1), []); title('FFT of IFFT Image'); % Log scale for spectrum
subplot(1,3,3);
imshow(reconstructed_img, []); title('Reconstructed Grayscale Image');

% Clip values to [0, 255] and convert back to uint8
reconstructed_img = uint8(max(min(reconstructed_img, 255), 0));

% Visualize the final reconstructed grayscale image
figure;
imshow(reconstructed_img);
title('Reconstructed Grayscale Image');
