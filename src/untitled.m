clc, clearvars, close all;

img = imread("../data/image.jpeg");
%img = imresize(img, 0.2);
gray_img = rgb2gray(img);


img(250:750,250:750) = 255;


imshow(img);