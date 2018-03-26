% TestScript Abberations
clc;clear;


% need wavelength, order, scale, and psf params. 

% based on star.wavelength(1) and star.wavelength(end), use wavelength and order to find a psf.

% abberation params wavelength dependent function(lambda)


load wave_coeff.mat
addpath('/Volumes/Software/Simulator/Classes')
scale = 1;
star = StarTest;
star.wavelength = star.wavelength;
% figure
% plot(star.wavelength,star.spectrum)

[rectangle, wavelength]=ConvolveOrder(star.wavelength,star.spectrum,wave_coeff(1,:,1),scale);
            