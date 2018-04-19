% TestScript Abberations
clc;clear; 

% need wavelength, order, scale, and psf params. 

% based on star.wavelength(1) and star.wavelength(end), use wavelength and order to find a psf.

% abberation params wavelength dependent function(lambda)


load wave_coeff.mat
addpath('/Volumes/Software/Simulator/Classes')
scale = 1;
star = StarTest;
pixSamp = 3;
R = 275e3;
dlam = 1/R/pixSamp;
step = dlam/scale;
interplambda = star.wavelength(1):step:star.wavelength(end);
interpflux = interp1(star.wavelength, star.spectrum, interplambda,'linear');


% figure
% plot(star.wavelength,star.spectrum)
wfe = ones(16,1);

% [APSF,order,wavelength]=ConvolveOrder(interplambda,interpflux./max(interpflux),wave_coeff(1,:,1),wfe,scale);


% indexing starts at 0 for Zernike #s. 0 is piston... etc. follows Wyant scheme 

% characteristic wavefront error in radians describing the amplitude (not RMS or P-V)
wfe(1) = 0;     % piston
wfe(2) = 0;     % distortion/tilt
wfe(3) = 0;     % -
wfe(4) = 0;     % defocus
wfe(5) = 10;     % Primary astigmatism
wfe(6) = 0;     % -
wfe(7) = 0;    % Primary coma
wfe(8) = 0;    % -
wfe(9) = 0;     % Spherical abb
wfe(10) = 0;    % Elliptical coma
wfe(11) = 0;    % -
wfe(12) = 0;    % Secondary astigmatism
wfe(13) = 0;    % -
wfe(14) = 0;    % Secondary coma
wfe(15) = 0;    % -
wfe(16) = 0;    % Secondary spherical abb

[APSF2,order2,wavelength2] = ConvolveOrder(interplambda,interpflux./max(interpflux),wave_coeff(1,:,1),wfe,scale);


% save PSF2 APSF FPgridx FPgridy
% 
% addpath('/Volumes/Software/ImageProcessing')
% run('FitSimulatedImage.m')