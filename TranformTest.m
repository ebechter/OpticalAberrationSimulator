%Description: Calculates the pupil based on the telescope
%parameters with a phase term. Also calculates the pupil based
%on the focal length of the fiber lens pair by backwards
%propagating. The fiber position error (from not centering the
%stage, and ADC offsets are incorporated in fiber pupil plane.
clear all
close all

Wavelength = 1;% microns
wfe = [0,0,0,0,0,0,0,0];

%Zernike phase map
W = 0; %
for z = 1:length(wfe)
    Z = ZernikeCalc(z,1); % use ZernikeCalc to produce normalized zernike surface
    Z = padarray(Z,[400,400]);% increase padding on output;
    Z = Z./(max(max(Z)));
    Ab = wfe(z)*Z;
    W = W+Ab;
end

%Constants (we should make this accept variables, or use
%defaults)
N = size(W,1);% variable 7/25/17 Sampling points
L = 50e-3;% length of grid in meters
dl = L/N; %Pupil plane grid spacing (meters)
D = 5.8e-3; % Pupil diameter

k = 2*pi./(1e-6*Wavelength); %wavenumber
alpha = 0.11; % secondary blocking fraction
f = 24.5e-3;% (20.8e-3) focal length of lens
d = f; %distance before the lens
gridparams = [dl,f,d];

%Pupil definition
[x1,y1] = meshgrid((-N/2+0.5:N/2-0.5)*dl); %grid in pupil plane
PPgrid(:,:,1) = x1;
PPgrid(:,:,2) = y1;


PupilOffset = 0; %1mm offset
Pupil_cenx = x1 + PupilOffset;
Pupil_ceny = y1 + PupilOffset;
Uin1 = circ(Pupil_cenx,Pupil_ceny,D);%outer circle
% Uin2 = circ(Pupil_cenx,Pupil_ceny,D*alpha);%inner circle
% pupil = Uin1-Uin2;%total pupil plance image

[gauss]=circ_gauss(Pupil_cenx,Pupil_ceny,[D/8,D/8],[0,0]);
gauss = gauss./(max(max(gauss)));
pupil = Uin1.*gauss;

% W =pupil.*W;
Pupil = pupil.*exp(1i*W(2:end,2:end)); %Complex Pupil plane with phase term

%Propogate the pupil plane to focal plane (electric fields) at
%all wavelengths
[x2,y2,PSF] = lens_in_front_ft(Pupil,1e-6*Wavelength,dl,f,d);
FPgridx = x2;
FPgridy = y2;