classdef StarTest < Spectra
    properties
        spType
        vmag
        rmag
        imag
        epsilon
        vsini
    end
    methods
        %constructor
        function obj = StarTest(type,mag,epsilon,vsini,rv,units)
            if nargin == 0
                obj.spType = 'M0V'; %spectral type
                obj.vmag = 10; % V magntiude
                obj.epsilon = 1; % not sure about this parameter
                obj.vsini= 2.5; %rotational velocity in km/s (Rot broad uses km/s)
                obj.rv = 0; %Set RV shift
                obj.spectrumUnits = 'counts';
            else
                obj.spType = type; %spectral type
                obj.vmag = mag; % V magntiude
                obj.epsilon = epsilon; % not sure about this parameter
                obj.vsini= vsini; %rotational velocity in km/s (Rot broad uses km/s)
                obj.rv = rv; %Set rv shift
                obj.spectrumUnits = units;
                
            end
            pathprefix = pwd;
            
            global allardpath
            global starfile
            starfile =  [pathprefix '/../Simulator/RefFiles/Star/fullpecautmamajek.xlsx'];
            allardpath = [pathprefix '/../Spectral_Catalogs/FAllard/CIFIST6b_trimmed_resampled/'];
            
            %Load the calibrated spectrum. %Scales to magnitude and calculates color to temp
            obj = loadSpectrum(obj); %target.Spectrum is in w/m^2/micron, target.Wavelength in microns
            %-----Stellar Effects-----%
            obj.spectrum = Star.rotBroad(obj.vsini,obj.epsilon,obj.spectrum,obj.wavelength); %Broaden spectrum and overwrite target.spectrum property
            %Convert Spectum to counts and fills target.Counts property
%             obj = energy2Counts(obj); %target.Counts is in counts/s/m^2/micron
            obj.wavelength = Star.vacShift(obj.wavelength);
%             obj.dsWavelength = Star.dopplerShift(obj.wavelength,obj.rv); %Shift Wavelength and assign DsWavelength
%             obj.counts = Star.energy2Counts(obj.dsWavelength,obj.spectrum);
        end
        
        %methods for preparing spectrum and collecting ancillary info
        function [obj] = loadSpectrum(obj)
            % DESCRIPTION: Match the spectral type to the obj spectral type and grab all
            % relevant parameters. Scale magnitude to apparent Vmag
            global allardpath
            
            [SpT,Teff,R_sun,Mv,M_J,VRc,VIc] = Spectral_data(); %load stellar reference table
            a = strcmp(SpT,[obj.spType]);%String to search for
            b = find(a,1);%find the location of the string
            
            %Relevent Parameter List
            Teff = Teff{b};
            Mv = Mv{b};
            M_J= M_J{b};
            R = R_sun{b};
            V_Rc=VRc{b};
            V_Ic=VIc{b};
            allard_file = [num2str(round(Teff/100)) '.BT-Settl.txt'];
            %Load fallard spectrum
            c = dlmread(strcat(allardpath,allard_file));
            fd = c(:,2);
            wl = c(:,1);
            
            %flux conversion at surface of the star
%             fallard_factor = -8; %a constant on fallards website
%             fd_new = 10.^(fd+fallard_factor);%conversion posted on fallard website
            joules = fd*1E-7;% conversion from ergs to joules
            jm = joules*1E4; %/cm^2 to /m^2
            jmu= jm*1E4; %/angstroms to /microns
            mu = wl*1E-4; %angstroms to microns
            
            %flux conversion earth ignoring the atmosphere
            R_solar = 6.957E8; % meters
            R = R*R_solar;
            d_psec = 10^(([obj.vmag]-Mv+5)/5); % solved from distance modulus m-M=-5+5*log10(d)
            d = d_psec*3.086e+16;%convert parsecs to distance in meters
            mj = M_J+5*log10(d_psec)-5;
            spfluxden= jmu*((R/d)^2); %flux at earth in W/m^2/um
            wavelength=mu; % wavelength in microns
            
            %Calculate magnitudes from colors
            R_mag = [obj.vmag]-V_Rc;
            I_mag = [obj.vmag]-V_Ic;
            %I_J = (I_mag-mj);
            
            %Display some data
            disp(['Vmag=',num2str([obj.vmag]),', Rmag=',num2str(R_mag),' Imag=',num2str(I_mag),...
                ' Jmag=',num2str(mj),' Stellar type= ', [obj.spType]])
            disp('-----')
            
            %assign properties to object
            obj.rmag = R_mag;
            obj.imag = I_mag;
            obj.spectrum = spfluxden;
            obj.wavelength = wavelength; 

        end
    end
     
    methods (Static)
    
            function [result] = rotBroad(vsini,epsilon,flux,wvl)
            
            if vsini==0
                result = flux;
                return
            end
            
            dwl = wvl(2) - wvl(1);
            effWvl = mean(wvl);
            
            % # The number of bins needed to create the broadening kernel
            binnHalf = floor(((vsini / 299792.458) * effWvl / dwl)) + 1;
            gwvl = ((0:4*binnHalf)-2*binnHalf)*dwl+effWvl;
            
            % Create the broadening Kernel
            dl = gwvl - effWvl;
            g = create_broad(dl,effWvl,dwl,vsini,epsilon);
            
            
            % Find and remove zeros
            % indi = find(g>0);
            g = g(g>0);
            
            result = conv(flux,g,'same')*dwl;  % not sure what the dwl is for
        end
    end
    
end
function [SpT,Teff,R_sun,Mv,M_J,VRc,VIc] = Spectral_data()
global starfile
% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\Andrew Bechter\Downloads\fullpecautmamajek.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Import the data
[~, ~, raw] = xlsread(starfile,'Sheet1');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]);

% Allocate imported array to column variable names
SpT = cellVectors(:,1);
Teff = cellVectors(:,2);
BCv = cellVectors(:,3);
R_sun = cellVectors(:,4);
Mv = cellVectors(:,5);
M_J = cellVectors(:,6);
logL = cellVectors(:,7);
BV = cellVectors(:,8);
BtVt = cellVectors(:,9);
UB = cellVectors(:,10);
VRc = cellVectors(:,11);
VIc = cellVectors(:,12);
VKs = cellVectors(:,13);
JH = cellVectors(:,14);
HK = cellVectors(:,15);
KsW1 = cellVectors(:,16);
Msun = cellVectors(:,17);
logA = cellVectors(:,18);
geby = cellVectors(:,19);
M_Ks = cellVectors(:,20);
Mbol = cellVectors(:,21);

% Clear temporary variables
% clear vars data raw cellVectors;

%Previously used radius fractions for M0-M9 not sure of the source
% 0.62
% 0.49
% 0.44
% 0.39
% 0.26
% 0.2
% 0.15
% 0.12
% 0.11
% 0.08

end
function result  = create_broad(dl,refwvl,dwl,vsini,eps)
%
%       Calculates the broadening profile.
%
%       Parameters
%       ----------
%       dl : array
%           'Delta wavelength': The distance to the reference point in
%           wavelength space [A].
%       refwvl : array
%           The reference wavelength [A].
%       dwl : float
%           The wavelength bin size [A].
%
%       Returns
%       -------
%       Broadening profile : array
%           The broadening profile according to Gray.
%

vc = vsini/299792.458;
dlmax = vc * refwvl;

c1 = 2*(1 - eps) / (pi * dlmax * (1. - eps/3));

c2 = eps / (2 * dlmax * (1 - eps/3));

result = zeros(1,length(dl));
x = dl./dlmax;
indi = find(abs(x) < 1.0);


result(indi) = c1*sqrt(1 - x(indi).^2) + c2*(1 - x(indi).^2);


%     # Correct the normalization for numeric accuracy
%     # The integral of the function is normalized, however, especially in the case
%     # of mild broadening (compared to the wavelength resolution), the discrete
%     # broadening profile may no longer be normalized, which leads to a shift of
%     # the output spectrum, if not accounted for.

result = result ./ (sum(sum(result)) * dwl);

end

