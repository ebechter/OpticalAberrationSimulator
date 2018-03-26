function [APSF,FPgridx,FPgridy] = ConvolveOrder(wavelength,spectrum,wave_coeff,scale)

% function [trim, wavelength] = ConvolveOrder(wavelength,spectrum,wave_coeff,scale)
            
            % trim each order beyond the edge of the detector
            
            ind1 = find(wavelength <= polyval(wave_coeff,-22.48),1,'last');
            ind2 = find(wavelength >= polyval(wave_coeff, 22.48),1,'first');
            
            wavelength = wavelength(ind1:ind2)';
            spectrum = spectrum(ind1:ind2)';
            
            % sampled at the high end. 3 pixels at red, smooth function to
            % 6 pixels at blue.
            
            horSamp = linspace(6,3,length(wavelength));
            
            vertSamp = median(horSamp);
%             pix_samp = pix_samp;
            % Custom convolution
            
            % Do the first loop iteration outside the loop. Need to
            % calculate dim first
            ii = 1;
            wfe = [0,0,0,0,0,0,0,0];
%             [PSF,center] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);
            [APSF,FPgridx,FPgridy] = makeAbbPsf(wfe,wavelength(1));%,[horzSamp,vertSamp]);
            
%             figure
%             imagesc(PSF)
            figure
            imagesc(APSF)
            
   return
            
%             how to resample?
            figure
            RAPSF = imresize(APSF,size(PSF));
            imagesc(RAPSF)
            
            
            dim = size(PSF,1);
            rectangle = zeros(dim,length(wavelength)+dim-1);
            rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
            
            for ii = 2:length(wavelength)
                
                [PSF,~] = MakePSF(scale,horSamp(ii),vert_samp);
                
                %                 Conv1 = conv2(spectrum(ii),PSF,'full');
                
                rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
                %                 rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+Conv1;
                
                
            end
            
            
            trim = rectangle(:,center(1):end-(size(PSF,2)-center(1)));
            
        end