function [APSF,trim,wavelength] = ConvolveOrder(wavelength,spectrum,wave_coeff,wfe,scale)

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

            % Custom convolution
            
            % Do the first loop iteration outside the loop. 
            ii = 1;
%             [PSF,center] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);

%           determine wfe to use. scale it off of wavelength... 
            wfeList{ii} = wfe*(wavelength(ii)/0.97);
            
            
            [APSF,FPgridx,FPgridy] = makeAbbPsf(wfeList{ii},wavelength(ii),[3 3],scale);
            center = [round(size(FPgridx,2)/2),round(size(FPgridy,1)/2)];
           
            dim = size(APSF,1);
            rectangle = zeros(dim,length(wavelength)+dim-1);
            rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+APSF.*spectrum(ii);
            
            for ii = 2:length(wavelength)
                wfeList{ii} = wfe*(wavelength(ii)/0.97);

%                 [PSF,~] = MakePSF(scale,horSamp(ii),vert_samp);

                [APSF,FPgridx,FPgridy] = makeAbbPsf(wfeList{ii},wavelength(ii),[horSamp(ii) vertSamp],scale);
                            
                %                 Conv1 = conv2(spectrum(ii),PSF,'full');
                
                rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+APSF.*spectrum(ii);
                %                 rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+Conv1;
                
                
            end
            
            
            trim = rectangle(:,center(1):end-(size(APSF,2)-center(1)));
            
        end