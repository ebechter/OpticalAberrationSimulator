        function [PSF,center] = MakePSF(scale,pix_samp,vert_samp)
            
            % dispersion direction
            fwhm = pix_samp*scale;
            sigmax=fwhm/(2*sqrt(2*log(2)));
            
            % cross dispersion direction
            fwhmy = vert_samp*scale;
            sigmay=fwhmy/(2*sqrt(2*log(2)));
            
            % trim according to bigger dimension
            Nsigmas = 9;
            val = Nsigmas * sigmay;
            %
            
            %             clip = 3/(2*sqrt(2*log(2)));
            %             val = 25*clip;
            
            oddGrid=2*round(val/2)-1; % not super important but lets go with an odd grid size every time
            
            [MatX,MatY]=meshgrid(1:1:oddGrid,1:1:oddGrid); % make the mesh grid
            
            center = [round(size(MatX,2)/2),round(size(MatY,1)/2)]; % center x and y in the (likely) non-integer center of grid ;
            % this is the most important line - the center of the PSF must fall in a pixel center
            % subtract the center in the x direction from the full convolution to get a "central" conv. result
            %....maybe not certain.
            
            PSF=circ_gauss(MatX,MatY,[sigmax,sigmay],center);
            
            
        end
