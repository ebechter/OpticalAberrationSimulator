% interpolate per order per aberration. 

nOrders = 10;
nAberrations = 5;

% abberArray = (36,9,3,10)

aberrArray=randn(36,9,3,10);


for trace = 1:3
        for abberation = 1:10

            aberrSamples = squeeze(aberrArray(:,:,trace,aberration));
            Vq = interp2(X,Y,aberrSamples,Xq,Yq)
           
        end
end
    
    
    aberrSamples = abberArray(:,nOrders);
    wavelengthSamples = wavelengthSampleArray(:,ii);
    wavelengthVector = wavelengthArray(:,ii);
    
    for jj = 1:nAberrations
        a{ii}(:,jj) = interp1(wavelengthSamples,aberrSamples,wavelengthVector);
    end
end