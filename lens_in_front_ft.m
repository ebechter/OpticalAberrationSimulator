function[x2,y2,Uout] = lens_in_front_ft(Uin,wvl,dl,f,d)

N = size(Uin,1);
k = 2*pi/wvl;
fX = (-N/2:1:N/2-1)/(N*dl);
[x2,y2]=meshgrid(wvl*f*fX);
clear('fX')

Uout = 1/(1i*wvl*f)...
    .*exp(1i*k/(2*f)*(1-d/f)*(x2.^2 +y2.^2))...
    .*ft2(Uin,dl);
end
