function img7 = FFTFilter(img)
    img2=fftshift(fft2(img,512,512));
    %img3=log(1+abs(img2));
    %figure,imshow(img3,[]);
%% Generating rectangular mask_Method 
    BW=ones(512,512);
    BW(1:250,250:258)=0;
    BW(260:512,250:258)=0;
    img5=(BW.*img2);
    %figure,imshow(BW,[]);
%% IFFT
    img6= ifft2(ifftshift(img5));
    img7=abs(img6(1:size(img,1),1:size(img,2)));
   figure,imshow(img7,[]);
end