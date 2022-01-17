%% BM-OCT Image Reconstruction.......
  %.....Parameter Initialization.... 
  clc;clear all;close all;
  NPixel = 2048;%
  lambda0 = 835E-9; % m
  bandwidth = 50E-9; % m
  lambdaMin = (lambda0 - bandwidth/2);
  lambdaMax = (lambda0 + bandwidth/2);
  dLambda = bandwidth / (NPixel-1); % m / Pixel
  lambda = lambdaMin:dLambda:lambdaMax; % m
  c0 = 300000000; % m/s
  omega0 =  2*pi*c0 / lambda0; % rad/s5
  omega = 2*pi*c0 ./ lambda; % r
  %....Brightness Control Parameters.... 
  Brightness=-10;
  Contrast=0.9;
  dB_Range=50;

 %% Raw File Directory
 RawdataDir = '.\BMScan_K_RawFilesCube';
 nImage = 384;
%% Read Raw Data
 for iImage = 1:nImage
    fid = fopen(fullfile(RawdataDir, sprintf('%05d.raw', iImage)), 'r');
    %fid = fopen(fullfile(RawdataDir,'BMresizedStack.raw'));
    bScanCameraRawDataStream = fread(fid, 'uint16');
    bScanCameraRawData(:,:,iImage) = reshape(bScanCameraRawDataStream, [2048,500]);
    %bScanCameraRawData(:,:,iImage) = reshape(bScanCameraRawDataStream, [256,256,3]);
     fclose(fid);
end

[~, ~, nSection]=size(bScanCameraRawData);

%% OCT Reconstruction Process
 for iSection = 1:nSection
    singleBscanCameraRawData = bScanCameraRawData(:, :, iSection);
    singleBscanCameraRawData1=double((singleBscanCameraRawData)-32768); 
    bScanCameraRawDataTranspose = singleBscanCameraRawData';
     
    %...Background Subtraction...........
    sizeData=size(bScanCameraRawDataTranspose);
    DCmatrix=ones(sizeData(1),1)*mean(bScanCameraRawDataTranspose);
    dcSubData = bScanCameraRawDataTranspose-DCmatrix;

    %...Resampling......................
    %...Load Spectrometer Dechirp File..
    preprocessing = load('mn_new.txt');
    aScanLength=1:1:sizeData(2);
    for noAscan=1:sizeData(1)
      dcSubData(noAscan,:) = interp1(aScanLength,dcSubData(noAscan,aScanLength),preprocessing(aScanLength),'spline','extrap');
    end
    
    %...Windowing the Data.............
    window = hann(sizeData(2)).*ones(sizeData(1),1)*4096;
    windowData=dcSubData.*window;
    
    %..Dispersion Compensation............
    a2 = 6000e-30;      
    a3 = 1000e-45; 
    
    %..Phase_correction for Dispersion compensation............
    PPdata=windowData;
    lambda = lambdaMin:dLambda:lambdaMax;
    omega = 2*pi*c0 ./ lambda;
    Taylor = -a2*((omega - omega0)).^2-a3*((omega - omega0)).^3;
    PhaseCorr = ones(sizeData(1),1)*Taylor;
    %DataPhaseCorr = DataHT.*exp(-1i*PhaseCorr);
    DataPhaseCorr = PPdata.*exp(-1i*PhaseCorr);
    DataPhaseCorr=DataPhaseCorr';
    %..IFFT......
    InvFT=abs(ifft(DataPhaseCorr,4096));
    ImageStack(:,:,iSection) = InvFT(50:1024,:);
 end

  %% Save the BM-OCT in Log Scale
 writeMode = 1;
 if writeMode == 1
     mkdir('BM-OCTReconstImage')
     for imageItr = 1:nImage
       Aline_max = max(max(ImageStack(:,:,imageItr)));
       OCT_log = (20 * log10(ImageStack(:,:,imageItr)./Aline_max)+ dB_Range)./dB_Range*65535;
       OCT_display = Contrast*(OCT_log) + Brightness;
       OCT_display(find(OCT_display >65535)) = 65535;
       OCT_display(find(OCT_display < 0)) = 0;
       OCT_display = uint16(OCT_display);
       OCT_Stack(:,:,imageItr) = OCT_display;
       basefilename = sprintf('NewReconstImage%03d.tif',imageItr);
       foldername='.\BM-OCTReconstImage\';
       destination = fullfile(foldername,basefilename);
       imwrite(OCT_display,destination);
     end
 end

  %% Angle Image Saving
 writeMode = 0;
 if writeMode == 1
     mkdir('BM-OCTAngleImage')
     for imageItr = 1:nImage
       Aline_max = max(max(angleStack(:,:,imageItr)));
       OCT_log = (20 * log10(angleStack(:,:,imageItr)./Aline_max)+ dB_Range)./dB_Range*65535;
       OCT_display = Contrast*(OCT_log) + Brightness;
       OCT_display(find(OCT_display >65535)) =65535;
       OCT_display(find(OCT_display < 0)) = 0;
       OCT_display = uint16(OCT_display);
       basefilename = sprintf('AngleImage%03d.tif',imageItr);
       foldername='.\BM-OCTAngleImage\';
       destination = fullfile(foldername,basefilename);
       imwrite(OCT_display,destination);
     end
 end
%%  Correlation based Axial Motion  Correction
volume_mcorr=corrAxialMotion(ImageStack);

%% Save the Axial Motion corrected Image in Log Scale
writeMode = 0;
if writeMode == 1
    mkdir('AxialMotionImage');
    for iImg=1:nImage
      motionCorrBscan=(volume_mcorr(:,:,iImg));
      maxIntensity=max(max(motionCorrBscan(:)));
      regImgLog=(20*log10(motionCorrBscan./maxIntensity)+dB_Range)./dB_Range*65535;
      regImgdB=Contrast*(regImgLog)+Brightness;
      regImgdB(find(regImgdB >65535)) =65535;
      regImgdB(find(regImgdB < 0)) = 0;
      regImgdB = uint16(regImgdB);
      basefilename = sprintf('motionCrtImage%03d.tif',iImg);
      foldername='.\AxialMotionImage\';
      destination = fullfile(foldername,basefilename);
      imwrite(regImgdB,destination);
    end
end

%% FFT  based Registration 
RegImage=RegFFT(volume_mcorr);

%% Save the Motion corrected Registered Image in Log Scale
writeMode = 0;
if writeMode == 1
    mkdir('RegisteredImage')
    for iImg=1:nImage
      motionCorrBscan=(RegImage(:,:,iImg));
      maxIntensity=max(max(motionCorrBscan(:)));
      regImgLog=(20*log10(motionCorrBscan./maxIntensity)+dB_Range)./dB_Range*65535;
      regImgdB=Contrast*(regImgLog)+Brightness;
      regImgdB(find(regImgdB >65535)) =65535;
      regImgdB(find(regImgdB < 0)) = 0;
      regImgdB = uint16(regImgdB); 
      basefilename = sprintf('motionCrtImage%03d.tif',iImg);
      foldername='.\RegisteredImage\';
      destination = fullfile(foldername,basefilename);
      imwrite(regImgdB,destination);
    end
end


%% Generate Speckle Variance OCT-A Images from Registered BM-OCT images
SVImage=Speckle_Variance(RegImage);

%% Save Registered Average Image
mkdir('avgImage1');
 ctr=1;
 for imgStack =1:3:384
 threeStack=RegImage(:,:,imgStack:imgStack+2);
 avgImg(:,:,ctr)=mean(threeStack,3);
 
 Aline_max = max(max(avgImg(:,:,ctr)));
 Avg_OCT_log = (20 * log10(avgImg(:,:,ctr)./Aline_max)+ dB_Range)./dB_Range*65535;
 Avg_OCT_display = Contrast*(Avg_OCT_log) + Brightness;
 Avg_OCT_display(find(Avg_OCT_display >65535)) =65535;
 Avg_OCT_display(find(Avg_OCT_display < 0)) = 0;
 
 Avg_OCT_display = uint16(Avg_OCT_display);
 basefilename = sprintf('AvgImg%03d.tif',ctr);
 foldername='.\avgImage1\';
 destination = fullfile(foldername,basefilename);
 imwrite(Avg_OCT_display,destination);
 ctr =ctr+1;
 end
%% Save OCT-A Image in Log Scale
Brightness = -10;
Contrast = 0.9;
dB_Range = 50;
writeMode = 0;
if writeMode == 0
    mkdir('OCTAImage1');
    for ii=1:size(SVImage,3)
      maxSVIntensity=max(max(SVImage(:,:,ii)));
      SVImgLog=(20*log10(SVImage(:,:,ii)./maxSVIntensity)+dB_Range)./dB_Range.*65535;
      %SVImgLog = SVImage(:,:,ii) ./ maxSVIntensity;
      SVImgdB=Contrast.*SVImgLog+Brightness;
      %SVImgdB = SVImgLog;
      SVImgdB(find(SVImgdB >65535)) = 65535;
      SVImgdB(find(SVImgdB < 0)) = 0;
      SVImgdB = uint16(SVImgdB);
      basefilename = sprintf('%03d.tif',ii);
      foldername='.\OCTAImage1\';
      destination = fullfile(foldername,basefilename);
      imwrite(SVImgdB,destination);
    end
end
%% Phase_Variance 
Dest_folder = '.\Result';
numBMscans = 3;
numFrames = 384;
strtFrame = 1;
lastFrame = ((numFrames/numBMscans - ceil((strtFrame-1)/numBMscans))*numBMscans)+(strtFrame-1);
cplxVol = RegImage(:,:,strtFrame:lastFrame);
clear volume_mcorr;
%cd(Dest_folder); 
delete *.mat;
%%  Main OCTA Process
%..Average OCT %..............
avgCplxVol = zeros([size(cplxVol,1) size(cplxVol,2) ((size(cplxVol,3)/numBMscans))]);
for I = 1:numBMscans:size(cplxVol,3)
    K = ((I-1)/numBMscans)+1;
    for J = 1:numBMscans
        cplxConjX     = cplxVol(:,:,I+(J-1)).*conj(cplxVol(:,:,I));
        bulkOffset = angle(sum(cplxConjX,1));
        avgCplxVol(:,:,K) = avgCplxVol(:,:,K) + (cplxVol(:,:,I+(J-1))...
            .*exp(-1j*repmat(bulkOffset, [size(cplxConjX,1) 1])));
    end
end

 avgOctVol_dB = 20.*log10(abs(avgCplxVol./(numBMscans)));


 %% OCT-A process %%
for I = 1:numBMscans:size(cplxVol,3)
    K = ((I-1)/numBMscans)+1;                       
    
    %..Speckle Variance & Complex Variance.........
    if numBMscans == 3
         Xconj_1 = cplxVol(:,:,I+1).*conj(cplxVol(:,:,I));
         Xconj_2 = cplxVol(:,:,I+2).*conj(cplxVol(:,:,I));            
     
        BulkOff_1 = repmat(angle(sum(Xconj_1,1)), [size(Xconj_1,1) 1]);
        BulkOff_2 = repmat(angle(sum(Xconj_2,1)), [size(Xconj_2,1) 1]);
        
        sv = var(cat(3,abs(cplxVol(:,:,I)),abs(cplxVol(:,:,I+1)),abs(cplxVol(:,:,I+2))),0,3);            
        cv = var(cat(3,(cplxVol(:,:,I)),(cplxVol(:,:,I+1).*exp(-1j*BulkOff_1)),(cplxVol(:,:,I+2).*exp(-1j*BulkOff_2))),0,3); 

    % Differential Averaging & Variance %
        Xconj_A   = cplxVol(:,:,I+1).*conj(cplxVol(:,:,I));
        Xconj_B   = cplxVol(:,:,I+2).*conj(cplxVol(:,:,I+1));
        Xconj_C   = cplxVol(:,:,I).*conj(cplxVol(:,:,I+2));
        BulkOff_A  = repmat(angle(sum(Xconj_A,1)), [size(Xconj_A,1) 1]);
        BulkOff_B  = repmat(angle(sum(Xconj_B,1)), [size(Xconj_B,1) 1]);
        BulkOff_C  = repmat(angle(sum(Xconj_C,1)), [size(Xconj_C,1) 1]);
        
        Dcplx_A = (cplxVol(:,:,I+1).*exp(-1j*BulkOff_A) - cplxVol(:,:,I));
        Dcplx_B = (cplxVol(:,:,I+2).*exp(-1j*BulkOff_B) - cplxVol(:,:,I+1));
        Dcplx_C = (cplxVol(:,:,I).*exp(-1j*BulkOff_C) - cplxVol(:,:,I+2));


    else
        Xconj_1 = cplxVol(:,:,I+1).*conj(cplxVol(:,:,I));
        Xconj_2 = cplxVol(:,:,I+2).*conj(cplxVol(:,:,I));
        Xconj_3 = cplxVol(:,:,I+3).*conj(cplxVol(:,:,I));            

        BulkOff_1 = repmat(angle(sum(Xconj_1,1)), [size(Xconj_1,1) 1]);
        BulkOff_2 = repmat(angle(sum(Xconj_2,1)), [size(Xconj_2,1) 1]);
        BulkOff_3 = repmat(angle(sum(Xconj_3,1)), [size(Xconj_3,1) 1]);

        sv = var(cat(3,abs(cplxVol(:,:,I)),abs(cplxVol(:,:,I+1)),abs(cplxVol(:,:,I+2)),abs(cplxVol(:,:,I+3))),0,3);            
        cv = var(cat(3,(cplxVol(:,:,I)),(cplxVol(:,:,I+1).*exp(-1j*BulkOff_1)),(cplxVol(:,:,I+2).*exp(-1j*BulkOff_2)),(cplxVol(:,:,I+3).*exp(-1j*BulkOff_3))),0,3); 


    end

     SV(:,:,K) = imadjust((sv-min(sv(:)))./(max(sv(:))-min(sv(:))));
     CV(:,:,K) = imadjust((cv-min(cv(:)))./(max(cv(:))-min(cv(:))));
     DA(:,:,K) = mean(cat(3,abs(Dcplx_A),abs(Dcplx_B),abs(Dcplx_C)),3);
     DV(:,:,K) = var(cat(3,abs(Dcplx_A),abs(Dcplx_B),abs(Dcplx_C)),0,3);   
end
avgSV = mov2Davg(SV, [3, 3]);
avgCV = mov2Davg(CV, [3, 3]);
avgDA = mov2Davg(DA, [3, 3]);
avgDV = mov2Davg(DV, [3, 3]);
%% OCT Intensity Image reconstruction
OCT_folder=strcat(Dest_folder,'\Intensity\');
mkdir(OCT_folder);
for ii=1:size(avgOctVol_dB,3)
    temp=uint16(avgOctVol_dB(:,:,ii));
    imwrite(temp,[OCT_folder,num2str(ii,'%03d'),'.tif']);
end


%% Speckle varaince OCTA
SV_folder=strcat(Dest_folder,'\SV\');
mkdir(SV_folder);
for ii=1:size(avgDA,3)
    temp=uint16(65536*avgSV(:,:,ii)/max(max(max(avgSV))));
    imwrite((temp),[SV_folder,num2str(ii,'%03d'),'.tif']);
end


%% Complex variance OCTA 
CV_folder=strcat(Dest_folder,'\CV\');
mkdir(CV_folder);
for ii=1:size(avgCV,3)
    temp=uint16(65535*avgCV(:,:,ii)/max(max(max(avgCV))));
    imwrite((temp),[CV_folder,num2str(ii,'%03d'),'.tif']);
end


%% Differential Averaging OCTA 
DA_folder=strcat(Dest_folder,'\DA\');
mkdir(DA_folder);
for ii=1:size(avgDA,3)
    temp=uint16(65535*avgDA(:,:,ii)/max(max(max(avgDA))));
    imwrite((temp),[DA_folder,num2str(ii,'%03d'),'.tif']);
end

%% Complex variance OCTA 
DV_folder=strcat(Dest_folder,'\DV\');
mkdir(DV_folder);
for ii=1:size(avgDV,3)
    temp=uint16(65535*avgDV(:,:,ii)/max(max(max(avgDV))));
    imwrite((temp),[DV_folder,num2str(ii,'%03d'),'.tif']);
end
%% Segmentation
inputDirectory = '.\avgImage1';
mkdir('.\avgImage1\Segmentation');
OutputDirectory = '.\avgImage1\Segmentation';
EightLayerSegmentation(inputDirectory, OutputDirectory,128);

%% Layer Mapping 
noPosition = 128;
n = input("Enter a number of Method");
switch n
    case 1
        OCTAImagedir ='.\Result\SV';
    case 2 
        OCTAImagedir ='.\Result\CV';
    case 3
        OCTAImagedir = '.\Result\DA';
    case 4
        OCTAImagedir = '.\Result\DV';
end
%OCTAImagedir ='C:\Users\AMD-PC-09\Desktop\Filter_OCTA';
SegmentedImagedir = 'E:\Matlab_OCTA\Result\Intensity\Segmentation';
%SegmentedImagedir = 'E:\PythonFiles\OCT_Angiography\Output\Layers';
[wholeInnerRetinaLayer,ChoridLayer,Layer1,Layer2,Layer3,Layer4] = LayerSegmentation(OCTAImagedir,SegmentedImagedir,noPosition);
%% Enface Projection
OCTAEnface = enfaceImage(Layer4);
%% FFT Filter
artifactRemovedEnface = FFTFilter(OCTAEnface);


%% Normalization
normEnface = artifactRemovedEnface ./ max(artifactRemovedEnface(:));
equalizedEnface = imsharpen(normEnface);
%% Denoising - Hessian Filter
%..preprocess..........
Ip = single(artifactRemovedEnface);
thr = prctile(Ip(Ip(:)>0),1) * 0.9;
Ip(Ip<=thr) = thr;
Ip = Ip - min(Ip(:));
Ip = Ip ./ max(Ip(:));
%% Compute enhancement for two different tau values
DenoisedImage = vesselness2D(Ip, 0.5:0.5:2.5, [1.5;1.5],1,true);

%figure(1),imshow(mean(RegImage(:,:,1:3),3),[]),title('Registered OCT Image');
figure(2),imshow(OCTAEnface),title('Enface OCT Image');
figure(3) ,imshow(equalizedEnface,[]),title('Enface OCT Angio Image');
figure(4),imshow(DenoisedImage),title('Denoised OCT Enface OCT Angio Image');

