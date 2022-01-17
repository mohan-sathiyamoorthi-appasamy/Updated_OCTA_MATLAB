function volume_regis=RegFFT(volume_log)
Repeatframe=3;
NumBScan=384;
usfac = 2;
volume_regis = zeros(size(volume_log));
IdxShift = zeros(2,NumBScan);
%Shift = zeros(2,NumBScan);
for framegroup=1:NumBScan/Repeatframe
frame1 = volume_log(:,:,(framegroup-1)*Repeatframe+1);
for num = 1:Repeatframe
frameNum = (framegroup-1) * Repeatframe + num;
frame2 = volume_log(:,:,frameNum);
a = dftregistration(fft2(frame1),fft2(frame2),usfac);
IdxShift(1,frameNum) = a(3);
IdxShift(2,frameNum) = a(4);
%Shift(:,frameNum) = [fix(IdxShift(3)); fix(IdxShift(4))];
frame = imresize(volume_log(:,:,frameNum),usfac);
frame = circshift(frame,round([IdxShift(1,frameNum)*usfac IdxShift(2,frameNum)*usfac]));
volume_regis(:,:,frameNum) = imresize(frame,1/usfac);
end
% framegroup
end
% mkdir('RegisteredImage')
% for i=1:192
% regImage=uint16(65535*volume_regis(:,:,i)/max(max(max(volume_regis))));
% basefilename = sprintf('RegImg%03d.tif',i);
% foldername='.\RegisteredImage\';
% destination = fullfile(foldername,basefilename);
% imwrite(regImage,destination);
 %end
end
