%...Average the Depth part (Enface Image)
function motionFreeAngio=enfaceImage(layer)
EnfaceOCTangio=mean(layer(:,:,:),1);
maxAvgVal=max(EnfaceOCTangio(:));
mulVal=65535/maxAvgVal;
EnfaceOCTangioImg=mulVal.*EnfaceOCTangio;
EnfaceOCTangioImg=squeeze(EnfaceOCTangioImg);
EnfaceOCTangioImg=uint16(EnfaceOCTangioImg);
EnfaceOCTangioImg=(EnfaceOCTangioImg');
motionFreeAngio=imresize(EnfaceOCTangioImg,[300 300]);
figure(1),imshow(motionFreeAngio,[]);
end