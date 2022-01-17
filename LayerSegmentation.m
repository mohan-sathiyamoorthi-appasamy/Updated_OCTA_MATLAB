function [wholeInnerRetinaLayer,ChoridLayer,Layer1,Layer2,Layer3,Layer4]=LayerSegmentation(OCTAImagedir,SegmentedImagedir,noPosition)
%function [wholeInnerRetinaLayer,ChoridLayer,Layer1,Layer2,Layer3,Layer4]=LayerSegmentation(noPosition)
% function wholeInnerRetinaLayer=LayerSegmentation(OCTAImagedir,SegmentedImagedir,noPosition)
% OCTAImagedir='C:\Users\Mohan\Desktop\OCTAReconstruction\SpeckleVariance';
% SegmentedImagedir='F:\Code\MATLAB\OCT-A\V2\Segmentation';

for noBscan=1:noPosition
   noBscan
matfilePath=[SegmentedImagedir '\' num2str(sprintf('%03d',noBscan)) '.mat'];
imgfilePath=[OCTAImagedir '\' num2str(sprintf('%03d',noBscan)) '.tif'];
layerPath=load(matfilePath);
layerPath = layerPath.bScan.Layers;
IlmRegion=layerPath(1,:);
GclRegion=layerPath(2,:);
IplRegion=layerPath(3,:);
InlRegion=layerPath(4,:);
OplRegion=layerPath(5,:);
RpeRegion=layerPath(8,:);

%..Read Image...
OCTImage=imread(imgfilePath);
[szImg2,szImg1]=size(OCTImage);
for i1 =1:szImg1-1
  for j1=IlmRegion(i1):OplRegion(i1)
     
    wholeInnerRetinaLayer(j1,i1,noBscan)=OCTImage(j1,i1);
  end
end

for i1 =1:szImg1
  for j1=RpeRegion(i1):szImg2
    ChoridLayer(j1,i1,noBscan)=OCTImage(j1,i1);
  end
end

for i1=1:szImg1
    for j1=IlmRegion(i1):GclRegion(i1)
        Layer1(j1,i1,noBscan)=OCTImage(j1,i1);
    end  
end

for i1=1:szImg1
        for j1=GclRegion(i1):IplRegion(i1)
            Layer2(j1,i1,noBscan)=OCTImage(j1,i1);
        end       
end

for i1=1:szImg1
        for j1=IplRegion(i1):InlRegion(i1)
            Layer3(j1,i1,noBscan)=OCTImage(j1,i1);
        end
        
end

for i2=1:szImg1
        for j2=InlRegion(i2):OplRegion(i2)
            Layer4(j2,i2,noBscan)=OCTImage(j2,i2);
        end
        
end


end

end