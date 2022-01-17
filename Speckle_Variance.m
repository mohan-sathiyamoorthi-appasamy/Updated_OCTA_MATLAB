function avgSV=Speckle_Variance(volume_mcorr)
bmFiles=double(volume_mcorr);
ctr=1;
for I = 1:3:384
%....Extracting Speckle Variance..................
sv = var(cat(3,abs(bmFiles(:,:,I)),abs(bmFiles(:,:,I+1)),abs(bmFiles(:,:,I+2))),0,3); 
%....Normalize the Speckle Variance.................
normalizeSV(:,:,ctr) = imadjust((sv-min(sv(:)))./(max(sv(:))-min(sv(:))));
ctr=ctr+1;
end

 for i = 1:size(normalizeSV,3)
   avgSV(:,:,i) = imguidedfilter(normalizeSV(:,:,i),'NeighborhoodSize', [3 3]);
 end 

end