%function corrAxialMotion(fn,volume,Dest_folder)
function volume_mcorr=corrAxialMotion(volume)
%     regImage=strcat(path,'regImage');
%     mkdir(regImage);
  %% Compute maximum cross-correlation (axial only)
    volume = double(volume);
    motionA = maxxcorrAx(20*log10(abs(volume))); 
    xaxis = 1:1:size(motionA,2);

    %% Set smoothing parameter
    p = polyfit(xaxis,motionA,2);
    f = polyval(p,xaxis);

    %% Compute motion correction parameters and do motion correction       
    disp_ind = motionA - f;
    m = size(volume,1);
    n = size(volume,3);
    topZero = max(disp_ind);
    botZero = abs(min(disp_ind));
    for k=1:n
        top = round(topZero-disp_ind(k));
        top_Stack(k) = top;
        volume_mcorr(top+1:top+m,:,k) = volume(:,:,k);
        disp(top);
    end
    clear volume;
    %% Crop
    cropOff = topZero+botZero;
    volume_mcorr(1:cropOff,:,:) = [];
    volume_mcorr(end-cropOff+1:end,:,:) = [];
%     for itr =1:n
%     regImgs=volume_mcorr(:,:,itr);
%     imwrite(regImgs,fullfile(regImage,sprintf('RegImgs%03d.tif',itr)));
%     end  
% cur_work_folder=pwd;
    %% Save
%     savepath = strrep(Dest_folder,'RAW DATA','Processed');
%     if exist(savepath)
%         savepath = savepath;
%     else
%         mkdir(savepath);
%     end
%     save(fullfile(savepath,[fn,'    _mcorr']), 'volume_mcorr', '-v7.3');
end