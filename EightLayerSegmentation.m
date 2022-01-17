% function layers = ThreeLayerSegExe(image, params)
function EightLayerSegmentation(inputDirectory, OutputDirectory,numBM3Image)    
    
% 
%inputDirectory = 'D:\OCT_Datas\DataSet_20191111\Reconstructed_data\LE_MaculaCube_C_VidyaFileRename\Reconstructed_image\Reconstruc_data';

%inputDirectory = 'D:\OCT_Datas\DataSet_20191111\20191111115751_LE_MaculaCube_A_Vidya\TifFilesCube';
%OutputDirectory = 'D:\OCT_Datas\DataSet_20191111\Reconstructed_data\LE_MaculaCube_C_VidyaFileRename\Reconstructed_image\Reconstruc_data\out';
%thicknessDirectory= ' D:\OCT_Datas\DataSet_20191111\TifFilesCube';
listOfFiles = dir([inputDirectory '/*.tif']);
nImage = length(listOfFiles);

doPlotRegionBoundaries = 1;
doShowImageSegmentation = 1;
doWriteOverlayImage = 1;

for iImageSlice = 1:numBM3Image/3
 imFileName = fullfile(inputDirectory, listOfFiles(iImageSlice).name);
  inputImage = imread (imFileName);%%oct engine
  % image=imresize(image,[500,250]);
  % image_1=imcrop(image_1,[1 30 950 500]);
   %image=imcrop(image,[1 11 512 500]);
%    
 %  imshow(image_1)
    % hold on
    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
%     if nargin < 2
%         params = [];

%     end
    
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    %
    % Validate the image is a 2D, grayscale image
    %
%     if isempty(image_1)
%         error('Image cannot be empty');
%     end
    
%     nDimensions = size(image,3);
%     
%     if nDimensions > 3
%         image(:,:,4:nDimensions) = [];
%     end
%     
%     if nDimensions == 3
%         image = rgb2gray(image);
%     end
%     
    %
    % Validate the parameters
    %
%     if isempty(params)
%         params = normal_getParameters();
%     end
    
    image_1 = double(inputImage);
    originalSize = size(image_1);
    
    
    %----------------------------------------------------------------------
    %  Make the registered black boundaries invalid
    %----------------------------------------------------------------------
    
    [image_1, invalidImage] = removeImageBorders(image_1, [0,255]);
    originalInvalidImage = invalidImage;
    
    
   
   
   
    MIN_WEIGHT=0.00001;

    
%    parameters = SegmentImageParameters;
    parameters.ALGORITHM_TYPE = 'normal';
    parameters.DEFAULT_LATERAL_RESOLUTION = 2;  % um/pixel
    parameters.DEFAULT_AXIAL_RESOLUTION = 2; % um/pixel
    parameters.X_RESOLUTION = 14.0;  % um/pixel
    parameters.Y_RESOLUTION = 5.8; % um/pixel
    
    % GetBwImageParameters
    parameters.getBwImageParams.X_FILTER_SIZE = 147.4;  % um
    parameters.getBwImageParams.Y_FILTER_SIZE = 73.7;  % um
    parameters.getBwImageParams.SIGMA = 11;  % um
    parameters.getBwImageParams.X_STREL_SIZE = 40.2;  % um;
    parameters.getBwImageParams.Y_STREL_SIZE = 20.1;  % um;
    parameters.getBwImageParams.MIN_CLUSTER_SIZE = 18000;  % um^2
    parameters.getBwImageParams.NUM_BORDERS = 4;
    
    % GetGraphCutParameters
   parameters.graphCutParams.NUM_LAYERS = 8;
   parameters.graphCutParams.MAX_NUM_LAYERS = 8;

%     parameters.graphCutParams.NUM_LAYERS = 3;
%     parameters.graphCutParams.MAX_NUM_LAYERS = 3;
    parameters.graphCutParams.SMOOTHING_CORRECTION = [0.05,0.05,0.1,0.1,0.05,0.05,0.05,0.05];
%     parameters.graphCutParams.LAYER_INDICES = [1,8,8,7];
%     parameters.graphCutParams.MATRIX_INDICES = [1,2,3,5];
    
   parameters.graphCutParams.LAYER_INDICES = [1,8,8,7,6,4,3,5,4,2];
   parameters.graphCutParams.MATRIX_INDICES = [1,2,3,5,5,5,4,4,5,6];
    
    % WeightingMatrixParameters
    parameters.graphCutParams.weightingMatrixParams.X_FILTER_SIZE = 40.2;  % um
    parameters.graphCutParams.weightingMatrixParams.Y_FILTER_SIZE = 20.1;  % um
    parameters.graphCutParams.weightingMatrixParams.SIGMA = 6;
    parameters.graphCutParams.weightingMatrixParams.EDGE_FILTER = [0.2 0.6 0.2; -0.2 -0.6 -0.2];
    parameters.graphCutParams.weightingMatrixParams.WEIGHT_RANGES =[];
    parameters.graphCutParams.weightingMatrixParams.WEIGHT_RANGES = { ...
        {[0,1]}, ...                        % 01: dark-light
        {[0,0.5],[0.5,1]}, ...              % 02: bright, distance
        {[0,0.4],[0.4,1]}, ...              % 03: light-dark, distance
        {[0,0.6],[0.6,1]}, ...              % 04: light-dark, distance
        {[0,0.5],[0.5,1]}, ...              % 05: dark-light distance
        {[0,0.5],[0.5,0.7],[0.7,1.1]}, ...  % 06: light-dark, dark, distance
    };

    % otherParams
    parameters.otherParams.EYE = [];
    parameters.otherParams.SCAN_ORIENTATION = [];
    
    
    %----------------------------------------------------------------------
    %  Resize the image to reduce computational complexity
    %----------------------------------------------------------------------
     LATERAL_RESOLUTION = 5.0;
     AXIAL_RESOLUTION = 6.0;
      
    % Rescale the lateral and axial resolutions
    xResizeScale = LATERAL_RESOLUTION / parameters.X_RESOLUTION;
    yResizeScale = AXIAL_RESOLUTION / parameters.Y_RESOLUTION;
    
    lateralRes = parameters.X_RESOLUTION;
    axialRes = parameters.Y_RESOLUTION;
    
    % Resize the image
    image = imresize(image_1,'scale', [yResizeScale, xResizeScale]);
    
    image = normalizeValues(image,0,255);
    resizedSize = size(image);
    
    % Resize the invalid image
    invalidImage = imresize(invalidImage,'scale', [yResizeScale, xResizeScale]);
    invalidImage(invalidImage <= 0.5) = 0;
    invalidImage(invalidImage > 0.5) = 1;
    invalidIndices = find(invalidImage);
  

    %----------------------------------------------------------------------
    %  Get th++++++++++++
    % 1111111111111111111111111111111111111111-+
    % rs
    %----------------------------------------------------------------------

    % Get a black and white image of the hyper-reflective bands
    [bwImage, borders] = normal_getBwImage( ...
        image, axialRes, lateralRes, invalidIndices, parameters.getBwImageParams);
    layers = borders;
     rpeTop = borders(3,:);
     rpeBottom = borders(4,:);
    
    
    %----------------------------------------------------------------------
    % Flatten the image based on bruchs
    %----------------------------------------------------------------------
    
     invalidImage = zeros(size(image));
     invalidImage(invalidIndices) = 1;
    
    [image,pixelShift,invalidIndices] = normal_flattenImage(image, rpeTop);
    rpeTop = flattenImage(rpeTop, -pixelShift);
    rpeBottom = flattenImage(rpeBottom, -pixelShift);
     
     invalidImage(invalidIndices) = 1;
     invalidIndices = find(invalidImage);

    
    %----------------------------------------------------------------------
    %  Segment layers
    %----------------------------------------------------------------------

    layers = normal_graphCut( ...
        image, ...
        axialRes, ...
        lateralRes, ...
        parameters.otherParams.EYE, ...
        rpeTop, ...
        rpeBottom, ...
        invalidIndices, ...
        1, ...
        parameters.graphCutParams);

    
    %----------------------------------------------------------------------
    %  Unflatten the layers if they were flattened
    %----------------------------------------------------------------------
    
     for iLayer = 1:size(layers,1)  
         layers(iLayer,:) = flattenImage(layers(iLayer,:), pixelShift);
     end
   

    %----------------------------------------------------------------------
    % resample layers
    %----------------------------------------------------------------------
    
     layers = resampleLayers(layers,resizedSize,originalSize);
    
    
    %----------------------------------------------------------------------
    %   Perform layer smoothing
    %----------------------------------------------------------------------
      
    layers = smoothLayers( ...
        layers, ...
        originalSize(1), ...
        parameters.graphCutParams.SMOOTHING_CORRECTION);
    
    
    %----------------------------------------------------------------------
    %   Ignore the colums of the cut where the image in valid
    %----------------------------------------------------------------------

    invalidIndices = find(originalInvalidImage);
    for iLayer = 1:size(layers,1)
        x = 1:originalSize(2);
        invalidInd = intersect(sub2ind(originalSize,layers(iLayer,:),x), invalidIndices);
        [yInvalid,xInvalid] = ind2sub(originalSize,invalidInd);
        xInvalid = unique(xInvalid);
        layers(iLayer,xInvalid) = NaN;
        layers(isnan(layers))=1;
    end
    %layers = layers +29;
    lineThickness = 2;
    %% Plot the layers
   
    for iLayer = 1:size(layers,1)  
        figure(1),
        imshow(inputImage,[]),hold on
        color = GetLayerColor(iLayer);
        h = plot(layers(iLayer,:),color,'LineWidth',lineThickness);
        hold on
    end
    
  %%  
%     ilm_layer=layers(1,:);
%     %ilm_layer(isnan(ilm_layer))=1;
%     rpe_layer=layers(3,:);
   % rpe_layer(isnan(rpe_layer))=1;
 
 l1=layers(1,:);
 l2=layers(2,:);
 l3=layers(3,:);
 l4=layers(4,:);
 l5=layers(5,:);
 l6=layers(6,:);
 l7=layers(7,:);
 l8=layers(8,:);
    
   % ilmToRpelayerThickness(iImageSlice,:)=rpe_layer-ilm_layer;
   %csvwrite(fullfile(OutputDirectory,['Image',num2str(iImageSlice) '.csv']),layers );
   csvwrite(fullfile(OutputDirectory,sprintf('Image%03d.csv',iImageSlice)),layers);
    %return;
%     
%   plot(layers(1,:))
%   plot(layers(2,:))
%   plot(layers(3,:))
%   plot(layers(4,:))
%   plot(layers(5,:))
%   plot(layers(6,:))
%   plot(layers(7,:))
%   plot(layers(8,:))

%sprintf()
%saveas(gcf, 'output', 'tif')
%saveas(gcf,['Image',num2str(iImageSlice)], 'tif');
%sw= imshow(image,[]);
%test = getimage(sw);
 
% ilmToRpelayerThickness(isnan(ilmToRpelayerThickness))=0;
%csvwrite(fullfile(OutputDirectory,'Thickness','ilmToRpelayerThickness_1.csv'),ilmToRpelayerThickness );   


% labelImage = zeros (size (image_1));
% nAscan = size (image_1, 2);
% for iAscan = 1:nAscan
%   if (~isnan(l1(iAscan)) & ~isnan(l2(iAscan)))
% 	labelImage(l1(iAscan):l2(iAscan), iAscan) = 1;
%   end
%   if (~isnan(l2(iAscan)) & ~isnan(l3(iAscan)))
%     labelImage(l2(iAscan):l3(iAscan), iAscan) = 2;
%   end
%   if (~isnan(l3(iAscan)) & ~isnan(l4(iAscan)))
% 	labelImage(l3(iAscan):l4(iAscan), iAscan) = 1;
%   end
%   if (~isnan(l4(iAscan)) & ~isnan(l5(iAscan)))
%     labelImage(l4(iAscan):l5(iAscan), iAscan) = 2;
%   end
%   if (~isnan(l5(iAscan)) & ~isnan(l6(iAscan)))
% 	labelImage(l5(iAscan):l6(iAscan), iAscan) = 1;
%   end
%   if (~isnan(l6(iAscan)) & ~isnan(l7(iAscan)))
%     labelImage(l6(iAscan):l7(iAscan), iAscan) = 2;
%   end
%   if (~isnan(l7(iAscan)) & ~isnan(l8(iAscan)))
% 	labelImage(l7(iAscan):l8(iAscan), iAscan) = 1;
%   end
%   
% end
% if (doShowImageSegmentation)
%   figure(2); imagesc (labelImage);
% end

%  segmentationFileName = fullfile (OutputDirectory, sprintf('%s_segmentation.tif', listOfFiles(iImageSlice).name(1:end-4)));
%  imwrite (uint8(labelImage), segmentationFileName);

% if (doWriteOverlayImage)
%   overlayDirectory =   mkdir (OutputDirectory, 'overlay');
%   overlayImage(:,:,1) = image_1 + ((labelImage == 1) * 60);
%   overlayImage(:,:,2) = image_1 + ((labelImage == 2) * 60);
%   overlayImage(:,:,3) = image_1;
%   overlayFileName = fullfile (OutputDirectory, 'overlay', sprintf('%s_segmentation.tif', listOfFiles(iImageSlice).name(1:end-4)));
%   imwrite (uint8(overlayImage), overlayFileName);
%   imshow(uint8(overlayImage));
% end



end
%end

% M_ILMV1=max(ilmToRpelayerThickness(:));
% ILMV_n=M_ILMV1/255;
% ILMV_img=uint8(ilmToRpelayerThickness./ILMV_n);
% ILMV_img1=ind2rgb(ILMV_img,parula(255));
% ILMV_img2=im2uint8(ILMV_img1);
% %imshow(ILMV_img2)
% outputFile3=fullfile(OutputDirectory,'Thickness','IlmtoRpeThicknessMap.tif');
% imwrite(ILMV_img2,outputFile3);
% 
% %Normative Data
% IlmtoRpe_1=fullfile(thicknessDirectory,'ilmToRpelayerThickness.csv');
% IlmtoRpe_1=csvread(IlmtoRpe_1);
% IlmtoRpe_2=ilmToRpelayerThickness-IlmtoRpe_1;
% M_ILMV2=max(IlmtoRpe_2(:));
% ILMV_n2=M_ILMV2/255;
% ILMV_img_1=uint8(IlmtoRpe_2./ILMV_n2);
% ILMV_img11=ind2rgb(ILMV_img_1,parula(255));
% ILMV_img12=im2uint8(ILMV_img11);
% %imshow(ILMV_img12)
% outputFile6=fullfile(OutputDirectory,'Thickness','IlmtoRpeDeviationMap.tif');
% imwrite(ILMV_img12,outputFile6);


function color = GetLayerColor(layerNumber)

    layerNumber = mod(layerNumber,5);
    
    switch (layerNumber)
        case {1}    % blue
            color = 'b';
        case {2}    % magenta
            color = 'm';
        case {3}    % cyan
            color = 'c';
        case {4}    % yellow
            color = 'y';
        case {0}    % green
            color = 'g';
        otherwise
            color = [];
    end
end

function [image,pixelShift,invalidIndices] = normal_flattenImage( ...
    image, line)

    %----------------------------------------------------------------------
    %   Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(line) || all(isnan(line))
        pixelShift = [];
        invalidIndices = [];
        return;
    end
    
    [imageHeight,imageWidth] = size(image);
    
    
    %----------------------------------------------------------------------
    %   Get the line to flatten
    %----------------------------------------------------------------------
    
    % Extrapolate missing values from the line
    x = 1:imageWidth;
    validInd = ~isnan(line);
    line = smooth(line,0.1);
    line = round(interp1(x(validInd),line(validInd),x,'nearest','extrap'));
    
    % Make sure line is not out of the image bounds
    line(line < 1) = 1;
    line(line > imageHeight) = imageHeight;
    
    
    %----------------------------------------------------------------------
    %   Flatten the image
    %----------------------------------------------------------------------
    
    % Flatten the image based on the line
    [image, pixelShift, invalidIndices] = flattenImage(image,line);
end


function [bwImage, borders] = normal_getBwImage( ...
    image, axialRes, lateralRes, invalidIndices, params)

    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    if nargin < 4
        invalidIndices = [];
    end
    
    if nargin < 5
        params = [];
    end
    
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(image)
        error('Image cannot be empty');
    end
    
    if size(image, 3) > 1
        error('Image must be 2D');
    end
    
    if isempty(params)
        params = normal_getParameters();
        params = params.getBwImageParams;
    end


    %----------------------------------------------------------------------
    %  Process the image
    %----------------------------------------------------------------------
    
    image = double(image);
    imageSize = size(image);
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
    %
    %  Smooth the image to reduce noise
    %
    xFilterSize = round(params.X_FILTER_SIZE / lateralRes);
    yFilterSize = round(params.Y_FILTER_SIZE / axialRes);
    filter = fspecial('gaussian',[yFilterSize,xFilterSize],params.SIGMA);
    
    smoothImage = blurImage(image,filter);
    
    
    %----------------------------------------------------------------------
    %  Generate a binary image of the hyper-reflective bands
    %----------------------------------------------------------------------
    
    %
    %  Find the edges of the image
    %         
    bwImage = [zeros(1,imageWidth);diff(smoothImage)];
    bwImage(bwImage <= 0.5) = 0;
    bwImage(bwImage > 0.5) = 1;   
    bwImage(invalidIndices) = 0;

    %
    %  Open any gaps in the clusters
    %
    xStrelSize = round(params.X_STREL_SIZE / lateralRes);
    yStrelSize = round(params.Y_STREL_SIZE / axialRes);
    structuringElement = strel('rectangle',[yStrelSize, xStrelSize]);
    bwImage = imopen(bwImage, structuringElement);

    %
    %  Remove all clusters smaller than a certain size
    %
    clusters = bwconncomp(bwImage);
    minClusterSize = params.MIN_CLUSTER_SIZE / (axialRes*lateralRes);
    
    for iCluster = 1:clusters.NumObjects
        clusterInd = clusters.PixelIdxList{iCluster};
        clusterSize = length(clusterInd);
        if clusterSize < minClusterSize
           bwImage(clusterInd) = 0;
        end
    end

    %
    %  Close any gaps in the clusters
    %
    bwImage = imclose(bwImage, structuringElement);
    %----------------------------------------------------------------------
    %  Get the borders of the two hyper-reflective bands
    %----------------------------------------------------------------------
    
    % Add an additional column to each side of the image
    bwImage = addColumns(bwImage, 1);
    newImageSize = size(bwImage);
    newImageWidth = newImageSize(2);
    
    % Get the weighting matrices
    weightingMatrices = bwWeightingMatrix(bwImage);

    
    % Get Borders from Weighting Matrix
    % Cut each border
    lines = NaN(params.NUM_BORDERS,newImageWidth);
    invalidIndices = [];
    numBorders = params.NUM_BORDERS;
    
    for iBorder = 1:numBorders

        % Exclude the previous borders from the region to cut
        if iBorder > 1
            x = 2:(newImageWidth-1);
           % y=load('C:\Users\Mohan\source\repos\CUDA_GraphCutSegmentation\TEST.txt');
           % y = round(y)+1;
           % y = y(1,x);
            y = lines(iBorder-1,x);
            removeInd = (y == 1 | y == imageHeight);
            x(removeInd) = [];
            y(removeInd) = [];
            invalidIndices = [invalidIndices, sub2ind(newImageSize, y, x)];
        end
        
        if iBorder < params.NUM_BORDERS
            yBottom = imageHeight*ones(1,newImageWidth);
        else  
            yBottom = nanmax(lines);
        end

        % Get the valid region to cut
        regionIndices = getRegion( ...
            newImageSize, ...
            ones(1,newImageWidth), ...
            yBottom, ...
            0, ...
            0, ...
            invalidIndices);

        % Cut the region to get the border
        lines(iBorder,:) = cutRegion( ...
            newImageSize, ...
            regionIndices, ...
            weightingMatrices{mod(iBorder,2)+1});
    end
    
    % Remove the added columns
    lines1 = lines;
    bwImage = bwImage(:,2:end-1);
    lines = lines(:,2:end-1);
    
    % Sort the lines in ascending order
    oddIndices = 1:2:params.NUM_BORDERS;
    evenIndices = 2:2:params.NUM_BORDERS;
    oddSortOrder = sortrows([nanmean(lines(oddIndices,:),2), (1:length(oddIndices))']);
    evenSortOrder = sortrows([nanmean(lines(evenIndices,:),2), (1:length(evenIndices))']);
    bottomBorders = lines(oddIndices(oddSortOrder(:,2)),:);
    topBorders = lines(evenIndices(evenSortOrder(:,2)),:);
    
    borders(oddIndices,:) = topBorders;
    borders(evenIndices,:) = bottomBorders;
    
    % Replace extrapolated points (those that do not lie along a
    % hyper-reflective band) with NaN
    for iBorder = 1:params.NUM_BORDERS
        border = borders(iBorder,:);
        if ~mod(iBorder,2)
            border = border - 1;
            border(border < 1) = 1;
        end
        ind = sub2ind(imageSize,border,1:imageWidth);
        [yStart,xStart] = ind2sub(imageSize, ind(find(bwImage(ind), 1, 'first')));
        [yEnd,xEnd] = ind2sub(imageSize, ind(find(bwImage(ind), 1, 'last')));
        borders(iBorder,:) = NaN(1,imageWidth);
        borders(iBorder,xStart:xEnd) = border(xStart:xEnd);
    end
end


function regionIndices = normal_getGraphCutRegion( ...
    image, ...
    layerNumber, ...
    axialRes, ...
    eye, ...
    rpeTop, ...
    rpeBottom, ...
    foveaParams, ...
    layers)
    
    eye=1;
    %----------------------------------------------------------------------
    %  Get the layer parameters if no previous imageLayers were input
    %----------------------------------------------------------------------
    
    if isempty(eye) || (eye ~= 1 && eye ~= 2)
        error('Eye must have a value of 1 or 2');
     end
    
    %----------------------------------------------------------------------
    %  Get the layer parameters if no previous imageLayers were input
    %----------------------------------------------------------------------

    imageSize = size(image);
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
    topLayerAddition = 0;
    bottomLayerAddition = 0;
    correctBottomLine = 0;
        
    switch (layerNumber)

        % vitreous-NFL
        case {1}            
            yTop = ones(1, imageWidth);
            yBottom = rpeTop;
            correctBottomLine = 1;

        % cut RPE-Choroid
        case {8}
            if all(isnan(layers(8,:))) 
                if all(isnan(rpeTop))
                    yTop = layers(1,:) + round(100.5/axialRes);
                else
                    yTop = rpeTop;
                end
                yBottom = rpeBottom;
            else
                yTop = layers(8,:) + round(13.4/axialRes);
                yBottom = layers(8,:) + round(67/axialRes);
            end

        % cut OS-RPE
        case {7}
            yTop = layers(8,:) - round(33.5/axialRes);
            yBottom = layers(8,:) - round(13.4/axialRes);

        % cut IS-OS
        case {6}
            yTop = layers(7,:) - round(67/axialRes);
            yBottom = layers(7,:) - round(20.1/axialRes);
            
        % cut 
        case {4}
            if all(isnan(layers(5,:)))
                difference = layers(6,:) - layers(1,:);
                yTop1 = layers(1,:) + round(difference/3);
                yTop2 = layers(6,:) - round(134/axialRes);
                yTop3 = layers(1,:) + round(13.4/axialRes);
                yTop = nanmax(nanmin(yTop1,yTop2),yTop3);
                yBottom = layers(6,:) - round(40.2/axialRes);
            else
                yTop = layers(3,:);
                yBottom = layers(5,:);
                correctBottomLine = 1;
            end

        % cut 
        case {3}
            yTop = layers(4,:) - round(46.9/axialRes);
            yBottom = layers(4,:);
            
        % cut 
        case {5}
            yTop = round(smooth(layers(4,:),0.1)');
            yBottom1 = yTop + round(100.5/axialRes);
            yBottom2 = layers(6,:) - round(13.4/axialRes);
            yBottom = nanmin(yBottom1,yBottom2);
            
            if ~isempty(foveaParams.Index)
                range = foveaParams.Range(1):foveaParams.Index;
                yFovea = layers(3,foveaParams.Index) - round(13.4/axialRes);
                thickness = 0:length(range)-1;
                thickness = round((yFovea-yTop(range(1)))/length(range) * thickness);
                yTop(range) = yTop(range(1)) + thickness;
                
                range = foveaParams.Index:foveaParams.Range(2);
                yFovea = layers(3,foveaParams.Index) - round(13.4/axialRes);
                thickness = 0:length(range)-1;
                thickness = round((yTop(range(end))-yFovea)/length(range) * thickness);
                yTop(range) = yFovea + thickness;
            end

        % cut 
        case {2}
            yTop = layers(1,:);
            yBottom = layers(3,:) - round(26.8/axialRes);

            if ~isempty(foveaParams.Index)
                range = foveaParams.Range(1):foveaParams.Range(2);
                yBottom(range) = layers(1,range) + round(20.1/axialRes);

                if eye == 1
                    thinRange = range(end):imageWidth;
                    thickRange = 1:range(1);
                elseif eye == 2
                    thinRange = 1:range(1);
                    thickRange = range(end):imageWidth;
                end
                thickness = 0:length(thickRange)-1;
                thickness = round(round(33.5/axialRes) / length(thickRange) * thickness);
                yTop(thickRange) = yTop(thickRange) + thickness;
                yBottom(thickRange) = layers(3,thickRange) - fliplr(thickness);
                yBottom(thinRange) = yTop(thinRange) + round(26.8/axialRes);
                
                correctBottomLine = 1;
            else
                range = 1:imageWidth;
                if eye == 1
                    range = fliplr(range);
                end
                thickness = 0:length(range)-1;
                thickness = round(round(13.4/axialRes) / length(range) * thickness);
                yTop(range) = yTop(range) + thickness;
            end
            
        otherwise
            error('Segmentation of layer %d is not supported', layerNumber);
    end

    
    %----------------------------------------------------------------------
    %  Handle cases where the input layer is NaN or empty.
    %----------------------------------------------------------------------
 
    nanIndices = isnan(yTop);   
    if sum(nanIndices) > 0
        yTop(nanIndices) = ones(1,sum(nanIndices));
    end
    
    nanIndices = isnan(yBottom);   
    if sum(nanIndices) > 0
        yBottom(nanIndices) = imageHeight*ones(1,sum(nanIndices));
    end


    %----------------------------------------------------------------------
    %  Make the region smaller with respect to the selected minimum 
    %  distance between lines
    %----------------------------------------------------------------------

    regionIndices = getRegion( ...
        imageSize, ...
        yTop, ...
        yBottom, ...
        topLayerAddition, ...
        bottomLayerAddition, ...
        [], ...
        correctBottomLine);
end





function object = addColumns(object, numColumns, fillValue)
    
    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    if nargin < 3
        fillValue = -1;
    end
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------

    if isempty(object)
        error('Layers cannot be empty');
    end

    %----------------------------------------------------------------------
    %  Add columns to each side, duplicating the rightmost and leftmost
    %  values
    %----------------------------------------------------------------------

    if fillValue == -1
        leftColumns = repmat(object(:,1), 1, numColumns);
        rightColumns = repmat(object(:,end), 1, numColumns);
    else
        leftColumns = fillValue*ones(size(object,1),numColumns);
        rightColumns = leftColumns;
    end
    object = [leftColumns, object, rightColumns];
end


function weightingMatrices = bwWeightingMatrix(image, matrixIndices, params)
    
    
    %----------------------------------------------------------------------
    %  Initialize missing parameters
    %----------------------------------------------------------------------
    
    if nargin < 2
        matrixIndices = [];
    end
    
    if nargin < 3
        params = [];
    end
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(params)
       % params = WeightingMatrixParameters;
    end    

    imageSize = size(image);
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
   
    %----------------------------------------------------------------------
    %  Generate a weighting matrix for the image
    %----------------------------------------------------------------------
    
    %
    %  Generate edges and points on the image
    %
    edges = createLattice(imageSize);
    
    %
    %  Separate the edges corresponding to the left-most and right-most
    %  columns that were artificially added to the image
    %
    maxIndex = imageHeight*imageWidth;
    leftColIndices = 1:(imageHeight-1);
    rightColIndices = ((imageHeight-1)*(imageWidth-2) + 1) : maxIndex;
    
    columnIndices = [leftColIndices, rightColIndices];
    imageIndices = setdiff(1:size(edges,1), columnIndices);
    
    columnEdges = edges(columnIndices,:);
    imageEdges = edges(imageIndices,:);
    
    
    %----------------------------------------------------------------------
    %  Calculate the weights based on the image gradient.  Lower weights 
    %  are assigned to areas with a higher gradient
    %----------------------------------------------------------------------
        
    %
    %  Create two edge maps (one for edges that transition from dark->light
    %  in the vertical direction, and one for edges transitioning from 
    %  light->dark).  
    %
    
   
    diffImage = diff([zeros(1,imageWidth); image]);
    
    
    lightDarkEdgeImage = double((diffImage > 0));
    darkLightEdgeImage = double(abs(diffImage < 0));
    
    ind = (lightDarkEdgeImage == 0) & (image == 1);
    lightDarkEdgeImage(ind) = -1;
    
    ind = (darkLightEdgeImage == 0) & (image == 1);
    darkLightEdgeImage(ind) = -1;
    
    %
    %  Calculate the gradient weights for each of the edge maps
    %    
    
    
    lightDarkGradientWeights = 2 - ...
        lightDarkEdgeImage(imageEdges(:,1)) - ...
        lightDarkEdgeImage(imageEdges(:,2));
    
    darkLightGradientWeights = 2 - ...
        darkLightEdgeImage(imageEdges(:,1)) - ...
        darkLightEdgeImage(imageEdges(:,2));
    
    
    %----------------------------------------------------------------------
    %  Create a weighting matrix. Rows represent the indices of the first
    %  node and columns represent the indices of the second node.  The 
    %  values represent the weights for the edge that is formed between
    %  the two nodes
    %----------------------------------------------------------------------
    
    % Define the matrices
    weights = { ...
        {lightDarkGradientWeights}, ...  % Dark-light
        {darkLightGradientWeights},...   % Light-dark
    };
    nMatrices = length(weights);  
    WEIGHT_RANGES=[];
    % Set default weight ranges if necessary   
    if isempty(WEIGHT_RANGES)
        WEIGHT_RANGES = {{[0,1]},{[0,1]}};
    end

    
    % Remove invalid matrix indices
    matrixIndices(matrixIndices < 1) = [];
    matrixIndices(matrixIndices > nMatrices) = [];
    if isempty(matrixIndices)
        matrixIndices = 1:nMatrices;
    end
    
    % Populate the matrices specified
    weightingMatrices = cell(nMatrices,1);
    matrixSize = maxIndex;
    
    for iMatrix = 1:length(matrixIndices)
        
        matrixIndex = matrixIndices(iMatrix);
        
        weightingMatrices{matrixIndex} = generateWeightingMatrix( ...
            matrixSize, ...
            imageEdges, ...                      
            weights{matrixIndex}, ...
            WEIGHT_RANGES{matrixIndex}, ...
            columnEdges, ...
            MIN_WEIGHT);
    end
end


function edges = createLattice(imageSize)

    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);

    %----------------------------------------------------------------------
    %  Get edge pairs
    %----------------------------------------------------------------------
    
    % Create a matrix containing the indices of each pixel
    indexMatrix = reshape(1:imageHeight*imageWidth, imageHeight, imageWidth);
    
    % Get (undirectional) vertical edges
    startNodes = indexMatrix(1:end-1,:);
    endNodes = indexMatrix(2:end,:);    
    edges = [startNodes(:), endNodes(:)];
    
    % Get (undirectional) horizontal edges
    startNodes = indexMatrix(:,1:end-1);
    endNodes = indexMatrix(:,2:end);
    edges = [edges; startNodes(:), endNodes(:)];
    
    % Get (undirectional) diagonal edges
    startNodes = indexMatrix(1:end-1,1:end-1);
    endNodes = indexMatrix(2:end,2:end);
    edges = [edges; startNodes(:), endNodes(:)];

    startNodes = indexMatrix(2:end,1:end-1);
    endNodes = indexMatrix(1:end-1,2:end);
    edges = [edges; startNodes(:), endNodes(:)];
end



function [cut,distance] = cutRegion( ...
    imageSize, ...
    regionIndices, ...
    weightingMatrix, ...
    coordinateIndices)

    %----------------------------------------------------------------------
    %  Initialize missing parameters
    %----------------------------------------------------------------------
    
    if nargin < 4
        coordinateIndices = [];
    end
    
    
    %----------------------------------------------------------------------
    %  Verify input parameters
    %----------------------------------------------------------------------

    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    regionIndices = sort(regionIndices);
    
    %
    %  Verify that all input coordinates are within the region 
    %  specified
    %
    if ~isempty(coordinateIndices)
        for iCoord = 1:length(coordinateIndices)
            if (0 == sum(regionIndices == coordinateIndices(iCoord)))
                error('Coordinates must be within the input region');
            end
        end
    end       

    %
    %  Make sure the region spans the entire width of the image
    %
    [y,x] = ind2sub(imageSize,regionIndices);   
    
               startIndex = regionIndices(1);
              temp= x
    endIndex = find(x == imageWidth,1,'first');
  
    if startIndex > imageHeight || isempty(endIndex)
        error('Region does not span the entire width of the image');
    end
    
    %
    %  Make sure the coordinate indices span the entire width of the image
    %
    if isempty(coordinateIndices) || coordinateIndices(1) > imageHeight
        coordinateIndices = [startIndex, coordinateIndices];
    end
    
    
    [y,x] = ind2sub(imageSize,coordinateIndices(end));
    
    if isempty(coordinateIndices) || x < imageWidth
        endIndex = regionIndices(endIndex);
        coordinateIndices = [coordinateIndices, endIndex];
    end
    
    nCoordinates = length(coordinateIndices);
    
    
    %----------------------------------------------------------------------
    %  Restrict the graph cut region to the regionIndices input
    %----------------------------------------------------------------------
        
    %
    %  Restrict the adjacency matrix based on the region indices
    %
    weightingMatrix = weightingMatrix(regionIndices, regionIndices);
        
    %    
    %  Generate lookup matrix, creating indices for the new region.
    %
    region = zeros(imageSize);
    region(regionIndices) = 1:length(regionIndices);
    
    
    %----------------------------------------------------------------------
    %  Calculate the best path that will connect the input coordinates
    %----------------------------------------------------------------------
    
    path = [];
    
    for iCoord = 1:nCoordinates - 1
        
        firstIndex = coordinateIndices(iCoord);
        secondIndex = coordinateIndices(iCoord + 1);
        
        %
        %  Find the shortest path two coordinates at a time
        %
    	[distance, pathSegment] = graphshortestpath( ...
            weightingMatrix, ...
            region(firstIndex), ...
            region(secondIndex), ...
            'Method', 'Dijkstra', ...
            'Directed', 0);
        
   % pathSegment=load('path.txt')';
   %   pathSegment = pathSegment + 1;
  %   pathSegment = flip(pathSegment);
        %  Add the paths together. Remove the repeated coordinate
        %
        if isempty(pathSegment)
            path = [];
            break;
        elseif iCoord ~= 1
            pathSegment(1) = [];
        end
        
        path = [path, pathSegment];
    end
    % path=load('C:\Users\Mohan\source\repos\CUDA_GraphCutSegmentation\TEST.txt');
     %path = path+1;
    [y,x] = ind2sub(imageSize,regionIndices(path));
    
    temp =x;
    %
    %  Since the path contains multiple points per column, take the first
    %  point in every column as the path
    %
    cut = NaN(1, imageWidth);
    
    if ~isempty(x) && ~isempty(y)
        for column = 1:imageWidth

            if column == 1
                index = find(x == column, 1, 'last'); 
            else
                index = find(x == column, 1, 'first'); 
            end

            cut(column) = y(index);
        end
    end
end


function matrix = generateWeightingMatrix( ...
    matrixSize, ...
    imageEdges, ...
    weights, ...
    weightRanges, ...
    columnEdges, ...
    columnWeight)

    %
    %  Set the column weights
    %
    columnWeights = columnWeight .* ones(length(columnEdges),1);    
    
    %
    %  Normalize the image weights and combine all weights together
    %
    imageWeights = columnWeight;
    
    for index = 1:length(weights)
        weights{index} = normalizeValues( ...
            weights{index}, weightRanges{index}(1), weightRanges{index}(2));
        
        imageWeights = imageWeights + weights{index};
    end 
                    
    %
    %  Combine the image and column weights, adding node pair duplicates
    %  such that the paths are bidirectional
    %
    totalEdges = [imageEdges; columnEdges];
    totalWeights = [imageWeights; columnWeights];
    
    
    
    %  Create a weighting matrix. Rows represent the indices of the first
    %  node and columns represent the indices of the second node.  The 
    %  values represent the weights for the edge that is formed between
    %  the two nodes   
    %  
    
%     U=load('C:\Users\Mohan\source\repos\Dijkstra_AdjacencyList\U.txt');
%     U = U +1;
%     V=load('C:\Users\Mohan\source\repos\Dijkstra_AdjacencyList\V.txt');
%     V = V+1;
%     totalWeights =load('C:\Users\Mohan\source\repos\Dijkstra_AdjacencyList\WeightValues.txt');
%     totalWeights = totalWeights';
    matrix = sparse( ...
         totalEdges(:,2), totalEdges(:,1), ...
        totalWeights, ...
        matrixSize, matrixSize);
end        


    function regionIndices = getRegion( ...
    imageSize, ...
    topLayer, ...
    bottomLayer, ...
    topAddition, ...
    bottomAddition, ...
    invalidIndices, ...
    correctBottomLine, ...
    spanImageWidth)
    
        
    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    if nargin < 4
        topAddition = 0;
    end
    if nargin < 5
        bottomAddition = 0;
    end
    if nargin < 6
        invalidIndices = [];
    end
    if nargin < 7
        correctBottomLine = 1;
    end
    if nargin < 8
        spanImageWidth = 1;
    end
        
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------

    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
    if length(topLayer) ~= length(bottomLayer)
        error('topLayer and bottomLayer must be of the same length');
    end
    
    if length(topLayer) ~= imageWidth
        error('The layers must have a length equal to the image width');
    end
    
    % Replace all layer values outside the image range
    topLayer(topLayer < 1) = 1;
    topLayer(topLayer > imageHeight) = imageHeight;
    bottomLayer(bottomLayer < 1) = 1;
    bottomLayer(bottomLayer > imageHeight) = imageHeight;
 
    %----------------------------------------------------------------------
    %  Expand the each layer boundary by the number of pixels to add,
    %  making sure to take care of any pixels that are out of bounds
    %----------------------------------------------------------------------
    
    bottomLayer = bottomLayer + bottomAddition;    
    topLayer = topLayer + topAddition;
    
    
    %----------------------------------------------------------------------
    %  Limit the layers by the invalid region
    %----------------------------------------------------------------------
    
    invalidImage = zeros(imageSize);
    invalidImage(invalidIndices) = 1;
    
    
    for iCol = 1:imageWidth
        topIndex = find(invalidImage(:,iCol) == 0, 1, 'first');
        bottomIndex = find(invalidImage(:,iCol) == 0, 1, 'last');
        if ~isempty(topIndex)
            topLayer(iCol) = max(topLayer(iCol),topIndex);
            bottomLayer(iCol) = min(bottomLayer(iCol),bottomIndex);
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Correct the appropriate line if it crosses the other line
    %----------------------------------------------------------------------
    
    difference = bottomLayer - topLayer;
    invalidInd = find(difference < 0);
    if ~isempty(invalidInd)
        if correctBottomLine
            bottomLayer(invalidInd) = topLayer(invalidInd);
        else
            topLayer(invalidInd) = bottomLayer(invalidInd);
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Get the indices of all pixels in between the two regions  
    %----------------------------------------------------------------------
    
    regionImage = zeros(imageSize);
    
     for iCol = 1:imageWidth
        
        % Close any vertical gaps that there may be in the region
        if (iCol < imageWidth)
            if topLayer(iCol) > bottomLayer(iCol+1)
                topLayer(iCol) = topLayer(iCol+1);
                bottomLayer(iCol+1) = bottomLayer(iCol);
                
            elseif bottomLayer(iCol) < topLayer(iCol+1)
                bottomLayer(iCol) = bottomLayer(iCol+1);
                topLayer(iCol+1) = topLayer(iCol);
            end
        end
        
        % Get the indices in the region
        yRegion = topLayer(iCol):bottomLayer(iCol);
        indices = round(sub2ind(imageSize, yRegion, iCol*ones(size(yRegion))));
        if ~isnan(indices)
            regionImage(indices) = 1;
            
        % Make sure the region extends across the width of the image
        elseif spanImageWidth
            regionImage(:,iCol) = 1;
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Take out any region indices that were specified as invalid
    %----------------------------------------------------------------------
      
    invalidImage = zeros(imageSize);
    invalidImage(invalidIndices) = 1;
    invalidImage = invalidImage & regionImage;
    
    % Remove invalid indices from the region
    regionImage(invalidIndices) = 0;
    
    % Make sure no columns are all NaN
    nanColumns = sum(regionImage) ==  0;
    regionImage(:,nanColumns) = invalidImage(:,nanColumns);
    nanColumns = sum(regionImage) ==  0;
    validImage = ones(imageSize);
    validImage(invalidIndices) = 0;
    regionImage(:,nanColumns) = validImage(:,nanColumns);
    nanColumns = sum(regionImage) ==  0;
    regionImage(:,nanColumns) = 1;
    
    % Get the indices    
    regionIndices = find(regionImage == 1);
end


function layers = normal_graphCut( ...
    image, ...
    axialRes, ...
    lateralRes, ...
    eye, ...
    rpeTop, ...
    rpeBottom, ...
    invalidIndices, ...
    recalculateWeights, ...
    params)
    
    
    %----------------------------------------------------------------------
    %  Persistent variables
    %----------------------------------------------------------------------
    
    persistent weightingMatrices;
    
    
    %----------------------------------------------------------------------
    %  Initialize missing parameters
    %----------------------------------------------------------------------
    
    if nargin < 5
        rpeTop = [];
    end    
    if nargin < 6
        rpeBottom = [];
    end    
    if nargin < 7
        invalidIndices = [];
    end
    if nargin < 8
        recalculateWeights = 1;
    end
    if nargin < 9
        params = [];
    end
    
    
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(params)
        params = normal_getParameters();
        params = params.graphCutParams;
    end
    
    numLayers = params.NUM_LAYERS;
    %numLayers=3;
    if numLayers < 1
        error('Number of layers must be greater than zero');
    end
    
    if params.MAX_NUM_LAYERS < numLayers
        error('Max number of layers cannot be less than the number of layers');
    end
    if numLayers < 1 || numLayers > params.MAX_NUM_LAYERS
        error('Layer indices must be greater than zero');
    end
    
    if isempty(image)
        error('Image cannot be empty');
    end

    
    %----------------------------------------------------------------------
    %  Normalize image
    %----------------------------------------------------------------------

    image = double(image);
    image = normalizeValues(image);  
    imageHeight = size(image,1);
    

    %------------------------------------------------------------------
    %  Add a column to each side of the images.  These are used to
    %  initialize the graph cut starting points
    %------------------------------------------------------------------

    image = addColumns(image, 1);

    if ~isempty(rpeTop)
        rpeTop = addColumns(rpeTop, 1);
    end

    if ~isempty(rpeBottom)
        rpeBottom = addColumns(rpeBottom, 1);
    end

    if ~isempty(invalidIndices)
        invalidIndices = invalidIndices + imageHeight;
    end
    
    
    %------------------------------------------------------------------
    %  Determine which layers to segment
    %------------------------------------------------------------------

    layerIndices = params.LAYER_INDICES;
    matrixIndices = params.MATRIX_INDICES;
    
    
    maxIndex = 0;
    uniqueLayerIndices = [];
    for index = 1:length(layerIndices)
        if ~sum(uniqueLayerIndices == layerIndices(index))
            uniqueLayerIndices = [uniqueLayerIndices, layerIndices(index)];
        end
    end
    
    for index = 1:numLayers
        lastIndex = find(layerIndices == uniqueLayerIndices(index),1,'last');
        if ~isempty(lastIndex) && lastIndex > maxIndex
            maxIndex = lastIndex;
        end
    end
    
    if maxIndex > 0
        layerIndices = layerIndices(1:maxIndex);
        if ~isempty(matrixIndices)
            matrixIndices = matrixIndices(1:maxIndex);
        end
    end
    
    
    %------------------------------------------------------------------
    %  Generate a weighting matrix for the image
    %------------------------------------------------------------------

    if recalculateWeights
        weightingMatrices = normal_weightingMatrix( ...
            image, ...
            axialRes, ...
            lateralRes, ...
            matrixIndices, ...
            invalidIndices, ...
            params.weightingMatrixParams);
    end
      
    
    %----------------------------------------------------------------------
    %  Segment each layer one at a time
    %----------------------------------------------------------------------
    
    imageSize = size(image);
    nLayers = max(layerIndices);
    layers = NaN(nLayers, imageSize(2));
    foveaParams = struct('Index', 0, 'Range', 0, 'Percentage', 0);
    
    % Loop through each layer and segment it
    for iLayer = 1:length(layerIndices)
        
        layerIndex = layerIndices(iLayer);
        
        % Get parameters for the current layer to segment
        regionIndices = normal_getGraphCutRegion( ...        
            image, ...
            layerIndex, ...
            axialRes, ...
            eye, ...
            rpeTop, ...
            rpeBottom, ...
            foveaParams, ...
            layers);
            
        % Assign the appropriate adjacency matrix
        matrixIndex = matrixIndices(iLayer);
        weightMatrix = weightingMatrices{matrixIndex};
        
        % Perform graph cut between each of the points, and combine them
        % together to create a continuous segmented line across the image
        cut = cutRegion(imageSize, regionIndices, weightMatrix);
        
        % See if there is a fovea present
        if layerIndex == 3
            ilm = layers(1,:);
            innerLayer = round(smooth(cut,0.1)');
            rpe = layers(6,:);
            foveaParams = locateFovea(imageSize, axialRes, lateralRes, ...
                ilm, innerLayer, rpe, invalidIndices);   
        end
        
        layers(layerIndex,:) = cut;
    end
	   

    %----------------------------------------------------------------------
    %  Make sure the layers do not cross
    %----------------------------------------------------------------------
    
    for iLayer = 1:nLayers
        if iLayer < nLayers
            topLayer = layers(iLayer, :);
            bottomLayer = layers(iLayer + 1, :);
            crossIndices = find((bottomLayer - topLayer) < 0);

            bottomLayer(crossIndices) = topLayer(crossIndices);
            layers(iLayer + 1,:) = bottomLayer;
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Remove the extra columns  and layers
    %----------------------------------------------------------------------
    
    layers = layers(:,2:end-1);
    
    % Remove extra layers that were not segmented or supposed to be
    % segmented
    layersToRemove = setdiff(1:nLayers,uniqueLayerIndices(1:numLayers));
    layers(layersToRemove,:) = [];
end


function foveaParams = locateFovea( ...
    imageSize, ...
    axialRes, ...
    lateralRes, ...
    ilm, ...
    innerLayer, ...
    rpe, ...
    invalidIndices)

    % Only look for the fovea in the middle of the image
    imageWidth = imageSize(2);
    margin = round(imageWidth/8);
    middleInd = zeros(1,imageWidth);
    middleInd(margin:imageWidth-margin) = 1;
    
    % Do not look at the fovea where the image is invalid
    x = 1:imageWidth;
    invalidInd = intersect(sub2ind(imageSize,[ilm,rpe],[x,x]), invalidIndices);
    [yInvalid,xInvalid] = ind2sub(imageSize,invalidInd);
    xInvalid = unique(xInvalid);
    middleInd(xInvalid) = 0;
            
    % Get the layer thickness
    thicknessOrig = innerLayer - ilm;
    thickness = smooth(thicknessOrig,0.1);
    
    % Locate where the locale max and min layer thicknesses are
    maxima = find( middleInd(2:end)' & diff( sign( diff([thickness(1); thickness]) ) ) < 0 );
    maxima(diff(maxima) == 1) = [];
    minima = find( middleInd(2:end)' & diff( sign( diff([thickness(1); thickness;]) ) ) > 0 );
    minima(diff(minima) == 1) = [];
    combined = [maxima, zeros(length(maxima),1); minima, ones(length(minima),1)];
    combined = sortrows(combined);
    
    % Look for the largest change in thickness (corresponding to the fovea
    % location)
    difference = abs(diff(thickness(combined(:,1))));
    locations = find(difference > round(33.5/axialRes));
    
    foveaParams.Index = [];
    foveaParams.Range = [];
    foveaParams.Percentage = [];
           
    for index = 1:length(locations)-1
        if (combined(locations(index),2) == 0 && combined(locations(index)+1,2) == 1) && ...
           (combined(locations(index+1),2) == 1 && combined(locations(index+1)+1,2) == 0)
    
            foveaRange = [combined(locations(index),1), combined(locations(index+1)+1,1)];
            foveaThick = thicknessOrig(foveaRange(1):foveaRange(2));
            foveaIndex = find(foveaThick == min(foveaThick));
            foveaIndex = foveaIndex(round(length(foveaIndex)/2)) + foveaRange(1) - 1;                   
            foveaPercentage = (innerLayer(foveaIndex) - ilm(foveaIndex)) / ...
                              (rpe(foveaIndex) - ilm(foveaIndex));
                          
            % Make the fovea range a constant width
            thickness = round(670/lateralRes);
            foveaRange(1) = max(1,foveaIndex-thickness);
            foveaRange(2) = min(imageWidth,foveaIndex+thickness);

            foveaParams.Index = foveaIndex;
            foveaParams.Range = foveaRange;
            foveaParams.Percentage = foveaPercentage;
            break;
        end
    end
end

function layers = smoothLayers( ...
    layers, ...
    imageHeight, ...
    smoothness, ...
    roundLayers)

    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    if nargin < 4
        roundLayers = 1;
    end

    
    %----------------------------------------------------------------------
    %  Verify input parameters
    %----------------------------------------------------------------------

    if isempty(smoothness)
        error('''smoothness'' cannot be empty');
    end
    
    nLayers = size(layers, 1);
    
    % If only one smoothness value was given, use that value for all layers
    if length(smoothness) == 1
        smoothness = smoothness * ones(1,nLayers);
        
    % If not enough smoothness values were given, do not smooth the layers
    % not given a value
    elseif length(smoothness) < nLayers
        temp = smoothness;
        smoothness = zeros(1,nLayers);
        smoothness(1:length(temp)) = temp;
    end
    
    
    %----------------------------------------------------------------------
    %  Smooth each layer individually
    %----------------------------------------------------------------------

    for iLayer = 1:nLayers

        layer = layers(iLayer,:);
        
        if smoothness(iLayer) > 0
            layer = smooth(layer, smoothness(iLayer), 'sgolay')';
            
            if roundLayers
                layer = round(layer);
            end
        end
        
        % Take care of out of bounds scenaries
        layer(layer < 1) = 1;
        layer(layer > imageHeight) = imageHeight;
        
        % Make sure the lower layers do not cross the layers above it
        if (iLayer > 1)
            crossIndices = (layer - layers(iLayer-1,:)) < 0;
            layer(crossIndices) = layers(iLayer-1, crossIndices);
        end
        
        layers(iLayer,:) = layer;
    end
end

function [image, borderImage, borderLines] = removeImageBorders( ...
    image, borderValues, replaceType, borderType)

    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------

    if nargin < 3 || isempty(replaceType)
        replaceType = 1;
    end
    if nargin < 4 || isempty(borderType)
        borderType = 0;
    end
        
    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(image)
        error('Image cannot be empty');
    end

    image = double(image);
    imageSize = size(image);
    imageWidth = imageSize(2);
    
    
    %----------------------------------------------------------------------
    %  Look for boundaries
    %----------------------------------------------------------------------
    
    borderLines = NaN(2,imageWidth);
    borderImage = false(imageSize);
    
    %
    %  First iteration searches through each column of the image, 
    %  second iteration searches through each row of the image,
    %  looking for the valid region (the region inside the borders)
    %
    for i = 1:2
        
        if borderType == 2
            continue;
        end
        
        %
        %  Search through each column of the image, looking for non-border
        %  regions
        %
        imageHeight = size(image,1);
        imageWidth = size(image,2);
    
        for iColumn = 1:imageWidth

            column = image(:,iColumn);
            borderColumn = zeros(size(column));

            % Locate where the border ends on the top of the image, and
            % where it starts on the bottom of the image
            for value = borderValues
                borderColumn(column == value) = 1;
            end

            startIndex = find(borderColumn ~= 1, 1, 'first');
            endIndex = find(borderColumn ~= 1, 1, 'last');

            % The full column is a border, so move on
            if isempty(startIndex) && isempty(endIndex)
                if i == 1
                    borderLines(1,iColumn) = imageHeight;
                    borderLines(2,iColumn) = 1;
                end
                continue;
            end                

            %--------------------------------------------------------------
            %  Address the top border
            %--------------------------------------------------------------
            
            % If no border was found at the top of the image, then move on
            if isempty(startIndex) || startIndex == 1
                if i == 1
                    borderLines(1,iColumn) = NaN;
                end
                
            % Otherwise a border was found
            else
                borderIndices = 1:startIndex-1;
                
                % Replace the border with the first non-border pixel
                if replaceType == 2
                    image(borderIndices,iColumn) = image(startIndex,iColumn);  
                    
                % Otherwise replace the border with the mirror image of the
                % non-border pixels
                else
                    startCol = startIndex;
                    endCol = startCol + length(borderIndices) - 1;
                    
                    if endCol <= imageHeight
                        image(borderIndices,iColumn) = flipud(image(startCol:endCol,iColumn));
                        
                    % During the first iteration, if there arent enough
                    % non-border pixels for mirror imaging, then save it
                    % for the next iteration
                    elseif i == 1
                        image(borderIndices,iColumn) = borderValues(1);
                    
                    % If there are not enough valid non-border pixels for
                    % the mirror image, then just duplicate the last pixel
                    else                        
                        replaceCol = image(startCol:imageHeight,iColumn);
                        missingLength = length(borderIndices) - length(replaceCol);
                        replaceCol = [replaceCol; replaceCol(end)*ones(missingLength,1)];
                        image(borderIndices,iColumn) = flipud(replaceCol);
                    end
                end
            
                % Assign the top border location for the given column
                if i == 1
                    borderLines(1,iColumn) = startIndex-1;
                end
                borderImage(1:startIndex-1,iColumn) = 1;
            end

            %--------------------------------------------------------------
            %  Address the bottom border
            %--------------------------------------------------------------
            
            % If no border was found at the bottom of the image, then move 
            % on
            if isempty(endIndex) || endIndex == imageHeight
                if i == 1
                    borderLines(2,iColumn) = NaN;
                end
                
            % Otherwise a border was found
            else
                borderIndices = endIndex+1:imageHeight;
                
                % Replace the border with the first non-border pixel
                if replaceType == 2
                    image(borderIndices,iColumn) = image(endIndex,iColumn);  
                    
                % Otherwise replace the border with the mirror image of the
                % non-border pixels
                else
                    endCol = endIndex;
                    startCol = endCol - length(borderIndices) + 1;

                    if startCol >= 1
                        image(borderIndices,iColumn) = flipud(image(startCol:endCol,iColumn));
                        
                    % During the first iteration, if there arent enough
                    % non-border pixels for mirror imaging, then save it
                    % for the next iteration
                    elseif i == 1
                        image(borderIndices,iColumn) = borderValues(1);
                    
                    % If there are not enough valid non-border pixels for
                    % the mirror image, then just duplicate the last pixel
                    else                        
                        replaceCol = image(1:endCol,iColumn);
                        missingLength = length(borderIndices) - length(replaceCol);
                        replaceCol = [replaceCol(1)*ones(missingLength,1);replaceCol];
                        image(borderIndices,iColumn) = flipud(replaceCol);
                    end
                end
            
                % Assign the bottom border location for the given column
                if i == 1
                    borderLines(2,iColumn) = endIndex+1;
                end
                borderImage(endIndex+1:end,iColumn) = 1;
            end
        end
                
        if borderType == 1
            break;
        end
        
        image = image';
        borderImage = borderImage';
    end
end


function values = normalizeValues(values, minValue, maxValue)

    %----------------------------------------------------------------------
    %  Initialize missing parameters
    %----------------------------------------------------------------------

    if nargin < 2
        minValue = 0;
    end
    
    if nargin < 3
        maxValue = 1;
    end

    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if minValue > maxValue
        error('''minValue'' must be <= ''maxValue''');
    end

    %----------------------------------------------------------------------
    %  Normalize the values
    %----------------------------------------------------------------------
    
    %
    %  Convert values to double for more accurate results
    %
    values = double(values);
    
    %
    %  If the new minima and maxima are the same value, then set the whole
    %  matrix to that value
    %
    if minValue == maxValue
        values(:) = minValue;
        return;
    end
    
    %
    %  Find old minima and maxima
    %
    oldMinValue = min(values(:));
    oldMaxValue = max(values(:));
    
    %
    %  If the new minima and maxima are same as the old ones, then just
    %  return the existing matrix
    %    
    if (oldMinValue == minValue && oldMaxValue == maxValue)
        return;
    end
    
    %
    %  If all old values are the same, just return the new minimum value
    %
    if (oldMinValue == oldMaxValue)
        values(:) = minValue;
        return;
    end
    
    %
    %  Perform normalization
    values = ((values - oldMinValue) ./ (oldMaxValue - oldMinValue) .* (maxValue - minValue)) + minValue;
end



function [flattenedImage, pixelShift, invalidIndices, flattenedLine] = ...
    flattenImage(image, lineToFlatten, referencePoint, fillType)

    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if nargin < 3 || isempty(referencePoint)
        referencePoint = 0;
    end
    
    if nargin < 4 || isempty(fillType)
        fillType = 1;
    end
    
    if isempty(image)
        error('Image cannot be empty');
    end
    
    if size(image, 3) > 1
        error('Image must be 1D or 2D');
    end
    
    if size(image, 2) == 1
        image = image';
    end
    
    % Do nothing if the line to flatten was not specified
    if isempty(lineToFlatten) || all(isnan(lineToFlatten))
        flattenedImage = image;
        pixelShift = [];
        invalidIndices = [];
        flattenedLine = [];
        return;
    end

    image = double(image);
    imageSize = size(image);
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
    if length(lineToFlatten) ~= imageWidth
        error('The line to flatten must be the same length as the image width');
    end
    
    % We are doing pixel shifts, so sub-pixel shifting is not allowed
    lineToFlatten = round(lineToFlatten);

    
    %----------------------------------------------------------------------
    %  Flatten a vector (line)
    %----------------------------------------------------------------------
    
    if imageHeight == 1
        
        % Shift the line so that it lies on the reference location        
        switch referencePoint
            case 0
                reference = round(mean(lineToFlatten));
            case 1
                reference = round(min(lineToFlatten));
            case 2
                reference = round(max(lineToFlatten));
            otherwise
                reference = round(mean(lineToFlatten) + referencePoint);
        end
        pixelShift = reference - lineToFlatten;
        flattenedImage = image + pixelShift;
        invalidIndices = [];
        flattenedLine =  reference * ones(imageSize);
        
        
    %----------------------------------------------------------------------
    %  Flatten a matrix (Bscan) based on the PRNL
    %----------------------------------------------------------------------
    
    else
        % Shift the image with respect reference location
        flattenedImage = image;
        switch referencePoint
            case 0
                reference = round(mean(lineToFlatten));
            case 1
                reference = round(min(lineToFlatten));
            case 2
                reference = round(max(lineToFlatten));
            otherwise
                reference = round(mean(lineToFlatten) + referencePoint);
        end
        pixelShift = reference - lineToFlatten;
        invalidImage = zeros(imageSize);

        % Loop through each column, shifting it up or down so that the line
        % to flatten lies on a flat line. 
        for index = 1:length(lineToFlatten)
            
            % Circular shift the column
            flattenedImage(1:end, index) = ...
                circshift(image(1:end, index), pixelShift(index));

            % If shifted down
            if pixelShift(index) > 0
                
                % Get locations where pixels have been shifted out. These
                % are invalid regions that need to be replaced
                invalidIndices = 1:pixelShift(index);
                
                switch(fillType)
                    case 0
                        mirroredColumn = zeros(1,length(invalidIndices));
                    case 1
                        % Get the mirror image of the valid regions of the column,
                        % which will be used to replace the invalid regions
                        mirroredColumn = padarray( ...
                            flattenedImage(invalidIndices(end) + 1:end,index), ...
                            length(invalidIndices), ...
                            'symmetric', 'pre');

                        mirroredColumn = mirroredColumn(1:length(invalidIndices));
                    case 2
                        mirroredColumn = (invalidIndices(end)+1)*ones(1,length(invalidIndices));
                end
                
            % If shifted up
            elseif pixelShift(index) < 0
                
                % Get locations where pixels have been shifted out. These
                % are invalid regions that need to be replaced
                invalidIndices = imageHeight + pixelShift(index) + 1:imageHeight;
               
                switch(fillType)
                    case 0
                        mirroredColumn = zeros(1,length(invalidIndices));
                        
                    case 1
                        % Get the mirror image of the valid regions of the column,
                        % which will be used to replace the invalid regions
                        mirroredColumn = padarray( ...
                            flattenedImage(1:invalidIndices(1)-1,index), ...
                            length(invalidIndices), ...
                            'symmetric', 'post');

                        mirroredColumn = mirroredColumn(end-length(invalidIndices)+1:end);
                        
                    case 2
                        mirroredColumn = (invalidIndices(1)-1)*ones(1,length(invalidIndices));
                end
                
            % If no shifting
            else
                invalidIndices = [];
                mirroredColumn = [];
            end
            
            % Replace the invalid indices with the mirror image of the
            % valid pixels. This is so that there is no artificial gradient
            % created when segmenting the image later on.
            flattenedImage(invalidIndices, index) = mirroredColumn;
            
            % Keep track of which indices on the image are invlaid 
            invalidImage(invalidIndices, index) = 1;
        end

        % Get the indices of the invalid regions
        invalidIndices = find(invalidImage == 1);
        
        % Get the resulting flattened line
        flattenedLine = lineToFlatten + pixelShift;
    end
end


function layers = resampleLayers(layers, originalSize, newSize)

    %----------------------------------------------------------------------
    % Calculate the scaling needed to upsample
    %----------------------------------------------------------------------

    scale = newSize ./ originalSize;
    
    if all(scale == 1)
        return;
    end
    
    
    %----------------------------------------------------------------------
    %  Upsample in the y direction
    %----------------------------------------------------------------------

    layers = round(layers * scale(1));
    
    
    %----------------------------------------------------------------------
    % Upsample each layer in the x direction using interpolation
    %----------------------------------------------------------------------
    
    nLayers = size(layers,1);
        
    if scale(2) < 1
        x = round(linspace(1,originalSize(2),newSize(2)));
        y = layers;
        layers = y(:,x);
    else
        newWidth = newSize(2);
        x = 1:newWidth;
        y = layers;

        layers = NaN(nLayers,newWidth);
        ind = round(1:scale(2):newWidth);
        ind(end) = newWidth;
        layers(:,ind) = y;

        % Loop through each layer
        for iLayer = 1:nLayers
            layer = layers(iLayer,:);

            if all(isnan(layer))
                continue;
            end

            % Interpolate
            validInd = ~isnan(layer);
            invalidInd = ~validInd;
            invalidInd(1:find(invalidInd == 0,1,'first')) = 0;
            invalidInd(find(invalidInd == 0,1,'last'):end) = 0;

            layer(invalidInd) = round(interp1( ...
                x(validInd),layer(validInd),x(invalidInd),'spline'));

            % Make sure the layers do not cross
            if iLayer > 1
                invInd = (layer - layers(iLayer-1,:)) < 0;
                layer(invInd) = layers(iLayer-1,invInd);
            end

            layers(iLayer,:) = layer;
        end
    end
end

function image = blurImage(image, filter)

    image = double(image);
    [imageHeight, imageWidth] = size(image);

    filter = double(filter);
    filterSize = size(filter);

    %
    %  Calculate the half widths of the filter
    %
    yHW = round(filterSize(1)/2);
    xHW = round(filterSize(2)/2);

    %
    %  Pad the borders to minimize convolution edge effects
    %
    image(1+yHW:imageHeight+yHW,1+xHW:imageWidth+xHW) = image;

    image(imageHeight+yHW+1:imageHeight+2*yHW,:) = image(imageHeight+yHW-1:-1:imageHeight,:);
    image(yHW:-1:1,:) = image(1+yHW+1:2*yHW+1,:);

    image(:,imageWidth+xHW+1:imageWidth+2*xHW) = image(:,imageWidth+xHW-1:-1:imageWidth);
    image(:,xHW:-1:1) = image(:,1+xHW+1:2*xHW+1);

    %
    %  Blur the image
    %
    image = conv2(image,filter,'same');

    %
    %  Restore the original size of the image
    %
    image = image(1+yHW:imageHeight+yHW,1+xHW:imageWidth+xHW);
end

function weightingMatrices = normal_weightingMatrix( ...f
    image, ...
    axialRes, ...
    lateralRes, ...
    matrixIndices, ...
    invalidIndices, ...
    params)

    %----------------------------------------------------------------------
    %  Initialize missing input parameters
    %----------------------------------------------------------------------
    
    if nargin < 4
        matrixIndices = [];
    end
    
    if nargin < 5
        invalidIndices = [];
    end
    
    if nargin < 6
        params = [];
    end

    %----------------------------------------------------------------------
    %  Validate input parameters
    %----------------------------------------------------------------------
    
    if isempty(matrixIndices)
        params = normal_getParameters();
        params = params.graphCutParams.weightingMatrixParams;
    end
    
    nMatrices = length(params.WEIGHT_RANGES);
    
    if ~isempty(matrixIndices)
        matrixIndices(matrixIndices < 1) = [];
        matrixIndices(matrixIndices > nMatrices) = [];
    end
    if isempty(matrixIndices)
        matrixIndices = 1:nMatrices;
    end   

    imageSize = size(image);
    imageHeight = imageSize(1);
    imageWidth = imageSize(2);
    
   
    %----------------------------------------------------------------------
    %  Separate the vertical weights on the leftmost and rightmost columns
    %  from the rest of the weights in the image
    %----------------------------------------------------------------------
    
    %
    %  Generate edges and points on the image
    %
    edges = createLattice(imageSize);
    
    %
    %  Separate the edges corresponding to the left-most and right-most
    %  columns that were artificially added to the image
    %
    maxIndex = imageHeight*imageWidth;
    leftColIndices = 1:(imageHeight-1);
    rightColIndices = ((imageHeight-1)*(imageWidth-1) + 1) : maxIndex;
    
    columnIndices = [leftColIndices, rightColIndices];
    imageIndices = setdiff(1:size(edges,1), columnIndices);
    
    columnEdges = edges(columnIndices,:);
    imageEdges = edges(imageIndices,:);
    
    
    %----------------------------------------------------------------------
    %  Calculate the weights based on the image gradient.  Lower weights 
    %  are assigned to areas with a higher gradient
    %----------------------------------------------------------------------
    
    %
    %  Filter the image
    %
    xFilterSize = round(params.X_FILTER_SIZE / lateralRes);
    yFilterSize = round(params.Y_FILTER_SIZE / axialRes); 
    filter = fspecial('gaussian',[yFilterSize,xFilterSize],params.SIGMA);   
    smoothImage = blurImage(image,filter);
    filter = fspecial('gaussian',[1,xFilterSize],params.SIGMA);   
    smoothImage2 = blurImage(image,filter);
    
    %
    %  Create two edge maps (one for edges that transition from dark->light
    %  in the vertical direction, and one for edges transitioning from 
    %  light->dark).  
    %
    lightDarkEdgeImage = (blurImage(smoothImage, -params.EDGE_FILTER) > 0) .* ...
                         blurImage(smoothImage, -params.EDGE_FILTER);
                     
    darkLightEdgeImage = (blurImage(smoothImage, params.EDGE_FILTER) > 0) .* ...
                         blurImage(smoothImage, params.EDGE_FILTER);       
    
    lightDarkEdgeImage2 = (blurImage(smoothImage2, -params.EDGE_FILTER) > 0) .* ...
                         blurImage(smoothImage2, -params.EDGE_FILTER);
                     
    darkLightEdgeImage2 = (blurImage(smoothImage2, params.EDGE_FILTER) > 0) .* ...
                         blurImage(smoothImage2, params.EDGE_FILTER);              
    
    % Make it difficult to cross the opposite gradient
    darkLightInd = (darkLightEdgeImage > 0);
    lightDarkInd = (lightDarkEdgeImage > 0);
    darkLightInd2 = (darkLightEdgeImage2 > 0);
    
    darkLightEdgeImage(lightDarkInd) = 0;
    lightDarkEdgeImage(darkLightInd) = 0;
    lightDarkEdgeImage2(darkLightInd2) = 0;
                     
    % Normalize the weights
    lightDarkEdgeImage = normalizeValues(lightDarkEdgeImage,0,1);  
    darkLightEdgeImage = normalizeValues(darkLightEdgeImage,0,1); 
    lightDarkEdgeImage2 = normalizeValues(lightDarkEdgeImage2,0,1);  

    lightDarkEdgeImage3 = lightDarkEdgeImage;
    
    % Only keep the strongest gradient, removing consecutive gradient
    % values
    for iCol = 1:imageWidth
        column = darkLightEdgeImage(:,iCol);
        maxima = find(diff(sign(diff([0;column;0]))) < 0);
        darkLightEdgeImage(:,iCol) = 0;
        darkLightEdgeImage(maxima,iCol) = column(maxima);
        
        column = lightDarkEdgeImage(:,iCol);
        maxima = find(diff(sign(diff([0;column;0]))) < 0);
        lightDarkEdgeImage(:,iCol) = 0;
        lightDarkEdgeImage(maxima,iCol) = column(maxima);
    end
    
    % Make it even more difficult to cross the opposite gradient
    darkLightEdgeImage(lightDarkInd) = -1;
    lightDarkEdgeImage(darkLightInd) = -1;
    
    % Set values in the invalid region to zero
    darkLightEdgeImage(invalidIndices) = 0;
    lightDarkEdgeImage(invalidIndices) = 0;
    lightDarkEdgeImage2(invalidIndices) = 0;
    lightDarkEdgeImage3(invalidIndices) = 0;
    
    %
    %  Calculate the gradient weights for each of the edge maps
    %        
    darkLightGradientWeights = 2 - ...
        darkLightEdgeImage(imageEdges(:,1)) - ...
        darkLightEdgeImage(imageEdges(:,2));
    
    lightDarkGradientWeights = 2 - ...
        lightDarkEdgeImage(imageEdges(:,1)) - ...
        lightDarkEdgeImage(imageEdges(:,2));
    
    lightDarkGradientWeights2 = 2 - ...
        lightDarkEdgeImage2(imageEdges(:,1)) - ...
        lightDarkEdgeImage2(imageEdges(:,2));
    
    lightDarkGradientWeights3 = 2 - ...
        lightDarkEdgeImage3(imageEdges(:,1)) - ...
        lightDarkEdgeImage3(imageEdges(:,2));
    
    
    %----------------------------------------------------------------------
    %  Calculate intensity weights
    %----------------------------------------------------------------------
    
    smoothImage(invalidIndices) = 0;
    
    brightIntensityWeights = - smoothImage(imageEdges(:,1)) ...
                             - smoothImage(imageEdges(:,2));
      
    darkIntensityWeights = smoothImage(imageEdges(:,1)) ...
                         + smoothImage(imageEdges(:,2));
    
                     
    %----------------------------------------------------------------------
    %  Calculate the geometric distances between pairs of points.  Lower
    %  weights go to pixel pairs that are closer together
    %----------------------------------------------------------------------
    
    [yFirstPoint, xFirstPoint] = ind2sub(imageSize, imageEdges(:,1));
    [ySecondPoint, xSecondPoint] = ind2sub(imageSize, imageEdges(:,2));
    
    distanceWeights = sqrt( ...
        (xFirstPoint - xSecondPoint).^2 + (yFirstPoint - ySecondPoint).^2);
     
    
    %----------------------------------------------------------------------
    %  Create a weighting matrix. Rows represent the indices of the first
    %  node and columns represent the indices of the second node.  The 
    %  values represent the weights for the edge that is formed between
    %  the two nodes
    %----------------------------------------------------------------------

    % Define the matrices    
    weights = { ...
        {darkLightGradientWeights}, ...                                          % 01: dark-light
        {brightIntensityWeights,distanceWeights}, ...                            % 02: bright, distance
        {lightDarkGradientWeights2,distanceWeights}, ...                         % 03: light-dark, distance
        {lightDarkGradientWeights3,distanceWeights}, ...                         % 04: light-dark, distance
        {darkLightGradientWeights,distanceWeights}, ...                          % 05: light-dark, dark, short weights
        {lightDarkGradientWeights,darkIntensityWeights,distanceWeights}, ...     % 06: light-dark, dark, distance
    };
    nMatrices = length(weights);   
    
    % Set default weight ranges if necessary 
    if isempty(params.WEIGHT_RANGES)
        params.WEIGHT_RANGES = { ...
            {[0,1]}, ...              % 01: dark-light
            {[0,1],[0,1]}, ...        % 02: bright, distance
            {[0,1],[0,1],[0,1]}, ...  % 03: light-dark, distance
            {[0,1],[0,1]}, ...        % 04: light-dark, distance
            {[0,1],[0,1]}, ...        % 05: light-dark, dark, short weights
            {[0,1],[0,1],[0,1]}, ...  % 06: light-dark, dark, distance
        };
    end
    
    % Remove invalid matrix indices
    matrixIndices(matrixIndices < 1) = [];
    matrixIndices(matrixIndices > nMatrices) = [];
    if isempty(matrixIndices)
        matrixIndices = 1:nMatrices;
    end
    
    % Populate the matrices specified
    weightingMatrices = cell(nMatrices,1);
    matrixSize = maxIndex;
    
    for iMatrix = 1:length(matrixIndices)
        
        matrixIndex = matrixIndices(iMatrix);
        
        weightingMatrices{matrixIndex} = generateWeightingMatrix( ...
            matrixSize, ...
            imageEdges, ...
            weights{matrixIndex}, ...
            params.WEIGHT_RANGES{matrixIndex}, ...
            columnEdges, ...
            MIN_WEIGHT);
    end
end
end