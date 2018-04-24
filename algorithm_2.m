%% Import 3D image from the file with resized images 

% IMPORTFILE('SampleXX_resized03')
% Imports data from the specified file
% SampleXX_resized03: is the file to read

% Import the file
newData1 = load('-mat','Sample34_resized03');

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

% Image processing

%FINDING AUTOMATICALLY THE NUMBER OF PEEKS ON THE HISTOGRAM

Ires = uint8(Ires);
[val,pos] = imhist(Ires);

for i = 1:5 %exclude the outliers in before 5
    val(i) = 0;
end

[~,id] = max(val);

for i = 1:(id+10) %exclude the big peek before idx+50
    val(i) = 0;
end

[valMax,posMax] = findpeaks(val,'MinPeakDist',10,'MinPeakWidth',5);
numThresh = numel(posMax)-1;

% if this automatic step is not working (for example with images where
% cortical bone and bone marrow have very similar contrast) the operator
% will be able to decide manually the number of threshold he desire.

%figure(1); plot(pos,val);

% BINARIZATION OF THE IMAGE
% Image processing to obtain a binarized image

thresh = multithresh(Ires,numThresh);
thresh = double(thresh);

T2 = thresh(numThresh)/255;

var = 'N';
while (var == 'N')
    
    % THRESHOLDING

    BWcortrab = imbinarize(Ires,T2);
    
    %CONTROL
    
    %showing the image obtained with thresholding technique (the operator can
    %decide if the binary image is good or not.
    figure(1); imshow(BWcortrab(:,:,5)); title('Binary image obtained with thresholding technique');

    % asking to the operator if the binarization has been successful
    fprintf('Has the binarization of the image been successfull? [Y/N]\n');
    var = input('','s');

    %if the binarization had some problems, the operator can insert manually the
    %number of threshold desired
    if (var == 'N')
        
        figure(2); imhist(Ires(:,:,5)); title('Histogram of the original image');
        fprintf('Insert manually the higher threshold T2 to binarize the image: (0<T2<255)\n');
        T2 = input('');
        
        if T2>255 
            T2 = 255;
        end
        
        if T2<0 
            T2 = 0;
        end
        
        T2 = T2/255;  
        
    end
end
close figure 1;

% CLOSING and OPENING techniques
%Also in this case the algorithm first apply opening and closing with a
%default radius r1 and r2. The operator can also decide to change this
%values.

%default values for r1 and r2
r1 = 3; %default value for opening
r2 = 3; %default value for closing

var = 'N';
while (var == 'N')

    % OPENING
    
    SEopening = strel('disk',r1);
    BWopening = imopen(BWcortrab,SEopening); 
    
    % CLOSING
    
    SEclosing = strel('disk',r2);
    BWclosing = imclose(BWopening,SEclosing);       
    
    % CONTROL
    
    % Showing the image obtained with opening and closing technique (the operator can
    % decide if the binary image is good or not)
    figure(2); imshow(BWclosing(:,:,5)); title('binary image obtained with opening and closing technique');

    % asking to the operator if OPENING and CLOSING techniques have been successful
    fprintf('Has the opening and closing technique been successfull? [Y/N]\n');
    var = input('','s');

    %if the opening and closing technique had some problems, the operator can insert manually the
    %values for R1 and R2.
    if (var == 'N')
        
         fprintf('Insert manually the radius (R1) for the OPENING technique :\n');
         r1 = input('');
         fprintf('Insert manually the radius (R2) for the CLOSING technique:\n');
         r2 = input('');
         
    end
end
close figure 2;

% FINDING CONNECTED ELEMENTS
% using bwconncomp() regionprops() to chose the connected element with more
% pixels: this operation make the algorithm more robust.

cc = bwconncomp(BWclosing);  
stats = regionprops3(cc, 'Volume'); 
idx = find( [stats.Volume] == max([stats.Volume]) ); 
BWconn = ismember(labelmatrix(cc), idx);  

% % SUBPLOT
% figure(2);
% subplot(2,3,1); imshow(array3d(:,:,5)); title('original image');
% subplot(2,3,2); imshow(BW(:,:,5)); title('binary image (BW)');
% subplot(2,3,3); imshow(BWopening(:,:,5)); title('BW after opening');
% subplot(2,3,4); imshow(BWclosing(:,:,5)); title('BW after closing');
% subplot(2,3,5); imshow(BWconn(:,:,5)); title('BW after finding connected elements ');
% 
% figure(3);
% imshowpair(array3d(:,:,5),BWconn(:,:,5))

BWcortrab = BWconn; % binary image of cortical and trabecular bone

% FILLING THE INTERNAL VOLUME OF THE BONE

% FIlling the internal volume of the bone is an useful operation for the next steps: 
% in fact if I overlap the original image with the binary image with all the holes filled,
% I can obtain an image with an uniform background (black) and the bone in different 
% shapes of gray.

BWmask = zeros(size(Ires,1),size(Ires,2),size(Ires,3)); 
%preallocating the size of the masked image

for i = 1:size(Ires,3)
    
% Load Mask
mask = BWcortrab(:,:,i);

% Fill holes
mask = imfill(mask, 'holes');

BWmask(:,:,i) = mask;

end

% MORPHOLOGICAL OPERATIONS ON THE BACKGROUND
% To make the algorithm much more robust, we need to do some morphological
% operations also on the background of the image.
% In this way, I can avoid the problems I have when the external boundaries
% of the image have holes or are not continuous.

background = ~BWmask; %inverting the binary image

%OPENING (on the background)

r = 7;
SEopening = strel('disk',r);
BWopening = imopen(background,SEopening); 

%CONNECTED ELEMENTS (on the background)

cc = bwconncomp(BWopening);  
stats = regionprops3(cc, 'Volume'); 
idx = find( [stats.Volume] == max([stats.Volume]) ); 
BWbackground = ismember(labelmatrix(cc), idx);  

BWmasked = ~BWbackground; %inverting again the image to obtain the original one

% overlapping the mask obtained before on the original image
Ibone = zeros(size(Ires,1),size(Ires,2),size(Ires,3)); %preallocate the size of the masked image

for i = 1:size(Ires,3)
    
% Load Mask
mask = BWmasked(:,:,i);

% Create masked image.
maskedImage = Ires(:,:,i);
maskedImage(~mask) = 0;

Ibone(:,:,i) = maskedImage;

end

% THRESHOLDING TECHIQUE TO SEPARATE BONE MARROW FROM CORTICAL AND TRABECULAR BONE.
% applying thresholding tecnique again on the masked image to separate bone
% marrow from the cortical and trabecular bone.
% now the number of threshold is 2 (fixed): I have to separate: -bone marrow
%                                                               -cortical+trabecular
%                                                               -background

% THRESHOLDING

T2 = T2*255;
BWbm = Ibone <= T2;    
BWbonemarrow = BWbm & BWmasked;
   
% IMAGE PROCESSING ON THE BONE MARROW IMAGE 

r1 = 3;
r2 = 3;

% OPENING
  
SEopening = strel('disk',r1);
BWopening = imopen(BWbonemarrow,SEopening); 
    
% CLOSING
    
SEclosing = strel('disk',r2);
BWclosing = imclose(BWopening,SEclosing);       
    
%CONNECTED ELEMENTS
% using bwconncomp() regionprops() to chose the connected element with more pixels

cc = bwconncomp(BWclosing);  
stats = regionprops3(cc, 'Volume'); 
idx = find( [stats.Volume] == max([stats.Volume]) ); 
BWconn = ismember(labelmatrix(cc), idx);  

BWbonemarrow=BWconn;

clearvars -except Ires BWcortrab Icortrab BWbonemarrow BWbackground;

% figure(4);
% subplot(2,2,1); imshow(array3d(:,:,5)); title('Original image');
% subplot(2,2,2); imshow(BW(:,:,5)); title('BW image (cortical and trabecular bone)');
% subplot(2,2,3); imshow(BWfinal(:,:,5)); title('Original overlapped with the mask');
% subplot(2,2,4); imshow(BWbm(:,:,5)); title('BW image of the bone marrow');
% figure(5);
% imshowpair(BW(:,:,5),BWbm(:,:,5));


% TECHNIQUES TO ELIMINATE TRABECULAR BONE FROM THE CORTICAL BONE
% there are mainly 3 ways to compute the convexHull of a 3D volume in
% matlab: -convhulln
%         -convexHull
%         -Alphashape
% Convhull is said to be the most robut and efficient way to do that
% (https://it.mathworks.com/help/matlab/math/computing-the-convex-hull.html)

%first of all I need trasform with the function ind2sub() the 3D
%matrix BWbm into 3 vectors x,y,z.

[x,y,z] = ind2sub(size(BWbonemarrow),find(BWbonemarrow == 1));  % trasform the 3D matrix in 3 vectors
z = double(z);

% % TECNIQUE 1 - Convhull
% 
% K = convhull(x,y,z);
% figure(5); trisurf(K,x,y,z);
% 
% TECNIQUE 2 - ConvexHull

% DT = delaunayTriangulation(x,y,z);
% [Hull,v] = convexHull(DT);
% figure(5); trisurf(Hull,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3));

% I can now compute the visualization of the volume with the function
% trisurf.

% TECNIQUE 3 - AlphaShape

% With the 3 vectors obtained before I can applay alpha-shape that creates a bounding
% volume that envelops the set of 3-D points.

shp = alphaShape(y,x,z,10); %obtaining alphaShape image
figure(6); Ishp = plot(shp); %plotting the image to show what obtained. Ishp is a struct type.

shpPatch = patch2struct(Ishp); %Converting Ishp (patch) into a struct that can be computed by polygon2voxel

VolumeSize = [size(BWcortrab,1) size(BWcortrab,2) size(BWcortrab,3)];
BWbonemarrowAlpha = polygon2voxel(shpPatch,VolumeSize,'none'); %creating a 3D array of the AlphaShape image

% polygon2voxel create the boundaries around alphaShape image: to obtain
% the final volume I need to fill these boundaries
BWbonemarrow = zeros(size(BWcortrab,1),size(BWcortrab,2),size(BWcortrab,3));
for i = 1:size(BWcortrab,3)
    
% Load Mask
mask = BWbonemarrowAlpha(:,:,i);
% Fill holes inside th
BWbonemarrow(:,:,i) = imfill(mask, 'holes');

end

BWbonemarrow = imbinarize(BWbonemarrow);

%overlap the mask created with the original image
Ibonemarrow = zeros(size(Ires,1),size(Ires,2),size(Ires,3)); %preallocate the size of the masked image

for i = 1:size(Ires,3)
    
% Create masked image.
maskedImage = Ires(:,:,i);
maskedImage(~BWbonemarrow(:,:,i)) = 0;

Ibonemarrow(:,:,i) = maskedImage;

end



% Creating the mask to eliminate trabecular bone.
maskCortical = ~imadd(BWbackground,BWbonemarrow);

%applying the mask to the original image
Icortical = zeros(size(Ires,1),size(Ires,2),size(Ires,3)); %preallocate the size of the masked image

for i = 1:size(Ires,3)
    
% Create masked image.
maskedImage = Ires(:,:,i);
maskedImage(~maskCortical(:,:,i)) = 0;

Icortical(:,:,i) = maskedImage;

end

clearvars -except Ires BWcortrab Ibonemarrow Ibone Icortical BWmask BWbonemarrow BWbonemarrowAlpha BWbackground maskCortical;

%% FIND ONLY TRABECULAS
BWmask = imbinarize(BWmask);
mask = imsubtract(BWmask,BWbonemarrow);

r1 = 3;
r2 = 3;

% OPENING
  
SEopening = strel('disk',r1);
BWopening = imopen(mask,SEopening); 
    
% CLOSING
    
SEclosing = strel('disk',r2);
BWclosing = imclose(BWopening,SEclosing);

BWtrabecular = zeros(size(Ires,1),size(Ires,2),size(Ires,3)); 
%preallocate the size of the masked image

for i = 1:size(Ires,3)
    
% Create masked image.
maskedImage = Ires(:,:,i);
maskedImage(~BWclosing(:,:,i)) = 0;

BWtrabecular(:,:,i) = maskedImage;

end
