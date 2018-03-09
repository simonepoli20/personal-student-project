%% importing multiples images in a 3D array
%WARNING: images must be resized by a factor ~0.3 otherwise my Matlab
%version cannot compute it

% Specify the folder where the files live.
myFolder = '/Users/simonepoli/Desktop/Corsi_ERASMUS/Personal student project/sample1';

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.jpg');
theFiles = dir(filePattern);

% Loading all the 2D images in a 3D matrix
% I'm working with resized images so I don't have anymore the necessity to resize the images

L = size( imread(theFiles(1).name) ); %Reading the image and getting the resize

array3d = zeros( L(1) , L(2) , length(theFiles) ); %3D matrix of 0s with the right size

for i=1:length(theFiles)
    
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);

    array3d(:,:,i) = rgb2gray( imread(theFiles(i).name) ); 

end

fprintf(1, 'Uploading ended\n');
clear i L;

%%resize the 3D array if needed
%
%k = 0.3;
%array3d = imresize3(array3d,k);

%% image processing

%FINDING NUMBER OF PEEKS ON THE HISTOGRAM

array3d = uint8(array3d);
[val,pos] = imhist(array3d);

for i = 1:50 %exclude the big peek before 50
  val(i) = 0;
end

[posMax,valMax] = findpeaks(val,'MinPeakDist',10,'MinPeakWidth',5);
numThresh = numel(posMax)-1;
% figure(1); plot(pos,val);


% Image processing 
% THRESHOLDING

% in this case I use the function imbinarize() to binarize a 3D volume with
% a given threshold. 

thresh = multithresh(array3d,numThresh);
thresh=double(thresh);

T = thresh(numThresh)/255;

BW = imbinarize(array3d,T);

%OPENING
[x,y,z] = ndgrid(-5:5);

r1 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWopening = imerode(BW,sphere);          %erosion with sphere radius = r1

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWopening = imdilate(BWopening,sphere);  %dilation with sphere radius = r1

% figure(1); imshowpair(BW(:,:,4),BWopening(:,:,4), 'montage'); title('opening effect'); %show the effect of opening to a single image.

% to visualize on the entire volume use: 
% volumeViewer(BWopening)

% CLOSING
[x,y,z] = ndgrid(-5:5);

r1 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWclosing = imdilate(BWopening,sphere);  %dilation with sphere radius = r1

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWclosing = imerode(BWclosing,sphere);          %erosion with sphere radius = r1

%figure(2); imshowpair(BW(:,:,4),BWclosing(:,:,4), 'montage'); title('opening effect'); 
%show the effect of closing after opening in a single image.

% to visualize on the entire volume use: 
% volumeViewer(BWclosing)

%CONNECTED ELEMENTS
% using bwconncomp() regionprops() to chose the connected element with more pixels

cc = bwconncomp(BWclosing);  
stats = regionprops3(cc, 'Volume'); 
idx = find( [stats.Volume] == max([stats.Volume]) ); 
BWconn = ismember(labelmatrix(cc), idx);  

% to visualize on the entire volume use: 
% volumeViewer(BW2)

%% subplot
%figure(1);
%subplot(2,3,1); imshow(array3d(:,:,5)); title('original_image');
%subplot(2,3,2); imshow(BW(:,:,5)); title('binarization with automatic threshold');
%subplot(2,3,3); imshow(BWopening(:,:,5)); title('opening');
%subplot(2,3,4); imshow(BWconn(:,:,5)); title('bwconncomp() and regionprops()');
%subplot(2,3,5); imshow(BWclosing(:,:,5)); title('closing');

BW = BWconn;

%figure(2);
%imshowpair(array3d(:,:,5),BW(:,:,5))

clearvars -except array3d BW;

%% second part of the project (after meeting 5/03)

% If I overlap the original image with the binary image (BW) obtained before
% (with all the holes filled), I can obtain an image with an uniform
% background (black) and the bone in different shape of gray.


BWmasked = zeros(size(array3d,1),size(array3d,2),size(array3d,3)); %preallocate the size of the masked image

BW1 = smooth3(BW,'gaussian',5); %smooth the 3D image in order to avoid holes in the surface of the image that could compromise
%the use of the function imfill(). 

for i = 1:size(array3d,3)
    
% Load Mask
mask = BW1(:,:,i);

% Fill holes
mask = imfill(mask, 'holes');

% Create masked image.
maskedImage = array3d(:,:,i);
maskedImage(~mask) = 0;

BWmasked(:,:,i) = maskedImage;

end

clearvars -except array3d BW BWmasked;

%% apply thresholding tecnique again on the masked image to separate bone marrow to che cortical+trabecular bone
%now the number of threshold is 2 (fixed): I have to separate: -bone marrow
%                                                      -cortical+trabecular
%                                                      -background


BWmasked = uint8(BWmasked);

thresh = multithresh(BWmasked,2);
%thresh=double(thresh);

T1 = thresh(1);
T2 = thresh(2);

BWbm = BWmasked >= T1 & BWmasked <= T2;
clearvars -except array3d BW BWmasked BWbm;


%% opening and closing for the bonemarrow
% OPENING
[x,y,z] = ndgrid(-5:5);
r1 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWopening = imerode(BWbm,sphere);          %erosion with sphere radius = r1

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWopening = imdilate(BWopening,sphere);          %dilation with sphere radius = r1

% CLOSING
[x,y,z] = ndgrid(-5:5);
r1 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWclosing = imdilate(BWopening,sphere);         %dilation with sphere radius = r1

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWclosing = imerode(BWclosing,sphere);          %erosion with sphere radius = r1

BWbm = BWclosing;

clearvars -except array3d BW BWmasked BWbm;

% %% 
% figure(1);
% subplot(2,2,1); imshow(array3d(:,:,5)); title('1 original image');
% subplot(2,2,2); imshow(BW(:,:,5)); title('2 binary image of the cortical+trabecular bone');
% subplot(2,2,3); imshow(BWmasked(:,:,5)); title('3 original image after the overlapping with the mask');
% subplot(2,2,4); imshow(BWbm(:,:,5)); title('4 binary image of the bone marrow');
% figure(2);
% imshowpair(BW(:,:,5),BWbm(:,:,5));
%% I need to apply some tecniques to eliminate trabecular bone from the cortical bone

%alpha shape

[x,y,z] = ind2sub(size(BWbm),find(BWbm == 1));
shp = alphaShape(x,y,z,inf);
figure(2); plot(shp);








