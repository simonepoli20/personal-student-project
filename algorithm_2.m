%% importing multiples images in an array

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

%% Load all the 2D images in a 3D matrix

%no more need to resize the images

for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  
  % such as reading it in as an image array with imread()
  I = imread(fullFileName);
  Ibw = I(:,:,1);
  
 
  % I use the function cat() to create the 3D matrix(:,:,numberImage)
  if k == 1
    array3d = Ibw; 
  else
    array3d = cat(length(theFiles), array3d, Ibw);
  end
  
end

fprintf(1, 'Uploading ended\n');

array3d = squeeze(array3d);


%% finding the peeks of the histogram to decide the number of threshold needed

[val,pos] = imhist(array3d);

for i = 1:50 %exclude the big peek before 50
  val(i) = 0;
end
[posMax,valMax] = peakfinder(val); 
numThresh = numel(posMax)-1;


%% Image processing (option 1) 
% in this case I use the function imbinarize() to binarize a 3D volume with
% a given threshold. Don't allow a double thresholding.

thresh = multithresh(array3d,numThresh);
thresh=double(thresh);

T = thresh(numThresh)/255;

BW = imbinarize(array3d,T);


%% Image processing (option 2) 
% the array is computed without dividing it slice by slice.

%Thresholding hysteresis: the following Matlab function contains the source
%code used for hysteresis thresholding for 3d images (or 2d).
%This hysteresis function performs a dual thresholding operation on a grayscale image (2D or 3D) 
%using two threshold values (lower and upper).

thresh = multithresh(array3d,numThresh);
thresh=double(thresh);

T2 = thresh(2)/255;
T3 = thresh(3)/255;

[tri,~]=hysteresis3d(array3d,T2,T3,26);

Tr = multithresh(tri,2);
Tr=double(Tr);

BW = tri >= Tr(2);


% BW_old = array3d >= T; %doesn't work. The image visualized with
% volumeViewer() is not the one expected. Don't know why
%comparing the 2 images obtained in the attempt 1 and 2, only a few voxels
%changed. I'm not able to note important differences.


%% image processing (option 3) 
% in this case the 3D array is computed slice by slice. It's the old
% version of the algorithm

% cycle on all the images of the 3D matrix
for k = 1 : size(array3d,3)
    
    %ally to every .jpeg in the array3d the function imadjust() to improve contrast
    array3d(:,:,k) = imadjust(array3d(:,:,k));
    
    % Threshold 
    thresh = multithresh(array3d(:,:,k),3);
    BWsingle = array3d(:,:,k)>=thresh(3);
    
    if k == 1
       BW = BWsingle; 
    else
       BW = cat(length(theFiles), BW, BWsingle);
    end
    
    fprintf(1, 'Iteration number %d\n', k);
end

BW = squeeze(BW);

%% opening for 3D image

[x,y,z] = ndgrid(-5:5);

r1 = 2;
r2 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWerosion = imerode(BW,sphere);          %erosion with sphere radius = r1

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r2;
BWopening = imdilate(BWerosion,sphere);  %dilation with sphere radius = r2

figure(1); imshowpair(BW(:,:,4),BWopening(:,:,4), 'montage'); title('opening effect'); %show the effect of opening to a single image.

% to visualize on the entire volume use: 
% volumeViewer(BWopening)

%% bwconncomp() regionprops()

cc = bwconncomp(BWopening); 
stats = regionprops3(cc, 'Volume'); 
idx = find( [stats.Volume] == max([stats.Volume]) ); 
BW2 = ismember(labelmatrix(cc), idx);  

% to visualize on the entire volume use: 
% volumeViewer(BW2)


%% closing for 3D image
[x,y,z] = ndgrid(-5:5);

r1 = 3;
r2 = 3;

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r1;
BWdilate = imdilate(BW2,sphere);  %dilation with sphere radius = r2

sphere = sqrt(x.^2 + y.^2 + z.^2) <= r2;
BWclosing = imerode(BWdilate,sphere);          %erosion with sphere radius = r1

figure(2); imshowpair(BWopening(:,:,4),BWclosing(:,:,4), 'montage'); title('opening effect'); 
%show the effect of closing after opening in a single image.

% to visualize on the entire volume use: 
% volumeViewer(BWclosing)

%% subplot

subplot(2,3,1); imshow(array3d(:,:,5)); title('original_image');
subplot(2,3,2); imshow(BW(:,:,5)); title('binarization with automatic threshold');
subplot(2,3,3); imshow(BWopening(:,:,5)); title('opening');
subplot(2,3,4); imshow(BW2(:,:,5)); title('bwconncomp() and regionprops()');
subplot(2,3,5); imshow(BWclosing(:,:,5)); title('closing');



