function [H, inliers, res, panaroma_img] = panaroma_stitching(IM1,IM2)

% IM1 = imread('..\data\part1\uttower\left.jpg');
% IM2 = imread('..\data\part1\uttower\right.jpg');
%  IM1 = imread('..\data\part1\hill\1.jpg');
%  IM2 = imread('..\data\part1\hill\2.jpg');
%  IM1 = imread('..\data\part1\ledge\1.jpg');
%  IM2 = imread('..\data\part1\ledge\2.jpg');
%  IM1 = imread('..\data\part1\pier\1.jpg');
%  IM2 = imread('..\data\part1\pier\2.jpg');
leftImage = IM1;
rightImage = IM2;

IM1G = rgb2gray(IM1);
IM1 = im2double(IM1G);
IM2G = rgb2gray(IM2);
IM2 = im2double(IM2G);
sigma = 3;    thresh = 0.05;
radius =3;      disp = 0;
sr = 5;         ploty = 0;
nsize = 5;
%% GETTING OUTPUT FROM HARRIS DETECTOR 
[cim1, r1, c1] = harris(IM1, sigma, thresh, radius, disp); 
[cim2, r2, c2] = harris(IM2, sigma, thresh, radius, disp);
% For Image Left - Excluding corners
n1 = 1;
for i = 1:size(r1)
    if((r1(i) > nsize) && (r1(i) < size(IM1,1)-nsize-1) && (c1(i) > nsize)  && (c1(i) < size(IM1,2)-nsize-1))
        n_r1(n1) = r1(i);
        n_c1(n1) = c1(i);
        n1 = n1+1;
    end
end
r1 = n_r1';c1 = n_c1';
% For Image Right - Excluding Corners
n2 = 1;
for i = 1:size(r2)
    if((r2(i) > nsize) && (r2(i) < size(IM2,1)-nsize-1)  && (c2(i) > nsize) && (c2(i) < size(IM2,2)-nsize-1))
        n_r2(n2) = r2(i);
        n_c2(n2) = c2(i);
        n2 = n2+1;
    end
end
r2 = n_r2';c2 = n_c2';
%% Building Desciptors for Key Points in both Images.
% Descriptors for Image Left im1
for i=1:length(r1)
    r = r1(i);
    c = c1(i);
    im_1= IM1(r-nsize:r+nsize,c-nsize:c+nsize);
    im1(i,:) = im_1(:).';
    
end
% Descriptors for Image Right im2
for i=1:length(r2)
    r = r2(i);
    c = c2(i);
    im_2= IM2(r-nsize:r+nsize,c-nsize:c+nsize);
    im_2 = im_2(:).';
    im2(i,:) = im_2(:).';
    
end
%% DIST2
% dist_2 is the pairwise distance between descriptors.
dist_2 = zeros(size(im1,1) * size(im2,1), 3);
index_i = 1;
for i = 1:size(im1,1)
    for j = 1:size(im2,1)
        distance = dist2(im1(i,:), im2(j,:));
        dist_2(index_i, :) = [distance, i, j];
        index_i = index_i + 1;
    end
end
%% FINDING PUTATIVE MATCHES
p_hundred = 100;
dist_2 = sortrows(dist_2,1);
for i=1:p_hundred 
    putative_matches(i,:,:) = dist_2(i,:,:);
end
l_img = zeros(size(putative_matches, 1), 2);
r_img = zeros(size(putative_matches, 1), 2);
for i = 1:size(putative_matches, 1)  
    l_img(i, :) = [c1(putative_matches(i, 2)), r1(putative_matches(i, 2))];
    r_img(i, :) = [c2(putative_matches(i, 3)), r2(putative_matches(i, 3))]; 
end

%% RANSAC
rand= 4;
iterations = 10000;
threshold_distance = 2;
pointsNum= size(l_img,1);
inlier_num= zeros(1,iterations);
homos = cell(1,iterations);
lp = l_img;
rp = r_img;
l_img = l_img';
r_img = r_img';
for i = 1:iterations
	random_numbers = randperm(pointsNum,rand);
	H1 = calculateHomo(l_img(:,random_numbers),r_img(:,random_numbers));
    size1 = size(l_img,2);
    trans_img = H1*[l_img;ones(1,size1)];
    trans_img = trans_img(1:2,:)./repmat(trans_img(3,:),2,1);
    dist = sum((r_img - trans_img).^2);
	inlier_i = find(dist < threshold_distance);
    inlier_num(i) = length(inlier_i);
	if (length(inlier_i) < rand)
        continue;
    end
    homos{i} = H1;
end
%% Calculating Inliers
   [inlier_max,indi] = max(inlier_num);
    H = homos{indi};
    size2 = size(l_img,2);
    trans_img = H*[l_img;ones(1,size2)];
    trans_img = trans_img(1:2,:)./repmat(trans_img(3,:),2,1);
    dist1 = sum((r_img - trans_img).^2,1);
    res = sum(dist1)/size2;
    inliers = find(dist1 < threshold_distance);
%% PLOTTING MATCHES
if (ploty == 1)
    plots = zeros(max(size(leftImage,1), size(rightImage,1)), size(leftImage,2)+size(rightImage,2), size(leftImage,3)); 
    plots(1:size(leftImage,1),1:size(leftImage,2),:) = leftImage;
    plots(1:size(rightImage,1),size(leftImage,2)+1 : size(leftImage,2)+size(rightImage,2),:) = rightImage;
    figure, 
    imshow(uint8(plots));
    hold on;
    N_i= size(inliers', 1);
    for i = 1:N_i
        inlier_i = inliers(i);
        l_inliers(i, :) = lp(inlier_i, :);
        r_inliers(i, :) = rp(inlier_i, :);
    end 

    x1 = l_inliers(:,2);
    y1 = l_inliers(:,1);
    x2 = r_inliers(:,2);
    y2 = r_inliers(:,1);
    for ind = 1:numel(x1)
       plot([y1(ind), y2(ind)+size(leftImage,2)], [x1(ind), x2(ind)],'Color', 'r', 'linewidth', 1 );
    end
end

%% WARPING ONE IMAGE ONTO OTHER
trans = maketform('projective', H');
I = [ 1, 0, 0; 0, 1, 0; 0, 0, 1];
IdentityForm = maketform('affine', I);
[wews,xdata,ydata]=imtransform(leftImage,trans,'bicubic');
xdata1=[min(1,xdata(1)) max(size(rightImage,2),xdata(2))];
ydata1=[min(1,ydata(1)) max(size(rightImage,1),ydata(2))];
left_img_t=imtransform(leftImage,trans,'XData',xdata1,'YData',ydata1);
right_img_t=imtransform(rightImage,IdentityForm,'XData',xdata1,'YData',ydata1);

left_img_t = im2double(left_img_t);
right_img_t = im2double(right_img_t);
for i =1:size(left_img_t,1)
    for j=1:size(left_img_t,2)
        if(left_img_t(i,j) == 0 || right_img_t(i,j) ==0)
            pixelR = max(left_img_t(i,j,1), right_img_t(i,j,1));
            pixelG = max(left_img_t(i,j,2), right_img_t(i,j,2));
            pixelB = max(left_img_t(i,j,3), right_img_t(i,j,3));
            stitched(i,j,1) = pixelR;
            stitched(i,j,2) = pixelG;
            stitched(i,j,3) = pixelB;
        end
        if(left_img_t(i,j) ~= 0 && right_img_t(i,j) ~=0)
            pixelR = (left_img_t(i,j,1) + right_img_t(i,j,1))/2;
            pixelG = (left_img_t(i,j,2) + right_img_t(i,j,2))/2;
            pixelB = (left_img_t(i,j,3) + right_img_t(i,j,3))/2;
            stitched(i,j,1) = pixelR;
            stitched(i,j,2) = pixelG;
            stitched(i,j,3) = pixelB;
        end
    end
end
panaroma_img = stitched;
imshow(stitched);

end
