function H = calculateHomo(l_img,r_img)
l_img1= l_img';
r_img1= r_img';
% r_img1 = [141, 131;480, 159 ;493, 630;64, 601];
% l_img1 = [318, 256;534, 372;316, 670;73, 473];
x1 = l_img1(:,1);
y1 = l_img1(:,2);
x2 = r_img1(:,1);
y2 = r_img1(:,2);
A = zeros(2*size(l_img1, 1), 9);
for ind=1:size(l_img1,1)
    row = 2*ind-1;
    val = [x1(ind), y1(ind), 1];
    A(row,1:3) = val;
    A(row,4:6) = [0,0,0];
    A(row,7:9) = -x2(ind) * val;
    A(row+1,1:3) = [0,0,0];
    A(row+1,4:6) = val;
    A(row+1,7:9) = -y2(ind) * val;
end
[~,~,V] = svd(A);
H = reshape(V(:,end), [3,3])';
H = H/H(3,3);
end