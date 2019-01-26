function F = getFundamentalMatrix(l_img,r_img)
x1 = l_img(:,1);
y1 = l_img(:,2);
x2 = r_img(:,1);
y2 = r_img(:,2);
A = zeros(size(l_img,1),9); 
for ind=1:size(l_img,1)
    A(ind,:) = [ (x1(ind))*(x2(ind)), (x1(ind))*(y2(ind)), x1(ind), (y1(ind))*(x2(ind)),...
        (y1(ind))*(y2(ind)), y1(ind), x2(ind), y2(ind),1];
end
[~,~,V] = svd(A);
f1 = reshape(V(:,end), [3,3])';
f1 = f1/f1(3,3);

[U,S,V] = svd(f1);
S(3,3) = 0;
F = U*S*V';
end