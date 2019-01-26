function F = fit_fundamental(matches)

N=100;
t_cm = 0;
cm = 0;
threshold = 0.003;
I1 = imread('../data/part2/house1.jpg');
I2 = imread('../data/part2/house2.jpg');
matches = load('../data/part2/house_matches.txt'); 
matches1(:,1) = matches(:,1);
matches1(:,2) = matches(:,2);
matches2(:,1) = matches(:,3);
matches2(:,2) = matches(:,4);
P = size(matches1,1);
for iteration=1:N
    index = randperm(size(matches1,1),8);
    % finding the fundamental matrix for the current iteration.
    t_fit = getFundamentalMatrix(matches2(index,:), matches1(index,:));
    t_cm = 0;
    % checking the number of temporart correct matches based on the
    % threshold value.
    for i = 1:P
        v1 = [matches2(i,:) 1] * t_fit; 
        v2 = v1 *[matches1(i,1); matches1(i,2); 1];
        res = abs(v2);
        if( res <= threshold )
            t_cm = t_cm + 1;
        end
    end
    % selecting the best fit if the no. of matches in the iteration are
    % equal to total no.of matches given.
    if (t_cm == P)
        cm = t_cm;
        b_Fit = t_fit;
        break;
    end
    % setting the best match if current matches are greater than prev
    % matches.
    if (t_cm > cm)
        cm = t_cm;
        b_Fit = t_fit;
    end
end
F = b_Fit;

end