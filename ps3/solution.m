%% Solution 1a
clear, clc;

pts3d_norm = dlmread('../../input/ps3/pts3d-norm.txt');
pts2d_norm_pic = dlmread('../../input/ps3/pts2d-norm-pic_a.txt');

A = zeros(40, 12);

for i = 1:20
    A((i-1)*2+1, 1:3) = pts3d_norm(i, :);
    A((i-1)*2+1, 4) = 1;
    A((i-1)*2+1, 9:11) = -pts2d_norm_pic(i, 1) * pts3d_norm(i, :);
    A((i-1)*2+1, 12) = -pts2d_norm_pic(i, 1);
    
    A(i*2, 5:7) = pts3d_norm(i, :);
    A(i*2, 8) = 1;
    A(i*2, 9:11) = -pts2d_norm_pic(i, 2) * pts3d_norm(i, :);
    A(i*2, 12) = -pts2d_norm_pic(i, 2);
end

[~,~,V] = svd(A' * A);

M = V(:, end);
M = reshape(M, 4, 3)';
disp('M is');
disp(M);    % scaled equivalent to the matrix in ps3

last_pt = M * [pts3d_norm(end, :), 1]';
last_pt = last_pt/last_pt(3);   % go to inhomogeneous
disp(['Projected point: ', num2str(last_pt(1:2)')]);
% disp(['Ground truth: ', num2str(pts2d_norm_pic(end, :))]);
disp(['Residual: ', num2str(norm(last_pt(1:2)' - pts2d_norm_pic(end, :)))]);

%% Solution 1b
clear;

pts3d = dlmread('../../input/ps3/pts3d.txt');
pts2d_pic_b = dlmread('../../input/ps3/pts2d-pic_b.txt');

all_norms = [];
num_rep = 10;
min_norm = realmax;

for k = [8, 12, 16]
    mean_norms = zeros(num_rep, 1);
    
    for n_rep = 1 : num_rep
        randnums = randperm(20);
        A = zeros(k*2, 12);

        % use k points to construct the matrix
        for i = 1 : k
            j = randnums(i);
            A((i-1)*2+1, 1:3) = pts3d(j, :);
            A((i-1)*2+1, 4) = 1;
            A((i-1)*2+1, 9:11) = -pts2d_pic_b(j, 1) * pts3d(j, :);
            A((i-1)*2+1, 12) = -pts2d_pic_b(j, 1);

            A(i*2, 5:7) = pts3d(j, :);
            A(i*2, 8) = 1;
            A(i*2, 9:11) = -pts2d_pic_b(j, 2) * pts3d(j, :);
            A(i*2, 12) = -pts2d_pic_b(j, 2);
        end

        [~,~,V] = svd(A' * A);

        M = V(:, end);
        M = reshape(M, 4, 3)';

        % 4 unused points
        last_pts = M * [pts3d(end-3 : end, :), ones(4, 1)]';
        last_pts = last_pts./repmat(last_pts(3, :), 3, 1);   % go to inhomogeneous

        % calculate residuals
        norms = zeros(1, 4);
        for i = 1 : 4
            norms(i) = norm(last_pts(1:2, i) - pts2d_pic_b(20-(4-i), :)');
        end
        mean_norm = mean(norms);

        mean_norms(n_rep) = mean_norm;
        
        % record the best M
        if mean_norm < min_norm
            min_norm = mean_norm;
            best_M = M;
        end
    end
    all_norms = cat(2, all_norms, mean_norms);
end

disp('Average residual:');
disp(all_norms);
disp('Best M:');
disp(best_M);

%% Solution 1c
Q = best_M(:, 1:3);
C = -inv(Q) * best_M(:, 4);

disp(['Location of the camera: ', num2str(C')]);

%% Solution 2a
clear;

pts2d_pic_a = dlmread('../../input/ps3/pts2d-pic_a.txt');
pts2d_pic_b = dlmread('../../input/ps3/pts2d-pic_b.txt');

A = zeros(20, 9);

for i = 1:20
    u = pts2d_pic_a(i, 1);
    v = pts2d_pic_a(i, 2);
    up = pts2d_pic_b(i, 1);
    vp = pts2d_pic_b(i, 2);
    A(i, :) = [up*u, up*v, up, vp*u, vp*v, vp, u, v, 1];
end

[~,~,V] = svd(A'*A);

F = V(:, end);
F = reshape(F, 3, 3)';

format long
disp('Estimated F:');
disp(F);
format short

%% Solution 2b
[U,S,V] = svd(F);
S(3,3) = 0;
F = U * S * V';

format long
disp('Estimated F of rank 2:');
disp(F);
format short

%% Solution 2c
im_a = imread('../../input/ps3/pic_a.jpg');
im_b = imread('../../input/ps3/pic_b.jpg');

% two sides of the image
l_L = cross([1, 1, 1], [1, size(im_a, 1), 1]);
l_R = cross([size(im_a, 2), 1, 1], [size(im_a, 2), size(im_a, 1), 1]);

figure(1);
imshow(im_a);
hold on;
for i = 1 : 20
    epi_line = F' * [pts2d_pic_b(i, 1); pts2d_pic_b(i, 2); 1];

    i_L = cross(epi_line, l_L);
    i_L = i_L/i_L(3);
    i_R = cross(epi_line, l_R);
    i_R = i_R/i_R(3);

    line([i_L(1), i_R(1)], [i_L(2), i_R(2)], 'Color','black');
end
hold off;

figure(2);
imshow(im_b);
hold on;
for i = 1 : 20
    epi_line = F * [pts2d_pic_a(i, 1); pts2d_pic_a(i, 2); 1];   % changed

    i_L = cross(epi_line, l_L);
    i_L = i_L/i_L(3);
    i_R = cross(epi_line, l_R);
    i_R = i_R/i_R(3);

    line([i_L(1), i_R(1)], [i_L(2), i_R(2)], 'Color','black');
end
hold off;

%% Solution 2d
mean_a = mean(pts2d_pic_a);
std_a = std(pts2d_pic_a);

T_a = diag([1./std_a 1]) * [1 0 -mean_a(1); 0 1 -mean_a(2); 0 0 1];

mean_b = mean(pts2d_pic_b);
std_b = std(pts2d_pic_b);

T_b = diag([1./std_b 1]) * [1 0 -mean_b(1); 0 1 -mean_b(2); 0 0 1];

% pts2d_pic_a_norm = T_a * cat(1, pts2d_pic_a', ones(1, 20));
% pts2d_pic_a_norm = pts2d_pic_a_norm ./ repmat(pts2d_pic_a_norm(3, :), 3, 1);
% pts2d_pic_a_norm = pts2d_pic_a_norm(1:2, :)';
% pts2d_pic_b_norm = T_b * cat(1, pts2d_pic_b', ones(1, 20));
% pts2d_pic_b_norm = pts2d_pic_b_norm ./ repmat(pts2d_pic_b_norm(3, :), 3, 1);
% pts2d_pic_b_norm = pts2d_pic_b_norm(1:2, :)';

% since T_a(3, 3) == T_b(3, 3) == 1, the above block can be simplified as:
pts2d_pic_a_norm = T_a(1:2, :) * cat(1, pts2d_pic_a', ones(1, 20));
pts2d_pic_a_norm = pts2d_pic_a_norm';
pts2d_pic_b_norm = T_b(1:2, :) * cat(1, pts2d_pic_b', ones(1, 20));
pts2d_pic_b_norm = pts2d_pic_b_norm';

% calculating F is the same as 2a and 2b
A = zeros(20, 9);

for i = 1:20
    u = pts2d_pic_a_norm(i, 1);
    v = pts2d_pic_a_norm(i, 2);
    up = pts2d_pic_b_norm(i, 1);
    vp = pts2d_pic_b_norm(i, 2);
    A(i, :) = [up*u, up*v, up, vp*u, vp*v, vp, u, v, 1];
end

[~,~,V] = svd(A'*A);

F_hat = V(:, end);
F_hat = reshape(F_hat, 3, 3)';

% enforce rank of 2
[U,S,V] = svd(F_hat);
S(3,3) = 0;
F_hat = U * S * V';

format long
disp('T_a is:');
disp(T_a);
disp('T_b is:');
disp(T_b);
disp('Normalized F is:');
disp(F_hat);
format short

%% Solution 2e
F_new = T_b' * F_hat * T_a;
format long
disp('New F is:');
disp(F_new);
format short

% calculating epipolar lines is similar to 2c
figure(3);
imshow(im_a);
hold on;
for i = 1 : 20
    % notice that do NOT use normalized points here
    epi_line = F_new' * [pts2d_pic_b(i, 1); pts2d_pic_b(i, 2); 1];

    i_L = cross(epi_line, l_L);
    i_L = i_L / i_L(3);
    i_R = cross(epi_line, l_R);
    i_R = i_R / i_R(3);

    line([i_L(1), i_R(1)], [i_L(2), i_R(2)], 'Color','black');
end
hold off;

figure(4);
imshow(im_b);
hold on;
for i = 1 : 20
    % notice that do NOT use normalized points here
    epi_line = F_new * [pts2d_pic_a(i, 1); pts2d_pic_a(i, 2); 1];

    i_L = cross(epi_line, l_L);
    i_L = i_L / i_L(3);
    i_R = cross(epi_line, l_R);
    i_R = i_R / i_R(3);

    line([i_L(1), i_R(1)], [i_L(2), i_R(2)], 'Color','black');
end
hold off;
