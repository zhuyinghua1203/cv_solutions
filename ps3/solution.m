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
%%

% 1b
clear, clc;

pts3d = dlmread('pts3d.txt');
pts2d_pic_b = dlmread('pts2d-pic_b.txt');


%%

% 2a
clear, clc;

pts2d_pic_a = dlmread('pts2d-pic_a.txt');
pts2d_pic_b = dlmread('pts2d-pic_b.txt');

F = zeros(20, 9);

for i = 1:20
    u = pts2d_pic_a(i, 1);
    v = pts2d_pic_a(i, 2);
    up = pts2d_pic_b(i, 1);
    vp = pts2d_pic_b(i, 2);
    F(i, :) = [up*u, up*v, up, vp*u, vp*v, vp, u, v, 1];
end

[U,S,V] = svd(F'*F);

f = V(:, end);
f = reshape(f, 3, 3);
f = f';
% disp(rank(f));

% 2b
[U,S,V] = svd(f);
S(3,3) = 0;
f = U*S*V';
% disp(rank(f));

% 2c
F = f;
line = F*[43, 203, 1]';


im_a=mean(imread('pic_a.jpg'), 3);
im_b=mean(imread('pic_b.jpg'), 3);

l_L = cross([1, 1, 1], [1, size(im_a, 1), 1]);
l_R = cross([size(im_a, 2), 1, 1], [size(im_a, 2), size(im_a, 1), 1]);



i_L = cross(line, l_L);
i_L = i_L/i_L(3)
i_R = cross(line, l_R);
i_R = i_R/i_R(3)

(i_R(2)-i_L(2))/(1072-1)*(22-1) + i_L(2) % should be 248





