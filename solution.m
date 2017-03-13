% 1a
clear, clc;

pts3d_norm = dlmread('pts3d-norm.txt');
pts2d_norm_pic = dlmread('pts2d-norm-pic_a.txt');

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

[U,S,V] = svd(A'*A);

m = V(:, end);
m = reshape(m, 4, 3);
m = m';

example = m * [pts3d_norm(end, :), 1]';
example = example/example(3);
disp(example)

% 1c
Q = m(:, 1:3);
C = -inv(Q)*m(:, 4);

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





