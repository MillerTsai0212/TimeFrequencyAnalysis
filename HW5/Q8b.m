clear all
x = double(imread('pepper.jpg'));
tic
y = wavedbc10(x);
toc
L = size(y, 2)/4;
x1L = y(:, 1:L); x1H1 = y(:, L+1:2*L); x1H2 = y(:, 2*L+1:3*L); x1H3 = y(:, 3*L+1:4*L);

z = iwavedbc10(x1L, x1H1, x1H2, x1H3);
figure(3);
colormap(gray(256));
image(z);

function y = iwavedbc10(x1L, x1H1, x1H2, x1H3)
g1 = [0.1601 -0.6038 0.7243 -0.1384 -0.2423 0.0322 0.0776 0.0062 -0.0126 -0.0033];
h1 = [-0.0033 -0.0126 0.0062 0.0776 0.0322 -0.2423 -0.1384 0.7243 -0.6038 0.1601];
h1T = transpose(h1); g1T = transpose(g1);

X1L = upsample(x1L, 2);
X1H1 = upsample(x1H1, 2);
x0 = conv2(X1L, g1T)+conv2(X1H1, h1T);

X1H2 = upsample(x1H2, 2);
X1H3 = upsample(x1H3, 2);
x1 = conv2(X1H2, g1T)+conv2(X1H3, h1T);

x0T = transpose(x0);
X0 = transpose(upsample(x0T, 2));
x1T = transpose(x1);
X1 = transpose(upsample(x1T, 2));
newx = conv2(X0, g1)+conv2(X1, h1); 
for i=1:10
    newx(1, :)=[]; newx(:, 1)=[];
end
L = size(newx, 1);
for i=1:9
    newx(L-8, :)=[]; newx(:, L-8)=[];
end
y = newx;
end
function y = wavedbc10(x);
%Step1: to set g[n] and h[n]
g = [0.0033 -0.0126 -0.0062 0.0776 -0.0322 -0.2423 0.1384 0.7243 0.6038 0.1601];
gT = transpose(g);
h = [0.1601 -0.6038 0.7243 -0.1384 -0.2423 0.0322 0.0776 0.0062 -0.0126 -0.0033];
hT = transpose(h);
%Step2: to compute v1L[m,n] and v1H[m,n]
xg = conv2(g, x);
xgT = transpose(xg);
v1L = transpose(downsample(xgT, 2));
xh=conv2(h, x);
xhT = transpose(xh);
v1H = transpose(downsample(xhT, 2));

%Step 3: to compute x1L[m,n] and x1H1[m,n]
v1Lg = conv2(gT, v1L);
x1L = downsample(v1Lg, 2);
v1Lh = conv2(hT, v1L);
x1H1 = downsample(v1Lh, 2);

%Step 4: to compute x1H2[m,n] and x1H3[m,n]
v1Hg = conv2(gT, v1H);
x1H2 = downsample(v1Hg, 2);
v1Hh = conv2(hT, v1H);
x1H3 = downsample(v1Hh, 2);
y = [x1L, x1H1, x1H2, x1H3];
end