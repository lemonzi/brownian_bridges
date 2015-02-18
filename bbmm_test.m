%% This is for playing with data coming from Google Takeout
% filtered with parse_location.js

disp('reading data...')

data = dlmread('locations.tsv', '\t');
d.x = data(:, 2:3);
d.t = data(:, 1);
d.s = data(:, 4) .^ 2;

%%
disp('building model...')
tic
p_xt = brownianb(d);
toc

%%
disp('simulating...')

% x = linspace(min(d.x(:,1)), max(d.x(:,1)), 15);
% y = linspace(min(d.x(:,2)), max(d.x(:,2)), 15);
x = linspace(4.1e8, 4.25e8, 100);
y = linspace(1e7,  max(d.x(:,2)), 100);
[xx,yy] = meshgrid(x,y);
n = length(x);
p = zeros(size(xx));

% expected occupation time is the integral over time of
% the probability of the particle being in a given point

dt = 12 * 3600 * 1000; % n-hour increments (2 years of data)
tt = min(d.t):dt:max(d.t);
% tt = linspace(min(d.t), max(d.t), 1e3);
% dt = tt(2) - tt(1);
dp = zeros(size(p,1), size(p,2), length(tt));
for i = 1:length(tt)
    t = tt(i);
    dp(:,:,i) = dt * arrayfun(@(x,y) p_xt([x,y], t), xx, yy);
    %imagesc(x,y,dp);
    %drawnow
    p = p + dp(:,:,i);
end
p = p / max(p(:)); % normalized occupation time

%%
disp('plot')

figure;
imagesc(y,x,p');
axis xy
drawnow

%%

pp = log(p(:));
m = [xx(:)/1e7, yy(:)/1e7, (pp+100)/100];
m = m(pp > -100, :);
dlmwrite('web/data/demo1.csv', m, 'precision', '%.8f');

%%

out.p = p;
out.x = x;
out.y = x;
saveData('google_bug_allx.mat', out);
