%% This is for playing with data coming from Movebank
% filtered with movebank.sh

disp('reading data...')

f = fopen('./movebank_sample.csv');
data = textscan(f, '%*d %s %f %f %*[^\n\r]', ...
    'CollectOutput', true, 'HeaderLines', 1, 'Delimiter', ',');
fclose(f);

[d.t,idx] = sort(datenum(data{1}));
d.x = data{2}(idx,:);
d.s = 1e-3;

%%
disp('building model...')
tic
p_xt = brownianb(d);
toc

%% USING BUILT-IN INTEGRAL
disp('simulating...')

x = linspace(min(d.x(:,1)), max(d.x(:,1)), 100);
y = linspace(min(d.x(:,2)), max(d.x(:,2)), 100);
% x = linspace(4.1e8, 4.25e8, 100);
% y = linspace(1e7,  max(d.x(:,2)), 100);
[xx,yy] = meshgrid(x,y);

% expected occupation time is the integral over time of
% the probability of the particle being in a given point

tlim = [min(d.t), min(d.t)+700];
p = arrayfun(@(x,y) integral(@(t)p_xt([x,y],t), tlim(1), tlim(2)), xx, yy);
p = p / max(p(:));

%% USING CUSTOM INTEGRATION METHOD
disp('simulating...')

x = linspace(min(d.x(:,1)), max(d.x(:,1)), 20);
y = linspace(min(d.x(:,2)), max(d.x(:,2)), 20);
% x = linspace(4.1e8, 4.25e8, 100);
% y = linspace(1e7,  max(d.x(:,2)), 100);
[xx,yy] = meshgrid(x,y);
n = length(x);
p = zeros(size(xx));

% expected occupation time is the integral over time of
% the probability of the particle being in a given point

dt = 24 / 24; 
tt = min(d.t):dt:min(d.t)+700;
% tt = linspace(min(d.t), max(d.t), 1e4);
dt = tt(2) - tt(1);
% dp = zeros(size(p,1), size(p,2), length(tt));
for i = 1:length(tt)
    t = tt(i);
    % dp(:,:,i) = dt * arrayfun(@(x,y) p_xt([x,y], t), xx, yy);
    dp = dt * arrayfun(@(x,y) p_xt([x,y], t), xx, yy);
    % imagesc(x,y,dp(:,:,i));
    % drawnow
    % p = p + dp(:,:,i);
    p = p + dp;
end
p = p / max(p(:)); % normalized occupation time

%%
disp('plot')

figure;
imagesc(y,x,log(p'));
axis xy
drawnow

%%

pp = log(p(:));
m = [xx(:), yy(:), (pp+100)/100];
m = m(pp > -100, :);
dlmwrite('web/data/demo2.csv', m, 'precision', '%.8f');

%%

out.p = p;
out.x = x;
out.y = x;
saveData('movebank.mat', out);
