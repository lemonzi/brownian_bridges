%% This is for playing with data coming from Google Takeout
% filtered with parse_location.js

disp('reading data...')

% Modify this with your own data
% Filter it with parse_location.js first
data = dlmread('samples/locations.tsv', '\t');

% Heavy down-sampling to speed up later integration
d.x = data(1:1000:end, 2:3);
d.t = data(1:1000:end, 1);
% d.s = data(1:1000:end, 4) .^ 2; % Too big
d.s = 1e3;

%%
disp('building model...')
tic
p_xt = brownianb(d);
toc

%%
disp('simulating...')

% x = linspace(min(d.x(:,1)), max(d.x(:,1)), 15);
% y = linspace(min(d.x(:,2)), max(d.x(:,2)), 15);

% In my case, this was roughly Catalonia
% x = linspace(4.1e8, 4.25e8, 100);
% y = linspace(1e7,  max(d.x(:,2)), 100);

% This is Barcelona
x = linspace(41.368830e7, 41.423108e7, 100);
y = linspace(2.127428e7,  2.207594e7, 100);

[xx,yy] = meshgrid(x,y);

tlim = [min(d.t), min(d.t)+7*24*3600*1000]; % 1/2 day
p = arrayfun(@(x,y) integral(@(t)p_xt([x,y],t), tlim(1), tlim(2)), xx, yy);
p = p / max(p(:));

%%
disp('plot')

figure;
imagesc(y,x,p');
axis xy
drawnow

%%

% pp = log(p(:));
% m = [xx(:)/1e7, yy(:)/1e7, (pp+100)/100];
% m = m(pp > -100, :);

p = (p-min(p(:))) / (max(p(:))-min(p(:)));
m = [xx(:)/1e7, yy(:)/1e7, p(:)];
m = m(p(:) > 0.3, :);

% Output file, draggable to the web interface
dlmwrite('web/data/demo1.csv', m, 'precision', '%.8f');
