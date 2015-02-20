%% This is for playing with data coming from Movebank
% filtered with movebank.sh

disp('reading data...')

% Modify this with your own data
f = fopen('samples/movebank_sample.csv');

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

%%
disp('simulating...')

x = linspace(min(d.x(:,1)), max(d.x(:,1)), 100);
y = linspace(min(d.x(:,2)), max(d.x(:,2)), 100);
[xx,yy] = meshgrid(x,y);

% expected occupation time is the integral over time of
% the probability of the particle being in a given point

tlim = [min(d.t), min(d.t)+700]; % 700 seconds
p = arrayfun(@(x,y) integral(@(t)p_xt([x,y],t), tlim(1), tlim(2)), xx, yy);
p = p / max(p(:));

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

% Output file, draggable to the web interface
dlmwrite('web/data/demo2.csv', m, 'precision', '%.8f');
