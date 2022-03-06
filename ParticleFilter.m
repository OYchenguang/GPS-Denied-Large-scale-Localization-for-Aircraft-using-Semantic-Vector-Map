clear all
close all


global gpsdata dxy imgx imgright  % 测量数据
global NUM Nps validn sigma sige  % 控制量
global S1 S2 id i idd                 % 误差
global Nmeas SID Nm ha
global score premup estmup estsscore curmup   
global fig
%% match,local odometry, gps
[imgx,score] = readdata();
% load('data.mat');
%% 测量模型

Nps = 15000;          % 粒子总数
validn = 3;          % 有效粒子
sigma = [100,100];   % 匹配误差
sige = 50;           % VO不确定度
fig = 0;             % 1小图，0大图

ha = tight_subplot(3,4,[.05 .025],[.1 .1],[.1 .1]);
initialization();
t=[];
for i = 2:NUM
   tic;
   propagation();
   z = measurement();
   resampling(z);
   t=[t,toc];
   
    % 观察 测量和预测
    axes(ha(i));
%     subplot(3,4,i,'Position',[0, 0, 1/3, 1/4] + [floor(i/4)/NUM, mod(i,4)/NUM, 0, 0 ]);
    k=0.001;
    
    plot(k*curmup(:,1),k*curmup(:,2),'b.');hold on;
    plot(k*estmup(:,1),k*estmup(:,2),'g.');hold on;
    %     plot(k*premup(:,1),k*premup(:,2),'c.');hold on;
    
    plot(k*premup(idd,1),k*premup(idd,2),'r.');hold on;
    plot(k*gpsdata(i,1),k*gpsdata(i,2),'o','Color',[0,0,0],'linewidth',2.5);hold on;
    grid on
    axis equal
    
    if fig==1
        xlim([-1.3 1])
        ylim([-1.3 1])
    else
        xlim([-30 75])
        ylim([-75 30])
    end
%     text(-15,20,['timestep ',num2str(i)],'FontSize',10);


   
%     xlim([-12000 2000])
%     ylim([-2000 1000])
%     xlim([-1500 800])
%     ylim([-1500 1000])
    title(['timestep ',num2str(i)])
%     title(['timestep ',num2str(i),': + prior, o posterior，* measurement'])
end
% legend('measurement','propagation','update','valid particle')

% set(ha,'XTickLabel',''); 
% set(ha,'YTickLabel','')


% hold on; 
% for i=2:NUM
%     sscore = score(i,:);
%     sscore = sscore(sscore~=0);
%     mup = imgx(i,1:2*length(sscore));
%     curmup = reshape(mup,[2,length(sscore)])';
%     plot(curmup(SID(i),1),curmup(SID(i),2),'*');
% end

%% 滤波结果
figure
plot(gpsdata(:,1),gpsdata(:,2),'*');
hold on; grid on;
% SID = SID +mean(gpsdata(2:end,:)) - mean(SID(2:end,:));
plot(SID(2:end,1),SID(2:end,2),'o')
plot(imgright(1:end,1),imgright(1:end,2),'+')
xlim([-1200 700])
ylim([-700 200]);
title('o filter，* gps, + closest measurement')
    
    
matcherr = gpsdata-imgright;
matcherr = sqrt(matcherr(:,1).^2 + matcherr(:,2).^2);

err = gpsdata - SID;
err = sqrt(err(:,1).^2 + err(:,2).^2);

filter = imgright- SID;
filter = sqrt(filter(:,1).^2 + filter(:,2).^2);

figure;
subplot(2,1,1);
plot(err);grid on;
title('filter error')
xlim([validn+1 12]);
ylim([0 200]);
subplot(2,1,2);
plot(matcherr);grid on;
title('match error')
xlim([1 12]);
ylim([0 200]);
% subplot(3,1,3);
% plot(filter);grid on;
% title('error between filter and match results')
% xlim([1 12]);
% ylim([0 50]);


figure;
subplot(2,1,1);
plot(err);grid on;
title('filter error')
xlim([validn+1 12]);
ylim([0 200]);
subplot(2,1,2);
plot(matcherr);grid on;
title('match error')
xlim([1 12]);
ylim([0 200]);
%%
figure;
ha = tight_subplot(4,1,[.05 .025],[.1 .1],[.1 .1]);
axes(ha(1));
plot(err); hold on
plot(err,'*');grid on;

title('position error of filter')
xlim([validn 12]);
ylim([0 200]);

axes(ha(2));
plot(matcherr);hold on
plot(matcherr,'*');grid on;
title('position error of matching')
xlim([validn 12]);
ylim([0 100]);

axes(ha(3));
plot(S2(:,1));hold on
plot(S2(:,1),'*');
title('S1'); grid on
plot(0.99*ones(1,12),'r--');
xlim([validn,12])
ylim([0.95,1])

axes(ha(4));
plot(log(S1/(sigma(1)*2*5)));hold on
plot(log(S1/(sigma(1)*2*5)),'*');
grid on
hold on
plot(zeros(1,12),'r--');
title('S2')
xlabel('timestep')
xlim([validn,12])

%%


figure;
subplot(3,1,1);
plot(err); hold on
plot(err,'*');grid on;

title('position error of filter')

xlim([validn 12]);
ylim([0 200]);

subplot(3,1,2);
plot(S2(:,1));hold on
plot(S2(:,1),'*');
plot(0.99*ones(1,12),'r--');
title('S1'); grid on

xlim([validn,12])
ylim([0.95,1])

subplot(3,1,3);
plot(log(S1/(sigma(1)*2*5)));hold on
plot(log(S1/(sigma(1)*2*5)),'*');
grid on
hold on
plot(zeros(1,12),'r--');
title('S2')
xlabel('timestep')
xlim([validn,12])

%%
figure
subplot(2,1,1);
plot(S2(:,1));hold on
plot(S2(:,1),'*');
title('S1'); grid on
xlabel('timestep')
plot(0.99*ones(1,12),'r--');
xlim([validn,12])
ylim([0.95,1])

subplot(2,1,2);
plot(log(S1/(sigma(1)*2*5)));hold on
plot(log(S1/(sigma(1)*2*5)),'*');
grid on
hold on
plot(zeros(1,12),'r--');
title('S2')
xlabel('timestep')
xlim([validn,12])
ylim([-1,5])

% %% 比较结果
% return
% close all
% for i=1:NUM
%    figure(i);
%    plot(gpsdata(i,1),gpsdata(i,2),'ro'); hold on
%     for j = 1:2:length(imgx(i,:))-1
%         if imgx(i,j)==0 || isnan(imgx(i,j))
%             break;
%         end
%         plot(imgx(i,j),imgx(i,j+1),'b+');hold on
%     end
%     xlim([-1500 800])
%     ylim([-1500 1000]);
%     grid
% end
 
%% 子函数
function Z = MeanValue(estmup, curmup, sigma, presscore, cursscore)
% 如果est在cur+sigma范围内，获得score值，否则score=0？或者衰减？
z = 0;
Z = [];
for i = 1:length(estmup)
    z = presscore(i)/2;
    for j=1:length(curmup)
%         d = abs(estmup(i,:) - curmup(j,:));
%         if(d(1)<sigma(1) && d(2)<sigma(2))
        d = norm(estmup(i,:) - curmup(j,:));
        if(d<sigma(1))
            z = max(z,cursscore(j));
        end
    end
    Z(end+1) = z; 
end
end

function [mup, score, Nmearnow] = resamp2(premup, sscore, Nmeas, sigma, Nps)
% 重新分配粒子
mup = [];
score = [];
Nmearnow = [];
idd = [];
sscore = sscore/sum(sscore); % 归一化
for i = 1:length(sscore)
    npi(i) = round(Nps * sscore(i));
    d = size(premup,2);
    for j=1:npi(i)
        mup(end+1,:) = premup(i,:);% + sigma(1:d).*(rand(1,d)-0.5);
        score(end+1,:) = sscore(i);
        Nmearnow(end+1,:) = Nmeas(i) + 1;
    end
end
score = score/sum(score); % 归一化
end

function [mup, score] = resamp(premup, sscore, sigma, Nps)
% 重新分配粒子
mup = [];
score = [];
sscore = sscore/sum(sscore); % 归一化
for i = 1:length(sscore)
    npi(i) = round(Nps * sscore(i));
    d = size(premup,2);
    for j=1:npi(i)
        mup(end+1,:) = premup(i,:);% + sigma(1:d).*(rand(1,d)-0.5);
        score(end+1,:) = sscore(i);
    end
end
score = score/sum(score); % 归一化
end

function [imgx,score] = readdata()

global NUM gpsdata dxy imgx id
%% 读取gps真值
file_path =  './test/';
img_path_list = dir(strcat(file_path,'*.jpg'));
NUM = 12;

gps = importdata('gps.txt');
for j=2:length(gps.textdata)
    t_gps(j) =  str2num(gps.textdata{j})/10e8;
end

for i=1:NUM
    t = str2num(img_path_list(i).name(1:15));
    dt = abs(t_gps  - t);
    id(i,:) = find(dt == min(dt));
end
for i=1:NUM
    gpsdata(i,:) = [str2num(gps.textdata{id(i),9}),str2num(gps.textdata{id(i),7})];
end

%% 经纬度转米
% [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid)
lon0 = gpsdata(1,1);
lat0 = gpsdata(1,2);

for i=1:NUM
    lat = gpsdata(i,2);
    lon = gpsdata(i,1);
    [x,y,~] = geodetic2ned(lat,lon,0,lat0,lon0,0,wgs84Ellipsoid);
    gpsdata(i,:)  = [y,x];
end

% figure;plot(gpsdata(:,1),gpsdata(:,2),'o'); hold on; plot(gpsdata(1,1),gpsdata(1,2),'*')

%% 相对位移估计值
for i=2:NUM
    dxy(i,:) = gpsdata(i,:) - gpsdata(i-1,:);
    dxy(i,:) = dxy(i,:) + norm(dxy(i,:))*0.05*(randn(1));
end

%% 读取匹配值
file_path =  '.\result\';    % 图像文件夹路径
file_path_list = dir(strcat(file_path,'*.txt'));%获取该文件夹中所
num = length(file_path_list);% 获取总数量

    result = importdata([file_path,file_path_list(1).name]);
    result(isnan(result))=0;
    imgdatac = rshape(result.');
for j=2:num
    result = importdata([file_path,file_path_list(j).name]);
    result(isnan(result))=0;
    imgdatac = [imgdatac, rshape(result.')];
end

% result2 = importdata('result22.txt');
% result3 = importdata('result23.txt');
% result4 = importdata('result24.txt');
% result1 = result1.';
% result2 = result2.';
% result3 = result3.';
% result4 = result4.';
% imgdata1 = rshape(result1);
% imgdata2 = rshape(result2);
% imgdata3 = rshape(result3);
% imgdata4 = rshape(result4);
% imgdatac = [imgdata1,imgdata2,imgdata3, imgdata4];

for i=1:size(imgdatac,1)
%     for j=1:size(imgdatac,2)
        t = imgdatac(i,:);
        id = t~=0;
        imgdata(i,1:sum(id)) = imgdatac(i,id);
%     end
end

%% 分解匹配值
imgx = [];
num = 1;
for i=1:NUM
    num = 1;
    for j = 2:3:length(imgdata(i,:))
        if imgdata(i,j)==0 || isnan(imgdata(i,j))
            break;
        end
        score(i,num) = imgdata(i,j-1);
        
        [x,y,~] = geodetic2ned(imgdata(i,j+1),imgdata(i,j),0,lat0,lon0,0,wgs84Ellipsoid);
        imgx(i,num*2-1:num*2) = [y,x]+[60,30]; % 粗略估计姿态误差
        num = num + 1;
    end
end
end

function imgdata = rshape(result1)
[w,h] = size(result1);
imgdata = [];
row = 0;
col = 1;
for i=1:w*h
    if(result1(i) == row)
        row = row + 1;
        col = 1;
    else
        imgdata(row,col) = result1(i);
        col = col+1;
    end
end
end

function propagation()
global dxy  sige
global premup estmup
global  i
  % t预测, 权值不变
%     estmup = premup + dxy(i,:) + norm(dxy(i,:))*0.01*(randn(size(premup)));
    estmup = premup + dxy(i,:) + sige*(randn(size(premup)));
end

function z = measurement()
global  gpsdata   imgright Nm
global score  estmup imgx curmup
global   sigma  i estsscore id
  % t测量
    sscore = score(i,:);
    sscore = sscore - max(sscore)*0.6;   % 大大加快滤波效果，就很神奇
    id = sscore>0;
    sscore = sscore(id);
    sscore = sscore/sum(sscore);
    cursscore = sscore;
    
    mup = imgx(i,:);
    mup = reshape(mup,[2,length(mup)/2])';
    curmup = mup(id,:);
    
    Nm(i) = length(curmup);
    
    % 计算匹配误差
    imgerr(i) = 10000;
    for j=1:length(curmup)
        d = norm(curmup(j,:) - gpsdata(i,:));
        if d<imgerr(i)
            imgerr(i) = d;
            imgright(i,:) = curmup(j,:);
        end
    end
        
    % 更新w 
    z = MeanValue(estmup, curmup, sigma, estsscore, cursscore);
end

function resampling(z)
global     SID S1 S2
global score premup estmup  Nmeas curmup
global Nps validn sigma  i estsscore id idd
    % 粒子 重采样
    estsscore = z(:)/sum(z);
%     [premup, estsscore] = resamp(estmup, estsscore, sigma, Nps);
    [premup, estsscore, Nmeas] = resamp2(estmup, estsscore, Nmeas, sigma, Nps);
    estsscore = estsscore/sum(estsscore);
   
           
    % 增加测量值粒子
    ssscore = score(i,id);
    for j=1:length(ssscore)
        if(ssscore(j) > 0.7)
            estsscore(end+1) = 2/Nps;
            premup(end+1,:) = curmup(j,:);
            Nmeas(end+1,:) = 1;
        end
    end
    
%     idd = Nmeas>validn;
    idd = Nmeas>=max(validn,max(Nmeas)-1);   
 
    if(i>validn) && sum(idd)>0
        S1(i) = sqrt((max(premup(idd,1)) - min(premup(idd,1)))^2 + (max(premup(idd,2)) - min(premup(idd,2)))^2); 
    else
        S1(i) = sqrt((max(premup(:,1)) - min(premup(:,1)))^2 + (max(premup(:,2)) - min(premup(:,2)))^2); 
    end
%     SID(i,:) = sum(premup(idd,:).*estsscore(idd).*Nmeas(idd))/sum(estsscore.*Nmeas);
    SID(i,:) = sum(premup(idd,:).*estsscore(idd))/sum(estsscore(idd));

    if(sum(idd)>1)
        S2(i,:) = [sum(Nmeas(idd).*estsscore(idd))/sum(estsscore.*Nmeas), sum(estsscore(idd))];
    end
    
end

function initialization()
global  gpsdata    S1  ha
global score premup  imgx Nmeas curmup
global Nps  sigma   estsscore id Nm fig

% t-1测量
sscore = score(1,:);
% sscore = sscore - 0.53;
id = sscore>0;
sscore = sscore(id);
sscore = sscore/sum(sscore);
presscore = sscore';

mup = imgx(1,1:2*length(sscore));
premup = reshape(mup,[2,length(sscore)])';
premup = premup(id,:);

% 重新分配粒子
[premup, estsscore] = resamp(premup, sscore, sigma, Nps);
mup = imgx(1,:);
mup = reshape(mup,[2,length(mup)/2])';
curmup = mup(id,:);
% 粒子连续观测到的次数,应该和estsscore，premup一一对应
Nmeas = ones(size(estsscore));
imgerr(1) = 1000;
for j=1:length(curmup)
    d = norm(curmup(j,:) - gpsdata(1,:));
    if d<imgerr(1)
        imgerr(1) = d;
        imgright(1,:) = curmup(j,:);
    end
end
S1(1) = sqrt((max(premup(:,1)) - min(premup(:,1)))^2 + (max(premup(:,2)) - min(premup(:,2)))^2);
Nm(1) = length(curmup);
    
% figure
% subplot(3,4,1)
axes(ha(1));
k=0.001;
plot(k*premup(:,1),k*premup(:,2),'g.');hold on;
plot(k*curmup(:,1),k*curmup(:,2),'b.');hold on;
% plot(k*gpsdata(1,1),k*gpsdata(1,2),'o','Color',[0,0,0]);hold on;
plot(k*gpsdata(1,1),k*gpsdata(1,2),'o','Color',[0,0,0],'linewidth',2.5);hold on;

grid on
axis equal

if fig==1
    xlim([-1.3 1])
    ylim([-1.3 1])
else
    xlim([-30 75])
    ylim([-75 30])
end

title(['timestep ',num2str(1)])
% text(-15,22,['timestep ',num2str(i)],'FontSize',10);
%     xlim([-12000 2000])
%     ylim([-2000 1000])
end

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);

    for ix = 1:Nw
        ii = ii+1;
%         ha(ii) = axes('Units','normalized', ...
%             'Position',[px py axw axh], ...
%             'XTickLabel','', ...
%             'YTickLabel','');
        ha(ii) = axes('Units','normalized', ...
'Position',[px py axw axh]);
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
end