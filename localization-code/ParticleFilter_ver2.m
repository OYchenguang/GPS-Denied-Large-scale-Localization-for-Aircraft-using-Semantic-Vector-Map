clear all
close all

global gpsdata imgdata imgdata_closest
global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1         % 测量更新结果
global u_m sig_m alpha_m u_mi
global x_v sig_v
global X S1 S2
global NUM fig ha Nps validn i imgerr
global gps_id truetraj  Nmea
global mu_registration



looperr = [];

for loop =  1:37
    loop

    %%
    if loop>1
        loopS1(loop-1,:) = S;
        loopN(loop-1) = mean(Nmea);
        loopD(loop-1,:) = [find(S<0,1),12-sum(S<0)];
        looperr(loop-1,:) = [mean(err(S<0)),mean(matcherr), err.'];
        loopX(loop-1) = mu_registration;
    end
    
    
    close all
    clearvars -except looperr  loop loopS1 loopD loopN loopX...
        gpsdata imgdata imgdata_closest  u_f sig_f alpha_f c_f u_f_1 sig_f_1...
        u_m sig_m alpha_m u_mi  x_v sig_v X S1 S2 NUM fig ha Nps validn i imgerr  gps_id truetraj  Nmea mu_registration
    
    
    %% match,local odometry, gps
    readdata();
    Nps = 50000;          % 粒子总数
    validn = 3;           % 有效粒子
    sig_m = [30,30];      % 匹配误差
    sig_v = [50,50];      % VIO不确定度
    fig = 1;              % 1小图�?0大图, -1不画�?
    mu_registration = 0.49+loop/100; % 匹配剔除
%     mu_registration = 0.6
    % ha = tight_subplot(3,4,[.05 .025],[.1 .1],[.1 .1]);
    if(fig~=-1)
        ha = tight_subplot(6,2,[.025 .03],[.04 .04],[.04 .04]);
    end
    Nmea = []; % 有效测量�?
    initialization();
    t=[];
    
    for i = 2:NUM
        tic;
        propagation();
        measurement();
        resampling();
        
        idd = c_f>=validn;
        
        if(i>validn) && sum(idd)>0
            S1(i) = sqrt((max(u_f_1(idd,1)) - min(u_f_1(idd,1)))^2 + (max(u_f_1(idd,2)) - min(u_f_1(idd,2)))^2);
        else
            S1(i) = sqrt((max(u_f_1(:,1)) - min(u_f_1(:,1)))^2 + (max(u_f_1(:,2)) - min(u_f_1(:,2)))^2);
        end
        
        if(sum(idd)>0)
            %         S2(i,:) = [sum(c_f(idd).*alpha_f(idd)./sig_f_1(idd,1))/sum(c_f.*alpha_f./sig_f_1(:,1)), sum(alpha_f(idd))];
            S2(i,:) = [sum(c_f(idd).*alpha_f(idd))/sum(c_f(c_f>1).*alpha_f(c_f>1)), sum(alpha_f(idd))];
            %         S2(i,:) = [sum(alpha_f(idd))/sum(alpha_f(c_f>2)), sum(alpha_f(idd))];
        end
        
        idd = c_f>=max(validn,max(c_f)-2);
        X(i,:) = [sum(u_f_1(idd,1).*alpha_f(idd)./sig_f_1(idd,1)),sum(u_f_1(idd,2).*alpha_f(idd)./sig_f_1(idd,1))]/sum(alpha_f(idd)./sig_f_1(idd,1));
        
        t=[t,toc];
        Nmea(i) = length(u_mi);
        if(fig~=-1)
            axes(ha(i));
            k=0.001;
            plot(k*u_mi(:,1),k*u_mi(:,2),'b.','linewidth',2);hold on;
            plot(k*gpsdata(i,1),k*gpsdata(i,2),'o','Color',[0,0,0],'linewidth',2.5);hold on;
            %     plot(k*u_f(:,1),k*u_f(:,2),'g.');hold on;
            plot(k*u_f_1(idd,1),k*u_f_1(idd,2),'r.');hold on;
            
            grid on
            axis equal
            if fig==1
                %         xlim([-1.3 1])
                %         xlim([-1.3 1])
                xlim([gpsdata(i,1)/1000-0.05 gpsdata(i,1)/1000+0.05])
                ylim([gpsdata(i,2)/1000-0.05 gpsdata(i,2)/1000+0.05])
            else
                xlim([-30 75])
                ylim([-75 30])
            end
            %     title(['timestep ',num2str(i)])
        end
    end
    if(fig~=-1)
        set(figure(1),'Position',[2800,-400,500,1400]);
    end
    % annotation('textbox',[0.09,0.88,0.1,0.1], 'LineStyle','-','LineWidth',0.1,'String',"Km")
    
    %% filtering result
    % 生成轨迹数据
    traj = [];
    load('dx.mat');
    for ii = validn:length(gps_id)-1
        traj(end+1,:) = X(ii,:);
        for jj = gps_id(ii):gps_id(ii+1)
            traj(end+1,:) = traj(end,:) + dx(ii,:,jj);
        end
    end
    if(fig~=-1)
        figure
        hold on;
        plot(traj(:,1),traj(:,2),'r.')
        plot(truetraj(:,1),truetraj(:,2),'--','Color',[0,0,0],'linewidth',1.5)
        
        quiver(truetraj(1,1),truetraj(1,2),truetraj(10,1)-truetraj(1,1),truetraj(10,2)-truetraj(1,2), 10,'Color',[0,0,0],'linewidth',1.5,'MaxHeadSize',15);
        plot(truetraj(1,1),truetraj(1,2),'x','Color',[0,0,0],'linewidth',2)
        
        
        plot(gpsdata(:,1),gpsdata(:,2),'o','Color',[0,0,0]);
        axis equal
        hold on; grid on;
        plot(X(validn:end,1),X(validn:end,2),'r*')
        plot(imgdata_closest(1:end,1),imgdata_closest(1:end,2),'b+')
        
        
        xlim([-1200 700])
        ylim([-700 200]);
        title('Trajectory in ENU')
        xlabel('East (m)')
        ylabel('North (m)')
        % title('o filter�?* gps, + closest measurement')
    end
    matcherr = gpsdata-imgdata_closest;
    matcherr = sqrt(matcherr(:,1).^2 + matcherr(:,2).^2);
    
    err = gpsdata - X;
    err = sqrt(err(:,1).^2 + err(:,2).^2);
    
    % filter = imgdata_closest- X;
    % filter = sqrt(filter(:,1).^2 + filter(:,2).^2);
    
%     mean(matcherr(2:end))
%     max(matcherr(2:end))
%     mean(err(6:end))
%     max(err(6:end))
    
    %% fig1
    if(fig~=-1)
        % validn = 1;
        figure;
        % ha = tight_subplot(4,1,[.05 .025],[.1 .1],[.1 .1]);
        % axes(ha(1));
        subplot(4,1,1);
        plot(err); hold on
        plot(err,'*');grid on;
        
        title('position error of filter (m)')
        % xlabel('timestep')
        xlim([validn 12]);
        ylim([0 max(0.1,max(err(7:end)))]);
        
        % axes(ha(2));
        subplot(4,1,2);
        plot(matcherr);hold on
        plot(matcherr,'*');grid on;
        title('position error of matching (m)')
        % xlabel('timestep')
        xlim([validn 12]);
        ylim([0 10]);
        
        % axes(ha(3));
        subplot(4,1,3);
        plot(S2(:,1));hold on
        plot(S2(:,1),'*');
        title('credible indicator S_1'); grid on
        % xlabel('timestep')
        plot(0.5*ones(1,12),'r--');
        xlim([validn,12])
        ylim([0.0,1])
        
        % axes(ha(4));
        subplot(4,1,4);
        plot(log(S1/(sig_m(1)+sig_v(1))*2*5));hold on
        plot(log(S1/(sig_m(1)+sig_v(1))*2*5),'*');
        
        grid on
        hold on
        plot(zeros(1,12),'r--');
        title('credible indicator S_2')
        xlabel('timestep')
        xlim([validn,12])
        
    end
    S = log(S1(:)/(sig_m(1)+sig_v(1))*2*5);
end

%%
if loop>1
figure;
yyaxis left
hold on;
err = looperr(:,1);
plot(loopX,(looperr(:,1)));
ylim([0 9]);
ylabel('mean position error of filtering (m)','FontName','Times New Roman','FontSize',14)
yyaxis right
plot(loopX,(looperr(:,2)));
grid on
xlabel('\mu')
ylabel('mean position error of registration (m)','FontName','Times New Roman','FontSize',14)
% ylabel('error/m')
% h = legend( 'mean position error (m) of filtering','mean position error (m) of registration');
% set(h,'FontName','Times New Roman','FontSize',12,'FontWeight','normal')

figure;
yyaxis left
plot(loopX,(loopD(:,2)));
hold on;
ylim([0 7]);
ylabel('number of convergence','FontName','Times New Roman','FontSize',14)
yyaxis right
plot(loopX,loopN);
grid on
xlabel('\mu')
ylabel('mean number of matching results','FontName','Times New Roman','FontSize',14)
% ylabel('error/m')
% h = legend('number of convergence','number of matching results');
% set(h,'FontName','Times New Roman','FontSize',12,'FontWeight','normal')
end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function readdata()

global NUM gpsdata x_v u_m alpha_m  truetraj gps_id
%% 读取gps真�??
file_path =  './test/';
img_path_list = dir(strcat(file_path,'*.jpg'));
NUM = 12;

% gps = importdata('gps.txt');
load('gps.mat')
for j=2:length(gps.textdata)
    t_gps(j) =  str2num(gps.textdata{j})/10e8;
end

for i=1:NUM
    t = str2num(img_path_list(i).name(1:15));
    dt = abs(t_gps  - t);
    gps_id(i,:) = find(dt == min(dt));
end
for i=1:NUM
    gpsdata(i,:) = [str2num(gps.textdata{gps_id(i),9}),str2num(gps.textdata{gps_id(i),7})];
end


%% 经纬度转�?
% [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid)
lon0 = gpsdata(1,1);
lat0 = gpsdata(1,2);

for i=1:NUM
    lat = gpsdata(i,2);
    lon = gpsdata(i,1);
    [x,y,~] = geodetic2ned(lat,lon,0,lat0,lon0,0,wgs84Ellipsoid);
    gpsdata(i,:)  = [y,x];
end

for i = 2:length(gps.textdata)
    truetraj(i-1,1:2) = [str2num(gps.textdata{i,9}),str2num(gps.textdata{i,7})];
    lat = truetraj(i-1,2);
    lon = truetraj(i-1,1);
    [x,y,~] = geodetic2ned(lat,lon,0,lat0,lon0,0,wgs84Ellipsoid);
    truetraj(i-1,:)  = [y,x];
end


% figure;plot(gpsdata(:,1),gpsdata(:,2),'o'); hold on; plot(gpsdata(1,1),gpsdata(1,2),'*')

%% 相对位移估计�?
load('x_v.mat');
%% 分解匹配�?
load('u_m.mat')
load('alpha_m.mat')
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
function initialization()
global  gpsdata S1 ha
global  Nps sigma estsscore Nmea fig
global  alpha_m u_m imgdata_closest
global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1 sig_m mu_registration

% t-1测量
alpha_mi = alpha_m(1,:);
alpha_mi = alpha_mi - mu_registration*max(alpha_mi);

u_mi = u_m(1,1:2*length(alpha_mi));
u_mi = reshape(u_mi,[2,length(alpha_mi)])';

id = alpha_mi>0;
u_mi = u_mi(id,:);
alpha_mi = alpha_mi(id);
alpha_mi = alpha_mi/sum(alpha_mi);

% 粒子连续观测到的次数,应该和estsscore，premup�?�?对应
imgerr(1) = 1000;
for j=1:size(u_mi,1)
    d = norm(u_mi(j,:) - gpsdata(1,:));
    if d<imgerr(1)
        imgerr(1) = d;
        imgdata_closest(1,:) = u_mi(j,:);
    end
end
S1(1) = sqrt((max(u_mi(:,1)) - min(u_mi(:,1)))^2 + (max(u_mi(:,2)) - min(u_mi(:,2)))^2);
Nmea(1) = size(u_mi,1);


u_f_1 = u_mi;
sig_f_1 = sig_m.*ones(size(u_mi));
alpha_f = alpha_mi;
c_f = ones(size(alpha_mi));

if(fig~=-1)
    % figure
    % subplot(3,4,1)
    axes(ha(1));
    k=0.001;
    plot(k*u_mi(:,1),k*u_mi(:,2),'b.','linewidth',2);hold on;
    % plot(k*gpsdata(1,1),k*gpsdata(1,2),'o','Color',[0,0,0]);hold on;
    plot(k*gpsdata(1,1),k*gpsdata(1,2),'o','Color',[0,0,0],'linewidth',2.5);hold on;
    
    grid on
    axis equal
    
    if fig==1
        %     xlim([-1.3 1])
        %     ylim([-1.3 1])
        xlim([gpsdata(1,1)/1000-0.05 gpsdata(1,1)/1000+0.05])
        ylim([gpsdata(1,2)/1000-0.05 gpsdata(1,2)/1000+0.05])
    else
        xlim([-30 75])
        ylim([-75 30])
    end
    set(gca,'xtick',-0.05:0.02:0.05)
    set(gca,'ytick',-0.05:0.02:0.05)
end
% title(['timestep ',num2str(1)])
% text(-15,22,['timestep ',num2str(i)],'FontSize',10);
%     xlim([-12000 2000])
%     ylim([-2000 1000])
end
function propagation()
global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1
global x_v sig_v
global i
% t预测, 权�?�不�?
u_f   = u_f_1   + x_v(i,:) + rand(size(u_f_1)).*sig_v(1,:)*0.05;
sig_f = sig_f_1 + sig_v(1,:);
alpha_f = alpha_f;
c_f = c_f;
end
function measurement()
global u_m sig_m alpha_m  u_mi
global i gpsdata imgdata_closest imgerr measure_id alpha_mi mu_registration
% t测量
alpha_mi = alpha_m(i,:);
alpha_mi = alpha_mi - max(alpha_mi)*mu_registration;   % 剔除多余匹配

u_mi = u_m(i,:);
u_mi = reshape(u_mi,[2,length(u_mi)/2])';

measure_id = alpha_mi>0;
alpha_mi = alpha_mi(measure_id);
%     alpha_mi = alpha_mi/sum(alpha_mi);
u_mi = u_mi(measure_id,:);

Nm(i) = length(u_mi);

% 计算匹配误差
imgerr(i) = 10000;
for j=1:size(u_mi,1)
    d = norm(u_mi(j,:) - gpsdata(i,:));
    if d<imgerr(i)
        imgerr(i) = d;
        imgdata_closest(i,:) = u_mi(j,:);
    end
end

% if i==11
%     i
% end
% 更新w
MeanValue(u_mi, alpha_mi);
end
function resampling()

global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1
global u_m sig_m alpha_m
global  X S1 S2
global Nps validn i idd
% 粒子 重采�?
%     alpha_f = alpha_f(:)/sum(alpha_f(:));
alpha_f = alpha_f(:)/1;
%     [u_f_1, alpha_f, c_f] = resamp2(u_f, alpha_f, c_f, 0, Nps);
resmp();
%     alpha_f = alpha_f/sum(alpha_f);

% 增加测量值粒�?
sscore = alpha_m(i,:);
for j=1:length(sscore)
    if(sscore(j) > 0.6)
        alpha_f(end+1) = 2/Nps;
        u_f_1(end+1,:) = u_m(i,j*2-1:j*2);
        c_f(end+1,:) = 1;
        sig_f_1(end+1,:) = sig_m;
    end
end

end
function resmp()
global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1  Nps
alpha_f_s = [];
u_f_s = [];
sig_f_s = [];
c_f_s = [];
for i = 1:length(alpha_f)
    npi(i) = round(Nps * alpha_f(i)/sum(alpha_f));
    %     if(npi(i)>=1 && c_f(i)>0)
    %         u_f_s(end+1,:) = u_f_1(i,:);
    %         alpha_f_s(end+1,:) = alpha_f(i)*npi(i);
    %         c_f_s(end+1,:) = c_f(i);
    %         sig_f_s(end+1,:) = sig_f_1(i,:);
    %     end
    Nps = max(100, Nps/1.2);
    for j=1:npi(i)
        u_f_s(end+1,:) = u_f_1(i,:);
        alpha_f_s(end+1,:) = alpha_f(i);
        c_f_s(end+1,:) = c_f(i);
        sig_f_s(end+1,:) = sig_f_1(i,:);
    end
    
end
alpha_f = alpha_f_s;
u_f_1 = u_f_s;
sig_f_1 = sig_f_s;
c_f = c_f_s;
end
function MeanValue(u_mi, alpha_mi)
global u_f sig_f alpha_f c_f
global u_f_1 sig_f_1
global u_m sig_m  Nps

for i = 1:size(u_f,1)
    s = 0;
    inline = 0;
    alpha_fd = 0;
    sig_fd = [0,0];
    u_fd = [0,0];
    
    for j=1:size(u_mi,1)
        d = norm(u_f(i,:) - u_mi(j,:));
        if(d<min(sig_m(1)+sig_f(1),100))
            % 测量更新
            %             u_f_1(i,:) = (u_f(i,:) * alpha_f(i) / sig_f(i,1) + u_mi(j,:) * alpha_mi(j) / sig_m(1,1))...
            %                     /( alpha_f(i) / sig_f(i,1) +  alpha_mi(j) / sig_m(1,1));
            if(alpha_mi(j)>s) % 保留范围内权值最大的
                s = alpha_mi(j);
                % 更新位置、方差�?�权重�?�次�?
                inline = 1;
                alpha_fd = alpha_mi(j);
                sig_fd = sig_m;
                u_fd =  u_mi(j,:);
            end
        end
    end
    
    % 权重怎么取很重要
    if(inline>0)
        u_f_1(i,:) = (u_f(i,:) * alpha_f(i) / sig_f(i,1) + u_fd * alpha_fd / sig_m(1,1))...
            /( alpha_f(i) / sig_f(i,1) + alpha_fd / sig_m(1,1));
        
        alpha_f(i) = alpha_f(i) + alpha_fd;
        c_f(i) = c_f(i) + 1;
        %         sig_f_1(i,:) = sig_f(i,:).* sig_f(i,:)./(sig_f(i,:) + sig_fd);
        %         alpha_f(i) = alpha_f(i)/length(alpha_mi)*length(alpha_f) + alpha_fd;
        %         alpha_f(i) = alpha_f(i)+alpha_fd;
        
        %          alpha_f(i) = (alpha_f(i) / sig_f(i,1) + alpha_fd / sig_m(1,1))...
        %                     /(  sig_f(i,1) +  sig_m(1,1));
        
        sig_f_1(i,:) = sig_fd.*sig_fd./sig_f_1(i,:); % sig_f(i,:).* sig_f(i,:)./(sig_f(i,:) + sig_fd);
    else
        u_f_1(i,:) = u_f(i,:);
        alpha_f(i) = alpha_f(i)/2;
        sig_f_1(i,:) = sig_f(i,:)*1.2;
        c_f(i) = max(c_f(i) - 0.5,0);
    end
end



% % 如果est在cur+sigma范围内，获得score值，否则score=0？或者衰减？
% z = 0;
% Z = [];
% for i = 1:length(estmup)
%     z = presscore(i)/3;
%     for j=1:length(curmup)
% %         d = abs(estmup(i,:) - curmup(j,:));
% %         if(d(1)<sigma(1) && d(2)<sigma(2))
%         d = norm(estmup(i,:) - curmup(j,:));
%         if(d<sigma(1))
%             z = max(z,cursscore(j));
%         end
%     end
%     Z(end+1) = z;
% end
end


