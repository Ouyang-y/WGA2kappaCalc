clear
close all
num = 13;
% Path = {
%     '.\test\WGA-1.tiff';
%     '.\test\WGA-2.tiff';
%     '.\test\WGA-3.tiff';
%     '.\test\WGA-4.tiff';
%     '.\test\WGA-5.tiff';
%     '.\test\WGA-6.tiff';
%     '.\test\WGA-7.tiff';
%     };
% printP = '';
Path = {
    '.\test\WGA-8.tiff';
    '.\test\WGA-9.tiff';
    '.\test\WGA-10.tiff';
    '.\test\WGA-11.tiff';
    '.\test\WGA-12.tiff';
    '.\test\WGA-13.tiff';
    '.\test\WGA-14.tiff';
    };
printP = '2';
L = 20:-3:1-1;dis = 11e-3;
CalculateNum = 20;
beta0 = 3.4;
kappa_Calculate = linspace(0.45,0.7,CalculateNum);
load("OSI_rainbow.mat")
%% set Kappa calculate in total
cal = 1;
% kappa_Calculate = linspace(0,max(Kappa_DC)*10,length(Kappa_DC)*500);
KappaCalculatePath = [pwd,filesep,'KappaCalculate.mat'];
if exist(KappaCalculatePath,"file")
    load(KappaCalculatePath);
end
if cal||~exist("kappa_Calculate0","var")||any(kappa_Calculate-kappa_Calculate0)
    Output = nan(num,length(kappa_Calculate),length(L));
    for Ltemp = 1:length(L)
        for tt = 1:length(kappa_Calculate)
            Wtemp = WGA_evaluate(num,1000,L(Ltemp),beta0,kappa_Calculate(tt));
            Output(:,tt,Ltemp) = Wtemp(:,end);
        end
    end
    kappa_Calculate0 = kappa_Calculate;
    save(KappaCalculatePath,"Output","kappa_Calculate0");
    cal = 0;
end
%%
warning('off','curvefit:prepareFittingData:sizeMismatch')
warning('off','curvefit:prepareFittingData:removingNaNAndInf')
Rsqure = nan(size(Path,2),CalculateNum);
ft_gs = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );opts.Display = 'Off';
opts.Lower = [0 -Inf 0];
opts.StartPoint = [3 91 0.8696440157339];
opts.Upper = [130 Inf Inf];
ft = fittype('A*exp((-(x-X0).^2/sigmaX2-(y-Y0).^2/sigmaY2)/2)', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "sigmaX2" "sigmaY2"]);
[res.A,res.X,res.Y,res.WX,res.WY] = deal(nan(size(Path,1),num));
for temp = 1:size(Path,1)
    pathtemp = Path{temp};
    im = imread(pathtemp);
    imgSize = size(im);
    imgRemoveBG = im(:,:,1);
    img.(['img0_',num2str(temp)])=imgRemoveBG;
    C = sum(imgRemoveBG);
    % findpeaks(C, 'MinPeakDistance', 50, 'MinPeakHeight', 200, 'MinPeakProminence', 50,'Threshold', 0.5)
    [peaks,locs]=findpeaks(C, 'MinPeakDistance', 30, 'MinPeakHeight', 200, 'MinPeakProminence', 50);
    % if strcmp('9.png',pathtemp),locs(end-2)=[];end
    minDis = min(diff(locs));
    locstemp = [1,locs,size(C,2)];
    while length(locstemp) < num+2
        diffLocs = diff(locstemp);
        [~,maxLocsNum] = max(diffLocs);
        if maxLocsNum == 1
            locstemp = [locstemp(1),locstemp(2)-minDis,locstemp(2:end)];
        else
            locstemp = [locstemp(1:maxLocsNum),locstemp(maxLocsNum)+minDis,locstemp(maxLocsNum+1:end)];
        end
    end
    locs = locstemp(2:end-1);
    diffLocs = diff(locs);
    A = ceil(diffLocs/2+locs(1:end-1));
    A = [locs(1)-floor(diffLocs(1)/2),A,locs(end)+ceil(diffLocs(end)/2)]; %#ok
    A(A<1)=1;A(A>size(C,2))=size(C,2);
    range = [A(1:end-1);A(2:end)];
    [restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY] = deal(nan(1,num));
    y = 1:imgSize(1);
    for tempSub = 1:num
        x = range(1,tempSub):range(2,tempSub);
        % if sum(sum(imgRemoveBG(y,x)>10))<10  %
        %     res.A(tempSub) = 0;
        if C(locs(tempSub))<1500  % 单线高斯拟合
            z = double(imgRemoveBG(y,locs(tempSub)));
            [xData, yData] = prepareCurveData(y,z);
            [fitGS, gof] = fit(xData,yData,ft_gs,opts);
            restemp.A(tempSub) = fitGS.a1;
            % restemp.X(tempSub) = locs(tempSub);
            % restemp.Y(tempSub) = fitGS.b1;
            % res.(['GS',num2str(temp),'_',num2str(tempSub)])=sfit(ft,fitGS.a1,locs(tempSub),fitGS.b1,fitGS.c1^2,fitGS.c1^2);
            continue;
        end
        % end
        % if C(locs(tempSub))<1000,res.A(tempSub) = 0;continue;end
        z = log(double(imgRemoveBG(y,x)));
        [fitresult,gof,gauss2] = gauss2fit(x,y,z);
        restemp.A(tempSub) = gauss2.A;
        restemp.X(tempSub) = gauss2.X0;
        restemp.Y(tempSub) = gauss2.Y0;
        restemp.WX(tempSub) = gauss2.sigmaX2;
        restemp.WY(tempSub) = gauss2.sigmaY2;
        % if strcmp(['GS',num2str(temp),'_',num2str(tempSub)],'GS1_10')
        %     XXXXX=1;
        % end
        res.(['GS',num2str(temp),'_',num2str(tempSub)])=gauss2;
    end
    % Normalize
    restemp.A = restemp.A/sum(restemp.A);
    res.A(temp,:)=restemp.A;
    res.X(temp,:)=locs;
    res.Y(temp,:)=restemp.Y;
    res.WX(temp,:)=restemp.WX;
    res.WY(temp,:)=restemp.WY;
    Outputtemp = Output(:,:,temp);
    IAll = repmat(restemp.A',1,size(Outputtemp,2));
    SSE = sum((IAll-Outputtemp).^2);
    SST = sum((IAll-mean(IAll).^2));
    Rsqure(temp,:) = 1-SSE./SST;

    % figure;clf
    % f=tiledlayout(3,1,"TileSpacing","tight");
    % nexttile
    % imagesc(imgRemoveBG);colormap(OSI_rainbow(2:end-1,:));
    % axis off
    % ax = nexttile;
    % hold on;
    % plot(C);
    % plot(locs,C(locs),'ro');
    % plot(range(1,:),C(range(1,:)),'ko');
    % plot(range(2,:),C(range(2,:)),'mo');
    % legend({'','center','begin','end'})
    % hold off;axis off;ax.XLim = [1,size(imgRemoveBG,2)];
    % nexttile
    % [restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY]=deal( ...
    %     res.A(temp,:),res.X(temp,:),res.Y(temp,:),res.WX(temp,:),res.WY(temp,:));
    % ser = isnan(restemp.Y);
    % numt = sum(ser);
    % restemp.Y(ser)=ones(1,numt)*mean(restemp.Y(~ser));
    % restemp.WX(ser)=ones(1,numt)*mean(restemp.WX(~ser));
    % restemp.WY(ser)=ones(1,numt)*mean(restemp.WY(~ser));
    % ExpImg = zeros(size(imgRemoveBG));
    % [X,Y]=meshgrid(1:size(ExpImg,2),1:size(ExpImg,1));
    % for tempSub = 1:num
    %     GausExp = sfit(ft,restemp.A(tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub));
    %     ExpImg = ExpImg + feval(GausExp,X,Y);
    % end
    % imagesc(ExpImg);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    % [~,name]=fileparts(pathtemp);
    % title(f,name)
    % exportgraphics(f,[name,'.tiff'],'Resolution',300)
end
%% kappa
r2 = sum(Rsqure);
plot(kappa_Calculate,r2/temp)
[aa,bb] = max(r2);
kappa = kappa_Calculate(bb);
fprintf('Kappa = %f',kappa);
%% Sim.
WGA = WGA_evaluate(num,(L(1)+1)*100,L(1)+1,beta0,kappa);
addCol = 10;
y = -addCol:addCol;
Y = repmat(y',[num (L(1)+1)*100]);
imgWGA = nan((addCol*2+1)*(num-1)+num,size(WGA,2));
I0 = repelem(WGA,(addCol*2+1),1);
Im0 = I0 .* exp(-Y.^2/(5)^2);
f = figure;
IMrgb = ind2rgb(int8(ceil(Im0*128)),OSI_rainbow);
imagesc(IMrgb)
f.Position = [20 20 1800 600];
% xlim([0,2100])
xticks(linspace(0,L(1)*100,11))
xticklabels(split(num2str(linspace(0,L(1),11))))
yticks(addCol+1:addCol*2+1:size(imgWGA,1))
yticklabels(split(num2str(num:-1:1)))
drawnow;ax = gca;ax.TickDir='out';ax.FontSize=22;ax.LineWidth=2;
xlabel('L(mm)');ylabel('num');
box off
exportgraphics(f,['Sim',printP,'.tiff'],'Resolution',300)
%% Sim. Cut
Lselect = L*100+1;
Sim = WGA(:,Lselect)';
SSE = sum((Sim'-res.A').^2);
SST = sum((Sim'-mean(Sim').^2));
Rsqure = 1-SSE./SST;
% f=figure;
% tiledlayout(5,2,"TileSpacing","tight")
% nexttile([1 2])
lineHalfWid = 2;
imgWGAt = IMrgb;
for temp = 1:size(L,2)
    Ls = Lselect(temp)-lineHalfWid:Lselect(temp)+lineHalfWid;
    imgWGAt(:,Ls,:)=ones(size(imgWGAt,1),size(Ls,2),3)*255;
end
f = figure;
IMrgb = ind2rgb(int8(ceil(Im0*128)),OSI_rainbow);
imagesc(imgWGAt)
f.Position = [20 20 1800 600];
% xlim([0,2100])
xticks(linspace(0,L(1)*100,11))
xticklabels(split(num2str(linspace(0,L(1),11))))
yticks(addCol+1:addCol*2+1:size(imgWGA,1))
yticklabels(split(num2str(num:-1:1)))
drawnow;ax = gca;ax.TickDir='out';ax.FontSize=22;ax.LineWidth=2;
xlabel('L(mm)');ylabel('num');
box off
exportgraphics(f,['SimCut',printP,'.tiff'],'Resolution',300)
%% Exp Sim
for temp = 1:size(Path,1)
    figure
    t = tiledlayout(4,1,"TileSpacing","none");
    [restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY]=deal( ...
        res.A(temp,:),res.X(temp,:),res.Y(temp,:),res.WX(temp,:),res.WY(temp,:));
    ser = isnan(restemp.Y);
    numt = sum(ser);
    restemp.Y(ser)=ones(1,numt)*mean(restemp.Y(~ser));
    restemp.WX(ser)=ones(1,numt)*mean(restemp.WX(~ser));
    restemp.WY(ser)=ones(1,numt)*mean(restemp.WY(~ser));
    ExpImg0 = img.(['img0_',num2str(temp)]);
    SimImg = zeros(size(ExpImg0));
    ExpImg = zeros(size(ExpImg0));
    [X,Y]=meshgrid(1:size(ExpImg0,2),1:size(ExpImg0,1));
    for tempSub = 1:num
        GausExp = sfit(ft,restemp.A(tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub));
        ExpImg = ExpImg + feval(GausExp,X,Y);
        GausSim = sfit(ft,Sim(temp,tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub));
        SimImg = SimImg + feval(GausSim,X,Y);
    end
    nexttile
    imagesc(ExpImg0);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    ax = nexttile;
    hold on;
    C = sum(ExpImg0);
    plot(C);
    plot(restemp.X,C(restemp.X),'ro');
    plot(range(1,:),C(range(1,:)),'ko');
    plot(range(2,:),C(range(2,:)),'mo');
    legend({'','center','begin','end'})
    hold off;axis off;ax.XLim = [1,size(imgRemoveBG,2)];ax.YLim(1)=-500;
    nexttile
    imagesc(ExpImg);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    nexttile
    imagesc(SimImg);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    pathtemp = Path{temp};
    [~,name]=fileparts(pathtemp);
    title(t,name)
    exportgraphics(t,[name,'.tiff'],'Resolution',300)
end
% nexttile
% add = 40;
% y = -addCol:addCol;
% Y = repmat(y',[num length(y)]);
% X = repmat(y,length(Y),1);
% imgWGAt = nan((addCol*2+1)*(num-1)+num,size(WGA,2));
% I0 = repelem(WGA(:,Lselect(1)),(addCol*2+1),(addCol*2+1));
% Im0 = I0 .* exp(-(Y.^2+X.^2)/(20)^2);
% imagesc(Im0')
% nexttile
% imagesc(img.(['img0_',num2str(1)]))
% %%
% for temp = 1:size(Path,1)
% figure(temp);clf
%     Sim = WGA_evaluate(num,1000,L(temp),beta0,kappa);
%
%
% end
%%

function [fitresult,gof,gauss2] = gauss2fit(x,y,z)
%gauss2fit(X, Y, Z) 二维高斯拟合
%  [f,g,gauss2_WG] = gauss2fit(1:size_WG(2),1:size_WG(1),log(double(mf_WG)));
%  二维高斯拟合
%
%  https://blog.csdn.net/u012366767/article/details/90743083
%
%  f(x,y) = A*exp(-(x-X0)^2/(2*sigmaX2)-(y-Y0)^2/(2*sigmaY2))
%
%  ln(f) = lnA - (x-X0)^2/(2*sigmaX2) - (y-Y0)^2/(2*sigmaY2)
%  ln(f) = p00 + p10*x + p01*y + p20*x^2 + p02*y^2;
%  p00 = lnA - X0^2/(2*sigmaX2) - Y0^2/(2*sigmaY2)
%  p10 = X0/sigmaX2
%  p01 = Y0/sigmaY2
%  p20 = -0.5/sigmaX2
%  p02 = -0.5/sigmaY2
%
%

% 利用多项式拟合
[xData,yData,zData] = prepareSurfaceData(x,y,z);
ft = fittype('poly22');
opts = fitoptions('Method','LinearLeastSquares');
opts.Lower = [-Inf 0 0 -Inf 0 -Inf];
opts.Upper = [Inf Inf Inf 0 0 0];
% opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
% opts.Upper = [Inf Inf Inf Inf Inf Inf];
[fitresult,gof] = fit([xData,yData],zData,ft,opts);
% 将多项式拟合结果转换为二维高斯结果
ft = fittype('A*exp((-(x-X0).^2/sigmaX2-(y-Y0).^2/sigmaY2)/2)', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "sigmaX2" "sigmaY2"]);
sigmaX2 = -0.5/fitresult.p20;
sigmaY2 = -0.5/fitresult.p02;
X0 = fitresult.p10*sigmaX2;
Y0 = fitresult.p01*sigmaY2;
A = exp(fitresult.p00+X0^2/(2*sigmaX2)+Y0^2/(2*sigmaY2));
gauss2 = sfit(ft,A,X0,Y0,sigmaX2,sigmaY2);
end

function I = WGA_evaluate(WGNum,EvalutionNum,L,beta,kappa)
%WGA_EVALUATE 此处显示有关此函数的摘要
%   此处显示详细说明
% WGNum=13;
% EvalutionNum=1000;
% L=8.5;beta=1;kappa=0.0001;
Ltemp = linspace(0,L,EvalutionNum);
dL = Ltemp(2);
H = diag(ones(1,WGNum)*beta)+diag(ones(1,WGNum-1)*kappa,-1)+diag(ones(1,WGNum-1)*kappa,1);
array = nan(WGNum,EvalutionNum);
array(:,1)=zeros(1,WGNum);
array(ceil(WGNum/2),1)=1;  %态矢
for t=2:EvalutionNum
    array(:,t)=(H*dL+2i*eye(WGNum))/(-H*dL+2i*eye(WGNum))*array(:,t-1);
end
I = mat2gray(array.*conj(array));
% imshow(I)
end

function img = removeBackground(img0)
[~, ~, z] = size(img0);
if z > 3
    img0(:,:,4:end) = [];
end
% 除去非灰度背景
background = logical(double(img0(:,:,1)) - double(img0(:,:,2)));
img = img0(:,:,1);
img(background) = 0;
end
