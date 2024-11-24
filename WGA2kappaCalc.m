clear
close all
mfilePath = mfilename("fullpath");
addpath([mfilePath(1:end-length(mfilename)),'\resource'])
load("OSI_rainbow.mat")
% Path = {                  % dis = 11μm
%     '.\test\WGA-1.tiff';  % L = 20-0.22mm  (0.22为抛光去除量，见.\test\wga_01.bmp
%     '.\test\WGA-2.tiff';  % L = 17
%     '.\test\WGA-3.tiff';  % L = 14
%     '.\test\WGA-4.tiff';  % L = 11
%     '.\test\WGA-5.tiff';  % L = 8
%     '.\test\WGA-6.tiff';  % L = 5
%     '.\test\WGA-7.tiff';  % L = 2
%     };
% printP = "";
Path = {                    % dis = 9μm
    % '.\test\WGA-8.tiff';    % L = 20
    % '.\test\WGA-9.tiff';    % L = 17
    % '.\test\WGA-10.tiff';   % L = 14
    % '.\test\WGA-11.tiff';   % L = 11
    % '.\test\WGA-12.tiff';   % L = 8
    '.\test\WGA-13.tiff';   % L = 5
    '.\test\WGA-14.tiff';   % L = 2
    };
printP = "2_";  % 输出图片名字前缀，可为空，不可没有该变量
num = 13;  % 阵列波导数
dL = -0.22;  % 从俯视图测量WGA与设计的误差
L = (5:-3:1)+dL;  % 对应path顺序的耦合长度
CalculateNum = 200;  % kappa计算数量
kappa_Calculate = linspace(0.55,0.65,CalculateNum);  % kappa计算范围
outputPicType = ".pdf";

outputPath = fileparts(fileparts(Path{1}))+"\output\";  % 输出文件夹 = Path{1}与test同级的output文件夹中
if ~exist(outputPath,"dir"),mkdir(outputPath),end
% load("OSI_rainbow.mat")
%% set Kappa calculate in total
KappaCalculatePath = outputPath + "KappaCalculate.mat";
if exist(KappaCalculatePath,"file"),load(KappaCalculatePath);end
fprintf("Kappa calculate Checking...\n")
if ~exist("kappa_Calculate0","var")||logical(CalculateNum-CalculateNum0)||...  % 计算过>计算数量>L数量>计算内容
        logical(length(L)-length(L0))||any([kappa_Calculate,num,L]-[kappa_Calculate0,num0,L0])
    Output = nan(num,length(kappa_Calculate),length(L));
    for Ltemp = 1:length(L)
        for tt = 1:length(kappa_Calculate)
            Wtemp = UWGA_evaluate(num,L(Ltemp),kappa_Calculate(tt));
            Output(:,tt,Ltemp) = Wtemp(:,end);
        end
    end
    [kappa_Calculate0,CalculateNum0,num0,L0]=deal(kappa_Calculate,CalculateNum,num,L);
    save(KappaCalculatePath,"Output","kappa_Calculate0","CalculateNum0","num0","L0");
    clearvars kappa_Calculate0 CalculateNum0 num0 L0
end
fprintf("Kappa calculate finished!\n")
%% Fit Calc
% 关闭拟合数目、Nan Warning提示
warning('off','curvefit:prepareFittingData:sizeMismatch')
warning('off','curvefit:prepareFittingData:removingNaNAndInf')
% 单线高斯拟合设置
ft_gs = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 5];  % 限制单线高斯拟合的宽度
opts.Robust = 'LAR';
% 二维高斯拟合设置
ft = fittype('A*exp((-(x-X0).^2/sigmaX2-(y-Y0).^2/sigmaY2)/2)', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "sigmaX2" "sigmaY2"]);
% 初始化结果存储
PathNum = size(Path,1);
Rsqure = nan(PathNum,CalculateNum);
[res.A,res.X,res.Y,res.WX,res.WY,res.rangeRight,res.rangeleft] = ...
    deal(nan(PathNum,num));
% 按图开始计算
for temp = 1:PathNum
    fprintf("%d/%d:%s --- Calculating...\n",temp,PathNum,Path{temp})
    % debuge用做断点
    % if tempSub==10
    %     XXXXX=1;
    % end
    im = imread(Path{temp});
    imgSize = size(im);
    imgRemoveBG = im(:,:,1);  % 灰度图只需要第一层
    img.(['img0_',num2str(temp)])=imgRemoveBG;
    C = sum(imgRemoveBG);  % 累加平滑曲线为多高斯峰，利于findpeaks
    [peaks,locs]=findpeaks(C, 'MinPeakDistance', imgSize(2)/num/1.3, 'MinPeakHeight', 200);
    % findpeaks不一定找全点，可以在此debug，运行:
    % findpeaks(C, 'MinPeakDistance', imgSize(2)\num/1.3, 'MinPeakHeight', 200)
    % 下面基于等间距分布，对遗漏点进行补充
    minDis = min(diff(locs));  % 以最低峰值间隔为基础进行加点
    locstemp = [1,locs,size(C,2)];  % 加入左边界1和右边界size(C,2)
    % 以最大间隔循环加点
    while length(locstemp) < num+2  % num+2包括左边界和右边界
        diffLocs = diff(locstemp);
        [~,maxLocsNum] = max(diffLocs);
        if maxLocsNum == 1
            locstemp = [locstemp(1),locstemp(2)-minDis,locstemp(2:end)];
        else
            locstemp = [locstemp(1:maxLocsNum),locstemp(maxLocsNum)+minDis,locstemp(maxLocsNum+1:end)];
        end
    end
    locs = locstemp(2:end-1);  % 移除左边界1和右边界size(C,2)
    % A - 各点边界
    diffLocs = diff(locs);
    A = ceil(diffLocs/2+locs(1:end-1));
    A = [locs(1)-floor(diffLocs(1)/2),A,locs(end)+ceil(diffLocs(end)/2)]; %#ok
    A(A<1)=1;A(A>size(C,2))=size(C,2);
    % range - 各点[上边界;下边界]
    range = [A(1:end-1);A(2:end)];
    % 初始化拟合结果
    [restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY] = deal(nan(1,num));
    y = 1:imgSize(1);
    for tempSub = 1:num
        x = range(1,tempSub):range(2,tempSub);
        if sum(sum(imgRemoveBG(y,x)>1))<10  % 大于1的数据点少于10个 - ccd读数误差
            restemp.A(tempSub) = 0;
            continue;
        end
        if C(locs(tempSub))<1200  % 峰值过低 - 单线高斯拟合
            z = double(imgRemoveBG(y,locs(tempSub)));
            [xData, yData] = prepareCurveData(y,z);
            [opts.StartPoint(1),opts.StartPoint(2)] = max(yData);
            opts.StartPoint(3) = length(yData)-sum(yData==0);
            [fitGS, gof] = fit(xData,yData,ft_gs,opts);
            restemp.A(tempSub) = fitGS.a1;
            restemp.Y(tempSub) = fitGS.b1;
            % res.(['GS',num2str(temp),'_',num2str(tempSub)])=sfit(ft,fitGS.a1,locs(tempSub),fitGS.b1,fitGS.c1^2,fitGS.c1^2);
            continue;
        end
        z = double(imgRemoveBG(y,x));
        [fitresult,gof,gauss2,~] = gauss2fit(x,y,z);
        [restemp.A(tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub)] = ...
            deal(gauss2.A,gauss2.X0,gauss2.Y0,gauss2.sigmaX2,gauss2.sigmaY2);
        % [res.rangeUp,res.rangeDown] = deal(range(1,:),)
        res.(['GS',num2str(temp),'_',num2str(tempSub)])=gauss2;
    end
    restemp.A = restemp.A/sum(restemp.A);  % Normalize A
    restemp.X(isnan(restemp.X))=locs(isnan(restemp.X));  % use locs pad nan X
    restemp.Y(isnan(restemp.Y))=deal(mean(restemp.Y(~isnan(restemp.Y))));  % use mean(Y) pad nan Y
    restemp.WX(isnan(restemp.WX))=deal(mean(restemp.WX(~isnan(restemp.WX))));  % use mean(WX) pad nan WX
    restemp.WY(isnan(restemp.WY))=deal(mean(restemp.WY(~isnan(restemp.WY))));  % use mean(WY) pad nan WY
    % save restemp > res
    [res.A(temp,:),res.X(temp,:),res.Y(temp,:),res.WX(temp,:),res.WY(temp,:)] = ...
        deal(restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY);
    [res.rangeRight(temp,:),res.rangeleft(temp,:)] = deal(range(1,:),range(2,:));
    % calc Rsqure
    Outputtemp = Output(:,:,temp);
    IAll = repmat(restemp.A',1,size(Outputtemp,2));
    SSE = sum((IAll-Outputtemp).^2);
    SST = sum((IAll-mean(IAll).^2));
    Rsqure(temp,:) = 1-SSE./SST;
end
%% fig1 - 计算绘制保存R方图
r2 = sum(Rsqure)/temp;  % max(Rsqure) > kappa
[aa,bb] = max(r2);
kappa = kappa_Calculate(bb);
fprintf('Kappa = %f,R^2 = %f\n',kappa,aa);
fprintf("Drawing Figures...\n")
f1 = figure(1);
plot(kappa_Calculate,r2)
title(sprintf('Kappa = %f,R^2 = %f',kappa,aa))
exportgraphics(f1,outputPath+printP+"Rsquare"+outputPicType,ContentType="vector")
%% fig2 - 计算绘制保存模拟(俯视)演化图
WGASim.L = L(1)-dL+1;
WGA = UWGA_evaluate(num,WGASim.L,kappa,"EvalutionNum",(WGASim.L)*100);  % 使用计算出的kappa进行演化
% 在num分布方向上进行高斯拓展(图片美化
addCol = 10;
y = -addCol:addCol;
Y = repmat(y',[num (WGASim.L)*100]);  % [num UWGA_evaluate计算数目]
imgWGA = nan((addCol*2+1)*(num-1)+num,size(WGA,2));  % 预分配图片尺寸
I0 = repelem(WGA,(addCol*2+1),1);  % 把计算演化结果WGA的num行数拓展到addCol上
Im0 = I0 .* exp(-Y.^2/(addCol/2)^2);  % 高斯化
f2 = figure(2);
IMrgb = ind2rgb(int8(ceil(Im0*128)),OSI_rainbow);
imagesc(IMrgb)
f2.Position = [20 80 1800 size(Y,2)+100];
xticks(linspace(0,WGASim.L*100,11))
xticklabels(split(num2str(linspace(0,WGASim.L,11))))
yticks(addCol+1:addCol*2+1:size(imgWGA,1))
yticklabels(split(num2str(num:-1:1)))
drawnow;ax = gca;ax.TickDir='out';ax.FontSize=22;ax.LineWidth=2;
xlabel('L(mm)');ylabel('num');
box off
exportgraphics(f2,outputPath+printP+"Sim"+outputPicType,ContentType="vector")
%% fig3 - 绘制保存模拟(俯视)演化图 带有Path对应长度标记
Lselect = (L-dL)*100;
Sim = WGA(:,Lselect)';
lineHalfWid = 1;
imgWGAt = IMrgb;
for temp = 1:PathNum
    Ls = Lselect(temp)-lineHalfWid:Lselect(temp)+lineHalfWid;
    imgWGAt(:,Ls,:)=ones(size(imgWGAt,1),size(Ls,2),3)*255;
end
f3 = figure(3);
IMrgb = ind2rgb(int8(ceil(Im0*128)),OSI_rainbow);
imagesc(imgWGAt)
f3.Position = [20 80 1800 size(Y,2)+100];
xticks(linspace(0,WGASim.L*100,11))
xticklabels(split(num2str(linspace(0,WGASim.L,11))))
yticks(addCol+1:addCol*2+1:size(imgWGA,1))
yticklabels(split(num2str(num:-1:1)))
drawnow;ax = gca;ax.TickDir='out';ax.FontSize=22;ax.LineWidth=2;
xlabel('L(mm)');ylabel('num');
box off
exportgraphics(f2,outputPath+printP+"SimCut"+outputPicType,ContentType="vector")
%% fig4 -- fig3+PathNum - 各L对应绘图[数据图;findpeaks图;拟合图;kappa演化图]
for temp = 1:PathNum
    figure
    t = tiledlayout(4,1,"TileSpacing","none");
    [restemp.A,restemp.X,restemp.Y,restemp.WX,restemp.WY]=deal( ...
        res.A(temp,:),res.X(temp,:),res.Y(temp,:),res.WX(temp,:),res.WY(temp,:));
    [range(1,:),range(2,:)] = deal(res.rangeRight(temp,:),res.rangeleft(temp,:));
    % 数据图
    nexttile
    ExpImg0 = img.(['img0_',num2str(temp)]);
    imagesc(ExpImg0);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    % findpeaks图
    ax = nexttile;hold on;
    C = sum(ExpImg0);
    plot(C);
    plot(round(restemp.X),C(round(restemp.X)),'ro');
    plot(range(1,:),C(range(1,:)),'ko');
    plot(range(2,:),C(range(2,:)),'mo');
    legend({'','center','begin','end'})
    hold off;axis off;ax.XLim = [1,size(C,2)];ax.YLim(1)=-500;
    % 拟合图;kappa演化图
    SimImg = zeros(size(ExpImg0));
    ExpImg = zeros(size(ExpImg0));
    [X,Y]=meshgrid(1:size(ExpImg0,2),1:size(ExpImg0,1));
    for tempSub = 1:num  % 使用二维高斯函数，复现各光斑拟合/演化结果，并加和
        GausExp = sfit(ft,restemp.A(tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub));
        ExpImg = ExpImg + feval(GausExp,X,Y);
        GausSim = sfit(ft,Sim(temp,tempSub),restemp.X(tempSub),restemp.Y(tempSub),restemp.WX(tempSub),restemp.WY(tempSub));
        SimImg = SimImg + feval(GausSim,X,Y);
    end
    nexttile
    imagesc(ExpImg);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    nexttile
    imagesc(SimImg);colormap(OSI_rainbow(2:end-1,:));grid off;box off;axis off
    [~,name]=fileparts(Path{temp});
    title(t,name)
    exportgraphics(t,outputPath+printP+string(name)+outputPicType,ContentType="vector")
end
fprintf("Done! Result saving in: %s\n",outputPath);
