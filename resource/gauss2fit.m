function [fitresult,gof,gauss2,WB] = gauss2fit(x,y,z)
% gauss2fit(X, Y, Z) 二维高斯拟合
%   二维高斯拟合
%
%   https://blog.csdn.net/u012366767/article/details/90743083
%
%   f(x,y) = A*exp(-(x-X0)^2/(2*sigmaX2)-(y-Y0)^2/(2*sigmaY2))
%
%   ln(f) = lnA - (x-X0)^2/(2*sigmaX2) - (y-Y0)^2/(2*sigmaY2)
%   ln(f) = p00 + p10*x + p01*y + p20*x^2 + p02*y^2;
%   p00 = lnA - X0^2/(2*sigmaX2) - Y0^2/(2*sigmaY2)
%   p10 = X0/sigmaX2
%   p01 = Y0/sigmaY2
%   p20 = -0.5/sigmaX2
%   p02 = -0.5/sigmaY2
%
% 以上拟合与Whole Beam Fit Eq仍有差距,可以依赖以下四行在matlab中运算
%   syms xxe y ye sX2 sY2 Al A2 Wx Wy AO
%   A1 = A*exp(-(x-x0)^2/2/sX2-(y-y0)^2/2/sY2)
%   A2 = A*exp(-2*(((x-x8)/(wx/2))^2-((y-y0)/(wy/2))^2))
%   solve(2*sX2==Wx^2/8)
%   最后可以得到WB
%
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[fitresult,gof,gauss2,WB] = gauss2fit(x,y,z);
%
% Params:
%   - x [required]  [vector] x轴
%   - y [required]  [vector] y轴
%   - z [required]  [vector] z值,double
%
% Return:
%   - fitresult 多项式拟合结果
%   - gof       拟合gof
%   - gauss2    二维高斯函数
%   - WB        Whole Beam Fit Eq
%
% Matlab Version: R2024b
%
% Author: oyy

% 利用多项式拟合
[xData,yData,zData] = prepareSurfaceData(x,y,log(z));
ft = fittype('poly22');
opts = fitoptions('Method','LinearLeastSquares');
opts.Lower = [-Inf 0 0 -Inf 0 -Inf];
opts.Upper = [Inf Inf Inf 0 0 0];
% opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
% opts.Upper = [Inf Inf Inf Inf Inf Inf];
[fitresult,gof] = fit([xData,yData],zData,ft,opts);
sigmaX2 = -0.5/fitresult.p20;
sigmaY2 = -0.5/fitresult.p02;
% 将多项式拟合结果转换为二维高斯结果
ft = fittype('A*exp((-(x-X0).^2/sigmaX2-(y-Y0).^2/sigmaY2)/2)', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "sigmaX2" "sigmaY2"]);
X0 = fitresult.p10*sigmaX2;
Y0 = fitresult.p01*sigmaY2;
A = exp(fitresult.p00+X0^2/(2*sigmaX2)+Y0^2/(2*sigmaY2));
gauss2 = sfit(ft,A,X0,Y0,sigmaX2,sigmaY2);
% 将多项式拟合结果转换为全光束拟合(Whole Beam Fit Equations)结果
WBft = fittype('A*exp(-2*(((x-X0)./(wX/2)).^2+((y-Y0)./(wY/2)).^2))', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "wX" "wY"]);
wX=sqrt(sigmaX2)*4;
wY=sqrt(sigmaY2)*4;
A = exp(fitresult.p00+X0^2/(2*sigmaX2)+Y0^2/(2*sigmaY2));
WB = sfit(WBft,A,X0,Y0,wX,wY);
end