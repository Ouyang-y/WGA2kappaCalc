function I = WGA_evaluate(WGNum,L,kappa,options)
% WGA_evaluate(WGNum,L,kappa,options) waveguide array evaluation
%   I = EdgeWGA_evaluate_expm(WGNum,L,kappa,options)
%   I = EdgeWGA_evaluate_expm(WGNum,L,kappa,"EdgeKappa",0.4,"beta",0)
%   I = EdgeWGA_evaluate_expm(13,8.5,0.5,"EdgeKappa",0.45,"beta",0,"EdgeBeta",0,"isEnd",false,"EvalutionNum",1000)
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[I] = EdgeWGA_evaluate(WGNum, L, kappa ...
%                   [, 'dl', 0 ...
%                    , 'EdgeKappa', 0 ...
%                    , 'beta', 0 ...
%                    , 'EdgeBeta', 0 ...
%                    , 'isEnd', false ...
%                    , 'EvalutionNum', 1000]);
%
% Params:
%   - WGNum         [required]  [integer; >1] 阵列波导
%   - L             [required]  [positive] 演化长度
%   - kappa         [required]  [numeric] 耦合系数(哈密顿量次对角线
%   - dl            [namevalue] [<=0] 演化长度误差
%   - EdgeKappa     [namevalue] [numeric] 边界耦合系数(哈密顿量次对角线第一个和最后一个
%   - beta          [namevalue] [numeric] 传播常数(哈密顿量主对角线
%   - EdgeBeta      [namevalue] [numeric] 边界传播常数(哈密顿量主对角线第一个和最后一个
%   - isEnd         [namevalue] [logical] 是否只计算L处的结果
%   - EvalutionNum  [namevalue] [integer; >2] 计算数目(划分L的数量
%
% Return:
%   - I 光强大小,行数对应波导数,列数对应演化长度
%
% Matlab Version: R2024b
%
% Author: oyy
arguments
    WGNum (1,1) {mustBeInteger,mustBeGreaterThan(WGNum,1)}
    L (1,1) {mustBePositive}
    kappa (1,1) {mustBeFinite}
    options.dl (1,1) {mustBeNonpositive} = 0
    options.EdgeKappa (1,1) {mustBeFinite} = kappa
    options.beta (1,1) {mustBeFinite} = 0
    options.EdgeBeta (1,1) {mustBeFinite} = 0
    options.isEnd (1,1) {mustBeNumericOrLogical} = false
    options.EvalutionNum (1,1) {mustBeInteger,mustBeGreaterThan(options.EvalutionNum,2)} = 1000
end
H = diag(ones(1,WGNum)*options.beta)+diag(ones(1,WGNum-1)*kappa,-1)+diag(ones(1,WGNum-1)*kappa,1);
[H(1,2),H(2,1),H(end-1,end),H(end,end-1)] = deal(options.EdgeKappa);
[H(1,1),H(end,end)] = deal(options.EdgeBeta);
array=zeros(WGNum,1);
array(ceil(WGNum/2))=1;  % 态矢
if options.isEnd
    I = mat2gray(abs(expm(-1i*L.*H)).^2.*array);
else
    I = nan(WGNum,options.EvalutionNum);
    L0 = linspace(0,L,options.EvalutionNum);
    for temp = 1:options.EvalutionNum
        I(:,temp) = abs(expm(-1i*(L0(temp)+options.dl).*H)).^2*array;
    end
end
end