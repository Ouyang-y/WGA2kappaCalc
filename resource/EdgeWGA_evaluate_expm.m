function I = EdgeWGA_evaluate_expm(WGNum,L,kappa,options)
% EdgeWGA_evaluate_expm(WGNum,L,kappa,EdgeKappa) uniformly spaced central incidence waveguide array evaluation
%   I = EdgeWGA_evaluate_expm(WGNum,L,kappa,options)
%   I = EdgeWGA_evaluate_expm(WGNum,L,kappa,"EdgeKappa",0.4,"beta",0)
%   I = EdgeWGA_evaluate_expm(13,8.5,0.5,"EdgeKappa",0.45,"beta",0,"EdgeBeta",0,"isEnd",false,"EvalutionNum",1000)
%
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[I] = EdgeWGA_evaluate(WGNum, L, kappa ...
%                   [, 'EdgeKappa', 0 ...
%                    , 'beta', 0 ...
%                    , 'EdgeBeta', 0 ...
%                    , 'isEnd', false ...
%                    , 'EvalutionNum', 1000]);
%
% Params:
%   - WGNum         [required]  [integer; >1] 阵列波导
%   - L             [required]  [positive] 演化长度
%   - kappa         [required]  [numeric] 耦合系数(哈密顿量次对角线
%   - EdgeKappa     [namevalue] [numeric] 边界耦合系数(哈密顿量次对角线第一个和最后一个
%   - beta          [namevalue] [numeric] 传播常数(哈密顿量主对角线
%   - EdgeBeta      [namevalue] [numeric] 边界传播常数(哈密顿量主对角线第一个和最后一个
%   - isEnd         [namevalue] [numeric] 是否只计算L处的结果
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
    options.EdgeKappa (1,1) {mustBeFinite} = kappa
    options.beta (1,1) {mustBeFinite} = 0
    options.EdgeBeta (1,1) {mustBeFinite} = 0
    options.isEnd (1,1) {mustBeNumericOrLogical} = false
    options.EvalutionNum (1,1) {mustBeInteger,mustBeGreaterThan(options.EvalutionNum,2)} = 1000
end
H = diag(ones(1,WGNum)*options.beta)+diag(ones(1,WGNum-1)*kappa,-1)+diag(ones(1,WGNum-1)*kappa,1);
[H(1,2),H(2,1),H(end-1,end),H(end,end-1)] = deal(options.EdgeKappa);
[H(1,1),H(end,end)] = deal(options.EdgeBeta);
array(:,1)=zeros(1,WGNum);
array(ceil(WGNum/2),1)=1;  % 态矢
if options.isEnd
    Ltemp = L;
else
    Ltemp = linspace(0,L,options.EvalutionNum);
    Ltemp = repmat(Ltemp,WGNum,1);
    array = repmat(array,1,options.EvalutionNum);
end
I = mat2gray(abs(expm(-1i*Ltemp.*H)).^2.*array);
end