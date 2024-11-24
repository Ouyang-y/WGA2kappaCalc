function I = UWGA_evaluate(WGNum,L,kappa,options)
% UWGA_evaluate(WGNum,EvalutionNum,L,beta,kappa) uniformly spaced central incidence waveguide array evaluation
%   I = UWGA_evaluate(WGNum,L,kappa,options)
%   I = UWGA_evaluate(WGNum,L,kappa,"EvalutionNum",1000,"beta",0)
%   I = UWGA_evaluate(13,8.5,0.5,"EvalutionNum",1000,"beta",0)
%
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[I] = UWGA_evaluate(WGNum, L, kappa ...
%                   [, 'EvalutionNum', 1000 ...
%                    , 'beta', 0]);
%
% Params:
%   - WGNum         [required]  [integer; >2] 阵列波导
%   - L             [required]  [positive] 演化长度
%   - kappa         [required]  [numeric] 耦合系数(哈密顿量次对角线
%   - EvalutionNum  [namevalue] [integer; >2] 计算数目(划分L的数量
%   - beta          [namevalue] [numeric] 传播常数(哈密顿量主对角线
%
% Return:
%   - I 光强矩阵,行数对应波导数
%
% Matlab Version: R2024b
%
% Author: oyy
arguments
    WGNum (1,1) {mustBeInteger,mustBeGreaterThan(WGNum,2)}
    L (1,1) {mustBePositive}
    kappa (1,1) {mustBeFinite}
    options.EvalutionNum (1,1) {mustBeInteger,mustBeGreaterThan(options.EvalutionNum,2)} = 1000
    options.beta (1,1) {mustBeFinite} = 0
end
Ltemp = linspace(0,L,options.EvalutionNum);
dL = Ltemp(2);
H = diag(ones(1,WGNum)*options.beta)+diag(ones(1,WGNum-1)*kappa,-1)+diag(ones(1,WGNum-1)*kappa,1);
array = nan(WGNum,options.EvalutionNum);
array(:,1)=zeros(1,WGNum);
array(ceil(WGNum/2),1)=1;  % 态矢
for t=2:options.EvalutionNum
    array(:,t)=(H*dL+2i*eye(WGNum))/(-H*dL+2i*eye(WGNum))*array(:,t-1);
end
I = mat2gray(array.*conj(array));
end