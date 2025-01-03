function Az = WGA_evaluate_expm_func(options,zn)
% WGA_evaluate_expm_func(options,zn) using in <lsqcurvefit>
%   options = [dl,kappa,EdgeKappa,beta,EdgeBeta,WGNum]
%   zn = [z;n]
%   z - evaluate length
%   n - num of WG
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[Az] = WGA_evaluate_expm_func(options,zn);
%
% Params:
%   - options   [required]  [size=1,6] [dl,kappa,EdgeKappa,beta,EdgeBeta,WGNum]
%   - zn        [required]  [size=:,2] 演化长度；波导序号
%
% Return:
%   - Az 在z处对应的几率幅
%
% Matlab Version: R2024b
%
% Author: oyy
options(6) = round(options(6));
H = diag(ones(1,options(6))*options(4))+diag(ones(1,options(6)-1)*options(2),-1)+diag(ones(1,options(6)-1)*options(2),1);
if options(3)
[H(1,2),H(2,1),H(end-1,end),H(end,end-1)] = deal(options(3));
end
[H(1,1),H(end,end)] = deal(options(5));
array=zeros(options(6),1);
array(ceil(options(6)/2))=1;  % 态矢
z = unique(zn(:,1),"stable");
Az = nan(size(zn,1),1);
for temp = 1:length(z)
Az(temp*options(6)-options(6)+1:temp*options(6)) = abs(expm(-1i*(z(temp)+options(1)).*H)).^2*array;
end
end