function [W H ]=Lee_SNMF(V,r,opts,lambda)

% for L1/2 sparse， ref lambda /5 and ref delta = 5; 
defopts=struct('maxiter',1000,'tol',1e-6,'trackit',20,'w',[],'h',[],'Vtrue',[]);
if ~exist('opts','var')
    opts=struct;
end
[MaxIt, tol, trackit,W,H,Vtrue]=scanparam(defopts,opts);
if isempty(Vtrue)
    Vtrue=V;
end
ts=clock;
eps=1e-9;

[M L]=size(V);
V=max(V,0);

  
%% Initialization;
if isempty(W)
    rand('seed',243)
%       rand('seed',24)
    W=rand(M,r);
end
if isempty(H)
    rand('seed',76)
%     rand('seed',156)
    H=rand(r,L);
end
AS = sum(H);
% H = H./repmat(AS,size(H,1),1);

%%%%%%%%%%%% sparsity
attend = zeros(1,M);
for i = 1:M
    attend(i)= (sqrt(L)- norm(V(i,:),1)/norm(V(i,:),2))/(sqrt(L-1));
end
% lambda = sum(attend)/sqrt(M);
% lambda = sum(attend)/sqrt(M)/1.5;
%% SECTION TITLE
% DESCRIPTIVE TEXT
%%%%%%%%%%%%sum to 1 or not
delta = 1.5e1; Vo = [V;delta*ones(1,L)];
%%%%%%%%%%%
% I=eye(r,r);
obj=inf;
objH = inf;
for it=1:MaxIt
    
    W0=W;
    H0 = H;
%    Cold = operatorA(W,H,V);
    W=W.*((V*H')./max(W*H*H',eps));
  
    WE = [W;delta*ones(1,r)];
    idx = find(H<=1e-4); PG = (max(H,eps)).^(-0.5); PG(idx) = eps;
    H=H.*((WE'*Vo)./max((WE'*WE)*H + 0.5*lambda*PG,eps)); %L1/2 Sparse
%    H=H.*((WE'*Vo)./max((WE'*WE)*H + lambda*(abs(H.^(-1))),eps)); %L1 Sparse
%    H=H.*((W'*V)./max((W'*W)*H,eps)); %L1 Sparse
%     H=H.*((WE'*Vo)./max((WE'*WE)*H,eps)); %L1 Sparse with sum to 1 constraint

% 
    if (trackit>0)&&(~rem(it,trackit))&&(it>MaxIt*0.2)  %mean trackit >0, it是其整数倍，并且大于0.2 
        obj(it)=norm(W-W0,'fro');  
%         objH = norm(H-H0,'fro');
        if obj(it)<tol %&&  objH<tol
            break;
        end    
    end
end


