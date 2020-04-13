function [S, A_est] = nnschalfprior3( X, A, Z, p, fname )
% nnschalfprior - Do hyperspectral unmixing using L1/2 method and Z regular
%
% SYNTAX:
% nnschalfprior( X, A p, fname )
%
% INPUT:    
% X       the non-negative data (in columns)
% A       initial matrix of A
% fname   name of file to write
% p       algorithm options:
%
% p.sources      number of components
% p.derta        value of derta
% p.lambda       value of lambda

dims = size(X,1);
samples = size(X,2);
sources = p.sources;
derta = p.derta;
lambda = p.lambda;
mu = p.mu;

times = 0;
%rand('state',sum(100*clock));
% Initializing A
%A = abs(randn(dims,sources));
%A = A./(ones(dims,1)*sqrt(sum(A.^2)));
A_est = A;

% Initializing S
S = p.S;
%S = hyperFcls(X,A);
%  S = (pinv(A_est' * A_est) * A_est' * X);
%  S(S<0) = 0;
%S = S./(ones(size(S,1),1)*sum(S));
%     S = abs(randn(sources,samples));
%     S = S./(ones(size(S,1),1)*sum(S));
% S = (pinv(A_est) * X);
% S(S<0) = 0;
% S = S./(ones(size(S,1),1)*sum(S));

% S = max(pinv(A_est' * A_est) * A_est' * X, 1e-10);
% S = S./(ones(size(S,1),1)*sum(S));

% X(X < 0) = 1e-19;
% S(S < 0) = 1e-19;
% A_est(A_est < 0) = 1e-19;
  

% These will store the history of the objective function
objhistory = [];
iterhistory = [];

% Loop indefinitely
iter = 0;
tt = 0;


% obj = 0.5 * norm((X-A_est*S),'fro') + lambda * sum(sum(sqrt(S)));
% obj = 0.5 * norm((X-A_est*S),'fro')^2 + lambda * sum(sum(S));

%     Xf = [X;derta * ones([1 samples])];
%     Af = [A_est;derta * ones([1 sources])]; 
    
obj = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S))) + mu * sum(sum((S*Z).^2));
% obj = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(S)) + mu * sum(sum((S*Z).^2));

% obj = 0.5 * norm((X-A_est*S),'fro') + lambda * norm(S,1);
% obj = 0.5 * norm((X-A_est*S),'fro');
%obj = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S)));

zz = Z*Z';
while 1,
%     if rem(iter,100)==0,
%         obj = 0.5*sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S)));
%     end
  % Save basis (only every 10th iteration)
       % Save objective function evolution 
        objhistory = [objhistory obj];
      iterhistory = [iterhistory iter];
   if rem(iter,100)==0,
      % Update activity measures
      activations = sqrt(mean(S'.^2));
%      fprintf(['\nSaving file: ' fname '...']);
%      save(fname,'p','A_est','S','objhistory','iterhistory','activations');
%      fprintf('DONE!\n');
   end  

         
  %Augment X & A to Xf & Af 
 

%  Xf = [X];
%  Af = [A_est];  

  % Update S and A with multiplicative steps
   %tmp11 = S.*(Af'*Xf)./(Af'*Af*S);
     %[m n] = find(S>1e-4);
%      tmp12 = S.*(Af'*Xf)./((Af'*Af*S) + lambda/2 ./ sqrt(S) + 1e-9);
%      tmp11 = S.*(Af'*Xf)./((Af'*Af*S) + 1e-9);
%      %S = S.*(Af'*Xf)./((Af'*Af*S));
%  %    S = S./(ones(size(S,1),1)*sum(S));
%      %S = S.*(Af'*Xf)./((Af'*Af*S) + lambda);
%     index1 = tmp11 < 1e-4;
%     index2 = tmp12 >= 1e-4;
%     tmp1(index1) = tmp11(index1);
%     tmp1(index2) = tmp12(index2);
%     tmp1 = reshape(tmp1,sources,samples);
%     S = tmp1;
    %S = S.*(Af'*Xf)./(Af'*Af*S + lambda/2 ./ sqrt(S));

    
    Xf = [X;derta * ones([1 samples])];
    Af = [A_est;derta * ones([1 sources])]; 
    
    S = S.*(Af'*Xf)./((Af'*Af*S) + 0.5 * lambda ./ sqrt(S) + 1e-19 + 2*mu * S*zz);
%     S = S.*(Af'*Xf)./((Af'*Af*S) + lambda + 1e-19 + 2*mu * S*zz);
    %S = S.*(Af'*Xf)./((Af'*Af*S));
    A_est = A_est.*(X*S')./(A_est*S*S' );
%     scaling = sqrt(sum(A_est.^2));  
%     A_est = A_est./(ones(size(A_est,1),1)*scaling);
%     S = S.*(scaling'*ones(1,size(S,2))); 
    
  %obj2 = sum(sum((X-A_est*S).^2))/samples + lambda * sum(sum(sqrt(S)))/samples;
    %S = S_tmp;
    %A_est = A_est_tmp;
  % Caculate new object function
 

     %objnew = 0.5 * norm((X-A_est*S),'fro') + lambda * sum(sum(sqrt(S)));
     %objnew = 0.5 * norm((X-A_est*S),'fro')^2 + lambda * sum(sum(S));
     
     objnew = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S))) +  mu * sum(sum((S*Z).^2));
%      objnew = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(S)) +  mu * sum(sum((S*Z).^2));
     
     %objnew = 0.5 * norm((X-A_est*S),'fro') + lambda * norm(S,1);
    %objnew = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S)));%0.5 * norm((X-A_est*S),'fro');
    %objnew = 0.5 * norm((X-A_est*S),'fro') + lambda * sum(sum(S));
    
    toltemp = abs(obj - objnew)/obj;
    if toltemp < 0.0001
        times = times + 1;
    end
    
    obj = objnew;
    error(iter + 1) = obj;
  if rem(iter,100)==0,
    display(['iter: ', num2str(iter), ' error: ', ...
    num2str(obj), ' tol: ', num2str(toltemp)]);
  end
  
  % Normalize columns of A (and scale rows of S correspondingly)
%   scaling = sqrt(sum(A.^2));  
%   A = A./(ones(size(A,1),1)*scaling);
%   S = S.*(scaling'*ones(1,size(S,2)));      
  
  % Update iteration counter
      %  error(iter) = 0.5 * sum(sum((X-A_est*S).^2))/samples + lambda * sum(sum(sqrt(S)))/samples;
%      errorA(iter) = norm(A(:,2) - Known);
  iter = iter+1;
  %if iter > 3000 || diff < 0.00001
%   if diff < 0.00001
%     times = times + 1;
%   end
  if iter > 3000 || times>=1%%
     %figure;plot(error);
    display(['finnay iter: ', num2str(iter), ' error: ', ...
     num2str(obj), ' tol: ', num2str(toltemp)]);
%      figure;plot(errorA);
%      fprintf(['\nSaving file: ' fname '...']);
%      %save(fname,'p','A_est','S','objhistory','iterhistory','activations');
%      save(fname,'p','A_est','S','objhistory','iterhistory');
%      fprintf('DONE!\n');
    break;
  end
end

