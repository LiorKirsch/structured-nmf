function [W,H,diff_record,time_record]=nmf_prob(parms,X,W_init,H_init)
%
% Probabilistic NFM interpretating X as samples from a multinomial
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Lars Kai Hansen, IMM-DTU (c) November 2005
%
maxiter = take_from_struct(parms, 'maxiter', 1000);
loglevel = take_from_struct(parms, 'loglevel', 1);
print_interval = take_from_struct(parms, 'print_interval', 100);
enforce_W_prob_norm = take_from_struct(parms, 'enforce_W_prob_norm', false);
record_scores = take_from_struct(parms, 'record_scores', false);


[N,M]=size(X);
[N,K]=size(W_init);
powers=1.5+(2.5-1.5)*((1:maxiter)-1)/(maxiter-1);
% INITIALIZE
[D,N]=size(X);
X_factor = (sum(sum(X)));
X_org = X;
X=X/X_factor;

W=W_init;
W=W./repmat(sum(W,1),D,1);
H=H_init;
H=H./repmat(sum(H,2),1,N);
P=ones(K,1);
P=P/sum(P);
W1=W;H1=H;

% use W*H to test for convergence
Xr_old = W*H;

diff_record =nan(1,maxiter);
time_record =nan(1,maxiter); tic;
for iter=1:maxiter,
    %E-step
    Qnorm=(W*diag(P))*H;

    for k=1:K,
        %E-step
        Q=(W(:,k)*(H(k,:)*P(k)))./(Qnorm+eps);
        XQ=X.*Q;
        %M-step W
        dummy=sum(XQ,2);
        W1(:,k)=dummy/(sum(dummy));
        dummy=sum(XQ,1);
        H1(k,:)=dummy/(sum(dummy));
    end

    W=W1;
    H=H1;

    if enforce_W_prob_norm
        W=W./(repmat(sum(W,2),1,size(W,2))+eps); % normalize rows to unit length
    end



    %%%%%%%%%%%%%%%%%%%%%%%
    % print to screen
    %%%%%%%%%%%%%%%%%%%%%%%
    if (rem(iter,print_interval)==0) & loglevel >0,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X_org,W*diag(sqrt(P))*X_factor*diag(sqrt(P))*H);
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
    if record_scores
        diff_record(iter) = nmf_euclidean_dist(X_org,W*diag(sqrt(P))*X_factor*diag(sqrt(P))*H);
        time_record(iter) = toc;
    end
end,

W=W*diag(sqrt(P))*X_factor;
H=diag(sqrt(P))*H;