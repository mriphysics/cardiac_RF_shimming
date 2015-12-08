%  Constrained optimisation computation for finding optimal shims within
%  constraints for in vivo experiments run at fixed pulse amplitude and for
%  fixed global SAR constraint.
%
% Created by Arian Beqiri, King's College London, December 2015.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [all_shims,Error,LSAR,GSAR] = CVX_b1_shim_LCurve(A,Amask,VOP,QG,targmask,idx,drmax,crs,Nc)

% Compute initial quad mode solution for starting MLS phase
mls_tmpsoln = sum(A,2);

% Preallocate variables
Error = zeros(length(crs),1); LSAR=Error; GSAR=Error;
all_shims = zeros(length(crs),Nc);

% Set initial target phase to that of quadrature map
z = exp(1i*(angle(mls_tmpsoln(idx))));

%%%%%%%%%%% MLS Iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for ii=1:length(crs)

    lsarlim = crs(ii);   

    mlsinit = ones(Nc,1);
    mls_shim = zeros(Nc,1);
    counter = 0;

    while (abs(norm(mls_shim - mlsinit)) > 0.1)

        mlsinit = mls_shim;

        %%%% CVX Optimisation %%%%%%%%%%%
        cvx_begin quiet
             variable y(Nc) complex;
             minimise(norm(Amask*y - targmask.*z));
             subject to           
                (y)'*QG*(y) <= 0.1;    % Global SAR Constraint
                SAR(y,VOP) <= (lsarlim + 0i);    % Local SAR Constraint
                max(abs(y)) <= drmax;
        cvx_end

        mls_shim = y;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mls_tmpsoln = A*mls_shim;  
        z = exp(1i*(angle(mls_tmpsoln(idx))));

        err = norm((abs(mls_tmpsoln(idx))-targmask))/norm(targmask);
        conv = abs(norm(mls_shim - mlsinit));

        fprintf('Convergence = %1.6f\tError = %1.4f\n',conv,err);
        counter = counter+1;
        
        % Solutions do not improve much for L-cuve solutions of ii>1 
        % after 1st MLS iteration so force only one iteration for these
        if (ii>1) && (counter>0);break;end

    end
    
    all_shims(ii,:) = mls_shim;
    Error(ii) = err;
    disp(['Local SAR constraint = ' num2str(crs(ii))])
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Local and Global SAR values from shims
parfor jj=1:length(crs)
    shim = squeeze(all_shims(jj,:)).';
    LSAR(jj) = abs(SAR_par(VOP,shim));
    GSAR(jj) = abs(shim'*QG*shim);
end

end