function sysmodel = SparseRegression(Y,Yp,U,dt, options_method, lambda)

Nstates = size(Y,1);
Ninputs = size(U,1);

% Sparse regression
switch options_method.sparsify
    case 'LASSO'
        % using the l1 norm to promote sparsity
        G = zeros(Nstates+Ninputs,Nstates);
        for i = 1:Nstates
            [G(:,i), FitInfo] = lasso([Y;U]',Yp(i,:)','Lambda',lambda); %,'Lambda',0.00000001
            %     lassoPlot(G{i},FitInfo);
        end
    case 'ILSTH'
        G = sparsifyDynamics([Y;U]',Yp',lambda,Nstates);
end

G = G';
A = G(1:Nstates,1:Nstates);
B = G(1:Nstates,Nstates+1:Nstates+Ninputs);
C = eye(Nstates,Nstates);
D = zeros(Nstates,1);
sysmodel = ss(A,B,C,D,dt);
end