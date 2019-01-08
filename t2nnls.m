function [T2Spectra, k, sfit] = t2nnls(signal, te, doPlot, T2Range, ...
                                       NT2Samples, regularization)
%
% T2 non-negative least squares algorithm for NMR relaxation  
% based on the reference: 
% Whittall and Mackay, J Magn Reson 84, 134-152 (1989)
% this implementation written by Galen Reed, 01/07/2019



T2 = logspace(log10(T2Range(1)), log10(T2Range(2)), NT2Samples);

T2DecayMatrix = exp( -te(:) * (1./T2(:).') );

%cnum = cond(T2DecayMatrix);
%disp(cnum);


if (regularization == 1)% not tested in forever

    lambdamax = find_lambdamax_l1_ls_nonneg(T2DecayMatrix', s);
    lambda = .001 * lambdamax;
    [k,status] = l1_ls_nonneg(T2DecayMatrix, signal, lambda, rel_tol,quiet);

elseif(regularization == 2)% broken
    
    lambdamax = norm(2*(T2DecayMatrix'*signal), 2);
    lambda = 5e-5* lambdamax;
 
    li = lambda * eye(NT2Samples);
    zp = zeros([NT2Samples 1]);
    size(li)
    size(NT2Samples)
    
    Acat = vertcat(NT2Samples, li);
    
    ycat = vertcat(s, zp);

    
  
    k = lsqnonneg(Acat, ycat);
    k = k(1:NT2);
    
    
elseif(regularization == 0) % works great
    
    k = lsqnonneg(T2DecayMatrix, signal);
     
end


sfit = T2DecayMatrix*k;
T2Spectra = k;

if(doPlot == 1)

    fontsize = 20;
    figure();
    subplot(2, 1, 1)
    semilogy(te, sfit, '-', te, signal, 'x');
    grid on;
    xlim([0 te(end)]);
    xlabel('TE [ms]', 'fontsize', fontsize);
    ylabel('signal', 'fontsize', fontsize);
    
    
    subplot(2, 1, 2);
    semilogx(T2, T2Spectra)
    grid on;
    xlim([T2(1) T2(end)]);
    xlabel('T_2 [ms]', 'fontsize', fontsize);
    ylabel('signal', 'fontsize', fontsize);
end





