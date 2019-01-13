function [T2Spectra, k, sfit] = t2nnls(signal, te, doPlot, T2Range, ...
                                       NT2Samples, regularizationType, regularizationLambda)
%
% T2 non-negative least squares algorithm for NMR relaxation  
% based on the reference: 
% Whittall and Mackay, J Magn Reson 84, 134-152 (1989)
% this implementation written by Galen Reed, 01/07/2019



T2 = logspace(log10(T2Range(1)), log10(T2Range(2)), NT2Samples);

T2DecayMatrix = exp( -te(:) * (1./T2(:).') );

%cnum = cond(T2DecayMatrix);
%disp(cnum);



lambda = regularizationLambda;

if (regularizationType == 1)% not tested in forever
  
    [k,status] = l1_ls_nonneg(T2DecayMatrix, signal, lambda, rel_tol,quiet);

elseif(regularizationType == 2)% works great
   
    lambdaDiag = lambda * eye(NT2Samples);
    zeroPad = zeros([NT2Samples 1]);
    Acat = vertcat(T2DecayMatrix, lambdaDiag);
    ycat = vertcat(signal, zeroPad);
    k = lsqnonneg(Acat, ycat);    
    
elseif(regularizationType == 0) % works great
    
    k = lsqnonneg(T2DecayMatrix, signal);
     
end


sfit = T2DecayMatrix*k;
T2Spectra = k;

if(doPlot == 1)

    fontsize = 20;

    figure();
    semilogy(te, sfit, '-', te, signal, 'x');
    grid on;
    xlim([0 te(end)]);
    xlabel('TE [s]', 'fontsize', fontsize);
    ylabel('signal', 'fontsize', fontsize);
    
    
    figure();
    semilogx(T2, T2Spectra)
    grid on;
    xlim([T2(1) T2(end)]);
    xlabel('T_2 [s]', 'fontsize', fontsize);
    ylabel('signal', 'fontsize', fontsize);
end





