%Truncated conjugate gradient method for solving min 0.5<s, Hs> - <b, s>
function s = tCG(H, b, Delta, kappa, theta, P_x, S_prec)
    N_max = 2*length(b);
    
    % The following recurrences for Prec-based norms and inner
    % products come from [CGT2000], pg. 205, first edition.
    % Below, P is the preconditioner.
    %
    % <eta_k,P*delta_k> = 
    %          beta_k-1 * ( <eta_k-1,P*delta_k-1> + alpha_k-1 |delta_k-1|^2_P )
    % |delta_k|^2_P = <r_k,z_k> + beta_k-1^2 |delta_k-1|^2_P
    %
    % Therefore, we need to keep track of:
    % 1)   |delta_k|^2_P
    % 2)   <eta_k,P*delta_k> = <eta_k,delta_k>_P
    % 3)   |eta_k  |^2_P
    %
    % Initial values are given by
    %    |delta_0|_P = <r,z>
    %    |eta_0|_P   = 0
    %    <eta_0,delta_0>_P = 0
    % because we take eta_0 = 0 (if useRand = false).
    
    %Initialize problem variables
    v_n = zeros(length(b), 1);
    r_n_minus = b;
    p_n = r_n_minus;

    %Precompute r_0 threshold term
    r_thresh = norm(r_n_minus)*min(norm(r_n_minus)^theta, kappa)

    s = v_n;
    if r_n_minus == 0
        return;
    end

    for n=1:N_max
        %Ensure p_n is in TxM
        p_n = P_x*p_n;
        
        H_p = H*p_n;
        pHp = p_n'*H_p;
        alpha_n = (r_n_minus'*r_n_minus) / pHp;
        v_n_plus = v_n + alpha_n * p_n;
        
        php_neg = pHp <= 0;
        tr_exceeded = norm(S_prec*v_n_plus) >= Delta;
        
        if php_neg
            disp('tCG pHp <= 0!')
        end
        
        if tr_exceeded
            disp(['tCG Trust region exceeded (Delta=' num2str(Delta) ')'])
        end
        
        if php_neg || tr_exceeded
            a = p_n'*p_n;
            b = 2*v_n'*p_n;
            c = norm(S_prec*v_n)^2 - Delta^2;
            
            %Set t equal to the + root of the quadratic
            t = max(roots([a b c]));
            v_n = v_n + t*p_n;

            s = v_n;
            
            disp(['tCG terminated in ' num2str(n) ' iterations with truncated solution.'])
            return;
        end

        v_n = v_n_plus;
        r_n = r_n_minus - alpha_n*H*p_n;
        
        if norm(r_n) <= r_thresh
            s = v_n;
            disp(['tCG converged in ' num2str(n) ' iterations.'])
            
            norm(r_n)
            return;
        end

        beta_n = (r_n'*r_n)/(r_n_minus'*r_n_minus);
        p_n = r_n + beta_n*p_n;

        r_n_minus = r_n;        
    end
    
    s = v_n;
    disp(['tCG terminated after ' num2str(n) ' iterations']);
end




