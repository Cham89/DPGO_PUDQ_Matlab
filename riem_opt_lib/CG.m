%Conjuage Gradient method for solving Hx=b
function s = CG(H, b)
    N_max = 2*length(b);
    
    %Initialize problem variables
    v_n = zeros(length(b),1);
    r_n_minus = b;
    p_n = r_n_minus;
    
    s = v_n;
    if r_n_minus == 0
        return;
    end
    
    for n=1:N_max
        H_p = H*p_n;
        pHp = p_n'*H_p;
        
        alpha_n = (r_n_minus'*r_n_minus)/pHp;
        v_n = v_n + alpha_n * p_n;
        r_n = r_n_minus - alpha_n * H_p;
        
        if norm(r_n) < 1e-10
            s = v_n;
            disp(['CG converged in ' num2str(n) ' iterations.'])
            return;
        end
        
        beta_n = (r_n'*r_n)/(r_n_minus'*r_n_minus);
        p_n = r_n + beta_n*p_n;
        
        r_n_minus = r_n;
    end
    
    s = v_n;
    disp(['CG terminated after ' num2str(n) ' iterations']);
end