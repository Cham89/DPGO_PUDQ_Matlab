function R = boumal_check_multirobot2(local_graph_r, varargin)
    p = inputParser;
    addParameter(p, 'num_dirs', 3);      
    addParameter(p, 't_min_exp', -8);    
    addParameter(p, 't_max_exp',  0);     
    addParameter(p, 'N_t', 400);            
    addParameter(p, 'make_plot', true);
    parse(p, varargin{:});
    num_dirs   = p.Results.num_dirs;
    t_min_exp  = p.Results.t_min_exp;
    t_max_exp  = p.Results.t_max_exp;
    N_t        = p.Results.N_t;
    make_plot  = p.Results.make_plot;

    num_vertices = length(local_graph_r.vertices);
    anchor_first = isfield(local_graph_r,'anchor_first') && local_graph_r.anchor_first && isfield(local_graph_r,'robot_id') && local_graph_r.robot_id == 1;

    if anchor_first
        movable = 5:(4*num_vertices);
    else
        movable = 1:(4*num_vertices);
    end

    X = G_get_X(local_graph_r);
    [rgrad, ~, F0, P_X] = rgn_gradhess_multi_J(local_graph_r); 
    P = P_X; 
    t = logspace(t_min_exp, t_max_exp, N_t);

    slopes = zeros(1, num_dirs);
    dir_ok = true(1, num_dirs);

    if make_plot, clf; end

    for d = 1:num_dirs
        v = zeros(4*num_vertices,1);
        v(movable) = randn(length(movable),1);
        v = P * v;
        nv = norm(v);
        if nv == 0
            dir_ok(d) = false; slopes(d) = NaN; continue;
        end
        v = v / nv;

        g = P * rgrad;
        gdotv = g' * v;

        Eg = zeros(1, N_t);
        for j = 1:N_t
            if exist('Exp_X_N','file') == 2
                Xnew = Exp_X_N(X, t(j) * v);
            else
                Xnew = X + t(j) * v;
            end
            G2 = G_set_X(local_graph_r, Xnew);
            [~, ~, Fnew] = rgn_gradhess_multi_J(G2);

            Eg(j) = Fnew - F0 - t(j) * gdotv;
        end

        mask = isfinite(Eg) & Eg > 0 & isfinite(t);
        if nnz(mask) < 20
            dir_ok(d) = false;
            slopes(d) = NaN;
            continue;
        end

        tt  = t(mask);
        EE  = Eg(mask);
        i1  = max(1, floor(numel(tt)/3));
        i2  = min(numel(tt), floor(2*numel(tt)/3));
        tt  = tt(i1:i2); EE = EE(i1:i2);

        pfit = polyfit(log10(tt), log10(EE), 1);
        slopes(d) = pfit(1);

        if make_plot
            hold on;
            set(gca,'xscale','log','yscale','log');
            plot(t, max(Eg, eps), 'DisplayName', sprintf('dir %d (slope %.2f)', d, slopes(d)));
        end
    end

    if make_plot
        yref = (t.^2) * 1;  
        plot(t, yref, 'k--', 'DisplayName', 'slope 2 ref');
        grid on; box on;
        title(sprintf('Boumal check — Robot %d (%d verts)', local_graph_r.robot_id, num_vertices));
        xlabel('t'); ylabel('E_g(t)');
        legend('Location','southwest');  
    end

    good  = dir_ok & isfinite(slopes);
    med_slope = median(slopes(good));
    mean_slope = mean(slopes(good));

    fprintf('\n=== BOUMAL CHECK for Robot %d ===\n', local_graph_r.robot_id);
    fprintf('F(x) = %.6e, tested %d/%d directions\n', F0, nnz(good), num_dirs);
    for d = 1:num_dirs
        if good(d)
            fprintf('  dir %d: slope = %.3f\n', d, slopes(d));
        else
            fprintf('  dir %d: (invalid)\n', d);
        end
    end
    fprintf('  median slope = %.3f, mean slope = %.3f\n', med_slope, mean_slope);

    if ~isempty(med_slope) && abs(med_slope - 2.0) < 0.35
        fprintf('Gradient consistent with O(t^2) (Boumal)\n');
    else
        fprintf('Bla there is a bug');
    end

    R = struct('slopes',slopes, 'median_slope',med_slope, 'mean_slope',mean_slope);
end
