function T = SE2_Exp_1(X)
   v = X(1:2, 3);
   alpha = X(2, 1);
   
   t = V_alpha(alpha)*v;
   T = [R_alpha(alpha), t; 0, 0 1];
end
