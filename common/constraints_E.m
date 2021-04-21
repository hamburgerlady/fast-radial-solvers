function eqs = constraints_E(E)

eqs = 2*(E*E.')*E - sum(diag(E*E.'))*E;
eqs = [eqs(:);det(E)];

