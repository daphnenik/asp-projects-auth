function h = GiannakisFormula(cum3,q,L)
    for k = 0:q
    h(k+1) = cum3(q+L+1,k+L+1)/cum3(q+L+1,L+1);
    end
end