
function construct_IRF(xfull,n)
    d = Dict{String,Any}()
    for i = 1:n
        d[(allinputs[i])] = xfull[:,i]
    end
    return d
end

function qfun(P,theta)
    eval = P["m"].*theta.^(-P["alph"])
end
function ffun(P,theta)
    eval = P["m"].*theta.^(1-P["alph"])
end