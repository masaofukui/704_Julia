
function plot_policy_functions(param, Bellman_result;a_idx=1:param.Na)
    @unpack ag,amin = param
    @unpack c_pol, a_pol = Bellman_result
    # Set the default font to Computer Modern
    default(fontfamily="Computer Modern")

    p_c = plot(ag[a_idx], c_pol[a_idx,[2,1]], 
        lw=5,
        label=["high income" "low income" ],
        linestyle=[:solid :dash])
    title!(L"Consumption Policy Functions, $c(a,y)$")
    xlabel!(L"Asset, $a$")

    p_a = plot(ag[a_idx], a_pol[a_idx,[2,1]], 
        lw=5,
        label=["high income" "low income" ],
        linestyle=[:solid :dash])
    plot!(ag[a_idx], ag[a_idx], 
        lw=2,
        label="45 degree line",
        linestyle=:dash,
        color=:grey)
    title!(L"Asset Policy Functions, $a'(a,y)$")
    xlabel!(L"Asset, $a$")

    p_pol = plot(p_c, p_a, layout=(1,2), size=(900, 500))
    display(p_pol)

    return p_pol
end