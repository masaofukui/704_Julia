
function wrapper_PE_policy_plot(param,param_changed; fig_save = 0, fig_name = "")
    @unpack beta,r,phi = param_changed
    if beta*(1+r) >= 1
        println("beta*(1+r) >= 1")
    end
    @assert param.amin - phi > - minimum(param.yg)./r || r <= 0
    Bellman_result = solve_policy_EGM(param, param_changed = param_changed)
    p_pol = plot_policy_functions(param, Bellman_result,param_changed;a_idx=1:50)
    
    return p_pol
end

function plot_policy_functions(param, Bellman_result,param_changed;a_idx=1:param.Na)
    @unpack ag,amin = param
    @unpack c_pol, a_pol = Bellman_result
    @unpack phi = param_changed
    ag_plot = ag .- phi;
    # Set the default font to Computer Modern
    default(fontfamily="Computer Modern",
    xgrid = :none)

    p_c = plot(ag_plot[a_idx], c_pol[a_idx,[2,1]], 
        lw=5,
        label=["high income" "low income" ],
        linestyle=[:solid :dash])
    title!(L"Consumption Policy Functions, $c(a,y)$")
    xlabel!(L"Asset, $a$")

    p_a = plot(ag_plot[a_idx], a_pol[a_idx,[2,1]] .-phi, 
        lw=5,
        label=["high income" "low income" ],
        linestyle=[:solid :dash])
    plot!(ag_plot[a_idx], ag_plot[a_idx], 
        lw=2,
        label="45 degree line",
        linestyle=:dash,
        color=:grey)
    title!(L"Asset Policy Functions, $a'(a,y)$")
    xlabel!(L"Asset, $a$")

    p_pol = plot(p_c, p_a, layout=(1,2), size=(900, 500))
    display(p_pol)
    if fig_save == 1
        savefig(p_pol, fig_dir*"policy_functions_"*fig_name*".pdf")
    end


    p_a_only = plot(ag_plot[a_idx], a_pol[a_idx,1] .-phi, 
    lw=5,
    label="low income",
    linestyle=:dash,
    color=default_colors[2])
    plot!(ag_plot[a_idx], ag_plot[a_idx], 
        lw=2,
        label="45 degree line",
        linestyle=:dash,
        color=:grey)
    title!(L"Asset Policy Functions, $a'(a,y)$")
    xlabel!(L"Asset, $a$")

    if fig_save == 1
        savefig(p_a_only, fig_dir*"policy_functions_a_only_"*fig_name*".pdf")
    end

    return p_pol
end


function plot_distribution(param, Bellman_result, ss_distribution; fig_save = 0, fig_name = "")
    @unpack ag, Ny, Na = param
    @unpack phi = param_changed
    @unpack c_pol = Bellman_result
    ag_plot = ag .- phi;

    default(; titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern", 
    legendfontfamily = "Computer Modern",
    titlefontsize = 20,
    xguidefontsize = 12,
    legendfontsize = 12,
    yguidefontsize = 12,
    xgrid = :none)



    c_distribution = copy(ss_distribution)
    for iy in 1:Ny
        dc = diff(c_pol[:,iy])
        dc = [dc[1]; dc]
        c_distribution[:,iy] = ss_distribution[:,iy]./dc
    end


    p1 = plot(ag_plot, c_pol[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],label=["high income" "low income"])
    title!("Consumption policy, \$c(a,y)\$")
    xlabel!("Asset, \$a\$")

    #plot wealth distribution 
    p2 = plot(ag_plot, ss_distribution[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],label=["high income" "low income"])
    title!("Wealth distribution, \$g(a)\$")
    xlabel!("Asset, \$a\$")


    p3 = plot( c_pol[:,[2,1]], c_distribution[:,[2,1]], 
        lw=5,linestyle = :solid)
        plot!( c_pol[:,[2,1]], c_distribution[:,[2,1]], 
        lw=5,linestyle = :dash)
    title!("Consumption Distribution")

    p_all = plot(p2,p1,layout = (2,1),size = (1000,700))
    if fig_save == 1
        savefig(p_all, fig_dir*"distribution_"*fig_name*".pdf")
    end
    return p_all
end


function plot_distribution_compare(param, Bellman_list, ss_distribution_list,var_input_vec; 
    fig_save = 0, fig_name = "", var_change = "phi")
    @unpack ag, Ny, Na = param
    @unpack c_pol = Bellman_result
    if var_change == "phi"
        ag_plot1 = ag .- var_input_vec[1];
        ag_plot2 = ag .- var_input_vec[2];
    else
        ag_plot1 = ag .- param.phi;
        ag_plot2 = ag .- param.phi;
    end

    c_pol1 = Bellman_list[1].c_pol;
    c_pol2 = Bellman_list[2].c_pol;
    ss_distribution1 = ss_distribution_list[1];
    ss_distribution2 = ss_distribution_list[2];


    default(; titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern", 
    legendfontfamily = "Computer Modern",
    titlefontsize = 20,
    xguidefontsize = 12,
    legendfontsize = 12,
    yguidefontsize = 12,
    xgrid = :none)


    c_distribution1 = copy(ss_distribution1)
    c_distribution2 = copy(ss_distribution2)
    for iy in 1:Ny
        dc1 = diff(c_pol1[:,iy])
        dc1 = [dc1[1]; dc1]
        c_distribution1[:,iy] = ss_distribution1[:,iy]./dc1
        dc2 = diff(c_pol2[:,iy])
        dc2 = [dc2[1]; dc2]
        c_distribution2[:,iy] = ss_distribution2[:,iy]./dc2
    end



    p1 = plot(ag_plot1, c_pol1[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],
        label=:none,
        alpha = 0.2,
        linecolor=[default_colors[1] default_colors[2]])
    plot!(ag_plot2, c_pol2[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],
        label=["high income" "low income"],
        linecolor=[default_colors[1] default_colors[2]])
    title!("Consumption policy, \$c(a,y)\$")
    xlabel!("Asset, \$a\$")

    #plot wealth distribution 
    p2 = plot(ag_plot1, ss_distribution1[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],
        label=:none,
        alpha = 0.2,
        linecolor=[default_colors[1] default_colors[2]])
    plot!(ag_plot2, ss_distribution2[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],
        label=["high income" "low income"],
        alpha = 1.0,
        linecolor=[default_colors[1] default_colors[2]])
    title!("Wealth distribution, \$g(a)\$")
    xlabel!("Asset, \$a\$")

    p3 = plot(c_pol1[:,[2,1]], c_distribution1[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],label=["high income" "low income"],alpha=0.1,color=default_colors[1:2])
    plot!(c_pol2[:,[2,1]], c_distribution2[:,[2,1]], 
        lw=5,
        linestyle=[:solid :dash],label=["high income" "low income"],alpha=1.0,color=default_colors[1:2])
    title!("Consumption distribution")

    p_all = plot(p2,p1,layout = (2,1),size = (1000,700))
    if fig_save == 1
        savefig(p_all, fig_dir*"compare_distribution_"*fig_name*".pdf")
    end
    return p_all

end



function plot_asset_demand_function(param, r_vec, A_vec; compare = false, A_vec_old = nothing)
    A_supply = zeros(length(r_vec))
    plt_asset = plot(r_vec,A_vec,
        lw=5,
        label="Asset Demand, \$A(r)\$",
        xlabel="Interest rate, \$r\$",
        color=default_colors[1]
    )
    if compare
        plot!(r_vec,A_vec_old,
            lw=5,
            label=:none,
            color=default_colors[1],
            alpha = 0.2
        )
    end
    plot!(r_vec,A_supply,
        lw=5,
        label="Asset Supply",
        color=default_colors[2],
        linestyle=:dash
    )
    vline!([1/param.beta-1], 
        label=:none,
        lw=3,
        linestyle=:dot,
        color=default_colors[3]
    )
    if !compare
        annotate!([(1/param.beta-1, minimum(A_vec)+0.02, 
            text("\$1/\\beta - 1\$", :black, :right, 10))])
    end
    display(plt_asset)
    if fig_save == 1
        if compare
            savefig(plt_asset, fig_dir*"asset_demand_function_compare.pdf")
        else
            savefig(plt_asset, fig_dir*"asset_demand_function.pdf")
        end
    end
end


function plot_power_law(param,dist_a,ag; fig_save = 0, fig_name = "",title_label ="",low_q=0.1,high_q=0.001)
    @unpack yg, Ny = param
    empirical_cdf_a = cumsum(dist_a)
    empirical_cdf_a = empirical_cdf_a./empirical_cdf_a[end]
    ccdf_a = 1.0.-empirical_cdf_a
    idx_a = ccdf_a .<= low_q .&& ccdf_a .>= high_q
    
    # Calculate linear fit
    x = log.(ag[idx_a])
    y = log.(ccdf_a[idx_a])
    A = [ones(length(x)) x] 
    coeffs = A \ y
    slope = coeffs[2]
    intercept = coeffs[1]
    fit_line = intercept .+ slope.*x
    
    # Plot data and fit
    plt_power_law = plot(x, y,
        lw=5,
        label=:none,
        color=default_colors[1]
    )
    plot!(x, fit_line,
        lw=3,
        label=:none,
        linestyle=:dash,
        color=default_colors[2]
    )
    find_x = findfirst(x .>= mean(x))
    range_y = maximum(y)-minimum(y)
    annotate!(x[find_x], y[find_x]+range_y/5,
        text("Slope = $(round(slope,digits=2))", :black, :left, 10, font("Computer Modern")))

    title!(title_label)
    ylabel!("Log CCDF")
    xlabel!("Log value")

    return plt_power_law
end