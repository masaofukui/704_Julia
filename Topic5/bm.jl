using Debugger
using SparseArrays
using Plots
using DelimitedFiles
using DataFrames
using LinearAlgebra
using JLD
using Optim
using NLsolve
using Statistics
using Printf
using LaTeXStrings
using Plots.PlotMeasures
using Distributions
using QuadGK
using ForwardDiff
lw = 3;
default(; titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern", 
    legendfontfamily = "Computer Modern",
    titlefontsize = 15,
    xguidefontsize = 12,
    legendfontsize = 12,
    yguidefontsize = 12,
    xgrid = :none)



function compute_homo_bm(fE)
        underw = 0.4;
        z = 0.999
        v = 1;
        q = 1;
        delta = 0.02;

        G(w) = min.( (1 .+ delta./fE).*(1  .- sqrt.((z .- w)./(z .-underw))),1.0)

        wg = (underw-0.1):0.001:(z);
        
        g = zeros(length(wg))
        a = ForwardDiff.jacobian(G, wg[wg .>= underw])
        gtemp = diag(a)
        G_init = G(wg);
        G_init[wg .<= underw] .= 0.0;
        g[wg .>= underw] .= gtemp

        lfun(w) = q.*v.*delta./( (delta + fE.*(1-G(w))).^2)
        lg = lfun.(wg)
        lg[wg .< underw] .= 0.0
        lg[G_init .>=1] .= 0.0
        total_l = sum(lg.*(wg[2]-wg[1]))
        lg = lg./total_l

        mug = (z .- wg)./wg;
        mug = vcat(mug)
        mug[wg .< underw] .= NaN
        mug[G_init .>=1] .= NaN

        Hg = delta.*G_init./(delta .+ fE.*(1 .- G_init))
        hg = diff(Hg)./diff(wg)

        lw = 3
        gplot = plot(wg,g,linewidth = lw,legend=:none)
        plot!(title="Wage Offer Distribution, \$g(w)\$")
        plot!(xlabel="Wage, \$w\$")

        
        lplot = plot(wg,lg,linewidth = lw,legend=:none)
                plot!(title="Firm Size Distribution, \$l(w)\$")
                plot!(xlabel="Wage, \$w\$")
     

        muplot = plot(wg,mug,linewidth = lw,legend=:none)
                plot!(title="Wage Markdown, \$(z-w)/w\$")
                plot!(xlabel="Wage, \$w\$")

        hplot = plot(wg[1:(end-1)],hg,linewidth = lw,legend=:none)
                plot!(title="Wage Distribution, \$h(w)\$")
                plot!(xlabel="Wage, \$w\$")


        plt = plot(gplot,hplot,muplot,layout=(1,3),size=(1000,300))
        plot!(margin = 4mm)
        display(plt)
        savefig("./figure/bm_homo"*string(fE)*".pdf")
end
compute_homo_bm(0.05)
compute_homo_bm(0.0001)
compute_homo_bm(10000)

function compute_bm(underw)
        q = 1;
        v = 1;
        σ = 0.4;
        delta = 0.02;
        μ = 0
        d = LogNormal(μ,σ)
        #d = Pareto(10,underw)
        J(x) = cdf(d,x)
        J_pdf(x) = pdf(d,x)

        zg = range(0,4,length=10000)
        wg = vcat( copy(zg))
        for i = eachindex(zg)
        z = zg[i]
        w_temp(tildez) = (delta.*(1-J.(underw)) + fE.*(1-J(z))).^2 ./
        ( (delta.*(1-J(underw)) + fE.*(1-J(tildez))).^2) 
        integ_temp = quadgk(w_temp,underw,z)[1]
        wg[i] = z - integ_temp;
        end


        wg[zg .<= underw] .= NaN
        Gz(z) = (J(z) - J(underw))./(1-J(underw))

        Gzg = Gz.(zg)
        Gzg[zg .<= underw] .= NaN
        Gzg[Gzg .>= 1.0] .= 1.0

        lzg = q.*v.*delta./( (delta .+ fE.*(1 .-Gzg)).^2)
        lzg[zg .<= underw] .= 0


        Hzg = delta.*Gzg./(delta .+ fE .*(1 .- Gzg))
        dHzg = diff(Hzg)./diff(zg);
        dwzg = diff(wg)./diff(zg)
        plot(zg,Hzg)
        hzg = diff(Hzg)./diff(wg)
        plot(wg,hzg)
        hzz = diff(Hzg)./diff(zg)
        hzz[isnan.(hzz)] .= 0.0
        J_pdf_z = J_pdf.(zg)

        return wg, lzg,Hzg,hzg,Gzg,hzz,zg,J_pdf_z
end


wg, lzg,Hzg,hzg,Gzg,hzz,zg,J_pdf_z = compute_bm(0.7)


zplot = plot(zg,J_pdf_z,linewidth = lw,legend=:none)
plot!(title="Productivity Distribution, \$J'(z)\$")
plot!(xlabel="Productivity, \$z\$")


wplot = plot(zg,wg,linewidth = lw,label="\$w(z)\$")
plot!(zg,zg,linestyle=:dash,linwidth = lw,label="\$z\$")
plot!(title="Wage, \$w(z)\$")
plot!(xlabel="Productivity, \$z\$")


muplot = plot(zg,(zg .- wg)./wg,linewidth = lw,label="w")
plot!(legend=:none)
plot!(title="Wage Markdown, \$(z-w)/w\$")
plot!(xlabel="Productivity, \$z\$")



lplot = plot(zg,lzg,linewidth = lw,legend=:none)
plot!(title="Firm Size, \$l(w(z))\$")
plot!(xlabel="Productivity, \$z\$")



plot(zplot,wplot,muplot,lplot,layout=(2,2),size=(1000,600))
plot!(margin = 8mm)
savefig("./figure/bm_hetero.pdf")


hplot = plot(wg,hzg,linewidth = lw,legend=:none)
plot!(title="Wage Distribution, \$h(w)\$")
plot!(xlabel="Wage, \$w\$")
ylims!(0.0,3.0)
savefig("./figure/h_hetero.pdf")
