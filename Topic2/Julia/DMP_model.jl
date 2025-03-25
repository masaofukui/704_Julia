
function Jd_def(P, wd,Jd_p1,Jd)
    eval =  P["z shock"] - wd + P["beta"].*(1-P["s"])*Jd_p1 - Jd
    return eval;
end
Jd_def_inputs = arg_name(Jd_def)

function Free_entry_def(P,Jd_p1,theta )
    eval =  P["beta"]*qfun(P,theta).*Jd_p1 - P["c"]
    return eval;
end
Free_entry_def_inputs = arg_name(Free_entry_def)


function wd_def(P,wd,theta )
    eval =   (1-P["gamma"]).*P["b"] + P["gamma"].*(P["z shock"] + theta.*P["c"]) -wd
    return eval;
end
wd_def_inputs = arg_name(wd_def)


function ud_def(P,ud,theta ,ud_m1)
    eval =  ud - ud_m1 - ( P["s"].*(1-ud_m1) - ffun(P,theta).*ud_m1 )
    return eval;
end
ud_def_inputs = arg_name(ud_def)


allinputs_string = "theta,Jd,wd,ud"




varargin = [
    Jd_def,Jd_def_inputs,
    Free_entry_def,Free_entry_def_inputs,
    wd_def,wd_def_inputs,
    ud_def,ud_def_inputs
    ]
    


references = process_inputs(allinputs_string,varargin)
n = length(references)

varargin_eq = varargin[1:2:end]
nvar = Int(length(varargin)/2);
allinputs = replace(allinputs_string,"\n" => "")
allinputs = split(allinputs,',');
n = length(allinputs);