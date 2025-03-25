using Parameters, Setfield

function process_inputs(allinputs_string,varargin)
    @assert (mod(length(varargin),2) == 0); # even number of inputs

    n = Int(length(varargin)/2);
    allinputs_string = replace(allinputs_string,"\n" => "")
    allinputs = split(allinputs_string,',');


    @assert (length(allinputs) == n); # of equations = # of unknowns
    references = Any[]

    for i = 1:n
        println("i=",i)
        # count the number of arguments
        nspec = first(methods(varargin[2*i-1])).nargs-2;
        specinputs = split(varargin[2*i],',');
        @assert (length(specinputs) == nspec); # inputs string should include # of function inputs except P
        
        refarray = zeros(Int64,nspec,2);
        
        for j = 1:nspec
            #println("j=",j)
            input = specinputs[j];
            pos = 0;
            inputname = input;
            if length(input) >= 3
                if cmp(input[(end-2):end],"_m1") == 0
                    pos = -1;
                    inputname = input[1:(end-3)];
                elseif cmp(input[(end-2):end],"_p1") == 0
                    pos = 1;
                    inputname = input[1:(end-3)];
                end
            end
            inputplace = findall(cmp.(inputname,allinputs) .== 0);
            @assert length(inputplace) == 1 "Error with "*inputname*" in "*varargin[2*i]; # should be exactly one match of input name in master list of inputs
            refarray[j,1] = inputplace[1];
            refarray[j,2] = pos;
        end
        push!(references,refarray);
    end
    
    return references
end



function system_solve(ss0,ss1,references,P,varargin_eq,T,Shock_Var,Shock_Path; trate = 0.06)
    
    n = length(references);
    @assert (n == length(varargin_eq));
    
    # start by assuming convergence from ss0 to ss1 within first Tconv
    mixes = ((1+trate).^(-(1:T)));
    xfull = mixes*ss0'+(1 .-mixes)*ss1';
    
    count = 0;
    maxcount = 100;
    
    Jdims = (T*n,T*n);

    all_eq =  replace(string(varargin_eq),"Function[" => "")
    all_eq =  replace(all_eq,"Any[" => "")
    all_eq =  replace(all_eq,"]" => "")
    all_eq =  replace(all_eq," " => "")
    all_eq = split(all_eq,',');

    
    while count < maxcount
        count = count+1
        J = zeros(T*n,T*n);
        
        # construct input matrix from xfull
        xfull_lag  = [ss0'; xfull[1:end-1,:]];
        xfull_lead = [xfull[2:end,:]; ss1'];
        
        input_mat = zeros(T,n,3);
        input_mat[:,:,1] = xfull_lag;
        input_mat[:,:,2] = xfull;
        input_mat[:,:,3] = xfull_lead;
        
        outputs,ders = all_computations(input_mat,references,P,varargin_eq,Shock_Var,Shock_Path);
        
        
        
        
        errors = maximum(abs.(outputs))
        error = maximum(errors)
        println("Error = ", error)

        if error < 1E-7;
        #if count == 100
            break;
        end

        
        # assemble Jacobian piece by piece!
        # 'i' is index of condition, 'j' is index of variable
        for i = 1:n
            for j = 1:n
                indices_lag = Base._sub2ind(Jdims,(((i-1)*T+2):(i*T)),(((j-1)*T+1):(j*T-1)));
                indices = Base._sub2ind(Jdims,(((i-1)*T+1):(i*T)),(((j-1)*T+1):(j*T)));
                indices_lead = Base._sub2ind(Jdims,(((i-1)*T+1):(i*T-1)),(((j-1)*T+2):(j*T)));
                
                J[indices_lag] = J[indices_lag] + ders[2:T,i,j,1];
                J[indices] = J[indices]+ ders[1:T,i,j,2];
                J[indices_lead] = J[indices_lead] + ders[1:(T-1),i,j,3];
            end
        end
        

        yvec = outputs[:];

        J = sparse(J)

        xupdate = - J \ yvec;

        xupdate = reshape(xupdate,T,n);
        xfull = xfull + xupdate;
    end
    @assert count < maxcount

    return xfull
end



function myeval(P,curfunc,allinputs_t)
    output = curfunc(P,allinputs_t...);
    return output
end


function all_computations(input_mat,references,P,varargin_eq,Shock_Var,Shock_Path)
    
    allfunctions = varargin_eq;
    
    h = 1e-8;
    T = size(input_mat,1);      # number of time periods doing this for
    n = length(allfunctions);
    
    # should have full list of functions, matching size of "references" AND
    # matching second dimension of input matrix (# of eqs = # of unknowns)
    @assert (n == length(references) && n == size(input_mat,2))
    @assert (size(input_mat,3) == 3) # should have one lag and one lead
    
    outputs = zeros(T,n);
    ders = zeros(T,n,n,3);
    
    # go through each function and find the derivative at these inputs!
    for i = 1:n
        curfunc = allfunctions[i];
        refarray = references[i];
        
        linearindices = Base._sub2ind((n, 3),refarray[:,1],refarray[:,2].+2);
        ninputs = length(linearindices);


        for t = 1:T
            temp = input_mat[t,:,:];
            temp = temp[:]
            allinputs_t = temp[linearindices]
            # shock!
            P[Shock_Var] = Shock_Path[t];
            # in the final period, set Q =0
            
            if t == T
                #P["final"] = 1;
                P["final"] = 0;
            else
                P["final"] = 0;
            end

            outputs[t,i] = myeval(P,curfunc,allinputs_t);

           
            # obtain gradient by going through each input and taking partial

            for j = 1:ninputs
                allinputs_high = copy(allinputs_t);
                allinputs_low  = copy(allinputs_t);
                allinputs_high[j] = allinputs_high[j] + h;
                allinputs_low[j] = allinputs_low[j] - h;
                
                derivatives = (myeval(P,curfunc,allinputs_high) .- myeval(P,curfunc,allinputs_low))/(2*h);
                
                ders[t,i,refarray[j,1],refarray[j,2]+2] = copy(derivatives);
            end
        end
        
        
       
    end

    return outputs,ders
end
    
    
function method_argnames(m::Method)
    argnames = ccall(:jl_uncompress_argnames, Vector{Symbol}, (Any,), m.slot_syms)
    isempty(argnames) && return argnames
    return argnames[1:m.nargs]
end

function arg_name(f)
    ms = collect(methods(f))
    argname_temp = method_argnames(last(ms))[2:end]
    argname_temp = string(argname_temp)
    argname = replace(argname_temp,"[:P," => "" )
    argname = replace(argname," :" => "" )
    argname = replace(argname,"]" => "" )
    return argname
end
