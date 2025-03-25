
function steady_state(guess,references_ss,P,varargin_eq_ss,Shock_Var,Shock_Path)
    n = length(varargin_eq_ss);
    @assert (n == length(guess) && n == length(references_ss));
    x = guess;
    count = 0;

    all_eq =  replace(string(varargin_eq_ss),"Function[" => "")
    all_eq =  replace(all_eq,"Any[" => "")
    all_eq =  replace(all_eq,"]" => "")
    all_eq =  replace(all_eq," " => "")
    all_eq = split(all_eq,',');


    maxcount = 40;
    while count < maxcount    # shouldn't take too many iterations!
        count = count + 1;
        input_mat = repeat(x',outer = [1 1 3]);
        outputs,ders = all_computations(input_mat,references_ss,P,varargin_eq_ss,Shock_Var,Shock_Path)
        
        

        J = dropdims( sum(ders,dims = 4), dims=(1,4))

     
        y = outputs[:];


        error = maximum(abs.(y));
        if error < 1E-10
            println("SS Found")
            break
        end

        x = x - J\y;
    end
    
    @assert (count < maxcount)
    return x
end

