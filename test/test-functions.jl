function test_eq(a::Any, b::Any)
    if (length(a) == 0) && (length(b) == 0) && a == b
        return true
    elseif length(a) != length(b)
        return false
    else
        return all(a .≈ b)
    end
end

function tpt_test(ftest::String, tpt_result::AbstractTPTHomogResult)
    testf = h5open(ftest)

    for fn in fieldnames(typeof(tpt_result))
        if fn == :sets
            for fns in fieldnames(typeof(tpt_result.sets))
                tpt = getfield(tpt_result.sets, fns)
                test = read(testf["tpt_homog/indices/$fns"])
        
                if !test_eq(tpt, test)
                    @show tpt
                    @show test
                    @error("$fns mismatch") 
                    return false
                end 
            end
        else
            tpt = getfield(tpt_result, fn)
            test = read(testf["tpt_homog/statistics/$fn"])
            
            if !test_eq(tpt, test) 
                @show tpt
                @show test
                @error("$fn mismatch") 
                return false
            end 
        end
    end

    close(testf)

    return true
end