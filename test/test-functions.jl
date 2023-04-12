function test_eq(a::Any, b::Any)
    if (length(a) == 0) && (length(b) == 0)
        return true
    elseif length(a) != length(b)
        return false
    else
        return all(a .≈ b)
    end
end

function tpt_test(ftest::String, tpt_homog::TPTHomog, tpt_res::TPTStats)
    testf = h5open(ftest)

    # inds
    for fn in fieldnames(typeof(tpt_homog.sets))
        tpt = getfield(tpt_homog.sets, fn)
        test = read(testf["tpt_homog/inds/$fn"])

        if !test_eq(tpt, test) 
            @error("$fn mismatch") 
            return false
        end 
    end

    # stats
    for fn in fieldnames(typeof(tpt_res))
        tpt = getfield(tpt_res, fn)
        test = read(testf["tpt_homog/stats/$fn"])
        
        if !test_eq(tpt, test) 
            @error("$fn mismatch") 
            return false
        end 
    end

    close(testf)

    return true
end