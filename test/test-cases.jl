function test_cases()
    # Ptest = [
    # 0.0374453  0.143588   0.122089  0.0933665  0.0647982   0.180506    0.164599   0.193608;
    # 0.140666   0.198493   0.224272  0.0121098  0.176694    0.0108653   0.0162701  0.2206298;
    # 0.166866   0.197242   0.127246  0.0308489  0.0287401   0.134896    0.194696   0.119465;
    # 0.0916178  0.0106214  0.136104  0.0747192  0.205926    0.195282    0.179782   0.1059476;
    # 0.111462   0.120019   0.186935  0.150579   0.0560114   0.00462604  0.184768   0.18559956;
    # 0.0417938  0.148309   0.16509   0.122614   0.111876    0.16153     0.150306   0.0984812;
    # 0.122551   0.0415795  0.235409  0.176576   0.00610855  0.152889    0.13524    0.12964695;
    # 0.0599739  0.159639   0.164662  0.0809561  0.207265    0.0321942   0.091017   0.2042928;
    # ]

    Ptest = [
        1/3 2/3 0 0 0 0 0 0;
        1/6 1/4 1/24 1/12 1/8 1/3 0 0;
        0 1/12 1/4 1/24 1/6 1/3 1/8 0;
        0 1/12 1/6 1/4 1/24 1/8 1/3 0;
        0 1/24 1/8 1/3 1/4 1/6 1/12 0;
        0 1/3 1/8 1/12 1/6 1/4 1/24 0;
        0 0 1/4 1/12 1/3 1/6 1/8 1/24;
        0 0 0 0 0 0 3/4 1/4;
    ]
    
    ### stationary

    Atest = [2]
    Btest = [7]


    t1_stat = [TPTHomog(Ptest, Atest, Btest), joinpath(@__DIR__, "TPT_8_no_int_stat.h5")]

    Atest = [2, 4]
    Btest = [7, 4]

    t2_stat = [TPTHomog(Ptest, Atest, Btest), joinpath(@__DIR__, "TPT_8_with_int_stat.h5")]

    test_cases_stationary = [
        t1_stat,
        t2_stat
    ]

    ### nonstationary

    Atest = [2]
    Btest = [7]

    t1_nonstat = [TPTHomog(Ptest, Atest, Btest), joinpath(@__DIR__, "TPT_8_no_int_nonstat.h5")]

    Atest = [2, 4]
    Btest = [7, 4]

    t2_nonstat = [TPTHomog(Ptest, Atest, Btest), joinpath(@__DIR__, "TPT_8_with_int_nonstat.h5")]

    test_cases_nonstationary = [
        t1_nonstat,
        t2_nonstat
    ]

    return test_cases_stationary, test_cases_nonstationary
end

function generate_tests()
    test_cases_stationary, test_cases_nonstationary = test_cases()

    for t in test_cases_stationary
        tpt_stat, fout = tpt_stationary_statistics(t[1]), t[2]
        parts = minimal_partitions(t[1], tpt_stat)
        tpt_write(fout, tpt_stat)
        tpt_write(fout, parts)
    end

    for t in test_cases_nonstationary
        tpt_nonstat, fout = tpt_nonstationary_statistics(t[1]), t[2]
        parts = minimal_partitions(t[1], tpt_nonstat)
        tpt_write(fout, tpt_nonstat)
        tpt_write(fout, parts)
    end

    return
end

### in src directory
# include("types.jl")
# include("tpt-homogeneous.jl")
# include("tpt-write.jl")
# include("partitions/partitions.jl")
# include("../test/test-cases.jl")
# generate_tests()