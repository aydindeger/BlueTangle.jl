using BenchmarkTools
function bench(a::BenchmarkTools.Trial, b::BenchmarkTools.Trial)
    t_ratio = median(a).time > median(b).time ? median(a).time / median(b).time : median(b).time / median(a).time
    m_ratio = median(a).memory > median(b).memory ? median(a).memory / median(b).memory : median(b).memory / median(a).memory
    
    rt = round(t_ratio, sigdigits=3)
    rm = round(m_ratio, sigdigits=3)
    
    if median(a).time < median(b).time
        ti=1
        println("First is faster: $(rt)x")
    else
        ti=2
        println("Second is faster: $(rt)x")
    end
    
    if median(a).memory < median(b).memory
        mi=1
        println("First uses less memory: $(rm)x")
    else
        mi=2
        println("Second uses less memory: $(rm)x")
    end
    
    return [(ti,rt), (mi,rm)]
end