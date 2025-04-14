function createMagProf(L::Int, profile)
    h = [1:length(L);]
    if profile == "parabolic"
        return h
    end
    if profile == "linear"
        return h
    end
    return h
end
# To be done
# Don't feel like doing it right now