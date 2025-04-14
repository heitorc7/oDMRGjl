using ITensors
using ITensorMPS

ITensors.space(::SiteType"TwoSpinHalf") = 4

#Define all states
function ITensors.state(::SiteType"TwoSpinHalf", st::AbstractString)
    if st == "UpUp"
        return 1
    elseif st == "UpDn"
        return 2
    elseif st == "DnUp"
        return 3
    elseif st == "DnDn"
        return 4
    end
    throw(ArgumentError("State $st not recognized for twoSpinHalf"))
end

function ITensors.op(::SiteType"TwoSpinHalf", s::Index, opname::AbstractString)
    sP = prime(s)
    Op = ITensor(dag(s), sP)
    # State number definitions (1=UpUp, 2=UpDn, 3=DnUp, 4=DnDn)
    if opname == "SzL"
        Op[s(1),sP(1)] = 0.5
        Op[s(2),sP(2)] = -0.5
        Op[s(3),sP(3)] = 0.5
        Op[s(4),sP(4)] = -0.5
    elseif opname == "SzR"
        Op[s(1),sP(1)] = 0.5
        Op[s(2),sP(2)] = 0.5
        Op[s(3),sP(3)] = -0.5
        Op[s(4),sP(4)] = -0.5
    elseif opname == "SxL"
        Op[s(2),sP(1)] = 0.5
        Op[s(1),sP(2)] = 0.5
        Op[s(4),sP(3)] = 0.5
        Op[s(3),sP(4)] = 0.5
    elseif opname == "SxR"
        Op[s(3),sP(1)] = 0.5
        Op[s(4),sP(2)] = 0.5
        Op[s(1),sP(3)] = 0.5
        Op[s(2),sP(4)] = 0.5
    elseif opname == "SyL"
        Op[s(2),sP(1)] = -0.5im
        Op[s(1),sP(2)] = 0.5im
        Op[s(4),sP(3)] = -0.5im
        Op[s(3),sP(4)] = 0.5im
    elseif opname == "SyR"
        Op[s(3),sP(1)] = 0.5im
        Op[s(4),sP(2)] = 0.5im
        Op[s(1),sP(3)] = -0.5im
        Op[s(2),sP(4)] = -0.5im
    elseif opname == "SmL"
        Op[s(1),sP(2)] = 1.0
        Op[s(3),sP(4)] = 1.0
    elseif opname == "SmR"
        Op[s(3),sP(1)] = 1.0
        Op[s(4),sP(2)] = 1.0
    elseif opname == "SpL"
        Op[s(2),sP(1)] = 1.0
        Op[s(4),sP(3)] = 1.0
    elseif opname == "SpR"
        Op[s(1),sP(3)] = 1.0
        Op[s(2),sP(4)] = 1.0
    elseif opname == "SmLSpR"
        Op[s(1),sP(4)] = 1.0
    elseif opname == "SpLSmR"
        Op[s(4),sP(1)] = 1.0
    elseif opname == "SpSmL"
        Op[s(1),sP(1)] = 1.0
        Op[s(3),sP(3)] = 1.0
    elseif opname == "SpSmR"
        Op[s(1),sP(1)] = 1.0
        Op[s(2),sP(2)] = 1.0
    elseif opname == "SmSpL"
        Op[s(2),sP(2)] = 1.0
        Op[s(4),sP(4)] = 1.0
    elseif opname == "SmSpR"
        Op[s(3),sP(3)] = 1.0
        Op[s(4),sP(4)] = 1.0
    else
        throw(ArgumentError("Operator $opname not defined for twoSpinHalf"))
    end
    return Op
end