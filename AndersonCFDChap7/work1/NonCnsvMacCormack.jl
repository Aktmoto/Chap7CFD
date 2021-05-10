module NonCnsvMacCormack
    #比熱比
    htSpecRatio = 1.4
    courant = 0.5

    mutable struct AverageDif
        ave_dρ::Vector{Float64}
        ave_dv::Vector{Float64}
        ave_dt::Vector{Float64}
        AverageDif(length) = new(zeros(length),zeros(length),zeros(length))
    end

    mutable struct PrimitiveValuesWithStepNum
        step:: Int64
        p::Vector{Float64}
        M::Vector{Float64}
        v::Vector{Float64}
        t::Vector{Float64}
        ρ::Vector{Float64}
        aveDif::AverageDif
        PrimitiveValuesWithStepNum(count,length) = new(count,zeros(length),zeros(length),zeros(length),zeros(length),zeros(length),AverageDif(length))
    end

    mutable struct Dif
        dρ::Vector{Float64}
        dv::Vector{Float64}
        dt::Vector{Float64}
    end
    mutable struct PredictedDif
        pre_dρ::Vector{Float64}
        pre_dv::Vector{Float64}
        pre_dt::Vector{Float64}
    end
    mutable struct PredictedValue
        v::Vector{Float64}
        t::Vector{Float64}
        ρ::Vector{Float64}
    end

    # countは今何回目かということ.
    function sub_superOneStep(deltax::Float64, count::Int64, length::Int64, ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64}, A::Vector{Float64},courant::Float64)
        next::PrimitiveValuesWithStepNum = PrimitiveValuesWithStepNum(count,length)
        dif::Dif = finiteDif(deltax, length,ρ,v,t,A)
        timeStep::Float64 = calcTimeStepWithCourant(length,deltax,v,t,courant)
        preValue::PredictedValue = predict(length,ρ,v,t,dif,timeStep)
        preDif::PredictedDif = correctDif(deltax,length,preValue,A)
        aveDif::AverageDif = averageDif(length,dif,preDif)
        calcNextTimeValue(length,next,aveDif,timeStep,ρ,v,t)
        sub_superCalcEdgeValue(length, next)
        sub_superCalcPressureAndMach(length, next)
        next.aveDif = aveDif
        return next
    end

    function sub_subOneStep(deltax::Float64, count::Int64, length::Int64, ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64}, A::Vector{Float64}, specificBackPressure::Float64)
        next::PrimitiveValuesWithStepNum = PrimitiveValuesWithStepNum(count,length)
        dif::Dif = finiteDif(deltax, length,ρ,v,t,A)
        timeStep::Float64 = calcTimeStep(length,deltax,v,t)
        preValue::PredictedValue = predict(length,ρ,v,t,dif,timeStep)
        preDif::PredictedDif = correctDif(deltax,length,preValue,A)
        aveDif::AverageDif = averageDif(length,dif,preDif)
        calcNextTimeValue(length,next,aveDif,timeStep,ρ,v,t)
        # 亜音速-亜音速は境界条件が異なる
        sub_subCalcEdgeValue(length, next,specificBackPressure)
        sub_subCalcPressureAndMach(length, next,specificBackPressure)
        return next
    end

    function finiteDif(deltax::Float64, length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64},A::Vector{Float64})::Dif
        x = FrontContEq(deltax, length, ρ,v,A)
        y = FrontNseq(deltax,length,ρ, v, t,)
        z = FrontEnergyEq(deltax, length,v,t,A)
        return Dif(x,y,z)
    end

    #前進差分
    function FrontContEq(deltax::Float64,length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, A::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length-1
             x[i] = (-ρ[i]*(v[i+1]-v[i])/deltax
            -ρ[i]v[i]*(log(A[i+1])-log(A[i]))/deltax
            -v[i]*(ρ[i+1]-ρ[i])/deltax)
        end
        return x
    end

    function FrontNseq(deltax::Float64,length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length-1
            x[i] = (-v[i]*(v[i+1]-v[i])/deltax
            -1/htSpecRatio*((t[i+1]-t[i])/deltax
            +(t[i]/ρ[i])*(ρ[i+1]-ρ[i])/deltax))
        end
        return x
    end

    function FrontEnergyEq(deltax::Float64,length::Int64,v::Vector{Float64}, t::Vector{Float64}, A::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length-1
            x[i] = (-v[i]*(t[i+1]-t[i])/deltax
            -(htSpecRatio-1)*t[i]*((v[i+1]-v[i])/deltax
            +v[i]*(log(A[i+1])
            -log(A[i]))/deltax))
        end
        return x
    end

    function predict(length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64}, dif::Dif, timeStep::Float64)::PredictedValue
        x::PredictedValue = PredictedValue(zeros(length),zeros(length),zeros(length))
        # ここインデックスの始点に注意. 
        for i in 1:1:length-1
            x.ρ[i] = (ρ[i] + dif.dρ[i]*timeStep)
            x.v[i] = (v[i] + dif.dv[i]*timeStep)
            x.t[i] = (t[i] + dif.dt[i]*timeStep)
        end
        return x
    end
    #後進差分
    function RearContEq(deltax::Float64,length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, A::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 2:1:length-1
            x[i] = (-ρ[i]*(v[i]-v[i-1])/deltax
            -ρ[i]v[i]*(log(A[i])
            -log(A[i-1]))/deltax
            -v[i]*(ρ[i]-ρ[i-1])/deltax)
        end
        return x
    end

    function RearNseq(deltax::Float64,length::Int64,ρ::Vector{Float64}, v::Vector{Float64}, t::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 2:1:length-1
            x[i] = (-v[i]*(v[i]-v[i-1])/deltax
            -1/htSpecRatio*((t[i]-t[i-1])/deltax
            +t[i]/ρ[i]*(ρ[i]-ρ[i-1])/deltax))
        end
        return x
    end

    function RearEnergyEq(deltax::Float64,length::Int64,v::Vector{Float64}, t::Vector{Float64}, A::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 2:1:length-1
            x[i] = (-v[i]*(t[i]-t[i-1])/deltax
            -(htSpecRatio-1)*t[i]*((v[i]-v[i-1])/deltax
            +v[i]*(log(A[i])-log(A[i-1]))/deltax))
        end
        return x
    end
    function correctDif(deltax::Float64,length::Int64,x::PredictedValue, A::Vector{Float64})::PredictedDif
        w = RearContEq(deltax, length, x.ρ,　x.v,　A)
        y = RearNseq(deltax,length,　x.ρ, x.v, x.t)
        z = RearEnergyEq(deltax, length, x.v, x.t,　A)
        return PredictedDif(w,y,z)
    end

    function averageDif(length::Int64, x::Dif, y::PredictedDif)::AverageDif
        aveDif::AverageDif = AverageDif(length)
         for i in 2:1:length-1
            aveDif.ave_dρ[i] = ((x.dρ[i]+y.pre_dρ[i])/2)
            aveDif.ave_dv[i] = ((x.dv[i]+y.pre_dv[i])/2)
            aveDif.ave_dt[i] = ((x.dt[i]+y.pre_dt[i])/2)
         end
         
         return aveDif
    end

    function calcTimeStep(length::Int64,deltax::Float64,v::Vector{Float64}, t::Vector{Float64})::Float64
        timeStep= Inf
        for i in 2:1:length-1
            tmp= courant*deltax/(sqrt(t[i])+v[i])
            timeStep = ifelse(tmp<timeStep, tmp, timeStep)
        end
        return timeStep
    end
    
    function calcTimeStepWithCourant(length::Int64,deltax::Float64,v::Vector{Float64}, t::Vector{Float64}, specificCourant::Float64)::Float64
        timeStep= Inf
        for i in 2:1:length-1
            tmp= specificCourant*deltax/(sqrt(t[i])+v[i])
            timeStep = ifelse(tmp<timeStep, tmp, timeStep)
        end
        return timeStep
    end
    function calcNextTimeValue(length::Int64,
        next::PrimitiveValuesWithStepNum, 
        aveDif::AverageDif,
        timeStep::Float64,
        ρ::Vector{Float64},
        v::Vector{Float64}, 
        t::Vector{Float64})
        for i in 2:1:length-1
            next.ρ[i] = (ρ[i]+aveDif.ave_dρ[i]*timeStep)
            next.v[i] = (v[i]+aveDif.ave_dv[i]*timeStep)
            next.t[i] = (t[i]+aveDif.ave_dt[i]*timeStep)
        end
    end
    #境界条件　亜音速-超音速
    function sub_superCalcEdgeValue(length::Int64, result::PrimitiveValuesWithStepNum)
        #流入側
        result.v[1] = (2*result.v[2]-result.v[3])
        result.ρ[1] = 1
        result.t[1] = 1
        #流出側
        result.v[length] = (2*result.v[length-1]-result.v[length-2])
        result.ρ[length] = (2*result.ρ[length-1]-result.ρ[length-2])
        result.t[length] = (2*result.t[length-1]-result.t[length-2])
    end

    function sub_superCalcPressureAndMach(length::Int64, result::PrimitiveValuesWithStepNum)
        for i in 1:1:length
            result.p[i] =(result.ρ[i]*result.t[i])
            result.M[i] =(result.v[i]/(sqrt(result.t[i])))
        end
    end
    #境界条件　亜音速-亜音速
    function sub_subCalcEdgeValue(length::Int64, result::PrimitiveValuesWithStepNum, specificBackPressure::Float64)
        #流入側
        result.v[1] = (2*result.v[2]-result.v[3])
        result.ρ[1] = 1
        result.t[1] = 1
        #流出側
        result.v[length] = (2*result.v[length-1]-result.v[length-2])
        result.ρ[length] = (2*result.ρ[length-1]-result.ρ[length-2])
        result.t[length] = specificBackPressure/ result.ρ[length]
    end

    function sub_subCalcPressureAndMach(length::Int64, result::PrimitiveValuesWithStepNum, specificBackPressure::Float64)
         for i in 1:1:length
            result.p[i] = ifelse(i!=length,(result.ρ[i]*result.t[i]),specificBackPressure)
            result.M[i] =(result.v[i]/(sqrt(result.t[i])))
        end
    end
end