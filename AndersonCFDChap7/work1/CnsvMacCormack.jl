module CnsvMacCormack
    include("./NonCnsvMacCormack.jl")
    import .NonCnsvMacCormack
    #比熱比
    htSpecRatio = 1.4

    #解ベクトルへの変換用の関数
    function encodeTo_u1(ρ::Vector{Float64},A::Vector{Float64})::Vector{Float64}
        return [ρ[i]*A[i] for i in 1:1:length(A)]
    end
    function encodeTo_u2(ρ::Vector{Float64},A::Vector{Float64},v::Vector{Float64})::Vector{Float64}
        return [ρ[i]*A[i]*v[i] for i in 1:1:length(A)]
    end
    function encodeTo_u3(ρ::Vector{Float64},A::Vector{Float64},v::Vector{Float64},t::Vector{Float64})::Vector{Float64}
        return [ρ[i]*A[i]*(t[i]/(htSpecRatio-1)+htSpecRatio/2*(v[i])^2) for i in 1:1:length(A)]
    end
    

    mutable struct DecodedPrimitiveValues
        p::Vector{Float64}
        M::Vector{Float64}
        v::Vector{Float64}
        t::Vector{Float64}
        ρ::Vector{Float64}
        DecodedPrimitiveValues(length) = new(zeros(length),zeros(length),zeros(length),zeros(length),zeros(length))
    end

    mutable struct SolutionVectorsWithStepAndPrim
        step::Int64
        u1::Vector{Float64}
        u2::Vector{Float64}
        u3::Vector{Float64}
        primitive::DecodedPrimitiveValues
        SolutionVectorsWithStepAndPrim(length) = new(0,zeros(length),zeros(length),zeros(length),DecodedPrimitiveValues(length))
    end


    mutable struct FluxAndSourceVector
        f1::Vector{Float64}
        f2::Vector{Float64}
        f3::Vector{Float64}
        j2::Vector{Float64}
        FluxAndSourceVector(length) = new(zeros(length),zeros(length),zeros(length),zeros(length))
    end

    mutable struct PredictedFluxVector
        pre_f1::Vector{Float64}
        pre_f2::Vector{Float64}
        pre_f3::Vector{Float64}
        PredictedFluxVector(length) = new(zeros(length),zeros(length),zeros(length))
    end

    mutable struct PredictedPrimitiveValue
        pre_ρ::Vector{Float64}
        pre_t::Vector{Float64}
        PredictedPrimitiveValue(length) = new(zeros(length),zeros(length))
    end

    mutable struct Dif
        du1::Vector{Float64}
        du2::Vector{Float64}
        du3::Vector{Float64}
        Dif(length) = new(zeros(length),zeros(length),zeros(length))
    end

    mutable struct PredictedDif
        pre_du1::Vector{Float64}
        pre_du2::Vector{Float64}
        pre_du3::Vector{Float64}
        PredictedDif(length) = new(zeros(length),zeros(length),zeros(length))
    end
    mutable struct AverageDif
        ave_du1::Vector{Float64}
        ave_du2::Vector{Float64}
        ave_du3::Vector{Float64}
        AverageDif(length) = new(zeros(length),zeros(length),zeros(length))
    end
    mutable struct PredictedSolutionVector
        pre_u1::Vector{Float64}
        pre_u2::Vector{Float64}
        pre_u3::Vector{Float64}
        PredictedSolutionVector(length) = new(zeros(length),zeros(length),zeros(length))
    end
    
    function sub_superOneStep(deltax::Float64, 
        count::Int64, 
        length::Int64, 
        ρ::Vector{Float64},#ソースベクターを求める際に必要
        t::Vector{Float64},
        v::Vector{Float64},
        u1::Vector{Float64}, 
        u2::Vector{Float64}, 
        u3::Vector{Float64}, 
        A::Vector{Float64})::SolutionVectorsWithStepAndPrim
        next::SolutionVectorsWithStepAndPrim = SolutionVectorsWithStepAndPrim(length)
        next.step = count
        
        fluxAndSource::FluxAndSourceVector = calcFluxAndSource(length,deltax,u1,u2,u3,ρ,t,A)
        timeStep::Float64 = NonCnsvMacCormack.calcTimeStep(length,deltax,v,t)
        #予測子
        dif::Dif = calcDif(length,deltax,fluxAndSource)
        preSolution::PredictedSolutionVector = calcPreSolution(length,u1,u2,u3,dif,timeStep)
        prePrimitive::PredictedPrimitiveValue = calcPrePrimitive(length,preSolution,A)
        preFluxAndSource::FluxAndSourceVector = calcFluxAndSource(length,deltax,preSolution.pre_u1,
        preSolution.pre_u2,
        preSolution.pre_u3,
        prePrimitive.pre_ρ,
        prePrimitive.pre_t,
        A)
        #修正子
        preDif::PredictedDif = calcPreDif(length,deltax,preFluxAndSource,prePrimitive,A)
        aveDif::AverageDif = calcAveDif(length,dif,preDif)
        calcNextTimeValue(length,next,u1,u2,u3,aveDif,timeStep)
        #境界条件
        sub_superCalcEdgeValue(length,next,A)
        CalcPrimitive(length,next,A)
        return next
    end

    function sub_subOneStepWithArtificialViscosity(deltax::Float64, 
        count::Int64, 
        length::Int64,
        specificPressure::Float64,
        ρ::Vector{Float64},#ソースベクターを求める際に必要
        t::Vector{Float64},
        v::Vector{Float64},
        p::Vector{Float64},
        u1::Vector{Float64}, 
        u2::Vector{Float64}, 
        u3::Vector{Float64}, 
        A::Vector{Float64},
        Cx::Float64)::SolutionVectorsWithStepAndPrim
        next::SolutionVectorsWithStepAndPrim = SolutionVectorsWithStepAndPrim(length)
        next.step = count
        
        fluxAndSource::FluxAndSourceVector = calcFluxAndSource(length,deltax,u1,u2,u3,ρ,t,A)
        timeStep::Float64 = NonCnsvMacCormack.calcTimeStep(length,deltax,v,t)
        #予測子
        dif::Dif = calcDif(length,deltax,fluxAndSource)
        preSolution::PredictedSolutionVector = calcPreSolutionWithArtificialViscosity(length,u1,u2,u3,dif,timeStep,Cx,p)
        prePrimitive::PredictedPrimitiveValue = calcPrePrimitive(length,preSolution,A)
        preFluxAndSource::FluxAndSourceVector = calcFluxAndSource(length,deltax,preSolution.pre_u1,
        preSolution.pre_u2,
        preSolution.pre_u3,
        prePrimitive.pre_ρ,
        prePrimitive.pre_t,
        A)
        #修正子
        preDif::PredictedDif = calcPreDif(length,deltax,preFluxAndSource,prePrimitive,A)
        aveDif::AverageDif = calcAveDif(length,dif,preDif)
        calcNextTimeValueWithArtificialViscosity(length,
        next,
        u1,
        u2,
        u3,
        aveDif,
        timeStep,
        Cx,
        prePrimitive.pre_ρ .* prePrimitive.pre_t,
        preSolution)
        #境界条件
        sub_subCalcEdgeValue(length,next,A,specificPressure)
        CalcPrimitive(length,next,A)
        return next
    end
    
    #region 7-5用
    function calcF1(u2::Vector{Float64})::Vector{Float64}
        return u2
    end
    function calcF2(length::Int64,u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length
            x[i] = (u2[i]^2/u1[i])+((htSpecRatio-1)/htSpecRatio)*(u3[i]-(htSpecRatio/2)*(u2[i]^2)/u1[i])
        end
        return x
    end
    function calcF3(length::Int64,u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length
            x[i] = htSpecRatio*u2[i]*u3[i]/u1[i]-htSpecRatio*(htSpecRatio-1)/2*(u2[i]^3/u1[i]^2)
        end
        return x
    end

    function calcJ2(length::Int64,deltax::Float64,ρ::Vector{Float64},t::Vector{Float64},A::Vector{Float64})::Vector{Float64}
        x::Vector{Float64} = zeros(length)
        for i in 1:1:length-1
            x[i] = (1/htSpecRatio)*ρ[i]*t[i]*((A[i+1])-(A[i]))/deltax
        end
        return x
            
    end

    function  calcFluxAndSource(length::Int64,
        deltax::Float64,
        u1::Vector{Float64},
        u2::Vector{Float64},
        u3::Vector{Float64},
        ρ::Vector{Float64},
        t::Vector{Float64},
        A::Vector{Float64})::FluxAndSourceVector
        fluxAndSource = FluxAndSourceVector(length)
        fluxAndSource.f1 = calcF1(u2)
        fluxAndSource.f2 = calcF2(length,u1,u2,u3)
        fluxAndSource.f3 = calcF3(length,u1,u2,u3)
        fluxAndSource.j2 = calcJ2(length,deltax,ρ,t,A)
        return fluxAndSource
    end

    function calcDif(length::Int64,deltax::Float64,fluxAndSource::FluxAndSourceVector)::Dif
        dif::Dif = Dif(length)
        f1 = fluxAndSource.f1
        f2 = fluxAndSource.f2
        f3 = fluxAndSource.f3
        j2 = fluxAndSource.j2
        for i in 1:1:length-1
            dif.du1[i] = -(f1[i+1]-f1[i])/deltax
            dif.du2[i] = -(f2[i+1]-f2[i])/deltax+j2[i]
            dif.du3[i] = -(f3[i+1]-f3[i])/deltax
        end
        return dif
    end

    function calcPreSolution(length::Int64,
        u1::Vector{Float64}, 
        u2::Vector{Float64},
        u3::Vector{Float64},
        dif::Dif,
        timeStep::Float64)::PredictedSolutionVector
        x::PredictedSolutionVector = PredictedSolutionVector(length)
        for i in 1:1:length-1
            x.pre_u1[i] = u1[i]+dif.du1[i]*timeStep
            x.pre_u2[i] = u2[i]+dif.du2[i]*timeStep
            x.pre_u3[i] = u3[i]+dif.du3[i]*timeStep
        end
        return x
    end
    function calcPrePrimitive(length::Int64, preSolution::PredictedSolutionVector,A::Vector{Float64})::PredictedPrimitiveValue
        x::PredictedPrimitiveValue = PredictedPrimitiveValue(length)
        u1 = preSolution.pre_u1
        u2 = preSolution.pre_u2
        u3 = preSolution.pre_u3
        for i in 1:1:length-1
            x.pre_ρ[i] = u1[i]/A[i]
            x.pre_t[i] = (htSpecRatio-1)*(u3[i]/u1[i]-htSpecRatio/2*(u2[i]/u1[i])^2)
        end
        return x
    end
    function calcPreDif(length::Int64,deltax::Float64,preFluxAndSource::FluxAndSourceVector,prePrimitive::PredictedPrimitiveValue,A::Vector{Float64})::PredictedDif
        preDif::PredictedDif = PredictedDif(length)
        f1 = preFluxAndSource.f1
        f2 = preFluxAndSource.f2
        f3 = preFluxAndSource.f3
        ρ = prePrimitive.pre_ρ
        t = prePrimitive.pre_t
        for i in 2:1:length-1
            preDif.pre_du1[i] = -(f1[i]-f1[i-1])/deltax
            preDif.pre_du2[i] = -(f2[i]-f2[i-1])/deltax+1/htSpecRatio*ρ[i]*t[i]*((A[i])-(A[i-1]))/deltax
            preDif.pre_du3[i] = -(f3[i]-f3[i-1])/deltax
        end
        return preDif
    end

    function calcAveDif(length::Int64,dif::Dif,preDif::PredictedDif)::AverageDif
        aveDif::AverageDif = AverageDif(length)
        for i in 2:1:length-1
            aveDif.ave_du1[i] = (dif.du1[i]+preDif.pre_du1[i])/2
            aveDif.ave_du2[i] = (dif.du2[i]+preDif.pre_du2[i])/2
            aveDif.ave_du3[i] = (dif.du3[i]+preDif.pre_du3[i])/2
        end
        return aveDif
    end

    function calcNextTimeValue(length::Int64,
        next::SolutionVectorsWithStepAndPrim,
        u1::Vector{Float64},
        u2::Vector{Float64},
        u3::Vector{Float64},
        aveDif::AverageDif,
        timeStep::Float64)
        primitive = next.primitive
        for i in 2:1:length-1
            next.u1[i] = (u1[i]+aveDif.ave_du1[i]*timeStep)
            next.u2[i] = (u2[i]+aveDif.ave_du2[i]*timeStep)
            next.u3[i] = (u3[i]+aveDif.ave_du3[i]*timeStep)
        end
    end

    function sub_superCalcEdgeValue(length::Int64,next::SolutionVectorsWithStepAndPrim,A::Vector{Float64})
        #流入側
        next.u1[1] = A[1] #ρ[1]=1のため
        next.u2[1] = (2*next.u2[2]-next.u2[3])
        next.u3[1] = next.u1[1]*(1/(htSpecRatio-1)+htSpecRatio/2*(next.u2[1]/next.u1[1])^2)#t[1]=1のため
        #流出側
        next.u1[length] = (2*next.u1[length-1]-next.u1[length-2])
        next.u2[length] = (2*next.u2[length-1]-next.u2[length-2])
        next.u3[length] = (2*next.u3[length-1]-next.u3[length-2])
    end

    function CalcPrimitive(length::Int64,next::SolutionVectorsWithStepAndPrim,A::Vector{Float64})
        primitive = next.primitive
        for i in 1:1:length
            primitive.ρ[i] = next.u1[i]/A[i]
            primitive.v[i] = next.u2[i]/next.u1[i]
            primitive.t[i] = (htSpecRatio-1)*(next.u3[i]/next.u1[i]-htSpecRatio/2*(primitive.v[i]^2))
            primitive.p[i] = primitive.ρ[i]*primitive.t[i]
            primitive.M[i] = primitive.v[i]/sqrt(primitive.t[i])
        end
    end
    #endregion

    #region 7-6用 基本的に人工粘性を追加するだけ

    function calcArtificialVsicosityOf(i::Int64,p::Vector{Float64},u::Vector{Float64},Cx::Float64)::Float64
        return Cx*(abs(p[i+1]-2*p[i]+p[i-1]))/(p[i+1]+2*p[i]+p[i-1])*(u[i+1]-2*u[i]+u[i-1])
    end
    function calcPreSolutionWithArtificialViscosity(length::Int64,
        u1::Vector{Float64}, 
        u2::Vector{Float64},
        u3::Vector{Float64},
        dif::Dif,
        timeStep::Float64,
        Cx::Float64,
        p::Vector{Float64})::PredictedSolutionVector
        x::PredictedSolutionVector = PredictedSolutionVector(length)
        for i in 2:1:length-1
            x.pre_u1[i] = u1[i]+dif.du1[i]*timeStep+calcArtificialVsicosityOf(i,p,u1,Cx)
            x.pre_u2[i] = u2[i]+dif.du2[i]*timeStep+calcArtificialVsicosityOf(i,p,u2,Cx)
            x.pre_u3[i] = u3[i]+dif.du3[i]*timeStep+calcArtificialVsicosityOf(i,p,u3,Cx)
        end
            #i=1,length-1も必要であるが, 人工粘性は境界ではindex範囲外になって使えないので, 境界ではそのままの値
            x.pre_u1[1] = u1[1]
            x.pre_u2[1] = u2[1]
            x.pre_u3[1] = u3[1]
            x.pre_u1[length] = u1[length]
            x.pre_u2[length] = u2[length]
            x.pre_u3[length] = u3[length]
        return x
    end

    function calcNextTimeValueWithArtificialViscosity(length::Int64,
        next::SolutionVectorsWithStepAndPrim,
        u1::Vector{Float64},
        u2::Vector{Float64},
        u3::Vector{Float64},
        aveDif::AverageDif,
        timeStep::Float64,
        Cx::Float64,
        p::Vector{Float64},
        pre_U::PredictedSolutionVector)
        primitive = next.primitive
        for i in 2:1:length-1
            next.u1[i] = (u1[i]+aveDif.ave_du1[i]*timeStep)+calcArtificialVsicosityOf(i,p,pre_U.pre_u1,Cx)
            next.u2[i] = (u2[i]+aveDif.ave_du2[i]*timeStep)+calcArtificialVsicosityOf(i,p,pre_U.pre_u2,Cx)
            next.u3[i] = (u3[i]+aveDif.ave_du3[i]*timeStep)+calcArtificialVsicosityOf(i,p,pre_U.pre_u3,Cx)
        end
    end

    function sub_subCalcEdgeValue(length::Int64,next::SolutionVectorsWithStepAndPrim,A::Vector{Float64},specificPressure)
        #流入側
        next.u1[1] = A[1] #ρ[1]=1のため
        next.u2[1] = (2*next.u2[2]-next.u2[3])
        next.u3[1] = next.u1[1]*(1/(htSpecRatio-1)+htSpecRatio/2*(next.u2[1]/next.u1[1])^2)#t[1]=1のため
        #流出側
        next.u1[length] = (2*next.u1[length-1]-next.u1[length-2])
        next.u2[length] = (2*next.u2[length-1]-next.u2[length-2])
        next.u3[length] = specificPressure*A[length]/(htSpecRatio-1)+htSpecRatio/2*next.u2[length]*next.primitive.v[length]
    end
end