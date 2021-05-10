#=
強引だが,  matplotlibでの日本語表示のために
C:\Users\(ユーザー名)\.julia\conda\3\pkgs\matplotlib-base-3.3.4-py38h49ac443_0\Lib\site-packages\matplotlib\mpl-data
にあるmatplotlibrc内の#font.family : sans-serifをfont.family : MS Gothicに変更する.
入れてあるパッケージ
GR v0.57.4    
IJulia v1.23.2
Plots v1.13.2 
PyCall v1.92.3
PyPlot v2.9.0
=#

#=
反射してpが振動するなら, 出口側圧力をがっつり固定しないで, 例えば, 固定圧力と線形化から出てきた圧力との平均を取ってなましてやる
@timeを用いて, 格子点の数が処理時間にどのくらい影響を与えるのかを調べる. 
https://github.com/Aktmoto/Chap7CFD
=#


include("./NonCnsvMacCormack.jl")
include("./CnsvMacCormack.jl")
include("./graphMaker.jl")
import .NonCnsvMacCormack
import .CnsvMacCormack
import .graphMaker

pointLength = 31
ductLength = 3
htSpecRatio = 1.4
deltax= ductLength/(pointLength-1)

#教科書のデータの桁と合わせて見やすくする.
function roundTo_3digits(x::Vector{Float64})::Vector{Float64}
    for i in 1:1:length(x)
        x[i] = round(x[i],digits=3)
    end
    return x
end

#region 7-3

#juliaの配列アクセスは1から
function area7_3(x::Vector{Float64})::Vector{Float64}
    a::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        a[i] = round(1+2.2(x[i]-1.5)^2,digits=3)
    end
    return a
end

function density7_3(x::Vector{Float64})::Vector{Float64}
    d::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        d[i] = round(1-0.3146*x[i],digits=3)
    end
    return d
end

function temperature7_3(x::Vector{Float64})::Vector{Float64}
    t::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        t[i] = round(1-0.2314*x[i],digits=3)
    end
    return t
end

function velocity7_3(x::Vector{Float64}, t::Vector{Float64})::Vector{Float64}
    v::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        v[i] = round((0.1+1.09*x[i])*sqrt(t[i]),digits=3)
    end
    return v
end


#配列はgeneratorで一気に作る.

x7_3 = [i for i in 0:deltax:ductLength]
ρ7_3 = density7_3(x7_3)
t7_3 = temperature7_3(x7_3)
v7_3 = velocity7_3(x7_3,t7_3)
A7_3 = area7_3(x7_3)

totalStep7_3 = 1400
function sub_superNozzleSolutionWith(courant::Float64)::Vector{NonCnsvMacCormack.PrimitiveValuesWithStepNum}
    one = NonCnsvMacCormack.sub_superOneStep(deltax,
                                1,
                                pointLength,
                                ρ7_3,
                                v7_3,
                                t7_3,
                                A7_3,
                                courant
                                )
    totalResult = Vector{NonCnsvMacCormack.PrimitiveValuesWithStepNum}()
    push!(totalResult,one)
    let
        prev = one
        for i in 2:1:totalStep7_3
            next = NonCnsvMacCormack.sub_superOneStep(deltax,
                                i,
                                pointLength,
                                prev.ρ,
                                prev.v,
                                prev.t,
                                A7_3,
                                courant
                                )
            push!(totalResult,next)
            prev = next
        end
    return totalResult
    end
end
totalResult7_3 = sub_superNozzleSolutionWith(0.5)
totalResult7_3_courant1 = sub_superNozzleSolutionWith(1.0)
totalResult7_3_limitCourant = sub_superNozzleSolutionWith(1.13)

#=
see = totalStep7_3
println("回数"*string(totalResult7_3[see].step))
println("圧力"*string(totalResult7_3[see].p))
println("速度"*string(totalResult7_3[see].v))
println("密度"*string(totalResult7_3[see].ρ))
println("温度"*string(totalResult7_3[see].t))
println("マッハ数"*string(totalResult7_3[see].M))
println("質量流量"*string([totalResult7_3[see].ρ[i]*totalResult7_3[see].v[i]*area7_3(x7_3)[i] for i in 1:pointLength]))
=#

#7-3 グラフ出力
rigid_M = [0.098,0.110,0.124,0.140,0.160,0.185,0.214,0.249,0.293,0.347,0.413,0.494,0.592,0.709,0.845,1.000,1.169,1.348,1.531,1.715,1.8896,2.071,2.240,2.402,2.557,2.706,2.848,2.983,3.114,3.239,3.359]
rigid_ρ = [(1+(htSpecRatio-1)/2*rigid_M[i]^2)^(-1/(htSpecRatio-1)) for i in 1:1:pointLength] 
rigid_p = [(1+(htSpecRatio-1)/2*rigid_M[i]^2)^(-htSpecRatio/(htSpecRatio-1)) for i in 1:1:pointLength] 
rigid_t = [(1+(htSpecRatio-1)/2*rigid_M[i]^2)^(-1) for i in 1:1:pointLength] 
graphMaker.makeTwoLineGraph(x7_3, totalResult7_3[totalStep7_3].M, rigid_M,"数値解析","厳密解","位置","マッハ数","7_3Mach")
graphMaker.makeTwoLineGraph(x7_3, totalResult7_3[totalStep7_3].ρ, rigid_ρ,"数値解析","厳密解","位置","密度","7_3Rho")
graphMaker.makeTwoLineGraph(x7_3, totalResult7_3[totalStep7_3].p, rigid_p,"数値解析","厳密解","位置","圧力","7_3Pressure")
graphMaker.makeTwoLineGraph(x7_3, totalResult7_3[totalStep7_3].t, rigid_t,"数値解析","厳密解","位置","温度","7_3Temperature")
graphMaker.makeMultiLineGraph(x7_3, "位置","圧力","7_3pressureConvergence",(totalResult7_3[1].p,"0Δt"),(totalResult7_3[400].p,"400Δt"),(totalResult7_3[800].p,"800Δt"),(totalResult7_3[1200].p,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_3, "位置","密度","7_3densityConvergence",(totalResult7_3[1].ρ,"0Δt"),(totalResult7_3[400].ρ,"400Δt"),(totalResult7_3[800].ρ,"800Δt"),(totalResult7_3[1200].ρ,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_3, "位置","温度","7_3temperatureConvergence",(totalResult7_3[1].t,"0Δt"),(totalResult7_3[400].t,"400Δt"),(totalResult7_3[800].t,"800Δt"),(totalResult7_3[1200].t,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_3, "位置","マッハ数","7_3machConvergence",(totalResult7_3[1].M,"0Δt"),(totalResult7_3[400].M,"400Δt"),(totalResult7_3[800].M,"800Δt"),(totalResult7_3[1200].M,"1200Δt"))

massFlow1 = totalResult7_3[1].ρ .*totalResult7_3[1].v .* A7_3
massFlow50 = totalResult7_3[50].ρ .*totalResult7_3[50].v .* A7_3
massFlow100 = totalResult7_3[100].ρ .*totalResult7_3[100].v .* A7_3
massFlow150 = totalResult7_3[150].ρ .*totalResult7_3[150].v .* A7_3
massFlow200 = totalResult7_3[200].ρ .*totalResult7_3[200].v .* A7_3
massFlow700 = totalResult7_3[700].ρ .*totalResult7_3[700].v .* A7_3
graphMaker.makeMultiLineGraph(x7_3,"位置","流量","7_3massFlow",(massFlow1,"0Δt"),(massFlow50,"50Δt"), (massFlow100,"100Δt"), (massFlow150,"150Δt"),(massFlow200,"200Δt"),(massFlow700,"700Δt"))
graphMaker.makeMultiLineLogGraph([i for i in 1:1:totalStep7_3],"ステップ数","残余","7_3decaying",([abs(totalResult7_3[i].aveDif.ave_dρ[16]) for i in 1:1:totalStep7_3],"|∂ρ/∂t|"),([abs(totalResult7_3[i].aveDif.ave_dv[16]) for i in 1:1:totalStep7_3],"|∂v/∂t|"))
massFlow700_courant_113 = totalResult7_3_limitCourant[700].ρ .*totalResult7_3_limitCourant[700].v .* A7_3
massFlow700_courant_1 = totalResult7_3_courant1[700].ρ .*totalResult7_3_courant1[700].v .* A7_3
graphMaker.makeMultiLineGraph(x7_3, "位置","流量","7_3courantComparison",(massFlow700,"C=0.5"),(massFlow700_courant_1,"C=1.0"),(massFlow700_courant_113,"C=1.13"))
#endregion

#region 7-4
specificBackPressure = 0.93

function area7_4(x::Vector{Float64})::Vector{Float64}
    a::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        a[i] = ifelse(x[i] <= 1.5,
         round(1+2.2*(x[i]-1.5)^2,digits=3),
         round(1+0.2223*(x[i]-1.5)^2,digits=3))
    end
    return a
end

function density7_4(x::Vector{Float64})::Vector{Float64}
    d::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        d[i] = round(1-0.023*x[i],digits=3)
    end
    return d
end

function temperature7_4(x::Vector{Float64})::Vector{Float64}
    t::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        t[i] = round(1-0.009333*x[i],digits=3)
    end
    return t
end

function velocity7_4(x::Vector{Float64}, t::Vector{Float64})::Vector{Float64}
    v::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        v[i] = round(0.05+0.11*x[i],digits=3)
    end
    return v
end


#配列はgeneratorで一気に作る.

x7_4 = [i for i in 0:deltax:ductLength]
ρ7_4 = density7_4(x7_4)
t7_4 = temperature7_4(x7_4)
v7_4 = velocity7_4(x7_4,t7_4)
A7_4 = area7_4(x7_4)


one7_4 = NonCnsvMacCormack.sub_subOneStep(deltax,
                            1,
                            pointLength,
                            ρ7_4,
                            v7_4,
                            t7_4,
                            A7_4,
                            specificBackPressure
                            )
totalResult7_4 = Vector{NonCnsvMacCormack.PrimitiveValuesWithStepNum}()
push!(totalResult7_4,one7_4)

totalStep7_4 = 1200
let
    prev = one7_4
    for i in 2:1:totalStep7_4
        next = NonCnsvMacCormack.sub_subOneStep(deltax,
                            i,
                            pointLength,
                            prev.ρ,
                            prev.v,
                            prev.t,
                            A7_4,
                            specificBackPressure
                            )
        push!(totalResult7_4,next)
        prev = next
    end
end

#=
see = totalStep7_4
println("回数"*string((totalResult7_4[see].step)))
println("圧力"*string(roundTo_3digits(totalResult7_4[see].p)))
println("速度"*string(roundTo_3digits(totalResult7_4[see].v)))
println("密度"*string(roundTo_3digits(totalResult7_4[see].ρ)))
println("温度"*string(roundTo_3digits(totalResult7_4[see].t)))
println("マッハ数"*string(roundTo_3digits(totalResult7_4[see].M)))
println("質量流量"*string(roundTo_3digits([totalResult7_4[see].ρ[i]*totalResult7_4[see].v[i]*area7_4(x7_4)[i] for i in 1:pointLength])))
=#

graphMaker.makeMultiLineGraph(x7_4, "位置","圧力","7_4pressureConvergence",(totalResult7_4[1].p,"0Δt"),(totalResult7_4[400].p,"400Δt"),(totalResult7_4[800].p,"800Δt"),(totalResult7_4[1200].p,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_4, "位置","密度","7_4densityConvergence",(totalResult7_4[1].ρ,"0Δt"),(totalResult7_4[400].ρ,"400Δt"),(totalResult7_4[800].ρ,"800Δt"),(totalResult7_4[1200].ρ,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_4, "位置","温度","7_4temperatureConvergence",(totalResult7_4[1].t,"0Δt"),(totalResult7_4[400].t,"400Δt"),(totalResult7_4[800].t,"800Δt"),(totalResult7_4[1200].t,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_4, "位置","マッハ数","7_4machConvergence",(totalResult7_4[1].M,"0Δt"),(totalResult7_4[400].M,"400Δt"),(totalResult7_4[800].M,"800Δt"),(totalResult7_4[1200].M,"1200Δt"))
#書くもの: 流量変化

#endregion

#region 7-5
specific_U2 = 0.59

function density7_5(x::Vector{Float64})::Vector{Float64}
    density::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        if x[i] <= 0.5
            density[i] = 1.0
        elseif 0.5 < x[i] <= 1.5
            density[i] = 1.0-0.366*(x[i]-0.5)
        else
            density[i] = 0.634-0.3879*(x[i]-1.5)
        end
    end
    return density
end

function temperature7_5(x::Vector{Float64})::Vector{Float64}
    temperature::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        if x[i] <= 0.5
            temperature[i] = 1.0
        elseif 0.5 < x[i] <= 1.5
            temperature[i] = 1.0-0.167*(x[i]-0.5)
        else
            temperature[i] = 0.833-0.3507*(x[i]-1.5)
        end
    end
    return temperature
end

function velocity7_5(specific_U2::Float64, area::Vector{Float64}, density::Vector{Float64})
    velocity = zeros(pointLength)
    for i in 1:1:pointLength
        velocity[i] = specific_U2/(area[i]*density[i])
    end
    return velocity
end


x7_5 = [i for i in 0:deltax:ductLength]
A7_5 = area7_3(x7_5)
ρ7_5 = density7_5(x7_5)
v7_5 = velocity7_5(specific_U2,A7_5,ρ7_5)
t7_5 = temperature7_5(x7_5)


one7_5 = CnsvMacCormack.sub_superOneStep(deltax,
                            1,
                            pointLength,
                            ρ7_5,
                            t7_5,
                            v7_5,
                            CnsvMacCormack.encodeTo_u1(ρ7_5,A7_5),
                            CnsvMacCormack.encodeTo_u2(ρ7_5,A7_5,v7_5),
                            CnsvMacCormack.encodeTo_u3(ρ7_5,A7_5,v7_5,t7_5),
                            A7_5
                            )
                            
totalResult7_5 = Vector{CnsvMacCormack.SolutionVectorsWithStepAndPrim}()
push!(totalResult7_5,one7_5)
totalStep7_5 = 1400

let
    prev = one7_5
    for i in 2:1:totalStep7_5
        next = CnsvMacCormack.sub_superOneStep(deltax,
                            i,
                            pointLength,
                            prev.primitive.ρ,
                            prev.primitive.t,
                            prev.primitive.v,
                            prev.u1,
                            prev.u2,
                            prev.u3,
                            A7_5)
        push!(totalResult7_5,next)
        prev = next
    end
end
see = totalStep7_5
#=
println("回数"*string((totalResult7_5[see].step)))
println("u1"*string(roundTo_3digits(totalResult7_5[see].u1)))
println("u2"*string(roundTo_3digits(totalResult7_5[see].u2)))
println("u3"*string(roundTo_3digits(totalResult7_5[see].u3)))
println("温度"*string(roundTo_3digits(totalResult7_5[see].primitive.t)))
println("密度"*string(roundTo_3digits(totalResult7_5[see].primitive.ρ)))
println("速度"*string(roundTo_3digits(totalResult7_5[see].primitive.v)))
println("マッハ"*string(roundTo_3digits(totalResult7_5[see].primitive.M)))
=#
graphMaker.makeMultiLineGraph(x7_5, "位置","圧力","7_5pressureConvergence",(totalResult7_5[1].primitive.p,"0Δt"),(totalResult7_5[400].primitive.p,"400Δt"),(totalResult7_5[800].primitive.p,"800Δt"),(totalResult7_5[1200].primitive.p,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_5, "位置","密度","7_5densityConvergence",(totalResult7_5[1].primitive.ρ,"0Δt"),(totalResult7_5[400].primitive.ρ,"400Δt"),(totalResult7_5[800].primitive.ρ,"800Δt"),(totalResult7_5[1200].primitive.ρ,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_5, "位置","温度","7_5temperatureConvergence",(totalResult7_5[1].primitive.t,"0Δt"),(totalResult7_5[400].primitive.t,"400Δt"),(totalResult7_5[800].primitive.t,"800Δt"),(totalResult7_5[1200].primitive.t,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_5, "位置","マッハ数","7_5machConvergence",(totalResult7_5[1].primitive.M,"0Δt"),(totalResult7_5[400].primitive.M,"400Δt"),(totalResult7_5[800].primitive.M,"800Δt"),(totalResult7_5[1200].primitive.M,"1200Δt"))
graphMaker.makeMultiLineGraph(x7_5,"位置","流量","7_5massFlowConvergence",(totalResult7_5[1].u2,"0Δt"),(totalResult7_5[50].u2,"50Δt"),(totalResult7_5[100].u2,"100Δt"),(totalResult7_5[150].u2,"150Δt"),(totalResult7_5[200].u2,"200Δt"),(totalResult7_5[700].u2,"700Δt"))
#正確値の表示のための設定
exact_u2 = [0.579 for i in 1:1:pointLength]
#保存非保存の比較
graphMaker.makeMultiLineGraph(x7_5, "位置","流量","7_5CnsvVsNon",(exact_u2,"理論値"), (totalResult7_5[totalStep7_5].u2,"保存形"),(totalResult7_3[totalStep7_3].ρ.*totalResult7_3[totalStep7_3].v .* A7_5,"非保存形"))

#endregion7-5

#region 7-6

specificPressure7_6 = 0.6784
specific_U2_7_6 = 0.59
Cx1 = 0.0
Cx2 = 0.1
Cx3 = 0.2
Cx4 = 0.3
totalStep7_6 = 1400

function density7_6(x::Vector{Float64})::Vector{Float64}
    density::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        if x[i] <= 0.5
            density[i] = 1.0
        elseif 0.5 < x[i] <= 1.5
            density[i] = 1.0-0.366*(x[i]-0.5)
        elseif 1.5 < x[i] <= 2.1
            density[i] = 0.634-0.702*(x[i]-1.5)
        else
            density[i] = 0.5892+0.10228*(x[i]-2.1)
        end
    end
    return density
end

function temperature7_6(x::Vector{Float64})::Vector{Float64}
    temperature::Vector{Float64} = zeros(pointLength)
    for i in 1:1:pointLength
        if x[i] <= 0.5
            temperature[i] = 1.0
        elseif 0.5 < x[i] <= 1.5
            temperature[i] = 1.0-0.167*(x[i]-0.5)
        elseif 1.5 < x[i] <= 2.1
            temperature[i] = 0.833-0.4908*(x[i]-1.5)
        else 
            temperature[i] = 0.93968+0.0622*(x[i]-2.1)
        end
    end
    return temperature
end

function velocity7_6(specific_U2::Float64, area::Vector{Float64}, density::Vector{Float64})
    velocity = zeros(pointLength)
    for i in 1:1:pointLength
        velocity[i] = specific_U2/(area[i]*density[i])
    end
    return velocity
end


x7_6 = [i for i in 0:deltax:ductLength]
A7_6 = area7_3(x7_6)
ρ7_6 = density7_6(x7_6)
v7_6 = velocity7_6(specific_U2_7_6,A7_6,ρ7_6)
t7_6 = temperature7_6(x7_6)
p7_6 = ρ7_6 .* t7_6

function sonicWaveSolutionWith(Cx::Float64)::Vector{CnsvMacCormack.SolutionVectorsWithStepAndPrim}
    one = CnsvMacCormack.sub_subOneStepWithArtificialViscosity(deltax,
                                1,
                                pointLength,
                                specificPressure7_6,
                                ρ7_6,
                                t7_6,
                                v7_6,
                                p7_6,
                                CnsvMacCormack.encodeTo_u1(ρ7_6,A7_6),
                                CnsvMacCormack.encodeTo_u2(ρ7_6,A7_6,v7_6),
                                CnsvMacCormack.encodeTo_u3(ρ7_6,A7_6,v7_6,t7_6),
                                A7_6,
                                Cx
                                )
                                
    totalResult = Vector{CnsvMacCormack.SolutionVectorsWithStepAndPrim}()
    push!(totalResult,one)

    let
        prev = one
        for i in 2:1:totalStep7_6
            next = CnsvMacCormack.sub_subOneStepWithArtificialViscosity(deltax,
                                i,
                                pointLength,
                                specificPressure7_6,
                                prev.primitive.ρ,
                                prev.primitive.t,
                                prev.primitive.v,
                                prev.primitive.p,
                                prev.u1,
                                prev.u2,
                                prev.u3,
                                A7_6,
                                Cx)
            push!(totalResult,next)
            prev = next
        end
    end
    return totalResult
end
totalResult7_6_Cx1 = sonicWaveSolutionWith(Cx1)
totalResult7_6_Cx2 = sonicWaveSolutionWith(Cx2)
totalResult7_6_Cx3 = sonicWaveSolutionWith(Cx3)
totalResult7_6_Cx4 = sonicWaveSolutionWith(Cx4)

#=
see = totalStep7_6
println("回数"*string((totalResult7_6[see].step)))
println("u1"*string(roundTo_3digits(totalResult7_6[see].u1)))
println("u2"*string(roundTo_3digits(totalResult7_6[see].u2)))
println("u3"*string(roundTo_3digits(totalResult7_6[see].u3)))
println("温度"*string(roundTo_3digits(totalResult7_6[see].primitive.t)))
println("密度"*string(roundTo_3digits(totalResult7_6[see].primitive.ρ)))
println("速度"*string(roundTo_3digits(totalResult7_6[see].primitive.v)))
println("圧力"*string(roundTo_3digits(totalResult7_6[see].primitive.p)))
println("マッハ"*string(roundTo_3digits(totalResult7_6[see].primitive.M)))
=#

graphMaker.makeMultiLineGraph(x7_6, "位置","圧力","7_6pressureBig",(totalResult7_6_Cx1[totalStep7_6].primitive.p,"Cx=0.0"),(totalResult7_6_Cx2[totalStep7_6].primitive.p,"Cx=0.1"),(totalResult7_6_Cx3[totalStep7_6].primitive.p,"Cx=0.2"),(totalResult7_6_Cx4[totalStep7_6].primitive.p,"Cx=0.3"))
graphMaker.makeMultiLineGraph(x7_6, "位置","流量","7_6massFlow",(totalResult7_6_Cx1[totalStep7_6].u2,"Cx=0.0"),(totalResult7_6_Cx2[totalStep7_6].u2,"Cx=0.1"),(totalResult7_6_Cx3[totalStep7_6].u2,"Cx=0.2"),(totalResult7_6_Cx4[totalStep7_6].u2,"Cx=0.3"))

#endregion 7-6

#region prob2
#endregion

#endregion

