module graphMaker
    using Plots
    using PyPlot
    using PyCall
    default(show=false)
    @pyimport matplotlib.animation as anim

    function makeTwoLineGraph(xValue::Vector{Float64},
         yValue1::Vector{Float64},
         yValue2::Vector{Float64},
         yLabel1::String,
         yLabel2::String,
         xAxisLabel::String,
         yAxisLabel::String,
         saveName::String
         )
        fig = figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(xValue, yValue1, label=yLabel1)
        ax.plot(xValue, yValue2, label=yLabel2)
        ax.set_xlim(xValue[1],xValue[end])
        ax.legend()
        ax.set_xlabel(xAxisLabel)
        ax.set_ylabel(yAxisLabel)
        fig.savefig("./graph/"*saveName*".png")
        fig.clf()
        close(fig)
    end
    function makeTwoLineGraph(xValue::Vector{Float64},
         yValue1::Vector{Float64},
         yValue2::Vector{Float64},
         yLabel1::String,
         yLabel2::String,
         xAxisLabel::String,
         yAxisLabel::String,
         saveName::String
         )
        fig = figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(xValue, yValue1, label=yLabel1)
        ax.plot(xValue, yValue2, label=yLabel2)
        ax.set_xlim(xValue[1],xValue[end])
        ax.legend()
        ax.set_xlabel(xAxisLabel)
        ax.set_ylabel(yAxisLabel)
        fig.savefig("./graph/"*saveName*".png")
        fig.clf()
        close(fig)
    end
    
    function makeMultiLineGraph(xValue::Vector{Float64},
         xlabel::String,
         ylabel::String,
         saveName::String,
         yValueWithLabel...
         )
        fig = figure()
        ax = fig.add_subplot(1, 1, 1)
        for i in 1:1:length(yValueWithLabel)
            ax.plot(xValue, yValueWithLabel[i][1], label=yValueWithLabel[i][2])
        end
        ax.set_xlim(xValue[1],xValue[end])
        ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.savefig("./graph/"*saveName*".png")
        fig.clf()
        close(fig)
    end
    function makeMultiLineLogGraph(xValue::Vector{Int64},
         xlabel::String,
         ylabel::String,
         saveName::String,
         yValueWithLabel...
         )
        fig = figure()
        ax = fig.add_subplot(1, 1, 1)
        for i in 1:1:length(yValueWithLabel)
            ax.plot(xValue, yValueWithLabel[i][1], label=yValueWithLabel[i][2])
        end
        ax.set_xlim(xValue[1],xValue[end])
        ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        fig.savefig("./graph/"*saveName*".png")
        fig.clf()
        close(fig)
    end
    #整理する時間がなかったのでここはもうめちゃくちゃ yValueWithLabelはSolutionVectorsWithStepAndPrim型
    function makeShockAnimation(xValue::Vector{Float64}, totalStep::Int64,name::String , yValueWithLabel)
        anim = Plots.Animation()
        anim2 = Plots.Animation()
        anim3 = Plots.Animation()
        tmp = 0.0
        for i=1:totalStep
            tmp += yValueWithLabel[i].timeStep
            plt = Plots.plot(xValue, yValueWithLabel[i].primitive.v,label="velocity",xlabel="x m",ylabel="u m/s",title="time: "*string(tmp )*" s")
            Plots.frame(anim, plt)

            plt2 = Plots.plot(xValue, yValueWithLabel[i].primitive.ρ,label="density",xlabel="x m",ylabel="ρ kg/m3",title="time: "*string(tmp )*" s")
            Plots.frame(anim2, plt2)

            plt3 = Plots.plot(xValue, yValueWithLabel[i].primitive.p,label="pressure",xlabel="x m",ylabel="p Pa",title="time: "*string(tmp )*" s")
            Plots.frame(anim3, plt3)
        end
        
        Plots.gif(anim, "./gif/"*name*"-velo"*".gif", fps = 30)
        Plots.gif(anim2, "./gif/"*name*"-rho"*".gif", fps = 30)
        Plots.gif(anim3, "./gif/"*name*"-pres"*".gif", fps = 30)
    end
end