module graphMaker
    using PyPlot
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
end