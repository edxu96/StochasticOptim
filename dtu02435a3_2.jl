# DTU02435 Decision Making under Uncertainty
# Assignment 3, Task 2
# Solve by heuristics and Benders Algorithm / Gurobi
# Edward J. Xu
# April 27th, 2019
push!(LOAD_PATH, "$(homedir())/Desktop/optim")
cd("$(homedir())/Desktop/optim")
using JuMP
using GLPKMathProgInterface
using Gurobi
using CSV
using DataFrames
using BlackBoxOptim


function processData()
    datf_numMax = CSV.read("MaxPlanes.csv", header = ["p", "num"], delim = ',')
    datf_cost = CSV.read("OperationCosts2.csv", header = true, delim = ',')
    datf_cap = CSV.read("PlaneCapacity.csv", header = ["p", "cap"], delim = ',')
    datf_route = CSV.read("DemandsRoute.csv", header = true, delim = ',')
    datf_costAdmin = CSV.read("CostAdquisition.csv", header = ["p", "cost"], delim = ',')
    datf_prob = CSV.read("Probability.csv", header = ["prob"], delim = ',')
    #
    datf_demand = CSV.read("Scenario.csv", header = true, delim = ',')
    mat4_demand = zeros(4, 5, 12, 10)
    for i = 1: 4
        for j = 1: 5
            for w = 1: 10
                mat4_demand[i, j, :, w] = datf_demand[w] * datf_route[i, j+1]
            end
        end
    end
    return (datf_numMax.num, datf_cost, datf_cap.cap, datf_costAdmin.cost, mat4_demand, datf_prob.prob)
end


function solve_sub(vec_yBar, datf_cost, mat_demand, vec_cap)
    model_sub = Model(solver = GurobiSolver(Presolve = 0, OutputFlag = 0, gurobi_env))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model_sub, mat3_z[1: 4, 1: 5, 1: 10] >= 0, Int)
    @variable(model_sub, mat_m[1: 4, 1: 5] >= 0)
    @objective(model_sub, Min, sum([datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j),
        (datf_cost[:p] .== p))][1] * mat3_z[i, j, p] + mat_m[i, j] * 1000] for p = 1: 10, i = 1: 4,
        j = 1: 5)[1])
    @constraint(model_sub, [i = 1: 4, j = 1: 5], mat_demand[i, j] -
        sum([mat3_z[i, j, p] * vec_cap[p]] for p = 1: 10)[1] <= mat_m[i, j])
    @constraint(model_sub, [p = 1: 10], sum([mat3_z[i, j, p]] for i = 1: 4, j = 1: 5)[1] <=
        vec_yBar[p])
    solve(model_sub)
    obj = getobjectivevalue(model_sub)
    vec_result_z = getvalue(mat3_z)
    vec_result_m = getvalue(mat_m)
    return (obj, vec_result_z, vec_result_m)
end


function solve_sub_all(vec_yBarRaw; datf_cost = datf_cost, vec_cap = vec_cap, vec_costAdmin = vec_costAdmin,
        mat4_demand = mat4_demand, vec_prob = vec_prob)
# function solve_sub_all(vec_yBarRaw)
    vec_yBar = round.(vec_yBarRaw)
    println("--------------------------------------------------------------------------------\n",
            "------------------------------ Solve Sub-Problem -------------------------------\n",
            "Input: vec_yBarRaw = $(vec_yBarRaw)\n",
            "Input: vec_yBar = $(vec_yBar)\n",
            "--------------------------------------------------------------------------------\n")
    mat_obj = zeros(10, 12)
    for w = 1: 10
        for t = 1: 12
            println("Input: w = $(w), t = $(t), step = $(12 * (w - 1) + t) / 120 ")
            (mat_obj[w, t], mat_result_z, mat_result_m) = solve_sub(
                vec_yBar, datf_cost, mat4_demand[:, :, t, w], vec_cap)
        end
    end
    obj = sum([vec_yBar[p] * vec_costAdmin[p] + vec_prob[w] * mat_obj[w, t]] for p = 1: 10, w = 1: 5, t = 1: 12)[1]
    println("--------------------------------------------------------------------------------\n",
            "Result: obj = $(obj)\n",
            "--------------------------------------------------------------------------------\n")
    return obj
end


function main()
    # Begin Optimization
    println("################################################################################\n",
            "############################ 1/2. Begin Optimization ###########################\n",
            "################################################################################\n")
    mod_master = bboptimize(solve_sub_all; SearchRange = [(0.0, 12.0), (0.0, 12.0), (0.0, 12.0), (0.0, 12.0),
        (0.0, 12.0), (0.0, 12.0), (0.0, 12.0), (0.0, 8.0), (0.0, 8.0), (0.0, 8.0)], NumDimensions = 10,
        MaxTime = false)
    println("################################################################################\n",
            "############################## 2/2. Nominal Ending #############################\n",
            "################################################################################\n")
    # Print the Result
    # println("Result: obj = $(obj).")
    # println("Result: vec_y = $(vec_result_y).")
    # println("Result: mat5_z = $(vec_result_z).")
    # println("Result: mat4_m = $(vec_result_m).")
    # return (vec_result_z, vec_result_m)
end

const gurobi_env = Gurobi.Env()
const (vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob) = processData()
main()
