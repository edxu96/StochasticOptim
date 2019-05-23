# DTU02435 Decision Making under Uncertainty
# Assignment 3, Task 2
# Out-of-Sample Test of the Plan
# Edward J. Xu
# May 5th, 2019
push!(LOAD_PATH, "$(homedir())/Desktop/StochasticAirplaneRent")
cd("$(homedir())/Desktop/StochasticAirplaneRent")
using JuMP
# using GLPKMathProgInterface
# using Gurobi
using CPLEX
using CSV
using DataFrames


function processData()
    datf_cost = CSV.read("Data/OperationCosts.csv", header = true, delim = ',')
    datf_cap = CSV.read("Data/PlaneCapacity.csv", header = ["p", "cap"], delim = ',')
    datf_route = CSV.read("Data/DemandsRoute.csv", header = true, delim = ',')
    datf_costAdmin = CSV.read("Data/CostAdquisition.csv", header = ["p", "cost"], delim = ',')
    datf_demand = CSV.read("Data/OutofSample.csv", header = true, delim = ',')
    # Calculate the demand matrix
    mat4_demand = zeros(4, 10, 12, 500)
    for i = 1: 4
        for j = 1: 10
            for w = 1: 500
                for t = 1: 12
                    mat4_demand[i, j, t, w] = datf_demand[t, w+1] * datf_route[i, j+1]
                end
            end
        end
    end
    return (datf_cost, datf_cap.cap, datf_costAdmin.cost, mat4_demand)
end


function optim(vec_y, datf_cost, vec_cap, vec_costAdmin, mat_demand)
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    # model = Model(solver = GurobiSolver(Presolve=0))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model, mat3_z[1:4, 1:10, 1:10] >= 0, Int)
    @variable(model, mat_m[1:4, 1:10] >= 0)
    @objective(model, Min, sum(datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j),
        (datf_cost[:p] .== p))][1] * mat3_z[i, j, p] for p = 1:10, i = 1:4, j = 1:10) + sum(mat_m) * 1000)
    @constraint(model, [i = 1:4, j = 1:10], mat_demand[i, j] - sum(mat3_z[i, j, p] * vec_cap[p] for p = 1:10)
        <= mat_m[i, j])
    @constraint(model, [p = 1:10], sum(mat3_z[i, j, p] for i = 1:4, j = 1:10) <= vec_y[p])
    solve(model)
    result_obj = getobjectivevalue(model)
    return result_obj
end


function testOutOfSample(vec_y, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    mat_result_obj = zeros(500, 12)
    vec_result_obj = zeros(500)
    ## 1,  Begin Optimization
    for w = 1:500
        for t = 1:12
            mat_result_obj[w, t] = optim(vec_y, datf_cost, vec_cap, vec_costAdmin, mat4_demand[:, :, t, w])
        end
        vec_result_obj[w] = sum(mat_result_obj[w, tt] for tt = 1:12) + sum(vec_y[p] * vec_costAdmin[p] for p = 1: 10)
        println("Result: obj = $(vec_result_obj[w]), when w = $(w).")
    end
    return vec_result_obj
end


function main()
    (datf_cost, vec_cap, vec_costAdmin, mat4_demand) = processData()
    println("################################################################################\n",
            "############### 1/4, Out-of-Sample Test of Deterministic Program ###############\n",
            "################################################################################")
    vec_y_1 = [0.0, 1.0, 2.0, 2.0, 0.0, 12.0, 5.0, 4.0, 8.0, 8.0]
    timeStart = time()
    vec_result_obj_1 = testOutOfSample(vec_y_1, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    println("Elapsed time: $(time()-timeStart) seconds.")
    CSV.write("result/out-of-sample_1.csv", DataFrame(vec_result_obj_1'), writeheader = false, delim = ',')
    println("################################################################################\n",
            "################# 2/4, Out-of-Sample Test of Stochastic Program ################\n",
            "################################################################################")
    vec_y_2 = [0.0, 7.0, 2.0, 3.0, 0.0, 12.0, 4.0, 4.0, 8.0, 8.0]
    timeStart = time()
    vec_result_obj_2 = testOutOfSample(vec_y_2, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    println("Elapsed time: $(time()-timeStart) seconds.")
    CSV.write("result/out-of-sample_2.csv", DataFrame(vec_result_obj_2'), writeheader = false, delim = ',')
    println("################################################################################\n",
            "################## 3/4, Out-of-Sample Test of Robust Program ###################\n",
            "################################################################################")
    vec_y_3 = [0.0, 3.0, 1.0, 2.0, 0.0, 12.0, 2.0, 8.0, 8.0, 8.0]
    timeStart = time()
    vec_result_obj_3 = testOutOfSample(vec_y_3, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    println("Elapsed time: $(time()-timeStart) seconds.")
    CSV.write("result/out-of-sample_3.csv", DataFrame(vec_result_obj_3'), writeheader = false, delim = ',')
    println("################################################################################\n",
            "################################### 4/4, End ###################################\n",
            "################################################################################")
end


main()
