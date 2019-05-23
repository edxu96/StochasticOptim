# DTU02435 Decision Making under Uncertainty
# Assignment 3, Task 2
# Edward J. Xu
# May 5th, 2019
push!(LOAD_PATH, "$(homedir())/Desktop/StochasticAirplaneRent")
cd("$(homedir())/Desktop/StochasticAirplaneRent")
using JuMP
# using GLPKMathProgInterface
# using Gurobi
using CPLEX
using CSV


function processData()
    datf_numMax = CSV.read("Data/MaxPlanes.csv", header = ["p", "num"], delim = ',')
    datf_cost = CSV.read("Data/OperationCosts.csv", header = true, delim = ',')
    datf_cap = CSV.read("Data/PlaneCapacity.csv", header = ["p", "cap"], delim = ',')
    datf_route = CSV.read("Data/DemandsRoute.csv", header = true, delim = ',')
    datf_costAdmin = CSV.read("Data/CostAdquisition.csv", header = ["p", "cost"], delim = ',')
    datf_prob = CSV.read("Data/Probability.csv", header = ["prob"], delim = ',')
    datf_demand = CSV.read("Data/Scenario.csv", header = true, delim = ',')
    datf_expect = CSV.read("Data/ExpectedValue.csv", header = ["value"], delim = ',')
    # Calculate the demand matrix
    mat4_demand = zeros(4, 10, 12, 10)
    for i = 1: 4
        for j = 1: 10
            for w = 1: 10
                for t = 1: 12
                    mat4_demand[i, j, t, w] = 1000 * datf_demand[t, w] * datf_route[i, j+1]
                end
            end
        end
    end
    mat3_demandExpect = zeros(4, 10, 12)
    for i = 1: 4
        for j = 1: 10
            for t = 1: 12
                mat3_demandExpect[i, j, t] = 1000 * datf_expect.value[t] * datf_route[i, j+1]
            end
        end
    end
    return (datf_numMax.num, datf_cost, datf_cap.cap, datf_costAdmin.cost, mat4_demand,
            datf_prob.prob, mat3_demandExpect)
end


function solvePrintResult(model, variable_y, variable_m)
    println("---------------------------- 1/2. Begin Optimization ---------------------------")
    timeStart = time()
    solve(model)
    println("Elapsed time: $(time()-timeStart) seconds.")
    ## 2,  Print the Result
    objResult = getobjectivevalue(model)
    vec_yResult = getvalue(variable_y)
    vec_mResult = getvalue(variable_m)
    println("Result: obj = $(objResult)\n",
            "Result: vec_y = $(vec_yResult)\n",
            "Result: total unfulfilled passenger demand = $(sum(vec_mResult))\n",
            "---------------------------------- 2/2. Ending ---------------------------------")
    return vec_yResult
end


function optim_determin(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat3_demandExpect)
    ## 1,  Begin Optimization
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    # model = Model(solver = GurobiSolver(Presolve=0))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model, vec_y[1: 10] >= 0, Int)
    @variable(model, mat4_z[1: 4, 1: 10, 1: 10, 1: 12] >= 0, Int)
    @variable(model, mat3_m[1: 4, 1: 10, 1: 12] >= 0)
    @objective(model, Min, sum(vec_y[p] * vec_costAdmin[p] for p = 1: 10) +
        sum(datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j), (datf_cost[:p] .== p))][1] *
        mat4_z[i, j, p, t] for p = 1: 10, i = 1: 4, j = 1: 10, t = 1: 12) +
        sum(mat3_m[i, j, t] * 1000 for i = 1: 4, j = 1: 10, t = 1: 12))
    @constraint(model, [p = 1: 10], vec_y[p] <= vec_numMax[p])
    @constraint(model, [i = 1: 4, j = 1: 10, t = 1: 12], mat3_demandExpect[i, j, t] -
        sum([mat4_z[i, j, p, t] * vec_cap[p]] for p = 1: 10)[1] <= mat3_m[i, j, t])
    @constraint(model, [p = 1: 10, t = 1: 12], sum(mat4_z[i, j, p, t] for i = 1: 4, j = 1: 10) <= vec_y[p])
    vec_yResult = solvePrintResult(model, vec_y, mat3_m)
    return vec_yResult
end


function optim(vec_y, datf_cost, vec_cap, vec_costAdmin, mat_demand)
    ## Most basic function to optimize the arrangement of plane for any month in any scenario, with fixed fleet.
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
    objResult = getobjectivevalue(model)
    vec_mResult = getvalue(mat_m)
    return (objResult, sum(vec_mResult))
end


function optim_determinStochas(vec_yResult, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    # Function to calculate the obj when first stage variable is fixed by solutions from deterministic program, using
    # stochastic model. The program is decomposed into programs for different months.
    mat_objResult = zeros(10, 12)
    vec_objResult = zeros(10)
    mat_mResult = zeros(10, 12)
    vec_mResult = zeros(10)
    ## 1,  Begin Optimization
    println("---------------------------- 1/2. Begin Calculation ----------------------------\n")
    timeStart = time()
    for w = 1:10
        for t = 1:12
            (mat_objResult[w, t], mat_mResult[w, t]) = optim(vec_yResult, datf_cost, vec_cap, vec_costAdmin,
                mat4_demand[:, :, t, w])
        end
        vec_objResult[w] = sum(mat_objResult[w, tt] for tt = 1:12) +
            sum(vec_yResult[p] * vec_costAdmin[p] for p = 1: 10)
        vec_mResult[w] = sum(mat_mResult[w, tt] for tt = 1:12)
    end
    objResult = sum(vec_objResult[w] * vec_prob[w] for w = 1:10)
    mResult = sum(vec_mResult)
    println("Elapsed time: $(time()-timeStart) seconds.\n",
            "Result: obj = $(objResult)\n",
            "Result: vec_y = $(vec_yResult)\n",
            "Result: total unfulfilled passenger demand = $(mResult)\n",
            "---------------------------------- 2/2. Ending ---------------------------------")
end


function optim_stochas(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    ## 1,  Begin Optimization
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    # model = Model(solver = GurobiSolver(Presolve=0))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model, vec_y[1: 10] >= 0, Int)
    @variable(model, mat5_z[1: 4, 1: 10, 1: 10, 1: 12, 1: 10] >= 0, Int)
    @variable(model, mat4_m[1: 4, 1: 10, 1: 12, 1: 10] >= 0)
    @objective(model, Min, sum(vec_y[p] * vec_costAdmin[p] for p = 1: 10) +
        sum((datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j), (datf_cost[:p] .== p))][1] *
        mat5_z[i, j, p, t, w] + mat4_m[i, j, t, w] * 1000) * vec_prob[w] for w = 1: 10, p = 1: 10, i = 1: 4,
        j = 1: 10, t = 1: 12))
    @constraint(model, [p = 1: 10], vec_y[p] <= vec_numMax[p])
    @constraint(model, [i = 1: 4, j = 1: 10, t = 1: 12, w = 1: 10], mat4_demand[i, j, t, w] -
        sum([mat5_z[i, j, p, t, w] * vec_cap[p]] for p = 1: 10)[1] <= mat4_m[i, j, t, w])
    @constraint(model, [p = 1: 10, t = 1: 12, w = 1: 10], sum(mat5_z[i, j, p, t, w] for i = 1: 4, j = 1: 10) <=
        vec_y[p])
    solvePrintResult(model, vec_y, mat4_m)
end


# function optim_stochas_2(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
#     ## 1,  Begin Optimization
#     model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
#     # model = Model(solver = GurobiSolver(Presolve=0))
#     # model = Model(solver = GLPKSolverMIP())
#     @variable(model, mat_y[1: 10, 1: 10] >= 0, Int)
#     @variable(model, mat5_z[1: 4, 1: 10, 1: 10, 1: 12, 1: 10] >= 0, Int)
#     @variable(model, mat4_m[1: 4, 1: 10, 1: 12, 1: 10] >= 0)
#     @objective(model, Min, sum(mat_y[p, w] * vec_costAdmin[p] +
#         (datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j), (datf_cost[:p] .== p))][1] *
#         mat5_z[i, j, p, t, w] + mat4_m[i, j, t, w] * 1000) * vec_prob[w] for w = 1: 10, p = 1: 10, i = 1: 4,
#         j = 1: 10, t = 1: 12))
#     @constraint(model, [p = 1:10, w = 1:10], mat_y[p, w] <= vec_numMax[p])
#     @constraint(model, [i = 1:4, j = 1:10, t = 1:12, w = 1:10], mat4_demand[i, j, t, w] -
#         sum([mat5_z[i, j, p, t, w] * vec_cap[p]] for p = 1:10)[1] <= mat4_m[i, j, t, w])
#     @constraint(model, [p = 1:10, t = 1: 12, w = 1: 10], sum(mat5_z[i, j, p, t, w] for i = 1: 4, j = 1: 10) <=
#         mat_y[p, w])
#     @constraint(model, [p = 1:10, w = 2:10], mat_y[p, 1] == mat_y[p, w])
#     solvePrintResult(model, mat_y, mat4_m)
# end


function optim_waitSeeSub(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat3_demand)
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    # model = Model(solver = GurobiSolver(Presolve=0))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model, vec_y[1:10] >= 0, Int)
    @variable(model, mat4_z[1:4, 1:10, 1:10, 1:12] >= 0, Int)
    @variable(model, mat3_m[1:4, 1:10, 1:12] >= 0)
    @objective(model, Min, sum(vec_y[p] * vec_costAdmin[p] for p = 1:10) +
        sum(datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j), (datf_cost[:p] .== p))][1] *
        mat4_z[i, j, p, t] + mat3_m[i, j, t] * 1000 for p = 1: 10, i = 1: 4, j = 1: 10, t = 1: 12))
    @constraint(model, [p = 1:10], vec_y[p] <= vec_numMax[p])
    @constraint(model, [i = 1:4, j = 1:10, t = 1:12], mat3_demand[i, j, t] -
        sum(mat4_z[i, j, p, t] * vec_cap[p] for p = 1:10) <= mat3_m[i, j, t])
    @constraint(model, [p = 1: 10, t = 1: 12], sum(mat4_z[i, j, p, t] for i = 1: 4, j = 1: 10) <= vec_y[p])
    solve(model)
    objResult = getobjectivevalue(model)
    vec_yResult = getvalue(vec_y)
    vec_mResult = getvalue(mat3_m)
    println("Result: obj = $(objResult)\n",
            "Result: vec_y = $(vec_yResult)\n",
            "Result: unfulfilled passenger demand = $(sum(vec_mResult))")
    return objResult
end


function optim_waitSee(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    ## 1,  Begin Optimization
    vec_objResult = zeros(10)
    for w = 1:10
        println("----- $(w)-th Sub-Problem -----")
        vec_objResult[w] = optim_waitSeeSub(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand[:,:,:,w])
    end
    objWeightedSum = sum(vec_prob[i] .* vec_objResult[i] for i = 1:10)
    println("----- Final Result -----\n",
            "Result: weighted sum of obj = $(objWeightedSum)")
end


function optim_robust(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    ## 1,  Begin Optimization
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0))
    # model = Model(solver = GurobiSolver(Presolve=0))
    # model = Model(solver = GLPKSolverMIP())
    @variable(model, vec_y[1: 10] >= 0, Int)
    @variable(model, mat5_z[1: 4, 1: 10, 1: 10, 1: 12, 1: 10] >= 0, Int)
    @variable(model, mat4_m[1: 4, 1: 10, 1: 12, 1: 10] >= 0)
    @variable(model, alpha >= 0)
    @objective(model, Min, sum(vec_y[p] * vec_costAdmin[p] for p = 1: 10) + alpha)
    @constraint(model, [w = 1: 10],
        sum((datf_cost.cost[.&((datf_cost[:i] .== i), (datf_cost[:j] .== j), (datf_cost[:p] .== p))][1] *
        mat5_z[i, j, p, t, w] + mat4_m[i, j, t, w] * 1000) for p = 1: 10, i = 1: 4, j = 1: 10, t = 1: 12) <= alpha)
    @constraint(model, [p = 1: 10], vec_y[p] <= vec_numMax[p])
    @constraint(model, [i = 1: 4, j = 1: 10, t = 1: 12, w = 1: 10], mat4_demand[i, j, t, w] -
        sum([mat5_z[i, j, p, t, w] * vec_cap[p]] for p = 1: 10)[1] <= mat4_m[i, j, t, w])
    @constraint(model, [p = 1: 10, t = 1: 12, w = 1: 10], sum(mat5_z[i, j, p, t, w] for i = 1: 4, j = 1: 10) <=
        vec_y[p])
    solvePrintResult(model, vec_y, mat4_m)
end


function main()
    (vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob, mat3_demandExpect) = processData()
    println("################################################################################\n",
            "######################## 1/5, Deterministic Programming ########################\n",
            "################################################################################")
    vec_yResult = optim_determin(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat3_demandExpect)
    # vec_yResult = [0.0, 1.0, 2.0, 2.0, 0.0, 12.0, 5.0, 4.0, 8.0, 8.0]
    # optim_determinStochas(vec_yResult, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    # println("################################################################################\n",
    #         "########################## 2/5, Stochastic Programming #########################\n",
    #         "################################################################################")
    # optim_stochas(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    # println("################################################################################\n",
    #         "########################### 3/5, Wait-See Programming ##########################\n",
    #         "################################################################################")
    # optim_waitSee(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand, vec_prob)
    # println("################################################################################\n",
    #         "########################### 4/5, Robust Programming ############################\n",
    #         "################################################################################")
    # optim_robust(vec_numMax, datf_cost, vec_cap, vec_costAdmin, mat4_demand)
    println("################################################################################\n",
            "################################### 5/5, End ###################################\n",
            "################################################################################")
end


main()
