export Poisson_Regression, predict, fit

type poissonRegression
    #array of weights
    w
    #w0 constnat
    w_0::Float64
    #learning rate
    eta::Float64
    #decay rate for the learning rate 
    alpha::Float64
    #number of samples parsed through (will be incremnted automatically)
    n::Int64
    #number of features
    d::Int64
end

#constructor for convenience
function Poisson_Regression(;eta::Float64=0.0001, alpha::Float64=1.00, d::Int64=0)
    if d != 0
        model = poissonRegression(zeros(d), 0.0, eta, alpha, 0, d)
    else
        model = poissonRegression(nothing, 0.0, eta, alpha, 0, d)
    end
    
    #return the model
    model
end

#prediction using n examles (nxd matrix)
function predict(model::poissonRegression, x::Array{Float64,2})
    linear_prediction = x*model.w+model.w_0
    prediction = exp(linear_prediction)
end

#prediction using 1xd array 
function predict(model::poissonRegression, x::Array{Float64,1})
    linear_prediction = x*model.w+model.w_0
    prediction = exp(linear_prediction)
end

#takes in nxd batch data as an input, conducts stochastic gradient descent
function fit(model::poissonRegression, y::Array{Float64, 1}, x::Array{Float64, 2})
    #checking if y and x match in size
    if length(y)!=size(x)[1]
        return nothing
    end 
    
    #initiating weight array if necessary
    if model.d == 0
        model.d = size(x)[2]
        model.w = zeros(model.d)
    end
    
    #updating info (right now its just the number of examples parsed)
    model.n += length(y)
    num_data = length(y)
    
    #making prediction
    prediction = predict(model, x)
    
    #updating w_0
    model.w_0 = model.w_0 + model.eta*(sum(y-prediction, 1)[1]*1.0/num_data)
    
    #updating w
    model.w = model.w + model.eta*((transpose(x)*(y-prediction))/num_data)
    
    model
end 