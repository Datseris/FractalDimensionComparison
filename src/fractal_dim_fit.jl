export linear_regression_fit_linalg
export linear_regression_fit_lsqfit
export logarithmic_corrected_fit_lsqfit
using Statistics: std
using Distributions: TDist, cdf
using FractalDimensions: linreg

#########################################################################################
# Simplest linear fit via linear algebra
#########################################################################################
function linear_regression_fit_linalg(x, y)
    a, s = linreg(x, y)
    n = length(y)
    # CI computed via https://stattrek.com/regression/slope-confidence-interval
    # standard error of slope
    df = max(n - 2, 1)
    yhat = @. a + s*x
    standard_error = sqrt((sum((y .- yhat).^2) ./ df)) / sqrt(sum((x .- mean(x)).^2))
    ci = 0.95 # 95% confidence interval
    α = 1 - ci
    pstar = 1 - α/2
    tdist = TDist(df)
    critical_value = quantile(tdist, pstar)
    margin_of_error = critical_value * standard_error
    s05 = s - margin_of_error
    s95 = s + margin_of_error
    return s, s05, s95
end

#########################################################################################
# Linear fit with GLM
#########################################################################################
# This code does exactly the same thing as above, because it also calls the
# `TDist` with same quantile processing etc. It is kept here as legacy

# import GLM

function linear_regression_fit_glm(x, y)
    # `ones` is used here to obtain the value of the intercept
    X = hcat(ones(length(x)), x)
    out = GLM.lm(X, y)
    o, s = GLM.coef(out)
    # We have try-block, because if inappropriate slope is detected
    # by `linear_region`, with too few points, the `confint` fails.
    try
        s05, s95 = GLM.confint(out)[2, :]
        return s, s05, s95
    catch err
        @warn "GLM fit failed with error: "*string(err)
        return s, NaN, NaN
    end
end

#########################################################################################
# Linear fit with LsqFit
#########################################################################################
import LsqFit

function linear_regression_fit_lsqfit(x,y, p0 = [0.0, 0.0])
    linearregression(x, p) = @. p[1] + p[2]*x
    fit = LsqFit.curve_fit(linearregression, x, y, p0)
    s = LsqFit.coef(fit)[2]
    s05, s95 = LsqFit.confidence_interval(fit, 0.05)[2]
    return s, s05, s95
end

#########################################################################################
# Logarithmically corrected fit with LsqFit (Sprott & Rownalds correction)
#########################################################################################
function logarithmic_corrected_fit_lsqfit(x, y, p0 = [0.0, 0.0, 0.0])
    corrected_fit(x, p) = @. p[1] + p[2]*x + p[3]*log(-x)
    # TODO: Simply reduce ε to start at 0...
    i = findlast(<(0), x)
    if (!isnothing(i) && i < length(x)÷2)
        @info "More than half of the data have log(ε) > 0 for corrected version. Returning standard linear..."
        return linear_regression_fit_linalg(x,y)
    elseif !isnothing(i)
        x2 = x[1:i]
        y2 = y[1:i]
    else # last case: all x are negative
        x2 = x
        y2 = y
    end
    if !all(x2 .< 0)
        @info "There still exist positive log(ε) for corrected version. Using standard linear..."
        return linear_regression_fit_linalg(x,y)
    end
    try
        fit = LsqFit.curve_fit(corrected_fit, x2, y2, p0)
        s = LsqFit.coef(fit)[2]
        s05, s95 = LsqFit.confidence_interval(fit, 0.05)[2]
        return s, s05, s95
    catch err
        @info "Corrected version errored with $(err). Using standard linear..."
        return linear_regression_fit_linalg(x,y)
    end
end