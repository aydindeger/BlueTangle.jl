using Plots

plot_measurement(mVector::Measurement)=plot_sample2([mVector],[""])

"""
    plot_measurement(mVector::Vector{Measurement}, labels::Vector{String} = [])

Plot the outcome probabilities for quantum measurements in a bar chart using Plots.jl.
Each measurement should include fields `bitstr` and `sample`, where `bitstr` is an array of possible outcomes
and `sample` is an array of corresponding probabilities.

# Arguments
- `mVector`: Vector of `Measurement` objects. Each `Measurement` should have fields `bitstr`, `sample`, and optionally `circuit_name`.
- `labels`: Optional vector of strings to label each set of measurements. If empty, labels are generated based on the index.

# Notes
This function limits to a maximum of five measurements for visual clarity.
"""
function plot_measurement(mVector::Vector{Measurement}, labels::Vector{String} = [""])
    # Limit the number of measurements to plot
    if length(mVector) > 5
        throw(ArgumentError("You can only plot up to five measurements."))
    end

    # Initialise the plot
    plt = plot(size=(800, 600), legend=:topright, framestyle=:box)

    # Define a colour palette
    # colors = [:blue, :orange, :green, :red, :purple, :brown]

    # Loop through each measurement to plot
    for (i, measurement) in enumerate(mVector)
        x_outcome = measurement.bitstr
        y_sample = measurement.sample

        # Generate or use provided labels
        if isempty(labels)
            label = "Measurement $(i)"
        else
            label = labels[i]
        end

        # Bar plot for each measurement
        bar!(x_outcome, y_sample, label=label, alpha=.7)#, width=0.8/i)
    end

    # Setting plot labels and title
    xlabel!("Outcomes")
    ylabel!("Probabilities")
    title!("Outcome Probabilities")

    # Display the plot
    display(plt)
end


