using VegaLite, DataFrames

# Load the MPG dataset
mpg = dataset("ggplot2", "mpg")

# Define the Vega-Lite specification for the plot
spec = @vlplot(
    data = mpg,
    mark = "point",
    encoding = {
        x = {"field" = "displ", "type" = "quantitative"},
        y = {"field" = "hwy", "type" = "quantitative"},
        color = {"field" = "class", "type" = "nominal"}
    }
)

# Embed the plot in the Markdown document using Vega-Embed
vega_embed(spec)