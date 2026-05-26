# ideal-day-oa
Shiny app to investigate personalised activity composition changes in knee osteoarthritis pain intensity and interference

Please see the below paper(s) for more details and context:

* "[My Best Day for My Knee OA: personalised 24-hour recommendations for people with knee osteoarthritis](https://tobedetermined.info)". MacIntyre et al. (2026). *Journal TBD.*


## Shiny app in action (no installation required)

Hosted at the web location:

* [https://arena2024.shinyapps.io/Best-Day-OA/](https://arena2024.shinyapps.io/Best-Day-OA/)


## Running the Shiny App locally (on your own computer, installation required)

The following steps will allow you to run the `ideal-day-oa` Shiny App on your own computer:

* Make sure you have [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/) installed.
* Download and unzip this Github repository (green `<> Code` button in top right -> `Download ZIP`).
* Double click the `ideal-day-oa.Rproj` file - this should open an RStudio session.
* (only has to be performed once) Make sure the prerequisite packages are installed by running the below command in the R console:
```r
install.packages(
  c(
    "shiny", "shinydashboard", "shinyjs", "bslib",
    "ggplot2", "ggthemes", "viridis" , "plotly"
    "forcats", "dplyr", "lubridate", "readr", "tidyr", "forcats", "purrr",
    "compositions", "expm", "mosaic", "knitr"
  )
)
```
* Now to finally run the Shiny App by running the below command in the R console:
```r
library(shiny)
runApp()
```



## Example screenshots

... To be added shortly ...





