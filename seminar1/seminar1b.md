# seminar_1b
Santina  
January 11, 2017  


# Basic inspection of my data

## smaller heading

### even smaller

I want to know how many rows are in my data using R. 

First I read my data: 

```r
prDat <- read.table("GSE4051_MINI.tsv", header = TRUE, row.names = 1)
```


Then I count how many rows are in my data: 

```r
nrow(prDat)
```

```
## [1] 39
```

I have 39 in my data. Notice how I used inline R code (you can see it by inspecting my .Rmd).

I also want to plot something just to show you how it works. 


```r
plot(cars)
```

![](seminar1b_files/figure-html/car_plot-1.png)<!-- -->

Becauese I named the above R chunk, the figure generated from it is not "unnamed..." but it's named as car_plot-1.png
