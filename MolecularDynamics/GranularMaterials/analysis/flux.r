source('lib.r')

# Store with key=angle
paths = list()
timeWindows = list()


# Hole width = 1.5d -------------------------------------------------------
fpathTemplate = '../dumps/hw=1.5d grain_h=20 g=2.0 damp=-10 dt=0.01/angle=%s/grain_count.csv'  # note angle is a string

# xlim <- range( c( x1, x2)) 
# ylim <- range( c( y1, y2)) 
xlim = c(0, 100)
ylim = c(0, 150)

angle = '15'
paths[[angle]] = sprintf(fpathTemplate, angle)
timeWindows[[angle]] = c(0, 55)

angle = '30'
paths[[angle]] = sprintf(fpathTemplate, angle)
timeWindows[[angle]] = c(0, 65)

angle = '45'
timeWindows[[angle]] = c(0, 80)
paths[[angle]] = sprintf(fpathTemplate, angle)

angle = '60'
timeWindows[[angle]] = c(0, 80)
paths[[angle]] = sprintf(fpathTemplate, angle)


# Hole width = 3d ---------------------------------------------------------
fpathTemplate = '../dumps/hw=3.0d grain_h=20 g=2.0 damp=-10 dt=0.01/angle=%s/grain_count.csv'

xlim = c(0, 35)
ylim = c(0, 155)

angle = '15'
paths[[angle]] = sprintf(fpathTemplate, angle)
timeWindows[[angle]] = c(0, 20)

angle = '30'
paths[[angle]] = sprintf(fpathTemplate, angle)
timeWindows[[angle]] = c(5, 24)

angle = '45'
timeWindows[[angle]] = c(5, 24)
paths[[angle]] = sprintf(fpathTemplate, angle)

angle = '60'
timeWindows[[angle]] = c(0, 100)
paths[[angle]] = sprintf(fpathTemplate, angle)



# Plot margins etc
# b = 3.5; t = 2
# l = b; r = t
# mgp is how much margin in titles/labels
par(mgp=c(1.5, 0.5, 0))#, mar=c(b, l, t, r))


getColor = function (angle, alpha) {
    # Don't specify alpha if you wan't opaque
    iAngle = as.integer(angle)
#     myColor = colors()[[iAngle]]
    if (!missing(alpha)) alpha = 90 * alpha  # specify as fraction
    return(rgb(90-iAngle, 0, iAngle, alpha, maxColorValue=90))
}

plot(0, 0, type ="n", xlim = xlim, ylim = ylim,
     xlab='Time', ylab='Cumulative grains through aperture', 
     main='Accumulation of grains vs time') 

ix = 1
legendColors = c()
slopeLabels = c()
for (ix in 1:length(paths)) {
    angle = names(paths)[[ix]]  # as string
    fpath = paths[[ix]]
    window = timeWindows[[ix]]
    
    data = read.csv(fpath, header=T)
    windowData = subset(data, subset= (time > window[[1]]) & (time < window[[2]]))
    maxTime = max(data$time)
    minTime = min(data$time)
    
    myColor = getColor(angle)
    legendColors = c(legendColors, myColor)
    lines(data, col=myColor, cex=0.2) # cex=0.1 is for small points
    
    cat(paste(angle, 'degrees:\n'))
    fullFit = lm(cumulative_grains_below_aperture ~ time, data)
    print(fullFit)
    windowFit = lm(cumulative_grains_below_aperture ~ time, windowData)
    print(windowFit)
#     abline(windowFit)  # adds the line to the plot
    
    # Set up linear fit line for plotting
    x = c(max(minTime, window[[1]]), min(maxTime, window[[2]]))
    m = windowFit$coefficients[[2]]
    slopeLabels = c(slopeLabels, sprintf('%.2f', m))
    b = windowFit$coefficients[[1]]
    yWindow = m*x + b
    lines(x, yWindow, col=getColor(angle, 0.3), lwd=2)  # lwd=line width
    # curve(m*x + b, add=T, from=window[[1]], to=window[[2]])  # same as above
    # text(10, 10, labels=sprintf('Slope = %.2f', m))
    # plot(fullFit)  # interesting; it shows a number of plots
}

legend('bottomright', legend=names(paths), title='Angle', 
       fill=legendColors, cex=0.8,
       bty='n', inset=0) # bty='n' is no box type (border)

legend('bottom', legend=slopeLabels, title='Slope', 
       fill=legendColors, cex=0.8,
       bty='n', inset=0) # bty='n' is no box type (border)

