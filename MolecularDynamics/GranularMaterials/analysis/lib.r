## savePath is if you want to save as an image file. Requires opening a window.
getSaveFunc = function(savePath, windowWidth=10, windowHeight=10) {
    f = function() {}
    if (is.null(savePath)) return(f)  # do nothing
    require('stringr')  # for file extension parsing
    windows(width=windowWidth, height=windowHeight)  # need to open new window in order to save
#     par(mgp=c(1.5, 0.5, 0))  # smaller margins axis title, labels, and line
    endFunc = function() dev.off()  # close window when done
    if (str_detect(savePath, ignore.case('svg$'))) {
        f = function() {
            dev.copy(svg, file=savePath)
            endFunc()
        }
    } else {  # default is png
        f = function () {
            savePlot(filename=savePath, type='png')
            endFunc()
        }
    }
    return(f)
    # fname=paste('images/', name, ' ', operationName, sep='')  # need to create 'images/' subfolder manually or R will error
}
