data = read.csv('anchor_accels.csv', header=T)
x = data
maxTime = max(x$time)
# x$ix = factor(x$ix)  # stringifies, FYI
# x$time = factor(x$time)
# Oddly large values, and for the middle anchors generally
x$time[x$value > 100]
x$ix[x$value > 100]

# numToTake = length(unique(x$ix)) * 10
# # steady = subset(x, subset= time > 15)  # arbitrary
# steady = tail(x, n=numToTake)
# 
# x = steady
# gs = split(x, x$ix)  # creates a list of dfs with keys=ixs (as strings e.g. "0")
# summary(gs[[1]]$value)
# sd(gs[[1]]$value)

b = 2.5; t = 1
l = b; r = t
par(mar=c(b, l, t, r))

for (x in list(data, subset(data, subset= time > maxTime-1))) {
    gByIX = split(x, x$ix)  # over >1 time value
    # collapse times into one statistic
    # v_min1qMedianMean3qMax = summary(gByIX[[1]]$value)
    statsPerIx = mapply(gByIX, FUN = function(df) summary(df$value))
    boxplot(statsPerIx)  # <-- What I want
    # rows are stat headings, so transpose to get those as col names and ix's as rows
    df = as.data.frame(t(statsPerIx))
}


# barplot(x$value[x$time == x$time[[10]]])
# Maybe color those left of the middle differently
