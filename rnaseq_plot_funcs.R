
plot_MA = function(logCounts, logFoldChange, FDR, FDR_thresh=0.05, logFC_thresh=2, xlab="logCounts", ylab="logFC", title="MA plot", pch='.') {

    plot(logCounts, logFoldChange, col=ifelse((FDR<FDR_thresh) & (abs(logFoldChange)>=logFC_thresh), "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_Volcano = function(logFoldChange, FDR, FDR_thresh=0.05, logFC_thresh=2,  xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch='.') {

   plot(logFoldChange, -1*log10(FDR), col=ifelse((FDR<FDR_thresh) & (abs(logFoldChange)>=logFC_thresh), "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}



#plot_MA = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch='.') {
#
#    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;
#
#}
#
#
#plot_Volcano = function(logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch='.') {
#
#   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
#
#}


plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    gridlayout = matrix(c(1:2),nrow=1,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(2,2), heights=c(1,1), respect=T) 

    # draw again, but use a smaller dot for data points
    plot_MA(logCounts, logFoldChange, FDR, pch='.');
    plot_Volcano(logFoldChange, FDR, pch='.');
    

    par(def.par)   
        
    
}


#plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {
#
#    def.par = par(no.readonly = TRUE) # save default, for resetting...
#
#    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
#    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
#
#    plot_MA(logCounts, logFoldChange, FDR);
#    plot_Volcano(logFoldChange, FDR);
#
#    # draw again, but use a smaller dot for data points
#    plot_MA(logCounts, logFoldChange, FDR, pch='.');
#    plot_Volcano(logFoldChange, FDR, pch='.');
#    
#
#    par(def.par)   
#        
#    
#}
