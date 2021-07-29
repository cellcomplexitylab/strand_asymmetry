RCircos.Draw.Chromosome.Ideogram <- function (ideo.pos=NULL, ideo.width=NULL)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    RCircos.Par  <- RCircos.Get.Plot.Parameters()
    
    if(is.null(ideo.pos)) ideo.pos <- RCircos.Par$chr.ideo.pos;
    if(is.null(ideo.width)) ideo.width <- RCircos.Par$chrom.width;
    
    #   Plot outlines for each chromosome
    #   =================================
    #   
    outerPos <- ideo.pos + ideo.width;
    innerPos <- ideo.pos;

    chromosomes <- unique(RCircos.Cyto$Chromosome);
    RCircos.Track.Outline(outerPos, innerPos, num.layers=1, chromosomes,
            track.colors=rep("gray80", length(chromosomes)));
#            track.colors=rep("white", length(chromosomes)));

# Do not add Giemsa bands after all...
#    #   Add chromosome bands (Giemsa stain positive only)
#    #   ================================================
#    #   
#    whiteBands <- which(RCircos.Cyto$BandColor == "white");
#   darkBands <- RCircos.Cyto; 
#   if(length(whiteBands)>0) darkBands <- darkBands[-whiteBands, ];
#
#    for(aBand in seq_len(nrow(darkBands)))
#    {   
#        aColor <- darkBands$BandColor[aBand];
#        aStart <- darkBands$StartPoint[aBand];
#        aEnd   <- darkBands$EndPoint[aBand];
#
#        posX <- c(RCircos.Pos[aStart:aEnd,1]*outerPos, 
#                RCircos.Pos[aEnd:aStart,1]*innerPos);
#        posY <- c(RCircos.Pos[aStart:aEnd,2]*outerPos, 
#                RCircos.Pos[aEnd:aStart,2]*innerPos);
#        polygon(posX, posY, col=aColor, border=NA);
#    }   
}


RCircos.Label.Chromosome.Names <- function (chr.name.pos=NULL)
{
    RCircos.Par  <- RCircos.Get.Plot.Parameters()
    if(is.null(chr.name.pos)) {
        chr.name.pos <- RCircos.Par$chr.name.pos;
    } else {
        if(chr.name.pos < RCircos.Par$track.in.start)
            stop("Chromosome name positions overlap with inside track.\n");
        if(chr.name.pos > RCircos.Par$track.out.start)
            message("May not plot data tracks outside of chromosomes.\n")
    }
   
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    chroms <- unique(RCircos.Cyto$Chromosome);

    rightSide <- nrow(RCircos.Pos)/2;
    for(aChr in seq_len(length(chroms)))
    {
        chrRows <- which(RCircos.Cyto$Chromosome == chroms[aChr]);
        chrStart <- min(RCircos.Cyto$StartPoint[chrRows]);
        chrEnd   <- max(RCircos.Cyto$EndPoint[chrRows]);
        mid <- round((chrEnd - chrStart + 1) / 2, digits=0) + chrStart;

        chrName <- sub(pattern="chr", replacement="", chroms[aChr]);
        text(RCircos.Pos[mid, 1]*RCircos.Par$chr.name.pos,
             RCircos.Pos[mid, 2]*RCircos.Par$chr.name.pos,
             label=chrName);
#             label=chrName, pos=ifelse(mid <= rightSide, 4, 2),
#             srt=RCircos.Pos$degree[mid]);
    }
}


Tile.Plot = function(tile.data, track.num, side, col="black", lwd=1) {
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    tile.data <- RCircos.Get.Single.Point.Positions(genomic.data=tile.data)
    # Pre-process data and arrange tiles in layers.
    the.layer <- 1
    the.chr <- tile.data[1, 1]
    start <- tile.data[1, 2]
    end <- tile.data[1, 3]
    tile.layers <- rep(1, nrow(tile.data))
    # The piece of code commented out below was used to put tiles on
    # different layers in case they overlap. In this version, the
    # tiles are put on the same layer.
#    for (a.row in 2:nrow(tile.data)) {
#        if (tile.data[a.row, 2] >= end) {
#            the.layer <- 1
#            start <- tile.data[a.row, 2]
#            end <- tile.data[a.row, 3]
#        }
#        else if (tile.data[a.row, 1] != the.chr) {
#            the.layer <- 1
#            the.chr <- tile.data[a.row, 1]
#            start <- tile.data[a.row, 2]
#            end <- tile.data[a.row, 3]
#        }
#        else {
#            the.layer <- the.layer + 1
#            if (tile.data[a.row, 3] > end) {
#                end <- tile.data[a.row, 3]
#            }
#        }
#        tile.layers[a.row] <- the.layer
#    }
    locations <- RCircos.Get.Track.Positions(side, track.num)
    out.pos <- locations[1]
    in.pos <- locations[2]
    layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers
    num.layers <- max(tile.layers)
    if (num.layers > RCircos.Par$max.layers) {
        if (side == "in") {
            in.pos <- out.pos - layer.height * num.layers
        }
        else {
            out.pos <- in.pos + layer.height * num.layers
        }
        cat(paste("Tiles plot will use more than one track.", 
            "Please select correct area for next track.\n"))
    }
    if (num.layers < RCircos.Par$max.layers) {
        layer.height <- RCircos.Par$track.height/num.layers
    }
    # tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color)
    # Do not plot outline, just the tiles.
    # RCircos.Track.Outline(out.pos, in.pos, num.layers)
    the.loc <- ncol(tile.data)
    tile.data$col = col
    for (a.row in 1:nrow(tile.data)) {
        tile.len <- tile.data[a.row, 3] - tile.data[a.row, 2]
        tile.range <- round(tile.len/RCircos.Par$base.per.unit/2, 
            digits = 0)
        start <- tile.data[a.row, the.loc] - tile.range
        end <- tile.data[a.row, the.loc] + tile.range
        layer.bot <- in.pos + layer.height * (tile.layers[a.row] - 
            1)
        layer.top <- layer.bot + layer.height * 0.8
        polygon.x <- c(RCircos.Pos[start:end, 1] * layer.top, 
            RCircos.Pos[end:start, 1] * layer.bot)
        polygon.y <- c(RCircos.Pos[start:end, 2] * layer.top, 
            RCircos.Pos[end:start, 2] * layer.bot)
        # polygon(polygon.x, polygon.y, col = tile.colors[a.row])
        polygon(polygon.x, polygon.y, border=tile.data$col[a.row], lwd=lwd)
    }
}

