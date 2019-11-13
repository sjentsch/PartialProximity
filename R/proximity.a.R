# =====================================================================================================================
# Algorithms to calculate proximities
# =====================================================================================================================
algPxm <- function(dtaMtx, mthNme, mthExp) {
    # =================================================================================================================
    # integer measures that implemented using R functions (cor, dist)
    # =================================================================================================================
    # integer - similarity - Pearson correlation
    if      (mthNme == "intCrr") {
        resMtx = cor(dtaMtx)
    }
    # integer - dissimilarity - Euclidian distance, squared Euclidian distance, Chebychev, Block, Minkowski, Customized
    else if (mthNme == "intEuc" || mthNme == "intSqE" || mthNme == "intChb" || mthNme == "intBlk" || mthNme == "intMnk" || mthNme == "intCst") {
       # distance measures from the R-function dist
       dstMth = ifelse(mthNme == "intEuc" || mthNme == "intSqE", "euclidian", ifelse(mthNme == "intChb", "maximum", ifelse(mthNme == "intBlk", "manhattan", ifelse(mthNme == "intMnk" || mthNme == "intCst", "minkowski", "error"))))
       dstExp = ifelse(mthNme == "intSqE", 2, ifelse(mthNme == "intCst", (mthExp[1] / mthExp[2]), 1))
       resMtx = as.matrix(dist(t(dtaMtx), upper=T, diag=T, method=dstMth, p = mthExp[1]) ^ dstExp)
    }
    # =================================================================================================================
    # integer, count and binary measures; implemented based upon:
    # www.ibm.com/support/knowledgecenter/SSLVMB_22.0.0/com.ibm.spss.statistics.algorithms/alg_proximities.htm
    # see www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/binmatch.htm for further possible measures
    # =================================================================================================================
    else {
        # create result matrix
        resMtx = array(0, c(ncol(dtaMtx), ncol(dtaMtx)))
        dimnames(resMtx) = list(colnames(dtaMtx), colnames(dtaMtx))

        # preparation: binary measures
        if (mthNme == "binShp")
            dtaMtx = scale(dtaMtx, center=T, scale=F)

        # for most calculations, the main diagonal doesn't need to be calculated
        # however, there are three measures that represent an exception from that
        # those methods are the same where the main diagonal is divided by two at
        #  the very end of this algorithm
        if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binRnR")
            minI = 1
        else
            minI = 2

        for (i in minI:ncol(dtaMtx)) {
            for (j in 1:(i - minI + 1)) {
                # =====================================================================================================
                # calculate match vector and t1 / t2 required for binary measures
                # =====================================================================================================
                if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binKc2" || mthNme == "binLmb" || mthNme == "binOch" || mthNme == "binPh4" || mthNme == "binPtD" || 
                    mthNme == "binSk4" || mthNme == "binSk5" || mthNme == "binSzD" || mthNme == "binVar" || mthNme == "binYlQ" || mthNme == "binYlY") {
                    mtcVec = c(as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] == 1)["TRUE"]), as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] != 1)["TRUE"]), as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] == 1)["TRUE"]), as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] != 1)["TRUE"]))
                    if (mthNme == "binAnD" || mthNme == "binLmb") {
                        t1 = max(mtcVec[1], mtcVec[2], na.rm=T) + max(mtcVec[3], mtcVec[4], na.rm=T) + max(mtcVec[1], mtcVec[3], na.rm=T) + max(mtcVec[2], mtcVec[4], na.rm=T)
                        t2 = max(sum(mtcVec[1], mtcVec[3], na.rm=T), sum(mtcVec[2], mtcVec[4], na.rm=T), na.rm=T) + max(sum(mtcVec[1], mtcVec[2], na.rm=T), sum(mtcVec[3], mtcVec[4], na.rm=T), na.rm=T)
                    }
                }

                # =====================================================================================================
                # measures for integer data
                # =====================================================================================================
                # integer - similarity - Cosine
                if      (mthNme == "intCos") {
                    resMtx[i, j] = crossprod(dtaMtx[, i], dtaMtx[, j]) / sqrt(crossprod(dtaMtx[, i]) * crossprod(dtaMtx[, j]))
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "int") {
#               }    
                # =====================================================================================================
                # measures for count data
                # =====================================================================================================
                # dissimilarities - counts
                    # DISSIMILARITIES
                    # Chi-square measure: disChi
                    # Phi-square measure: disPhi
                else if (mthNme == "cntChi") {
                    resMtx[i, j] = NA # to implement
                }
                else if (mthNme == "cntPhi") {
                    resMtx[i, j] = NA # to implement
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "cnt") {
#               }
                # =====================================================================================================
                # measures for binary data
                # =====================================================================================================
                # binary - similarity - Russel and Rao: binRnR
                    # Simple matching: binSmM
                    # Jaccard: binJcc
                    # Dice: binDic
                    # Rogers and Tanimoto: binRnT
                    # Sokal and Sneath 1 - 5: binSk1, binSk2, binSk3, binSk4, binSk5)
                    # Kulczynski 1 - 2: binKc1, binKc2
                    # Hamann: binHmn
                    # Lambda: binLmb
                    # Anderberg's D: binAnD
                    # Yule's Y & Q: binYlY, binYlQ
                    # Ochiai: binOch
                    # Phi 4-point correlation: binPh4
                    # Dispersion: binDsp

                # binary - dissimilarity
                    # DISSIMILARITIES
                    # Euclidian distance: binEuc - BEUCLID
                    # Squared Euclidian distance: binSqE - BSEUCLID
                    # Size difference: binSzD - SIZE
                    # Pattern difference: binPtD - PATTERN
                    # Variance: binVar - VARIANCE
                    # Shape: binShp - BSHAPE
                    # Lance and Williams: binLnW - BLWMN

                else if (mthNme == "binAnD") {
                    resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx))
                }
                else if (mthNme == "binDic") {
                    resMtx[i, j] = (as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) * 2) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                else if (mthNme == "binDsp") {
                    resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - ifelse(is.na(mtcVec[2]) || is.na(mtcVec[3]), 0, (mtcVec[2] * mtcVec[3]))) / (nrow(dtaMtx) ^ 2)
                }
                else if (mthNme == "binEuc") {
                    resMtx[i, j] = sqrt(as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"]))
                }
                else if (mthNme == "binJcc") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"])
                }
                else if (mthNme == "binHmm") {
                    resMtx[i, j] = (as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"])) / nrow(dtaMtx)
                }
                else if (mthNme == "binLmb") {
                    resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx) - t2)
                }
                else if (mthNme == "binKc1") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                else if (mthNme == "binKc2") {
                    resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3])) / 2
                }
                else if (mthNme == "binLnW") {
                    resMtx[i, j] = sum(abs(dtaMtx[, i] - dtaMtx[, j])) / sum(dtaMtx[, i] + dtaMtx[, j])
                }
                else if (mthNme == "binOch") {
                    resMtx[i, j] = sqrt((mtcVec[1] / (mtcVec[1] + mtcVec[2])) * (mtcVec[1] / (mtcVec[1] + mtcVec[3])))
                }
                else if (mthNme == "binPh4") {
                    resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - (mtcVec[2] * mtcVec[3])) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
                }
                else if (mthNme == "binPtD") {
                    resMtx[i, j] = (mtcVec[2] * mtcVec[3]) / (nrow(dtaMtx) ^ 2)
                }
                else if (mthNme == "binRnR") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
                }
                else if (mthNme == "binRnT") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + 2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
                }
                else if (mthNme == "binShp") {
                    resMtx[i, j] = sum((dtaMtx[, i] - dtaMtx[, j]) ^ 2) / nrow(dtaMtx)
                }
                else if (mthNme == "binSk1") {
                    resMtx[i, j] = (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"])) / (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
                }
                else if (mthNme == "binSk2") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (2 * as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                else if (mthNme == "binSk3") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"])
                }
                else if (mthNme == "binSk4") {
                    resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3]) + mtcVec[4] / (mtcVec[2] + mtcVec[4]) + mtcVec[4] / (mtcVec[3] + mtcVec[4])) / 4
                }
                else if (mthNme == "binSk5") {
                    resMtx[i, j] = (mtcVec[1] * mtcVec[4]) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
                }
                else if (mthNme == "binSmM") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
                }
                else if (mthNme == "binSqE") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"])
                }
                else if (mthNme == "binSzD") {
                    resMtx[i, j] = (mtcVec[2] - mtcVec[3]) ^ 2 / (nrow(dtaMtx) ^ 2)
                }
                else if (mthNme == "binVar") {
                    resMtx[i, j] = (mtcVec[2] + mtcVec[3]) / (4 * nrow(dtaMtx))
                }
                else if (mthNme == "binYlQ") {
                    resMtx[i, j] = (mtcVec[1] * mtcVec[4] - mtcVec[2] * mtcVec[3]) / (mtcVec[1] * mtcVec[4] + mtcVec[2] * mtcVec[3])
                }
                else if (mthNme == "binYlY") {
                    resMtx[i, j] = (sqrt(mtcVec[1] * mtcVec[4]) - sqrt(mtcVec[2] * mtcVec[3])) / (sqrt(mtcVec[1] * mtcVec[4]) + sqrt(mtcVec[2] * mtcVec[3]))
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "bin") {
#               }
            }
        }

        # =============================================================================================================
        # transpose and add the matrix (only the bottom triangle is calculated, th upper is "mirrored")
        # =============================================================================================================
        resMtx = resMtx + t(resMtx)
        # handle the main diagonal
        if      (mthNme == "binEuc" || mthNme == "binLnW" || mthNme == "binPtD" || mthNme == "binShp" || mthNme == "binSqE" || mthNme == "binSzD" || mthNme == "binVar")
            diag(resMtx) = 0
        else if (mthNme == "intCos" || mthNme == "cntChi" || mthNme == "cntChi" || mthNme == "binDic" || mthNme == "binJcc" || mthNme == "binHmm" || mthNme == "binKc2" || 
                 mthNme == "binLmb" || mthNme == "binOch" || mthNme == "binPh4" || mthNme == "binRnT" || mthNme == "binSmM" || mthNme == "binSk1" || mthNme == "binSk2" ||
                 mthNme == "binSk4" || mthNme == "binSk5" || mthNme == "binYlY" || mthNme == "binYlQ")
            diag(resMtx) = 1
        else if (mthNme == "binKc1" || mthNme == "binSk3")
            diag(resMtx) = NA
        else if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binRnR")
            diag(resMtx) = diag(resMtx) / 2
    
    }
    # =================================================================================================================
    # end: own implementations
    # =================================================================================================================

    # =================================================================================================================
    # return matrix with results
    # =================================================================================================================
    return(resMtx)
}

xfmDta <- function(dtaMtx, xfmMth, xfmDir) {
    xfmDim = ifelse(xfmDir == "xfmVar", 2, ifelse(xfmDir == "xfmSbj", 1, NA))

    # z-scores - Z - correct
    if      (xfmMth == "xfmZsc")
        xfmMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)),  "-"), xfmDim, as.vector(apply(dtaMtx, xfmDim, sd)),          "/")
    # range -1 to 1 - RANGE - over subjects correct, over variables possibly wrong in SPSS
    else if (xfmMth == "xfmRNP")
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
    # range 0 to 1 - RESCALE - over subjects correct, over variables possibly wrong in SPSS
    else if (xfmMth == "xfmRZP")
        xfmMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, min)),   "-"), xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
    # maximum magnitude of 1 - MAX - correct
    else if (xfmMth == "xfmMag")
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, max)),  "/")
    # mean of 1 - MEAN - correct
    else if (xfmMth == "xfmAvr")
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)), "/")
    # standard deviation of 1 - SD - correct
    else if (xfmMth == "xfmStd")
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, sd)),   "/")

    return(xfmMtx)
}


xfmRes <- function(resMtx, blnXfR) {
    if (blnXfmR[1])
        resMtx = abs(resMtx)
    if (blnXfmR[2])
        resMtx[row(resMtx) != col(resMtx)] = -resMtx[row(resMtx) != col(resMtx)]
    if (blnXfmR[3])
        resMtx[row(resMtx) != col(resMtx)] = (resMtx[row(resMtx) != col(resMtx)] - min(resMtx[row(resMtx) != col(resMtx)])) / diff(range(resMtx[row(resMtx) != col(resMtx)]))

    return(resMtx)
}

nmePxm <- function(mthNme) {
    # measures for interval data
    if      (mthNme == "intBlk")
        resNme = "Block"
    else if (mthNme == "intChb")
        resNme = "Chebychev"
    else if (mthNme == "intCst")
        resNme = "Customized (based upon Minkowski; Power: [IP], Root: [IR])"
    else if (mthNme == "intCos")
        resNme = "Cosine"
    else if (mthNme == "intCrr")
        resNme = "Pearson correlation"
    else if (mthNme == "intEuc")
        resNme = "Euclidian distance"
    else if (mthNme == "intMnk")
        resNme = "Minkowski (Power: [IP])"
    else if (mthNme == "intSqE")
        resNme = "Squared Euclidian distance"
    # count / frequency measures
    else if (mthNme == "cntChi")
        resNme = "Chi-squared measure"
    else if (mthNme == "cntPhi")
        resNme = "Phi-squared measure"
    # binary measures
    else if (mthNme == "binAnD")
        resNme = "Anderberg's D (present: [BP], absent: [BA])"
    else if (mthNme == "binDic")
        resNme = "Dice (present: [BP], absent: [BA])"
    else if (mthNme == "binDsp")
        resNme = "Dispersion (present: [BP], absent: [BA])"
    else if (mthNme == "binEuc")
        resNme = "(Binary) Euclidian distance (present: [BP], absent: [BA])"
    else if (mthNme == "binHmn")
        resNme = "Hamann (present: [BP], absent: [BA])"
    else if (mthNme == "binJcc")
        resNme = "Jaccard (present: [BP], absent: [BA])"
    else if (mthNme == "binKc1")
        resNme = "Kulczynski 1 (present: [BP], absent: [BA])"
    else if (mthNme == "binKc2")
        resNme = "Kulczynski 2 (present: [BP], absent: [BA])"
    else if (mthNme == "binLmb")
        resNme = "Lambda (present: [BP], absent: [BA])"
    else if (mthNme == "binLnW")
        resNme = "Lance and Williams (present: [BP], absent: [BA])"
    else if (mthNme == "binOch")
        resNme = "Ochiai (present: [BP], absent: [BA])"
    else if (mthNme == "binPh4")
        resNme = "Phi 4-point correlation (present: [BP], absent: [BA])"
    else if (mthNme == "binPtD")
        resNme = "Pattern difference (present: [BP], absent: [BA])"
    else if (mthNme == "binRnR")
        resNme = "Russel and Rao (present: [BP], absent: [BA])"
    else if (mthNme == "binRnT")
        resNme = "Rogers and Tanimoto (present: [BP], absent: [BA])"
    else if (mthNme == "binShp")
        resNme = "Shape (present: [BP], absent: [BA])"
    else if (mthNme == "binSk1")
        resNme = "Sokal and Sneath 1 (present: [BP], absent: [BA])"
    else if (mthNme == "binSk2")
        resNme = "Sokal and Sneath 2 (present: [BP], absent: [BA])"
    else if (mthNme == "binSk3")
        resNme = "Sokal and Sneath 3 (present: [BP], absent: [BA])"
    else if (mthNme == "binSk4")
        resNme = "Sokal and Sneath 4 (present: [BP], absent: [BA])"
    else if (mthNme == "binSk5")
        resNme = "Sokal and Sneath 5 (present: [BP], absent: [BA])"
    else if (mthNme == "binSmM")
        resNme = "Simple matching (present: [BP], absent: [BA])"
    else if (mthNme == "binSqE")
        resNme = "(Binary) squared Euclidian distance (present: [BP], absent: [BA])"
    else if (mthNme == "binSzD")
        resNme = "Size difference (present: [BP], absent: [BA])"
    else if (mthNme == "binVar")
        resNme = "(Binary) variance (present: [BP], absent: [BA])"
    else if (mthNme == "binYlQ")
        resNme = "Yule's Q (present: [BP], absent: [BA])"
    else if (mthNme == "binYlY")
        resNme = "Yule's Y (present: [BP], absent: [BA])"
#   else if (mthNme == "")
#       resNme = ""
    else
        resNme = ""

    return(resNme)
}

spsPxm <- function(mthNme) {
    # measures for interval data
    if      (mthNme == "intBlk")
        spssMt = "BLOCK"
    else if (mthNme == "intChb")
        spssMt = "CHEBYCHEV"
    else if (mthNme == "intCst")
        spssMt = "POWER([IP], [IR])"
    else if (mthNme == "intCos")
        spssMt = "COSINE"
    else if (mthNme == "intCrr")
        spssMt = "CORRELATION"
    else if (mthNme == "intEuc")
        spssMt = "EUCLID"
    else if (mthNme == "intMnk")
        spssMt = "MINKOWSKI([IP])"
    else if (mthNme == "intSqE")
        spssMt = "SEUCLID"
    # count / frequency measures
    else if (mthNme == "cntChi")
        spssMt = "CHISQ"
    else if (mthNme == "cntPhi")
        spssMt = "PH2"
    # binary measures
    else if (mthNme == "binAnD")
        spssMt = "D([BP], [BA])"
    else if (mthNme == "binDic")
        spssMt = "DICE([BP], [BA])"
    else if (mthNme == "binDsp")
        spssMt = "DISPER([BP], [BA])"
    else if (mthNme == "binEuc")
        spssMt = "BEUCLID([BP], [BA])"
    else if (mthNme == "binHmn")
        spssMt = "HAMANN([BP], [BA])"
    else if (mthNme == "binJcc")
        spssMt = "JACCARD([BP], [BA])"
    else if (mthNme == "binKc1")
        spssMt = "K1([BP], [BA])"
    else if (mthNme == "binKc2")
        spssMt = "K2([BP], [BA])"
    else if (mthNme == "binLmb")
        spssMt = "LAMBDA([BP], [BA])"
    else if (mthNme == "binLnW")
        spssMt = "BLWMN([BP], [BA])"
    else if (mthNme == "binOch")
        spssMt = "OCHIAI([BP], [BA])"
    else if (mthNme == "binPh4")
        spssMt = "PHI([BP], [BA])"
    else if (mthNme == "binPtD")
        spssMt = "PATTERN([BP], [BA])"
    else if (mthNme == "binRnR")
        spssMt = "RR([BP], [BA])"
    else if (mthNme == "binRnT")
        spssMt = "RT([BP], [BA])"
    else if (mthNme == "binShp")
        spssMt = "BSHAPE([BP], [BA])"
    else if (mthNme == "binSk1")
        spssMt = "SS1([BP], [BA])"
    else if (mthNme == "binSk2")
        spssMt = "SS2([BP], [BA])"
    else if (mthNme == "binSk3")
        spssMt = "SS3([BP], [BA])"
    else if (mthNme == "binSk4")
        spssMt = "SS4([BP], [BA])"
    else if (mthNme == "binSk5")
        spssMt = "SS5([BP], [BA])"
    else if (mthNme == "binSmM")
        spssMt = "SM([BP], [BA])"
    else if (mthNme == "binSqE")
        spssMt = "BSEUCLID([BP], [BA])"
    else if (mthNme == "binSzD")
        spssMt = "SIZE([BP], [BA])"
    else if (mthNme == "binVar")
        spssMt = "VARIANCE([BP], [BA])"
    else if (mthNme == "binYlQ")
        spssMt = "Q([BP], [BA])"
    else if (mthNme == "binYlY")
        spssMt = "Y([BP], [BA])"
#   else if (mthNme == "")
#       spssMt = ""
    else
        spssMt = ""

    return(spssMt)
}
