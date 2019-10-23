
# This file is a generated template, your changes will not be overwritten

proximityClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "proximityClass",
    inherit = proximityBase,
    private = list(
        # ====================================================
        .init = function() {
            # get variables
            pxmVar = self$options$get('vars')
            pxmLbl = self$options$get('label')
            if (length(pxmVar) > 1) {
print('init')
dtaRaw = self$options$.getData()
print(str(dtaRaw))
saveRDS(dtaRaw, file='/home/sjentsch/Downloads/Trial1.rds')
                dtaMtx = self$data[, pxmVar]
saveRDS(dtaMtx, file='/home/sjentsch/Downloads/Trial2.rds')
print(self$data)
print(str(dtaMtx))
                sbjLbl = self$data[, pxmLbl]
                blnVld = !apply(is.na(dtaMtx), 1, any)
print(blnVld)
                dtaMtx = dtaMtx[blnVld, ]
                sbjLbl = sbjLbl[blnVld, ]
                numIvN = length(blnVld[!blnVld])
                # for binary data, ensure that these only contain the two categories defined by lvlMsr == lvlBin
                if (self$options$get('lvlMsr') == 'lvlBin') {
print(c(self$options$get('binAbs'), self$options$get('binPrs')))
print(str(dtaMtx))
print(as.matrix(dtaMtx) == as.numeric(self$options$get('binPrs')))
                    blnBnC = apply(dtaMtx == as.numeric(self$options$get('binPrs')) | dtaMtx == as.numeric(self$options$get('binAbs')), 1, all)
                    dtaMtx = dtaMtx[blnBnC, ]
                    sbjLbl = sbjLbl[blnBnC, ]
                    numIvB = length(blnBnC[!blnBnC])
                }
                numSbj = dim(dtaMtx)[1]
                numVar = dim(dtaMtx)[2]

                # initialize tables

# implement case processing summary as notes
            }
        },
        # ====================================================
        .run = function() {
            # get variables
            pxmVar = self$options$get('vars')
            if (length(pxmVar) > 1) {
print('run')
print(self$data)
                dtaMtx = self$data[, pxmVar]
print(dtaMtx)
                blnVld = !apply(is.na(dtaMtx), 1, any)
                dtaMtx = dtaMtx[blnVld, ]
                # for binary data, ensure that these only contain the two categories defined by lvlMsr == lvlBin
                if (self$options$get('lvlMsr') == 'lvlBin') {
                    blnBnC = apply(dtaMtx == as.numeric(self$options$get('binPrs')) | dtaMtx == as.numeric(self$options$get('binAbs')), 1, all)
                    dtaMtx = dtaMtx[blnBnC, ]
# replace with present and absent
                }
                numSbj = dim(dtaMtx)[1]
                numVar = dim(dtaMtx)[2]
                clcDis = (self$options$get('disSim') == 'clcDis')
                clcSim = (self$options$get('disSim') == 'clcSim')
print(dim(dtaMtx))
print(clcSim)
print(clcDis)

                # transform data (if necessary)
                # The Transform Values group allows you to standardize data values for either cases or variables before computing proximities.
                # These transformations are not applicable to binary data. Available standardization methods are z scores, range –1 to 1,
                # range 0 to 1, maximum magnitude of 1, mean of 1, and standard deviation of 1.
                if (self$options$get('xfmMth') != 'xfmNon')
                    dtaMtx = xfmDta(dtaMtx, self$options$get('xfmMth'), self$options$get('xfmDir'))
                # transpose the matrix if similarities between subject are to be calculated
                if (self$options$get('btwDir') == 'btwSbj')
                    dtaMtx = t(dtaMtx)

                # calculate proximity measures
                # similarities and dissimilarites - interval
                if      (self$options$get('lvlMsr') == 'lvlInt') {
print("lvlInt")
                    # SIMILARITIES
                    # Pearson correlation: simCrr - CORRELATION - correct
                    # Cosine: simCos - COSINE - correct
                    # DISSIMILARITIES
                    # Euclidian distance:  disEuc - EUCLID
                    # Squared Euclidian distance: disSqE - SEUCLID
                    # Chebychev: disChb - CHEBYCHEV
                    # Block: disBlk - BLOCK
                    # Minkowski: disMnk - MINKOWSKI(p)
                    # Customized: disCst - POWER(p [disInP], r [disInR])
                    resMtx = algInt(dtaMtx,
                                    ifelse(clcDis, substr(self$options$get('intDis'), 4, 6), ifelse(clcSim, substr(self$options$get('intSim'), 4, 6), '')), # mthInt
                                    ifelse(clcDis, c(self$options$get('intPwr'), self$options$get('intRot')), c(0, 0)))                                     # expInt
                }
                # similarities and dissimilarities - binary
                # see https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/binmatch.htm for further possible measures
                else if (self$options$get('lvlMsr') == 'lvlBin') {
print("lvlBin")
                    # SIMILARITES
                    # Russel and Rao: simRnR
                    # Simple matching: simSmM
                    # Jaccard: simJcc
                    # Dice: simDic
                    # Rogers and Tanimoto: simRnT
                    # Sokal and Sneath 1 - 5: simSk1, simSk2, simSk3, simSk4, simSk5)
                    # Kulczynski 1 - 2: simKc1, simKc2
                    # Hamann: simHmn
                    # Lambda: simLmb
                    # Anderberg's D: simAnD
                    # Yule's Y & Q: simYlY, simYlQ
                    # Ochiai: simOch
                    # Phi 4-point correlation: simPh4
                    # Dispersion: simDsp
                    # DISSIMILARITIES
                    # Euclidian distance: disEuc - BEUCLID
                    # Squared Euclidian distance: disSqE - BSEUCLID
                    # Size difference: disSzD - SIZE
                    # Pattern difference: disPtD - PATTERN
                    # Variance: disVar - VARIANCE
                    # Shape: disShp - BSHAPE
                    # Lance and Williams: disLnW - BLWMN
                   resMtx = algBin(dtaMtx,
                                   ifelse(clcDis, substr(self$options$get('binDis'), 4, 6), ifelse(clcSim, substr(self$options$get('binSim'), 4, 6), ''))) # mthBin
                }
                # dissimilarities - counts
                else if (self$options$get('lvlMsr') == 'lvlCnt') {
                    # DISSIMILARITIES
                    # Chi-square measure: disChi
                    # Phi-square measure: disPhi
print("lvlCnt")
                    resMtx = algBin(dtaMtx,
                                    ifelse(clcDis, substr(self$options$get('cntDis'), 4, 6), ''))                                                           # mthCnt
                }

                # The Transform Measures group allows you to transform the values generated by the distance measure. They are applied after
                # the distance measure has been computed. Available options are absolute values, change sign, and rescale to 0–1 range.
                blnXfR = c(self$options$get('xfmAbs'), self$options$get('xfmInv'), self$options$get('xfmRsc'))
                if (any(blnXfR))
                    resMtx = xfmRes(resMtx, blnXfR)

                # assign results
                if      (self$options$get('btwDir') == 'btwVar') {
                    print(resMtx)
                }
                else if (self$options$get('btwDir') == 'btwSbj') {
    #                   print(resMtx)
                }
            }
        }
    )
)

algInt <- function(dtaMtx, mthNme, mthExp) {
    if      (mthNme == 'Crr') {
        resMtx = cor(dtaMtx)
    }
    else if (mthNme == 'Cos') {
        resMtx = array(0, c(ncol(dtaMtx), ncol(dtaMtx)))
        dimnames(resMtx) = list(colnames(dtaMtx), colnames(dtaMtx))
        for (i in 2:ncol(dtaMtx)) {
            for (j in 1:(i - 1)) {
                resMtx[i, j] = crossprod(dtaMtx[, i], dtaMtx[, j]) / sqrt(crossprod(dtaMtx[, i]) * crossprod(dtaMtx[, j]))
            }
        }
        resMtx = resMtx + t(resMtx)
        diag(resMtx) = 1
    }
    else {
       # distance measures from the R-function dist
       dstMth = ifelse(mthNme == 'Euc' || mthNme == 'SqE', 'euclidian', ifelse(mthNme == 'Chb', 'maximum', ifelse(mthNme == 'Blk', 'manhattan', ifelse(mthNme == 'Mnk' || mthNme == 'Cst', 'minkowski', 'error'))))
       dstExp = ifelse(mthNme == 'SqE', 2, ifelse(mthNme == 'Cst', (mthExp[1] / mthExp[2]), 1))
       resMtx = dist(t(dtaMtx), upper=T, diag=T, method=dstMth, p = mthExp[1]) ^ dstExp
    }

    return(resMtx)
}


algBin <- function(dtaMtx, mthNme) {
    resMtx = array(0, c(ncol(dtaMtx), ncol(dtaMtx)))
    dimnames(resMtx) = list(colnames(dtaMtx), colnames(dtaMtx))

    if (mthNme == "Shp")
        dtaMtx = scale(dtaMtx, center=T, scale=F)

    if (mthNme == "RnR" || mthNme == "AnD" || mthNme == "Dsp")
       minI = 1
    else
       minI = 2

    for (i in minI:ncol(dtaMtx)) {
        for (j in 1:(i - minI + 1)) {
            if (mthNme == "AnD" || mthNme == "Dsp" || mthNme == "Kc2" || mthNme == "Lmb" || mthNme == "Och" || mthNme == "Ph4" || mthNme == "PtD" || mthNme == "Sk4" || mthNme == "Sk5" || mthNme == "SzD" || mthNme == "Var" || mthNme == "YlQ" || mthNme == "YlY") {
                mtcVec = c(as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] == 1)['TRUE']), as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] != 1)['TRUE']), as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] == 1)['TRUE']), as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] != 1)['TRUE']))
                if (mthNme == "AnD" || mthNme == "Lmb") {
                    t1 = max(mtcVec[1], mtcVec[2], na.rm=T) + max(mtcVec[3], mtcVec[4], na.rm=T) + max(mtcVec[1], mtcVec[3], na.rm=T) + max(mtcVec[2], mtcVec[4], na.rm=T)
                    t2 = max(sum(mtcVec[1], mtcVec[3], na.rm=T), sum(mtcVec[2], mtcVec[4], na.rm=T), na.rm=T) + max(sum(mtcVec[1], mtcVec[2], na.rm=T), sum(mtcVec[3], mtcVec[4], na.rm=T), na.rm=T)
                }
            }

            else if (mthNme == "AnD") {
                resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx))
            }
            else if (mthNme == "Dic") {
                resMtx[i, j] = (as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) * 2) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
            }
            else if (mthNme == "Dsp") {
                resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - ifelse(is.na(mtcVec[2]) || is.na(mtcVec[3]), 0, (mtcVec[2] * mtcVec[3]))) / (nrow(dtaMtx) ^ 2)
            }
            else if (mthNme == "Euc") {
                resMtx[i, j] = sqrt(as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"]))
            }
            else if (mthNme == "Jcc") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"])
            }
            else if (mthNme == "Hmm") {
                resMtx[i, j] = (as.vector(table(dtaMtx[, i] == dtaMtx[, j])['TRUE']) - as.vector(table(dtaMtx[, i] != dtaMtx[, j])['TRUE'])) / nrow(dtaMtx)
            }
            else if (mthNme == "Lmb") {
                resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx) - t2)
            }
            else if (mthNme == "Kc1") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])['TRUE']) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])['TRUE']))
            }
            else if (mthNme == "Kc2") {
                resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3])) / 2
            }
            else if (mthNme == "LnW") {
                resMtx[i, j] = sum(abs(dtaMtx[, i] - dtaMtx[, j])) / sum(dtaMtx[, i] + dtaMtx[, j])
            }
            else if (mthNme == "Och") {
                resMtx[i, j] = sqrt((mtcVec[1] / (mtcVec[1] + mtcVec[2])) * (mtcVec[1] / (mtcVec[1] + mtcVec[3])))
            }
            else if (mthNme == "Ph4") {
                resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - (mtcVec[2] * mtcVec[3])) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
            }
            else if (mthNme == "PtD") {
                resMtx[i, j] = (mtcVec[2] * mtcVec[3]) / (nrow(dtaMtx) ^ 2)
            }
            else if (mthNme == "RnR") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
            }
            else if (mthNme == "RnT") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + 2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
            }
            else if (mthNme == "Shp") {
                resMtx[i, j] = sum((dtaMtx[, i] - dtaMtx[, j]) ^ 2) / nrow(dtaMtx)
            }
            else if (mthNme == "Sk1") {
                resMtx[i, j] = (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"])) / (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
            }
            else if (mthNme == "Sk2") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (2 * as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
            }
            else if (mthNme == "Sk3") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"])
            }
            else if (mthNme == "Sk4") {
                resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3]) + mtcVec[4] / (mtcVec[2] + mtcVec[4]) + mtcVec[4] / (mtcVec[3] + mtcVec[4])) / 4
            }
            else if (mthNme == "Sk5") {
                resMtx[i, j] = (mtcVec[1] * mtcVec[4]) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
            }
            else if (mthNme == "SmM") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
            }
            else if (mthNme == "SqE") {
                resMtx[i, j] = as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"])
            }
            else if (mthNme == "SzD") {
                resMtx[i, j] = (mtcVec[2] - mtcVec[3]) ^ 2 / (nrow(dtaMtx) ^ 2)
            }
            else if (mthNme == "Var") {
                resMtx[i, j] = (mtcVec[2] + mtcVec[3]) / (4 * nrow(dtaMtx))
            }
            else if (mthNme == "YlQ") {
                resMtx[i, j] = (mtcVec[1] * mtcVec[4] - mtcVec[2] * mtcVec[3]) / (mtcVec[1] * mtcVec[4] + mtcVec[2] * mtcVec[3])
            }
            else if (mthNme == "YlY") {
                resMtx[i, j] = (sqrt(mtcVec[1] * mtcVec[4]) - sqrt(mtcVec[2] * mtcVec[3])) / (sqrt(mtcVec[1] * mtcVec[4]) + sqrt(mtcVec[2] * mtcVec[3]))
            }
#           else if (mthNme == "") {
#           }
        }
    }

    resMtx = resMtx + t(resMtx)
    if      (mthNme == "Euc" || mthNme == "LnW" || mthNme == "PtD" || mthNme == "Shp" || mthNme == "SqE" || mthNme == "SzD" || mthNme == "Var")
        diag(resMtx) = 0
    else if (mthNme == "Dic" || mthNme == "Jcc" || mthNme == "Hmm" || mthNme == "Kc2" || mthNme == "Lmb" || mthNme == "Och" || mthNme == "Ph4" || mthNme == "RnT" ||
             mthNme == "SmM" || mthNme == "Sk1" || mthNme == "Sk2" || mthNme == "Sk4" || mthNme == "Sk5" || mthNme == "YlY" || mthNme == "YlQ")
        diag(resMtx) = 1
    else if (mthNme == "Kc1" || mthNme == "Sk3")
        diag(resMtx) = NA
    else if (mthNme == "AnD" || mthNme == "Dsp" || mthNme == "RnR")
        diag(resMtx) = diag(resMtx) / 2

    return(resMtx)
}


algCnt <- function(dtaMtx, mthNme) {
    resMtx = array(0, c(ncol(dtaMtx), ncol(dtaMtx)))
    dimnames(resMtx) = list(colnames(dtaMtx), colnames(dtaMtx))

    for (i in 2:ncol(dtaMtx)) {
        for (j in 1:(i - 1)) {
            if      (mthNme == 'Chi')
                resMtx = matrix(rep(0, nVar ^ 2), nrow = nVar) # to implement
            else if (mthNme == 'Phi')
                resMtx = matrix(rep(0, nVar ^ 2), nrow = nVar) # to implement
        }
    }

    resMtx = resMtx + t(resMtx)
    diag(resMtx) = 1

    return(resMtx)
}


xfmDta <- function(dtaMtx, xfmMth, xfmDir) {
    xfmDim = ifelse(xfmDir == 'xfmVar', 2, ifelse(xfmDir == 'xfmSbj', 1, NA))

    # z-scores - Z - correct
    if      (xfmMth == 'xfmZsc')
        xfmMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)),  "-"), xfmDim, as.vector(apply(dtaMtx, xfmDim, sd)),          "/")
    # range -1 to 1 - RANGE - over subjects correct, over variables possibly wrong in SPSS
    else if (xfmMth == 'xfmRNP')
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
    # range 0 to 1 - RESCALE - over subjects correct, over variables possibly wrong in SPSS
    else if (xfmMth == 'xfmRZP')
        xfmMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, min)),   "-"), xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
    # maximum magnitude of 1 - MAX - correct
    else if (xfmMth == 'xfmMag')
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, max)),  "/")
    # mean of 1 - MEAN - correct
    else if (xfmMth == 'xfmAvr')
        xfmMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)), "/")
    # standard deviation of 1 - SD - correct
    else if (xfmMth == 'xfmStd')
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
