library(ordinal)

stopifnot(isTRUE(all.equal(ordinal:::pgumbelR(-10:10), ordinal:::pgumbel(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::pgumbel2R(-10:10), ordinal:::pgumbel2(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::pAOR(-10:10, 1), ordinal:::pAO(-10:10, 1))))
stopifnot(isTRUE(all.equal(ordinal:::plgammaR(-10:10, 1), ordinal:::plgamma(-10:10, 1))))

stopifnot(isTRUE(all.equal(ordinal:::dgumbelR(-10:10), ordinal:::dgumbel(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::dgumbel2R(-10:10), ordinal:::dgumbel2(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::dAOR(-10:10, 1), ordinal:::dAO(-10:10, 1))))
stopifnot(isTRUE(all.equal(ordinal:::dlgammaR(-10:10, 1), ordinal:::dlgamma(-10:10, 1))))

stopifnot(isTRUE(all.equal(ordinal:::ggumbelR(-10:10), ordinal:::ggumbel(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::ggumbel2R(-10:10), ordinal:::ggumbel2(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::glogisR(-10:10), ordinal:::glogis(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::gnormR(-10:10), ordinal:::gnorm(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::gcauchyR(-10:10), ordinal:::gcauchy(-10:10))))
stopifnot(isTRUE(all.equal(ordinal:::gAOR(-10:10, 1), ordinal:::gAO(-10:10, 1))))
stopifnot(isTRUE(all.equal(ordinal:::glgammaR(-10:10, 1), ordinal:::glgamma(-10:10, 1))))

## Do we need a test suite for PFUN, DFUN and GFUN as well?
