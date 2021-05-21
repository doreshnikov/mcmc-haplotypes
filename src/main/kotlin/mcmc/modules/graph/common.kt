package mcmc.modules.graph

import utils.Randseed

typealias PathId = Int

fun Int.Companion.randomId(): PathId =
    Randseed.INSTANCE.randint()

typealias Pathway = List<Int>