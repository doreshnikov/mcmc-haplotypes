package mcmc.modules.graph

import graph.AlignedDeBruijnGraph
import graph.components.Edge
import mcmc.Model
import utils.Randseed
import kotlin.math.*

class PathsOverlay(
    entity: AlignedDeBruijnGraph,
    private val distributionConfig: DistributionConfig,
    val alpha: Double
//    private val paths: MutableList<Path>
) : Model<AlignedDeBruijnGraph, PathsOverlay.PathsDelta, List<Pair<String, Double>>>(entity) {

    private val weights = mutableMapOf<Edge, Double>()
    private val coverage = mutableMapOf<Edge, Double>()

    val paths = mutableMapOf<PathId, Path>()
    private val pathsBackRefs = mutableMapOf<Edge, MutableSet<PathId>>()

    private val common = mutableMapOf<PathId, MutableMap<PathId, Double>>()

    private fun newPathId(): PathId {
        var id = PathId.randomId()
        while (id in paths.keys) {
            id = PathId.randomId()
        }
        return id
    }

    inner class Path(
        val vertices: Pathway,
        var initialWeight: Double,
        val id: PathId = newPathId()
    ) {

        var weight = 0.0
            private set
        val edges = edges(vertices)
//        val impact = MutableList(vertices.size - 1) { 0.0 }

        var leftSeparationPoints = HashMap<Int, Int>()
        var rightSeparationPoints = HashMap<Int, Int>()

        init {
            vertices.forEachIndexed { i, v ->
                if (this@PathsOverlay.entity.edges[v].size > 1) {
                    rightSeparationPoints[v] = i
                }
                if (this@PathsOverlay.entity.reversedEdges[v].size > 1) {
                    leftSeparationPoints[v] = i
                }
            }
        }

        fun commit() {
            edges.forEach { e ->
                pathsBackRefs.getOrPut(e) { hashSetOf() }.add(id)
            }
            paths[id] = this
            if (initialWeight != 0.0) {
                plusAssign(initialWeight)
            }
            common[id] = mutableMapOf()
            edges.forEach { e ->
                pathsBackRefs[e]?.forEach {
                    if (it != id) {
                        val com = common.getValue(id).getOrPut(it) { 0.0 } + weights.getValue(e)
                        common[id]!![it] = com
                        common[it]!![id] = com
                    }
                }
            }
        }

        fun erase() {
            check(weight == 0.0) { "Can't delete paths with non-zero weights" }
            common.remove(id)
            edges.forEach { e ->
                pathsBackRefs[e]!!.remove(id)
                pathsBackRefs[e]!!.forEach {
                    common[it]!!.remove(id)
                }
                if (pathsBackRefs[e]!!.isEmpty()) {
                    pathsBackRefs.remove(e)
                }
            }
            paths.remove(id)
        }

        operator fun plusAssign(w: Double) {
            weight += w
//            edges.forEachIndexed { i, e ->
//                val c = coverage.compute(e) { _, c -> (c ?: 0.0) + w }!!
//                pathsBackRefs[e]?.forEach { pid ->
//                    val path = paths.getValue(pid)
//                    path.impact[i] = min(c - weights.getValue(e), path.weight)
//                }
//            }
        }

        operator fun minusAssign(w: Double) {
            plusAssign(-w)
        }

        fun collectDNA(): String {
            val sb = StringBuilder()
            sb.append(entity[vertices[0]].data.dropLast(1))
            vertices.forEach { v ->
                sb.append(entity[v].data.last())
            }
            return sb.toString()
        }

    }

    init {
        entity.indices.forEach { v ->
            entity.edges[v].forEach { u ->
                val locator = u.locator
                weights[locator] = u.weight
            }
        }

        fun dfs(v: Int): MutableList<Int> {
            if (entity.edges[v].isEmpty()) {
                return mutableListOf(v)
            }
            val next = entity.edges[v]
            val u = next[Randseed.INSTANCE.scoreIntegerEnum(next.map { it.weight }).sample()].target
            return dfs(u).also { it.add(v) }
        }

        val start = entity.indexOfFirst { it.position == 0 }
        Path(dfs(start).reversed(), 1.0).commit()
    }

    private fun checkIfPresent(vertices: Pathway): PathId? {
        val path = Path(vertices, 0.0)
        // possibly can be replaced with stream reduce, but the last attempt was painful
        val fittingIds = HashSet<PathId>(paths.keys)
        for (e in path.edges) {
            val followingPaths = pathsBackRefs[e]
            followingPaths?.let { fittingIds.retainAll(it) } ?: fittingIds.clear()
            if (fittingIds.isEmpty()) {
                break
            }
        }
        return fittingIds.firstOrNull()
    }

    private val distribution = Randseed.INSTANCE.poisson(distributionConfig.penaltyLambda)

    private fun edgeError(edge: Edge): Double {
        return coverage.getOrDefault(edge, 0.0) - weights.getValue(edge)
    }

    private fun logCoverageLikelihoodDelta(edge: Edge, delta: Double): Double {
        // e = c - w => -(c + d - w)^2 instead of -(c - w)^2 adds -d (2e + d)
        return -delta * (2 * edgeError(edge) + delta)
    }

    private fun logPathCountLikelihoodDelta(pathCountDelta: Int): Double {
        // 0 -> 0.0
        // else -> dist.logProb(paths.size + pathCountDelta) - dist.logProb(paths.size)
        // log(l^k e^-l / k!) = k logl - l - logk!
        val fact = paths.size.toDouble() + (1 + pathCountDelta) / 2
        return pathCountDelta * (ln(distributionConfig.penaltyLambda) - ln(fact))
    }

    abstract inner class PathsDelta : Delta()

    inner class New(
        private val sourcePath: Path,
        private val targetPathway: Pathway,
        private val delta: Double,
        override val logJumpDensity: Double
    ) : PathsDelta() {
        override val logLikelihoodDelta: Double
            get() {
                val targetEdges = edges(targetPathway)
                val edgesDelta = targetEdges.indices.sumByDouble { i ->
                    if (targetEdges[i] != sourcePath.edges[i]) {
                        logCoverageLikelihoodDelta(sourcePath.edges[i], -delta) +
                                logCoverageLikelihoodDelta(targetEdges[i], delta)
                    } else 0.0
                }
                return edgesDelta + logPathCountLikelihoodDelta(1)
            }

        override fun accept() {
            sourcePath -= delta
            Path(targetPathway, delta).commit()
        }
    }

    inner class Del(
        private val sourcePath: Path,
        private val targetPath: Path,
        override val logJumpDensity: Double
    ) : PathsDelta() {
        override val logLikelihoodDelta: Double
            get() {
                val delta = sourcePath.weight
                val edgesDelta = targetPath.edges.indices.sumByDouble { i ->
                    if (targetPath.edges[i] != sourcePath.edges[i]) {
                        logCoverageLikelihoodDelta(sourcePath.edges[i], -delta) +
                                logCoverageLikelihoodDelta(targetPath.edges[i], delta)
                    } else 0.0
                }
                return edgesDelta + logPathCountLikelihoodDelta(-1)
            }

        override fun accept() {
            targetPath += sourcePath.weight
            sourcePath -= sourcePath.weight
            sourcePath.erase()
        }
    }

    inner class Transfer(
        private val sourcePath: Path,
        private val targetPath: Path,
        private val delta: Double,
        override val logJumpDensity: Double
    ) : PathsDelta() {
        override val logLikelihoodDelta: Double
            get() {
                val edgesDelta = targetPath.edges.indices.sumByDouble { i ->
                    if (targetPath.edges[i] != sourcePath.edges[i]) {
                        logCoverageLikelihoodDelta(sourcePath.edges[i], -delta) +
                                logCoverageLikelihoodDelta(targetPath.edges[i], delta)
                    } else 0.0
                }
                return edgesDelta
            }

        override fun accept() {
            sourcePath -= delta
            targetPath += delta
        }
    }

    inner class Fail : PathsDelta() {
        override val logJumpDensity: Double = 0.0
        override val logLikelihoodDelta: Double = 0.0

        override fun accept() {}
    }

    override fun logLikelihood(): Double {
        val coverageLikelihood = weights.keys.sumByDouble { e ->
            -(coverage.getOrDefault(e, 0.0) - weights.getValue(e)).pow(2)
        }
        return coverageLikelihood + distribution.logProbability(paths.size)
    }

    companion object {
        // This all is also probably unusable now xD
        // kmp

        private fun transferOptionLogLikelihood(expectedDelta: Double, impactScore: Double): Double = when {
            impactScore == 0.0 -> -0.5
            sign(expectedDelta) != sign(impactScore) -> -1.0
            else -> min(1.0, impactScore / expectedDelta)
        }

        private fun impactScore(impact: Double): Double = when (impact) {
            0.0 -> 0.0
            else -> sign(impact) * ln(1.0 + abs(impact))
        }

        private fun scoreToDelta(score: Double): Double =
            sign(score) * (exp(abs(score)) - 1.0)

        fun ln(value: Int): Double = ln(value.toDouble())
        fun ln(value: Long): Double = ln(value.toDouble())
    }

    private val sourceSelector = Randseed.INSTANCE.probIntegerEnum(
        entity.sources.map { entity[it].pathCountRight }
    )
    private val uniformWeightSelector = Randseed.INSTANCE.uniformReal(0.0, 1.0)
    private val normalWeightSelector = Randseed.INSTANCE.normal(0.0, 1.0)

    private fun selectUniformWeight(from: Double, to: Double): Double {
        val w = uniformWeightSelector.sample()
        return from + w * (to - from)
    }

    private fun selectNormalWeight(mean: Double, sd: Double): Double {
        val w = normalWeightSelector.sample()
        return w * sd + mean
    }

    fun proposeCandidateNew(pNew: Double, pDel: Double): PathsDelta {
        @Deprecated("messes all densities")
        fun selectTargetPathway(): Pathway {
            val source = entity.sources[sourceSelector.sample()]
            val pathway = mutableListOf(source)

            while (!entity.isRightmost(pathway.last())) {
                val v = pathway.last()
                val next = entity[v].rightSelector.sample()
                pathway.add(entity.edges[v][next].target)
            }

            return pathway
        }

        fun selectClosePathway(path: Path): Pair<Pathway, Double> {
            val right = Randseed.INSTANCE.randint() % 2
            val next = MutableList(entity.size) { -1 }
            val score = MutableList(entity.size) { 0.0 }
            val prob = MutableList(entity.size) { 0.0 }
            var cur = -1
            val pathway = mutableListOf<Int>()
            var g = 0.0
            val common: Int

            if (right % 2 == 0) {
                val sep = path.rightSeparationPoints.entries.toList()
                val point = sep[Randseed.INSTANCE.uniformInteger(0, sep.size - 1).sample()]
                g -= ln(sep.size)
                pathway.addAll(path.vertices.take(point.value))
                common = point.value

                fun walk(v: Int, depth: Int) {
                    if (next[v] > -1) return
                    for (e in entity.edges[v]) {
                        walk(e.target, depth - 1)
                    }
                    if (entity.edges[v].isEmpty()) return

                    val opts = entity.edges[v].map { exp(score[it.target] / depth) }
                    val select = Randseed.INSTANCE.scoreIntegerEnum(opts)
                    val idx = select.sample()
                    val nedge = entity.edges[v][idx]
                    next[v] = nedge.target
                    score[v] = score[next[v]] - edgeError(nedge.locator)
                    prob[v] = select.logProbability(idx)
                }
                walk(point.key, path.edges.size - point.value)
                cur = point.key
                while (cur != -1) {
                    pathway.add(cur)
                    g += prob[cur]
                    cur = next[cur]
                }
            } else {
                val sep = path.leftSeparationPoints.entries.toList()
                val point = sep[Randseed.INSTANCE.uniformInteger(0, sep.size - 1).sample()]
                g -= ln(sep.size)
                common = path.edges.size - point.value

                fun walk(v: Int, depth: Int) {
                    if (next[v] > -1) return
                    for (e in entity.reversedEdges[v]) {
                        walk(e.target, depth - 1)
                    }
                    if (entity.reversedEdges[v].isEmpty()) return

                    val opts = entity.reversedEdges[v].map { exp(score[it.target] / depth) }
                    val select = Randseed.INSTANCE.scoreIntegerEnum(opts)
                    val idx = select.sample()
                    val nedge = entity.reversedEdges[v][idx]
                    next[v] = nedge.target
                    score[v] = score[next[v]] - edgeError(nedge.locator.reversed)
                    prob[v] = select.logProbability(idx)
                }
                walk(point.key, point.value)
                cur = point.key
                while (cur != -1) {
                    pathway.add(cur)
                    g += prob[cur]
                    cur = next[cur]
                }
                pathway.reverse()
                pathway.addAll(path.vertices.drop(point.value + 1))
            }
            return pathway to g / (path.edges.size - common) + ln(0.5)
        }

        var gForward = 0.0
        val pathIds = paths.keys.filter {
            paths.getValue(it).weight >= 2 * alpha
        }.toMutableList()
        if (pathIds.isEmpty()) {
            return Fail()
        }
        val dropPathIds = Randseed.INSTANCE.randint() % 2 == 0
        if (dropPathIds) {
            Randseed.INSTANCE.shuffle(pathIds)
            val ok = 5000000 / entity.originLength
            if (ok < pathIds.size) {
                pathIds.dropLast(pathIds.size - ok)
            }
        }

        val ok = 5000000 / pathIds.size
        val weights = pathIds.map { pid ->
            val path = paths.getValue(pid)
            val delta = if (dropPathIds && ok < path.edges.size) {
                var accum = 0.0
                val dist = Randseed.INSTANCE.uniformInteger(0, path.edges.size - 1)
                repeat(ok) {
                    val idx = dist.sample()
                    accum -= edgeError(path.edges[idx])
                }
                accum / ok
            } else {
                -path.edges.sumByDouble { e -> edgeError(e) } / path.edges.size
            }
            if (delta > 0) min(delta, path.weight) else delta
        }
        val pathSelector = Randseed.INSTANCE.scoreIntegerEnum(weights.map(::exp))
        val idx = pathSelector.sample()
        gForward += pathSelector.logProbability(idx)
        val sourcePath = paths.getValue(pathIds[idx])

//        var targetPathway = selectTargetPathway()
        var (targetPathway, g) = selectClosePathway(sourcePath)
        var retries = 5
        while (checkIfPresent(targetPathway) != null) {
//            targetPathway = selectTargetPathway()
            val res = selectClosePathway(sourcePath)
            targetPathway = res.first
            g = res.second
            retries--
            if (retries == 0) {
                return Fail()
            }
        }
        gForward += g

//        val delta = Randseed.INSTANCE.uniformReal(alpha, sourcePath.weight - alpha).sample()
        val wS = selectNormalWeight(0.0, 1.0)
        var delta = wS + weights[idx]
        var gW = normalWeightSelector.logDensity(wS)
        if (delta < alpha) {
            gW = ln(normalWeightSelector.cumulativeProbability(alpha - weights[idx]))
            delta = alpha
        } else if (delta > sourcePath.weight - alpha) {
            gW = ln(1.0 - normalWeightSelector.cumulativeProbability(sourcePath.weight - alpha - weights[idx]))
            delta = sourcePath.weight - alpha
        }
        gForward += gW

        var gBackward = 0.0
        gBackward = gForward // todo there's still the same problem as before
        val gTotal = gBackward - gForward + ln(pDel) - ln(pNew)

        return New(
            sourcePath, targetPathway, delta,
            selectNormalWeight(gTotal / 2, abs(gTotal))
        )
    }

    fun proposeCandidateDel(pDel: Double, pNew: Double): PathsDelta {
        val pathIds = paths.keys.toList()
        val weights = pathIds.map { pid ->
            val path = paths.getValue(pid)
            val delta = path.edges.sumByDouble { e -> edgeError(e) } / path.edges.size
            delta / path.weight
        }
        val pathSelector = Randseed.INSTANCE.scoreIntegerEnum(weights)

        var gForward = 0.0
        val idx = pathSelector.sample()
        val sourcePath = paths.getValue(pathIds[idx])
        gForward += pathSelector.logProbability(idx)

        val targetWeights = pathIds.map {
            common[sourcePath.id]!![it] ?: 0.0
        }
        val targetPathSelector = Randseed.INSTANCE.scoreIntegerEnum(targetWeights)
        val targetIdx = targetPathSelector.sample()
        val targetPathId = pathIds[targetIdx]
        val targetPath = paths.getValue(targetPathId)
        gForward += targetPathSelector.logProbability(targetIdx)

        var gBackward = 0.0
        gBackward = gForward // todo same problem - density of New << Del
        val gTotal = gBackward - gForward + ln(pNew) - ln(pDel)

        return Del(
            sourcePath, targetPath,
            selectNormalWeight(gTotal / 2, abs(gTotal))
        )
    }

    fun proposeCandidateTransfer(): PathsDelta {
        val weights = paths.mapValues {
            val path = it.value
            val delta = -path.edges.sumByDouble { e -> edgeError(e) } / path.edges.size
            if (delta > 0) min(delta, path.weight) else delta
        }

        val sourcePathIds = paths.keys.filter {
            paths.getValue(it).weight > alpha
        }
        if (sourcePathIds.isEmpty()) return Fail()
        val sourceWeights = sourcePathIds.map { weights.getValue(it) }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourceWeights.map(::exp))
        val sourcePath = paths.getValue(sourcePathIds[sourcePathSelector.sample()])

        val targetPathIds = paths.keys.toList()
        val targetWeights = targetPathIds.map { weights.getValue(it) }
        val pathSelector = Randseed.INSTANCE.scoreIntegerEnum(targetWeights.map { exp(-it) })
        var targetPathId = targetPathIds[pathSelector.sample()]
        while (targetPathId == sourcePath.id) {
            targetPathId = targetPathIds[pathSelector.sample()]
        }
        val targetPath = paths.getValue(targetPathId)


//        val delta = Randseed.INSTANCE.uniformReal(0.0, sourcePath.weight - alpha).sample()
        val expect = (weights.getValue(targetPathId) - weights.getValue(sourcePath.id)) / 2
        val sd = max(alpha, (weights.getValue(targetPathId) + weights.getValue(sourcePath.id)) / 2)
        val wSelector = Randseed.INSTANCE.normal(expect, sd)
        var delta = wSelector.sample()
        if (delta < 0) {
            delta = min(sd, (sourcePath.weight - alpha) / 2)
        }
        if (delta > sourcePath.weight - alpha) {
            delta = sourcePath.weight - alpha
        }
//        var newPotentialSources = 0
//        if (targetPath.weight <= alpha && targetPath.weight + delta > alpha) {
//            newPotentialSources += 1
//        }
//        if (sourcePath.weight - delta <= alpha) {
//            newPotentialSources -= 1
//        }

        return Transfer(
            sourcePath, targetPath, delta,
            ln(sourcePath.weight - alpha) - ln(targetPath.weight + delta - alpha)
        ) // todo this g/g also needs tuning ofc
    }

    private val proposer = Randseed.INSTANCE.uniformReal(0.0, 1.0)

    override fun proposeCandidate(): Delta {
        fun pDel(pathCount: Int): Double {
            return if (pathCount > 1) 0.5 / distributionConfig.penaltyLambda else 0.0
        }

        fun pNew(pathCount: Int): Double {
//            if (pathCount > 20) return 0.0
            return if (pathCount > 1) 0.5 / max(pathCount, 3) else 1.0
        }

        val pNewCurrent = pNew(paths.size)
        val pDelCurrent = pDel(paths.size)

        val option = proposer.sample()
        val result = when {
            option <= pNewCurrent ->
                proposeCandidateNew(pNewCurrent, pDel(paths.size + 1))
            option <= pNewCurrent + pDelCurrent ->
                proposeCandidateDel(pDelCurrent, pNew(paths.size - 1))
            else ->
                proposeCandidateTransfer()
        }
        return if (result is Fail) proposeCandidate() else result
    }

    override fun extractResult(): List<Pair<String, Double>> {
        return paths.map { it.value.collectDNA() to it.value.weight }
    }

}