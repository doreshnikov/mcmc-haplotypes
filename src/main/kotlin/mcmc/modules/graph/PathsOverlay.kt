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
) : Model<AlignedDeBruijnGraph, PathsOverlay.PathsDelta>(entity) {

    private val weights = mutableMapOf<Edge, Double>()
    private val coverage = mutableMapOf<Edge, Double>()

    private val paths = mutableMapOf<PathId, Path>()
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
        private val initialWeight: Double,
        val id: PathId = newPathId()
    ) {

        var weight = 0.0
            private set
        val edges = edges(vertices)
        val impact = MutableList(vertices.size - 1) { 0.0 }

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
        }

        fun erase() {
            check(weight == 0.0) { "Can't delete paths with non-zero weights" }
            edges.forEach { e ->
                pathsBackRefs[e]?.remove(id)
            }
            paths.remove(id)
        }

        operator fun plusAssign(w: Double) {
            weight += w
            edges.forEachIndexed { i, e ->
                val c = coverage.compute(e) { _, c -> (c ?: 0.0) + w }!!
                pathsBackRefs[e]?.forEach { pid ->
                    val path = paths.getValue(pid)
                    path.impact[i] = min(c - weights.getValue(e), path.weight)
                }
            }
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
            val options = entity.edges[v].map { Pair(it.target, it.weight) }
            val u = Randseed.INSTANCE.probEnum(options).sample()
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

    private fun logCoverageLikelihoodDelta(edge: Edge, delta: Double): Double {
        val error = coverage.getOrDefault(edge, 0.0) - weights.getValue(edge)
        // e = c - w => -(c + d - w)^2 instead of -(c - w)^2 adds -d (2e + d)
        return -delta * (2 * error + delta)
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

    override fun logLikelihood(): Double {
        val coverageLikelihood = weights.keys.sumByDouble { e ->
            -(coverage.getOrDefault(e, 0.0) - weights.getValue(e)).pow(2)
        }
        return coverageLikelihood + distribution.logProbability(paths.size)
    }

    fun collectResult(): List<Pair<String, Double>> {
        return paths.values.map { path -> path.collectDNA() to path.weight }
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

    private val vertexSelector = Randseed.INSTANCE.probIntegerEnum(
        entity.sources.map { entity[it].pathCountRight.toDouble() }
    )
    private val totalPaths = entity.sources.sumOf { entity[it].pathCountRight }
    private val weightSelector = Randseed.INSTANCE.uniformReal(0.0, 1.0)

    private fun selectWeight(from: Double, to: Double): Double {
        val w = weightSelector.sample()
        return from + w * (to - from)
    }

    fun proposeCandidateNew(pNew: Double, pDel: Double): New {
        fun selectTargetPathway(): Pathway {
            val source = entity.sources[vertexSelector.sample()]
            val pathway = mutableListOf(source)

            while (!entity.isRightmost(pathway.last())) {
                val v = pathway.last()
                // todo can incorporate distributions into vertices to avoid re-initializing every time
                val next = entity[v].rightSelector.sample()
                pathway.add(entity.edges[v][next].target)
            }

            return pathway
        }

        val pathIds = paths.keys.filter {
            paths.getValue(it).weight >= 2 * alpha
        }
        val pathSelector = Randseed.INSTANCE.uniformInteger(0, pathIds.size - 1)
        val sourcePath = paths.getValue(pathIds[pathSelector.sample()])

        var targetPathway = selectTargetPathway()
        while (checkIfPresent(targetPathway) != null) {
            targetPathway = selectTargetPathway()
        }
        val delta = Randseed.INSTANCE.uniformReal(alpha, sourcePath.weight - alpha).sample()

        return New(
            sourcePath, targetPathway, delta,
            ln(pathIds.size) + ln(totalPaths - paths.size) + ln(sourcePath.weight - 2 * alpha) -
                    ln(paths.size) - ln(paths.size + 1) +
                    ln(pDel) - ln(pNew)
        )
    }

    fun proposeCandidateDel(pDel: Double, pNew: Double): Del {
        val pathIds = paths.keys.toList()
        val pathSelector = Randseed.INSTANCE.uniformInteger(0, pathIds.size - 1)

        val sourcePath = paths.getValue(pathIds[pathSelector.sample()])
        var targetPathId = pathIds[pathSelector.sample()]
        while (targetPathId == sourcePath.id) {
            targetPathId = pathIds[pathSelector.sample()]
        }
        val targetPath = paths.getValue(targetPathId)

        val backwardsPathIds = paths.keys.filter {
            paths.getValue(it).weight >= 2 * alpha && it != sourcePath.id
        }
        return Del(
            sourcePath, targetPath,
            ln(paths.size) + ln(paths.size - 1) -
                    ln(backwardsPathIds.size) - ln(totalPaths - paths.size - 1) -
                    ln(sourcePath.weight + targetPath.weight - 2 * alpha) +
                    ln(pNew) - ln(pDel)
        )
    }

    fun proposeCandidateTransfer(): Transfer {
        val sourcePathIds = paths.keys.filter {
            paths.getValue(it).weight > alpha
        }
        val sourcePathSelector = Randseed.INSTANCE.uniformInteger(0, sourcePathIds.size - 1)
        val sourcePath = paths.getValue(sourcePathIds[sourcePathSelector.sample()])

        val pathIds = paths.keys.toList()
        val pathSelector = Randseed.INSTANCE.uniformInteger(0, pathIds.size - 1)
        var targetPathId = pathIds[pathSelector.sample()]
        while (targetPathId == sourcePath.id) {
            targetPathId = pathIds[pathSelector.sample()]
        }
        val targetPath = paths.getValue(targetPathId)

        val delta = Randseed.INSTANCE.uniformReal(0.0, sourcePath.weight - alpha).sample()
        var newPotentialSources = 0
        if (targetPath.weight <= alpha && targetPath.weight + delta > alpha) {
            newPotentialSources += 1
        }
        if (sourcePath.weight - delta <= alpha) {
            newPotentialSources -= 1
        }

        return Transfer(
            sourcePath, targetPath, delta,
            ln(sourcePathIds.size) + ln(sourcePath.weight - alpha) -
                    ln(sourcePathIds.size + newPotentialSources) - ln(targetPath.weight + delta - alpha)
        )
    }

    override fun proposeCandidate(): Delta {
        fun pDel(pathCount: Int): Double {
            return if (pathCount > 1) 1.0 / distributionConfig.penaltyLambda else 0.0
        }

        fun pNew(pathCount: Int): Double {
            return if (pathCount > 1) 1.0 / (paths.size + distributionConfig.penaltyLambda) else 1.0
        }

        val pNewCurrent = pNew(paths.size)
        val pDelCurrent = pDel(paths.size)

        val option = Randseed.INSTANCE.uniformReal(0.0, 1.0).sample()
        return when {
            option <= pNewCurrent ->
                proposeCandidateNew(pNewCurrent, pDel(paths.size + 1))
            option <= pNewCurrent + pDelCurrent ->
                proposeCandidateDel(pDelCurrent, pNew(paths.size - 1))
            else ->
                proposeCandidateTransfer()
        }
    }

    @Deprecated("Old version")
    fun proposeCandidateThisDoesntWork(): PathsDelta {
        val averageImpacts = paths.mapValues { (_, p) ->
            p.impact.sumByDouble(::impactScore) / p.impact.size
        }
        val entries = averageImpacts.entries
        val pid = Randseed.INSTANCE.scoreEnum(entries.map { it.key }, entries.map { it.value }).sample()
        val sourcePath = paths.getValue(pid)

        val prefixes = mutableListOf(0.0)
        for (i in sourcePath.impact) {
            prefixes.add(prefixes.last() + i)
        }
        val suffixes = prefixes.map { prefixes.last() - it }

        val prefixAverages = (1 until prefixes.size - 1).map { i -> prefixes[i] / i }
        val suffixAverages = (1 until suffixes.size - 1).map { i -> suffixes[i] / (suffixes.size - i) }
        var splitPoint = Randseed.INSTANCE.scoreIntegerEnum((prefixAverages + suffixAverages).map(::abs)).sample()

        if (splitPoint >= prefixAverages.size) {
            splitPoint -= prefixAverages.size
            val delta = min(scoreToDelta(suffixAverages[splitPoint]), sourcePath.weight)
            splitPoint++

            fun dfs(v: Int): Pair<Double, MutableList<Int>> {
                val options = entity.edges[v].map { e ->
                    val tail = dfs(e.target)
                    val impact = min(coverage.getOrDefault(e.locator, 0.0) - e.weight, e.weight)
                    impactScore(impact) + tail.first to tail.second
                }
                if (options.isEmpty()) {
                    return 0.0 to mutableListOf(v)
                }
                return Randseed.INSTANCE.scoreEnum(
                    options,
                    options.map { exp(transferOptionLogLikelihood(delta, it.first)) }
                ).sample().also { it.second.add(v) }
            }

            val tail = dfs(sourcePath.vertices[splitPoint])
            tail.second.reverse()
//            return Transfer(
//                sourcePath,
//                sourcePath.vertices.take(splitPoint) + tail.second,
//                sign(delta) * Randseed.INSTANCE.uniformReal(
//                    0.0, min(abs(delta), abs(scoreToDelta(tail.first)))
//                ).sample()
//            )
        } else {
            val delta = min(scoreToDelta(prefixAverages[splitPoint]), sourcePath.weight)
            splitPoint++

            fun dfs(v: Int): Pair<Double, MutableList<Int>> {
                val options = entity.reversedEdges[v].map { e ->
                    val head = dfs(e.target)
                    val impact = min(coverage.getOrDefault(e.locator.reversed, 0.0) - e.weight, e.weight)
                    impactScore(impact) + head.first to head.second
                }
                if (options.isEmpty()) {
                    return 0.0 to mutableListOf(v)
                }
                return Randseed.INSTANCE.scoreEnum(
                    options,
                    options.map { exp(transferOptionLogLikelihood(delta, it.first)) }
                ).sample().also { it.second.add(v) }
            }

            val head = dfs(sourcePath.vertices[splitPoint])
//            return Transfer(
//                sourcePath,
//                head.second + sourcePath.vertices.drop(splitPoint + 1),
//                sign(delta) * Randseed.INSTANCE.uniformReal(
//                    0.0, min(abs(delta), abs(scoreToDelta(head.first)))
//                ).sample()
//            )
        }
        error("This function is deprecated")
    }

}