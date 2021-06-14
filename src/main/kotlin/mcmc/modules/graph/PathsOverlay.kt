package mcmc.modules.graph

import graph.AlignedDeBruijnGraph
import graph.components.Edge
import mcmc.Model
import utils.Randseed
import java.io.File
import kotlin.math.*

class PathsOverlay(
    entity: AlignedDeBruijnGraph,
    private val distributionConfig: DistributionConfig,
    val alpha: Double,
    val maxHandleLength: Int
//    private val paths: MutableList<Path>
) : Model<AlignedDeBruijnGraph, PathsOverlay.PathsDelta, List<Pair<String, Double>>>(entity) {

    private val VERBOSE = true
    val SHIFT_NEW = -25.95
    val SHIFT_DELETE = 87.64 // thats for 500.45
    val SCALE_NEW = 30.0
    val SCALE_DELETE = 80.0
//    val SHIFT_NEW = -32.81
//    val SHIFT_DELETE = 200.88 // 500.40
//    val SCALE_NEW = 30.0
//    val SCALE_DELETE = 500.0
//    val SHIFT_NEW = -2.0
//    val SHIFT_DELETE = 2.5
//    val SCALE_NEW = 3.0
//    val SCALE_DELETE = 4.0
    val log = File("stats").bufferedWriter()

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
        val totalEdgesWeight: Double
        var totalCoverage = 0.0
            private set
        private val verticesSet = HashSet(vertices)

        val excess: Double
            get() = min(weight, (totalCoverage - totalEdgesWeight) / edges.size)

        val edges = edges(vertices)

        val leftSeparationPoints = HashMap<Int, Int>()
        val rightSeparationPoints = HashMap<Int, Int>()

        val allHandles = mutableListOf<Pathway>()
        val handlesGraph = ArrayList<MutableSet<Int>>()

        init {
            vertices.forEachIndexed { i, v ->
                if (this@PathsOverlay.entity.edges[v].size > 1) {
                    rightSeparationPoints[v] = i
                }
                if (this@PathsOverlay.entity.reversedEdges[v].size > 1) {
                    leftSeparationPoints[v] = i
                }
            }

            totalEdgesWeight = edges.sumOf { e -> weights.getValue(e) }
            totalCoverage = edges.sumOf { e -> coverage.getOrDefault(e, 0.0) }
        }

        fun collectHandles() {
            if (allHandles.isNotEmpty()) {
                return
            }
            val blacklist = mutableSetOf<Int>()

            fun walkLeft(stack: MutableList<Int>) {
                val v = stack.last()
                if (stack.size > maxHandleLength + 1) {
                    return
                }
                if (entity[v].position == 0) {
                    allHandles.add(stack.reversed())
                }
                for (e in entity.reversedEdges[v]) {
                    if (e.target !in verticesSet) {
                        stack.add(e.target)
                        walkLeft(stack)
                        stack.removeLast()
                    }
                }
            }

            fun walkRight(stack: MutableList<Int>): Int {
                if (stack.size > maxHandleLength + 1) {
                    return 0
                }
                val v = stack.last()
                if (stack.size > 1 && (entity.isRightmost(v) || v in verticesSet)) {
                    if ((entity.isRightmost(v) && v !in verticesSet) || stack.size > 2) {
                        allHandles.add(ArrayList(stack))
                        return 1
                    }
                    return 0
                } else {
                    var found = 0
                    for (e in entity.edges[v]) {
                        if (e.target !in blacklist) {
                            stack.add(e.target)
                            val go = walkRight(stack)
                            stack.removeLast()

                            if (go == 0) {
                                blacklist.add(e.target)
                            }
                            found += go
                        }
                    }
                    return found
                }
            }

            vertices.forEachIndexed { idx, v ->
                blacklist.clear()
                if (idx <= maxHandleLength && entity.reversedEdges[v].isNotEmpty()) {
                    walkLeft(mutableListOf(v))
                }
                if (entity.edges[v].isNotEmpty()) {
                    walkRight(mutableListOf(v))
                }
            }
//            if (VERBOSE) {
            if (false) {
                println(allHandles.size)
                allHandles.groupBy { it.size }.forEach { (t, u) ->
                    println("- size $t : $u")
                }
            }
            allHandles.sortBy { entity[it.first()].position }

            handlesGraph.addAll(allHandles.map { mutableSetOf() })
            var r = 0
            allHandles.forEachIndexed { idx, handle ->
                while (r < allHandles.size &&
                    entity[allHandles[r].first()].position <= entity[handle.last()].position
                ) {
                    r++
                }
                val start = entity[handle.first()].position
                for (i in 0 until r) {
                    if (entity[allHandles[i].last()].position >= start && i != idx) {
                        handlesGraph[idx].add(i)
                    }
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
                pathsBackRefs.getValue(e).remove(id)
                pathsBackRefs.getValue(e).forEach {
                    common[it]!!.remove(id)
                }
                if (pathsBackRefs[e].isNullOrEmpty()) {
                    pathsBackRefs.remove(e)
                }
            }
            paths.remove(id)
        }

        operator fun plusAssign(w: Double) {
            weight += w
            // TODO calculate common edges and do it faster
            edges.forEach { e ->
                coverage.compute(e) { _, v -> (v ?: 0.0) + w }
                pathsBackRefs[e]?.forEach { pid ->
                    val path = paths.getValue(pid)
                    path.totalCoverage += w
                }
            }
        }

        operator fun minusAssign(w: Double) {
            plusAssign(-w)
        }

        fun edgeExcess(idx: Int): Double {
            return min(weight, edgeError(edges[idx]))
        }

        fun collectDNA(): String {
            val sb = StringBuilder()
            sb.append(entity[vertices[0]].data.dropLast(1))
            vertices.forEach { v ->
                sb.append(entity[v].data.last())
            }
            return sb.toString()
        }

        fun differenceWith(pid: PathId): Double {
            return totalEdgesWeight + paths.getValue(pid).totalEdgesWeight - 2 * common.getValue(id)
                .getOrDefault(pid, 0.0)
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
        return when (pathCountDelta) {
            0 -> 0.0
            else -> distribution.logProbability(paths.size + pathCountDelta) - distribution.logProbability(paths.size)
        }
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
            Path(targetPathway, delta).apply {
//                collectHandles()
                commit()
            }
        }
    }

    inner class Move(
        private val sourcePath: Path,
        private val targetPathway: Pathway,
        override val logJumpDensity: Double
    ) : PathsDelta() {
        override val logLikelihoodDelta: Double
            get() {
                val targetEdges = edges(targetPathway)
                return targetEdges.indices.sumByDouble { i ->
                    if (targetEdges[i] != sourcePath.edges[i]) {
                        logCoverageLikelihoodDelta(sourcePath.edges[i], -sourcePath.weight) +
                                logCoverageLikelihoodDelta(targetEdges[i], sourcePath.weight)
                    } else 0.0
                }
            }

        override fun accept() {
            Path(targetPathway, sourcePath.weight).apply {
//                collectHandles()
                commit()
            }
            sourcePath -= sourcePath.weight
            sourcePath.erase()
        }
    }

    inner class Delete(
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

        fun soft(x: Double): Double {
            return 1 + sign(x) * sqrt(abs(x))
        }
    }

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

    fun proposeCandidateNewLast(pNew: Double, pDelete: Double): PathsDelta {
        val sourcePathIds = paths.keys
            .filter { paths.getValue(it).weight >= 2 * alpha }
        val sourcePathScores = sourcePathIds
            .map { soft(paths.getValue(it).excess) }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourcePathScores)
        val sourceIdx = sourcePathSelector.sample()
        val sourcePath = paths.getValue(sourcePathIds[sourceIdx])

        sourcePath.collectHandles()
        val handleGraph = sourcePath.handlesGraph
        val handles = sourcePath.allHandles
        if (handles.isEmpty()) {
            return Fail()
        }

        val selectedHandles = mutableListOf<Int>()
        val selector = Randseed.INSTANCE.uniformReal(0.0, 1.0)
        val t = handleGraph.count { it.isNotEmpty() }
        var xi = (handles.size - t / 2).toDouble() / 2.0 * handles.size / handles.sumOf { it.size - 1 }
        xi = max(xi, 1.0 / handles.minOf { it.size - 2 })
        handles.forEachIndexed { idx, handle ->
            if (selector.sample() < 1.0 / xi / (handle.size - 1)) {
                selectedHandles.add(idx)
            }
        }
        val newGraph = MutableList<MutableSet<Int>>(selectedHandles.size) { mutableSetOf() }
        var edgesCount = 0
        selectedHandles.forEachIndexed { i, h1 ->
            selectedHandles.forEachIndexed { j, h2 ->
                if (h1 in handleGraph[h2]) {
                    newGraph[j].add(i)
                    edgesCount++
                }
            }
        }
        if (selectedHandles.isEmpty()) {
            return Fail()
        }
        val removed = MutableList(selectedHandles.size) { false }
        val indexSelector = Randseed.INSTANCE.uniformInteger(0, selectedHandles.size - 1)
        while (edgesCount > 0) {
            val i = indexSelector.sample()
            if (removed[i]) {
                continue
            }
            removed[i] = true
            for (j in newGraph[i]) {
                newGraph[j].remove(i)
            }
            edgesCount -= 2 * newGraph[i].size
            newGraph[i].clear()
        }

        val resultHandles = selectedHandles.filterIndexed { i, _ -> !removed[i] }
        val resultHandlesNeighbors = resultHandles.flatMap { handleGraph[it] }.toSet()
        val targetPathway = mutableListOf<Int>()
        resultHandles.forEach {
            val handle = handles[it]
            targetPathway += sourcePath.vertices.take(entity[handle.first()].position).drop(targetPathway.size)
            targetPathway += handle
        }
        targetPathway += sourcePath.vertices.drop(targetPathway.size)
        if (checkIfPresent(targetPathway) != null) {
            return Fail()
        }

        var pathsDifference = 0.0
        val targetEdges = edges(targetPathway)
        val w0 = targetEdges.mapIndexed { i, edge ->
            if (sourcePath.edges[i] != edge) {
                pathsDifference += weights.getValue(sourcePath.edges[i]) + weights.getValue(edge)
                sourcePath.edgeExcess(i) - edgeError(edge)
            } else 0.0
        }.sum() / 2.0 / sourcePath.edges.size
        var mean = (sourcePath.weight + w0 * 3.0) / 5.0
        if (mean < 0.0) {
            mean = alpha
        }
        var std = max(mean - alpha, sourcePath.weight - mean - alpha)
        if (std < alpha) {
            std = alpha
        }
        val weightSelector = Randseed.INSTANCE.normal(mean, std)
        val weight = weightSelector.sample()
        val correctedWeight = when {
            weight < alpha -> alpha
            weight > sourcePath.weight - alpha -> sourcePath.weight - alpha
            else -> weight
        }

        val gForward = ln(pNew) +
                sourcePathSelector.logProbability(sourceIdx) +
                resultHandles.sumOf { -ln(xi) - ln(handles[it].size - 1) } +
                handles.indices.toHashSet().minus(resultHandles).sumOf {
                    val d = 1.0 / xi / (handles[it].size - 1)
                    ln(if (it in resultHandlesNeighbors) 1.0 - d / 2.0 else 1.0 - d)
                } +
                weightSelector.logDensity(weight) - weightSelector.logDensity(mean)
        val targetExcess = min(
            correctedWeight,
            targetEdges.sumOf { edgeError(it) } / (targetPathway.size - 1)
        )
        val totalRevDifference = paths.values
            .filter { it.id != sourcePath.id }
            .sumOf { 1.0 / sqrt(it.differenceWith(sourcePath.id)) } +
                1.0 / sqrt(pathsDifference)
        val gBackward = ln(pDelete) +
                ln(soft(targetExcess)) - ln(paths.values.sumOf { soft(it.excess) }) +
                -ln(pathsDifference) / 2.0 - ln(totalRevDifference)
        val gDiv = gBackward - gForward
        if (VERBOSE)
            log.write("New: $gDiv\n")

        return New(
            sourcePath, targetPathway, correctedWeight,
            (gDiv + SHIFT_NEW) / SCALE_NEW
        )
    }

    fun proposeCandidateMoveLast(): PathsDelta {
        val sourcePathIds = paths.keys.toList()
        val sourcePathScores = sourcePathIds
            .map { soft(paths.getValue(it).excess) }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourcePathScores)
        val sourceIdx = sourcePathSelector.sample()
        val sourcePath = paths.getValue(sourcePathIds[sourceIdx])

        sourcePath.collectHandles()
        val handles = sourcePath.allHandles
        if (handles.isEmpty()) {
            return Fail()
        }
        val handleScores = handles.map { handle ->
            var excessMetric = -edges(handle).sumOf { e -> edgeError(e) }
            for (pos in entity[handle.first()].position until entity[handle.last()].position) {
                excessMetric += sourcePath.edgeExcess(pos)
            }
            soft(excessMetric / (handle.size - 1) / 2.0) / (handle.size - 1).toDouble()
        }
        val handleSelector = Randseed.INSTANCE.scoreIntegerEnum(handleScores)
        val handleIdx = handleSelector.sample()
        val handle = handles[handleIdx]

        val targetPathway = mutableListOf<Int>()
        targetPathway += sourcePath.vertices.take(entity[handle.first()].position)
        targetPathway += handle
        targetPathway += sourcePath.vertices.drop(targetPathway.size)
        if (checkIfPresent(targetPathway) != null) {
            return Fail()
        }

        var backHandleScore = edges(handle).sumOf { min(sourcePath.weight, edgeError(it)) } / (handle.size - 1)
        for (pos in entity[handle.first()].position until entity[handle.last()].position) {
            backHandleScore -= edgeError(sourcePath.edges[pos])
        }
        backHandleScore /= (handle.size) * 2.0
        return Move(
            sourcePath, targetPathway,
            ln(soft(backHandleScore) / 2.0 / (handle.size - 1)) - ln(handleScores[handleIdx])
        )
    }

    fun proposeCandidateDeleteLast(pDelete: Double, pNew: Double): PathsDelta {
        val sourcePathIds = paths.keys.toList()
        val sourcePathScores = sourcePathIds
            .map { soft(paths.getValue(it).excess) }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourcePathScores)
        val sourceIdx = sourcePathSelector.sample()
        val sourcePath = paths.getValue(sourcePathIds[sourceIdx])

        val targetPathIds = paths.keys
            .filter { it != sourcePath.id }
            .filter { common.getValue(sourcePath.id).getOrDefault(it, 0.0) > 0.0 }
        if (targetPathIds.isEmpty()) {
            return Fail()
        }
        val targetPathScores = targetPathIds
            .map { soft(-paths.getValue(it).excess) / sourcePath.differenceWith(it) }
        val targetPathSelector = Randseed.INSTANCE.scoreIntegerEnum(targetPathScores)
        val targetIdx = targetPathSelector.sample()
        val targetPath = paths.getValue(targetPathIds[targetIdx])

        val gForward = ln(pDelete) +
                sourcePathSelector.logProbability(sourceIdx) +
                targetPathSelector.logProbability(targetIdx)

        targetPath.collectHandles()
        val handles = targetPath.allHandles
        val t = targetPath.handlesGraph.count { it.isNotEmpty() }
        var xi = (handles.size - t / 2).toDouble() / 2.0 * handles.size / handles.sumOf { it.size - 1 }
        xi = max(xi, 1.0 / handles.minOf { it.size - 2 })
        val covered = MutableList(targetPath.edges.size) { false }
        val resultHandles = handles.indices.filter { i ->
            val handle = handles[i]
            if (edges(handle).count { pathsBackRefs[it]?.contains(sourcePath.id) == true } == handle.size - 1) {
                for (pos in entity[handle.first()].position until entity[handle.last()].position) {
                    covered[pos] = true
                }
                true
            } else false
        }
        val extraHandlesSizes = mutableListOf<Int>()
        var idx = 0
        while (idx < targetPath.edges.size) {
            while (idx < targetPath.edges.size && covered[idx]) {
                idx++
            }
            val start = idx
            while (idx < targetPath.edges.size && !covered[idx]) {
                idx++
            }
            if (idx < start) {
                extraHandlesSizes.add(idx - start)
            }
        }
        val resultHandlesNeighbors = handles.indices.filter { i ->
            val handle = handles[i]
            var intersects = false
            for (pos in entity[handle.first()].position until entity[handle.last()].position) {
                if (covered[pos]) {
                    intersects = true
                    break
                }
            }
            intersects
        }
        val gBackward = ln(pNew) +
                sourcePathSelector.logProbability(sourcePathIds.indexOf(targetPath.id)) +
                resultHandles.sumOf { -ln(xi) - ln(handles[it].size - 1) } +
                extraHandlesSizes.sumOf { -ln(xi) - ln(it) } +
                handles.indices.toHashSet().minus(resultHandles).sumOf {
                    val d = 1.0 / xi / (handles[it].size - 1)
                    ln(if (it in resultHandlesNeighbors) 1.0 - d / 2.0 else 1.0 - d)
                }

        val gDiv = gBackward - gForward
        if (VERBOSE)
            log.write("Delete: $gDiv\n")

        return Delete(
            sourcePath, targetPath,
            (gDiv + SHIFT_DELETE) / SCALE_DELETE
        )
    }

    fun proposeCandidateTransferLast(): PathsDelta {
        val sourcePathIds = paths.keys.toList()
        val sourcePathScores = sourcePathIds
            .map { soft(paths.getValue(it).excess) }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourcePathScores)
        val sourceIdx = sourcePathSelector.sample()
        val sourcePath = paths.getValue(sourcePathIds[sourceIdx])

        val targetPathIds = paths.keys
            .filter { it != sourcePath.id }
        val targetPathScores = targetPathIds
            .map { soft(-paths.getValue(it).excess) / sourcePath.differenceWith(it) }
        val targetPathSelector = Randseed.INSTANCE.scoreIntegerEnum(targetPathScores)
        val targetIdx = targetPathSelector.sample()
        val targetPath = paths.getValue(targetPathIds[targetIdx])

        val w0 = targetPath.edges.mapIndexed { i, edge ->
            if (sourcePath.edges[i] != edge) {
                sourcePath.edgeExcess(i) - edgeError(edge)
            } else 0.0
        }.sum() / 2.0 / sourcePath.edges.size
        var mean = (sourcePath.weight + w0 * 3.0) / 5.0
        if (mean < 0.0) {
            mean = alpha
        }
        val weightSelector = Randseed.INSTANCE.normal(mean, max(mean, sourcePath.weight - mean - alpha))
        val weight = weightSelector.sample()
        val correctedWeight = when {
            weight < 0 -> 0.0
            weight > sourcePath.weight - alpha -> sourcePath.weight - alpha
            else -> weight
        }

        val gForward = sourcePathSelector.logProbability(sourceIdx) +
                targetPathSelector.logProbability(targetIdx) +
                weightSelector.logDensity(weight) - weightSelector.logDensity(mean)
        val w0Back = sourcePath.edges.mapIndexed { i, edge ->
            if (targetPath.edges[i] != edge) {
                targetPath.edgeExcess(i) + 2 * correctedWeight - edgeError(edge)
            } else 0.0
        }.sum() / 2.0 / sourcePath.edges.size
        var meanBack = (targetPath.weight + correctedWeight + w0Back * 3.0) / 5.0
        if (meanBack < 0.0) {
            meanBack = alpha
        }
        val weightSelectorBack = Randseed.INSTANCE.normal(meanBack, max(meanBack, sourcePath.weight - meanBack - alpha))
        val targetPathIdsBack = paths.keys
            .filter { it != targetPath.id }
        val targetPathScoresBack = targetPathIdsBack
            .map { soft(-paths.getValue(it).excess) / targetPath.differenceWith(it) }
        val targetPathSelectorBack = Randseed.INSTANCE.scoreIntegerEnum(targetPathScoresBack)
        val gBackward = sourcePathSelector.logProbability(sourcePathIds.indexOf(targetPath.id)) +
                targetPathSelectorBack.logProbability(targetPathIdsBack.indexOf(sourcePath.id)) +
                weightSelectorBack.logDensity(correctedWeight) - weightSelectorBack.logDensity(meanBack)

        val gDiv = gBackward - gForward
        if (VERBOSE)
            log.write("Transfer: $gDiv\n")
        return Transfer(
            sourcePath, targetPath, correctedWeight,
            gDiv
        )
    }

    fun proposeCandidateCommon(): PathsDelta {
        var gForward = 0.0
        var gBackward = 0.0

        val sourcePathIds = paths.keys.toList()
        val sourcePathScores = sourcePathIds.map { pid ->
            exp(paths.getValue(pid).excess)
        }
        val sourcePathSelector = Randseed.INSTANCE.scoreIntegerEnum(sourcePathScores)
        val sourceIdx = sourcePathSelector.sample()
        gForward += sourcePathSelector.logProbability(sourceIdx)
        val sourcePath = paths.getValue(sourcePathIds[sourceIdx])

        sourcePath.collectHandles()
        val handleBuckets = sourcePath.allHandles.groupBy { it.size }
        val handleLengths = handleBuckets.keys.toList()
        val handleLengthScores = handleLengths.map { len ->
            handleBuckets.getValue(len).size.toDouble() / len
        }
        val handleLengthSelector = Randseed.INSTANCE.scoreIntegerEnum(handleLengthScores)
        val handleLengthIdx = handleLengthSelector.sample()
        gForward += handleLengthSelector.logProbability(handleLengthIdx)
        val handleLength = handleLengths[handleLengthIdx]

        val handles = handleBuckets.getValue(handleLength)
        gForward -= ln(handles.size)
        val handle = handles[Randseed.INSTANCE.uniformInteger(0, handles.size - 1).sample()]

        val handleStart = entity[handle[0]].position
        val targetPathway = sourcePath.vertices.take(handleStart) +
                handle +
                sourcePath.vertices.drop(handleStart + handle.size)

        var delta = Randseed.INSTANCE.uniformReal(0.0, sourcePath.weight).sample()
        val targetPathId = checkIfPresent(targetPathway)
        gBackward = gForward // todo fix
        if (targetPathId == null) {
            if (delta > sourcePath.weight - alpha) {
                return Move(
                    sourcePath, targetPathway,
                    gBackward - gForward
                )
            } else {
                if (delta < alpha) {
                    delta = alpha
                }
                return New(
                    sourcePath, targetPathway, delta,
                    gBackward - gForward
                )
            }
        } else {
            return Fail()
            if (delta > sourcePath.weight - alpha) {
                return Delete(
                    sourcePath, paths.getValue(targetPathId),
                    gBackward - gForward
                )
            } else {
                return Transfer(
                    sourcePath, paths.getValue(targetPathId), delta,
                    gBackward - gForward
                )
            }
        }
    }

//    fun proposeCandidateNew(pNew: Double, pDel: Double): PathsDelta {
//        @Deprecated("messes all densities")
//        fun selectTargetPathway(): Pathway {
//            val source = entity.sources[sourceSelector.sample()]
//            val pathway = mutableListOf(source)
//
//            while (!entity.isRightmost(pathway.last())) {
//                val v = pathway.last()
//                val next = entity[v].rightSelector.sample()
//                pathway.add(entity.edges[v][next].target)
//            }
//
//            return pathway
//        }
//
//        @Deprecated("doesnt rely on handles")
//        fun selectClosePathway(path: Path): Pair<Pathway, Double> {
//            val right = Randseed.INSTANCE.randint() % 2
//            val next = MutableList(entity.size) { -1 }
//            val score = MutableList(entity.size) { 0.0 }
//            val prob = MutableList(entity.size) { 0.0 }
//            var cur = -1
//            val pathway = mutableListOf<Int>()
//            var g = 0.0
//            val common: Int
//
//            if (right % 2 == 0) {
//                val sep = path.rightSeparationPoints.entries.toList()
//                val point = sep[Randseed.INSTANCE.uniformInteger(0, sep.size - 1).sample()]
//                g -= ln(sep.size)
//                pathway.addAll(path.vertices.take(point.value))
//                common = point.value
//
//                fun walk(v: Int, depth: Int) {
//                    if (next[v] > -1) return
//                    for (e in entity.edges[v]) {
//                        walk(e.target, depth - 1)
//                    }
//                    if (entity.edges[v].isEmpty()) return
//
//                    val opts = entity.edges[v].map { exp(score[it.target] / depth) }
//                    val select = Randseed.INSTANCE.scoreIntegerEnum(opts)
//                    val idx = select.sample()
//                    val nedge = entity.edges[v][idx]
//                    next[v] = nedge.target
//                    score[v] = score[next[v]] - edgeError(nedge.locator)
//                    prob[v] = select.logProbability(idx)
//                }
//                walk(point.key, path.edges.size - point.value)
//                cur = point.key
//                while (cur != -1) {
//                    pathway.add(cur)
//                    g += prob[cur]
//                    cur = next[cur]
//                }
//            } else {
//                val sep = path.leftSeparationPoints.entries.toList()
//                val point = sep[Randseed.INSTANCE.uniformInteger(0, sep.size - 1).sample()]
//                g -= ln(sep.size)
//                common = path.edges.size - point.value
//
//                fun walk(v: Int, depth: Int) {
//                    if (next[v] > -1) return
//                    for (e in entity.reversedEdges[v]) {
//                        walk(e.target, depth - 1)
//                    }
//                    if (entity.reversedEdges[v].isEmpty()) return
//
//                    val opts = entity.reversedEdges[v].map { exp(score[it.target] / depth) }
//                    val select = Randseed.INSTANCE.scoreIntegerEnum(opts)
//                    val idx = select.sample()
//                    val nedge = entity.reversedEdges[v][idx]
//                    next[v] = nedge.target
//                    score[v] = score[next[v]] - edgeError(nedge.locator.reversed)
//                    prob[v] = select.logProbability(idx)
//                }
//                walk(point.key, point.value)
//                cur = point.key
//                while (cur != -1) {
//                    pathway.add(cur)
//                    g += prob[cur]
//                    cur = next[cur]
//                }
//                pathway.reverse()
//                pathway.addAll(path.vertices.drop(point.value + 1))
//            }
//            return pathway to g / (path.edges.size - common) + ln(0.5)
//        }
//
//        var gForward = 0.0
//        val pathIds = paths.keys.filter {
//            paths.getValue(it).weight >= 2 * alpha
//        }.toMutableList()
//        if (pathIds.isEmpty()) {
//            return Fail()
//        }
//        val dropPathIds = Randseed.INSTANCE.randint() % 2 == 0
//        if (dropPathIds) {
//            Randseed.INSTANCE.shuffle(pathIds)
//            val ok = 5000000 / entity.originLength
//            if (ok < pathIds.size) {
//                pathIds.dropLast(pathIds.size - ok)
//            }
//        }
//
//        val ok = 5000000 / pathIds.size
//        val weights = pathIds.map { pid ->
//            val path = paths.getValue(pid)
//            val delta = if (dropPathIds && ok < path.edges.size) {
//                var accum = 0.0
//                val dist = Randseed.INSTANCE.uniformInteger(0, path.edges.size - 1)
//                repeat(ok) {
//                    val idx = dist.sample()
//                    accum -= edgeError(path.edges[idx])
//                }
//                accum / ok
//            } else {
//                -path.edges.sumByDouble { e -> edgeError(e) } / path.edges.size
//            }
//            if (delta > 0) min(delta, path.weight) else delta
//        }
//        val pathSelector = Randseed.INSTANCE.scoreIntegerEnum(weights.map(::exp))
//        val idx = pathSelector.sample()
//        gForward += pathSelector.logProbability(idx)
//        val sourcePath = paths.getValue(pathIds[idx])
//
////        var targetPathway = selectTargetPathway()
//        var (targetPathway, g) = selectClosePathway(sourcePath)
//        var retries = 5
//        while (checkIfPresent(targetPathway) != null) {
////            targetPathway = selectTargetPathway()
//            val res = selectClosePathway(sourcePath)
//            targetPathway = res.first
//            g = res.second
//            retries--
//            if (retries == 0) {
//                return Fail()
//            }
//        }
//        gForward += g
//
////        val delta = Randseed.INSTANCE.uniformReal(alpha, sourcePath.weight - alpha).sample()
//        val wS = selectNormalWeight(0.0, 1.0)
//        var delta = wS + weights[idx]
//        var gW = normalWeightSelector.logDensity(wS)
//        if (delta < alpha) {
//            gW = ln(normalWeightSelector.cumulativeProbability(alpha - weights[idx]))
//            delta = alpha
//        } else if (delta > sourcePath.weight - alpha) {
//            gW = ln(1.0 - normalWeightSelector.cumulativeProbability(sourcePath.weight - alpha - weights[idx]))
//            delta = sourcePath.weight - alpha
//        }
//        gForward += gW
//
//        var gBackward = 0.0
//        gBackward = gForward // todo there's still the same problem as before
//        val gTotal = gBackward - gForward + ln(pDel) - ln(pNew)
//
//        return New(
//            sourcePath, targetPathway, delta,
//            selectNormalWeight(gTotal / 2, abs(gTotal))
//        )
//    }

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

        return Delete(
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
        val lam = distributionConfig.penaltyLambda
        fun pNew(pathCount: Int): Double {
            return if (pathCount > 1) 1.0 / max(6.0, lam + pathCount) else 0.5
        }
        fun pMove(pathCount: Int): Double {
            return if (pathCount > 1) 1.0 / max(6.0, lam + lam / pathCount.toDouble()) else 0.5
        }
        fun pDel(pathCount: Int): Double {
            return if (pathCount > 1) pathCount.toDouble() / lam / (lam + pathCount) else 0.0
        }

        val pNewCurrent = pNew(paths.size)
        val pMoveCurrent = pMove(paths.size)
        val pDelCurrent = pDel(paths.size)

        val option = proposer.sample()
        val result = when {
            option <= pNewCurrent ->
//                proposeCandidateNew(pNewCurrent, pDel(paths.size + 1))
//                proposeCandidateCommon()
                proposeCandidateNewLast(pNewCurrent, pDel(paths.size + 1))
            option <= pNewCurrent + pMoveCurrent ->
                proposeCandidateMoveLast()
            option <= pNewCurrent + pMoveCurrent + pDelCurrent ->
//                proposeCandidateDel(pDelCurrent, pNew(paths.size - 1))
                proposeCandidateDeleteLast(pDelCurrent, pNew(paths.size - 1))
            else ->
//                proposeCandidateTransfer()
                proposeCandidateTransferLast()
        }
        return if (result is Fail) proposeCandidate() else result
    }

    override fun extractResult(): List<Pair<String, Double>> {
        return paths.map { it.value.collectDNA() to it.value.weight }
    }

}