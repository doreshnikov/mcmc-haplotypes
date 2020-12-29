package mcmc.modules.graph

import graph.AlignedDeBruijnGraph
import graph.components.Edge
import mcmc.Model

import org.apache.commons.math3.distribution.PoissonDistribution
import kotlin.math.ln
import kotlin.math.pow

abstract class PathsOverlay(
    entity: AlignedDeBruijnGraph,
    private val distributionConfig: DistributionConfig,
    private val paths: MutableList<Path>
) : Model<AlignedDeBruijnGraph, PathsOverlay.PathsDelta>(entity) {

    private val weights = hashMapOf<Edge, Double>()
    private val coverage = hashMapOf<Edge, Double>()

    private val distribution = PoissonDistribution(distributionConfig.penaltyLambda)

    init {
        entity.indices.forEach { v ->
            entity.edges[v].forEach { u ->
                val locator = u.locator
                weights[locator] = u.weight
            }
        }
        paths.forEach { layPath(coverage, it) }
    }

    private fun logCoverageError(weight: Double, cover: Double): Double {
        return -(weight - cover).pow(2)
    }

    private fun layPath(localCoverage: HashMap<Edge, Double>, path: Path) {
        for (i in 0 until path.size - 1) {
            val locator = path.edge(i)
            check(locator in weights) { "Can't lay path on non-existent edge" }
            localCoverage.compute(locator) { _, cover -> (cover ?: 0.0) + path.weight }
        }
    }

    private fun erasePath(pathId: Int) {
        layPath(coverage, -paths[pathId])
        paths.removeAt(pathId)
    }

    private fun pathId(path: Path): Int {
        return paths.indexOfFirst { localPath ->
            localPath.size == path.size && localPath.indices.all { localPath[it] == path[it] }
        }
    }

    private fun transferWeight(fromPathId: Int, toPathId: Int, weight: Double): Boolean {
        check(weight <= paths[fromPathId].weight) { "Can't have negative weights in coverage" }
        layPath(coverage, paths[fromPathId] withWeight -weight)
        paths[fromPathId].weight -= weight
        layPath(coverage, paths[toPathId] withWeight weight)
        paths[toPathId].weight += weight

        if (paths[fromPathId].weight == 0.0) {
            paths.removeAt(fromPathId)
            return true
        }
        return false
    }

    private fun logCoverageDeltaLikelihood(coverageDelta: HashMap<Edge, Double>) =
        coverageDelta.entries.sumByDouble { (locator, delta) ->
            val weight = weights.getValue(locator)
            val cover = coverage.getOrDefault(locator, 0.0)
            // logCoverageError(weight, cover + delta) - logCoverageError(weight, cover)
            // -(w - c - d)^2 + (w - c)^2 = d (2w - 2c - d)
            delta * (2 * weight - 2 * cover - delta)
        }

    private fun logPathCountDeltaLikelihood(pathAdded: Boolean, pathRemoved: Boolean): Double =
        when ((if (pathAdded) 1 else 0) + (if (pathRemoved) -1 else 0)) {
            0 -> 0.0
            // else -> distribution.logProbability(paths.size + pathCountDelta) - distribution.logProbability(paths.size)
            // log(l^k e^-l / k!) = k logl - l - logk!
            1 -> ln(distributionConfig.penaltyLambda) - ln((paths.size + 1).toDouble())
            else -> -ln(distributionConfig.penaltyLambda) + ln(paths.size.toDouble())
        }

    abstract inner class PathsDelta : Delta() {
        inner class Transfer(
            private val fromPathId: Int,
            private val toPath: Path
        ) : PathsDelta() {
            private val pathAdded = pathId(toPath) == -1
            private val pathRemoved = paths[fromPathId].weight == toPath.weight

            override fun logLikelihoodDelta(): Double {
                val coverageDelta = hashMapOf<Edge, Double>()
                layPath(coverageDelta, paths[fromPathId] withWeight -toPath.weight)
                layPath(coverageDelta, toPath)
                return logCoverageDeltaLikelihood(coverageDelta) + logPathCountDeltaLikelihood(pathAdded, pathRemoved)
            }

            override fun accept() {
                val toPathId = if (pathAdded) {
                    paths.add(toPath withWeight 0.0)
                    paths.size - 1
                } else {
                    pathId(toPath)
                }
                transferWeight(fromPathId, toPathId, toPath.weight)
            }
        }

        inner class Interchange(
            private val pathId1: Int,
            private val pathId2: Int,
            crossingPoint: Int,
            weightTransferRatio: Double
        ) : PathsDelta() {
            private val path1 = paths[pathId1]
            private val path2 = paths[pathId2]

            private val newPath1 = Path(
                path1.subList(0, crossingPoint) + path2.subList(crossingPoint, path2.size),
                path1.weight * weightTransferRatio + path2.weight * (1 - weightTransferRatio)
            )
            private val newPath2 = Path(
                path2.subList(0, crossingPoint) + path1.subList(crossingPoint, path1.size),
                path2.weight * weightTransferRatio + path1.weight * (1 - weightTransferRatio)
            )

            init {
                check(paths[pathId1][crossingPoint] == paths[pathId2][crossingPoint]) {
                    "Paths should have same vertex on the crossing point"
                }
            }

            override fun logLikelihoodDelta(): Double {
                val coverageDelta = hashMapOf<Edge, Double>()
                layPath(coverageDelta, -path1)
                layPath(coverageDelta, -path2)
                layPath(coverageDelta, newPath1)
                layPath(coverageDelta, newPath2)
                return logCoverageDeltaLikelihood(coverageDelta)
            }

            override fun accept() {
                erasePath(pathId1)
                erasePath(pathId2)
                newPath1.let {
                    layPath(coverage, it)
                    paths.add(it)
                }
                newPath2.let {
                    layPath(coverage, it)
                    paths.add(it)
                }
            }
        }
    }

    override fun logLikelihood(): Double {
        val coverageLikelihood = weights.keys.sumByDouble { locator ->
            logCoverageError(weights.getValue(locator), coverage.getOrDefault(locator, 0.0))
        }
        return coverageLikelihood + distribution.logProbability(paths.size)
    }

    override fun proposeCandidate(): PathsDelta {
        TODO("Not yet implemented")
    }

}