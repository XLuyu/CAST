package edu.nus

object Util {
    val complementary = mapOf('A' to 'T', 'C' to 'G', 'G' to 'C', 'T' to 'A')
    fun PearsonCorrelationSimilarity(rawA: List<Double>, rawB: List<Double>): Double {
        val AB = rawA.indices.filterNot { rawA[it].isNaN() || rawB[it].isNaN() }
        if (AB.size < 2) return 0.0 //throw Exception("Attempt to compute Pearson between two NaN-rich vectors!")
        val A = rawA.slice(AB)
        val B = rawB.slice(AB)
        val avgA = A.average()
        val avgB = B.average()
        val biasA = A.map { it - avgA }
        val biasB = B.map { it - avgB }
        val numerator = A.indices.sumByDouble { biasA[it] * biasB[it] }
        val dominator = Math.sqrt(biasA.sumByDouble { it * it } * biasB.sumByDouble { it * it })
        return if (dominator == 0.0) 0.0 else numerator / dominator
    }
}

class Genotype(acgt_: Array<Int>, badCount: Int = 0, lowerbound: Double, upperBound:Double) {
    var sum = acgt_.sum()
    var acgt = if (sum != 0) acgt_.map { it / sum.toDouble() } else acgt_.map { 0.0 }
    val isReliable = badCount <= 0.1 * sum && lowerbound < sum && sum <= upperBound
    val known = 0 < sum
    fun distance(other: Genotype): Double {
        if (!this.isReliable || !other.isReliable) throw Exception("[Error] attempt to compute unreliable genotypes' distance!!")
        if (!known || !other.known) return Double.NaN
        return acgt.indices.sumByDouble { Math.abs(acgt[it] - other.acgt[it]) } / 2
    }
}

class GenotypeVector(vector: List<Genotype>, depth: List<Double>): ArrayList<Genotype>(vector) {
    val isReliable = all { it.isReliable } && count { it.sum == 0 } <= size / 2  && isHeterogeneous(depth)
    private inline fun selfOuterProduct(calcValue:(Genotype, Genotype) -> Double) =
        map { genotypeA ->
            map { genotypeB ->
                calcValue(genotypeA, genotypeB)
            }.toTypedArray()
        }.toTypedArray()
    fun toCMatrix() = selfOuterProduct { A, B -> if ( A.distance(B).isNaN()) 0.0 else 1.0}
    fun toDMatrix() = selfOuterProduct { A, B -> A.distance(B).let { if (it < 0.1 || it.isNaN()) 0.0 else it }}
    fun toDistMatrix() = selfOuterProduct { A, B -> A.distance(B).let { if (it < 0.1) 0.0 else it }}
    private fun isHeterogeneousPrecheck(): Boolean {
        val nonZero = filter { it.known }
        return (0..4).sumByDouble { i->
            val column = nonZero.map { it.acgt[i] }
            if (column.isNotEmpty()) column.max()!!-column.min()!! else 0.0
        } >= 0.4
    }
    private fun isHeterogeneous(depth: List<Double>) = isHeterogeneousPrecheck() && toDistMatrix().let { m ->
        val fold = indices.map { this[it].sum / depth[it] }
        m.none { x -> // to check if heterogeneity is not from mutation
            val sx = x.filter { !it.isNaN() }.sorted()
            sx.isNotEmpty() && sx.last() - sx[1] < 0.2
        }  &&
        m.map { Math.abs(Util.PearsonCorrelationSimilarity(it.toList(), fold))}.average() < 0.8
    }
}

class Matrix(var data:Array<Array<Double>>){
    constructor(size: Int) : this(Array(size) { Array(size){0.0} }) // new zero matrix
    fun containsZero() = data.any { row -> row.any { it==0.0 }}
    fun containsNaN() = data.any { row -> row.any { it.isNaN() }}
    private inline fun mapEachCell(f: (Int, Int) -> Double): Matrix {
        val new = Matrix(data.size)
        for (i in data.indices)
            for (j in data[i].indices)
                new.data[i][j] = f(i, j)
        return new
    }
    private inline fun updateEachCell(f: (Int, Int) -> Double) {
        for (i in data.indices)
            for (j in data[i].indices)
                data[i][j] = f(i, j)
    }
    operator fun times(factor: Double) = mapEachCell { i, j -> data[i][j] * factor }
    operator fun div(other: Matrix) = mapEachCell { i, j -> data[i][j] / other.data[i][j]}
    operator fun plus(other: Matrix) = mapEachCell { i, j -> data[i][j] + other.data[i][j]}
    operator fun plusAssign(other: Matrix) = updateEachCell { i, j -> data[i][j] + other.data[i][j]}
    fun updateBy(other:Matrix) = updateEachCell { i, j -> if (other.data[i][j].isNaN()) data[i][j] else other.data[i][j] }
    fun fill(v:Double) = updateEachCell { _, _ -> v }
}