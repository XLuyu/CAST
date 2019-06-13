package edu.nus

import kotlin.math.roundToInt

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

class Genotype(acgt_: List<Double>, badCount: Int = 0) { //TODO: distinguish reliable/known
    var sum = acgt_.sum().roundToInt()
    var acgt = if (sum != 0) acgt_.map { it / sum } else acgt_
    val isReliable = badCount <= 0.1 * sum && acgt[4] < 0.1
    fun checksum(lowerbound: Double) {
        if (sum <= lowerbound) sum = 0
    }

    fun distance(other: Genotype): Double {
        if (!this.isReliable || !other.isReliable) throw Exception("[Error] attempt to compute unreliable genotypes' distance!!")
        if (this.sum == 0 || other.sum == 0) return Double.NaN
        return acgt.indices.sumByDouble { Math.abs(acgt[it] - other.acgt[it]) } / 2
    }
}

class GenotypeVector(val vector: List<Genotype>) {
    fun toDistMatrix(): Array<DoubleArray> =
        vector.map { i ->
            vector.map { j ->
                val d = i.distance(j)
                if (d < 0.1) 0.0 else d
            }.toDoubleArray()
        }.toTypedArray()

    val isReliable = vector.all { it.isReliable } && vector.count { it.sum == 0 } <= vector.size / 2
    private fun isHeterogeneousPrecheck(): Boolean {
        val nonZero = vector.filter { it.sum > 0 }
        return (0..4).sumByDouble { i->
            val column = nonZero.map { it.acgt[i] }
            if (column.isNotEmpty()) column.max()!!-column.min()!! else 0.0
        } >= 0.4
    }

    fun isHeterogeneous() = isHeterogeneousPrecheck() && toDistMatrix().none { x ->
            val sx = x.filter { !it.isNaN() }.sorted()
            sx.isNotEmpty() && sx.last() - sx[1] < 0.2
        } // to check if heterogeneity is not from mutation

    fun consistentWithDepth(depth: List<Double>): Boolean {
        val fold = vector.indices.map { vector[it].sum / depth[it] }
        try {
            return fold.all { it <= 1.8 }
                    && toDistMatrix().map { Math.abs(Util.PearsonCorrelationSimilarity(it.toList(), fold))}.average() < 0.8
        } catch (e:Exception){
            println(fold.joinToString())
            vector.forEach { println(it.sum) }
            throw e
        }
    }
}

class Matrix(var data:Array<DoubleArray>){
    constructor(size: Int) : this(Array<DoubleArray>(size) { DoubleArray(size){0.0} })
    fun timesUpdate(factor:Double) {
        for (i in data.indices)
            for (j in data[i].indices)
                data[i][j] *= factor
    }
    operator fun times(factor:Double): Matrix  {
        val rt = this.copy()
        for (i in rt.data.indices)
            for (j in rt.data[i].indices)
                rt.data[i][j] *= factor
        return rt
    }
    operator fun plusAssign(other:Matrix)  {
        for (i in data.indices)
            for (j in data[i].indices)
                data[i][j] += other.data[i][j]
    }
    fun copy(): Matrix  {
        val other = Matrix(data.size)
        for ( i in data.indices)
            for (j in data[i].indices)
            other.data[i][j] = data[i][j]
        return other
    }
    fun normalized(): Matrix  {
        val max = data.map{ it.max()!! }
        if (!max.contains(0.0)) {
            for ( i in data.indices)
                for (j in data[i].indices)
                    data[i][j] /= max[i]
        }
        return this
    }
    fun updateBy(other:Matrix) {
        for (i in data.indices)
            for (j in data[i].indices) if (!other.data[i][j].isNaN())
                data[i][j] = other.data[i][j]
    }
    fun fill(v:Double) {
        for (i in data.indices)
            for (j in data[i].indices)
                data[i][j] = v
    }
    fun containsNaN() = data.any { row -> row.any { it.isNaN() }}
    fun printline(hint:String) {
        println("======$hint")
        for ( i in data){
            i.forEach{ print("$it%.2f\t")}
            println()
        }
    }
}
