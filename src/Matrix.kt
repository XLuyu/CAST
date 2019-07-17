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
    val known = lowerbound < sum

    fun distance(other: Genotype): Double {
        if (!this.isReliable || !other.isReliable) throw Exception("[Error] attempt to compute unreliable genotypes' distance!!")
        if (!known || !other.known) return Double.NaN
        return acgt.indices.sumByDouble { Math.abs(acgt[it] - other.acgt[it]) } / 2
    }
    operator fun times(other:Genotype) = (0..3).map {  acgt[it]*other.acgt[it] }.sum()
}

class GenotypeVector(val vector: List<Genotype>, depth: List<Double>) {
    operator fun get(i: Int) = vector[i]
    fun toSqrMatrix(): Array<DoubleArray> = this.vector.any { !it.known }.let { z ->
        vector.map { i ->
            vector.map { j -> if (!z) i * j else 0.0 }.toDoubleArray()
        }.toTypedArray()
    }
    fun toCMatrix(): Array<DoubleArray> =
        vector.map { i ->
            vector.map { j ->
                if ( i.distance(j).isNaN()) 0.0 else 1.0
            }.toDoubleArray()
        }.toTypedArray()
    fun toDMatrix(): Array<DoubleArray> =
        vector.map { i ->
            vector.map { j ->
                val d = i.distance(j)
                if (d < 0.1 || d.isNaN()) 0.0 else d
            }.toDoubleArray()
        }.toTypedArray()
    fun toDistMatrix(): Array<DoubleArray> =
        vector.map { i ->
            vector.map { j ->
                val d = i.distance(j)
                if (d < 0.1) 0.0 else d
            }.toDoubleArray()
        }.toTypedArray()

    val isReliable = vector.all { it.isReliable }
            && vector.count { it.sum == 0 } <= vector.size / 2
            && isHeterogeneous(depth)
    private fun isHeterogeneousPrecheck(): Boolean {
        val nonZero = vector.filter { it.known }
        return (0..4).sumByDouble { i->
            val column = nonZero.map { it.acgt[i] }
            if (column.isNotEmpty()) column.max()!!-column.min()!! else 0.0
        } >= 0.4
    }

    private fun isHeterogeneous(depth: List<Double>) = isHeterogeneousPrecheck() && toDistMatrix().let { m ->
        val fold = vector.indices.map { vector[it].sum / depth[it] }
        m.none { x -> // to check if heterogeneity is not from mutation
            val sx = x.filter { !it.isNaN() }.sorted()
            sx.isNotEmpty() && sx.last() - sx[1] < 0.2
        }  &&
        m.map { Math.abs(Util.PearsonCorrelationSimilarity(it.toList(), fold))}.average() < 0.8
    }
}

class Matrix(var data:Array<DoubleArray>){
    constructor(size: Int) : this(Array<DoubleArray>(size) { DoubleArray(size){0.0} })
    fun diagonalZero() = data.indices.any {data[it][it]==0.0}
    fun containsZero() = data.any { row -> row.any { it.isNaN() }}
    fun toCosineMatrix() =
        data.indices.map { i ->
            data.indices.map { j -> if (data[i][i] == 0.0 || data[j][j] == 0.0) 0.0 else data[i][j] / Math.sqrt(data[i][i] * data[j][j]) }.toDoubleArray()
        }.toTypedArray()
    fun cosineClustering(): List<Int> {
        val pairs = ArrayList<Triple<Int, Int, Double>>()
        for (i in data.indices) if (data[i][i] != 0.0) {
            for (j in data.indices) if (data[j][j] != 0.0) {
                val s = data[i][j] / Math.sqrt(data[i][i] * data[j][j])
                pairs.add(Triple(i, j, s))
            }
        }
        pairs.sortByDescending { it.third }
        val UFset = UnionFind(data.size)
        for ((a, b, s) in pairs) {
            if (s<0.9) break
            UFset.union(a, b)
        }
        return data.indices.map { UFset.find(it) }
    }
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
    operator fun div(other:Matrix): Matrix  {
        val rt = this.copy()
        for (i in rt.data.indices)
            for (j in rt.data[i].indices)
                rt.data[i][j] /= other.data[i][j]
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
}
class UnionFind(size: Int): HashMap<Int, Int>(size){
    var setNum = size
    fun find(element: Int): Int{
        val father = this[element]
        if (father==null) this[element] = element
        if (father==null || father==element) return element
        val root = find(father)
        this[element] = root
        return root
    }
    fun union(A: Int, B: Int): Boolean {
        val a = find(A)
        val b = find(B)
        if (a<b) this[b] = a
        if (b<a) this[a] = b
        if (a!=b) setNum--
        return a!=b
    }
}