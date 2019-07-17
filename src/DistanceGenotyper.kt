package edu.nus

import org.tc33.jheatchart.HeatChart
import java.io.File
import java.lang.Math.*
import kotlin.random.Random
import kotlin.collections.ArrayList

data class SegmentEnd(val pos: Int, var extPos: Int, val profile: Matrix, val threshold: Double)
class DistanceGenotyper(private val bamFiles:List<String>, private val outDir: File, private val heatmap:Boolean=false){
    private var hapLen = 100000L
    private fun diffSum(a: Array<DoubleArray>, b: Array<DoubleArray>) = a.indices.sumByDouble { i ->
            a.indices.sumByDouble { j -> Math.abs(a[i][j]-b[i][j])}
        }
    private fun getGenotypeVectorOnChrom(bamScanner: BamFilesParallelScanner): ArrayList<Pair<Int, GenotypeVector>> {
        val contig = bamScanner.getChrom()
        val contiglen = bamScanner.header.getSequence(contig).sequenceLength
        val sites = ArrayList<Pair<Int, GenotypeVector>> (contiglen / 1000)
        val bwseq = StringBuilder()
        do {
            val gv = bamScanner.get()
            if (gv.isReliable) {
                bwseq.append("$contig\t${bamScanner.pos-1}\t${bamScanner.pos}\n")
                sites.add(Pair(bamScanner.pos, gv))
            }
        } while (bamScanner.nextPosition())
//        File(outDir,"SNP_$contig.bed").writeText(bwseq.toString())
        return sites
    }
    private fun scanLeftToRight(contig: String, sites: ArrayList<Pair<Int, GenotypeVector>>): ArrayList<Matrix?> {
        val dist = Matrix(bamFiles.size)
        val count = Matrix(bamFiles.size)
        var lastsite = 0
        val snapshot = ArrayList<Matrix?>(sites.size)
        for ((pos, gv) in sites) {
            while (sites[lastsite].first+hapLen<pos) {
                dist += (Matrix(sites[lastsite].second.toDMatrix())*(-1.0))
                count += (Matrix(sites[lastsite].second.toCMatrix())*(-1.0))
                lastsite++
            }
            dist += Matrix(gv.toDMatrix())
            count += Matrix(gv.toCMatrix())
            snapshot.add(if (count.containsZero()) null else dist/count)
        }
        val unknown = bamFiles.indices.filter { count.data[it][it]==0.0 }
        if (unknown.isNotEmpty()) {
            println("\n[Info] No reliable genotype on $contig from: ")
            unknown.forEach { println(bamFiles[it]) }
        }
        return snapshot
    }
    private fun scanSampling(snapshot: ArrayList<Matrix?>): Double {
        val random = Random(0)
        fun Random.shuffle(n: Int): Array<Int> {
            val order = Array(n){it}
            order.indices.reversed().forEach { i ->
                val j = this.nextInt(i+1)
                val temp = order[i]
                order[i] = order[j]
                order[j] = temp
            }
            return order
        }
        val repeatPerPair = Math.ceil(10000.0/snapshot.count { it!=null }).toInt()
        val diffsum = mutableListOf<Double>()
        for (mat in snapshot.filterNotNull()) {
            for (i in 1..repeatPerPair) {
                val order = random.shuffle(bamFiles.size)
                val shuffle = Array(order.size) { mat.data[order[it]] }
                diffsum.add(diffSum(shuffle, mat.data))
            }
        }
        return diffsum.min()!!
    }
    private fun scanGV(sites: ArrayList<Pair<Int, GenotypeVector>>, snapshot: ArrayList<Matrix?>, contiglen: Int): ArrayList<Pair<SegmentEnd, SegmentEnd>?> {
        val dist = Matrix(bamFiles.size)
        val count = Matrix(bamFiles.size)
        var lastsite = sites.size - 1
        var lastSegEnd = sites.size - 1
        val segment = ArrayList<Pair<SegmentEnd,SegmentEnd>?>()
        val threshold = scanSampling(snapshot)
        for (i in sites.indices.reversed()) {
            val (pos, gv) = sites[i]
            while (pos+hapLen<sites[lastsite].first) {
                dist += (Matrix(sites[lastsite].second.toDMatrix())*(-1.0))
                count += (Matrix(sites[lastsite].second.toCMatrix())*(-1.0))
                lastsite--
            }
            dist += Matrix(gv.toDMatrix())
            count += Matrix(gv.toCMatrix())
            if (!count.containsZero() && 0 < i && snapshot[i - 1] != null && hapLen<sites[i-1].first && pos+hapLen<contiglen) {
                val avg = dist/count
                val diff = diffSum(snapshot[i - 1]!!.data, avg.data)
                if (diff > threshold ) {
                    segment.add(Pair(SegmentEnd(sites[i].first, sites[i - 1].first, avg, threshold), SegmentEnd(
                            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
                                if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first,snapshot[lastSegEnd]!!, threshold)
                    ))
                    lastSegEnd = i - 1
                }
            }
        }
        val seg = Pair(SegmentEnd(1, 1, dist/count, threshold), SegmentEnd(
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first,
            snapshot[lastSegEnd]!!, threshold))
        segment.add(seg)
        return segment
    }
    private fun getSegmentOfChrom(bamScanner: BamFilesParallelScanner): List<Triple<String, Matrix, Double>> {
        val contig = bamScanner.getChrom()
        val contiglen = bamScanner.header.getSequence(contig).sequenceLength
        val sites = getGenotypeVectorOnChrom(bamScanner)
        if (sites.size < 1) return  mutableListOf()
        val snapshot = scanLeftToRight(contig, sites)
        if (snapshot.last()==null) return mutableListOf()
        val segment = scanGV(sites, snapshot, contiglen)
        // merge non-solid segments (100% flexible end)
        segment.reverse()
        if (segment.size<=10)
            segment.forEach { println("[${it!!.first.extPos},${it!!.first.pos}][${it!!.second.pos},${it!!.second.extPos}]") }
        for (i in 1 until segment.size - 1) {
            val (j, k, l) = Triple(segment[i - 1], segment[i], segment[i + 1])
            if (k!!.first.pos==k.second.pos) {
                l!!.first.extPos = k.first.extPos
                j!!.second.extPos = k.second.extPos
                segment[i] = j
                segment[i - 1] = null
            }
        }
        println("\n$contig[${contiglen} bp / ${sites.size} SNP] -> ${segment.count{it!=null}} segments (${segment.size} before adjacent merge)")
        val named_segment = segment.filterNotNull().flatMap{x -> listOf(
            Triple("+$contig(${x.first.extPos},${x.first.pos})(${x.second.pos},${x.second.extPos})", x.first.profile, x.first.threshold),
            Triple("-$contig(${x.first.extPos},${x.first.pos})(${x.second.pos},${x.second.extPos})", x.second.profile, x.second.threshold))}
        if (heatmap) {
            val heatmapDir = File(outDir, "heatmap")
            if (!heatmapDir.exists()) heatmapDir.mkdirs()
            for ((name, mat, _) in named_segment) {
                HeatChart(mat.data).saveToFile(File(heatmapDir,"${name.substring(1)}_${if (name[0]=='+') "L" else "R"}.png"))
                val temp = Matrix(bamFiles.size)
                for (i in temp.data.indices)
                    for (j in temp.data.indices)
                        temp.data[i][j] = Util.PearsonCorrelationSimilarity(mat.data[i].toList(), mat.data[j].toList())
                HeatChart(temp.data).saveToFile(File(heatmapDir,"${name.substring(1)}_${if (name[0]=='+') "L" else "R"}pearson.png"))
//                mat.data.forEach { row -> println(row.joinToString(separator = ","){"%.2f".format(it)}) }
            }
        }
//        exitProcess(1)
        return named_segment
    }
    private fun getPairwiseSegment(bamScanner: BamFilesParallelScanner): ArrayList<Triple<String, String, Double>> {
        val contigsGV = mutableListOf<Triple<String, Matrix, Double>>()
        while (bamScanner.hasNext()) {
            contigsGV.addAll(getSegmentOfChrom(bamScanner))
        }
        val candidate = ArrayList<Triple<String, String, Double>>()
        for ((name, mat, t) in contigsGV) {
            val contig = name.substring(1)
            contigsGV.filter{ it.first.substring(1) != contig}
                .map{Triple(it.first, diffSum(mat.data, it.second.data), max(t, it.third))}
                .filter{ it.second < it.third }
                .forEach{candidate.add(Triple(name, it.first, it.second))}
        }
        return candidate
    }
    fun run(): File{
        val reportFile = File(outDir, "report")
        if (reportFile.exists()) reportFile.delete()
        val reportWriter = reportFile.bufferedWriter()
        val coverage = DepthDetector(bamFiles).getCoverage
        println(coverage.map { it.first }.joinToString("|","\nLower:"))
        println(coverage.map { it.second.toInt() }.joinToString("|","\nDepth:"))
        println(coverage.map { it.third }.joinToString("|","\nUpper:"))
        val bamScanner = BamFilesParallelScanner(bamFiles, coverage)
        hapLen = bamScanner.header.referenceLength/1000
        val finalLink = getPairwiseSegment(bamScanner)
        finalLink.filter { it.first < it.second }.forEach { (k, v, ds) -> reportWriter.write("$k\t$v\t$ds\n") }
        reportWriter.close()
        return reportFile
    }
}