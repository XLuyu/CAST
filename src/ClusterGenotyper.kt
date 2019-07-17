package edu.nus

import org.tc33.jheatchart.HeatChart
import java.io.File
import org.jgrapht.graph.DefaultWeightedEdge
import org.jgrapht.graph.SimpleWeightedGraph
import org.jgrapht.alg.matching.MaximumWeightBipartiteMatching
import kotlin.system.exitProcess

typealias ListList = Pair<List<Int>, List<Int>>
class ClusterGenotyper(private val bamFiles:List<String>, private val outDir: File, private val heatmap:Boolean=false){
    private var hapLen = 100000L
    private fun getOneChromGVFromScanner(bamScanner: BamFilesParallelScanner): List<Pair<String, ListList?>> {
        // scan this contig to get all snp sites
        val contig = bamScanner.getChrom()
        val contiglen = bamScanner.header.getSequence(contig).sequenceLength
        val sites = ArrayList<Pair<Int, GenotypeVector>> (contiglen / 1000)
        val bwseq = StringBuilder()
        do {
            val gv = bamScanner.get()
            if (gv.isReliable) {
                bwseq.append("$contig\t${bamScanner.pos-1}\t${bamScanner.pos}\t${gv.vector.indices.filter { !gv[it].known }.joinToString(separator = "_")}\n")
                sites.add(Pair(bamScanner.pos, gv))
            }
        } while (bamScanner.nextPosition())
        File(outDir,"SNP_$contig.bed").writeText(bwseq.toString())
        if (sites.size < 1) return  mutableListOf<Pair<String, ListList>>()

        val snapshot = ArrayList<List<Int>?>(sites.size)
        // left to right scan
        var lastsite = 0
        var xy = Matrix(bamFiles.size)
        val indices = 0 until bamFiles.size
        for ((pos, gv) in sites) {
            while (sites[lastsite].first+hapLen<pos) {
                xy += (Matrix(sites[lastsite].second.toSqrMatrix())*(-1.0))
                lastsite++
            }
            val delta = Matrix(gv.toSqrMatrix())
            xy += delta
            snapshot.add(if (xy.diagonalZero()) null else xy.cosineClustering())
        }
        val unknown = indices.filter { xy.data[it][it]==0.0 }
        if (unknown.isNotEmpty()) {
            println("\n[Info] No reliable genotype on $contig from: ")
            unknown.forEach { println(bamFiles[it]) }
            return ArrayList<Pair<String, ListList>>()
        }
        // right to left scan
        val heatmapDir = File(outDir, "Inheretance").bufferedWriter()
        lastsite = sites.size - 1
        xy = Matrix(bamFiles.size)
        var lastSegEnd = sites.size - 1
        val segment = ArrayList<Pair<Array<Int>, ListList>?>()
        for (i in sites.indices.reversed()) {
            val (pos, gv) = sites[i]
            while (sites[lastsite].first>pos+hapLen) {
                xy += (Matrix(sites[lastsite].second.toSqrMatrix())*(-1.0))
                lastsite--
            }
            xy += Matrix(gv.toSqrMatrix())
            if (!xy.diagonalZero() && 0 < i && snapshot[i - 1] != null) {
                val cvector = xy.cosineClustering()
                val reombcount = recomb(snapshot[i - 1]!!, cvector) // network flow
                heatmapDir.write("$contig:${sites[i-1].first}-$pos $reombcount\n")
                heatmapDir.write(cvector.joinToString(separator = " ", postfix = "\n"))
                heatmapDir.write(snapshot[i - 1]!!.joinToString(separator = " ", postfix = "\n"))
                if (reombcount > 2) { // define threshold later
                    val seg = Pair(arrayOf(sites[i - 1].first, sites[i].first,
                        if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
                        if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first),
                        Pair(cvector, snapshot[lastSegEnd]!!))
                    segment.add(seg)
                    lastSegEnd = i - 1
//                    if (pos==13186942){
                        HeatChart(xy.data).saveToFile(File(outDir,"${pos}_data.png"))
                        HeatChart(xy.toCosineMatrix()).saveToFile(File(outDir,"${pos}_cos.png"))
//                        for (i in xy.data) println(i.joinToString { "%.2f".format(it) })
//                    }
                }
            }
        }
        val seg = Pair(arrayOf(1, 1,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first),
            Pair(xy.cosineClustering(), snapshot[lastSegEnd]!!))
        segment.add(seg)
        // merge non-solid segments (100% flexible end)
        segment.reverse()
        for (i in 1 until segment.size - 1) {
            val (j, k, l) = Triple(segment[i - 1], segment[i], segment[i + 1])
            if (k!!.first[1]==k.first[2]) {
                l!!.first[0] = k.first[0]
                j!!.first[3] = k.first[3]
                segment[i] = j
                segment[i - 1] = null
            }
        }
        println("\n$contig[${sites.size} SNPs] -> ${segment.count{it!=null}} segments (${segment.size} before adjacent merge)")
        val named_segment = segment.filterNotNull().map{x ->
            Pair("$contig(${x.first[0]},${x.first[1]})(${x.first[2]},${x.first[3]})", x.second)}
//        if (heatmap) {
//            val heatmapDir = File(outDir, "Inheretance").bufferedWriter()
//            for ((name, headListAndTailList) in named_segment) {
//                heatmapDir.write(name+"\n")
//                heatmapDir.write(headListAndTailList.first.joinToString(separator = " ", postfix = "\n"))
//                heatmapDir.write(headListAndTailList.second.joinToString(separator = " ", postfix = "\n"))
//            }
//        }
        heatmapDir.close()
        exitProcess(1)
        return named_segment
    }
    private fun getAllGVFromScanner(bamScanner: BamFilesParallelScanner): MutableMap<String, ListList> {
        val GVmap = mutableMapOf<String, ListList>()
        while (bamScanner.hasNext()) {
            val CGV = getOneChromGVFromScanner(bamScanner)
            for ((contig, headMatrixAndTailMatrix) in CGV)
                GVmap[contig] = headMatrixAndTailMatrix!!
        }
        return GVmap
    }
    private fun pairwiseMutualBest(contigsGV: List<Pair<String, List<Int>>>): ArrayList<Pair<String, Pair<String, Int>>> {
        val candidate = ArrayList<Pair<String, Pair<String, Int>>>()
        for ((k, v) in contigsGV) {
            val contig = k.substring(1)
            val similarities = contigsGV.filter{ it.first.substring(1) != contig}.map{Pair(it.first, recomb(v, it.second))} //network
            similarities.filter{it.second < 3}.forEach{candidate.add(Pair(k, Pair(it.first,it.second)))} //define threshold
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
        val contigsTwoEndGV = getAllGVFromScanner(bamScanner)
        val contigsOneEndGV = contigsTwoEndGV .flatMap { (contig, v) -> listOf(Pair("+$contig", v.first), Pair("-$contig", v.second)) }
        val finalLink = pairwiseMutualBest(contigsOneEndGV)
        for ((k, v) in finalLink) if (k < v.first) {
            reportWriter.write("$k\t${v.first}\t${v.second}\n")
        }
        reportWriter.close()
        return reportFile
    }
}
fun recomb(A: List<Int>, B: List<Int>): Int{
    val g = SimpleWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge::class.java)
    val AB = (A zip B).map { "A${it.first}B${it.second}" }
    val a = A.map { "A$it" }.toSet()
    val b = B.map { "B$it" }.toSet()
    a.forEach { g.addVertex(it) }
    b.forEach { g.addVertex(it) }
    for (i in a)
        for (j in b){
            g.addEdge(i,j)
            g.setEdgeWeight(i,j, AB.count{ it==i+j }.toDouble())
        }
    return AB.size - MaximumWeightBipartiteMatching(g, a, b).matching.weight.toInt()
}