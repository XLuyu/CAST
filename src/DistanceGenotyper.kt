package edu.nus

import java.io.File
import java.lang.Math.*
import kotlin.random.Random
import kotlin.collections.ArrayList
import edu.nus.Util.shuffle

data class SegmentEnd(val pos: Int, var extPos: Int, val profile: Matrix, val threshold: Double)
class HaplotypeWindow(private val sites:ArrayList<Pair<Int, GenotypeVector>>, matrixSize:Int){
    private val dist = Matrix(matrixSize)
    private val count = Matrix(matrixSize)
    private var head = 0
    private var tail = 0
    private fun add(gv: GenotypeVector){
        dist += Matrix(gv.toDMatrix())
        count += Matrix(gv.toCMatrix())
    }
    private fun remove(gv: GenotypeVector){
        dist += (Matrix(gv.toDMatrix()) * -1.0)
        count += (Matrix(gv.toCMatrix()) * -1.0)
    }
    fun shiftTo(left: Int, right: Int): Matrix? { // left and right are inclusive
        while (tail<sites.size && sites[tail].first<=right) add(sites[tail++].second)
        while (head<sites.size && sites[head].first<left) remove(sites[head++].second)
        while (0<head && left<=sites[head-1].first) add(sites[--head].second)
        while (0<tail && right<sites[tail-1].first) remove(sites[--tail].second)
        return if (count.containsZero()) null else dist/count
    }
    fun absentee() = count.data.indices.filter { count.data[it][it]==0.0 }
}
class DistanceGenotyper(private val bamFiles:List<String>, private val outDir: File, private val heatmap:Boolean=false){
    private var hapLen = 100000
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
        val haplotype = HaplotypeWindow(sites, bamFiles.size)
        val snapshot = ArrayList<Matrix?>(sites.size)
        for ((pos, _) in sites) {
            val avg = haplotype.shiftTo(pos-hapLen, pos)
            snapshot.add(avg)
        }
        val unknown = haplotype.absentee()
        if (unknown.isNotEmpty()) {
            println("\n[Info] No reliable genotype on $contig from: ")
            unknown.forEach { println("\t"+bamFiles[it]) }
        }
        return snapshot
    }
    private fun scanSampling(snapshot: ArrayList<Matrix?>): Double {
        val random = Random(0)
        val repeatPerPair = Math.ceil(10000.0/snapshot.count { it!=null }).toInt()
        val diffsum = mutableListOf<Double>()
        for (mat in snapshot.filterNotNull()) {
            for (i in 1..repeatPerPair) {
                val order = random.shuffle(bamFiles.size)
                val shuffle = Array(order.size) { mat.data[order[it]] }
                diffsum.add(mat.absDiffSum(shuffle))
            }
        }
        return diffsum.min()!!
    }
    private fun scanGV(sites: ArrayList<Pair<Int, GenotypeVector>>, snapshot: ArrayList<Matrix?>, contiglen: Int): ArrayList<Pair<SegmentEnd, SegmentEnd>?> {
        val haplotype = HaplotypeWindow(sites, bamFiles.size)
        val threshold = scanSampling(snapshot)
        val segment = ArrayList<Pair<SegmentEnd,SegmentEnd>?>()
        var lastSegEnd = SegmentEnd(1, 1, haplotype.shiftTo(sites[0].first, sites[0].first+hapLen)!!, threshold)
        for (i in 1 until sites.size) {
            val pos = sites[i].first
            val avg = haplotype.shiftTo(pos, pos+hapLen)
            if (avg==null || snapshot[i-1]==null || sites[i-1].first<hapLen || contiglen-hapLen+1<pos) continue
            val diff = snapshot[i-1]!!.absDiffSum(avg)
            if (diff > threshold ) {
                segment.add(Pair(lastSegEnd, SegmentEnd(sites[i-1].first, sites[i].first, snapshot[i-1]!!, threshold)))
                lastSegEnd = SegmentEnd(sites[i].first, sites[i - 1].first, avg, threshold)
            }
        }
        segment.add(Pair(lastSegEnd, SegmentEnd(contiglen, contiglen, snapshot.last()!!, threshold)))
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
        println()
        if (segment.size<=10)
            segment.forEach { println("\t[${it!!.first.extPos},${it.first.pos}][${it.second.pos},${it.second.extPos}]") }
        for (i in 1 until segment.size - 1) {
            val (j, k, l) = Triple(segment[i - 1], segment[i], segment[i + 1])
            if (k!!.first.pos==k.second.pos) {
                l!!.first.extPos = k.first.extPos
                j!!.second.extPos = k.second.extPos
                segment[i] = j
                segment[i - 1] = null
            }
        }
        println("$contig[$contiglen bp / ${sites.size} SNP] -> ${segment.count{it!=null}} segments (${segment.size} before adjacent merge)")
        val named_segment = segment.filterNotNull().flatMap{x -> listOf(
            Triple("+$contig(${x.first.extPos},${x.first.pos})(${x.second.pos},${x.second.extPos})", x.first.profile, x.first.threshold),
            Triple("-$contig(${x.first.extPos},${x.first.pos})(${x.second.pos},${x.second.extPos})", x.second.profile, x.second.threshold))}
        if (heatmap) {
            val heatmapDir = File(outDir, "heatmap")
            if (!heatmapDir.exists()) heatmapDir.mkdirs()
            for ((name, mat, _) in named_segment) {
                mat.saveHeatMap(File(heatmapDir,"${name.substring(1)}_${if (name[0]=='+') "L" else "R"}.png"))
                val temp = Matrix(bamFiles.size)
                for (i in temp.data.indices)
                    for (j in temp.data.indices)
                        temp.data[i][j] = Util.PearsonCorrelationSimilarity(mat.data[i].toList(), mat.data[j].toList())
                temp.saveHeatMap(File(heatmapDir,"${name.substring(1)}_${if (name[0]=='+') "L" else "R"}pearson.png"))
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
                .map{Triple(it.first, mat.absDiffSum(it.second), max(t, it.third))}
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
        hapLen = (bamScanner.header.referenceLength/1000).toInt()
        val finalLink = getPairwiseSegment(bamScanner)
        finalLink.filter { it.first < it.second }.forEach { (k, v, ds) -> reportWriter.write("$k\t$v\t$ds\n") }
        reportWriter.close()
        return reportFile
    }
}