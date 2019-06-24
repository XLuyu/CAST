package edu.nus

import org.tc33.jheatchart.HeatChart
import java.io.File

class Genotyper(private val bamFiles:List<String>, private val outDir: File, private val heatmap:Boolean=false){
    private val decay = 0.99995
    private fun support(a: Array<DoubleArray>, b: Array<DoubleArray>) =
        a.indices.map { i ->
            val q = a[i].indices.filter{it != i}.map{a[i][it]}
            val p = b[i].indices.filter{it != i}.map{b[i][it]}
            Math.abs(Util.PearsonCorrelationSimilarity(q, p))
        }.sum() / a.size
    private fun getOneChromGVFromScanner(bamScanner: BamFilesParallelScanner): List<Pair<String, MatMat?>> {
        // scan this contig to get all snp sites
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
//        reliablePlot.writeText(bwseq.toString())
        if (sites.size < 1) return  mutableListOf<Pair<String, MatMat>>()
        println("\n[$contig] SNP: ${sites.size} sites")
        // train matrix
        val rightScanSum = Matrix(bamFiles.size)
        val leftScanSum = Matrix(bamFiles.size)
        val matrix = Matrix(bamFiles.size)
        val snapshot = ArrayList<Matrix?>(sites.size)
        // left to right scan
        var lastpos = 0
        matrix.fill(Double.NaN)
        for ((pos, gv) in sites) {
            matrix.updateBy(Matrix(gv.toDistMatrix()))
            if (!matrix.containsNaN()) {
                val decayN = Math.pow(decay, (pos - lastpos).toDouble())
                leftScanSum.timesUpdate(decayN)
                leftScanSum += (matrix * ((1 - decayN) / (1 - decay)))
                snapshot.add(leftScanSum.copy())
            } else
                snapshot.add(null)
            lastpos = pos
        }
        if (snapshot.last() == null) {
            println("\n[Info] No reliable genotype on $contig from: ")
            for (i in matrix.data.indices) if (matrix.data[i].all {it.isNaN()}) println(bamFiles[i])
            return ArrayList<Pair<String, MatMat>>()
        }
        // right to left scan
        lastpos = contiglen + 1
        matrix.fill(Double.NaN)
        var lastSegEnd = sites.size - 1
        val segment = ArrayList<Pair<Array<Int>, MatMat>?>()
        for (i in sites.indices.reversed()) {
            val (pos, gv) = sites[i]
            matrix.updateBy(Matrix(gv.toDistMatrix()))
            if (!matrix.containsNaN()) {
                val decayN = Math.pow(decay, (lastpos - pos).toDouble())
                rightScanSum.timesUpdate(decayN)
                rightScanSum += (matrix * ((1 - decayN) / (1 - decay)))
                if (0 < i && snapshot[i - 1] != null) {
                    val matrixPearson = support(snapshot[i - 1]!!.data, rightScanSum.data)
//          FileUtils.write(conf.geneticPlot,f"$matrixPearson\t${sites(i).first-sites(i-1).first}\n",true)
                    if (matrixPearson < 0.3) {
                        val seg = Pair(arrayOf(sites[i - 1].first, sites[i].first,
                            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
                            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first),
                            Pair(rightScanSum.copy(), snapshot[lastSegEnd]!!))
                        segment.add(seg)
                        lastSegEnd = i - 1
                    }
                }
            }
            lastpos = pos
        }
        val seg = Pair(arrayOf(1, 1,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first),
            Pair(rightScanSum.copy(), snapshot[lastSegEnd]!!))
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
        val named_segment = segment.filterNotNull().map{x ->
            Pair("$contig(${x.first[0]},${x.first[1]})(${x.first[2]},${x.first[3]})", x.second)}
        if (heatmap) {
            val heatmapDir = File(outDir, "heatmap")
            if (!heatmapDir.exists()) heatmapDir.mkdirs()
            for ((name, headMatrixAndTailMatrix) in named_segment) {
                HeatChart(headMatrixAndTailMatrix.first.data).saveToFile(File(heatmapDir,"${name}_L.png"))
                HeatChart(headMatrixAndTailMatrix.second.data).saveToFile(File(heatmapDir,"${name}_R.png"))
            }
        }
        return named_segment
    }
    private fun getAllGVFromScanner(bamScanner: BamFilesParallelScanner): MutableMap<String, MatMat> {
        val GVmap = mutableMapOf<String, MatMat>()
        while (bamScanner.hasNext()) {
            val CGV = getOneChromGVFromScanner(bamScanner)
            for ((contig, headMatrixAndTailMatrix) in CGV)
                GVmap[contig] = headMatrixAndTailMatrix!!
        }
        return GVmap
    }
    private fun pairwiseMutualBest(contigsGV: List<Pair<String, Matrix>>): ArrayList<Pair<String, Pair<String, Double>>> {
        val candidate = ArrayList<Pair<String, Pair<String, Double>>>()
        for ((k, v) in contigsGV) {
            val contig = k.substring(1)
            val similarities = contigsGV.filter{ it.first.substring(1) != contig}.map{Pair(it.first, support(v.data, it.second.data))}
            similarities.filter{it.second > 0.6}.forEach{candidate.add(Pair(k, Pair(it.first,it.second)))}
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