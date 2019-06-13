package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.file
import org.tc33.jheatchart.HeatChart
import java.io.File

class CAST(private val args: Array<String>) : CliktCommand() {
    private val output by option("-o", help="output directory [CAST_output]").default("CAST_output")
    private val heatmap by option(help="output distance matrix heatmap [None]").flag()
    private val draftFile by argument(help="draft assembly (usually parent) to be improved").file(true)
    private val bamFiles by argument(help=".bam files (usually progeny) to provide genetic information, one file for each individual")
        .multiple(true)
    private val coverage by lazy{ DepthDetector(bamFiles).getCoverage }
    private val bamScanner by lazy{ BamFilesParallelScanner(bamFiles) }
    private val decay = 0.99995

    private val outDir by lazy { File(output) }
    private val heatmapDir by lazy { File(outDir, "heatmap") }
    private val reportFile by lazy { File(outDir, "report") }
    private val reliablePlot by lazy { File(outDir, "reliable.bed") }

    private fun support(a: Array<DoubleArray>, b: Array<DoubleArray>) =
        a.indices.map { i ->
            val q = a[i].indices.filter{it != i}.map{a[i][it]}
            val p = b[i].indices.filter{it != i}.map{b[i][it]}
            Math.abs(Util.PearsonCorrelationSimilarity(q, p))
        }.sum() / a.size
    private fun getOneChromGVFromScanner(): List<Pair<String, Pair<Matrix, Matrix>?>> {
        // scan this contig to get all snp sites
        val contig = bamScanner.getChrom()
        val contiglen = bamScanner.header.getSequence(contig).sequenceLength
        var sites = ArrayList<Pair<Int, GenotypeVector>> (contiglen / 1000)
        val bwseq = StringBuilder()
        do {
            val gv = GenotypeVector(bamScanner.get())
            gv.vector.indices.forEach {i -> gv.vector[i].checksum(coverage[i]*0.1)}
            if (gv.isReliable && gv.isHeterogeneous()) {
                bwseq.append("$contig\t${bamScanner.pos-1}\t${bamScanner.pos}\n")
                sites.add(Pair(bamScanner.pos, gv))
//                if (contig=="scf7180000001717|quiver")
//                    HeatChart(gv.toDistMatrix()).saveToFile(File(heatmapDir,"${contig}_${bamScanner.pos}.png"))
            }
        } while (bamScanner.nextPosition())
        reliablePlot.writeText(bwseq.toString())
        sites = sites.filterTo(ArrayList()){ (_, gv) -> gv.consistentWithDepth(coverage) }
        if (sites.size < 1) {
            return  mutableListOf<Pair<String, Pair<Matrix, Matrix>>>()
        } else {
            println("\n[$contig] SNP: ${sites.size} sites")
        }
        // train matrix
        val rightScanSum = Matrix(bamFiles.size)
        val leftScanSum = Matrix(bamFiles.size)
        val matrix = Matrix(bamFiles.size)
        val snapshot = ArrayList<Matrix?>(sites.size)
        val segment = ArrayList<Pair<Pair<Pair<Int, Int>, Pair<Int, Int>>, Pair<Matrix, Matrix>?>>() //TODO: type alias
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
            println("\n[$contig] SNP: some sample is missing across whole contig")
            for (i in matrix.data.indices)if (matrix.data[i].all {it.isNaN()}) print("sample $i is missing\n")
            return ArrayList<Pair<String, Pair<Matrix, Matrix>>>()
        }
        // right to left scan
        lastpos = contiglen + 1
        matrix.fill(Double.NaN)
        var lastSegEnd = sites.size - 1
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
                        val seg = Pair(Pair(Pair(sites[i - 1].first, sites[i].first), Pair(
                        if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
                        if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first)),
                        Pair(rightScanSum.copy(), snapshot[lastSegEnd]!!))
                        segment.add(seg)
                        lastSegEnd = i - 1
                    }
                }
            }
            lastpos = pos
        }
        val seg = Pair(Pair(Pair(1, 1), Pair(
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd].first,
            if (lastSegEnd == sites.size - 1) contiglen else sites[lastSegEnd + 1].first)),
        Pair(rightScanSum.copy(), snapshot[lastSegEnd]!!))
        segment.add(seg)
        // merge non-solid segments (100% flexible end)
        for (i in 1 until segment.size - 1) {
            val (j, k, l) = Triple(segment[i - 1], segment[i], segment[i + 1])
            if (j.first.first == k.first.second && k.first.first == l.first.second) {
                segment[i + 1] = Pair(Pair(l.first.first, k.first.second), Pair(l.second!!.first, k.second!!.second)) //TODO: Bug check!
                segment[i] = Pair(k.first,null)
            }
        }
        val named_segment = segment.filter{it != null}.map{x -> Pair(
            "$contig(${x.first.first.first},${x.first.first.second})(${x.first.second.first},${x.first.second.second})", x.second)} // TODO: Bug check!
        if (heatmap) {
            if (!heatmapDir.exists()) heatmapDir.mkdirs()
            for ((name, headMatrixAndTailMatrix) in named_segment) {
                HeatChart(headMatrixAndTailMatrix!!.first.data).saveToFile(File(heatmapDir,"${name}_L.png"))
                HeatChart(headMatrixAndTailMatrix.second.data).saveToFile(File(heatmapDir,"${name}_R.png"))
            }
        }
        return named_segment
    }
    private fun getAllGVFromScanner(): MutableMap<String, Pair<Matrix, Matrix>> {
        val GVmap = mutableMapOf<String, Pair<Matrix, Matrix>>()
        while (bamScanner.hasNext()) {
            val CGV = getOneChromGVFromScanner()
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

    override fun run() {
        outDir.mkdirs()
        if (reportFile.exists()) reportFile.delete()
        val reportWriter = reportFile.bufferedWriter()
        reportWriter.write(args.joinToString(" ","cmd: ","\n"))
        println(coverage.map { it.toInt() }.joinToString("|","\nDepth:"))
        val contigsTwoEndGV = getAllGVFromScanner()
        val contigsOneEndGV = contigsTwoEndGV .flatMap { (contig, v) -> listOf(Pair("+$contig", v.first), Pair("-$contig", v.second)) }
        val finalLink = pairwiseMutualBest(contigsOneEndGV)
        for ((k, v) in finalLink) if (k < v.first) {
            reportWriter.write("$k\t${v.first}\t${v.second}\n")
        }
        reportWriter.close()
        //======
        val stringer = Stringer(draftFile, outDir)
        stringer.loadEdgeFromFile(reportFile.absolutePath)
        stringer.correctAndScaffoldFasta()
    }
}
fun main(args: Array<String>) = CAST(args).main(if (args.isEmpty()) arrayOf("--help") else args)